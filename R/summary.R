#' Extract the logLik of the BV_chain at any position(s)
#'
#' @param object An object of s3 class BV_chain
#' @param chain_range Which part of the chain should be used
#' @param raw_beta Should the raw beta (or the marginal beta) be calculated
#' @param ... more arguments
#'
#' @return A list containing the estimated values of all the parameters in the chain.
#'
#' @exportS3Method summary BV_chains
#' @exportS3Method summary BV_chain
#'
summary.BV_chains = function(object,
                             chain_range = NULL,
                             ...) {
  start_pos = object[[1]]$start_pos
  end_pos = object[[1]]$end_pos
  if (is.null(chain_range)) {
    chain_range =  max(ceiling(0.2 * end_pos), start_pos):end_pos
  }
  chain_range = chain_range - start_pos + 1
  data = object[[1]]$data
  ###########
  beta_Mat = sapply(object, function(BV_chain) {
    chain = BV_chain$chain
    beta_Mat = chain$beta_chain[, chain_range]
    delta_Mat = chain$delta_chain[, chain_range]
    sigma_beta_Mat = chain$sigma_beta_chain[, chain_range]
    beta_marginalized_Mat =
      beta_Mat * delta_Mat * sigma_beta_Mat
    return(beta_marginalized_Mat)
  },
  simplify = "array")
  beta_Mat = matrix(beta_Mat, nrow(beta_Mat))
  ###########
  delta_Mat = sapply(object, function(BV_chain) {
    chain = BV_chain$chain
    delta_Mat = chain$delta_chain[, chain_range]
    return(delta_Mat)
  },
  simplify = "array")
  delta_Mat = matrix(delta_Mat, nrow(delta_Mat))

  beta_table = matrix(0, nrow(beta_Mat), 5)
  beta_table[, 1] = rowMeans(beta_Mat)
  beta_table[, 2] = rowMeans(delta_Mat)
  apply(beta_Mat, 1, function(x) {
    quantile(x, c(0.05, 0.95), names = FALSE)
  }) %>%
    t() ->
    beta_table[, 3:4]
  apply(beta_Mat, 1, function(x) {
    mean(x <= 0)
  }) ->
    beta_table[, 5]
  rownames(beta_table) = rownames(object[[1]]$data$L)
  colnames(beta_table) = c('beta', 'delta', 'lower95', 'upper95', 'p')
  res = list(chain_range = chain_range,
             beta_table = beta_table)
  ##############################
  coef_lst = coef(object, chain_range)
  l = l_calcu(coef_lst$beta_marginalized, coef_lst$Alpha, data)
  ls = sapply(object, function(BV_chain) {
    coef_lst = coef(BV_chain, chain_range)
    l = l_calcu(coef_lst$beta_marginalized,
                coef_lst$Alpha, data)
    l
  })
  l_Mat = sapply(object, function(BV_chain) {
    logLik(BV_chain, chain_range)
  })
  D1 = -2 * l
  D2 = -2 * mean(l_Mat)
  DIC = 2*D2 - D1
  res$lik = list(all_chains = l,
               each_chain = ls,
               trace = l_Mat,
               DIC = DIC)
  ##############################
  if (object[[1]]$enrichment) {
    gammaG_Mat = sapply(object, function(BV_chain) {
      chain = BV_chain$chain
      gammaG_Mat = chain$gammaG_chain[, chain_range]
      sigma_gammaG_Mat = chain$sigma_gammaG_chain[, chain_range]
      gammaG_marginalized_Mat =
        gammaG_Mat * sigma_gammaG_Mat
      return(gammaG_marginalized_Mat)
    },
    simplify = "array")
    gammaG_Mat = matrix(gammaG_Mat, nrow(gammaG_Mat))

    #########
    gammaG_table = matrix(0, nrow(gammaG_Mat), 4)
    gammaG_table[, 1] = rowMeans(gammaG_Mat)
    apply(gammaG_Mat, 1, function(x) {
      quantile(x, c(0.05, 0.95), names = FALSE)
    }) %>%
      t() ->
      gammaG_table[, 2:3]
    apply(gammaG_Mat, 1, function(x) {
      mean(x <= 0)
    }) ->
      gammaG_table[, 4]
    rownames(gammaG_table) = colnames(object[[1]]$data$G)
    colnames(gammaG_table) = c('gammaG', 'lower95', 'upper95', 'p')
    res$gammaG_table = gammaG_table
  }

  return(res)
}
##########################################
summary.BV_chain = function(object,
                            chain_range = NULL,
                            ...) {
  BV_chains = as.BV_chains(list(object))
  res = summary(BV_chains, chain_range)
  return(res)
}
