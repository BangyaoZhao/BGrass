#' Extract the logLik of the BV_chain at any position(s)
#'
#' @param object An object of s3 class BV_chain
#' @param chain_range Which part of the chain should be used
#' @param raw_beta Should the raw beta (or the marginal beta) be calculated
#' @param ... more arguments
#'
#' @return A list containing the estimated values of all the parameters in the chain.
#'
#' @exportS3Method coef BV_chain
#' @exportS3Method coef BV_chains
#' @export as.BV_chains
#'

coef.BV_chain = function(object,
                         chain_range = NULL,
                         ...) {
  list2env(object, envir = environment())
  if (is.null(chain_range)) {
    chain_range =  max(ceiling(0.2 * end_pos), start_pos):end_pos
  }
  chain_range = chain_range[chain_range %in% (start_pos:end_pos)]
  chain_range = chain_range + 1 - start_pos


  coef_lst = list()
  var_names = gsub('_chain', '', names(chain))
  for (var in var_names) {
    var_chain = chain[[paste0(var, '_chain')]]
    if (is.null(dim(var_chain))) {
      var_esti = mean(var_chain[chain_range])
    } else if (length(dim(var_chain)) == 2) {
      var_esti = rowMeans(var_chain[, chain_range, drop = FALSE])
    } else if (length(dim(var_chain)) == 3) {
      var_esti = rowMeans(var_chain[, , chain_range, drop = FALSE], dims = 2)
    }

    coef_lst[[var]] = var_esti
  }

  # marginalized beta
  sigma_beta_chain = chain$sigma_beta_chain[, chain_range, drop = FALSE]
  delta_chain = chain$delta_chain[, chain_range, drop = FALSE]
  beta_chain = chain$beta_chain[, chain_range, drop = FALSE]

  beta_marginalized = sapply(1:nrow(beta_chain), function(j) {
    # browser()
    betaj_chain = beta_chain[j, ] *
      sigma_beta_chain[j, ] * delta_chain[j, ]
    mean(betaj_chain)
  })
  coef_lst[['beta_marginalized']] = beta_marginalized

  # marginalized gammaG
  if (enrichment) {
    gammaG_chain = chain$gammaG_chain[, chain_range, drop = FALSE]
    sigma_gammaG_chain = chain$sigma_gammaG_chain[, chain_range, drop = FALSE]
    gammaG_marginalized = rowMeans(gammaG_chain * sigma_gammaG_chain)
    coef_lst[['gammaG_marginalized']] = gammaG_marginalized
  }

  return(coef_lst)
}
#####################################
as.BV_chains = function(BV_chains) {
  class(BV_chains) = 'BV_chains'
  return(BV_chains)
}
#####################################
coef.BV_chains = function(object,
                          chain_range = NULL,
                          ...) {
  n_chains = length(object)
  coef_lsts = lapply(object, function(BV_chain) {
    coef_lst = coef(BV_chain, chain_range = chain_range)
    return(coef_lst)
  })

  coef_lst = lst_mean(coef_lsts)
  return(coef_lst)
}
#####################################
lst_mean = function(lsts) {
  n_lsts = length(lsts)
  mean_lst = lsts[[1]]
  for (var in names(mean_lst)) {
    for (j in 1:n_lsts) {
      if (j > 1) {
        mean_lst[[var]] =
          mean_lst[[var]] + lsts[[j]][[var]]
      }
    }
    mean_lst[[var]] = mean_lst[[var]] / n_lsts
  }
  return(mean_lst)
}
