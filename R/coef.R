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

  chain = chain[chain_range]
  chain = lapply(chain, one_ite_summ)
  coeflst = lst_mean(chain)

  return(coeflst)
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
    if (var == 'beta_conditional') {
      varMat = sapply(lsts, function(x)
        x[[var]])
      mean_lst[[var]] = rowMeans(varMat,na.rm = TRUE)
    } else {
      for (j in 1:n_lsts) {
        mean_lst[[var]] =
          mean_lst[[var]] * ((j - 1) / j) +
          lsts[[j]][[var]] / j
      }
    }
  }
  return(mean_lst)
}
