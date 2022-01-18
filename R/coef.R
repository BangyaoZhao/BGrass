#' Extract the logLik of the BV_chain at any position(s)
#'
#' @param object An object of s3 class BV_chain
#' @param chain_range Which part of the chain should be used, by default it will use the later 80\%
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
#' Combine multiple BV chains
#'
#' @param BV_chains A list of objects of s3 class BV_chain
#'
#' @return An object of s3 class BV_chains
#'
as.BV_chains = function(BV_chains) {
  class(BV_chains) = 'BV_chains'
  return(BV_chains)
}
#####################################
coef.BV_chains = function(object,
                          chain_range = NULL,
                          ...) {
  list2env(object[[1]], envir = environment())
  if (is.null(chain_range)) {
    chain_range =  max(ceiling(0.2 * end_pos), start_pos):end_pos
  }
  chain_range = chain_range[chain_range %in% (start_pos:end_pos)]
  chain_range = chain_range + 1 - start_pos

  chain = list()
  for (BV_chain in object) {
    chain = c(chain, BV_chain$chain[chain_range])
  }
  chain = lapply(chain, one_ite_summ)
  coef_lst = lst_mean(chain)
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
      mean_lst[[var]] = rowMeans(varMat, na.rm = TRUE)
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
