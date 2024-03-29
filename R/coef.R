#' Extract the logLik of the BGrass_chain at any position(s)
#'
#' @param object An object of s3 class BGrass_chain
#' @param chain_range Which part of the chain should be used, by default it will use the later 80\%
#' @param ... more arguments
#'
#' @return A list containing the estimated values of all the parameters in the chain.
#'
#' @exportS3Method coef BGrass_chain
#' @exportS3Method coef BGrass_chains
#' @export as.BGrass_chains
#'

coef.BGrass_chain = function(object,
                         chain_range = NULL,
                         ...) {
  list2env(object, envir = current_env())
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
#' Combine multiple BGrass chains
#'
#' @param BGrass_chains A list of objects of s3 class BGrass_chain
#'
#' @return An object of s3 class BGrass_chains
#'
as.BGrass_chains = function(BGrass_chains) {
  class(BGrass_chains) = 'BGrass_chains'
  return(BGrass_chains)
}
#####################################
coef.BGrass_chains = function(object,
                          chain_range = NULL,
                          ...) {
  list2env(object[[1]], envir = current_env())
  if (is.null(chain_range)) {
    chain_range =  max(ceiling(0.2 * end_pos), start_pos):end_pos
  }
  chain_range = chain_range[chain_range %in% (start_pos:end_pos)]
  chain_range = chain_range + 1 - start_pos

  chain = list()
  for (BGrass_chain in object) {
    chain = c(chain, BGrass_chain$chain[chain_range])
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
