#' Extract the logLik of the BV_chain at any position(s)
#'
#' @param object An object of s3 class BV_chain
#' @param chain_range The range of the chain to be evaluated
#' @param ... more arguments
#'
#' @return a vector of length \code{length(chain_range)} representing the
#' corresponding log likelihood.
#'
#' @exportS3Method logLik BV_chain
#'

logLik.BV_chain = function(object,
                           chain_range = object$start_pos:object$end_pos,
                           ...) {
  list2env(object, envir = environment())
  chain_range = chain_range - start_pos + 1
  #browser()
  l_chain = sapply(chain, function(coeflst) {
    coeflst = one_ite_summ(coeflst)
    list2env(coeflst, envir = environment())
    l_calcu(beta_marginalized, Alpha, data)
  })
  return(l_chain)
}
##########################################
l_calcu = function(beta, Alpha, data) {
  list2env(data, envir = environment())
  P = 1 / (1 + exp(-(
    tcrossprod(X, Alpha) + tcrossprod(V, beta)
  )))
  l = sum(A * log(P) + (nn - A) * log(1 - P))
  return(l)
}
##########################################
one_ite_summ = function(coeflst) {
  list2env(coeflst, envir = environment())
  beta_marginalized = beta * sigma_beta * delta
  beta_conditional = beta * sigma_beta
  beta_conditional[!as.logical(delta)] = NA

  coeflst$beta_marginalized = beta_marginalized
  coeflst$beta_conditional = beta_conditional

  if ('gammaG' %in% names(coeflst)) {
    gammaG_marginalized = gammaG * sigma_gammaG
    coeflst$gammaG_marginalized = gammaG_marginalized
  }
  return(coeflst)
}
