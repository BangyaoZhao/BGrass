#' Extract the logLik of the BGrass_chain at any position(s)
#'
#' @param object An object of s3 class BGrass_chain
#' @param chain_range The range of the chain to be evaluated
#' @param ... more arguments
#'
#' @return a vector of length \code{length(chain_range)} representing the
#' corresponding log likelihood.
#'
#' @exportS3Method logLik BGrass_chain
#'

logLik.BGrass_chain = function(object,
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
full_l_calcu = function(data, coeflst, hyperparameters) {
  list2env(coeflst, envir = environment())
  list2env(data, envir = environment())
  list2env(hyperparameters, envir = environment())
  J = length(beta)
  p1 = ncol(Alpha)
  P = 1 / (1 + exp(-(
    tcrossprod(X, Alpha) + tcrossprod(V, beta_marginalized)
  )))
  l = sum(A * log(P) + (nn - A) * log(1 - P))

  l = l - sum(t(Alpha ^ 2) / sigma_alpha_2) / 2 -
    sum(log(sigma_alpha_2)) * J / 2
  l = l - sum(b_alpha / sigma_alpha_2) -
    (a_alpha + 1) * sum(log(sigma_alpha_2))
  l = l - t(beta) %*% L %*% beta / 2
  l = l - J * sum(log(tau_2)) / 2 - sum(sigma_beta ^ 2 / tau_2) / 2
  l = l - (a_tau + 1) * log(tau_2) + b_tau / tau_2
  if (pi_delta < 1) {
    l = l + sum(delta) * log(pi_delta) + (J - sum(delta)) * log(1 - pi_delta)
  }
  l = as.numeric(l)
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
