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
  data = aggregate_data(object$data)
  list2env(object$chain, envir = environment())
  chain_range = chain_range - object$start_pos + 1
  #browser()
  l_chain = sapply(chain_range, function(i) {
    l_calcu(delta_chain[, i] * sigma_beta_chain[, i] * beta_chain[, i],
            Alpha_chain[, , i],
            data)
  })
  return(l_chain)
}

l_calcu = function(beta, Alpha, data) {
  list2env(data, envir = environment())
  P = 1 / (1 + exp(-(
    tcrossprod(X, Alpha) + tcrossprod(V, beta)
  )))
  l = sum(A * log(P) + (nn - A) * log(1 - P))
  return(l)
}
