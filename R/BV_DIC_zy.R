BV_DIC <- function(object,
                   chain_range = object$start_pos:object$end_pos,
                   ...){
  list2env(object, envir = environment())
  list2env(data, envir = environment())
  chain_range = chain_range - start_pos + 1
  chain = chain[chain_range]
  mcmc_beta <-  sapply(chain, function(x)x$beta)
  mcmc_delta <- sapply(chain, function(x)x$delta)
  mcmc_sigma_beta <- sapply(chain, function(x)x$sigma_beta)
  mcmc_alpha <- lapply(chain, function(x)x$Alpha)
  
  #mcmc_loglik
  mcmc_loglik = t(sapply(1:ncol(mcmc_beta), function(i) {
    beta = mcmc_delta[,i] * mcmc_sigma_beta[,i] * mcmc_beta[,i]
    # P = 1 / (1 + exp(-(
    #   tcrossprod(X, mcmc_alpha[[i]]) + tcrossprod(V, beta)
    # )))
    # l = colSums(A * log(P) + (nn - A) * log(1 - P))
    # return(l)
    
    # alternate way of loglik
    # for individual data
    lik_temp = tcrossprod(X, mcmc_alpha[[i]]) + tcrossprod(V, beta)
    loglik = colSums(A * lik_temp) - colSums(nn*log1pexp(lik_temp))
    return(loglik)
  }))

  
  DIC_mat  = -4*mean(apply(mcmc_loglik,1, function(s) sum(s))) +
    2*sum(colMeans(mcmc_loglik))
  
  return(DIC_mat)
}
