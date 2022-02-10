#' Update the BGrass chain
#'
#' @param object An object of s3 class BGrass_chain
#' @param n_ite The number of new MCMC draws to add to the currrent chain
#' @param max_length The maximum length of the chain to be stored
#' @param ... more arguments
#'
#' @return The updated BGrass_chain
#'
#' @importFrom progress progress_bar
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom rlang current_env
#' @exportS3Method update BGrass_chain

update.BGrass_chain = function(object,
                           n_ite = 1,
                           max_length = Inf ,
                           ...) {
  list2env(object, envir = current_env())
  # extend the logl
  newlogl = rep(0,n_ite*n_thin)

  # extract the latest parameters
  chain_length = end_pos - start_pos + 1
  coeflst = chain[[chain_length]]

  # resize the chain to appropriate length
  n_drop = chain_length + n_ite - max_length
  n_drop = max(0, n_drop)
  n_drop = min(chain_length, n_drop)
  if (n_drop < chain_length) {
    chain = chain[(n_drop + 1):chain_length]
  } else {
    chain = list()
  }

  end_pos = end_pos + n_ite
  start_pos = max(end_pos - max_length + 1, start_pos)
  # Update
  isRStudio <- Sys.getenv("RSTUDIO") == "1"
  # select appropriate progression bar
  if (!isRStudio) {
    pb <- txtProgressBar(min = 0,
                         max = n_ite,
                         style = 3)
  } else {
    pb <-
      progress_bar$new(format = "  updating [:bar] :percent eta: :eta",
                       total = n_ite)
  }
  for (j in 1:n_ite) {
    for (i in 1:n_thin) {
      coeflst = update_one_time(coeflst, object)
      coeflst = one_ite_summ(coeflst)
      newlogl[(j - 1) * n_thin + i] =
        full_l_calcu(data,
                     coeflst,
                     hyperparameters)
    }
    pos = end_pos - start_pos  + 1 + j - n_ite
    if (pos >= 1) {
      chain[[pos]] = coeflst
    }
    # update progression bar

    if (!isRStudio) {
      setTxtProgressBar(pb, j)
    } else {
      pb$tick()
    }
  }
  if (!isRStudio) {
    close(pb)
  }
  # store the updated chain
  object$logl = c(logl,newlogl)
  object$chain = chain
  object$end_pos = end_pos
  object$start_pos = start_pos
  return(object)
}
##################################
update_one_time = function(coeflst, object) {
  list2env(coeflst, envir = current_env())
  list2env(object, envir = current_env())
  list2env(hyperparameters, envir = current_env())
  list2env(data, envir = current_env())
  J = length(beta)
  p1 = ncol(Alpha)
  new_coeflst = list()
  Omega = draw_Omega(data,
                     Alpha,
                     beta * delta * sigma_beta,
                     aggregation)
  if (enrichment) {
    delta = draw_delta(data,
                       Alpha,
                       beta * sigma_beta,
                       Omega,
                       0.5,
                       gammaG * sigma_gammaG)
    omega = pgdraw(1, G %*% (gammaG * sigma_gammaG))
    gammaG = draw_gammaG(data,
                         omega,
                         gammaG,
                         sigma_gammaG,
                         delta)
    sigma_gammaG = draw_sigma_gammaG(data,
                                     omega,
                                     gammaG,
                                     sigma_gammaG,
                                     delta,
                                     tau_gammaG_2)
    tau_gammaG_2 = 1 / rgamma(1,
                              a_gammaG + ncol(G) / 2,
                              b_gammaG + sum(sigma_gammaG ^ 2) / 2)

    new_coeflst$omega = omega
    new_coeflst$gammaG = gammaG
    new_coeflst$sigma_gammaG = sigma_gammaG
    new_coeflst$tau_gammaG_2 = tau_gammaG_2
  } else {
    delta = draw_delta(data, Alpha, beta * sigma_beta, Omega, pi_delta)
  }
  beta = draw_beta(beta,
                   data,
                   Alpha,
                   Omega,
                   sigma_beta,
                   delta)
  sigma_alpha_2 = 1 / rgamma(ncol(X),
                             a_alpha + J / 2,
                             b_alpha + colSums(Alpha ^ 2) / 2)
  tau_2 = 1 / rgamma(1, a_tau + J / 2, b_tau + sum(sigma_beta ^ 2) / 2)
  Alpha = draw_Alpha(data, beta * delta * sigma_beta,  Omega,  sigma_alpha_2)
  sigma_beta = draw_sigma_beta(data,
                               beta,
                               Alpha,
                               Omega,
                               delta,
                               tau_2)

  new_coeflst$delta = delta
  new_coeflst$beta = beta
  new_coeflst$sigma_beta = sigma_beta
  new_coeflst$sigma_alpha_2 = sigma_alpha_2
  new_coeflst$tau_2 = tau_2
  new_coeflst$Alpha = Alpha

  return(new_coeflst)
}
