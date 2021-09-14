#' Update the BV chain
#'
#' @param object An object of s3 class BV_chain
#' @param n_ite The number of new MCMC draws to add to the currrent chain
#' @param max_length The maximum length of the chain to be stored
#' @param ... more arguments
#'
#' @return The updated BV_chain
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @exportS3Method update BV_chain

update.BV_chain = function(object,
                           n_ite = 1,
                           max_length = Inf ,
                           ...) {
  list2env(object, envir = environment())

  # extract the latest parameters
  coeflst = coef(object, chain_range = end_pos)
  coeflst$Omega = chain$Omega

  # resize the chain to appropriate length
  chain_length = end_pos - start_pos + 1
  n_drop = chain_length + n_ite - max_length
  n_drop = max(0, n_drop)
  n_drop = min(chain_length, n_drop)
  n_add = min(n_ite, max_length)
  chain = resize_chain(chain, n_drop, n_add)

  chain_length = chain_length + n_add - n_drop
  end_pos = end_pos + n_ite
  start_pos = end_pos - chain_length + 1
  # Update
  pb <- txtProgressBar(min = 0, max = n_ite, style = 3)
  for (j in 1:n_ite) {
    coeflst = update_one_time(coeflst, object)
    pos = end_pos - start_pos  + 1 + j - n_ite
    if (pos >= 1) {
      chain = fill_in(coeflst, chain, pos)
    }
    # update progression bar
    setTxtProgressBar(pb, j)
  }
  close(pb)
  # store the updated chain

  object$chain = chain
  object$end_pos = end_pos
  object$start_pos = start_pos
  return(object)
}
##################################
resize_chain = function(chain, n_drop, n_add) {
  chain_names = names(chain)[grepl('_chain', names(chain))]
  chain_length = ncol(chain[['beta_chain']])
  if (n_drop - n_add >= chain_length) {
    chain_range = numeric(0)
    n_add = 0
  } else if (n_drop < chain_length) {
    chain_range = (n_drop + 1):(chain_length + min(n_add, 0))
    n_add = max(n_add, 0)
  } else {
    chain_range = numeric(0)
    n_add = n_add - n_drop + chain_length
  }
  for (chain_name in chain_names) {
    chain_dim = dim(chain[[chain_name]])
    if (is.null(chain_dim)) {
      chain[[chain_name]] = chain[[chain_name]][chain_range]
      chain[[chain_name]] = c(chain[[chain_name]], rep(0, n_add))
    } else if (length(chain_dim) == 2) {
      chain[[chain_name]] =
        chain[[chain_name]][, chain_range, drop = FALSE]
      chain[[chain_name]] = cbind(chain[[chain_name]],
                                  matrix(0, nrow(chain[[chain_name]]), n_add))
    } else if (length(chain_dim) == 3) {
      chain[[chain_name]] =
        chain[[chain_name]][, , chain_range, drop = FALSE]
      dim1 = dim(chain[[chain_name]])
      n0 = prod(dim1[1:2]) * n_add
      chain[[chain_name]] = c(chain[[chain_name]], rep(0, n0))
      chain[[chain_name]] = array(chain[[chain_name]],
                                  dim = c(dim1[1:2], dim1[3] + n_add))
    }
  }
  return(chain)
}
##################################
fill_in = function(coeflst, chain, pos) {
  var_names = gsub('_chain', '', names(chain))
  for (var_name in var_names) {
    chain_name = paste0(var_name, '_chain')
    chain_dim = dim(chain[[chain_name]])
    if (is.null(chain_dim)) {
      chain[[chain_name]][pos] = coeflst[[var_name]]
    } else if (length(chain_dim) == 2) {
      chain[[chain_name]][, pos] = coeflst[[var_name]]
    } else if (length(chain_dim) == 3) {
      chain[[chain_name]][, , pos] = coeflst[[var_name]]
    }
  }
  return(chain)
}
##################################
update_one_time = function(coeflst, object) {
  list2env(coeflst, envir = environment())
  list2env(object, envir = environment())
  list2env(hyperparameters, envir = environment())
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
                   delta,
                   epsilon)
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
