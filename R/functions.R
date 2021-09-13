#' @importFrom stats rbinom rgamma rnorm
#' @importFrom utils data globalVariables
#' @importFrom BayesLogit rpg

draw_Omega = function(data, Alpha, beta, aggregation) {
  list2env(data, envir = environment())
  b = rep(nn, ncol(A))
  c = as.vector(tcrossprod(X, Alpha) + V %*% t(beta))
  if (aggregation) {
    Omega = rpg(length(b), b, c)
  } else {
    Omega = pgdraw(b, c)
  }
  Omega = matrix(Omega, nrow(X), nrow(Alpha))
  return(Omega)
}

# draw_Alpha = function(data, beta,  Omega,  sigma_alpha_2) {
#   #browser()
#   list2env(data, envir = environment())
#   Alpha = matrix(0, length(beta), ncol(X))
#   for (j in 1:length(beta)) {
#     Sigma_alpha_j = crossprod(X, Omega[, j] * X) + diag(ncol(X)) / sigma_alpha_2
#     Sigma_alpha_j = solve(Sigma_alpha_j)
#
#     mu_alpha_j = crossprod(X, A[, j] - nn / 2 - beta[j] * Omega[, j] * V)
#     mu_alpha_j = Sigma_alpha_j %*% mu_alpha_j
#
#     Alpha[j,] = rmvnorm(1, mean = mu_alpha_j, sigma = Sigma_alpha_j)
#   }
#   return(Alpha)
# }

draw_Alpha = function(data, beta,  Omega,  sigma_alpha_2) {
  #browser()
  list2env(data, envir = environment())
  Alpha = matrix(0, length(beta), ncol(X))
  for (j in 1:length(beta)) {
    pre_alpha_j = crossprod(X, Omega[, j] * X) + diag(ncol(X)) / sigma_alpha_2
    mu_alpha_j = crossprod(X, A[, j] - nn / 2 - beta[j] * Omega[, j] * V)

    Alpha[j, ] = rmvnorm_chol(mu_alpha_j, pre_alpha_j)
  }
  return(Alpha)
}

draw_beta = function(beta,
                     data,
                     Alpha,
                     Omega,
                     sigma_beta,
                     delta,
                     epsilon) {
  list2env(data, envir = environment())

  beta_old = beta
  Lbeta = L %*% beta
  for (j in 1:length(beta)) {
    sigma_betaj_2 = sigma_beta[j] ^ 2 * delta[j] *
      sum(Omega[, j] * V ^ 2) + (L[j, j] + epsilon)
    sigma_betaj_2 = 1 / sigma_betaj_2

    mu_betaj = sigma_beta[j] * delta[j] *
      sum((A[, j] - nn / 2 - (X %*% Alpha[j, ]) * Omega[, j]) * V) -
      (sum(Lbeta[-j]) - sum(L[-j, j] * beta[j]))
    # should be equal to sum(L[-j,-j] %*% beta[-j])
    mu_betaj = sigma_betaj_2 * mu_betaj

    beta[j] = rnorm(1, mu_betaj, sqrt(sigma_betaj_2))
    Lbeta = Lbeta + L[, j] * (beta[j] - beta_old[j])
  }
  bound = 100
  beta[abs(beta) > bound] = bound * sign(beta[abs(beta) > bound])
  return(beta)
}

draw_sigma_beta = function(data,
                           beta,
                           Alpha,
                           Omega,
                           delta,
                           tau_2) {
  list2env(data, envir = environment())
  sigma_beta = rep(0, length(beta))

  for (j in 1:length(beta)) {
    sigmaj_2 = beta[j] ^ 2 * delta[j] * sum(Omega[, j] * V ^ 2) + 1 / tau_2
    sigmaj_2 = 1 / sigmaj_2

    muj = beta[j] * delta[j] *
      sum((A[, j] - nn / 2 - (X %*% Alpha[j,]) * Omega[, j]) * V)
    muj = sigmaj_2 * muj

    sigma_beta[j] = rtruncnorm(1,
                               a = 0,
                               mean = muj,
                               sd = sqrt(sigmaj_2))
  }
  return(sigma_beta)
}

draw_delta = function(data, Alpha, beta, Omega, pi_delta, gammaG = NULL) {
  list2env(data, envir = environment())
  delta = rep(0, length(beta))
  for (j in 1:length(beta)) {
    D = beta[j] * sum((A[, j] - nn / 2 - (X %*% Alpha[j,]) * Omega[, j]) * V) -
      0.5 * beta[j] ^ 2 * sum(Omega[, j] * V ^ 2)
    if (!is.null(gammaG)) {
      D = D + sum(G[j,] * gammaG)
    }
    if (pi_delta == 1) {
      p = 1
    } else {
      p = 1 / (1 + (1 - pi_delta) / pi_delta * exp(-D))
    }
    delta[j] = rbinom(1, 1, p)
  }
  return(delta)
}

draw_gammaG = function(data,
                       omega,
                       gammaG,
                       sigma_gammaG,
                       delta) {
  list2env(data, envir = environment())
  for (k in 1:ncol(G)) {
    pre = sum(omega * G[, k] ^ 2) * sigma_gammaG[k] ^ 2 + 1
    mu = sum((delta - 0.5) * G[, k]) -
      sum(gammaG[-k] * sigma_gammaG[-k] * colSums(omega * G[, k] * G[, -k]))
    mu = mu * sigma_gammaG[k] / pre
    gammaG[k] = rnorm(1, mu, 1 / sqrt(pre))
  }
  return(gammaG)
}

draw_sigma_gammaG = function(data,
                             omega,
                             gammaG,
                             sigma_gammaG,
                             delta,
                             tau_gammaG_2) {
  list2env(data, envir = environment())
  for (k in 1:ncol(G)) {
    pre = sum(omega * G[, k] ^ 2) * gammaG[k] ^ 2 + 1 / tau_gammaG_2
    mu = sum((delta - 0.5) * G[, k]) -
      sum(gammaG[-k] * sigma_gammaG[-k] * colSums(omega * G[, k] * G[,-k]))
    mu = mu * gammaG[k] / pre
    sigma_gammaG[k] = rtruncnorm(1,
                                 a = 0,
                                 mean = mu,
                                 sd = 1 / sqrt(pre))
  }

  return(sigma_gammaG)
}

rmvnorm_chol = function(mu, preMat) {
  R = chol(preMat)
  b = solve(t(R), mu)
  Z = rnorm(length(b))
  X = solve(R, Z + b)
  return(X)
}
