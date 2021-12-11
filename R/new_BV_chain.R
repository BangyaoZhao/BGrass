#' Initialize a new BV chain
#'
#' @param A A 0-1 matrix with n rows and J columns (J AEs).
#' The i-j th element of A indicates whether the j th AE is reported in report i.
#' @param V A 0-1 vector of length n indicating the type pf vaccine in each report.
#' @param X The covariate matrix of dim n by p + 1 (p covariates and a column of intercept)
#' @param L The J by J normalized graph Laplacian matrix of AEs.
#' You can get it using the function \code{getLaplacian}.
#' @param G The J by K matrix indicating the groups each AE comes from. if \code{!enrichment},
#' this argument should be set to \code{NULL}.
#' You can get it using the function \code{getG}.
#' @param nn If \code{A} is already aggregated, indicate the number of
#' reports corresponding the each row of \code{A}. The function will
#' assume the data is not aggregated, so \code{nn = rep(1, nrow(A))} by default.
#' @param epsilon A hyperparameter, controlling the level of information borrowing.
#' @param a_alpha A hyperparameter
#' @param b_alpha A hyperparameter
#' @param a_tau A hyperparameter
#' @param b_tau A hyperparameter
#' @param pi_delta A hyperparameter controlling delta, only useful when not doing enrichment.
#' If no variable selection is wanted, set \code{pi_delta = 1}.
#' @param sigma_gammaG_2 Hyperparameter of gammaG, only useful when
#' \code{enrichment = TRUE}
#' @param aggregation Whether to perform aggregation.
#' This will usually speed up the function without changing the result.
#' @param enrichment Whether to perform enrichment.
#'
#' @return \code{new_BV_chain} returns an object of S3 class \code{'BV_chain'} that is ready for method \code{update}
#'
#' \code{help('update.BV_chain')}
#'
#' \code{help('logLik.BV_chain')}
#'
#' @importFrom truncnorm rtruncnorm etruncnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom pgdraw pgdraw
#' @importClassesFrom Matrix dgCMatrix
#'
#' @export new_BV_chain

new_BV_chain = function(A,
                        V,
                        X,
                        L,
                        G = NULL,
                        nn = rep(1, nrow(A)),
                        epsilon = 0.2,
                        a_alpha = 0.1,
                        b_alpha = 0.1,
                        a_tau = 0.5,
                        b_tau = 0.5,
                        a_gammaG = 0.5,
                        b_gammaG = 0.5,
                        pi_delta = 0.5,
                        n_thin = 1,
                        aggregation = TRUE,
                        enrichment = FALSE) {
  L = scaleL(L, epsilon)

  data = list(
    A = A,
    V = V,
    X = X,
    L = L,
    nn = nn
  )
  if (aggregation) {
    data = aggregate_data(data)
  }


  J = ncol(A)
  if (enrichment) {
    gammaG = rep(0, ncol(G))
    tau_gammaG_2 = a_gammaG / (1+b_gammaG)
    sigma_gammaG = abs(rnorm(ncol(G),sd = sqrt(tau_gammaG_2)))

    omega = pgdraw(1, G %*% (gammaG * sigma_gammaG))
    pi_delta = 1 / (1 + exp(-G %*% gammaG))
    delta = rbinom(J, 1, pi_delta)
  } else {
    delta = rbinom(J, 1, pi_delta)
  }
  tau_2 = b_tau / (a_tau + 1)
  sigma_beta = rep(1, J)
  sigma_alpha_2 = rep(b_alpha / (a_alpha + 1), ncol(X))
  glm_coefs = get_glm_coefs(data)
  beta = glm_coefs$beta_glm
  Alpha = glm_coefs$Alpha_glm
  Omega = draw_Omega(data, Alpha, beta * delta * sigma_beta, aggregation)

  hyperparameters = list(
    a_alpha = a_alpha,
    b_alpha = b_alpha,
    a_tau = a_tau,
    b_tau = b_tau
  )
  if (enrichment) {
    data$G = G
    hyperparameters$a_gammaG = a_gammaG
    hyperparameters$b_gammaG = b_gammaG
  } else {
    hyperparameters$pi_delta = pi_delta
  }

  coef_lst = list(
    delta = delta,
    tau_2 = tau_2,
    sigma_beta = sigma_beta,
    sigma_alpha_2 = sigma_alpha_2,
    beta = beta,
    Alpha = Alpha
  )
  if (enrichment) {
    coef_lst$gammaG = gammaG
    coef_lst$sigma_gammaG= sigma_gammaG
    coef_lst$tau_gammaG_2 = tau_gammaG_2
    coef_lst$omega = omega
  }

  coeflst = one_ite_summ(coef_lst)
  BV_chain = list(
    data = data,
    hyperparameters = hyperparameters,
    chain = list(coef_lst),
    start_pos = 1,
    end_pos = 1,
    aggregation = aggregation,
    enrichment = enrichment
  )

  class(BV_chain) = 'BV_chain'
  BV_chain$logl = logLik(BV_chain)
  BV_chain$n_thin = n_thin

  return(BV_chain)
}

# ------------
aggregate_data = function(data) {
  list2env(data, envir = environment())
  keys = sapply(1:length(V), function(i) {
    key = c(X[i,], V[i])
    return(paste0(key, collapse = '_'))
  })
  indi_sets = split(1:nrow(X), keys)
  names(indi_sets) = NULL
  new_n = length(indi_sets)
  new_nn = sapply(indi_sets, function(indi_set) {
    return(sum(nn[indi_set]))
  })
  new_A = matrix(0, new_n, ncol(A))
  new_X = matrix(0, new_n, ncol(X))
  new_V = rep(0, new_n)
  for (i in 1:new_n) {
    indi_set = indi_sets[[i]]
    new_A[i,] = colSums(A[indi_set, , drop = FALSE])
    new_X[i,] = X[indi_set[1],]
    new_V[i] = V[indi_set[1]]
  }
  data$A = new_A
  data$X = new_X
  data$V = new_V
  data$nn = new_nn
  return(data)
}

# --------------------

get_glm_coefs = function(data,
                         cutoff = 10) {
  list2env(data, envir = environment())
  # run glm models
  Alpha_glm = matrix(0, ncol(A), ncol(X) + 1)
  suppressWarnings({
    for (j in 1:ncol(A)) {
      Aj = A[, j]
      model = glm(Aj / nn ~ 0 + X + V,
                  family = binomial(),
                  weights = nn)
      alpha_glm = coef(model)
      alpha_glm = pmin(alpha_glm,cutoff)
      alpha_glm = pmax(alpha_glm,-cutoff)
      Alpha_glm[j,] =  alpha_glm
    }
  })

  rownames(Alpha_glm) = NULL
  beta_glm = Alpha_glm[, ncol(Alpha_glm)]
  Alpha_glm = Alpha_glm[,-ncol(Alpha_glm)]
  colnames(Alpha_glm) = colnames(X)

  return(list(beta_glm = beta_glm,
              Alpha_glm = Alpha_glm))
}

# --------------------
scaleL = function(L, epsilon) {
  J = ncol(L)
  if (epsilon == Inf) {
    corMat_inv = diag(J)
  } else {
    PreMat = L + epsilon * diag(J)
    W = sqrt(diag(solve(PreMat)))
    corMat_inv = (W %*% t(W)) * PreMat
  }
  return(corMat_inv)
}
