% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/new_BGrass_chain.R
\name{new_BGrass_chain}
\alias{new_BGrass_chain}
\title{Initialize a new BGrass chain}
\usage{
new_BGrass_chain(
  A,
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
  enrichment = FALSE
)
}
\arguments{
\item{A}{A 0-1 matrix with n rows and J columns (J AEs).
The i-j th element of A indicates whether the j th AE is reported in report i.}

\item{V}{A 0-1 vector of length n indicating the type pf vaccine in each report.}

\item{X}{The covariate matrix of dim n by p + 1 (p covariates and a column of intercept)}

\item{L}{The J by J normalized graph Laplacian matrix of AEs.
You can get it using the function \code{getLaplacian}.}

\item{G}{The J by K matrix indicating the groups each AE comes from. if \code{!enrichment},
this argument should be set to \code{NULL}.
You can get it using the function \code{getG}. (It is not useful right now)}

\item{nn}{If \code{A} is already aggregated, indicate the number of
reports corresponding the each row of \code{A}. The function will
assume the data is not aggregated, so \code{nn = rep(1, nrow(A))} by default.}

\item{epsilon}{A hyperparameter, controlling the level of information borrowing.}

\item{a_alpha}{A hyperparameter}

\item{b_alpha}{A hyperparameter}

\item{a_tau}{A hyperparameter}

\item{b_tau}{A hyperparameter}

\item{a_gammaG}{For enrichment, not useful right now.}

\item{b_gammaG}{For enrichment, not useful right now.}

\item{pi_delta}{A hyperparameter controlling delta, only useful when not doing enrichment.
If no variable selection is wanted, set \code{pi_delta = 1}.}

\item{n_thin}{Thinning parameter. The number of draws for one valid MCMC sample.}

\item{aggregation}{Whether to perform aggregation.
This will usually speed up the function without changing the result.}

\item{enrichment}{Whether to perform enrichment. This is under development, do not use this.}
}
\value{
\code{new_BGrass_chain} returns an object of S3 class \code{'BGrass_chain'} that is ready for method \code{update}

\code{help('update.BGrass_chain')}

\code{help('logLik.BGrass_chain')}
}
\description{
Initialize a new BGrass chain
}
