% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_qcs.R
\name{.sim_ok_cells}
\alias{.sim_ok_cells}
\title{Simulates healthy cells}
\usage{
.sim_ok_cells(
  z,
  samples,
  mu_0,
  alpha_j,
  beta_k,
  Sigma_k,
  delta_jk,
  N,
  D,
  K,
  J,
  df
)
}
\arguments{
\item{z}{latent true cell components}

\item{D}{number of dimensions}

\item{K}{number of components}

\item{J}{number of samples}
}
\value{
?
}
\description{
Assumes the following distribution:
x | z = k ~ MVN( mu_0 + alpha_j + beta_k + delta_jk, Sigma_k)
z | J = j ~ p_jk
}
\keyword{internal}
