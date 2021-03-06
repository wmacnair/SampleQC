% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_qcs.R
\name{.sim_outliers}
\alias{.sim_outliers}
\title{Simulates outliers from ok cells}
\usage{
.sim_outliers(x_ok, samples, outliers, out_j, J)
}
\arguments{
\item{x_ok}{healthy cells}

\item{samples}{list of sample membership}

\item{out_j}{parameters for sample outliers}

\item{J}{number of samples}
}
\value{
\code{matrix} of cell QC metrics including outliers
}
\description{
Outliers are assumed to be a combination of loss of non-mitochondrial
counts and features. Each sample has a different proportion of outlier
cells, and a different proportion of log counts lost.
}
\keyword{internal}
