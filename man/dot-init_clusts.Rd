% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_sampleqc.R
\name{.init_clusts}
\alias{.init_clusts}
\title{Initializes clusters}
\usage{
.init_clusts(x, groups, J, K, N)
}
\arguments{
\item{x}{Matrix of QC values}

\item{groups}{sample_id labels as integers}

\item{J}{Number of samples}

\item{K}{Number of requested clusters}

\item{N}{Number of observations}
}
\value{
list, containing lots of cell outlier information
}
\description{
Initializes clusters
}
\keyword{internal}
