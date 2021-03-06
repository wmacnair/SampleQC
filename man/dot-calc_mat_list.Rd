% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_pairwise_mmds.R
\name{.calc_mat_list}
\alias{.calc_mat_list}
\title{Splits qc_dt into list of matrices, one for each sample}
\usage{
.calc_mat_list(qc_dt, qc_names, sample_list)
}
\arguments{
\item{qc_dt}{Data.table of QC metrics for all cells and samples}

\item{qc_names}{List of metrics to actually use for calculating}

\item{sample_list}{List of samples to iterate through}
}
\value{
List of matrices
}
\description{
Splits qc_dt into list of matrices, one for each sample
}
\keyword{internal}
