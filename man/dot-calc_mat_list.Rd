% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SampleQC_samples.R
\name{.calc_mat_list}
\alias{.calc_mat_list}
\title{Splits qc_dt into list of matrices, one for each sample}
\usage{
.calc_mat_list(qc_dt, qc_names, sample_list, centre_samples, scale_samples)
}
\arguments{
\item{qc_dt}{Data.table of QC metrics for all cells and samples}

\item{qc_names}{List of metrics to actually use for calculating}

\item{sample_list}{List of samples to iterate through}

\item{centre_samples, scale_samples}{Should we centre or scale the values
within each sample before calculating MMD?}
}
\value{
List of matrices
}
\description{
Splits qc_dt into list of matrices, one for each sample
}
\keyword{internal}