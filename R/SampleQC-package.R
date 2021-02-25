#' SampleQC: robust multivariate, multi-celltype, multi-sample quality control for single cell RNA-seq
#'
#' `SampleQC` is an R package for robust multivariate, multi-
#' celltype, multi-sample quality control for single cell RNA-seq. Unimodal
#' approaches to QC, such as `scater`, may preferentially remove healthy cells 
#' from celltypes with extreme QC statistics, introducing bias. `SampleQC` 
#' removes this bias by fitting a Gaussian mixture model. The fit is done 
#' across multiple celltypes simultaneously, and using robust estimators
#' allows it to find the distributions in the (expected) presence of outliers.
#' 
#' @docType package
#' @name SampleQC
#' @useDynLib SampleQC, .registration = TRUE
NULL
#> NULL
