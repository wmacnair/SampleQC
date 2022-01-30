# make_qc_dt.R
# Make nice data.table of QC metrics

#' Checks specified QC metrics and makes data.table for input to 
#' \code{calc_pairwise_mmds}.
#' 
#' Takes a \code{data.frame} of raw QC metrics, and makes a nice neat 
#' \code{data.table} output that can be used in \pkg{SampleQC}. For example, 
#' users with a \pkg{SingleCellExperiment} object \code{sce} may first run 
#' \code{scater::calculateQCMetrics}, then call \code{make_qc_dt(colData(sce))}.
#' We work with \code{data.frame}/\code{data.table} objects to have the most 
#' flexible possible approach (and to save work on e.g. keeping up with changes
#' to dependencies like \pkg{SingleCellExperiment} and \pkg{Seurat}).
#' 
#' This code also calculates some sample-level statistics, e.g. median log 
#' library size per sample, and adds columns with binned values for these.
#' 
#' @param qc_df data.frame object containing calculated QC metrics
#' @param sample_var which column of qc_df has sample labels? (e.g. sample, group, 
#' batch, library)
#' @param qc_names list of qc_names that need to be extracted
#' @param annot_vars list of user-specified sample-level annotations
#' @importFrom assertthat assert_that
#' @importFrom data.table setcolorder
#' @return qc_dt, a data.table containing the sample variable plus qc metrics
#' @export
make_qc_dt <- function(qc_df, sample_var = 'sample_id', 
  qc_names = c('log_counts', 'log_feats', 'logit_mito'), annot_vars = NULL) {

  # some checks
  if ( 'DFrame' %in% class(qc_df) )
    qc_df      = as.data.frame(qc_df)
  assert_that( is.data.frame(qc_df), msg = "qc_df must be a data.frame" )

  assert_that( sample_var %in% colnames(qc_df),
    msg = sprintf("%s is listed as variable for samples but is not in data.frame", 
      sample_var))

  reserved_ns  = c('sample_id', 'group_id', 'cell_id')
  assert_that( length(intersect(annot_vars, reserved_ns)) == 0,
    msg = paste0("The following variable names are reserved and cannot be used ", 
      "as annot_vars:\n", paste(reserved_ns, collapse = ", ")))

  assert_that( all(annot_vars %in% names(qc_df)),
    msg = sprintf("the following variables are listed in annot_vars but not in qc_df:\n%s", 
      paste(setdiff(annot_vars, names(qc_df)), collapse = ", ")))

  # set up qc_dt
  qc_dt   = .init_qc_dt(qc_df, sample_var)

  # add known metrics
  if ('log_counts' %in% qc_names) {
    qc_dt   = .add_log_counts(qc_dt, qc_df)
  }
  if ('log_feats' %in% qc_names) {
    qc_dt   = .add_log_feats(qc_dt, qc_df)
  }
  if ('logit_mito' %in% qc_names) {
    qc_dt   = .add_logit_mito(qc_dt, qc_df)
  }
  if ('splice_ratio' %in% qc_names) {
    qc_dt   = .add_splice_ratio(qc_dt, qc_df)
  }

  # add unknown metrics
  qc_dt   = .add_unknown_metrics(qc_dt, qc_df, qc_names)

  # add some useful annotations
  qc_dt   = .add_qc_annots(qc_dt)

  # add specified annotation variables
  qc_dt   = .add_annot_vars(qc_dt, qc_df, annot_vars)
  
  # put in nice order
  setcolorder(qc_dt, c('cell_id', 'sample_id', qc_names))

  # double-check everything is ok
  .check_qc_dt(qc_dt, qc_names, annot_vars)

  return(qc_dt)
}

#' Initializes qc_dt object
#'
#' @param qc_df input data
#' @param sample_var which column of df has sample labels? (e.g. sample, group, 
#' batch, library)
#' @importFrom assertthat assert_that
#' @keywords internal
.init_qc_dt <- function(qc_df, sample_var) {
  # add cell identifiers
  if ('cell_id' %in% colnames(qc_df)) {
    qc_dt   = data.table(cell_id = qc_df$cell_id)
  } else if ( !is.null(rownames(qc_df)) ) {
    qc_dt   = data.table(cell_id = rownames(qc_df))
  } else {
    stop("input data.frame must have either rownames or 'cell_id' as a column")
  }
  assert_that( length(unique(qc_dt$cell_id)) == nrow(qc_dt),
    msg = "cell identifiers are not unique")

  # add sample identifiers
  qc_dt[, sample_id := qc_df[[sample_var]] ]

  # check no missing values or NAs
  assert_that( all(!is.na(qc_dt$cell_id)), msg = "missing values in cell_id")
  assert_that( all(!is.na(qc_dt$sample_id)), msg = "missing values in sample_id")

  return(qc_dt)
}

#' Add log counts to qc_dt
#'
#' @param qc_dt data.table of QC metrics
#' @param qc_df input data
#' @importFrom data.table ":="
#' @importFrom assertthat assert_that
#' @keywords internal
.add_log_counts <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)
  valid_ns  = c('log_counts', 'total', 'sum', 'nCount_RNA')

  # check which are present
  here_ns   = vapply(valid_ns, function(v) v %in% df_names, logical(1))
  assert_that( sum(here_ns) >= 1,
    msg = paste0(
      "no valid column present for log_counts\n", 
      paste0("valid columns are: ", paste(valid_ns, collapse = ", "))
      ))
  to_use    = valid_ns[here_ns][[1]]

  # add values
  if (to_use %in% 'log_counts') {
    qc_dt[, log_counts := qc_df[[ to_use ]] ]

  } else if (to_use %in% c('total', 'sum', 'nCount_RNA')) {
    assert_that( all(qc_df[[ to_use ]] > 0) )
    qc_dt[, log_counts := log10(qc_df[[ to_use ]]) ]

  } else {
    stop("log_counts requested but required variables not present")

  }

  # do some checks
  assert_that( "log_counts" %in% names(qc_dt) ) 
  assert_that( !any(is.na(qc_dt$log_counts)), 
    msg = "some log_counts values are NA")
  assert_that( !any(is.infinite(qc_dt$log_counts)), 
    msg = "some log_counts values are infinite")
  assert_that( all(qc_dt$log_counts >= 0), 
    msg = "some log_counts values are <= 0")

  return(qc_dt)
}

#' Add log feats to qc_dt
#'
#' @param qc_dt data.table of QC metrics
#' @param qc_df input data
#' @importFrom data.table ":="
#' @importFrom assertthat assert_that
#' @keywords internal
.add_log_feats <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)
  valid_ns  = c('log_feats', 'detected', 'nFeature_RNA')

  # check which are present
  here_ns   = vapply(valid_ns, function(v) v %in% df_names, logical(1))
  assert_that( sum(here_ns) >= 1,
    msg = paste0(
      "no valid column present for log_feats\n", 
      paste0("valid columns are: ", paste(valid_ns, collapse = ", "))
      ))
  to_use    = valid_ns[here_ns][[1]]

  # add values
  if (to_use %in% 'log_feats') {
    qc_dt[, log_feats := qc_df[[ to_use ]] ]

  } else if (to_use %in% c('detected', 'nFeature_RNA')) {
    assert_that( all(qc_df[[ to_use ]] > 0) )
    qc_dt[, log_feats := log10(qc_df[[ to_use ]]) ]

  } else {
    stop("log_feats requested but required variables not present")

  }

  # do some checks
  assert_that( "log_feats" %in% names(qc_dt) ) 
  assert_that( !any(is.na(qc_dt$log_feats)), 
    msg = "some log_feats values are NA")
  assert_that( !any(is.infinite(qc_dt$log_feats)), 
    msg = "some log_feats values are infinite")
  assert_that( all(qc_dt$log_feats >= 0), 
    msg = "some log_feats values are <= 0")

  return(qc_dt)
}

#' Add logit-transformed mitochondrial proportions to qc_dt
#'
#' @param qc_dt data.table of QC metrics
#' @param qc_df input data
#' @importFrom data.table ":="
#' @importFrom assertthat assert_that
#' @keywords internal
.add_logit_mito <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)

  # add logit-transformed mitochondrial proportion to qc_dt
  if ('logit_mito' %in% df_names) {
    qc_dt[, logit_mito  := qc_df$logit_mito ]

  } else if ( ('subsets_mito_sum' %in% df_names) & ('total' %in% df_names) ) {
    qc_dt[, logit_mito  := qlogis( (qc_df$subsets_mito_sum + 1) / (qc_df$total + 2) ) ]

  } else if ( ('subsets_mt_sum' %in% df_names) & ('total' %in% df_names) ) {
    qc_dt[, logit_mito  := qlogis( (qc_df$subsets_mt_sum + 1) / (qc_df$total + 2) ) ]

  } else if ( ('percent.mt' %in% df_names) & ('nCount_RNA' %in% df_names) ) {
    total_counts  = qc_df$nCount_RNA
    mt_counts     = qc_df$nCount_RNA * qc_df$percent.mt / 100
    assert_that( all(abs(mt_counts - round(mt_counts, 0)) < 1e-10) )
    qc_dt[, logit_mito  := qlogis( (mt_counts + 1) / (total_counts + 2) ) ]

  } else if ( ('mito_prop' %in% df_names) & ('log_counts' %in% df_names) ) {
    total_counts  = 10^qc_df$log_counts
    mt_counts     = qc_df$mito_prop * total_counts
    assert_that( all(abs(mt_counts - round(mt_counts, 0)) < 1e-8) )
    qc_dt[, logit_mito  := qlogis( (mt_counts + 1) / (total_counts + 2) ) ]

  } else {
    stop("logit_mito requested but required variables not present")
  }

  # do some checks
  assert_that( "logit_mito" %in% names(qc_dt) ) 
  assert_that( !any(is.na(qc_dt$logit_mito)), 
    msg = "some logit_mito values are NA")
  assert_that( !any(is.infinite(qc_dt$logit_mito)), 
    msg = "some logit_mito values are infinite")

  return(qc_dt)
}

#' Add log splice ratio to qc_dt
#'
#' @param qc_dt data.table of QC metrics
#' @param qc_df input data
#' @importFrom data.table ":="
#' @importFrom assertthat assert_that
#' @keywords internal
.add_splice_ratio <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)

  # add logit-transformed mitochondrial proportion to qc_dt
  if ('splice_ratio' %in% df_names) {
    qc_dt[, splice_ratio  := qc_df$splice_ratio ]

  } else if ( ('total_spliced' %in% df_names) & ('total_unspliced' %in% df_names) ) {
    qc_dt[, splice_ratio  := qlogis( (qc_df$total_spliced + 1) / (qc_df$total_unspliced + 1) ) ]

  } else {
    stop("logit_mito requested but required variables not present")

  }

  # do some checks
  assert_that( "splice_ratio" %in% names(qc_dt) ) 
  assert_that( !any(is.na(qc_dt$logit_mito)), 
    msg = "some logit_mito values are NA")
  assert_that( !any(is.infinite(qc_dt$logit_mito)), 
    msg = "some logit_mito values are infinite")

  return(qc_dt)
}

#' Shows the list of QC metrics that for which \pkg{SampleQC} currently has 
#' specific functionality. \pkg{SampleQC} will happily use metrics that aren't
#' in this list, however for those in this list it can plot a couple of extra 
#' things.
#' 
#' @return character vector of QC metrics that SampleQC knows about
#' @export
list_known_metrics <- function() {
  return(c('log_counts', 'log_feats', 'logit_mito', 'splice_ratio'))
}

#' Adds metrics that SampleQC doesn't have specific functions for
#'
#' @param qc_dt data.table of QC metrics
#' @param qc_df data.frame object containing calculated QC metrics
#' @param qc_names list of qc_names that need to be extracted
#' @importFrom assertthat assert_that
#' @importFrom data.table set
#' @keywords internal
.add_unknown_metrics <- function(qc_dt, qc_df, qc_names)  {
  # anything to add?
  to_add  = setdiff(qc_names, list_known_metrics())
  if ( length(to_add) == 0 )
    return(qc_dt)

  # add them
  message("adding the following metrics that are not known to `SampleQC`:")
  message(paste(to_add, collapse = ", "))
  for (v in to_add) {
    assert_that( v %in% names(qc_df), msg = paste0(v, " missing from qc_df"))
    set(qc_dt, i = NULL, v, qc_df[[v]])
    assert_that( !any(is.na(qc_dt$v)), msg = paste0("NA values for ", v))
    assert_that( !any(is.infinite(qc_dt$v)), msg = paste0("infinite values for ", v))
  }

  return(qc_dt)
}

#' Adds sample-level annotations for each known QC metric
#'
#' @param qc_dt data.table
#' @importFrom data.table ":="
#' @return qc_dt, a data.table containing the sample variable plus qc metrics
#' @keywords internal
.add_qc_annots <- function(qc_dt) {
  # add annotations for sample size
  qc_dt[, log_N  := log10(.N), by='sample_id']

  # and factor version
  N_cuts    = c(1,100,200,400,1000,2000,4000,10000,20000,40000,Inf)
  N_labs    = paste0('<=', N_cuts[-1])
  qc_dt[, N_cat  := factor(
    cut(10^log_N, breaks = N_cuts, labels = N_labs), 
    levels = N_labs), by = 'sample_id']

  # add annotations relating to library sizes
  if ('log_counts' %in% names(qc_dt) ) {
    # add median log counts per sample
    qc_dt[, med_counts  := median(log_counts), by='sample_id']

    # put mito level into categories
    counts_cuts = c(1,100,300,1000,3000,10000,30000, Inf)
    counts_labs = paste0('<=', counts_cuts[-1])
    qc_dt[, counts_cat  := factor(
      cut(10^med_counts, breaks = counts_cuts, labels = counts_labs), 
      levels = counts_labs), by = 'sample_id']
  }

  # add annotations relating to features
  if ('log_feats' %in% names(qc_dt) ) {
    # add median log feats per sample
    qc_dt[, med_feats   := median(log_feats), by='sample_id']

    # put mito level into categories
    feats_cuts = c(1,100,300,1000,3000,10000,30000, Inf)
    feats_labs = paste0('<=', feats_cuts[-1])
    qc_dt[, feats_cat  := factor(
      cut(10^med_feats, breaks = feats_cuts, labels = feats_labs), 
      levels = feats_labs), by = 'sample_id']
  }

  # add annotations relating to mitochondrial proportions
  if ('logit_mito' %in% names(qc_dt) ) {
    # add median mito proportion
    qc_dt[, med_mito  := median(plogis(logit_mito)), by='sample_id']

    # put mito level into categories
    mito_cuts   = c(0,0.01,0.05,0.1,0.2,0.5,1)
    mito_labs   = paste0('<=', mito_cuts[-1])
    qc_dt[, mito_cat  := factor(
      cut(med_mito, breaks = mito_cuts, labels = mito_labs), 
      levels = mito_labs), by = 'sample_id']
  }

  # add annotations relating to mitochondrial proportions
  if ('splice_ratio' %in% names(qc_dt) ) {
    # add median mito proportion
    qc_dt[, med_splice  := median(plogis(splice_ratio)), by='sample_id']

    # put mito level into categories
    splice_cuts = c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1)
    splice_labs = paste0('<=', splice_cuts[-1])
    qc_dt[, splice_cat  := factor(
      cut(med_splice, breaks = splice_cuts, labels = splice_labs), 
      levels = splice_labs), by = 'sample_id']
  }

  return(qc_dt)
}

#' Adds user-specified annotation variables
#'
#' @param qc_dt data.table
#' @param qc_df data.frame object containing calculated QC metrics
#' @importFrom magrittr "%>%"
#' @importFrom data.table ":=" ".N"
#' @return qc_dt, a data.table containing the sample variable plus qc metrics
#' @keywords internal
.add_annot_vars <- function(qc_dt, qc_df, annot_vars) {
  # check all present
  assert_that( all(annot_vars %in% names(qc_df)) )
  # add them
  for (v in annot_vars) 
    qc_dt[[v]] = qc_df[[v]]

  # check that they all sample level
  for (v in annot_vars) {
    check_dt  = qc_dt[, c('sample_id', v), with = FALSE] %>%
      .[, .N, by = c('sample_id', v) ]
    assert_that( nrow(check_dt) == length(unique(check_dt$sample_id)),
      msg = paste0("annotation variable ", v, " has more than one value per\n", 
        "sample (should be sample-level only)"))
  }

  return(qc_dt)
}


#' Checks that output is ok
#'
#' @param qc_dt data.table
#' @param qc_names list of qc_names that need to be extracted
#' @param annot_vars list of annotation variables
#' @importFrom assertthat assert_that
#' @keywords internal
.check_qc_dt <- function(qc_dt, qc_names, annot_vars) {
  # unpack
  col_names   = colnames(qc_dt)

  # check specific names
  if ('log_counts' %in% col_names)
    assert_that( all(qc_dt$log_counts >= 0) )
  if ('log_feats' %in% col_names)
    assert_that( all(qc_dt$log_feats >= 0) )
  if ('logit_mito' %in% col_names)
    assert_that( all(is.finite(qc_dt$logit_mito)) )
  if ('splice_ratio' %in% col_names)
    assert_that( all(is.finite(qc_dt$splice_ratio)) )

  # check qc metrics and annotations for NAs
  for (n in qc_names) {
    assert_that( all(!is.na(qc_dt[[n]])) )
  }
  annots_auto   = c(
    "med_counts", "counts_cat", 
    "med_feats", "feats_cat", 
    "med_mito", "mito_cat", 
    "med_splice", "splice_cat", 
    "log_N", "N_cat")
  for (n in c(annots_auto, annot_vars)) {
    if ( n %in% names(qc_dt) )
      assert_that( all(!is.na(qc_dt[[n]])),
        msg = paste0('NA present in an annotation variable, ', n) )
  }
}
