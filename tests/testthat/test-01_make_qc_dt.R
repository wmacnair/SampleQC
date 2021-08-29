context("Preparing data.table of QC metrics")
# pkg_dir   = '/home/will/work/packages/SampleQC'
# devtools::document(pkg_dir); devtools::test(pkg_dir)
# devtools::document(pkg_dir); testthat::test_file(file.path(pkg_dir, 'tests/testthat/test-01_make_qc_dt.R'))

################
# set up
################

suppressPackageStartupMessages({
  library('data.table')
  library('SingleCellExperiment')
  library('Seurat')
  library('assertthat')
})
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)

# generate toy data.frame
qc_df     = .make_toy_qc_df()

# generate toy sce and Seurat
sce       = .make_toy_sce()
seu       = CreateSeuratObject(counts(sce), 
  meta.data = as.data.frame(colData(sce)))

# add QC metrics to both
sce       = scater::addPerCellQC(sce, 
  subsets = list(mito = grep("mt-", rownames(sce))))
seu[["percent.mt"]] = PercentageFeatureSet(seu, pattern = "^mt-")

# check
assert_that( ncol(sce) == ncol(seu), nrow(sce) == nrow(seu) )

# some variables we need
qc_names  = c('log_counts', 'log_feats', 'logit_mito')
bad_names = c('dummy1', 'dummy2', 'dummy3')


################
# tests
################

test_that("check make_qc_dt works in basic form", {
  # get right type of output
  expect_is(
    make_qc_dt(qc_df, sample_var = 'sample_id', annot_vars = 'annot_1'),
    'data.table')
  expect_is(
    make_qc_dt(colData(sce), sample_var = 'sample_id', annot_vars = 'annot_1'),
    'data.table')
  expect_is(
    make_qc_dt(seu@meta.data, sample_var = 'sample_id', annot_vars = 'annot_1'),
    'data.table')
})

test_that("check inputs (data.frame)", {
  
  # error when bad sample_id variable requested
  expect_error(make_qc_dt(qc_df, sample_var = 'bad_var'))

  # error when bad sample_id variable requested
  expect_error(make_qc_dt(qc_df, annot_vars = 'sample_id'))
  expect_error(make_qc_dt(qc_df, annot_vars = 'group_id'))
  expect_error(make_qc_dt(qc_df, annot_vars = 'cell_id'))
  
  # error when bad annotations are requested
  expect_error(make_qc_dt(qc_df, sample_var = 'sample_id', 
    annot_vars = 'bad_var'))

  # error when variables requested that are not present
  expect_error(make_qc_dt(qc_df, sample_var = 'sample_id', qc_names = 'log_splice'))

  # error when some variables are missing
  tmp_df    = copy(qc_df)
  tmp_df$log_counts = NULL
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id', 
    qc_names = c('log_counts', 'log_feats')))
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id', 
    qc_names = c('logit_mito', 'log_feats')))

  # error when some entries are missing
  tmp_df    = copy(qc_df)
  tmp_df$log_counts[1] = NA
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id'))

  # error when some entries are wrong
  tmp_df    = copy(qc_df)
  tmp_df$log_counts = -tmp_df$log_counts
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id'))
})

test_that("check inputs (sce)", {
  # error when bad sample_id variable requested
  expect_error(make_qc_dt(colData(sce), sample_var = 'bad_var'))
  
  # error when bad annotations are requested
  expect_error(make_qc_dt(colData(sce), sample_var = 'sample_id', 
    annot_vars = 'bad_var'))

  # error when variables requested that are not present
  expect_error(make_qc_dt(colData(sce), sample_var = 'sample_id', 
    qc_names = 'log_splice'))

  # error when some variables are missing
  tmp_df    = copy(colData(sce))
  tmp_df$total  = NULL
  tmp_df$sum    = NULL
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id', 
    qc_names = c('log_counts', 'log_feats')))
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id', 
    qc_names = c('logit_mito', 'log_feats')))

  # error when some entries are missing
  tmp_df    = copy(colData(sce))
  tmp_df$total[1] = NA
  tmp_df$sum[1]   = NA
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id'))

  # error when some entries are wrong
  tmp_df    = copy(colData(sce))
  tmp_df$total  = -tmp_df$total
  tmp_df$sum    = -tmp_df$sum
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id'))
})

test_that("check inputs (Seurat)", {
  # error when bad sample_id variable requested
  expect_error(make_qc_dt(seu@meta.data, sample_var = 'bad_var'))
  
  # error when bad annotations are requested
  expect_error(make_qc_dt(seu@meta.data, sample_var = 'sample_id', 
    annot_vars = 'bad_var'))

  # error when variables requested that are not present
  expect_error(make_qc_dt(seu@meta.data, sample_var = 'sample_id', qc_names = 'log_splice'))

  # error when some variables are missing
  tmp_df            = copy(seu@meta.data)
  tmp_df$nCount_RNA = NULL
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id', 
    qc_names = c('log_counts', 'log_feats')))
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id', 
    qc_names = c('logit_mito', 'log_feats')))

  # error when some entries are missing
  tmp_df                = copy(seu@meta.data)
  tmp_df$nCount_RNA[1]  = NA
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id'))

  # error when some entries are wrong
  tmp_df            = copy(seu@meta.data)
  tmp_df$nCount_RNA = -tmp_df$nCount_RNA
  expect_error(make_qc_dt(tmp_df, sample_var = 'sample_id'))
})
