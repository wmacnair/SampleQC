context("Utils")
# pkg_dir   = '/home/will/work/packages/SampleQC'
# Rcpp::compileAttributes(pkg_dir); devtools::document(pkg_dir); devtools::test(pkg_dir)
# testthat::test_file(file.path(pkg_dir, 'tests/testthat/test-05_utils.R'))

################
# set up
################

suppressPackageStartupMessages({
  library('data.table')
  library('SingleCellExperiment')
})
seed        = as.numeric(format(Sys.time(), "%s"))
set.seed(seed)

# generate toy dataset
qc_df       = .make_toy_qc_df()
qc_names    = c('log_counts', 'log_feats', 'logit_mito')
qc_dt       = make_qc_dt(qc_df, sample_var = 'sample_id', qc_names = qc_names, 
  annot_vars = 'annot_1')

# run mmds
annot_disc  = c('annot_1')
qc_obj      = calc_pairwise_mmds(qc_dt, qc_names, subsample=20, n_times=5, n_cores=1)
mmd_only    = calc_pairwise_mmds(qc_dt, qc_names, subsample=20, n_times=5, n_cores=1)

# define K_list
K_list      = rep(1, metadata(qc_obj)$n_groups)
qc_obj      = fit_sampleqc(qc_obj, K_list=K_list)

################
# tests
################

test_that("get_n_groups works", {
  # should work with defaults
  expect_type(get_n_groups(qc_obj), 'integer')
  expect_true( metadata(qc_obj)$n_groups == get_n_groups(qc_obj) )

  # should fail on irrelevant things
  expect_error(get_n_groups(qc_dt))
})

test_that("get_outliers works", {
  # if fit_sampleqc not yet run, should throw error
  expect_is(get_outliers(qc_obj), 'data.table')
  expect_named(get_outliers(qc_obj), c('sample_id', 'cell_id', 'out_cell', 'out_sample', 'out_cluster', 'outlier'))

  # check all sample_ids present
  outliers_dt   = get_outliers(qc_obj)
  expect_setequal(outliers_dt$cell_id, qc_dt$cell_id)
  expect_setequal(outliers_dt$sample_id, qc_dt$sample_id)
  expect_setequal(outliers_dt$sample_id, qc_obj$sample_id)

  # if fit_sampleqc not yet run, should throw error
  expect_error(get_outliers(mmd_only))
})

