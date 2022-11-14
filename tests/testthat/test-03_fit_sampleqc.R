context("Fitting SampleQC model")
# pkg_dir   = '/home/will/work/packages/SampleQC'
# devtools::document(pkg_dir); devtools::test(pkg_dir)
# testthat::test_file(file.path(pkg_dir, 'tests/testthat/test-03_fit_sampleqc.R'))

################
# set up
################

suppressPackageStartupMessages({
  library('data.table')
  library('SingleCellExperiment')
})
seed      = as.numeric(format(Sys.time(), "%s"))
set.seed(seed)

# generate toy dataset
sims_list = simulate_qcs(n_cells = 5e3)
qc_names  = sims_list$expt_params$qc_names
qc_dt     = make_qc_dt(sims_list$qc_out, sample_var = 'sample_id',
  qc_names = qc_names)

# run mmds
annots_disc = c('annot_1')
suppressMessages({
  qc_obj    = calc_pairwise_mmds(qc_dt, qc_names, 
    subsample=20, n_times=5, n_cores=1)
})

# define K_list
K_list    = rep(1, get_n_groups(qc_obj))

# run for 2-sample edge case
qc_dt_2   = qc_dt[ sample_id %in% c("sample01", "sample02") ]
suppressMessages({
  qc_obj_2  = calc_pairwise_mmds(qc_dt_2, qc_names, 
    subsample=20, n_times=5, n_cores=1)
})


################
# tests
################

test_that("check fitting function works", {
  # does it work ok with defaults?
  expect_is(fit_sampleqc(qc_obj, K_all=1), 'SingleCellExperiment')
})

test_that("check fitting function work with 2 samples", {
  # does it work ok with defaults?
  expect_is(fit_sampleqc(qc_obj_2, K_all=1), 'SingleCellExperiment')
})

test_that("parameter specifications for one vs multiple sample clusters are correct", {
  # K_all, K_list specified at the same time, or neither
  expect_error(fit_sampleqc(qc_obj, K_all=1, K_list=K_list))
  expect_error(fit_sampleqc(qc_obj))
  expect_error(fit_sampleqc(qc_obj, K_all=NULL, K_list=NULL))

  # do they both work individually?
  expect_is(fit_sampleqc(qc_obj, K_all=1), 'SingleCellExperiment')
  expect_is(fit_sampleqc(qc_obj, K_list=K_list), 'SingleCellExperiment')

  # non-integer specifications
  expect_error(fit_sampleqc(qc_obj, K_all=1.4))
  K_list_tmp    = K_list
  K_list_tmp[[1]] = 1.3
  expect_error(fit_sampleqc(qc_obj, K_list=K_list_tmp))

  # non-negative specifications
  expect_error(fit_sampleqc(qc_obj, K_all=0))
  expect_error(fit_sampleqc(qc_obj, K_all=-1))
  K_list_tmp    = K_list
  K_list_tmp[[1]] = -1
  expect_error(fit_sampleqc(qc_obj, K_list=c(1, -1, 2, 1)))

  # names of output lists
})

test_that("seeds are replicable", {
  K_list_2      = rep(2, length(K_list))
  set.seed(4321)
  suppressMessages({fit_1 = fit_sampleqc(qc_obj, K_list = K_list_2, n_cores = 1)})
  set.seed(4321)
  suppressMessages({fit_2 = fit_sampleqc(qc_obj, K_list = K_list_2, n_cores = 1)})
  expect_equal(fit_1, fit_2)
})

test_that("z and mahalanobis distances agree", {
})
