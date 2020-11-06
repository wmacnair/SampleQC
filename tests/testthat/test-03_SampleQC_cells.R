context("Cell-level functions")
# pkg_dir     = '/home/will/work/packages/SampleQC'
# devtools::document(pkg_dir); devtools::test(pkg_dir)

################
# set up
################

suppressPackageStartupMessages({
    library('data.table')
    library('SingleCellExperiment')
})
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)

# generate toy dataset
qc_df       = .make_toy_qc_df()
qc_names    = c('log_counts', 'log_feats', 'logit_mito')
qc_dt       = make_qc_dt(qc_df, qc_names)

# run mmds
annot_disc  = c('annot_1')
qc_obj      = calculate_sample_to_sample_MMDs(qc_dt, qc_names, subsample=20, n_times=5, n_cores=1)

# define K_list
K_list      = rep(1, get_n_groups(qc_obj))

################
# tests
################

test_that("cell functions work", {
    # does it work ok with defaults?
    expect_is(fit_sampleQC(qc_obj, K_all=1), 'SingleCellExperiment')
})

test_that("parameter specifications for one vs multiple sample clusters are correct", {
    # K_all, K_list specified at the same time, or neither
    expect_error(fit_sampleQC(qc_obj, K_all=1, K_list=K_list))
    expect_error(fit_sampleQC(qc_obj))
    expect_error(fit_sampleQC(qc_obj, K_all=NULL, K_list=NULL))

    # do they both work individually?
    expect_is(fit_sampleQC(qc_obj, K_all=1), 'SingleCellExperiment')
    expect_is(fit_sampleQC(qc_obj, K_list=K_list), 'SingleCellExperiment')

    # non-integer specifications
    expect_error(fit_sampleQC(qc_obj, K_all=1.4))
    K_list_tmp      = K_list
    K_list_tmp[[1]] = 1.3
    expect_error(fit_sampleQC(qc_obj, K_list=K_list_tmp))

    # non-negative specifications
    expect_error(fit_sampleQC(qc_obj, K_all=0))
    expect_error(fit_sampleQC(qc_obj, K_all=-1))
    K_list_tmp      = K_list
    K_list_tmp[[1]] = -1
    expect_error(fit_sampleQC(qc_obj, K_list=c(1, -1, 2, 1)))

    # names of output lists
})

test_that("z and mahalanobis distances agree", {
})
