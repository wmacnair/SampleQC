context("Preparing data.table of QC metrics")
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

# generate toy datasets
sce         = .make_toy_sce()
qc_df       = .make_toy_qc_df()

# some variables we need
qc_names    = c('log_counts', 'log_feats', 'logit_mito')
bad_names   = c('dummy1', 'dummy2', 'dummy3')


################
# tests
################

test_that("check make_qc_dt works", {
    # get right type of output
    expect_is(make_qc_dt(sce), 'data.table')
    expect_is(make_qc_dt(qc_df), 'data.table')
})

test_that("check arguments to make_qc_dt", {
    # errors when wrong thing used as input
    expect_error(make_qc_dt('dummy_string'))
    expect_error(make_qc_dt(1:10))

    # errors when qc names are missing
    expect_error(make_qc_dt(qc_df, bad_names))
    expect_warning(make_qc_dt(qc_df, qc_names=c('log_feats')))
    expect_warning(make_qc_dt(sce, qc_names=c('log_feats')))

    # errors when logit_mito is requested but log_counts not present
    qc_df_tmp    = copy(qc_df)
    qc_df_tmp$log_counts    = NULL
    expect_error(make_qc_dt(qc_df_tmp, qc_names=c('logit_mito', 'log_feats')))
    
    # errors when some entries are missing
    qc_df_tmp    = copy(qc_df)
    qc_df_tmp$log_counts[1] = NA
    expect_error(make_qc_dt(qc_df_tmp, qc_names))

    # errors when some entries are wrong
    qc_df_tmp    = copy(qc_df)
    qc_df_tmp$log_counts = -qc_df_tmp$log_counts
    expect_error(make_qc_dt(qc_df_tmp, qc_names))
})

test_that("sce should work as expected", {
    # get right type of output
    sce_tmp  = copy(sce)
    rownames(sce_tmp)[1] = 'not-mt-gene'
    expect_warning(make_qc_dt(sce_tmp))
})
