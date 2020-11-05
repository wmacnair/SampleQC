context("Sample-level functions")
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

# specify dummy annotation
annot_disc  = c('annot_1')

################
# tests
################

test_that("does sample function work?", {
    # do they work ok?
    expect_is(calculate_sample_to_sample_MMDs(qc_dt, qc_names, subsample=20, n_times=5, n_cores=1), 'SingleCellExperiment')
})

test_that("automatic handling of annot_disc, annot_continuous", {
    # run default
    qc_obj    = calculate_sample_to_sample_MMDs(qc_dt, qc_names, 
        subsample=20, n_times=5, n_cores=1)

    # get right type of output
    expect_equal(metadata(qc_obj)$annots$disc, c('group_id', 'N_cat', 'mito_cat', 'counts_cat'))
    expect_equal(metadata(qc_obj)$annots$cont, c('log_N', 'med_mito', 'med_counts'))

    # run default
    qc_obj    = calculate_sample_to_sample_MMDs(qc_dt, qc_names, 
        annots_disc=annot_disc, subsample=20, n_times=5, n_cores=1)

    # get right type of output
    expect_equal(metadata(qc_obj)$annots$disc, c('group_id', 'annot_1', 'N_cat', 'mito_cat', 'counts_cat'))
    expect_equal(metadata(qc_obj)$annots$cont, c('log_N', 'med_mito', 'med_counts'))
})

test_that("MMD seeds should be replicable", {
    # related to Simone's comment
})