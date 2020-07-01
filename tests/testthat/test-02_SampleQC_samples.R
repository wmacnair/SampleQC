context("Sample-level functions")
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

annot_discrete  = c('annot_1')

################
# tests
################

test_that("do sample functions work?", {
    # do they work ok?
    expect_is(calculate_sample_to_sample_MMDs(qc_dt, qc_names, subsample=20, n_times=5, n_cores=1), 'list')
    mmd_list    = calculate_sample_to_sample_MMDs(qc_dt, qc_names, subsample=20, n_times=5, n_cores=1)
    expect_is(embed_sample_to_sample_MMDs(mmd_list, qc_dt, n_nhbrs=5), 'list')
})

test_that("automatic handling of annot_discrete, annot_continuous", {
    # run default
    mmd_list    = calculate_sample_to_sample_MMDs(qc_dt, qc_names, subsample=20, n_times=5, n_cores=1)
    mmd_list    = embed_sample_to_sample_MMDs(mmd_list, qc_dt, n_nhbrs=5)

    # get right type of output
    expect_equal(mmd_list$annot_discrete, c('QC_clust', 'N_cat', 'mito_cat', 'counts_cat'))
    expect_equal(mmd_list$annot_cont, c('log_N', 'med_mito', 'med_counts'))

    # run with annotations specified
    mmd_list    = calculate_sample_to_sample_MMDs(qc_dt, qc_names, subsample=20, n_times=5, n_cores=1)
    mmd_list    = embed_sample_to_sample_MMDs(mmd_list, qc_dt, annot_discrete=annot_discrete, n_nhbrs=5)

    # get right type of output
    expect_equal(mmd_list$annot_discrete, c('annot_1', 'QC_clust', 'N_cat', 'mito_cat', 'counts_cat'))
    expect_equal(mmd_list$annot_cont, c('log_N', 'med_mito', 'med_counts'))
})

test_that("MMD seeds should be replicable", {
    # related to Simone's comment
})