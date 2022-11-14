context("Pairwise MMDs")
# pkg_dir     = '/home/will/work/packages/SampleQC'
# devtools::document(pkg_dir); devtools::test(pkg_dir)
# testthat::test_file(file.path(pkg_dir, 'tests/testthat/test-02_calc_pairwise_mmds.R'))

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

# prep function
qc_names    = c('log_counts', 'log_feats', 'logit_mito')
qc_dt       = make_qc_dt(qc_df, sample_var = 'sample_id', 
    qc_names, annot_vars = 'annot_1')

# make corner case
qc_dt_2     = qc_dt[ sample_id %in% c("sample01", "sample02") ]

# specify dummy annotation
annots_disc = c('annot_1')

################
# tests
################

test_that("does sample function work?", {
    # do they work ok?
    expect_is(calc_pairwise_mmds(qc_dt, qc_names, 
        subsample = 20, n_times = 5, n_cores = 1), 
      'SingleCellExperiment')
})

test_that("does sample function work for just 2 samples?", {
    # do they work ok?
    expect_is(calc_pairwise_mmds(qc_dt_2, qc_names, 
        subsample = 20, n_times = 5, n_cores = 1), 
      'SingleCellExperiment')
})

test_that("automatic handling of annots_disc, annot_continuous", {
    # run default
    suppressMessages({
      qc_obj    = calc_pairwise_mmds(qc_dt, qc_names, 
        subsample=20, n_times=5, n_cores=1)
    })

    # get right type of output
    expect_equal(metadata(qc_obj)$annots$disc, 
      c('group_id', 'N_cat', 'counts_cat', 'feats_cat', 'mito_cat'))
    expect_equal(metadata(qc_obj)$annots$cont, 
      c('log_N', 'med_counts', 'med_feats', 'med_mito'))

    # run default
    suppressMessages({
      qc_obj    = calc_pairwise_mmds(qc_dt, qc_names, 
        annots_disc = annots_disc, subsample = 20, n_times = 5, n_cores = 1)
    })

    # get right type of output
    expect_equal(metadata(qc_obj)$annots$disc, 
      c('group_id', 'annot_1', 'N_cat', 'counts_cat', 'feats_cat', 'mito_cat'))
    expect_equal(metadata(qc_obj)$annots$cont, 
      c('log_N', 'med_counts', 'med_feats', 'med_mito'))
})

test_that("MMD seeds should be replicable", {
    # see Simone's comment
    sample_list     = sort(unique(qc_dt$sample_id))
    n_samples       = length(sample_list)
    centre_samples  = TRUE
    scale_samples   = FALSE
    n_times         = 10
    n_cores         = 2
    sigma           = 3
    subsample       = 100

    # split qc metric values into one matrix per sample
    mat_list    = .calc_mat_list(qc_dt, qc_names, sample_list)

    # do it once
    invisible({
        mmds_1      = .calc_mmd_mat(sample_list, mat_list,
            n_times, subsample, sigma, n_cores)
        mmds_2      = .calc_mmd_mat(sample_list, mat_list,
            n_times, subsample, sigma, n_cores)
    })
    expect_equal(mmds_1, mmds_2)
})
