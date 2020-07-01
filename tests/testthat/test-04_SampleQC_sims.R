context("Simulations")
# devtools::document(pkg_dir); devtools::test(pkg_dir)
# devtools::document(pkg_dir); testthat::test_file(file.path(pkg_dir, 'tests/testthat/test-04_SampleQC_sims.R'))

################
# set up
################

suppressPackageStartupMessages({
    library('data.table')
    library('SingleCellExperiment')
})
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)

################
# tests
################

test_that("check .sim_sample_group works", {
    # define parameters
    N           = 1e4
    J           = 10
    K           = 4
    D           = 3
    mu_0        = matrix(c(3, 3.5, -2), nrow=1)
    p_out_0     = 0.2
    p_loss_0    = 0.9

    # draw component parameters
    beta_k      = .draw_beta_k(D, K)
    Sigma_k     = .draw_Sigma_k(D, K)

    # subdivide cells into samples
    samples     = .draw_samples(J, N)

    # draw mean, covariance parameters
    alpha_j     = .draw_alpha_j(D, J)
    delta_jk    = .draw_delta_jk(D, J, K)

    # draw cluster membership parameters
    dir_0       = .draw_dir_0(K)
    p_jk        = .draw_p_jk(J, K, dir_0)
    z           = .draw_z(samples, p_jk, N, J, K)

    # draw outlier parameters
    out_j       = .draw_out_j(J, p_out_0, p_loss_0)
    outliers    = .draw_outliers(samples, out_j, N)

    # simulate all healthy cells
    x_ok        = .sim_ok_cells(z, samples,
        mu_0, alpha_j, beta_k, Sigma_k, delta_jk, 
        N, D, K, J)

    # adjust values for outliers
    x_out       = .sim_outliers(x_ok, samples, outliers, out_j, J)

    # do samples have correct shapes, values?
    expect_equal(length(samples), N)
    expect_equal(min(samples), 1)
    expect_equal(max(samples), J)

    # do means, covariances have correct shapes, values?
    expect_equal(dim(mu_0), c(1,D))
    expect_equal(dim(alpha_j), c(J,D))
    expect_equal(dim(beta_k), c(K,D))
    expect_equal(dim(Sigma_k), c(D,D,K ))
    expect_equal(dim(delta_jk), c(J*K,D))

    # alpha_j, beta_k, delta_jk column means should be 0
    expect_true( all.equal(colMeans(alpha_j), rep(0, D)) )
    expect_true( all.equal(colMeans(beta_k), rep(0, D)) )
    expect_true( all.equal(colMeans(delta_jk), rep(0, D)) )

    # do cluster membership params have correct shapes?
    expect_equal(length(dir_0), K)
    expect_equal(dim(p_jk), c(J,K))
    expect_equal(length(z), N)
    expect_equal(min(z), 1)
    expect_equal(max(z), K)

    # p_jk rows should sum to 1
    expect_true( all.equal(rowSums(p_jk), rep(1, J)) )

    # do outlier params have correct shapes?
    expect_equal(dim(out_j), c(J,2))
    expect_equal(length(outliers), N)

    # do outlier params have correct values?
    expect_gt(min(as.vector(out_j)), 0)
    expect_lt(max(as.vector(out_j)), 1)
    expect_setequal(outliers, c(0,1))

    # do simulated QC matrices have correct shapes?
    expect_equal(dim(x_ok), c(N,D))
    expect_equal(dim(x_out), c(N,D))

    # do simulated QC matrices have correct values?
    x_diff      = x_ok - x_out
    changed     = as.integer(rowSums(x_diff^2) > 0)
    expect_gte(min(as.vector(x_diff[, 1:2])), 0)
    expect_lte(min(as.vector(x_diff[, 3])), 0)
    expect_true(all.equal(outliers, changed))

    # check whole thing
    sims_group  = .sim_sample_group(
        N, D, J, K,
        mu_0, beta_k, Sigma_k,
        p_out_0, p_loss_0
        )
    expect_type(sims_group, 'list')
})

test_that("check sim_experiment works", {
    # specify default parameters
    n_groups    = 4
    n_cells     = 1e5
    cells_p_s   = 2000
    D           = 3
    K           = 4

    # label groups
    sample_groups   = sprintf('QC%01d', n_groups)

    # decide sizes of groups
    group_prob  = rdirichlet(1, rep(1, n_groups))
    Ns          = rmultinom(1, n_cells, prob=group_prob) %>% as.vector

    # how many samples in each? (J)
    Js          = rpois(n_groups, Ns/cells_p_s) + 1

    # generate mu_0 for each
    mu_0s       = .draw_mu_0s(D, n_groups)

    # generate common components
    beta_k      = .draw_beta_k(D, K)
    Sigma_k     = .draw_Sigma_k(D, K)

    # select which for each group
    sel_ks      = .draw_sel_ks(K, n_groups)

    # generate p_out params for each
    p_out_0s    = .draw_p_out_0s(n_groups)
    p_loss_0s   = .draw_p_loss_0s(n_groups)

    # do mu_0s have correct shapes, values?
    expect_equal(dim(mu_0s), c(n_groups,D))

    # do means, covariances have correct shapes, values?
    expect_equal(dim(beta_k), c(K,D))
    expect_equal(dim(Sigma_k), c(D,D,K ))

    # beta_k column means should be 0
    expect_true( all.equal(colMeans(beta_k), rep(0, D)) )

    # check whole thing
    qc_names    = c('log_counts', 'log_feats', 'logit_mito')
    expt_params = .draw_expt_params(n_groups, n_cells, cells_p_s, 
        D, K, qc_names)
    expect_type(expt_params, 'list')
})

test_that("check sim_experiment works", {
    # default parameters
    n_groups    = 4
    n_cells     = 1e5
    cells_p_s   = 2000
    D           = 3
    K           = 4
    qc_names    = c('log_counts', 'log_feats', 'logit_mito')

    # draw experiment-level parameters
    expt_params = .draw_expt_params(n_groups, n_cells, cells_p_s, 
        D, K, qc_names)

    # do for each
    group_sims  = lapply(
        1:n_groups,
        function(ii) {
            # extract relevant values
            N_ii        = expt_params$Ns[[ii]]
            J_ii        = expt_params$Js[[ii]]
            sel_k       = expt_params$sel_ks[ii, ]
            K_ii        = sum(sel_k)
            mu_0_ii     = expt_params$mu_0s[ii, ]

            # extract relevant components
            beta_k_ii   = expt_params$beta_k[sel_k, , drop=FALSE]
            Sigma_k_ii  = expt_params$Sigma_k[, , sel_k, drop=FALSE]

            # extract relevant outlier parameters
            p_out_0_ii  = expt_params$p_out_0s[[ii]]
            p_loss_0_ii = expt_params$p_loss_0s[[ii]]

            # simulate
            sims_ii     = .sim_sample_group(
                N_ii, D, J_ii, K_ii,
                mu_0_ii, 
                beta_k_ii, Sigma_k_ii,
                p_out_0_ii, p_loss_0_ii
                )

            return(sims_ii)
        }) %>% setNames(expt_params$sample_groups)

    # join together
    sims_list   = .process_outputs(group_sims, expt_params)

    # check whole thing
    sims_list   = sim_experiment(n_groups, n_cells, cells_p_s, D, K)
    expect_type(expt_params, 'list')
})

test_that("seeds are replicable", {
    # maybe not necessary, since not using Rcpp
})

test_that("z sampling works", {
})
