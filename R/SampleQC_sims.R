# SampleQC: robust multivariate, multi-celltype, multi-sample quality control 
# for single cell RNA-seq
# SampleQC_sims.R
# Complex simulations of good and bad quality cells

# devtools::load_all('~/work/packages/SampleQC')
# devtools::document('~/work/packages/SampleQC')

#' Simulates QC metrics for whole experiment
#' 
#' @description 
#' Detailed simulation of QC metric distribution across whole experiment, 
#' assuming multiple celltypes, which are shared across different groups of 
#' samples ('sample groups'). Parameters are derived from real data and are 
#' hopefully reasonably realistic.
#' 
#' @details
#' The simulation first randomly generates \emph{K} celltypes (= mixture 
#' components), each with its own mean and covariance matrix. It also 
#' generates parameters for \emph{n_groups}, including: how many cells in each 
#' group; how many samples; which celltypes are present in that group; 
#' hyperparameters for outliers in that group.
#' 
#' Given these generated quantities, a set of QC metric vectors representing 
#' cells is drawn for each sample group, and combined into \emph{qc_dt}. Each 
#' row in this \code{data.table} represents a cell, and has annotations showing
#' 
#' @param n_groups How many sample groups?
#' @param n_cells How many cells in total to simulate
#' @param cells_p_s Cells per sample, average value (i.e. default is that each 
#' sample has 2000 cells on average)
#' @param D How many QC metrics do you want
#' @param K How many mixture components should there be in total? (see Details)
#' 
#' @return list with multiple entries:
#' - qc_ok: \code{data.table} of cell QC metrics \emph{before} outlier 
#' perturbation
#' - qc_out: \code{data.table} of cell QC metrics \emph{after} outlier 
#' perturbation
#' - x_ok, x_out: matrices of values in qc_ok, qc_out respectively
#' - groups: vector of true sample groups
#' - samples: vector of sample_ids
#' - z: vector of true celltype / mixture component values
#' - outliers: vector of outlier status (0/1)
#' - group_sims: detailed list of simulation outputs for each sample group
#' - expt_params: list of whole experiment-level parameters (e.g. celltype 
#' means and covariance matrices)
#' 
#' @export
sim_experiment <- function(n_groups=4, n_cells=1e5, cells_p_s=2000, D=3, K=4, 
    qc_names=c('log_counts', 'log_feats', 'logit_mito')) {
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
            theta_0_ii  = expt_params$theta_0s[[ii]]
            p_loss_0_ii = expt_params$p_loss_0s[[ii]]

            # simulate
            sims_ii     = .sim_sample_group(
                N_ii, D, J_ii, K_ii,
                mu_0_ii, 
                beta_k_ii, Sigma_k_ii,
                p_out_0_ii, theta_0_ii, p_loss_0_ii
                )

            return(sims_ii)
        }) %>% setNames(expt_params$sample_groups)

    # join together
    sims_list   = .process_outputs(group_sims, expt_params)

    return(sims_list)
}

#' Draws experiment-level parameters
#' 
#' @param n_groups Number of sample groups, each of which is composed of 
#' multiple samples
#' @param n_cells Total number of cells in experiment
#' @param cells_p_s Average number of cells per sample (not per 
#' sample group)
#' @param D Number of dimensions
#' @param K Number of components
#' @param qc_names Names for qc metrics
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#'
#' @return list of parameters
#' @keywords internal
.draw_expt_params <- function(n_groups, n_cells, cells_p_s, D, K, qc_names) {
    # label groups
    group_names = sprintf('QC%01d', seq_len(n_groups))

    # check qc_names
    assert_that( length(qc_names) == D, 
        msg='length of qc_names must be equal to D')

    # decide sizes of groups
    dirichlet_a = 10
    group_prob  = rdirichlet(1, rep(dirichlet_a, n_groups))
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
    theta_0s    = .draw_theta_0s(n_groups)
    p_loss_0s   = .draw_p_loss_0s(n_groups)

    # bundle together
    expt_params = list(
        n_groups        = n_groups,
        n_cells         = n_cells,
        cells_p_s       = cells_p_s,
        sample_groups   = group_names,
        D               = D,
        qc_names        = qc_names,
        K               = K,
        Ns              = Ns,
        Js              = Js,
        mu_0s           = mu_0s,
        beta_k          = beta_k,
        Sigma_k         = Sigma_k,
        sel_ks          = sel_ks,
        p_out_0s        = p_out_0s,
        theta_0s        = theta_0s,
        p_loss_0s       = p_loss_0s
        )

    return(expt_params)
}

#' Draws centre of sample group
#' 
#' @param D number of dimensions
#' @param n_groups number of sample groups
#' 
#' @importFrom mvtnorm rmvnorm
#' @importFrom magrittr "%>%"
#' @importFrom assertthat assert_that
#'
#' @return vector of D columns
#' @keywords internal
.draw_mu_0s <- function(D, n_groups) {
    # alpha_j ~ MVN
    mu_0_0          = matrix(
        c(log10(3000), log10(1500), qlogis(0.03)),
        nrow=1)

    # cor_12          = 0.98
    # cor_13          = -0.7
    # cor_23          = -0.8
    # L_mat           = matrix(0, nrow=D, ncol=D)
    # L_mat[1,1]      = 1
    # L_mat[2,1]      = cor_12
    # L_mat[3,1]      = cor_13
    # L_mat[2,2]      = sqrt(1 - cor_12^2)
    # L_mat[2,2]      = sqrt(1 - cor_12^2)
    # L_mat[3,2]      = (cor_23 - L_mat[3,1] * L_mat[2,1]) / L_mat[2,2]
    # L_mat[3,3]      = sqrt(1 - L_mat[3,1]^2 - L_mat[3,2]^2)
    # L_mat %*% t(L_mat)

    corr_mat        = diag(D)
    corr_mat[1,2]   = 0.98
    corr_mat[1,3]   = -0.7
    corr_mat[2,3]   = -0.8
    corr_mat        = corr_mat + t(corr_mat)
    diag(corr_mat)  = 1
    # check is positive def
    assert_that( all(eigen(corr_mat)$values>0) )

    sd_0            = c(0.3, 0.4, 1.5)
    Sigma           = diag(sd_0) %*% corr_mat %*% diag(sd_0)
    assert_that( all(eigen(Sigma)$values>0) )

    # simulate mu_0 values
    mu_0s     = rmvnorm(n_groups, mean=rep(0, D), sigma=Sigma) %>%
        sweep(2, mu_0_0, '+')

    # check
    assert_that(
        ncol(mu_0s) == D, 
        nrow(mu_0s) == n_groups, 
        msg = "mu_0s has wrong number of dimensions"
        )

    return(mu_0s)
}

#' Draws random parameters for mixture model components
#' 
#' @param D number of dimensions
#' @param K number of components
#' 
#' @importFrom assertthat assert_that
#' @importFrom mvtnorm rmvnorm
#' @importFrom data.table data.table
#' @importFrom magrittr "%>%"
#'
#' @return ?
#' @keywords internal
.draw_beta_k <- function(D, K) {
    # beta_k ~ MVN, based on real data
    beta_0          = matrix(c(0, 0, 0), nrow=1)
    corr_mat        = diag(D)
    corr_mat[1,2]   = 0.9
    corr_mat[1,3]   = -0.5
    corr_mat[2,3]   = -0.5
    corr_mat        = corr_mat + t(corr_mat)
    diag(corr_mat)  = 1

    # check is valid correlatin matrix
    assert_that( all(eigen(corr_mat)$values > 0) )

    # make into sigma matrxi
    beta_sd_0       = c(0.4, 0.3, 1.5)
    Sigma           = diag(beta_sd_0) %*% corr_mat %*% diag(beta_sd_0)

    # simulate beta_k, centre
    beta_k          = rmvnorm(K, mean=beta_0, sigma=Sigma) %>%
        sweep(2, colMeans(.), "-")

    # put ks in correct order
    k_order         = order(beta_k[, 1])
    beta_k          = beta_k[ k_order, ]

    # check outputs
    assert_that( ncol(beta_k) == D, 
        msg="beta_k has wrong number of dimensions")
    assert_that( nrow(beta_k) == K, 
        msg="beta_k has wrong number of rows")
    assert_that( all.equal(colMeans(beta_k), rep(0,D)), 
        msg="beta_k not centred")
    assert_that( all(diff(beta_k[,1])>0), 
        msg="beta_k not ordered")

    return(beta_k)
}

#' Draws random parameters for mixture model covariances
#' 
#' @param D number of dimensions
#' @param K number of components
#' 
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table
#' @importFrom magrittr "%>%"
#'
#' @return ?
#' @keywords internal
.draw_Sigma_k <- function(D, K) {
    # define scale function for Wishart distn
    # (based on real data)
    V_orig  = matrix(
        c(
            2.27e-02,  1.66e-02,  1.20e-03,
            1.66e-02,  1.66e-02, -1.26e-02,
            1.20e-03, -1.26e-02,  2.93e-01
        ), nrow=3)
    sigma_V = diag(V_orig) %>% sqrt

    corr_mat        = diag(D)
    corr_mat[1,2]   = 0.99
    corr_mat[1,3]   = -0.6
    corr_mat[2,3]   = -0.7
    corr_mat        = corr_mat + t(corr_mat)
    diag(corr_mat)  = 1
    # check is positive def
    assert_that( all(eigen(corr_mat)$values>0) )
    V       = diag(sigma_V) %*% corr_mat %*% diag(sigma_V)

    # draw covariance matrices
    df      = 10
    Sigma_k = rWishart(K, df, V/df)

    # Sigma_k = vapply(seq_len(K), function(k) {
    #     # do as Wishart
    #     # construct correlation matrix
    #     corr_mat        = diag(D)
    #     corr_mat[1,2]   = 0.99
    #     corr_mat[2,3]   = 0.3
    #     corr_mat        = corr_mat + t(corr_mat)
    #     diag(corr_mat)  = 1

    #     # turn into covariance matrix
    #     sd_vec          = rep(k/10, D)
    #     sigma_tmp       = diag(sd_vec) %*% corr_mat %*% diag(sd_vec)

    #     return(sigma_tmp)
    # }, array(0, c(D,D)) ) %>% array(dim=c(D, D, K))

    # check outputs
    assert_that( all(dim(Sigma_k) == c(D,D,K)), 
        msg="Sigma_k has wrong dimensions")

    return(Sigma_k)
}

#' Determines which components are observed in each sample group
#' 
#' @param K number of components
#' @param n_groups number of sample groups
#' 
#' @importFrom magrittr "%>%"
#' @importFrom assertthat assert_that
#'
#' @return matrix of 0s and 1s, where rows correspond to sample groups, and
#' columns correspond to components.
#' @keywords internal
.draw_sel_ks <- function(K, n_groups) {
    # keep drawing until every sample group has at least one component
    sel_ks  = matrix(0, ncol=K, nrow=n_groups)

    # generate things
    is_empty    = TRUE
    is_duped    = TRUE
    while( is_empty | is_duped ) {
        # generate random components
        sel_ks  = vapply(
            1:n_groups, 
            function(ii) rbinom(K, 1, 0.5) == 1,
            logical(K)) %>% t

        # check if any empty
        is_empty    = any(rowSums(sel_ks) == 0)
        # check if all different
        is_duped    = any( duplicated(sel_ks, MARGIN=1) )
    }

    # checks
    assert_that(
        nrow(sel_ks) == n_groups, 
        ncol(sel_ks) == K, 
        msg = "sel_ks has wrong number of dimensions"
        )
    assert_that(
        all(rowSums(sel_ks) >= 1), 
        msg = "sel_ks has some rows with no components allocated"
        )

    return(sel_ks)
}

#' Sample group-level probabilities of cells being outliers
#' 
#' @param n_groups number of sample groups
#' 
#' @importFrom assertthat assert_that
#'
#' @return vector of probabilities, one for each sample group
#' @keywords internal
.draw_p_out_0s <- function(n_groups) {
    # sample startpoints
    p_out_00    = 0.08
    p_out_sd    = 0.3
    p_out_0s    = plogis(qlogis(p_out_00) + p_out_sd * rnorm(n_groups))

    # checks
    assert_that(
        length(p_out_0s) == n_groups, 
        msg = "p_out_0s entries have wrong length"
        )
    assert_that(
        all(p_out_0s > 0), 
        all(p_out_0s < 1), 
        msg = "p_out_0s entries should be between 0 and 1"
        )

    return(p_out_0s)
}

#' Draws group-level probabilities of cells being outliers
#' 
#' @param n_groups number of sample groups
#' 
#' @importFrom assertthat assert_that
#'
#' @return vector of theta values for beta-binomial distributions, one for each sample group
#' @keywords internal
.draw_theta_0s <- function(n_groups) {
    # sample startpoints
    theta_00    = 4
    theta_sd    = 0.5
    theta_0s    = exp( theta_00 + theta_sd * rnorm(n_groups) )

    # checks
    assert_that(
        length(theta_0s) == n_groups, 
        msg = "theta_0s entries have wrong length"
        )

    return(theta_0s)
}

#' Draws group-level probabilities of reads being lost in outlier cells
#' 
#' @param n_groups number of sample groups
#' 
#' @importFrom assertthat assert_that
#'
#' @return vector of probabilities, one for each sample group
#' @keywords internal
.draw_p_loss_0s <- function(n_groups) {
    # sample startpoints
    p_loss_0s   = plogis(qlogis(0.5) + 0.5 * rnorm(n_groups))

    # checks
    assert_that(
        length(p_loss_0s) == n_groups, 
        msg = "p_loss_0s entries have wrong length"
        )
    assert_that(
        all(p_loss_0s > 0), 
        all(p_loss_0s < 1), 
        msg = "p_loss_0s entries should be between 0 and 1"
        )

    return(p_loss_0s)
}

#' Simulates QC metrics for one group of samples
#' 
#' @param D number of dimensions
#' @param J number of groups
#' @param K number of components
#' 
#' @section Details:
#' 
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table
#' @importFrom magrittr "%>%"
#'
#' @return list with all parameters for this group
#' @keywords internal
.sim_sample_group <- function(N, D, J, K, 
    mu_0, beta_k, Sigma_k, p_out_0, theta_0, p_loss_0) {

    # subdivide cells into samples
    samples     = .draw_samples(J, N)

    # draw mean, covariance parameters
    # mu_0        = .draw_mu_0(D)
    alpha_j     = .draw_alpha_j(D, J)
    delta_jk    = .draw_delta_jk(D, J, K)

    # draw cluster membership parameters
    dir_0       = .draw_dir_0(K)
    p_jk        = .draw_p_jk(J, K, dir_0)
    z           = .draw_z(samples, p_jk, N, J, K)

    # draw outlier parameters
    out_j       = .draw_out_j(J, p_out_0, theta_0, p_loss_0)
    outliers    = .draw_outliers(samples, out_j, N)

    # simulate all healthy cells
    x_ok        = .sim_ok_cells(z, samples,
        mu_0, alpha_j, beta_k, Sigma_k, delta_jk, 
        N, D, K, J)

    # adjust values for outliers
    x_out       = .sim_outliers(x_ok, samples, outliers, out_j, J)

    # make output list
    data_list = list(
        D           = D, 
        J           = J, 
        K           = K, 
        N           = N, 
        dir_0       = dir_0,
        mu_0        = mu_0,
        alpha_j     = alpha_j,
        beta_k      = beta_k,
        Sigma_k     = Sigma_k,
        delta_jk    = delta_jk,
        p_jk        = p_jk, 
        out_j       = out_j, 
        z           = z,
        samples     = samples, 
        outliers    = outliers, 
        x_ok        = x_ok, 
        x_out       = x_out
        )

    return(data_list)
}

#' Randomly splits up N cells into J samples with different sizes
#' 
#' @param J number of samples
#' @param N number of cells
#' 
#' @importFrom assertthat assert_that
#'
#' @return vector of samples
#' @keywords internal
.draw_samples <- function(J, N) {
    # generate samples
    j_weights   = exp(rnorm(J)/2)
    n_js        = as.vector(rmultinom(1, N, prob=j_weights))
    j_vals      = seq_len(J)
    samples     = rep(j_vals, times=n_js)

    # do some checks
    assert_that( max(samples) == J )
    assert_that( min(samples) == 1 )
    assert_that( length(samples) == N )
    assert_that( all.equal(sort(unique(samples)), 1:J) )

    return(samples)
}

#' Determines Dirichlet parameter for p_jk draws
#' 
#' @param K number of components
#' 
#' @importFrom assertthat assert_that
#'
#' @return vector of K columns
#' @keywords internal
.draw_dir_0 <- function(K) {
    alpha       = 5
    dir_0       = rep(alpha, K)

    assert_that( length(dir_0) == K )

    return(dir_0)
}

#' Draws random parameters for sample shifts
#' 
#' @param D number of dimensions
#' @param J number of samples
#' 
#' @importFrom assertthat assert_that
#' @importFrom mvtnorm rmvnorm
#' @importFrom data.table data.table
#' @importFrom magrittr "%>%"
#'
#' @return matrix of J rows by D columns
#' @keywords internal
.draw_alpha_j <- function(D, J) {
    # alpha_j ~ MVN, based on real data
    alpha_0         = matrix(c(0, 0, 0), nrow=1)
    corr_mat        = diag(D)
    corr_mat[1,2]   = 0.9
    corr_mat[1,3]   = -0.4
    corr_mat[2,3]   = -0.6
    corr_mat        = corr_mat + t(corr_mat)
    diag(corr_mat)  = 1
    alpha_j_sd      = c(0.25, 0.25, 0.6)
    Sigma           = diag(alpha_j_sd) %*% corr_mat %*% diag(alpha_j_sd)

    # simulate alpha_j, centre

    alpha_j         = rmvnorm(J, mean=alpha_0, sigma=Sigma) %>%
        sweep(2, colMeans(.), "-")

    # check outputs
    assert_that( ncol(alpha_j) == D, 
        msg="alpha_j has wrong number of dimensions")
    assert_that( nrow(alpha_j) == J, 
        msg="alpha_j has wrong number of rows")
    assert_that( all.equal(colMeans(alpha_j), rep(0,D)), 
        msg="alpha_j not centred")

    return(alpha_j)
}

#' Draws random parameters for mixture model components
#' 
#' @param D number of components
#' @param J number of groups
#' @param K number of components
#' 
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table
#' @importFrom magrittr "%>%"
#'
#' @return ?
#' @keywords internal
.draw_delta_jk <- function(D, J, K) {
    if (TRUE) {
        # delta_jk ~ MVN
        delta_0         = matrix(c(0, 0, 0), nrow=1)
        corr_mat        = diag(D)
        corr_mat[1,2]   = 0.95
        corr_mat[1,3]   = -0.6
        corr_mat[2,3]   = -0.7
        corr_mat        = corr_mat + t(corr_mat)
        diag(corr_mat)  = 1

        assert_that( all(eigen(corr_mat)$values > 0) )

        delta_jk_sd     = c(0.05, 0.05, 0.05)
        Sigma           = diag(delta_jk_sd) %*% corr_mat %*% diag(delta_jk_sd)

        # simulate delta_j, centre
        delta_jk        = rmvnorm(J*K, mean=delta_0, sigma=Sigma) %>%
            sweep(2, colMeans(.), "-")        
    } else {
        # just do zeros
        delta_jk        = matrix(0, nrow=J*K, ncol=D)
    }

    # check outputs
    assert_that( ncol(delta_jk) == D, 
        msg="delta_jk has wrong number of dimensions")
    assert_that( nrow(delta_jk) == J*K, 
        msg="delta_jk has wrong number of rows")
    assert_that( all.equal(colMeans(delta_jk), rep(0,D)), 
        msg="delta_jk not centred")

    return(delta_jk)
}

#' Draws random parameters for mixing parameters
#' 
#' @param J number of samples
#' @param K number of components
#' @param dir_0 parameters for Dirichlet distribution
#' 
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table
#' @importFrom gtools rdirichlet
#'
#' @return ?
#' @keywords internal
.draw_p_jk <- function(J, K, dir_0) {
    # for each group, sample p_jk
    p_jk        = matrix(NA, J, K)
    k_vals      = seq_len(K)
    for (j in seq_len(J)) {
        p_jk[j, ]   = rdirichlet(1, dir_0)
    }

    return(p_jk)
}

#' Draws latent true component for each cell
#' 
#' @param samples sample indices
#' @param p_jk component probabilities
#' @param N number of cells
#' @param J number of samples
#' @param K number of components
#' 
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table
#'
#' @return ?
#' @keywords internal
.draw_z <- function(samples, p_jk, N, J, K) {
    z           = matrix(NA, N, 1)
    k_vals      = seq_len(K)
    for (j in seq_len(J)) {
        j_idx       = samples == j
        z[j_idx, ]  = sample(k_vals, 
            size=sum(j_idx), replace=TRUE, 
            prob=p_jk[j, ])
    }
    return(z)
}

#' Draws random parameters to determine how to do outliers
#' 
#' @param J number of samples
#' 
#' @importFrom magrittr "%>%" set_colnames
#'
#' @return \code{matrix} with rows corresponding to samples
#' col 1 is \code{p_out}, proportion of each sample which is an outlier
#' col 2 is \code{p_lost}, mean proportion of non-mito counts lost in outliers
#' @keywords internal
.draw_out_j <- function(J, p_out_0, theta_0, p_loss_0) {
    # allow outlier fraction to vary by sample (p_out)
    p_out       = rbeta(J, shape1=theta_0 * p_out_0, theta_0 * (1-p_out_0))

    # allow outlier effect to vary by sample (p_lost)
    p_loss      = plogis( rnorm(J)*0.5 + qlogis(p_loss_0) )
    
    # put together
    out_j       = matrix(c(p_out, p_loss), ncol=2) %>%
        set_colnames(c('p_out', 'p_loss'))

    return(out_j)
}

#' Which cells are outliers?
#' 
#' @param samples sample indices
#' @param out_j probabilities of being outliers
#' @param N number of cells
#' 
#' @importFrom magrittr "%>%" set_colnames
#'
#' @return vector of outlier status
#' @keywords internal
.draw_outliers <- function(samples, out_j, N) {
    # extract outlier proportions
    p_out       = out_j[, 'p_out']
    outliers    = rbinom(N, rep(1,N), p_out[samples])

    return(outliers)
}

#' Simulates healthy cells
#' 
#' Assumes the following distribution:
#' x | z = k ~ MVN( mu_0 + alpha_j + beta_k + delta_jk, Sigma_k)
#' z | J = j ~ p_jk
#' 
#' @param z latent true cell components
#' @param D number of dimensions
#' @param K number of components
#' @param J number of samples
#' 
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table
#' @importFrom mvtnorm rmvnorm
#' @importFrom magrittr "%>%"
#'
#' @return ?
#' @keywords internal
.sim_ok_cells <- function(z, samples, 
    mu_0, alpha_j, beta_k, Sigma_k, delta_jk, 
    N, D, K, J) {

    # sample x
    x           = matrix(nrow=N, ncol=D)
    for (k in seq_len(K)) {
        # which entries, how many?
        k_idx           = z[, 1] == k
        n_k             = sum(k_idx)

        # what groups here?
        sample_idx      = samples[k_idx]

        # prep mu
        mu_0_mat        = matrix(rep(mu_0, n_k), ncol=D, byrow=TRUE)
        alpha_j_mat     = alpha_j[sample_idx, ]
        beta_k_mat      = matrix(rep(beta_k[k, ], n_k), ncol=D, byrow=TRUE)
        delta_idx       = (k-1)*J + sample_idx
        delta_jk_mat    = delta_jk[delta_idx, ]

        # prep sigma
        sigma       = Sigma_k[, ,k ]

        # do draw
        x[k_idx, ]  = mu_0_mat + 
            alpha_j_mat + beta_k_mat + delta_jk_mat +
            rmvnorm(sum(k_idx), mean=rep(0, D), sigma=sigma)
    }

    return(x)
}

#' Simulates outliers from ok cells
#' 
#' Outliers are assumed to be a combination of loss of non-mitochondrial 
#' counts and features. Each sample has a different proportion of outlier 
#' cells, and a different proportion of log counts lost.
#' 
#' @param x_ok healthy cells
#' @param samples list of sample membership
#' @param out_j parameters for sample outliers
#' @param J number of samples
#' 
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table
#' @importFrom magrittr "%>%"
#'
#' @return \code{matrix} of cell QC metrics including outliers
#' @keywords internal
.sim_outliers <- function(x_ok, samples, outliers, out_j, J) {
    # sample x
    x_out   = copy(x_ok)
    for (j in seq_len(J)) {
        # which group, cells?
        j_idx   = (samples == j) & (outliers == 1)
        x_tmp   = x_ok[ j_idx, , drop=FALSE]
        n_j     = nrow(x_tmp)

        if (n_j > 0) {
            # extract values
            total_old   = round(10^x_tmp[, 1],0)
            feats_old   = round(10^x_tmp[, 2],0)
            mt_prop_old = plogis(x_tmp[, 3])
            non_mt_old  = round(total_old * (1-mt_prop_old), 0)
            mt_old      = total_old - non_mt_old

            # downsample
            p_tmp       = out_j[j, 'p_loss']
            non_mt_new  = rbinom(n_j, non_mt_old, 1-p_tmp) + 1
            feats_new   = rbinom(n_j, feats_old, 1-p_tmp) + 1
            total_new   = non_mt_new + mt_old + 1
            mt_prop_new = (mt_old + 1) / total_new

            # update values
            x_new           = copy(x_tmp)
            x_new[, 1]      = log10(total_new)
            x_new[, 2]      = log10(feats_new)
            x_new[, 3]      = qlogis(mt_prop_new)
            x_out[j_idx, ]  = x_new
        }
    }

    # check no infinite values
    assert_that( sum(is.infinite(as.vector(x_out))) == 0,
        msg     = 'infinite values in x_out')

    return(x_out)
}

#' Gathers results simulated for individual sample groups together
#' 
#' Gathers results simulated for individual sample groups together into
#' results for a whole experiment. Results for individual sample groups are 
#' returned. 
#' 
#' @param group_sims List of results for individual sample groups
#' @param expt_params True parameters used to generate group_sims
#' 
#' @importFrom data.table data.table
#' @importFrom magrittr "%>%"
#' @importFrom magrittr set_colnames
#' @importFrom assertthat assert_that
#'
#' @return \code{list} of outputs
#' @keywords internal
.process_outputs <- function(group_sims, expt_params) {
    # define some useful things
    groups_idx      = seq_len(expt_params$n_groups)
    sample_groups   = expt_params$sample_groups
    comps           = seq_len(expt_params$K)

    # join x matrices together
    x_ok    = do.call(rbind, lapply(
        group_sims, function(s) s$x_ok)) %>% 
        set_colnames(expt_params$qc_names)
    x_out   = do.call(rbind, lapply(
        group_sims, function(s) s$x_out)) %>% 
        set_colnames(expt_params$qc_names)

    # define experiment-level sample group labels
    groups          = rep(sample_groups, times=expt_params$Ns)

    # define experiment-level sample labels, join
    Js_cumul        = cumsum(expt_params$Js)
    Js_cumul        = c(0, Js_cumul)
    labels_samples  = lapply(groups_idx, 
        function(jj) {
            jj_idx  = (Js_cumul[[jj]]+1):Js_cumul[[jj+1]]
            return(sprintf('sample%02d', jj_idx))
        }) %>% setNames(sample_groups)
    samples         = do.call(c, lapply(
        sample_groups, 
        function(s) labels_samples[[s]][group_sims[[s]]$samples]
        ))

    # define experiment-level components, join
    comp_labels     = lapply(groups_idx, 
        function(jj) {
            # which components used in this sample group?
            comps_used  = which( expt_params$sel_ks[jj, ] == 1 )
            comps_jj    = comps_used[ group_sims[[jj]]$z ]

            return(comps_jj)
        }) %>% unlist

    # join all outliers
    outliers        = do.call(c, lapply(
        sample_groups, function(s) group_sims[[s]]$outliers
        ))

    # put into a data.table which can be used as input to SampleQC
    id_pattern  = sprintf('cell%%0%dd', ceiling(log10(expt_params$n_cells)))
    cell_ids    = sprintf(id_pattern, seq.int(expt_params$n_cells))
    qc_ok       = data.table(cell_id=cell_ids, x_ok, sample_id=samples)
    qc_out      = data.table(cell_id=cell_ids, x_out, sample_id=samples)

    # put together
    sims_list = list(
        qc_ok       = qc_ok,
        qc_out      = qc_out,
        x_ok        = x_ok,
        x_out       = x_out,
        groups      = groups,
        samples     = samples,
        z           = comp_labels,
        outliers    = outliers,
        group_sims  = group_sims,
        expt_params = expt_params
        )

    return(sims_list)
}

#' Generate fake qc_df object for testing
#'
#' @importFrom data.table data.table rbindlist ":="
#' @return data.frame
#' @keywords internal
.make_toy_qc_df <- function() {
    # generate sample means
    J               = 10
    N_per_sample    = as.integer(exp(rnorm(J, log(200), 1)))
    mu_log_counts   = rnorm(J, log(4000), 0.1)
    mu_log_feats    = mu_log_counts - 0.5 + rnorm(J, 0, 0.1)
    mu_logit_mito   = rnorm(J, -3, 1)

    # generate cell data
    df_list = lapply(seq_len(J), function(j) {
        # unpack
        N           = N_per_sample[[j]]
        mu_counts   = mu_log_counts[[j]]
        mu_feats    = mu_log_feats[[j]]
        mu_mito     = mu_logit_mito[[j]]

        # generate random data
        df  = data.frame(
            sample_id   = sprintf('sample%02d', j),
            annot_1     = sprintf('annot%02d', j),
            log_counts  = rnorm(N, mu_counts, 1)
            )
        df$log_feats    = (df$log_counts - mu_counts) + mu_feats + rnorm(N, 0, 0.01)
        df$mito_prop    = plogis(rnorm(N, mu_mito, 0.1))

        return(df)
    })
    qc_df   = do.call(rbind, df_list)

    # add cell_id
    qc_df$cell_id = sprintf('cell%04d', seq_len(nrow(qc_df)))
    n_cols  = ncol(qc_df)
    qc_df   = qc_df[, c(n_cols, seq_len((n_cols-1)))]

    return(qc_df)
}

#' Generate fake sce object for testing
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @return sce
#' @keywords internal
.make_toy_sce <- function() {
    # generate sample_ids
    J               = 10
    N_per_sample    = as.integer(exp(rnorm(J, log(200), 1)))
    sample_ids      = rep(sprintf('sample%02d', seq_len(J)), times=N_per_sample)
    annot_1         = rep(sprintf('annot%02d', seq_len(J)), times=N_per_sample)

    # generate genes
    n_genes         = 100
    gene_names      = sprintf('gene%03d', seq_len(n_genes))
    n_mt            = 13
    mt_genes        = sprintf('mt-%02d', seq_len(n_mt))
    gene_names[seq_len(n_mt)]  = mt_genes

    # make count matrix
    n_cells         = length(sample_ids)
    counts_mat      = matrix(rpois(n_cells*n_genes, lambda=10), ncol=n_cells)

    # make column data
    cols_df         = data.frame(
        cell_id     = sprintf('cell%04d', seq_len(n_cells)),
        sample_id   = factor(sample_ids),
        annot_1     = factor(annot_1)
        )

    # make sce object
    sce         = SingleCellExperiment(
        list(counts=counts_mat),
        colData = cols_df
        )
    colnames(sce)   = cols_df$cell_id
    rownames(sce)   = gene_names

    return(sce)
}
