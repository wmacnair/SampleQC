# SampleQC: robust multivariate, multi-celltype, multi-sample quality control 
# for single cell RNA-seq
# SampleQC_cells.R
# Code to identify outlier cells within sample groups identified by code
# in SampleQC_samples. Includes plots of outputs, and plots of fitted 
# statistics

# devtools::load_all('~/work/packages/SampleQC')
# devtools::document('~/work/packages/SampleQC')

#' Fits `SampleQC` model to one cluster of samples
#' 
#' @description 
#' Running \link{calculate_sample_to_sample_MMDs} identifies groups 
#' of samples with similar distributions of QC metrics ('sample groups'). 
#' \link{SampleQC} then fits a multivariate Gaussian mixture model to 
#' each sample group, and uses these distributions to determine outlier cells
#' without celltype bias.
#'
#' @details
#' Determining the number of components to use in a Gaussian mixture model is 
#' a hard problem that \link{SampleQC} does not attempt to solve. 
#' However, \link{SampleQC} provides diagnostic plots that assist the 
#' user to make this decision themself. 
#' 
#' We recommend the following workflow. The first step is always to identify 
#' sample groups (and derive embeddings) by calling 
#' \link{calculate_sample_to_sample_MMDs}. Next, we recommend calling 
#' \link{fit_sampleQC} with \code{K_list=rep(1, get_n_groups(qc_obj)))},
#' then rendering outputs with \code{make_SampleQC_report(qc_obj, save_dir, 
#' 'test')}. This fits the simplest possible model to each sample group (i.e. 
#' one component), and renders a report with biaxial plots. 
#' 
#' The user can then use the biaxial plots to determine how many components 
#' are appropriate for each sample group. So if 
#' 
#' @param qc_obj Output from calculate_sample_to_sample_MMDs
#' @param K_all,K_list How many QC celltypes do we expect? Exactly one of K_all
#' and K_list should be specified. If the user wants to fit the same model to 
#' all samples they should specify `K_all` as an integer. If the user wants to 
#' fit a different model to each of the sample groups identified by 
#' `calculate_sample_to_sample_MMDs`, they should specify `K_list` as an 
#' integer vector with the same number of entries as sample groups. See Details.
#' @param n_cores How many cores to use; n_cores=1 runs in serial. The default 
#' depends on whether K_all or K_list is used: if K_all is set, n_cores=1; if 
#' K_list is set, n_cores=length(K_list).
#' @param alpha Chi squared threshold to define outliers
#' @param em_iters Maximum number of EM iterations
#' @param mcd_alpha,mcd_iters Parameters for robust estimation of celltype 
#' means and covariances. \emph{mcd_alpha} is the proportion of the data which 
#' the robust covariance estimator attempts to fit; adjust this to be lower 
#' when more outliers are present. \emph{mcd_iters} is the maximum number of 
#' iterations allowed for estimating the covariance matrix, and should only 
#' rarely be hit.
#' @param method Should the fitting be robust, or maximum likelihood 
#' estimation? In general, robust should be used.
#' @param bp_seed random seed for BiocParallel workers
#' @param track Track values of parameters during fitting (used for debugging)
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom mclust mclustBIC
#' @importFrom data.table setkey copy ":="
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam bpstart bpstop
#'
#' @return list, containing MMD values and sample clusters based on MMD values
#' @export
fit_sampleQC <- function(qc_obj, K_all=NULL, K_list=NULL, n_cores=NULL, 
    alpha=0.01, em_iters=50, mcd_alpha=0.5, mcd_iters=50, 
    method=c('robust', 'mle'), bp_seed=22, track=FALSE) {
    # check some inputs
    .check_is_qc_obj(qc_obj)
    method      = match.arg(method)
    fit_params  = .check_fit_params(qc_obj, K_all, K_list, n_cores, 
        alpha, em_iters, mcd_alpha, mcd_iters, method, track)

    # decide whether to fit to each sample group individually
    if (fit_params$do_list) {
        # split by sample groups
        df_list     = split.data.frame(colData(qc_obj), 
            colData(qc_obj)$group_id)

        # put in correct order
        df_list     = df_list[ metadata(qc_obj)$group_list ]

        # define parallelization
        if (fit_params$n_cores == 1) {
            bpparam = SerialParam()
        } else {
            bpparam = MulticoreParam(workers = fit_params$n_cores, 
                RNGseed = bp_seed)
        }

        # fit each sample group
        bpstart(bpparam)
        fit_list    = bplapply(
            seq_along(df_list), 
            function(i)
                .fit_one_sampleQC(df_list[[i]], K_list[[i]], fit_params, seed = i), 
                BPPARAM = bpparam)
        names(fit_list) = metadata(qc_obj)$group_list
        bpstart(bpstop)
    } else {
        # fit one model to all samples
        df          = colData(qc_obj)
        fit_list    = list(.fit_one_sampleQC(df, K_all, fit_params, seed = 1))
        names(fit_list) = 'All'
    }

    # reassemble
    qc_obj      = .assemble_fit_results(qc_obj, fit_list, fit_params)

    return(qc_obj)
}

#' Checks that the parameters specified are ok
#' 
#' @importFrom assertthat assert_that is.flag is.count
#' @importFrom S4Vectors metadata
#' 
#' @keywords internal
.check_fit_params <- function(qc_obj, K_all, K_list, n_cores, 
    alpha, em_iters, mcd_alpha, mcd_iters, method, track) {
    # check specification of type of run is ok
    if (!is.null(K_all) & !is.null(K_list))
        stop('only one of `K_all` and `K_list` should be specified at one 
            time; see ?fit_sampleQC')
    do_list   = !is.null(K_list)

    # check other parameters are ok
    assert_that(
        em_iters == as.integer(em_iters), em_iters > 0,
        msg     = 'em_iters must be an integer greater than 0')
    assert_that(
        mcd_iters == as.integer(mcd_iters), mcd_iters > 0,
        msg     = 'mcd_iters must be an integer greater than 0')
    assert_that(
        alpha > 0, alpha < 1, 
        msg     = 'alpha must be a values between 0 and 1')
    assert_that(
        mcd_alpha > 0, mcd_alpha <= 1, 
        msg     = 'mcd_alpha must be a values between 0 and 1')

    # check specification of number of components is ok for both options
    if (do_list) {
        n_groups    = metadata(qc_obj)$n_groups
        assert_that(
            all.equal(K_list, as.integer(K_list)), 
            all(K_list > 0),
            msg     = 'K_list must be a vector of integers greater than 0')
        assert_that(
            length(K_list) == n_groups, 
            msg     = paste0('K_list must have the same length as the number of\n',
                        'clusters in the mmd_clusts element of qc_obj'))
    } else {
        assert_that(
            K_all == as.integer(K_all), 
            K_all > 0, 
            msg     = 'K_all must be an integer greater than 0')
    }
    # check inputs
    if (do_list) {
        if (is.null(n_cores)) {
            n_cores     = length(K_list)
        } else {
            assert_that( is.count(n_cores) )
        }
    } else {
        n_cores     = 1
    }

    assert_that( is.flag(track),
        msg     = 'track must be either TRUE or FALSE')

    # assemble
    fit_params  = list(
        K_list      = K_list,
        K_all       = K_all,
        n_cores     = n_cores,
        method      = method,
        alpha       = alpha,
        em_iters    = em_iters,
        mcd_iters   = mcd_iters,
        mcd_alpha   = mcd_alpha,
        track       = track,
        do_list     = do_list
        )
    
    return(fit_params)
}

#' Fits `SampleQC` model to one cluster of samples
#' 
#' @param df DataFrame with lots of nested things
#' @param K How many QC celltypes are there?
#' @param fit_params List of parameters specifying the model
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom mclust mclustBIC summaryMclustBIC
#' @importFrom data.table setkey
#'
#' @return list, containing lots of cell outlier information
#' @keywords internal
.fit_one_sampleQC <- function(df, K=1, fit_params, seed) {
    # unpack inputs
    alpha       = fit_params$alpha
    em_iters    = fit_params$em_iters
    mcd_alpha   = fit_params$mcd_alpha
    mcd_iters   = fit_params$mcd_iters
    method      = fit_params$method
    track       = fit_params$track

    # extract x values
    x           = do.call(rbind, lapply(df$qc_metrics, as.matrix))
    rownames(x) = unlist(df$cell_id)
    qc_names    = colnames(x)

    # extract group info
    sample_list = df$sample_id
    sample_ids  = rep.int(
        df$sample_id, 
        vapply(df$cell_id, length, numeric(1))
        )
    groups      = as.integer(factor(sample_ids))
    groups_0    = groups - 1

    # extract info for cpp code
    D           = ncol(x)
    J           = length(sample_list)
    N           = nrow(x)

    # initialize clusters
    init_clusts = .init_clusts(x, groups, J, K, N)
    
    # run EM code
    if (method == 'robust') {
        message('running robust EM algorithm')
        fit_obj = fit_sampleQC_robust_cpp(
            x, init_clusts$init_z, groups_0, D, J, K, N, 
            em_iters, mcd_alpha, mcd_iters, track
            )

    } else if (method == 'mle') {
        message('running EM algorithm')
        fit_obj = fit_sampleQC_mle_cpp(
            x, init_clusts$init_gamma, groups_0, 
            D, J, K, N, em_iters)
    }

    # adjust mu_0, alpha_j and beta_k to have zero means
    alpha_means     = colMeans(fit_obj$alpha_j)
    beta_means      = colMeans(fit_obj$beta_k)
    fit_obj$mu_0    = as.vector(fit_obj$mu_0) + alpha_means + beta_means
    fit_obj$alpha_j = sweep(fit_obj$alpha_j, 2, alpha_means, "-")
    fit_obj$beta_k  = sweep(fit_obj$beta_k, 2, beta_means, "-")

    # process results
    fit_obj$sample_list = sample_list
    fit_obj$sample_id   = sample_ids
    fit_obj$qc_names    = qc_names

    # tidy some stuff
    fit_obj$like_1      = as.vector(fit_obj$like_1)
    fit_obj$like_2      = as.vector(fit_obj$like_2)

    # extract outliers
    fit_obj$outliers_dt = .calc_outliers_dt(fit_obj, x, groups, alpha)

    return(fit_obj)
}

#' Initializes clusters
#' 
#' @param x Matrix of QC values
#' @param groups sample_id labels as integers
#' @param J Number of samples
#' @param K Number of requested clusters
#' @param N Number of observations
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom mclust mclustBIC summaryMclustBIC
#'
#' @return list, containing lots of cell outlier information
#' @keywords internal
.init_clusts <- function(x, groups, J, K, N) {
    # K = 1 is special case
    if (K == 1) {
        init_z      = rep(0,N)
        init_gamma  = matrix(rep(0,N), ncol=1)
        return(list(init_gamma = init_gamma, init_z = init_z))
    }

    # centre x around group medians, for sensible GMM initialization
    grp_medians = vapply(
        seq_len(J), function(j) apply(x[groups == j, ], 2, median),
        numeric(ncol(x))) %>% t
    x_centred   = x - grp_medians[groups, ]

    # try to initialize clusters, multiple times if necessary
    n_tries     = 5
    n           = 1
    not_done    = TRUE
    message('fitting GMM to subset for initialization')
    while (n <= n_tries & not_done) {
        # take sample
        sample_idx  = sample(N, min(N, 1e4))
        x_sample    = x_centred[sample_idx, ]

        # fit mclust model
        model_spec  = "EVE"
        mclust_obj  = mclustBIC(x_sample, G = K, modelNames = model_spec,
            verbose = FALSE )

        # fit clusters
        mclust_init = summaryMclustBIC(mclust_obj, x_centred,
            G = K, modelNames = model_spec)
        K_mclust    = length(unique(mclust_init$classification))
        not_done    = K_mclust != K
        n           = n + 1
    }

    if (K_mclust == K) {
        message('initializing cluster membership')
        init_z      = mclust_init$classification - 1
        init_gamma  = mclust_init$z
    } else {
        warning(strwrap(prefix = " ", initial = "",
            "Mclust struggled to initialize, using random initialization instead.
            This may indicate that the selected value of K is too high."))        
        fake_data   = exp(matrix(rnorm(N*K), ncol = K))
        init_gamma  = fake_data / rowSums(fake_data)
        rownames(init_gamma) = rownames(x)
        init_z      = apply(init_gamma, 1, which.max) - 1
    }
    assert_that(all.equal(rowSums(init_gamma), rep(1, N), check.names = FALSE))
    assert_that(length(unique(init_z)) == K)

    return(list(init_gamma = init_gamma, init_z = init_z))
}

#' DEPRECATED: Calculate proportions of cells within each sample allocated 
#' to each 
#' 
#' @param fit_obj Most of output from .fit_one_sampleQC
#' @param x Matrix of QC values
#' @param groups Sample labels for every cell
#' @param alpha Chi squared threshold defining outliers
#' 
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table rbindlist dcast
#' @importFrom matrixStats rowMins
#'
#' @return data.table containing fitted parameters
#' @keywords internal
.calc_p_jk_proportions <- function(fit_obj, x, groups, alpha) {
    # unpack
    mu_0        = fit_obj$mu_0
    alpha_js    = fit_obj$alpha_j
    beta_ks     = fit_obj$beta_k
    sigma_ks    = fit_obj$sigma_k
    K           = fit_obj$K
    sample_list = fit_obj$sample_list

    # calc maha distance threshold for relevant chi-squared distn
    ndims       = length(fit_obj$qc_names)
    maha_cut    = qchisq(1 - alpha, df=ndims)

    # loop through groups
    groups_list = sort(unique(groups))
    message('calculating outliers for ', length(groups_list), ' groups:')
    outliers_dt = lapply(
        groups_list, function(j) {
            message('.', appendLF=FALSE)
            if (j%%20 == 0)
                message()
            # restrict to this group
            x_j     = x[ groups == j, ]
            n_j     = nrow(x_j)
            alpha_j = alpha_js[j, ]

            # calculate mahalanobis distances for each cell from each cluster
            mu_j    = mu_0 + alpha_j
            # mahas   = vapply(
            #     seq_len(K), 
            #     function(k)
            #         mvnfast::maha(x_j, mu_j + beta_ks[k,], sigma_ks[,,k]), 
            #         numeric(n_j)
            #         )
            mahas   = vapply(
                seq_len(K), 
                function(k) 
                    mahalanobis(x_j, mu_j + beta_ks[k,], sigma_ks[,,k]),
                numeric(n_j)
                )

            # find smallest, also label outliers
            x_cs    = paste0('C', apply(mahas, 1, which.min))
            out_idx = rowMins(mahas) > maha_cut
            x_cs[ out_idx ] = 'outlier'

            return(data.table(sample_id=sample_list[[j]], cluster=x_cs))
        }) %>% rbindlist
    message()
    outs_counts     = dcast(
        outliers_dt, sample_id ~ cluster, 
        value.var='cluster', fun.aggregate=length
        )
    p_jk_props      = as.matrix(outs_counts[, -'sample_id', with=FALSE])
    rownames(p_jk_props) = outs_counts$sample_id

    return(p_jk_props)
}

#' Calculate outliers for each cell
#' 
#' @param fit_obj Most of output from .fit_one_sampleQC
#' @param x Matrix of QC values
#' @param groups Sample labels for every cell
#' @param alpha Chi squared threshold defining outliers
#' 
#' @importFrom magrittr "%>%" set_colnames
#' @importFrom data.table data.table ":="
#' @importFrom matrixStats rowMins
#'
#' @return data.table containing fitted parameters
#' @keywords internal
.calc_outliers_dt <- function(fit_obj, x, groups, alpha) {
    # unpack
    mu_0        = fit_obj$mu_0
    alpha_js    = fit_obj$alpha_j
    beta_ks     = fit_obj$beta_k
    sigma_ks    = fit_obj$sigma_k
    K           = fit_obj$K
    sample_list = fit_obj$sample_list

    # calc maha distance threshold for relevant chi-squared distn
    maha_cut    = qchisq(1 - alpha, df=fit_obj$D)

    # centre x by alpha_js
    x_centred   = sweep(x, 2, mu_0, "-") - alpha_js[groups, ]
    n_x         = nrow(x)

    # calculate maha distance for each cluster component
    maha_mat    = vapply(
        seq_len(K), 
        function(k)
            mahalanobis(x_centred, beta_ks[k, ], sigma_ks[,,k]), 
        numeric(n_x)
        )
    maha_mins   = rowMins(maha_mat)

    # put together into data.table
    maha_names  = paste0('maha_', seq_len(K))
    outliers_dt = data.table(maha_mat) %>%
        set_colnames(maha_names) %>%
        .[, cell_id     := rownames(x) ] %>%
        .[, sample_id   := sample_list[groups] ] %>%
        .[, outlier     := maha_mins > maha_cut ]

    return(outliers_dt)
}

#' Puts results from fitting models into SingleCellExperiment format
#' 
#' @importFrom S4Vectors metadata "metadata<-"
#' @importFrom SummarizedExperiment "colData<-"
#' @importFrom assertthat assert_that
#' 
#' @return qc_obj with things added
#' @keywords internal
.assemble_fit_results <- function(qc_obj, fit_list, fit_params) {
    # put everything in metadata as lazy storage
    metadata(qc_obj)$fit_list   = fit_list
    metadata(qc_obj)$fit_params = fit_params

    # set up dummy columns
    sample_list = rownames(colData(qc_obj))
    n_samples   = nrow(colData(qc_obj))
    D           = metadata(qc_obj)$D
    colData(qc_obj)$K           = numeric(n_samples)
    colData(qc_obj)$mu_0        = 
        rep(list(matrix(NA, nrow=1, ncol=D)), n_samples) %>% 
        setNames(sample_list)
    colData(qc_obj)$alpha_j     = 
        rep(list(matrix(NA, nrow=1, ncol=D)), n_samples) %>% 
        setNames(sample_list)
    colData(qc_obj)$beta_k      = 
        rep(list(matrix(NA, nrow=1, ncol=D)), n_samples) %>% 
        setNames(sample_list)
    colData(qc_obj)$sigma_k     = 
        rep(list(matrix(NA, nrow=D, ncol=D)), n_samples) %>% 
        setNames(sample_list)
    colData(qc_obj)$p_jk        = 
        rep(list(matrix(NA, nrow=1, ncol=1)), n_samples) %>% 
        setNames(sample_list)
    colData(qc_obj)$z           = 
        rep(list(factor(NA)), n_samples) %>% 
        setNames(sample_list)
    colData(qc_obj)$outlier     = 
        rep(list(factor(NA)), n_samples) %>% 
        setNames(sample_list)

    # TODO: convert warning to error
    tryCatch({
        # loop through fitted models
        for (ii in seq_along(fit_list)) {
            # which samples to update?
            if ( fit_params$do_list == TRUE ) {
                # if we fit to each group, do this group
                g       = names(fit_list)[[ii]]
                g_idx   = qc_obj$group_id == g
    
                # get this fit_obj
                fit_obj = fit_list[[g]]

            } else {
                # otherwise update everything
                g_idx   = rep(TRUE, ncol(qc_obj))

                # get this fit_obj
                fit_obj = fit_list[[1]]
            }
            # how many here?
            n_ii    = sum(g_idx)

            # add K
            colData(qc_obj)$K[g_idx] = 
                rep(fit_obj$K, n_ii)
            # add mu_0
            colData(qc_obj)$mu_0[g_idx] = 
                rep(list(matrix(fit_obj$mu_0, nrow=1)), n_ii)
            # add alpha_j
            colData(qc_obj)$alpha_j[g_idx] = 
                asplit(fit_obj$alpha_j, 1) %>% 
                lapply(matrix, nrow=1)
            # add beta_k
            colData(qc_obj)$beta_k[g_idx] = 
                rep(list(fit_obj$beta_k), n_ii)
            # add sigma_k
            colData(qc_obj)$sigma_k[g_idx] = 
                rep(list(fit_obj$sigma_k), n_ii)
            # add p_jk
            colData(qc_obj)$p_jk[g_idx] = 
                asplit(fit_obj$p_jk, 1)
            # add z
            colData(qc_obj)$z[g_idx] = 
                split(as.vector(fit_obj$z), fit_obj$sample_id)
            # add outliers
            colData(qc_obj)$outlier[g_idx] = 
                split(fit_obj$outliers_dt, fit_obj$sample_id)
        }
    }, warning = function(w) 
        stop('problem in assembling results')
    )

    # check it worked ok
    assert_that( .check_is_qc_obj(qc_obj) == 'fit', 
        msg="fit didn't work somehow" )

    return(qc_obj)
}

#' Plots inferred sample-level statistics over QC biaxials
#'
#' @param qc_obj Output from function fit_sampleQC
#' @param sel_sample Which sample to plot?
#' @param qc_names List of metrics to be plotted (first is used for y axis).
#' Default is to use the qc_names specified in metadata(qc_obj)$qc_names; this
#' option is included to allow for "zooming in" on a subset, or changing the 
#' order in which they are plotted.
#' @param alpha_cut Chi-squared quantile to determine how large the ellipses
#' should be.
#' 
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr "%>%" set_colnames
#' @importFrom data.table copy ":=" setnames rbindlist data.table melt
#' @importFrom ggplot2 ggplot aes geom_bin2d geom_path geom_point
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_shape_manual
#' @importFrom ggplot2 scale_linetype_manual scale_fill_distiller
#' @importFrom ggplot2 coord_cartesian facet_grid theme_bw theme labs
#' @importFrom scales pretty_breaks
#' 
#' @return ggplot object
#' @export
plot_fit_over_biaxials_one_sample <- function(qc_obj, sel_sample, 
    qc_names=NULL, alpha_cut=0.01) {
    # can we do this?
    assert_that(
        .check_is_qc_obj(qc_obj) == 'fit',
        msg='SampleQC model must be fit (via `fit_sampleQC`) before calling 
        this function')

    # expect to use qc_names specified in qc_obj; qc_names parameter here
    # allows zooming in on subset of qc metrics used
    if (is.null(qc_names))
        qc_names    = metadata(qc_obj)$qc_names

    # get data for cells
    qc_1        = qc_names[[1]]
    qc_not_1    = qc_names[-1]
    points_dt   = copy(colData(qc_obj)[[sel_sample, 'qc_metrics']]) %>%
        .[, cell_id := colData(qc_obj)$cell_id[[sel_sample]] ] %>%
        data.table::melt(
            id      = c('cell_id', qc_1), 
            measure = qc_not_1,
            variable.name='other_qc', value.name='qc_x'
            ) %>%
        setnames(qc_1, 'qc_y')

    # extract ellipses
    ellipses_dt = lapply(
        2:length(qc_names),
        function(jj) 
            .calc_ellipses_dt(
                qc_obj, sel_sample, 1, jj, 
                qc_names[[jj]], alpha=alpha_cut
                )
        ) %>% rbindlist

    # extract means
    mu_0        = colData(qc_obj)[[sel_sample, 'mu_0']]
    alpha_j     = colData(qc_obj)[[sel_sample, 'alpha_j']]
    beta_k      = colData(qc_obj)[[sel_sample, 'beta_k']]
    sigma_k     = colData(qc_obj)[[sel_sample, 'sigma_k']]
    means_dt    = sweep(beta_k, 2, mu_0 + alpha_j, '+') %>%
        data.table %>%
        set_colnames(qc_names) %>%
        .[, component := 1:.N] %>%
        data.table::melt(id=c('component', qc_1), 
            variable.name='other_qc', value.name='qc_x') %>%
        setnames(qc_1, 'qc_y')

    # vars in nice order
    points_dt[, other_qc := factor(other_qc, levels=qc_not_1)]
    ellipses_dt[, other_qc := factor(other_qc, levels=qc_not_1)]
    means_dt[, other_qc := factor(other_qc, levels=qc_not_1)]

    # base range on points, not ellipses
    y_range = c(
        floor(min(points_dt$qc_y)*2)/2,
        ceiling(max(points_dt$qc_y)*2)/2
        )

    # plot
    g = ggplot() + 
        aes(y=qc_y, x=qc_x) +
        geom_bin2d( data=points_dt ) +
        geom_path( data=ellipses_dt, aes(group=component) ) +
        geom_point( data=means_dt, size=3) +
        scale_x_continuous( breaks=pretty_breaks() ) + 
        scale_y_continuous( breaks=pretty_breaks() ) +
        scale_shape_manual( values=c(1, 16) ) +
        scale_linetype_manual( values=c('dashed', 'solid') ) +
        scale_fill_distiller( palette='RdBu', trans='log10' ) +
        coord_cartesian( ylim=y_range ) +
        facet_grid( . ~ other_qc, scales='free_x' ) +
        theme_bw() + 
        # theme( aspect.ratio=1 ) +
        # labs( y=qc_1, x='other QC variable')
        labs( y=qc_1, x=NULL)

    return(g)
}

#' Calculates path for an ellipse based on 2D mu and sigma values
#'
#' @param qc_obj Output from function fit_sampleQC
#' @param sel_sample Which sample?
#' @param dim1,dim2 Dimensions to use
#' @param dim2_name Name of second dimension
#' @param alpha Chi-squared quantile to use for ellipse sizes
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table rbindlist ":="
#' 
#' @return data.table object
#' @keywords internal
.calc_ellipses_dt <- function(qc_obj, sel_sample, dim1, dim2, 
    dim2_name, alpha=0.01) {
    # extract data from qc_obj
    K           = colData(qc_obj)[sel_sample, 'K']
    mu_0        = colData(qc_obj)[[sel_sample, 'mu_0']]
    alpha_j     = colData(qc_obj)[[sel_sample, 'alpha_j']]
    beta_k      = colData(qc_obj)[[sel_sample, 'beta_k']]
    sigma_k     = colData(qc_obj)[[sel_sample, 'sigma_k']]

    # add together into means
    mus         = sweep(beta_k, 2, mu_0 + alpha_j, '+')

    # # check that we have valid values
    # this_m  = mus[[ sel_sample ]]
    # if ( is.null(this_m) )
    #     return(NULL)

    # specify dimensions to use
    sel_idx     = c(dim1, dim2)

    # go through samples and components
    ellipses_dt = lapply(
        seq_len(K),
        function(k) .calc_one_ellipse(
            mus[k, sel_idx], 
            sigma_k[sel_idx, sel_idx, k], 
            alpha
            ) %>% .[, component := k ]
        ) %>% rbindlist %>% .[, other_qc := dim2_name ]

    return(ellipses_dt)
}

#' Calculates path for an ellipse based on 2D mu and sigma values
#'
#' @param mu mean vector
#' @param sigma covariance matrix
#' @param alpha contour threshold
#' 
#' @importFrom magrittr "%>%"
#' @importFrom mixtools ellipse
#' @importFrom data.table data.table setnames ":="
#' 
#' @return data.table object
#' @keywords internal
.calc_one_ellipse <- function(mu, sigma, alpha) {
    ell = ellipse(mu, sigma, alpha=alpha, draw=FALSE)
    dt  = data.table(ell) %>%
        setnames(., names(.), c('qc_y', 'qc_x'))
    return(dt)
}

#' Calculates path for an ellipse based on 2D mu and sigma values
#'
#' @param qc_obj Output from function fit_sampleQC
#' @param sel_sample Selected sample
#' @param outliers_dt Optional alternative outlier specification; see function
#' \code{get_outliers}.
#'
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom data.table ":=" setnames rbindlist dcast copy melt
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 facet_grid theme_bw theme labs
#' @importFrom scales pretty_breaks
#' 
#' @return ggplot object
#' @export
plot_outliers_one_sample <- function(qc_obj, sel_sample, outliers_dt=NULL) {
    # unpack, arrange
    qc_names    = metadata(qc_obj)$qc_names
    qc_1        = qc_names[[1]]
    qc_not_1    = qc_names[-1]
    qc_sel      = copy(colData(qc_obj)[[sel_sample, 'qc_metrics']])
    # assert_that( all( qc_sel$cell_id == outlier_sel$cell_id ) )

    # check which outliers to use
    if (is.null(outliers_dt)) {
        # join together
        plot_dt     = cbind(
            cell_id = qc_obj$cell_id[[sel_sample]],
            qc_sel,
            outlier = qc_obj$outlier[[sel_sample]]$outlier
            ) %>%
            melt(
                id      = c('cell_id', 'outlier', qc_1), 
                measure = qc_not_1,
                variable.name='qc_x_name', value.name='qc_x'
                ) %>%
            setnames(qc_1, 'qc_y')        
    } else {
        # join together
        plot_dt     = cbind(
            cell_id = qc_obj$cell_id[[sel_sample]],
            qc_sel
            ) %>% 
            merge(outliers_dt[, .(cell_id, outlier)], by='cell_id', 
                all.x=TRUE, all.y=FALSE) %>%
            melt(
                id      = c('cell_id', 'outlier', qc_1), 
                measure = qc_not_1,
                variable.name='qc_x_name', value.name='qc_x'
                ) %>%
            setnames(qc_1, 'qc_y')        
    }

    # plot
    g   = ggplot() + 
        aes(y=qc_y, x=qc_x, colour=outlier) +
        geom_point(data=plot_dt[outlier == FALSE], size=1, colour='grey' ) +
        geom_point(data=plot_dt[outlier == TRUE],  size=1, colour='red' ) +
        scale_x_continuous( breaks=pretty_breaks() ) + 
        scale_y_continuous( breaks=pretty_breaks() ) +
        facet_grid( . ~ qc_x_name, scales='free_x' ) +
        theme_bw() + 
        # theme( aspect.ratio=1 ) +
        labs( y=qc_names[[1]], x='other QC variable', colour='outlier?' )

    return(g)
}

#' Plots inferred sample-level statistics, with empirical log-likelihood
#'
#' @param qc_obj Output from function fit_sampleQC
#' @param sel_group Which sample group to plot
#' @param df degrees of freedom for multivariate t
#' 
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table ":=" setnames melt
#' @importFrom forcats fct_reorder
#' @importFrom MASS cov.trob
#' @importFrom mvnfast dmvt
#' @importFrom ggplot2 ggplot aes geom_tile geom_text
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 theme_bw theme labs element_blank
#' @importFrom scales pretty_breaks
#' @importFrom patchwork plot_layout
#' 
#' @return ggplot object
#' @export
plot_alpha_js_likelihoods <- function(qc_obj, sel_group, df=5) {
    # unpack
    qc_names    = metadata(qc_obj)$qc_names
    g_idx       = colData(qc_obj)$group_id == sel_group
    samples     = colData(qc_obj)$sample_id[g_idx]
    mu_0        = do.call(rbind, colData(qc_obj)$mu_0[g_idx])
    alpha_j     = do.call(rbind, colData(qc_obj)[g_idx, 'alpha_j'])

    # some checks
    assert_that( all( dim(mu_0) == dim(alpha_j) ),
        msg='mu_0 and alpha_j should have the same dimensions')

    # extract values
    alphas_dt   = data.table(mu_0 + alpha_j) %>%
        set_colnames(qc_names) %>% 
        .[, sample_id := samples ] %>%
        .[, sample_id := fct_reorder(sample_id, -alpha_j[, 1]) ] %>%
        melt(
            id='sample_id', 
            variable.name='qc_metric', value.name='qc_value'
            ) %>%
        .[, scaled_qc := scale(qc_value), by='qc_metric' ]

    # max value
    max_val     = ceiling(10*max(abs(alphas_dt$scaled_qc)))/10

    # make plot
    g_alphas    = ggplot(alphas_dt) +
        aes(
            x=qc_metric, y=sample_id, fill=scaled_qc,
            label=sprintf("%0.2f", round(qc_value, digits = 2))
            ) +
        geom_tile() + geom_text() +
        scale_fill_distiller(
            palette='RdBu', limits=c(-max_val,max_val),
            breaks=pretty_breaks()
            ) +
        labs(
            x='QC metric', y='sample ID', 
            fill='scaled\nQC metric\nvalue\n'
            ) +
        theme_bw()

    # fit robust multivariate t to alpha values, extract likelihoods
    fit_robt    = cov.trob(alpha_j, cor=TRUE, nu=df)
    llikes      = dmvt(alpha_j, 
        mu=fit_robt$center, sigma=fit_robt$cov, 
        df=df, log=TRUE)
    likes_dt    = data.table(
        loglike     = llikes,
        sample_id   = factor(samples, levels=levels(alphas_dt$sample_id))
        )

    # plot likelihoods
    g_likes     = ggplot(likes_dt) +
        aes( x="", y=sample_id, fill=loglike ) +
        geom_tile() +
        scale_fill_distiller(
            palette='PiYG', direction=1,
            breaks=pretty_breaks()
            ) +
        labs( x=NULL, y=NULL, fill='empirical\nlog-likelihood\n' ) +
        theme_bw() + theme( axis.text.y=element_blank() )

    # assemble via patchwork
    g = g_alphas + g_likes + plot_layout( widths=c(4, 1) )

    return(g)
}

#' Plots inferred sample-level statistics, with empirical log-likelihood
#'
#' @param qc_obj Output from function fit_sampleQC
#' @param sel_group Which sample group to plot
#' @param df Degrees of freedom for multivariate t
#' @param qc_idx Which QC metrics to plot?
#' @param pc_idx Which PCA components to plot?
#' 
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr "%>%" set_colnames
#' @importFrom data.table data.table ":="
#' @importFrom MASS cov.trob
#' @importFrom mvnfast dmvt
#' @importFrom mixtools ellipse
#' @importFrom ggplot2 ggplot aes_string geom_point geom_path
#' @importFrom ggplot2 scale_fill_distiller theme_bw theme labs
#' @importFrom scales pretty_breaks
#' @importFrom patchwork plot_layout
#' 
#' @return ggplot object
#' @export
plot_alpha_js <- function(qc_obj, sel_group, df=5, qc_idx=c(1,2), pc_idx=c(1,2)) {
    # unpack
    qc_names    = metadata(qc_obj)$qc_names
    D           = metadata(qc_obj)$D
    g_idx       = colData(qc_obj)$group_id == sel_group
    samples     = colData(qc_obj)$sample_id[g_idx]
    mu_0        = do.call(rbind, colData(qc_obj)$mu_0[g_idx])
    alpha_j     = do.call(rbind, colData(qc_obj)[g_idx, 'alpha_j'])

    # fit robust multivariate t to alpha values, extract likelihoods
    fit_robt    = cov.trob(alpha_j, cor=TRUE, nu=df)
    mu          = fit_robt$center
    cov_mat     = fit_robt$cov
    llikes      = dmvt(alpha_j, mu=mu, sigma=cov_mat, df=df, log=TRUE)

    # get PCs from covariance matrix
    cov_eigen   = eigen(cov_mat)
    pc_scores   = sweep(alpha_j, 2, fit_robt$center, "-") %*% 
        cov_eigen$vectors
    std_devs    = sqrt(cov_eigen$values)

    # put together into dt
    plot_dt     = data.table(pc_scores) %>%
        set_colnames(paste0('PC', seq_len(D))) %>%
        .[, (qc_names) := asplit(alpha_j, 2) ] %>%
        .[, sample_id := samples] %>%
        .[, loglike   := llikes ]

    # make ellipse
    alpha       = 0.05
    ellipse_std = data.table(
        ellipse(mu[qc_idx], cov_mat[qc_idx, qc_idx], alpha=alpha, draw=FALSE)
        ) %>% set_colnames(qc_names[qc_idx])
    ellipse_pca = data.table(
        ellipse(rep(0,2), diag(std_devs[pc_idx]^2), alpha=alpha, draw=FALSE)
        ) %>% set_colnames(paste0('PC', pc_idx))

    # plot
    g_std = ggplot() +
        aes_string( x=qc_names[qc_idx[1]], y=qc_names[qc_idx[2]] ) +
        geom_point(
            data=plot_dt, 
            aes(fill=loglike), 
            colour='black', shape=21, size=3
            ) +
        geom_path( data=ellipse_std ) +
        scale_fill_distiller(
            palette='PiYG', direction=1,
            breaks=pretty_breaks()
            ) +
        labs( fill='empirical\nlog-likelihood\n' ) +
        theme_bw() + theme(aspect.ratio=1)
    g_pca = ggplot() +
        aes_string(
            x=paste0('PC', pc_idx[[1]]),
            y=paste0('PC', pc_idx[[2]])
            ) +
        geom_point(
            data=plot_dt, aes(fill=loglike), 
            colour='black', shape=21, size=3 ) +
        geom_path( data=ellipse_pca ) +
        scale_fill_distiller(
            palette='PiYG', direction=1,
            breaks=pretty_breaks()
            ) +
        labs( fill='empirical\nlog-likelihood\n' ) +
        # coord_cartesian( xlim=c(-1, 1), ylim=c(-1, 1) ) +
        theme_bw() + theme(aspect.ratio=1)

    # assemble via patchwork
    g = g_std + g_pca

    return(g)
}

#' Plots estimates of cluster splits by sample, with KL divergence from mean
#'
#' @param qc_obj Output from function fit_sampleQC
#' @param sel_group Which sample group to plot
#' 
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table ":=" setnames melt copy merge.data.table
#' @importFrom forcats fct_reorder
#' @importFrom MASS cov.trob
#' @importFrom mvnfast dmvt
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_distiller
#' @importFrom ggplot2 theme_bw theme labs element_blank
#' @importFrom scales pretty_breaks
#' @importFrom patchwork plot_layout
#' 
#' @return ggplot object
#' @export
plot_beta_ks <- function(qc_obj, sel_group) {
    # unpack
    g_idx       = colData(qc_obj)$group_id == sel_group
    alpha_j     = do.call(rbind, colData(qc_obj)$alpha_j[g_idx])
    sample_list = colData(qc_obj)$sample_id[ g_idx ]
    sample_ord  = sample_list[order(-alpha_j[,1])]
    z           = do.call(c, colData(qc_obj)$z[g_idx])
    outliers    = do.call(rbind, colData(qc_obj)$outlier[g_idx])

    # some checks
    assert_that( nrow(outliers) == length(z),
        msg="outliers and z don't match up :/" )

    # put dt together
    outliers_dt = copy(outliers) %>%
        .[, sample_id := factor(sample_id, levels=sample_ord)] %>%
        .[, z := paste0(sel_group, '.', z) ] %>%
        .[, cluster := z] %>%
        .[ outlier == TRUE, cluster := 'out']

    # calculate proportions
    p_jk_dt     = outliers_dt[, .N, by=c('sample_id', 'z')] %>% 
        .[, sum_j   := sum(N), by='sample_id' ] %>%
        .[, p_jk    := N / sum_j ] %>%
        .[, p_jk_label := sprintf("%0.2f", round(p_jk, digits = 2)) ]
    p_k_dt      = outliers_dt[, .N, by='z'] %>% 
        .[, total := sum(N) ] %>%
        .[, p_k := N / total ] %>%
        setorder('z')
    p_k         = p_k_dt$p_k

    # make into matrix
    p_jk_mat    = dcast(p_jk_dt, sample_id ~ z, value.var='p_jk') %>%
        .[, -'sample_id'] %>% 
        as.matrix

    # define KL divergence fn
    kl_fn   <- function(Q) 
        sum(p_k * (log2(p_k) - log2(Q)))

    # make data.table with kl vals in
    kl_dt       = data.table(
        sample_id   = factor(sample_list, levels=sample_ord),
        kl          = apply(p_jk_mat, 1, kl_fn)
        )

    # extract p_jk values
    diffs_dt    = data.table(sweep(p_jk_mat, 2, p_k, "-")) %>% 
        .[, sample_id := factor(sample_list, levels=sample_ord) ] %>%
        melt(
            id='sample_id', variable.name='z', 
            value.name='cluster_delta'
            )
    p_jk_dt     = data.table::merge.data.table(
        diffs_dt, p_jk_dt, by=c('sample_id', 'z')
        )

    # find limits for deltas
    delta_lim   = ceiling(max(abs(diffs_dt$cluster_delta))*10)/10

    # extract outlier proportions
    outs_dt     = outliers_dt[, .(outs=sum(outlier), total=.N), 
        by=c('sample_id')] %>%
        .[, out_prop := outs / total ] %>%
        .[, out_label := sprintf("%0.3f", round(out_prop, digits = 3))]

    # make plot
    g_pjks  = ggplot(p_jk_dt) +
        aes(x=z, y=sample_id, fill=cluster_delta, label=p_jk_label) +
        geom_tile() + geom_text() +
        scale_fill_distiller(
            palette='RdBu', breaks=pretty_breaks(),
            limits=c(-delta_lim, delta_lim)
            ) +
        labs(
            x='GMM component', y='sample ID', 
            fill='difference\nfrom mean\ncluster propn.\n'
            ) +
        theme_bw() + theme( legend.position='bottom', legend.direction='vertical' )

    # make plot
    g_outs  = ggplot(outs_dt) +
        aes( x="", y=sample_id, fill=out_prop, label=out_label ) +
        geom_tile() + geom_text() +
        scale_fill_distiller(
            palette='BuPu', direction=1,
            breaks=pretty_breaks()
            ) +
        expand_limits( fill=0 ) +
        labs( x=NULL, y=NULL, fill='propn.\ncells\nexcluded' ) +
        theme_bw() + 
        theme(
            legend.position='bottom', legend.direction='vertical', 
            axis.text.y=element_blank()
            )

    # plot likelihoods
    g_kl    = ggplot(kl_dt) +
        aes( x="", y=sample_id, fill=kl ) +
        geom_tile() +
        scale_fill_distiller(palette='PiYG', direction=-1,
            breaks=pretty_breaks()
            ) +
        expand_limits( fill=0 ) +
        labs(x=NULL, y=NULL,
            fill='KL divergence\nfrom mean\nproportion\n') +
        theme_bw() + 
        theme(legend.position='bottom', legend.direction='vertical',
            axis.text.y=element_blank())

    # assemble via patchwork
    g = g_pjks + g_outs + g_kl + plot_layout( widths=c(5, 2, 2) )

    return(g)
}

