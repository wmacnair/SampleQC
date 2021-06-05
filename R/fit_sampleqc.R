#' Fits `SampleQC` model to one cluster of samples
#' 
#' @description 
#' Running \link{calc_pairwise_mmds} identifies groups 
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
#' \link{calc_pairwise_mmds}. Next, we recommend calling 
#' \link{fit_sampleqc} with \code{K_list=rep(1, get_n_groups(qc_obj)))},
#' then rendering outputs with \code{make_sampleqc_report(qc_obj, save_dir, 
#' 'test')}. This fits the simplest possible model to each sample group (i.e. 
#' one component), and renders a report with biaxial plots. 
#' 
#' The user can then use the biaxial plots to determine how many components 
#' are appropriate for each sample group. So if 
#' 
#' @param qc_obj Output from calc_pairwise_mmds
#' @param K_all,K_list How many QC celltypes do we expect? Exactly one of K_all
#' and K_list should be specified. If the user wants to fit the same model to 
#' all samples they should specify `K_all` as an integer. If the user wants to 
#' fit a different model to each of the sample groups identified by 
#' `calc_pairwise_mmds`, they should specify `K_list` as an 
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
fit_sampleqc <- function(qc_obj, K_all=NULL, K_list=NULL, n_cores=NULL, 
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
                .fit_one_sampleqc(df_list[[i]], K_list[[i]], fit_params, seed = i), 
                BPPARAM = bpparam)
        names(fit_list) = metadata(qc_obj)$group_list
        bpstart(bpstop)
    } else {
        # fit one model to all samples
        df          = colData(qc_obj)
        fit_list    = list(.fit_one_sampleqc(df, K_all, fit_params, seed = 1))
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
            time; see ?fit_sampleqc')
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
.fit_one_sampleqc <- function(df, K=1, fit_params, seed) {
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
        fit_obj = fit_sampleqc_robust_cpp(
            x, init_clusts$init_z, groups_0, D, J, K, N, 
            em_iters, mcd_alpha, mcd_iters, track
            )

    } else if (method == 'mle') {
        message('running EM algorithm')
        fit_obj = fit_sampleqc_mle_cpp(
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
#' @param fit_obj Most of output from .fit_one_sampleqc
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
#' @param fit_obj Most of output from .fit_one_sampleqc
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
