# SampleQC: robust multivariate, multi-celltype, multi-sample quality control 
# for single cell RNA-seq
#' devtools::load_all('~/work/packages/BayesQC')
#' devtools::document('~/work/packages/BayesQC')

#' SampleQC_cells.R
#' Code to identify outlier cells within sample groups identified by code
#' in SampleQC_samples. Includes plots of outputs, and plots of fitted 
#' statistics

#' Fits `SampleQC` model to one cluster of samples
#' 
#' @param mmd_list Outputs from calculate_sample_to_sample_MMDs
#' @param qc_dt Data.table of QC metrics for all cells and samples
#' @param qc_names List of metrics to actually use for calculating sample-to-
#' sample distances
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
#' means and covariances
#' @param method Which of various implemented options should be used?
#' 
#' @section Details:
#'
#' [How to define K]
#'
#' [Start with `K_all=1`]
#'
#' [Add more details of K_all vs K_list?]
#'
#' [What difference do mcd parameters make?]
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%" set_rownames set_colnames
#' @importFrom mclust mclustBIC
#' @importFrom data.table setkey copy ":="
#'
#' @return list, containing MMD values and sample clusters based on MMD values
#' @export
fit_sampleQC <- function(mmd_list, qc_dt, qc_names=c('log_counts', 'log_feats',
    'logit_mito'), K_all=NULL, K_list=NULL, n_cores, alpha=0.01, em_iters=50, 
    mcd_alpha=0.5, mcd_iters=50, method=c('robust', 'mle')) {
    # check some inputs
    method      = match.arg(method)

    # check specification of type of run is ok
    if (!is.null(K_all) & !is.null(K_list))
        stop('only one of `K_all` and `K_list` should be specified at one 
            time; see ?fit_sampleQC')
    list_flag   = !is.null(K_list)

    # check specification of number of components is ok for both options
    if (list_flag) {
        N_clusts    = length(unique(mmd_list$mmd_clusts))
        assert_that(
            all.equal(K_list, as.integer(K_list)), 
            all(K_list > 0),
            msg     = 'K_list must be a vector of integers greater than 0')
        assert_that(
            length(K_list) == N_clusts, 
            msg     = 'K_list must have the same length as the number of 
            clusters in the mmd_clusts element of mmd_list')
    } else {
        assert_that(
            K_all == as.integer(K_all), 
            K_all > 0, 
            msg     = 'K_all must be an integer greater than 0')
    }
    # check inputs
    if (list_flag) {
        if (missing(n_cores))
            n_cores     = length(K_list)        
    } else {
        n_cores     = 1
    }

    if (list_flag) {
        # add clusters to qc_dt
        clusts_dt   = data.table(
            sample_id   = rownames(mmd_list$mmd_mat), 
            QC_clust    = paste0('QC',mmd_list$mmd_clusts)
            )
        qc_all      = data.table:::merge.data.table(clusts_dt, qc_dt,
            by='sample_id')

        # split data
        qc_dt_list  = split(qc_all, qc_all$QC_clust)
    } else {
        # don't split, but put into one list
        qc_all      = copy(qc_dt) %>% .[, QC_clust := 'All' ]
        qc_dt_list  = list(All=qc_all)
        K_list      = K_all
    }

    # fit each one on different core
    em_list         = mclapply(seq_along(qc_dt_list), 
        function(i) .fit_one_sampleQC(
            qc_dt_list[[i]], qc_names, K_list[[i]], 
            alpha, em_iters, mcd_alpha, mcd_iters, method
            ), mc.cores=n_cores)
    # em_list      = lapply(seq_along(qc_dt_list), 
    #     function(i) .fit_one_sampleQC(qc_dt_list[[i]], qc_names, 
    # K_list[[i]], alpha, em_iters, mcd_alpha, mcd_iters, method))

    names(em_list)   = names(qc_dt_list)
    
    return(em_list)
}

#' Fits `SampleQC` model to one cluster of samples
#' 
#' @param qc_dt Data.table of QC metrics for all cells and samples
#' @param qc_names List of metrics to actually use for calculating 
#' sample-to-sample distances
#' @param K How many QC celltypes are there?
#' @param alpha Chi squared threshold to define outliers
#' @param em_iters Maximum number of EM iterations
#' @param mcd_alpha,mcd_iters Parameters for robust estimation of celltype 
#' means and covariances
#' @param method Which of various implemented options should be used?
#' 
#' @section Details:
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%" set_rownames set_colnames
#' @importFrom mclust mclustBIC summaryMclustBIC
#' @importFrom data.table setkey
#'
#' @return list, containing lots of cell outlier information
#' @keywords internal
.fit_one_sampleQC <- function(qc_dt, qc_names=c('log_counts', 'log_feats', 
    'logit_mito'), K=1, alpha=0.01, em_iters=50, mcd_alpha=0.5, mcd_iters=50, 
    method=c('robust', 'mle')) {
    # checks of inputs
    # assert_that( is.integer(K) & K>0 )
    # assert_that( is.integer(em_iters) & em_iters>0 )
    method      = match.arg(method)
    assert_that( all(qc_names %in% names(qc_dt)) )

    # extract x values
    x           = as.matrix(qc_dt[, qc_names, with=FALSE])
    rownames(x) = qc_dt$cell_id

    # extract group info
    sample_list = sort(unique(qc_dt$sample_id))
    groups      = as.integer(factor(qc_dt$sample_id))
    groups_0    = groups - 1

    # extract info for cpp code
    D           = length(qc_names)
    J           = length(sample_list)
    N           = nrow(x)

    # centre x around group medians, for sensible GMM initialization
    grp_medians = vapply(
        seq_len(J), function(j) apply(x[groups == j, ], 2, median),
        numeric(ncol(x))) %>% t
    x_centred   = x - grp_medians[groups, ]

    # initialize via GMM over whole dataset
    message('fitting GMM to subset for initialization')
    model_spec  = "VVV"
    sample_idx  = sample(N, min(N, 1e4))
    x_sample    = x_centred[sample_idx, ]
    mclust_obj  = mclustBIC(
        x_sample, G=K, modelNames=model_spec, 
        verbose=FALSE
        )
    
    # run EM code
    if (method == 'robust') {
        if (K == 1) {
            init_z      = rep(0,N)
        } else {
            message('initializing cluster membership')
            init_z      = summaryMclustBIC(mclust_obj, x_centred,
                G=K, modelNames=model_spec)$classification-1
        }
        message('running robust EM algorithm')
        em_obj      = fit_sampleQC_robust_cpp(
            x, init_z, groups_0, D, J, K, N, 
            em_iters, mcd_alpha, mcd_iters
            )

    } else if (method == 'mle') {
        if (K == 1) {
            init_gamma  = matrix(rep(0,N), ncol=1)
        } else {
            message('initializing cluster membership')
            init_gamma  = summaryMclustBIC(mclust_obj, x_centred, 
                G=K, modelNames=model_spec)$z
        }

        message('running EM algorithm')
        em_obj     = fit_sampleQC_mle_cpp(
            x, init_gamma, groups_0, 
            D, J, K, N, em_iters)
    }

    # adjust mu_0, alpha_j and beta_k to have zero means
    alpha_means     = colMeans(em_obj$alpha_j)
    beta_means      = colMeans(em_obj$beta_k)
    em_obj$mu_0     = as.vector(em_obj$mu_0) + alpha_means + beta_means
    em_obj$alpha_j  = sweep(em_obj$alpha_j, 2, alpha_means, "-")
    em_obj$beta_k   = sweep(em_obj$beta_k, 2, beta_means, "-")

    # process results
    em_obj$sample_list  = sample_list
    em_obj$qc_names     = qc_names
    for (n in c("llike", "hyper_dt", "metadata"))
        em_obj[[n]]     = "blank"
    em_obj$sample_id    = .make_mu_jk_dt(em_obj)

    # tidy some stuff
    em_obj$like_1   = as.vector(em_obj$like_1)
    em_obj$like_2   = as.vector(em_obj$like_2)

    # extract mu and sigma values
    setkey(em_obj$sample_id, 'qc_metric')
    mu_list     = lapply(
        em_obj$sample_list, function(sel_sample) {
            lapply(seq_len(em_obj$K), function(k) {
                temp_dt     = em_obj$sample_id[ (sample_id == sel_sample) & 
                    (component == k) ]
                temp_dt     = temp_dt[qc_names]
                vals        = temp_dt$value %>% setNames(qc_names)
                return(vals)
            }) }) %>% setNames(sample_list)
    sigma_list  = lapply(
        em_obj$sample_list, function(s) {
            lapply(seq_len(em_obj$K), function(k) 
                em_obj$sigma_k[ , , k ] %>% 
                    set_rownames(qc_names) %>% 
                    set_colnames(qc_names))
        }) %>% setNames(sample_list)
    em_obj$mu_list      = mu_list
    em_obj$sigma_list   = sigma_list

    # extract outliers
    em_obj$p_jk_outliers    = .calc_p_jk_proportions(em_obj, x, groups, alpha)
    em_obj$outliers_dt      = .calc_outliers_dt(em_obj, x, groups, alpha)
    em_obj$outliers_plot    = .calc_outliers_plot(qc_dt, em_obj)

    return(em_obj)
}

#' Extract parameters from fitted em_obj
#' 
#' @param em_obj Most of output from .fit_one_sampleQC
#' 
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table ":="
#'
#' @return data.table containing fitted parameters
#' @keywords internal
.make_mu_jk_dt <- function(em_obj) {
    qc_names    = em_obj$qc_names
    mu_0    = data.table(
        qc_metric   = qc_names, 
        mu_0        = as.vector(em_obj$mu_0)
        )
    alphas  = data.table(
        qc_metric   = rep(qc_names, each=em_obj$J), 
        sample_id   = rep(em_obj$sample_list, times=em_obj$D), 
        alpha_j     = as.vector(em_obj$alpha_j)
        )
    betas   = data.table(
        qc_metric   = rep(qc_names, each=em_obj$K), 
        component   = rep(seq_len(em_obj$K), times=em_obj$D), 
        beta_k      = as.vector(em_obj$beta_k)
        )
    mu_jk_dt    = data.table:::merge.data.table(alphas, betas,
        by=c('qc_metric'), allow.cartesian=TRUE) %>% 
        data.table:::merge.data.table(mu_0,
            by=c('qc_metric'), allow.cartesian=TRUE) %>% 
        .[, value := mu_0 + alpha_j + beta_k ] %>%
        .[, stat := "mean" ]
    return(mu_jk_dt)
}

#' Calculate proportions of cells within each sample allocated to each 
#' 
#' @param em_obj Most of output from .fit_one_sampleQC
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
.calc_p_jk_proportions <- function(em_obj, x, groups, alpha) {
    # unpack
    mu_0        = em_obj$mu_0
    alpha_js    = em_obj$alpha_j
    beta_ks     = em_obj$beta_k
    sigma_ks    = em_obj$sigma_k
    K           = em_obj$K
    sample_list = em_obj$sample_list

    # calc maha distance threshold for relevant chi-squared distn
    ndims       = length(em_obj$qc_names)
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
#' @param em_obj Most of output from .fit_one_sampleQC
#' @param x Matrix of QC values
#' @param groups Sample labels for every cell
#' @param alpha Chi squared threshold defining outliers
#' 
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table setnames ":="
#' @importFrom matrixStats rowMins
#'
#' @return data.table containing fitted parameters
#' @keywords internal
.calc_outliers_dt <- function(em_obj, x, groups, alpha) {
    # unpack
    mu_0        = em_obj$mu_0
    alpha_js    = em_obj$alpha_j
    beta_ks     = em_obj$beta_k
    sigma_ks    = em_obj$sigma_k
    K           = em_obj$K
    sample_list = em_obj$sample_list

    # calc maha distance threshold for relevant chi-squared distn
    ndims       = length(em_obj$qc_names)
    maha_cut    = qchisq(1 - alpha, df=ndims)

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
        setnames(names(.), maha_names) %>%
        .[, cell_id     := rownames(x) ] %>%
        .[, sample_id   := sample_list[groups] ] %>%
        .[, outlier     := maha_mins > maha_cut ]

    return(outliers_dt)
}

#' Make data.table for plotting outliers over QC metrics
#' 
#' @param qc_dt data.table of QC metrics
#' @param em_obj Most of output from .fit_one_sampleQC
#' 
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table
#' @importFrom data.table melt
#' @importFrom data.table setnames
#'
#' @return data.table containing fitted parameters
#' @keywords internal
.calc_outliers_plot <- function(qc_dt, em_obj) {
    outliers_dt = em_obj$outliers_dt
    qc_names    = em_obj$qc_names
    outliers_plot   = outliers_dt[, list(cell_id, sample_id, outlier)] %>%
        data.table:::merge.data.table(
            qc_dt[, c('cell_id', qc_names), with=FALSE], on='cell_id'
            ) %>%
        melt(
            id      = c('cell_id', 'sample_id', 'outlier', qc_names[[1]]),
            measure = qc_names[-1], 
            variable.name='qc_x_name', value.name='qc_x'
            ) %>%
        setnames(qc_names[[1]], 'qc_y')
    return(outliers_plot)
}

#' Plots inferred sample-level statistics over QC biaxials
#'
#' @param qc_dt data.table of QC metrics
#' @param fit_obj fit_obj object from function fit_qc_model
#' @param fit_list list of fit_obj objects (for comparing multiple models)
#' @param sel_sample which sample to plot?
#' @param qc_names list of metrics to be plotted (first is used for y axis)
#' @param alpha_cut contour value for ellipses
#' @importFrom magrittr "%>%"
#' @importFrom data.table ":=" setnames rbindlist dcast
#' @importFrom ggplot2 ggplot aes geom_bin2d geom_path geom_point
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_shape_manual
#' @importFrom ggplot2 scale_linetype_manual scale_fill_distiller
#' @importFrom ggplot2 coord_cartesian facet_grid theme_bw theme labs
#' @importFrom scales pretty_breaks
#' @return ggplot object
#' @export
plot_fit_over_biaxials_one_sample <- function(fit_obj, qc_dt, sel_sample, 
    qc_names=NULL, alpha_cut=0.01) {
    # if qc_names not specified, use those in sample
    if (is.null(qc_names))
        qc_names    = fit_obj$qc_names

    # get points
    qc_1        = qc_names[[1]]
    qc_not_1    = qc_names[-1]
    points_dt   = qc_dt[ sample_id == sel_sample ] %>%
        melt(
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
                fit_obj, sel_sample, 1, jj, 
                qc_names[[jj]], alpha=alpha_cut
                )
        ) %>% rbindlist

    # extract means
    means_dt    = fit_obj$sample_id[
        (stat == 'mean') & (sample_id == sel_sample)
        ] %>%
        dcast(component ~ qc_metric, value.var='value') %>%
        melt(id=c('component', qc_1), 
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
        geom_bin2d(data=points_dt) +
        geom_path(
            data    = ellipses_dt[ sample_id == sel_sample ],
            aes(group=component)
            ) +
        geom_point(data=means_dt, size=3) +
        scale_x_continuous( breaks=pretty_breaks() ) + 
        scale_y_continuous( breaks=pretty_breaks() ) +
        coord_cartesian( ylim=y_range ) +
        scale_shape_manual( values=c(1, 16) ) +
        scale_linetype_manual( values=c('dashed', 'solid') ) +
        scale_fill_distiller( palette='RdBu', trans='log10' ) +
        facet_grid( . ~ other_qc, scales='free_x' ) +
        theme_bw() + 
        # theme( aspect.ratio=1 ) +
        # labs( y=qc_1, x='other QC variable')
        labs( y=qc_1, x=NULL)

    return(g)
}

#' Calculates path for an ellipse based on 2D mu and sigma values
#'
#' @param fit_obj fit_obj object from function fit_qc_model
#' @param dim1 dimension to use
#' @param dim2 dimension to use
#' @param alpha contour threshold
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table rbindlist ":="
#' @return data.table object
#' @keywords internal
.calc_ellipses_dt <- function(fit_obj, sel_sample, dim1, dim2, 
    dim2_name, alpha=0.01) {
    # extract data from fit_obj
    mus     = fit_obj$mu_list
    sigmas  = fit_obj$sigma_list
    n_comps = length(mus[[1]])

    # check that we have valid values
    this_m  = mus[[ sel_sample ]]
    if ( is.null(this_m) )
        return(NULL)

    # specify dimensions to use
    sel_idx = c(dim1, dim2)

    # go through samples and components
    ellipses_dt = lapply(
        seq_len(n_comps),
        function(c) .calc_one_ellipse(
            as.vector(mus[[sel_sample]][[c]][sel_idx]), 
            sigmas[[sel_sample]][[c]][sel_idx, sel_idx], 
            alpha, c, sel_sample
            )
        ) %>% rbindlist %>% .[, other_qc := dim2_name ]

    return(ellipses_dt)
}

#' Calculates path for an ellipse based on 2D mu and sigma values
#'
#' @param mu mean vector
#' @param sigma covariance matrix
#' @param alpha contour threshold
#' @param c component
#' @param s sample
#' @importFrom magrittr "%>%"
#' @importFrom mixtools ellipse
#' @importFrom data.table data.table setnames ":="
#' @return data.table object
#' @keywords internal
.calc_one_ellipse <- function(mu, sigma, alpha, c, s) {
    ell = ellipse(mu, sigma, alpha=alpha, draw=FALSE)
    dt  = data.table(ell) %>%
        setnames(., names(.), c('qc_y', 'qc_x')) %>%
        .[, component  := c] %>%
        .[, sample_id  := s]
    return(dt)
}

#' Calculates path for an ellipse based on 2D mu and sigma values
#'
#' @param em_obj Output from .fit_one_sampleQC
#' @param s Selected sample
#'
#' @importFrom magrittr "%>%"
#' @importFrom data.table ":=" setnames rbindlist dcast
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 facet_grid theme_bw theme labs
#' @importFrom scales pretty_breaks
#' @return ggplot object
#' @export
plot_outliers_one_sample <- function(em_obj, s) {
    # unpack
    qc_names    = em_obj$qc_names
    dt          = em_obj$outliers_plot[ sample_id == s ]

    # plot
    g   = ggplot() + 
        aes(y=qc_y, x=qc_x, colour=outlier) +
        geom_point(data=dt[outlier == FALSE], size=2, colour='grey' ) +
        geom_point(data=dt[outlier == TRUE],  size=2, colour='red' ) +
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
#' @param fit_obj fit_obj object from function fit_qc_model
#' @param df degrees of freedom for multivariate t
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
#' @return ggplot object
#' @export
plot_alpha_js_likelihoods <- function(fit_obj, df=5) {
    # extract values
    alphas_dt   = data.table(sweep(fit_obj$alpha_j, 2, fit_obj$mu_0, "+")) %>%
        setnames(names(.), fit_obj$qc_names) %>% 
        .[, sample_id := fit_obj$sample_list ] %>%
        .[, sample_id := fct_reorder(sample_id, -fit_obj$alpha_j[, 1]) ] %>%
        melt(
            id='sample_id', 
            variable.name='qc_metric', value.name='qc_value'
            ) %>%
        .[, scaled_qc := scale(qc_value), by='qc_metric' ]

    # max value
    max_val  = ceiling(10*max(abs(alphas_dt$scaled_qc)))/10

    # make plot
    g_alphas    = ggplot(alphas_dt) +
        aes(
            x=qc_metric, y=sample_id, fill=scaled_qc,
            label=sprintf("%0.2f", round(qc_value, digits = 2))
            ) +
        geom_tile() +
        geom_text() +
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
    alpha_mat   = fit_obj$alpha_j
    fit_robt    = cov.trob(alpha_mat, cor=TRUE, nu=df)
    mu          = fit_robt$center
    cov_mat     = fit_robt$cov
    llikes      = dmvt(alpha_mat, mu=mu, sigma=cov_mat, df=df, log=TRUE)
    likes_dt    = data.table(
        loglike     = llikes,
        sample_id   = factor(
            fit_obj$sample_list,
            levels  = levels(alphas_dt$sample_id)
            )
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
#' @param fit_obj fit_obj object from function fit_qc_model
#' @param df degrees of freedom for multivariate t
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table ":=" setnames melt
#' @importFrom forcats fct_reorder
#' @importFrom MASS cov.trob
#' @importFrom mixtools ellipse
#' @importFrom mvnfast dmvt
#' @importFrom ggplot2 ggplot aes_string geom_point geom_path
#' @importFrom ggplot2 scale_fill_distiller theme_bw theme labs
#' @importFrom scales pretty_breaks
#' @importFrom patchwork plot_layout
#' @return ggplot object
#' @export
plot_alpha_js <- function(fit_obj, df=5, qc_idx=c(1,2), pc_idx=c(1,2)) {
    # fit robust multivariate t to alpha values, extract likelihoods
    alpha_mat   = fit_obj$alpha_j
    fit_robt    = cov.trob(alpha_mat, cor=TRUE, nu=df)
    mu          = fit_robt$center
    cov_mat     = fit_robt$cov
    llikes      = dmvt(alpha_mat, mu=mu, sigma=cov_mat, df=df, log=TRUE)

    # get PCs from covariance matrix
    cov_eigen   = eigen(cov_mat)
    pc_scores   = sweep(alpha_mat, 2, mu, "-") %*% cov_eigen$vectors
    std_devs    = sqrt(cov_eigen$values)

    # put together into dt
    qc_names    = fit_obj$qc_names
    n_qcs       = length(qc_names)
    plot_dt     = data.table(pc_scores) %>%
        setnames(names(.), paste0('PC', seq_len(n_qcs))) %>%
        .[, (qc_names) := asplit(alpha_mat, 2) ] %>%
        .[, sample_id := fit_obj$sample_list] %>%
        .[, loglike   := llikes ]

    # make ellipse
    alpha       = 0.05
    ellipse_std = data.table(
        ellipse(mu[qc_idx], cov_mat[qc_idx, qc_idx], alpha=alpha, draw=FALSE)
        ) %>% setnames(names(.), qc_names[qc_idx])
    ellipse_pca = data.table(
        ellipse(rep(0,2), diag(std_devs[pc_idx]^2), alpha=alpha, draw=FALSE)
        ) %>% setnames(names(.), paste0('PC', pc_idx))

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
#' @param fit_obj fit_obj object from function fit_qc_model
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table ":=" setnames melt
#' @importFrom forcats fct_reorder
#' @importFrom MASS cov.trob
#' @importFrom mvnfast dmvt
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_distiller
#' @importFrom ggplot2 theme_bw theme labs element_blank
#' @importFrom scales pretty_breaks
#' @importFrom patchwork plot_layout
#' @return ggplot object
#' @export
plot_beta_ks <- function(fit_obj) {
    # unpack
    if (is.null(fit_obj$p_k)) {
        p_k     = colMeans(fit_obj$p_jk)
    } else {
        p_k     = as.vector(fit_obj$p_k)
    }
    p_jk        = fit_obj$p_jk
    sample_list = fit_obj$sample_list

    # define KL divergence fn
    kl_fn   <- function(Q) 
        sum(p_k * (log2(p_k) - log2(Q)))

    # make data.table with kl vals in
    sample_ord  = sample_list[order(-fit_obj$alpha_j[,1])]
    kl_dt       = data.table(
        sample_id = factor(fit_obj$sample_list, levels=sample_ord),
        kl = apply(fit_obj$p_jk, 1, kl_fn)
        )

    # extract p_jk values
    p_jks_dt    = data.table(fit_obj$p_jk) %>% 
        setnames(names(.), paste0('C', seq_len(fit_obj$K))) %>% 
        .[, sample_id := factor(sample_list, levels=sample_ord) ] %>%
        melt(
            id='sample_id', variable.name='component',
            value.name='cluster_prop'
            )
    diffs_dt    = data.table(sweep(p_jk, 2, p_k, "-")) %>% 
        setnames(names(.), paste0('C', seq_len(fit_obj$K))) %>% 
        .[, sample_id := factor(sample_list, levels=sample_ord) ] %>%
        melt(
            id='sample_id', variable.name='component', 
            value.name='cluster_delta'
            )
    p_jks_dt    = diffs_dt[ p_jks_dt, on=c('sample_id', 'component')]

    # find limits for deltas
    delta_lim   = ceiling(max(abs(diffs_dt$cluster_delta))*10)/10

    # extract outlier proportions
    outliers_dt = data.table(
        sample_id   = factor(sample_list, levels=sample_ord),
        out_prop    = fit_obj$p_jk_outliers[, 'outlier'] / 
            rowSums(fit_obj$p_jk_outliers)
        )

    # make plot
    g_pjks  = ggplot(p_jks_dt) +
        aes(
            x=component, y=sample_id, fill=cluster_delta, 
            label=sprintf("%0.2f", round(cluster_prop, digits = 2))
            ) +
        geom_tile() +
        geom_text() +
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
    g_outs  = ggplot(outliers_dt) +
        aes(
            x="", y=sample_id, fill=out_prop, 
            label=sprintf("%0.3f", round(out_prop, digits = 3))
            ) +
        geom_tile() +
        geom_text() +
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
        scale_fill_distiller(
            palette='PiYG', direction=-1, 
            breaks=pretty_breaks()
            ) +
        expand_limits( fill=0 ) +
        labs(
            x=NULL, y=NULL, 
            fill='KL divergence\nfrom mean\nproportion\n'
            ) +
        theme_bw() + 
        theme(
            legend.position='bottom', legend.direction='vertical', 
            axis.text.y=element_blank()
            )

    # assemble via patchwork
    g = g_pjks + g_outs + g_kl + plot_layout( widths=c(5, 2, 2) )

    return(g)
}

