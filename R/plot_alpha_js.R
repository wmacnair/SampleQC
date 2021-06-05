#' plot_alpha_js
#'
#' Plots inferred sample-level statistics, with empirical log-likelihood
#'
#' @param qc_obj Output from function fit_sampleqc
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
