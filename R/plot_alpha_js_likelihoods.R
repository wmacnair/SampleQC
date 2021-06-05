#' plot_alpha_js_likelihoods
#'
#' Plots inferred sample-level statistics, with empirical log-likelihood
#'
#' @param qc_obj Output from function fit_sampleqc
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
