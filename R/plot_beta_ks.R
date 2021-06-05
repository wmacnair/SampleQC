#' plot_beta_ks
#'
#' Plots estimates of cluster splits by sample, with KL divergence from mean
#'
#' @param qc_obj Output from function fit_sampleqc
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
