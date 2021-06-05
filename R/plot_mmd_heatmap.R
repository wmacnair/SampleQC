#' plot_mmd_heatmap
#'
#' Plots heatmap of sample-sample distances (as measured by MMD)
#' 
#' @param qc_obj Outputs from function calc_pairwise_mmds
#' 
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table ":=" melt as.data.table
#' @importFrom scales pretty_breaks
#' @importFrom ggplot2 ggplot aes geom_tile
#' @importFrom ggplot2 scale_fill_distiller expand_limits
#' @importFrom ggplot2 theme_bw theme element_blank element_text labs
#'
#' @return ggplot object
#' @export
plot_mmd_heatmap <- function(qc_obj) {
    # some checks
    .check_is_qc_obj(qc_obj)

    # get stuff we need
    mmd_mat     = assay(qc_obj, 'mmd')
    n_times     = metadata(qc_obj)$mmd_params$n_times
    
    # do hierarchical clustering, use this to order the samples
    hclust_obj  = hclust(as.dist(mmd_mat), method='complete')
    sample_ord  = hclust_obj$labels[hclust_obj$order]

    # make mmd_dt 
    mmd_dt      = as.data.table(mmd_mat, keep.rownames = 'sample_i') %>% 
        melt(id = 'sample_i', variable.name = 'sample_j', 
            value.name='mmd_mean') %>%
        .[, sample_i := factor(sample_i, levels=sample_ord) ] %>%
        .[, sample_j := factor(sample_j, levels=sample_ord) ] %>%
        .[ sample_i != sample_j ]

    # plot
    g = ggplot(mmd_dt) +
        aes( x=sample_i, y=sample_j, fill=mmd_mean ) +
        geom_tile() +
        scale_fill_distiller(
            palette='PiYG', direction=-1, 
            breaks=pretty_breaks()
            ) + 
        expand_limits( fill=0 ) +
        theme_bw() + theme( 
            aspect.ratio    = 1,
            axis.text.x     = element_blank(), 
            axis.text.y     = element_text(size=6) ) +
        labs( fill=sprintf('mean MMD\n(mean used\n%d subsamples)', n_times))

    return(g)
}
