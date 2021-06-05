#' plot_outliers
#'
#' Calculates path for an ellipse based on 2D mu and sigma values
#'
#' @param qc_obj Output from function fit_sampleqc
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
plot_outliers <- function(qc_obj, sel_sample, outliers_dt=NULL) {
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
