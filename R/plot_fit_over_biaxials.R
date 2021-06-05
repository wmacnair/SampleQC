#' plot_fit_over_biaxials
#'
#' Plots inferred sample-level statistics over QC biaxials
#'
#' @param qc_obj Output from function fit_sampleqc
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
plot_fit_over_biaxials <- function(qc_obj, sel_sample, 
    qc_names=NULL, alpha_cut=0.01) {
    # can we do this?
    assert_that(
        .check_is_qc_obj(qc_obj) == 'fit',
        msg='SampleQC model must be fit (via `fit_sampleqc`) before calling 
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
#' @param qc_obj Output from function fit_sampleqc
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
