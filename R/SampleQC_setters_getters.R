# SampleQC: robust multivariate, multi-celltype, multi-sample quality control 
# for single cell RNA-seq
#' devtools::load_all('~/work/packages/BayesQC')
#' devtools::document('~/work/packages/BayesQC')

#' SampleQC_setters_getters.R
#' Utility functions to define and extract values from a \code{qc_obj}.

#' Extracts number of groups
#' 
#' @param qc_obj Output from \code{calculate_sample_to_sample_MMD}
#' 
#' @importFrom assertthat assert_that
#' 
#' @return number of groups identified by \code{calculate_sample_to_sample_MMD}
#' @export
get_n_groups <- function(qc_obj) {
    # check ok
    assert_that( .check_is_qc_obj(qc_obj) %in% c('mmd', 'fit'),
        msg='not a SampleQC object')

    # extract n_groups
    n_groups    = metadata(qc_obj)$n_groups

    return(n_groups)
}

#' Extracts list of outliers
#' 
#' @param qc_obj Output from fit_sampleQC
#' @param exc_groups List of sample groups to exclude
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom data.table rbindlist
#' 
#' @return \code{data.table} containing cell_id, sample_id and whether outlier 
#' or not
#' @export
get_outliers <- function(qc_obj, exc_groups=NULL, exc_clusters=NULL) {
    # check ok
    assert_that( .check_is_qc_obj(qc_obj) == 'fit',
        msg='fit_sampleQC must be run before outliers can be returned')
    assert_that( all(exc_groups %in% metadata(qc_obj)$group_list),
        msg='invalid sample groups specified for exc_groups')
    assert_that( all(names(exc_clusters) %in% metadata(qc_obj)$group_list),
        msg='invalid sample groups specified for exc_clusters')

    # extract cell outliers
    outliers_dt     = colData(qc_obj)$outlier %>% 
        lapply(function(dt) dt[, .(sample_id, cell_id, out_cell=outlier)]) %>%
        rbindlist

    # exclude specified samples
    exc_samples     = qc_obj$sample_id[qc_obj$group_id %in% exc_groups]
    outliers_dt[, out_sample    := sample_id %in% exc_samples]

    # exclude specified clusters
    outliers_dt[, out_cluster := FALSE]
    for (ii in seq_along(exc_clusters)) {
        # which group?
        g_ii        = names(exc_clusters)[[ii]]
        z_ii        = exc_clusters[[ii]]

        # which cells to exclude?
        z_vec       = do.call(c, qc_obj$z[qc_obj$group_id == g_ii])
        cell_ids    = do.call(c, qc_obj$cell_id[qc_obj$group_id == g_ii])
        out_cells   = cell_ids[ z_vec %in% z_ii ]
        outliers_dt[ cell_id %in% out_cells , out_cluster := TRUE  ]
    }

    # combine
    outliers_dt[, outlier   := out_cell | out_sample | out_cluster ]

    return(outliers_dt)
}
