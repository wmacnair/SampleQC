# Miscellaneous useful functions 

# excellent set of colours copied from the CATALYST package, here:
# https://bioconductor.org/packages/release/bioc/html/CATALYST.html
.CLUSTER_COLS = c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", 
    "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", 
    "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4", 
    "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", 
    "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
    )

#' Extracts number of groups
#' 
#' @param qc_obj Output from \code{calc_pairwise_mmds}
#' 
#' @importFrom assertthat assert_that
#' 
#' @return number of groups identified by \code{calc_pairwise_mmds}
#' @export
get_n_groups <- function(qc_obj) {
    # check ok
    assert_that( .check_is_qc_obj(qc_obj) %in% c('mmd', 'fit'),
        msg='not a SampleQC object')

    # extract n_groups
    n_groups    = metadata(qc_obj)$n_groups

    return(n_groups)
}

#' Checks whether a given sce object is a SampleQC object
#' 
#' @param qc_obj SampleQC object
#' 
#' @importFrom assertthat assert_that
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assays
#' @return string saying progress
#' 
#' @keywords internal
.check_is_qc_obj <- function(qc_obj) {
    # check is an sce
    assert_that(
        is(qc_obj, "SingleCellExperiment"),
        msg='not SampleQC object: must be SingleCellExperiment'
        )

    # check that correct assay names are present
    assert_that(
    all( c('mmd', 'mmd_adj') %in% names(assays(qc_obj)) ),
    msg='not SampleQC object: "mmd" and "mmd_adj", must be present in assays'
        )
    # check that base colData names are present
    cd_names    = names(colData(qc_obj))
    base_names  = c('sample_id', 'cell_id', 'qc_metrics')
    assert_that(
        all( base_names %in% cd_names ),
        msg=paste0('not SampleQC object: "', 
            paste(setdiff(base_names, cd_names), collapse='", "'), 
            '" missing from colData')
        )
    # check that base annotations are present
    base_disc   = c('N_cat', 'mito_cat')
    names_disc  = names(colData(qc_obj)$annot_disc[[1]])
    assert_that(
        all( base_disc %in% names_disc ),
        msg=paste0('not SampleQC object: "', 
            paste(setdiff(base_disc, names_disc), collapse='", "'), 
            '" missing from annotations')
        )    
    base_cont   = c('log_N', 'med_mito', 'med_counts')
    names_cont  = names(colData(qc_obj)$annot_cont[[1]])
    assert_that(
        all( base_cont %in% names_cont ),
        msg=paste0('not SampleQC object: "', 
            paste(setdiff(base_cont, names_cont), collapse='", "'), 
            '" missing from annotations')
        )    
    # check that MMD colData names are present
    mmd_names   = c('group_id')
    assert_that(
        all( mmd_names %in% cd_names ),
        msg=paste0('not SampleQC object: "', 
            paste(setdiff(mmd_names, cd_names), collapse='", "'), 
            '" missing from colData')
        )
    # check that dims of colData match assay dims
    assert_that(
        ncol(qc_obj)==nrow(qc_obj),
        msg=paste0('not SampleQC object: assays should be square')
        )
    # check that nested colData entries have consistent sizes
    assert_that(
        ncol(qc_obj)==nrow(qc_obj),
        msg=paste0('not SampleQC object: assays should be square')
        )
    # check metadata
    meta_names  = names(metadata(qc_obj))
    meta_needed = c("qc_names", "D", "n_groups", "group_list", "mmd_params", 
        "annots")
    assert_that(
        all( meta_needed %in% meta_names ),
        msg=paste0('metadata for qc_obj incorrect')
        )
    assert_that(
        all(colnames(colData(qc_obj)$annot_disc[[1]]) == 
            metadata(qc_obj)$annots$disc),
        msg=paste0('names of annot_disc in colData(qc_obj) should match 
            metadata(qc_obj)$annots$disc')
        )
    assert_that(
        all(colnames(colData(qc_obj)$annot_cont[[1]]) == 
            metadata(qc_obj)$annots$cont),
        msg=paste0('names of annot_cont in colData(qc_obj) should match 
            metadata(qc_obj)$annots$cont')
        )

    # if ok so far, definitely return mmd
    obj_type    = 'mmd'

    # are fit colData names also present?
    fit_names   = c('K', 'mu_0', 'alpha_j', 
        'beta_k', 'sigma_k', 'p_jk', 'z', 'outlier')
    if ( !all(fit_names %in% cd_names) )
        return(obj_type)

    # are fit meta things present?
    fit_metas   = c('fit_list', 'fit_params')
    if ( !all(fit_metas %in% meta_names) )
        return(obj_type)

    # if we got through to here, assume that we fit the thing
    obj_type    = 'fit'
    
    return(obj_type)
}

#' Generate fake qc_df object for testing
#'
#' @importFrom data.table data.table rbindlist ":="
#' @return data.frame
#' @keywords internal
.make_toy_qc_df <- function() {
    # generate sample means
    J               = 20
    N_per_sample    = as.integer(exp(rnorm(J, log(200), 1)))
    mu_log_counts   = rnorm(J, log(4000), 0.1)
    mu_log_feats    = mu_log_counts - 0.5 + rnorm(J, 0, 0.1)
    mu_logit_mito   = rnorm(J, -3, 1)

    # generate cell data
    df_list = lapply(seq_len(J), function(j) {
        # unpack
        N           = N_per_sample[[j]]
        mu_counts   = mu_log_counts[[j]]
        mu_feats    = mu_log_feats[[j]]
        mu_mito     = mu_logit_mito[[j]]

        # generate random data
        df  = data.frame(
            sample_id   = sprintf('sample%02d', j),
            annot_1     = sprintf('annot%02d', j),
            log_counts  = rnorm(N, mu_counts, 1)
            )
        df$log_feats    = (df$log_counts - mu_counts) + mu_feats + rnorm(N, 0, 0.01)
        df$mito_prop    = plogis(rnorm(N, mu_mito, 0.1))

        return(df)
    })
    qc_df   = do.call(rbind, df_list)

    # add cell_id
    qc_df$cell_id = sprintf('cell%04d', seq_len(nrow(qc_df)))
    n_cols  = ncol(qc_df)
    qc_df   = qc_df[, c(n_cols, seq_len((n_cols-1)))]

    return(qc_df)
}

#' Generate fake sce object for testing
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @return sce
#' @keywords internal
.make_toy_sce <- function() {
    # generate sample_ids
    J               = 10
    N_per_sample    = as.integer(exp(rnorm(J, log(200), 1)))
    sample_ids      = rep(sprintf('sample%02d', seq_len(J)), times=N_per_sample)
    annot_1         = rep(sprintf('annot%02d', seq_len(J)), times=N_per_sample)

    # generate genes
    n_genes         = 100
    gene_names      = sprintf('gene%03d', seq_len(n_genes))
    n_mt            = 13
    mt_genes        = sprintf('mt-%02d', seq_len(n_mt))
    gene_names[seq_len(n_mt)]  = mt_genes

    # make count matrix
    n_cells         = length(sample_ids)
    counts_mat      = matrix(rpois(n_cells*n_genes, lambda=10), ncol=n_cells)

    # make column data
    cols_df         = data.frame(
        cell_id     = sprintf('cell%04d', seq_len(n_cells)),
        sample_id   = factor(sample_ids),
        annot_1     = factor(annot_1)
        )

    # make sce object
    sce         = SingleCellExperiment(
        list(counts=counts_mat),
        colData = cols_df
        )
    colnames(sce)   = cols_df$cell_id
    rownames(sce)   = gene_names

    return(sce)
}
