# SampleQC: robust multivariate, multi-celltype, multi-sample quality control 
# for single cell RNA-seq
# devtools::load_all('~/work/packages/SampleQC')
# devtools::document('~/work/packages/SampleQC')
# debug(calculate_sample_to_sample_MMDs); qc_obj = calculate_sample_to_sample_MMDs(qc_dt, qc_names=qc_names, annots_disc=annots_disc, n_cores=16)

# SampleQC_samples.R
# Code to estimate multivariate dissimilarities between samples, and use
# this to identify groups of clusters whose cell QC should be fit 
# together.

# excellent set of colours copied from the CATALYST package, here:
# https://bioconductor.org/packages/release/bioc/html/CATALYST.html
.CLUSTER_COLS = c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", 
    "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", 
    "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4", 
    "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", 
    "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
    )

#' Calculates MMD distances between all samples in qc_dt object, then and 
#' embeds them with some annotations.
#' 
#' @param qc_dt Data.table of QC metrics for all cells and samples
#' @param qc_names List of metrics to actually use for calculating 
#' sample-to-sample distances
#' @param sigma Scale for MMD kernel (defaults to length of qc_names)
#' @param n_cores How many cores to parallelize over
#' @param subsample Should we downsample the number of cells per sample to 
#' this number
#' @param n_times How many times do we sample MMD between each pair of samples?
#' (if subsampled, MMD is a stochastic value)
#' @param centre_samples,scale_samples Should we centre or scale the values 
#' within each sample before calculating MMD?
#' 
#' @section Details:
#' [Maximum mean discrepancy (MMD) ]
#' [Why to centre and not scale]
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table setkey
#' @importFrom igraph graph_from_edgelist graph_from_adjacency_matrix
#' @importFrom igraph as_adjacency_matrix
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim<-
#' @return slightly funky `sce` object
#' @export
calculate_sample_to_sample_MMDs <- function(qc_dt, qc_names=c('log_counts', 
    'log_feats', 'logit_mito'), annots_disc=NULL, annots_cont=NULL, 
    n_cores=4, sigma, subsample=200, n_times=10,
    centre_samples=TRUE, scale_samples=FALSE) {
    # check some inputs
    if( missing(sigma) )
        sigma   = length(qc_names)
    assert_that( length(sigma) == 1 )
    assert_that( sigma > 0 )

    message('sample-level parts of SampleQC:')
    # get list of samples to compare
    sample_list = sort(unique(qc_dt$sample_id))
    n_samples   = length(sample_list)

    # split qc metric values into one matrix per sample
    mat_list    = .calc_mat_list(qc_dt, qc_names, sample_list, 
        centre_samples, scale_samples)

    # calculate mmds
    mmd_mat     = .calc_mmd_mat(sample_list, mat_list, 
        n_times, subsample, sigma, n_cores)

    # cluster on basis of mmd_mat
    mmd_graph   = .make_mmd_graph(mmd_mat, n_nhbrs=5)
    group_ids   = .cluster_mmd_graph(mmd_graph)

    # do embeddings
    mds_mat     = .embed_mds(mmd_mat)
    n_nhbrs     = 5
    umap_mat    = .embed_umap(mmd_mat, n_nhbrs)

    # define automatic names
    message('  adding annotation variables')
    auto_disc   = c('N_cat', 'mito_cat', 'counts_cat')
    auto_cont   = c('log_N', 'med_mito', 'med_counts')
    annots_disc = .process_annot_vars(qc_dt, annots_disc, auto_disc)
    annots_cont = .process_annot_vars(qc_dt, annots_cont, auto_cont)
    disc_dt     = qc_dt[, c('sample_id', annots_disc), with=FALSE] %>% 
        unique %>% setkey('sample_id') %>% .[ sample_list ] %>%
        .[, group_id := group_ids ] %>% setcolorder('group_id')
    annots_disc = c('group_id', annots_disc)
    cont_dt     = qc_dt[, c('sample_id', annots_cont), with=FALSE] %>% 
        unique %>% setkey('sample_id') %>% .[ sample_list ]

    # assemble outputs
    message('  constructing SampleQC object')
    assay_list  = list(
        mmd     = mmd_mat,
        mmd_adj = as_adjacency_matrix(mmd_graph)
        )

    # assemble values for samples
    col_DF      = DataFrame(
        sample_id   = sample_list,
        cell_id     = I(split(qc_dt$cell_id, qc_dt$sample_id)),
        qc_metrics  = I(split(qc_dt[, qc_names, with=FALSE], qc_dt$sample_id)),
        group_id    = group_ids
        )

    # add annotations
    assert_that( all(disc_dt$sample_id == col_DF$sample_id) )
    assert_that( all(cont_dt$sample_id == col_DF$sample_id) )
    col_DF$annot_disc   = I(split(disc_dt[, -'sample_id'], disc_dt$sample_id))
    col_DF$annot_cont   = I(split(cont_dt[, -'sample_id'], cont_dt$sample_id))

    # record MMD parameters
    mmd_params  = list(
        sigma           = sigma, 
        subsample       = subsample, 
        n_times         = n_times, 
        centre_samples  = centre_samples, 
        scale_samples   = scale_samples
        )

    # put together into one object
    n_groups    = length(unique(group_ids))
    qc_obj      = SingleCellExperiment(
        assays      = assay_list,
        colData     = col_DF,
        metadata    = list(
            qc_names    = qc_names,
            D           = length(qc_names),
            n_groups    = n_groups,
            group_list  = paste0('SG', seq.int(n_groups)),
            mmd_params  = mmd_params,
            annots      = list(
                disc = annots_disc,
                cont = annots_cont
                )
            )
        )

    # add dim-reductions
    reducedDim(qc_obj, 'MDS')   <- mds_mat
    reducedDim(qc_obj, 'UMAP')  <- umap_mat

    # check it worked ok
    .check_is_qc_obj(qc_obj)

    return(qc_obj)
}

#' Splits qc_dt into list of matrices, one for each sample
#' 
#' @param qc_dt Data.table of QC metrics for all cells and samples
#' @param qc_names List of metrics to actually use for calculating 
#' @param sample_list List of samples to iterate through
#' @param centre_samples,scale_samples Should we centre or scale the values 
#' within each sample before calculating MMD?
#'
#' @importFrom magrittr "%>%"
#' 
#' @return List of matrices
#' 
#' @keywords internal
.calc_mat_list <- function(qc_dt, qc_names, sample_list, centre_samples, scale_samples) {
    # scale overall matrix
    qc_mat_all  = apply(as.matrix(qc_dt[ , qc_names, with=FALSE ]), 2, scale)

    # make list of individual matrices, scaled if necessary
    sample_ids  = qc_dt$sample_id
    mat_list    = lapply(sample_list, function(s) {
        qc_mat  = qc_mat_all[ sample_ids == s, ]
        qc_mat  = apply(qc_mat, 2, scale, 
            center=centre_samples, scale=scale_samples)
        return(qc_mat)
    }) %>% setNames(sample_list)

    return(mat_list)
}

#' Calculates mean MMD distances between all pairs of samples
#' 
#' @param sample_list List of samples to iterate through
#' @param mat_list List of matrices of QC metrics
#' @param n_times How many times do we sample MMD between each pair of samples?
#' (if subsampled, MMD is a stochastic value)
#' @param subsample Should we downsample the number of cells per sample to 
#' this number
#' @param sigma Scale for MMD kernel (defaults to length of qc_names)
#' @param n_cores How many cores to parallelize over
#'
#' @importFrom parallel mclapply
#' @importFrom data.table data.table dcast
#' @importFrom magrittr "%>%" set_rownames
#' 
#' @return matrix of MMD values
#' @keywords internal
.calc_mmd_mat <- function(sample_list, mat_list, n_times, subsample, sigma, n_cores) {
    # define combns to do
    ij_grid     = expand.grid(
        sample_i=factor(sample_list), 
        sample_j=factor(sample_list)
        )
    ij_grid     = ij_grid[ as.integer(ij_grid$sample_i) < 
        as.integer(ij_grid$sample_j), ]

    # iterate through them
    n_combns    = nrow(ij_grid)
    message('  calculating ', n_combns, 
        ' sample-sample MMDs (. = 20, one row ~= 1000):\n  ', appendLF=FALSE)
    mmd_means   = mclapply(
        seq_len(n_combns), function(r) {
            if( (r%%20) == 0 )
                message(".", appendLF=FALSE)
            if( (r%%1000) == 0 )
                message("\n  ", appendLF=FALSE)
            # prep
            sample_i    = as.character(ij_grid[r, 'sample_i'])
            sample_j    = as.character(ij_grid[r, 'sample_j'])
            mat_i       = mat_list[[sample_i]]
            mat_j       = mat_list[[sample_j]]

            # calc MMD
            mmd_vec     = vapply(
                seq_len(n_times), 
                function(z) return(.dist_fn(mat_i, mat_j, subsample, sigma)), 
                numeric(1))
            return(mean(mmd_vec))
        }, mc.cores=n_cores)
    message("\n", appendLF=FALSE)

    # put into data.table, add reverse versions
    mmd_dt      = data.table(ij_grid, mmd_mean=unlist(mmd_means))
    mmd_dt      = rbind(
        mmd_dt, 
        mmd_dt[, list(sample_i=sample_j, sample_j=sample_i, mmd_mean)]
        )

    # make matrix
    mmd_mat     = dcast(
        mmd_dt, sample_i ~ sample_j, 
        value.var='mmd_mean', fill=0) %>% 
        .[, -'sample_i', with=FALSE] %>% 
        as.matrix %>%
        set_rownames(., colnames(.))

    return(mmd_mat)
}

#' Wrapper for cpp MMD distance function
#' 
#' @param m_i matrix i
#' @param m_j matrix j
#' @param subsample Number to subsample down to
#' @param sigma Number to subsample down to
#'
#' @importFrom assertthat assert_that
#' 
#' @return MMD between subsamples of input matrices
#' @keywords internal
.dist_fn <- function(m_i, m_j, subsample, sigma) {
    # check matrices are the same size
    assert_that( ncol(m_i) == ncol(m_j) )

    # get values for sample i, subsample if necessary
    n_i     = nrow(m_i)
    if (n_i > subsample) {
        m_i   = m_i[ sample(n_i, subsample), ]
    }
    # get values for sample j, subsample if necessary
    n_j     = nrow(m_j)
    if (n_j > subsample) {
        m_j   = m_j[ sample(n_j, subsample), ]
    }
    return(abs(mmd_fn(m_i, m_j, sigma)))
}

#' Build igraph object from MMD matrix
#' 
#' @param mmd_mat
#' @param n_nhbrs How many neighbours in the NN graph
#'
#' @importFrom magrittr "%>%"
#' @importFrom igraph graph_from_edgelist
#' 
#' @return igraph object
#' @keywords internal
.make_mmd_graph <- function(mmd_mat, n_nhbrs=5) {
    # turn into nearest neighbours graph
    edge_list   = lapply(seq_len(nrow(mmd_mat)), 
        function(i) {
            row     = mmd_mat[i, ]
            nbrs    = order(row)[2:(n_nhbrs+1)]
            return(data.frame(i=i, j=nbrs))
        }) %>% rbindlist
    graph_obj   = graph_from_edgelist(as.matrix(edge_list), directed=FALSE)
    # graph_obj   = graph_from_adjacency_matrix(abs(mmd_mat), mode='undirected', weighted=TRUE)
    return(graph_obj)
}

#' Cluster igraph object using Louvain clustering
#' 
#' @param mmd_graph
#'
#' @importFrom igraph cluster_louvain
#' @importFrom magrittr "%>%"
#' 
#' @return vector of group_ids for all samples
#' @keywords internal
.cluster_mmd_graph <- function(mmd_graph) {
    message('  clustering samples using MMD values')
    # do clustering
    cluster_obj = cluster_louvain(mmd_graph)
    group_ids   = cluster_obj$membership

    # rename in order
    group_ns    = table(group_ids)
    new_ids     = seq.int(length(unique(group_ids))) %>%
        setNames(names(sort(-group_ns)))
    group_ids   = new_ids[ as.character(group_ids) ] %>% unname
    group_ids   = paste0('SG', group_ids) %>% factor

    return(group_ids)
}

#' MDS embedding of MMD distances
#' 
#' @param mmd_mat matrix of MMD distances
#' 
#' @importFrom scales rescale
#' 
#' @return matrix of embedding
#' @export
.embed_mds <- function(mmd_mat) {
    message('  calculating MDS embedding')
    # embed with multidimensional scaling
    k           = 2
    mds_mat     = cmdscale(mmd_mat, k)
    mds_mat     = apply(mds_mat, 2, rescale, to=c(0.05, 0.95))
    colnames(mds_mat)   = paste0('MDS', 1:k)

    return(mds_mat)
}

#' Calculates non-linear embeddings of MMD distances calculated by 
#' calculate_sample_to_sample_MMDs
#' 
#' @param mmd_mat MMD distances
#' @param n_nhbrs How many neighbours for UMAP nearest neighbours graph
#' 
#' @importFrom uwot umap
#' @importFrom scales rescale
#' @importFrom magrittr "%>%"
#'
#' @return matrix of embedding
#' @keyword internal
.embed_umap <- function(mmd_mat, n_nhbrs) {
    message('  calculating UMAP embedding')

    # convert mmd mat to nearest neighbours
    # nn_method   = list(
    #     idx     = apply(
    #         mmd_mat, 1, function(row) order(row)[2:(n_nhbrs+1)]
    #         ) %>% t,
    #     dist    = apply(
    #         mmd_mat, 1, function(row) sort(row)[2:(n_nhbrs+1)]
    #         ) %>% t
    #     )
    # umap_mat    = umap(NULL, nn_method=nn_method, min_dist=0.001)
    umap_mat    = umap(as.dist(mmd_mat), min_dist=0.001)

    # tidy up
    umap_mat    = apply(umap_mat, 2, rescale, to=c(0.05, 0.95))
    colnames(umap_mat)   = paste0('UMAP', 1:ncol(umap_mat))

    return(umap_mat)
}

#' Checks annotation variables
#' 
#' @param qc_dt
#' @param specified Requested annotation variables
#' @param auto Automatically generated annotation variables
#' 
#' @return character vector
#' @keyword internal
.process_annot_vars <- function(qc_dt, specified, auto) {
    # warn about missing annotations
    dt_names        = colnames(qc_dt)
    missing_vals    = setdiff(specified, dt_names)
    if (length(missing_vals) > 0)
        warning("These variables were requested as annotations, but 
            aren't present in qc_dt:\n", paste(missing_vals, collapse=', '))

    # add automatic names
    annot_vars      = c(specified, setdiff(auto, specified))
    annot_vars      = intersect(annot_vars, dt_names)

    return(annot_vars)
}

#' Plots non-linear embeddings of MMD distances calculated by 
#' embed_sample_to_sample_MMDs
#' 
#' @param qc_obj SampleQC object
#' @param var_type Is this annotation discrete or continuous?
#' @param sel_embed Which type of embedding to use
#' 
#' @importFrom assertthat assert_that
#' @importFrom S4Vectors metadata
#' @return NULL
#' @export
plot_embeddings <- function(qc_obj, var_type=c('discrete', 'continuous'), 
    sel_embed=c('MDS', 'UMAP')) {
    # check inputs
    .check_is_qc_obj(qc_obj)
    var_type    = match.arg(var_type)
    sel_embed   = match.arg(sel_embed)

    if (var_type == 'discrete') {
        annot_vars  = metadata(qc_obj)$annots$disc
    } else {
        annot_vars  = metadata(qc_obj)$annots$cont
    }
    for (v in annot_vars) {
        cat('### ', v, '\n')
        g = .plot_one_embedding(qc_obj, v, var_type, sel_embed)
        print(g)
        cat('\n\n')
    }
}

#' Plots non-linear embeddings of MMD distances calculated by 
#' embed_sample_to_sample_MMDs
#' 
#' @param qc_obj SampleQC object
#' @param v Which annotation variable to plot
#' @param var_type Is this annotation discrete or continuous?
#' @param sel_embed Which type of embedding to use
#' 
#' @importFrom assertthat assert_that
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales pretty_breaks
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 ggplot aes_string geom_point geom_text
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_fill_brewer
#' @importFrom ggplot2 scale_fill_distiller scale_fill_manual
#' @importFrom ggplot2 theme theme_bw labs
#'
#' @return ggplot object
#' @keywords internal
.plot_one_embedding <- function(qc_obj, v, var_type, sel_embed) {
    # unpack
    .check_is_qc_obj(qc_obj)    

    # checks
    assert_that( var_type %in% c('discrete', 'continuous') )
    assert_that( sel_embed %in% c('MDS', 'UMAP') )

    # make plot_dt
    plot_dt     = data.table(
        sample_id   = colData(qc_obj)$sample_id,
        reducedDim(qc_obj, sel_embed)
        )
    # add annotations
    if ( var_type=='discrete' ) {
        tmp_dt      = rbindlist(colData(qc_obj)$annot_disc)
    } else {
        tmp_dt      = rbindlist(colData(qc_obj)$annot_cont)
    }
    plot_dt     = cbind(plot_dt, tmp_dt)
    assert_that( v %in% names(plot_dt) )

    # restrict to chosen embedding
    if (var_type == 'discrete') {
        n_cols      = length(unique(plot_dt[[v]]))
        if (v %in% c('mito_cat', 'N_cat', 'counts_cat', 'totals_cat')) {
            if (n_cols < 3)
                col_vals    = rev(brewer.pal(n_cols, 'PiYG'))[c(1,3)]
            else
                col_vals    = rev(brewer.pal(n_cols, 'PiYG'))
        } else {
            col_vals    = .CLUSTER_COLS
        }
        # for N_cat, reverse order
        if ( v %in% c('N_cat', 'counts_cat', 'totals_cat') )
            col_vals    = rev(col_vals)

        g = ggplot(plot_dt) + 
            aes_string(
                x=paste0(sel_embed, '1'), y=paste0(sel_embed, '2'), 
                label='sample_id', fill=paste0("factor(", v, ")")
                ) + 
            geom_text(
                hjust=0.5, vjust=1.5, colour='black', 
                size=2, alpha=0.8
                ) + 
            geom_point( size=3, alpha=0.8, shape=21, colour='black' ) +
            scale_x_continuous( breaks=pretty_breaks() ) + 
            scale_y_continuous( breaks=pretty_breaks() )
        if (n_cols > 2) {
            g = g + scale_fill_manual( values=col_vals[seq_len(n_cols)] )
        } else {
            g = g + scale_fill_brewer( palette='Set1' )
        }
        g = g + 
            theme_bw() + theme( aspect.ratio=1, axis.text=element_blank() ) + 
            labs( fill=v )

    } else if (var_type == 'continuous') {
        # for some, reverse order
        direction   = -1
        if ( v %in% c('log_N', 'med_counts') )
            direction   = 1

        g = ggplot(plot_dt) + 
            aes_string( x=paste0(sel_embed, '1'), y=paste0(sel_embed, '2'), 
                label='sample_id', fill=v ) + 
            geom_text(
                hjust=0.5, vjust=1.5, 
                colour='black', size=2, alpha=0.8
                ) +
            geom_point( size=3, alpha=0.8, shape=21, colour='black' ) +
            scale_fill_distiller( palette='RdBu', breaks=pretty_breaks(), 
                direction=direction ) + 
            scale_x_continuous( breaks=pretty_breaks() ) + 
            scale_y_continuous( breaks=pretty_breaks() ) +
            theme_bw() + theme( aspect.ratio=1, axis.text=element_blank() )
    }

    # do nice label
    label_list  = list(
        group_id    = 'sample\ngroup', 
        N_cat       = 'number\nof cells', 
        mito_cat    = 'median\nmito propn.', 
        counts_cat  = 'median\nlog10(counts)',
        log_N       = 'log10(number\nof cells)', 
        med_mito    = 'median\nmito propn.', 
        med_counts  = 'median\nlog10(counts)'
    )
    g = g + labs( fill=label_list[[v]] )

    return(g)
}

#' Plots heatmap of sample-sample distances (as measured by MMD)
#' 
#' @param qc_obj Outputs from function calculate_sample_to_sample_MMDs
#' 
#' @importFrom assertthat assert_that
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table ":="
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
    mmd_dt      = reshape2:::melt.array(mmd_mat, 
        varnames=c('sample_i', 'sample_j'), value.name='mmd_mean') %>%
        data.table %>%
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
