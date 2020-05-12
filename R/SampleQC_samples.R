# SampleQC: robust multivariate, multi-celltype, multi-sample quality control for single cell RNA-seq
# devtools::load_all('~/work/packages/BayesQC')
# devtools::document('~/work/packages/BayesQC')

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

#' Wrapper for cpp MMD distance function
#' 
#' @param m_i matrix i
#' @param m_j matrix j
#' @param subsample Number to subsample down to
#' @param sigma Number to subsample down to
#'
#' @importFrom assertthat assert_that
#' @return MMD between subsamples of input matrices
#' @keywords internal
.dist_fn <- function(m_i, m_j, subsample, sigma) {
    # check matrices are the same size
    assert_that( ncol(m_i)==ncol(m_j) )

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

#' Calculates MMD distances between all samples in qc_dt object
#' 
#' @param qc_dt Data.table of QC metrics for all cells and samples
#' @param qc_names List of metrics to actually use for calculating sample-to-sample distances
#' @param sigma Scale for MMD kernel (defaults to length of qc_names)
#' @param n_cores How many cores to parallelize over
#' @param subsample Should we downsample the number of cells per sample to this number
#' @param n_times How many times do we sample MMD between each pair of samples? (if subsampled, MMD is a stochastic value)
#' @param centre_samples,scale_samples Should we centre or scale the values within each sample before calculating MMD?
#' 
#' @section Details
#' [Maximum mean discrepancy (MMD) ]
#' [Why to centre and not scale]
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom parallel mclapply
#' @importFrom data.table data.table
#' @importFrom data.table dcast
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_louvain
#' @return list, containing MMD values and sample clusters based on MMD values
#' @export
calculate_sample_to_sample_MMDs <- function(qc_dt, qc_names, sigma, n_cores=16, subsample=1000, n_times=10, centre_samples=TRUE, scale_samples=FALSE) {
    # check some inputs
    if( missing(sigma) )
        sigma   = length(qc_names)
    assert_that( length(sigma)==1 )
    assert_that( sigma>0 )

    # get list of samples to compare
    sample_list = sort(unique(qc_dt$sample_id))
    n_samples   = length(sample_list)

    # scale overall matrix
    qc_mat_all  = apply(as.matrix(qc_dt[ , qc_names, with=FALSE ]), 2, scale)

    # make list of individual matrices, scaled if necessary
    sample_ids  = qc_dt$sample_id
    mat_list    = lapply(sample_list, function(s) {
        qc_mat  = qc_mat_all[ sample_ids==s, ]
        qc_mat  = apply(qc_mat, 2, scale, center=centre_samples, scale=scale_samples)
        return(qc_mat)
    }) %>% setNames(sample_list)

    # define combns to do
    ij_grid     = expand.grid(sample_i=factor(sample_list), sample_j=factor(sample_list))
    ij_grid     = ij_grid[ as.integer(ij_grid$sample_i) < as.integer(ij_grid$sample_j), ]

    # iterate through them
    n_combns    = nrow(ij_grid)
    message(sprintf('calculating %d sample-sample MMDs (. = 20, one row ~= 1000):', n_combns, n_combns/20))
    mmd_means   = mclapply(
        1:n_combns, function(r) {
            if( (r%%20)==0 )
                message(".", appendLF=FALSE)
            if( (r%%1000)==0 )
                message("\n", appendLF=FALSE)
            # prep
            sample_i    = as.character(ij_grid[r, 'sample_i'])
            sample_j    = as.character(ij_grid[r, 'sample_j'])
            mat_i       = mat_list[[sample_i]]
            mat_j       = mat_list[[sample_j]]

            # calc MMD
            mmd_vec     = sapply(1:n_times, function(dummy) return(.dist_fn(mat_i, mat_j, subsample, sigma)) )
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
    mmd_mat     = dcast(mmd_dt, sample_i ~ sample_j, value.var='mmd_mean', fill=0) %>% 
        .[, -'sample_i', with=FALSE] %>% 
        as.matrix %>%
        set_rownames(., colnames(.))

    # turn into nearest neighbours graph
    message('   clustering samples using MMD values')
    n_nhbrs     = 5
    edge_list   = lapply(1:nrow(mmd_mat), 
        function(i) {
            row     = mmd_mat[i, ]
            nbrs    = order(row)[2:(n_nhbrs+1)]
            return(data.frame(i=i, j=nbrs))
        }) %>% rbindlist
    graph_obj   = graph_from_edgelist(as.matrix(edge_list), directed=FALSE)
    # graph_obj   = graph_from_adjacency_matrix(abs(mmd_mat), mode='undirected', weighted=TRUE)

    # do clustering
    cluster_obj = cluster_louvain(graph_obj)
    mmd_clusts  = cluster_obj$membership
    # mmd_clusts  = paste0('QC', graph_obj$membership)

    # assemble outputs
    mmd_list    = list(
        mmd_dt      = mmd_dt,
        mmd_mat     = mmd_mat, 
        mmd_clusts  = mmd_clusts,
        params      = list(sigma=sigma, subsample=subsample, n_times=n_times, centre_samples=centre_samples, scale_samples=scale_samples)
        )

    return(mmd_list)
}

#' Calculates non-linear embeddings of MMD distances calculated by 
#' calculate_sample_to_sample_MMDs
#' 
#' @param mmd_list Outputs from function calculate_sample_to_sample_MMDs
#' @param qc_dt Data.table of QC metrics for all cells and samples
#' @param annot_discrete List of discrete-valued variables to include for plotting
#' @param annot_cont List of continuous-valued variables to include for plotting
#' @param n_nhbrs How many neighbours in nearest neighbours graph used as UMAP input
#' 
#' @section Details
#' [Maximum mean discrepancy (MMD) ]
#' [Why to centre and not scale]
#' 
#' @importFrom assertthat assert_that
#' @importFrom uwot umap
#' @importFrom scales rescale
#' @importFrom magrittr "%>%"
#' @importFrom data.table data.table
#' @importFrom data.table setnames
#' @importFrom data.table ":="
#' @importFrom forcats fct_infreq
#' @importFrom data.table setorder
#'
#' @return data.table containing selected 2D embedding, plus annotations per sample
#' @export
embed_sample_to_sample_MMDs <- function(mmd_list, qc_dt, annot_discrete=NULL, annot_cont=NULL, n_nhbrs=5) {
    # unpack 
    assert_that( !is.null(mmd_list$mmd_mat) )
    assert_that( !is.null(mmd_list$mmd_clusts) )
    mmd_mat     = mmd_list$mmd_mat
    mmd_clusts  = mmd_list$mmd_clusts

    # calculate selected embedding
    k           = 2

    # embed with multidimensional scaling
    mds_mat     = cmdscale(mmd_mat, k)

    # convert mmd mat to nearest neighbours
    nn_method   = list(
        idx     = apply(mmd_mat, 1, function(row) order(row)[2:(n_nhbrs+1)]) %>% t,
        dist    = apply(mmd_mat, 1, function(row) sort(row)[2:(n_nhbrs+1)]) %>% t
        )
    umap_mat    = umap(NULL, nn_method=nn_method, min_dist=0.001)

    # scale nicely
    mds_mat     = apply(mds_mat, 2, rescale, to=c(0.05, 0.95))
    umap_mat    = apply(umap_mat, 2, rescale, to=c(0.05, 0.95))

    # annotate with values from qc_dt
    samples     = colnames(mmd_mat)
    n_samples   = length(samples)
    embed_dt    = data.table(rbind(mds_mat, umap_mat)) %>%
        setnames(names(.), paste0('dim', 1:k)) %>%
        .[, embedding := rep(c('MDS', 'UMAP'), each=n_samples)] %>%
        .[, sample_id := rep(samples, 2)] %>%
        .[, QC_clust := fct_infreq(factor(rep(mmd_clusts, 2)))]
    levels(embed_dt$QC_clust)   = paste0('QC', 1:max(mmd_clusts))

    # warn about missing annotations
    dt_names        = c(colnames(qc_dt), 'QC_clust')
    missing_disc    = setdiff(annot_discrete, dt_names)
    if (length(missing_disc)>0)
        warning("the following variables were requested as discrete annotations, but aren't present in qc_dt:\n", paste(missing_disc, collapse=', '))

    missing_cont    = setdiff(annot_cont, dt_names)
    if (length(missing_cont)>0)
        warning("the following variables were requested as continuous annotations, but aren't present in qc_dt:\n", paste(missing_cont, collapse=', '))

    # define automatic names
    auto_discrete   = c('QC_clust', 'N_cat', 'mito_cat', 'counts_cat')
    auto_cont       = c('log_N', 'med_mito', 'med_counts')

    # add automatic names
    annot_discrete  = c(annot_discrete, setdiff(auto_discrete, annot_discrete))
    annot_discrete  = intersect(annot_discrete, dt_names)
    annot_cont      = c(annot_cont, setdiff(auto_cont, annot_cont))
    annot_cont      = intersect(annot_cont, dt_names)

    # get annotations
    annots_dt   = qc_dt[, c('sample_id', setdiff(annot_discrete, 'QC_clust'), annot_cont), with=FALSE] %>%
        unique
    embed_dt    = data.table:::merge.data.table(annots_dt, embed_dt, by='sample_id') %>%
        setorder('embedding', 'sample_id')

    # add this stuff to mmd_list
    mmd_list$embed_dt       = embed_dt
    mmd_list$annot_discrete = annot_discrete
    mmd_list$annot_cont     = annot_cont

    return(mmd_list)
}

#' Plots non-linear embeddings of MMD distances calculated by 
#' embed_sample_to_sample_MMDs
#' 
#' @param embed_dt Output from embed_sample_to_sample_MMDs
#' @param annot_vars Which set of annotation variables to plot
#' @param var_type Is this annotation discrete or continuous?
#' @param sel_embed Which type of embedding to use
#' 
#' @importFrom assertthat assert_that
#' @return NULL
#' @export
plot_embeddings <- function(mmd_list, var_type, sel_embed) {
    assert_that( var_type %in% c('discrete', 'continuous') )
    if (var_type=='discrete') {
        annot_vars  = mmd_list$annot_discrete
    } else {
        annot_vars  = mmd_list$annot_cont
    }
    for (v in annot_vars) {
        cat('### ', v, '\n')
        g = .plot_one_embedding(mmd_list, v, var_type, sel_embed)
        print(g)
        cat('\n\n')
    }
}

#' Plots non-linear embeddings of MMD distances calculated by 
#' embed_sample_to_sample_MMDs
#' 
#' @param embed_dt Output from embed_sample_to_sample_MMDs
#' @param annot_var Which annotation variable to plot
#' @param var_type Is this annotation discrete or continuous?
#' @param sel_embed Which type of embedding to use
#' 
#' @importFrom assertthat assert_that
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales pretty_breaks
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#'
#' @return ggplot object
#' @keywords internal
.plot_one_embedding <- function(mmd_list, v, var_type, sel_embed) {
    # unpack
    embed_dt    = mmd_list$embed_dt

    # checks
    assert_that( var_type %in% c('discrete', 'continuous') )
    assert_that( v %in% names(embed_dt) )
    assert_that( sel_embed %in% c('MDS', 'UMAP') )

    # restrict to chosen embedding
    dt  = embed_dt[ embedding==sel_embed ]
    if (var_type=='discrete') {
        n_cols      = length(unique(dt[[v]]))
        if (v %in% c('mito_cat', 'N_cat')) {
            col_vals    = rev(brewer.pal(n_cols, 'PiYG'))
        } else {
            col_vals    = .CLUSTER_COLS
        }
        g = ggplot(dt) + 
            aes_string( x='dim1', y='dim2', label='sample_id', fill=paste0("factor(", v, ")")) + 
            geom_text( hjust=0.5, vjust=1.5, colour='black', size=2, alpha=0.8 ) + 
            geom_point( size=3, alpha=0.8, shape=21, colour='black' ) +
            scale_x_continuous( breaks=pretty_breaks() ) + scale_y_continuous( breaks=pretty_breaks() )
        if (n_cols > 2) {
            g = g + scale_fill_manual( values=col_vals[1:n_cols] )
        } else {
            g = g + scale_fill_brewer( palette='Set1' )
        }
        g = g + 
            theme_bw() + theme( aspect.ratio=1, axis.text=element_blank() ) + 
            labs( fill=v )
    } else if (var_type=='continuous') {
        g = ggplot(dt) + 
            aes_string( x='dim1', y='dim2', label='sample_id', fill=v) + 
            geom_text( hjust=0.5, vjust=1.5, colour='black', size=2, alpha=0.8 ) +
            geom_point( size=3, alpha=0.8, shape=21, colour='black' ) +
            scale_fill_distiller( palette='RdBu', breaks=pretty_breaks() ) + 
            scale_x_continuous( breaks=pretty_breaks() ) + scale_y_continuous( breaks=pretty_breaks() ) +
            theme_bw() + theme( aspect.ratio=1, axis.text=element_blank() )
    }
    g = g + labs( x=paste0(sel_embed, '1'), y=paste0(sel_embed, '2') )

    return(g)
}

#' Calculates non-linear embeddings of MMD distances calculated by 
#' calculate_sample_to_sample_MMDs
#' 
#' @param mmd_list Outputs from function calculate_sample_to_sample_MMDs
#' @param qc_dt Data.table of QC metrics for all cells and samples
#' @param annot_discrete List of discrete-valued variables to include for plotting
#' @param annot_cont List of continuous-valued variables to include for plotting
#' @param n_nhbrs How many neighbours in nearest neighbours graph used as UMAP input
#' 
#' @section Details
#' [Maximum mean discrepancy (MMD) ]
#' [Why to centre and not scale]
#' 
#' @importFrom assertthat assert_that
#' @importFrom data.table copy
#' @importFrom data.table ":="
#' @importFrom scales pretty_breaks
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_distiller
#' @importFrom ggplot2 expand_limits
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 labs
#'
#' @return ggplot object
#' @export
plot_mmd_heatmap <- function(mmd_list) {
    # unpack 
    assert_that( !is.null(mmd_list$mmd_mat) )
    assert_that( !is.null(mmd_list$mmd_dt) )
    assert_that( !is.null(mmd_list$params) )
    mmd_dt      = mmd_list$mmd_dt
    mmd_mat     = mmd_list$mmd_mat
    n_times     = mmd_list$params$n_times
    
    # do hierarchical clustering, use this to order the samples
    hclust_obj  = hclust(as.dist(mmd_mat), method='complete')
    sample_ord  = hclust_obj$labels[hclust_obj$order]
    plot_dt     = copy(mmd_dt)
    plot_dt[, sample_i := factor(sample_i, levels=sample_ord) ] %>%
        .[, sample_j := factor(sample_j, levels=sample_ord) ] %>%
        .[ sample_i != sample_j ]

    # plot
    g = ggplot(plot_dt) +
        aes( x=sample_i, y=sample_j, fill=mmd_mean ) +
        geom_tile() +
        scale_fill_distiller( palette='PiYG', direction=-1, breaks=pretty_breaks() ) + 
        expand_limits( fill=0 ) +
        theme_bw() + theme( 
            aspect.ratio    = 1,
            axis.text.x     = element_blank(), 
            axis.text.y     = element_text(size=6) ) +
        labs( fill=sprintf('mean MMD\n(over %d samples)', n_times))

    return(g)
}
