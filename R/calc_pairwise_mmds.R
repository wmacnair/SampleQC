#' calc_pairwise_mmds
#'
#' Calculates sample-to-sample MMD distances, clusters and embeds them
#' 
#' @description 
#' Maximum mean discrepancy (MMD) is a method for measuring distances between 
#' multivariate distributions. \emph{SampleQC} uses this to determine which 
#' samples have similar distributions, and should therefore have the same 
#' model fit to them. In addition
#'
#' @details 
#' The most important variable to consider here is \emph{qc_names}, which 
#' specifies the set of QC metrics that will be used both here and in later 
#' steps. The default is to use log counts, log features and the logit 
#' transformation of the mitochondrial proportion (under this transformation, 
#' the distribution is reasonably modelled as a gaussian mixture model)
#'
#' @param qc_dt data.table of QC metrics for all cells and samples
#' @param qc_names List of metrics to actually use for calculating 
#' sample-to-sample distances
#' @param annots_disc,annots_cont Which annotating variables should be added 
#' to the SampleQC object? These should be entries in qc_dt. annots_disc 
#' indicates discrete variables; annots_cont continuous variables.
#' @param one_group_only A TRUE/FALSE variable specifying whether samples 
#' should be first clustered into groups, or assumed to form one group. The 
#' default is to cluster when there more than 10 samples, and not cluster when 
#' there are fewer than ten samples.
#' @param n_cores How many cores to parallelize over?
#' @param sigma Scale for MMD kernel (defaults to length of qc_names)
#' @param subsample Should we downsample the number of cells per sample to 
#' this number
#' @param n_times How many times do we sample MMD between each pair of samples?
#' (if subsampled, MMD is a stochastic value)
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
calc_pairwise_mmds <- function(qc_dt, qc_names=c('log_counts', 
    'log_feats', 'logit_mito'), annots_disc=NULL, annots_cont=NULL, 
    one_group_only = NULL, n_cores=4, sigma=length(qc_names), 
    subsample=100, n_times=20, seed = 22) {
    # check some inputs
    if( missing(sigma) )
        sigma   = length(qc_names)
    assert_that( length(sigma) == 1 )
    assert_that( sigma > 0 )

    # get list of samples to compare
    sample_list = sort(unique(qc_dt$sample_id))
    n_samples   = length(sample_list)
    if( is.null(one_group_only) )
        one_group_only  = n_samples <= 10
    assert_that( is.flag(one_group_only) )

    message('sample-level parts of SampleQC:')

    # split qc metric values into one matrix per sample
    mat_list    = .calc_mat_list(qc_dt, qc_names, sample_list)

    # calculate mmds
    mmd_mat     = .calc_mmd_mat(sample_list, mat_list, 
        n_times, subsample, sigma, n_cores)

    # cluster on basis of mmd_mat
    set.seed(seed)
    mmd_graph   = .make_mmd_graph(mmd_mat, n_nhbrs=5)
    group_ids   = .cluster_mmd_graph(mmd_graph, one_group_only)

    # do embeddings
    mds_mat     = .embed_mds(mmd_mat)
    n_nhbrs     = 5
    umap_mat    = .embed_umap(mmd_mat, n_nhbrs)

    # define automatic names
    message('  adding annotation variables')
    auto_disc   = c('N_cat', 'counts_cat', 'feats_cat', 'mito_cat', 'splice_cat')
    auto_cont   = c('log_N', 'med_counts', 'med_feats', 'med_mito', 'med_splice')
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
        n_times         = n_times
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
#'
#' @importFrom magrittr "%>%"
#' 
#' @return List of matrices
#' 
#' @keywords internal
.calc_mat_list <- function(qc_dt, qc_names, sample_list) {
    # do PCA on overall matrix
    qc_mat      = as.matrix(qc_dt[ , qc_names, with=FALSE ])
    pca_obj     = prcomp(qc_mat, center = TRUE, scale. = TRUE)
    pca_qc      = pca_obj$x

    # scale
    pca_qc      = apply(pca_qc, 2, scale)

    # make list of individual matrices, scaled if necessary
    sample_ids  = qc_dt$sample_id
    mat_list    = lapply(sample_list, function(s) {
        qc_mat  = pca_qc[ sample_ids == s, ]
        qc_mat  = apply(qc_mat, 2, scale, center = TRUE, scale = FALSE)
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
#' @param bp_seed random seed for BiocParallel workers
#'
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam bpstart bpstop
#' @importFrom data.table data.table dcast
#' @importFrom magrittr "%>%" set_rownames
#' 
#' @return matrix of MMD values
#' @keywords internal
.calc_mmd_mat <- function(sample_list, mat_list, n_times, subsample, sigma, n_cores, bp_seed=22) {
    # define combns to do
    ij_grid     = expand.grid(
        sample_i=factor(sample_list), 
        sample_j=factor(sample_list)
        )
    ij_grid     = ij_grid[ as.integer(ij_grid$sample_i) < 
        as.integer(ij_grid$sample_j), ]

    # define parallelization
    if (n_cores == 1) {
        bpparam = SerialParam()
    } else {
        bpparam = MulticoreParam(workers = n_cores, RNGseed = bp_seed, 
            progressbar = TRUE)
    }
    # iterate through them
    n_combns    = nrow(ij_grid)
    message('  calculating ', n_combns, 
        ' sample-sample MMDs:')
    bpstart(bpparam)
    mmd_means   = bplapply(
        seq_len(n_combns), function(r) {
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
        }, BPPARAM = bpparam)
    bpstop(bpparam)

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
#' @param mmd_mat matrix of MMD distances
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
#' @importFrom igraph gorder cluster_louvain
#' @importFrom magrittr "%>%"
#' @importFrom stringr str_match
#' 
#' @return vector of group_ids for all samples
#' @keywords internal
.cluster_mmd_graph <- function(mmd_graph, one_group_only) {
    message('  clustering samples using MMD values')
    # if not clustering, return one group for every sample
    if (one_group_only == TRUE) {
        group_ids   = rep('SG1', gorder(mmd_graph)) %>% factor
        return(group_ids)
    }

    # otherwise do clustering
    cluster_obj = cluster_louvain(mmd_graph)
    group_ids   = cluster_obj$membership

    # rename in order
    group_ns    = table(group_ids)
    new_ids     = seq.int(length(unique(group_ids))) %>%
        setNames(names(sort(-group_ns)))
    group_ids   = new_ids[ as.character(group_ids) ] %>% unname
    group_ids   = paste0('SG', group_ids) %>% factor

    # put levels in nice order
    new_levels  = levels(group_ids) %>%
        str_match('[0-9]+') %>%
        as.integer %>%
        sort %>%
        paste0('SG', .)
    levels(group_ids) = new_levels

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
#' calc_pairwise_mmds
#' 
#' @param mmd_mat MMD distances
#' @param n_nhbrs How many neighbours for UMAP nearest neighbours graph
#' 
#' @importFrom uwot umap
#' @importFrom scales rescale
#' @importFrom magrittr "%>%"
#'
#' @return matrix of embedding
#' @keywords internal
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
    umap_mat    = umap(as.dist(mmd_mat), 
        n_neighbors=n_nhbrs, min_dist=0.001)

    # tidy up
    umap_mat    = apply(umap_mat, 2, rescale, to=c(0.05, 0.95))
    colnames(umap_mat)   = paste0('UMAP', 1:ncol(umap_mat))

    return(umap_mat)
}

#' Checks annotation variables
#' 
#' @param qc_dt \code{data.table} of QC metrics
#' @param specified Requested annotation variables
#' @param auto Automatically generated annotation variables
#' 
#' @return character vector
#' @keywords internal
.process_annot_vars <- function(qc_dt, specified, auto) {
    # warn about missing annotations
    dt_names        = colnames(qc_dt)
    missing_vals    = setdiff(specified, dt_names)
    if (length(missing_vals) > 0)
        warning("These variables were requested as annotations, but 
            aren't present in qc_dt:\n", paste(missing_vals, collapse=', '))

    # check sample_id wasn't specified
    if ('sample_id' %in% specified)
        warning(paste0("'sample_id' was specified as an annotation. This", 
            "variable is reserved for determining samples, so was removed ",
            "from the annotations."))
    specified       = setdiff(specified, 'sample_id')

    # add automatic names
    annot_vars      = c(specified, setdiff(auto, specified))
    annot_vars      = intersect(annot_vars, dt_names)

    return(annot_vars)
}
