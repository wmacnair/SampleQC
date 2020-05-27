# functions for hierarchy of mixtures
generate_data_alpha_j_beta_k <- function(J, K, N, mu_0, alpha_j, beta_k, sigma_k) {
    # define default values
    if (missing(mu_0))
        mu_0    = 1
    if (missing(alpha_j))
        alpha_j = rnorm(J)
    if (missing(beta_k))
        beta_k  = seq_len(K)
    if (missing(sigma_k))
        sigma_k = (seq_len(K)) / 10
    beta_k      = beta_k - mean(beta_k)

    # sample J different groups
    j_weights   = exp(rnorm(J)/2)
    n_js        = as.vector(rmultinom(1, N, prob=j_weights))
    j_vals      = seq_len(J)
    groups      = rep(j_vals, times=n_js)

    # for each group, sample p_jk, z, x
    p_jk        = matrix(NA, J, K)
    z           = matrix(NA, N, 1)
    k_vals      = seq_len(K)
    for (j in seq_len(J)) {
        # get which samples are here
        j_idx       = groups == j

        # sample pi, z
        this_p_jk   = gtools::rdirichlet(1, rep(5, K))
        p_jk[j, ]   = this_p_jk
        z[j_idx, ]  = sample(k_vals, size=sum(j_idx), replace=TRUE, prob=this_p_jk)
    }

    # sample x
    x           = mu_0 + alpha_j[groups] + beta_k[z] + rnorm(N, mean=0, sd=1) * sigma_k[z]
    x           = matrix(x, ncol=1)

    # make output list
    data_list = list(
        J       = J, 
        K       = K, 
        N       = N, 
        x       = x, 
        groups  = groups, 
        mu_0    = mu_0,
        alpha_j = alpha_j,
        beta_k  = beta_k,
        sigma_k = sigma_k,
        p_jk    = p_jk, 
        z       = z
        )
    return (data_list)
}

generate_mvn_alpha_j_beta_k <- function(D, J, K, N, mu_0, alpha_j, beta_k, sigma_k) {
    # define default values
    if (missing(mu_0))
        mu_0    = rep(1, D)
    if (missing(alpha_j))
        alpha_j = matrix(rnorm(J*D)/10, ncol=D)
    if (missing(beta_k))
        beta_k  = matrix(rnorm(K*D)/2, ncol=D)
    if (missing(sigma_k))
        sigma_k = vapply(seq_len(K), function(k) {
            corr_mat        = diag(D)
            sd_vec          = rep(k/10, D)
            sigma_tmp       = corr_mat %>% sweep(1, sd_vec, "*") %>% sweep(2, sd_vec, "*")
            return(sigma_tmp)
        }, matrix(D)) %>% array(dim=c(D, D, K))
    beta_k      = sweep(beta_k, 2, colMeans(beta_k), "-")

    # put ks in correct order
    k_order = order(beta_k[, 1])
    beta_k  = beta_k[ k_order, ]
    sigma_k = sigma_k[ , , k_order ]

    # sample J different groups
    j_weights   = exp(rnorm(J)/2)
    n_js        = as.vector(rmultinom(1, N, prob=j_weights))
    j_vals      = seq_len(J)
    groups      = rep(j_vals, times=n_js)

    # for each group, sample p_jk, z, x
    p_jk        = matrix(NA, J, K)
    z           = matrix(NA, N, 1)
    k_vals      = seq_len(K)
    for (j in seq_len(J)) {
        # get which samples are here
        j_idx       = groups == j

        # sample pi, z
        this_p_jk   = gtools::rdirichlet(1, rep(5, K))
        p_jk[j, ]   = this_p_jk
        z[j_idx, ]  = sample(k_vals, size=sum(j_idx), replace=TRUE, prob=this_p_jk)
    }

    # sample x
    x           = matrix(nrow=N, ncol=D)
    for (k in seq_len(K)) {
        # which entries, how many?
        k_idx       = z[, 1] == k
        n_k         = sum(k_idx)

        # prep mu
        mu_0_mat    = matrix(rep(mu_0, n_k), ncol=D, byrow=TRUE)
        alpha_j_mat = alpha_j[groups[k_idx], ]
        beta_k_mat  = matrix(rep(beta_k[k, ], n_k), ncol=D, byrow=TRUE)

        # prep sigma
        sigma_k_    = sigma_k[, ,k ]

        # do draw
        x[k_idx, ]  = mu_0_mat + alpha_j_mat + beta_k_mat + mvtnorm::rmvnorm(sum(k_idx), mean=rep(0, D), sigma=sigma_k_)
    }

    # make output list
    data_list = list(
        D       = D, 
        J       = J, 
        K       = K, 
        N       = N, 
        x       = x, 
        groups  = groups, 
        mu_0    = mu_0,
        alpha_j = alpha_j,
        beta_k  = beta_k,
        sigma_k = sigma_k,
        p_jk    = p_jk, 
        z       = z
        )
    return (data_list)
}

make_truth_alpha_j_beta_k <- function(data_list) {
    truth   = data.table(
        var     = c(
            paste0('mu_0'),
            paste0('alpha_j[', seq_len(data_list$J), ']'),
            paste0('beta_k[', seq_len(data_list$K), ']'),
            paste0('sigma_k[', seq_len(data_list$K), ']'),
            as.vector(outer(seq_len(data_list$J), seq_len(data_list$K), function(j, k) paste0('p_jk[', j, ",", k, ']')))
            ),
        value   = c(
            data_list$mu_0, 
            data_list$alpha_j, 
            data_list$beta_k, 
            data_list$sigma_k,
            as.vector(data_list$p_j)
            )
        )
    # truth[ , var    := factor(var) ] %>%
    #     .[ , var_type   := factor(str_match(var, '^(.+)\\[')[,2]) ] %>%
    #     .[ , comp       := factor(str_match(var, '([0-9]+)\\]')[,2]) ] %>%
    #     .[ , group      := factor(str_match(var, 'j_\\[([0-9]+),')[,2]) ]
    truth   = .annotate_variables(truth)
    return (truth)
}

make_em_dt <- function(data_list, em_list) {
    em_dt   = data.table(
        var     = c(
            paste0('mu_0'),
            paste0('alpha_j[', seq_len(data_list$J), ']'),
            paste0('beta_k[', seq_len(data_list$K), ']'),
            paste0('sigma_k[', seq_len(data_list$K), ']'),
            paste0('p_k[', seq_len(data_list$K), ']'),
            as.vector(outer(seq_len(data_list$J), seq_len(data_list$K), function(j, k) paste0('p_jk[', j, ",", k, ']')))
            ),
        value   = c(
            em_list$mu_0, 
            em_list$alpha_j, 
            em_list$beta_k, 
            em_list$sigma_k,
            em_list$p_k,
            as.vector(em_list$p_jk)
            )
        )
    em_dt   = .annotate_variables(em_dt)
    return (em_dt)
}

calc_means <- function(dt) {
    mu_0        = dt[var == 'mu_0']$value
    alphas      = dt[ var_type == 'alpha_j', list(group, alpha_j=value, dummy=1) ]
    betas       = dt[ var_type == 'beta_k', list(comp, beta_k=value, dummy=1) ]
    means_dt = merge(alphas, betas, by='dummy', allow.cartesian=TRUE) %>% .[, dummy := NULL ]
    means_dt[, mu_jk := mu_0 + alpha_j + beta_k ]
    return(means_dt)
}

.annotate_variables_mvn <- function(dt) {
    dt[ , var    := factor(var) ] %>%
    .[ , var_type   := factor(str_match(var, '^(.+)[\\[$]')[,2]) ] %>%
    .[ ,            comp    := factor(str_match(var, '_k\\[([0-9]+),[0-9]+\\]')[,2]) ] %>%
    .[is.na(comp),  comp    := factor(str_match(var, 'jk\\[[0-9]+,([0-9]+)\\]')[,2]) ] %>%
    .[ ,            dim     := factor(str_match(var, '_0\\[([0-9]+)\\]')[,2]) ] %>%
    .[is.na(dim),   dim     := factor(str_match(var, '_j\\[[0-9]+,([0-9]+)\\]')[,2]) ] %>%
    .[is.na(dim),   dim     := factor(str_match(var, '_k\\[[0-9]+,([0-9]+)\\]')[,2]) ] %>%
    .[ ,            group   := factor(str_match(var, 'j\\[([0-9]+)')[,2]) ] %>%
    .[is.na(group), group   := factor(str_match(var, 'jk\\[([0-9]+),')[,2]) ]
    return(dt)
}

extract_params_mvn <- function(data_list, param_list) {
    dt   = data.table(
        var     = c(
            paste0('mu_0[', seq_len(data_list$D), ']'),
            as.vector(outer(seq_len(data_list$J), seq_len(data_list$D), function(j, d) paste0('alpha_j[', j,',', d, ']'))),
            as.vector(outer(seq_len(data_list$K), seq_len(data_list$D), function(k, d) paste0('beta_k[', k,',', d, ']'))),
            as.vector(outer(seq_len(data_list$D^2), seq_len(data_list$K),  function(k, d) paste0('sigma_k[', k,',', d, ']'))),
            as.vector(outer(seq_len(data_list$J), seq_len(data_list$K), function(j, k) paste0('p_jk[', j, ",", k, ']')))
            ),
        value   = c(
            param_list$mu_0,
            as.vector(param_list$alpha_j),
            as.vector(param_list$beta_k),
            as.vector(param_list$sigma_k),
            as.vector(param_list$p_j)
            )
        )
    # dt[ , var    := factor(var) ] %>%
    #     .[ , var_type   := factor(str_match(var, '^(.+)\\[')[,2]) ] %>%
    #     .[ , comp       := factor(str_match(var, '([0-9]+)\\]')[,2]) ] %>%
    #     .[ , group      := factor(str_match(var, 'j_\\[([0-9]+),')[,2]) ]
    dt   = .annotate_variables_mvn(dt)
    return (dt)
}

.annotate_variables <- function(dt) {
    dt[ , var    := factor(var) ] %>%
    .[ , var_type   := factor(str_match(var, '^(.+)[\\[$]')[,2]) ] %>%
    .[ ,            comp := factor(str_match(var, 'k\\[([0-9]+)\\]')[,2]) ] %>%
    .[is.na(comp),  comp := factor(str_match(var, 'jk\\[[0-9]+,([0-9]+)\\]')[,2]) ] %>%
    .[ ,            group := factor(str_match(var, 'j\\[([0-9]+)')[,2]) ] %>%
    .[is.na(group), group := factor(str_match(var, 'jk\\[([0-9]+),')[,2]) ]
    return(dt)
}

calc_means_mvn <- function(dt) {
    mu_0        = dt[ var_type == 'mu_0', list(dim, mu_0=value, dummy=1) ]
    alphas      = dt[ var_type == 'alpha_j', list(dim, group, alpha_j=value, dummy=1) ]
    betas       = dt[ var_type == 'beta_k', list(dim, comp, beta_k=value, dummy=1) ]
    means_dt    = merge(
        alphas, betas, 
        by=c('dim'), allow.cartesian=TRUE) %>% 
        merge(
            mu_0, 
            by=c('dim'), allow.cartesian=TRUE) %>% 
        .[, dim := paste0('dim', dim) ] %>%
        .[, mu_jk := mu_0 + alpha_j + beta_k ] %>% 
        dcast( group + comp ~ dim, value.var='mu_jk' ) 
    return(means_dt)
}

#' Generate fake qc_df object for testing
#'
#' @importFrom data.table data.table rbindlist ":="
#' @return data.frame
#' @keyword internal
.make_toy_qc_df <- function() {
    # generate sample means
    J               = 10
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
#' @keyword internal
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

