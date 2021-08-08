#' plot_embeddings
#'
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
        med_counts  = 'median\nlog10(counts)',
        med_splice  = 'median\nlog2 splice ratio'
    )
    g = g + labs( fill=label_list[[v]] )

    return(g)
}
