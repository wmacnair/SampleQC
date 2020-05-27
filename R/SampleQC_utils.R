# SampleQC: robust multivariate, multi-celltype, multi-sample quality control 
# for single cell RNA-seq
# devtools::load_all('~/work/packages/BayesQC')
# devtools::document('~/work/packages/BayesQC')

# SampleQC_utils.R
# Miscellaneous useful functions 


#' Extracts a data.table of QC metrics, either by calculating them from a 
#' SingleCellExperiment object, or copying them from a data.frame that the
#' user has already made. Also does a little bit of checking of the inputs.
#' 
#' WARNING: This bit of code definitely needs more testing / thought / 
#' improvements... Please let me know if it falls over!
#' @param x either SCE, or data.frame object containing calculated QC metrics
#' @param qc_names list of qc_names that need to be extracted
#' @return qc_dt, a data.table containing the sample variable plus qc metrics
#' @export
make_qc_dt <- function(x, qc_names=c('log_counts', 'log_feats', 'logit_mito')) {
    UseMethod("make_qc_dt")
}

#' Checks that input is ok, puts it into expected format
#'
#' @param x data.frame object containing calculated QC metrics
#' @param qc_names list of qc_names that need to be extracted
#' @importFrom assertthat assert_that
#' @importFrom data.table setcolorder
#' @return qc_dt, a data.table containing the sample variable plus qc metrics
#' @export
make_qc_dt.data.frame <- function(x, qc_names=c('log_counts', 'log_feats', 
    'logit_mito')) {
    # check that cell_id and sample_id are present
    df_names    = colnames(x)
    if ( !('cell_id' %in% df_names) )
        stop('cell_id must be column of data.frame')
    if ( !('sample_id' %in% df_names) )
        stop('sample_id must be column of data.frame')
    if ( length(unique(x$cell_id)) < nrow(x) )
        stop('cell_id values must be unique')

    # put into data.table
    qc_dt           = data.table(x)

    # check which names are present
    non_mito_names  = setdiff(qc_names, c('logit_mito', 'mito_prop'))
    missed_ns       = setdiff(non_mito_names, df_names)
    if (length(missed_ns) > 0)
        stop('the following qc metrics are missing from the input data.frame
            \n:', paste(missed_ns, collapse=', '))

    if ( ( ('logit_mito' %in% qc_names) & !('logit_mito' %in% df_names) ) & 
        ( !('mito_prop' %in% df_names) | !('log_counts' %in% df_names) )
        )
        stop("'logit_mito' found in qc_names but neither 'mito_prop' or 
            'logit_mito' is present in the input data.frame")

    if ('log_counts' %in% qc_names) {
        assert_that(
            all(qc_dt$log_counts >= 0), 
            msg='All log_counts values should be >= 0.'
            )
    } else {
        warning("'log_counts' is not present as a column - are you sure you 
            want to do QC without log_counts?")        
    }

    # add logit_mito if necessary
    if ( ('logit_mito' %in% qc_names) & !('logit_mito' %in% df_names) ) {
        # check that we have log_counts
        assert_that('log_counts' %in% df_names)

        # do calculations
        mito_props          = qc_dt$mito_prop
        assert_that(all(mito_props >= 0) & all(mito_props <= 1))
        total_counts        = 10^qc_dt$log_counts
        qc_dt$logit_mito    = qlogis( (total_counts*mito_props +1) / 
            (total_counts+2) )
    }

    # check logit_mito values
    assert_that(all(is.finite(qc_dt$logit_mito)))

    # add some useful annotations
    qc_dt       = .add_annotations(qc_dt)
    
    # put in nice order
    setcolorder(qc_dt, c('cell_id', 'sample_id', qc_names))

    # check values
    .check_qc_dt(qc_dt, qc_names)

    return(qc_dt)
}

#' Checks that input is ok, puts it into expected format
#'
#' @param x SingleCellExperiment or data.frame (or data.table, some other class
#' inheriting data.frame) containing calculated QC metrics
#' @param qc_names list of qc_names that need to be extracted
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom Matrix colSums
#' @importFrom stringr str_detect
#' @importFrom data.table data.table
#' @importFrom data.table ":=" as.data.table setcolorder
#' @return qc_dt, a data.table containing the sample variable plus qc metrics
#' @export
make_qc_dt.SingleCellExperiment <- function(x, qc_names=c('log_counts', 
    'log_feats', 'logit_mito')) {
    # some checks before starting
    assert_that( is(x, "SingleCellExperiment") )
    assert_that( !is.null(colnames(x)) )
    assert_that( !is.null(x$sample_id) )
    assert_that( length(unique(colnames(x))) == ncol(x) )

    # remove proteins if necessary
    type_col        = rowData(x)$Type
    if (!is.null(type_col)) {
        warning("column 'Type' found in `rowData`; assuming that this is a 
            multi-modal dataset and restricting to rows labelled as 'Gene 
            Expression'")
        warning("to avoid this message, or if this is not what you want to do,
            please make a copy of your sce with only the rows you are 
            interested in, and no 'Type' column in `rowData`")
        x         = x[ type_col == 'Gene Expression', ]
    }

    # warning
    if (!('log_counts' %in% qc_names))
        warning("'log_counts' is not specified as in qc_names - are you sure 
            you want to do QC without log_counts?")

    # restrict to cells >= 1 count
    total_counts    = Matrix::colSums(counts(x))
    if (any(total_counts == 0))
        stop("at least some cells have 0 reads; please remove these before 
            running `make_qc_dt`")

    # define output dt
    qc_dt = data.table(
        cell_id     = colnames(x),
        sample_id   = x$sample_id
        )

    # add specified QC metrics one by one
    if ('log_counts' %in% qc_names)
        qc_dt[, log_counts  := log10(total_counts) ]
    if ('log_feats' %in% qc_names) {
        total_feats     = (counts(x) > 0) %>% Matrix::colSums(.)
        qc_dt[, log_feats   := log10(total_feats) ]
    }
    if ('logit_mito' %in% qc_names) {
        # prepare to calc stats
        idx_mt          = grepl("^mt-", rownames(x), ignore.case = TRUE)
        if (sum(idx_mt) < 13)
            warning('found ', sum(idx_mt), ' mitochondrial genes, which is less
             than the expected 13 human genes; you may want to check your 
             data')

        # calculate stats by hand
        non_mt_counts   = x[ !idx_mt, ] %>% counts %>% Matrix::colSums(.)
        mt_counts       = x[ idx_mt, ] %>% counts %>% Matrix::colSums(.)

        # add all of them
        qc_dt[, log_mito    := log10(mt_counts + 1) ]
        qc_dt[, log_non_mt  := log10(non_mt_counts + 1) ]
        qc_dt[, logit_mito  := qlogis( (mt_counts+1) / (total_counts+2) ) ]
    }

    # add some useful annotations
    qc_dt       = .add_annotations(qc_dt)

    # add any annotations from coldata
    other_vars  = setdiff(names(colData(x)), colnames(qc_dt))
    qc_dt       = cbind(qc_dt, as.data.table(colData(x)[other_vars]))

    # check values
    .check_qc_dt(qc_dt, qc_names)

    return(qc_dt)
}

#' Checks that input is ok, puts it into expected format
#'
#' @param qc_dt data.table
#' @importFrom assertthat assert_that
#' @keywords internal
.check_qc_dt <- function(qc_dt, qc_names) {
    # unpack
    col_names   = colnames(qc_dt)

    # check specific names
    if ('log_counts' %in% col_names)
        assert_that( all(qc_dt$log_counts >= 0) )
    if ('log_feats' %in% col_names)
        assert_that( all(qc_dt$log_feats >= 0) )
    if ('logit_mito' %in% col_names)
        assert_that( all(is.finite(qc_dt$logit_mito)) )

    # check all names
    for (n in qc_names) {
        assert_that( all(!is.na(qc_dt[[n]])) )
    }
}

#' Checks that input is ok, puts it into expected format
#'
#' @param sce SingleCellExperiment or data.frame (or data.table, some other 
#' class inheriting data.frame) containing calculated QC metrics
#' @param qc_names list of qc_names that need to be extracted
#' @importFrom data.table ":="
#' @return qc_dt, a data.table containing the sample variable plus qc metrics
#' @keywords internal
.add_annotations <- function(qc_dt) {
    # add annotations relating to mitochondrial proportions
    if ('log_counts' %in% names(qc_dt) ) {
        # add median log counts per sample
        qc_dt[, med_counts  := median(log_counts), by='sample_id']

        # put mito level into categories
        counts_cuts = c(1,100,300,1000,3000,10000,30000)
        counts_labs = paste0('<=', counts_cuts[-1])
        qc_dt[, counts_cat  := factor(
            cut(10^med_counts, breaks=counts_cuts, labels=counts_labs), 
            levels=counts_labs), by='sample_id']
    }

    # add annotations relating to mitochondrial proportions
    if ('logit_mito' %in% names(qc_dt) ) {
        # add median mito proportion
        qc_dt[, med_mito    := median(plogis(logit_mito)), by='sample_id']

        # put mito level into categories
        mito_cuts   = c(0,0.01,0.05,0.1,0.2,0.5,1)
        mito_labs   = paste0('<=', mito_cuts[-1])
        qc_dt[, mito_cat    := factor(
            cut(med_mito, breaks=mito_cuts, labels=mito_labs), 
            levels=mito_labs), by='sample_id']
    }

    # add annotations for sample size
    qc_dt[, log_N    := log10(.N), by='sample_id']

    # and factor version
    N_cuts      = c(1,100,200,400,1000,2000,4000,10000,20000,40000)
    N_labs      = paste0('<=', N_cuts[-1])
    qc_dt[, N_cat    := factor(
        cut(10^log_N, breaks=N_cuts, labels=N_labs), 
        levels=rev(N_labs)), by=sample_id]

    return(qc_dt)
}

#' Renders SampleQC report into specified format.
#' 
#' @param mmd_list Outputs from calculate_sample_to_sample_MMDs
#' @param em_list Outputs from fit_sampleQC
#' @param save_dir Directory to save into (must exist)
#' @param proj_name Name to use to in file
#' 
#' @section Details
#' 
#' @importFrom assertthat assert_that
#' @importFrom rmarkdown render
#' @import BiocStyle
#' @import patchwork
#'
#' @return NULL
#' @export
make_SampleQC_report <- function(mmd_list, em_list, save_dir, proj_name) {
    # checks
    assert_that( dir.exists(save_dir) )
    # start_dir   = getwd()
    # setwd(save_dir)

    # define output file
    report_file = paste0('SampleQC_report_', proj_name, '.html')

    # render
    tmplt_file  = system.file("Rmd", "SampleQC_report_template.Rmd", 
        package="SampleQC")
    message('rendering ')
    tryCatch({
        render(tmplt_file, output_file=report_file, output_dir=save_dir)
    }, error=function(cond) {
        message("Something went wrong with rendering your report :/")
        message("Here's the original error message:")
        message(cond)
        message()
        return(NA)
    })

    # # go back to original directory
    # setwd(start_dir)
}
