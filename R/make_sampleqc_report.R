#' make_sampleqc_report
#'
#' Renders SampleQC report into specified format.
#' 
#' @param qc_obj Outputs from fit_sampleqc
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
make_sampleqc_report <- function(qc_obj, save_dir, proj_name) {
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
        render(tmplt_file, output_file=report_file, output_dir=save_dir, quiet=TRUE)
    }, error=function(cond) {
        message("Something went wrong with rendering your report :/")
        message("Here's the original error message:")
        message(cond)
        message()
        return(NA)
    })
}
