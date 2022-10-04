#' Creation of an HTML file with results from model calibration
#'
#' @description html_calibration creates an HTML file that summarizes all outputs
#' from model calibration, evaluation, and selection.
#'
#' @param path directory where the HTML file will be written; current directory
#' by default.
#' @param file.name (character) name of the HTML file without extension (e.g.,
#' "calibration_results")
#'
#' @return
#' An HTML file summarizing results from model calibration, evaluation, and
#' selection.
#'
#' @details
#' This function is used along with the functions \code{\link{kuenm_ceval}}
#' \code{\link{kuenm_cal_swd}}.
#'
#' @export
#' @importFrom rmarkdown render
#' @importFrom knitr kable
#'
#' @examples
#' path <- getwd() # directory with outputs of the kuenm_ceval function
#' name <- "evaluation_results"
#'
#' \dontrun{
#' html_calibration(path = path, file.name = name)
#' }

html_calibration <- function(path = getwd(), file.name) {

  if (missing(file.name)) {
    stop("Argument 'file.name' must be defined")
  }

  file.name <- paste0(path, "/", file.name, ".Rmd")
  suppressMessages(
    file.copy(from = system.file("extdata", "Rmd_calibration.Rmd",
                                 package = "kuenm"), to = file.name)
  )

  rmarkdown::render(file.name, "html_document", quiet = TRUE)
  unlink(file.name)
}
