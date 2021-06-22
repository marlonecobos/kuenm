#' Creation of an HTML file with results from model calibration
#'
#' @description html_eval creates an HTML file that summarizes all outputs from
#' the model calibration and selection processes.
#'
#' @param path directory in which the HTML file will be written; default = current directory.
#' @param file.name (character) name of the HTML file.
#'
#' @return An HTML file summarizing results from the model calibration and selection process.
#'
#' @details This function is used along with the \code{\link{kuenm_ceval}} function.
#' @export
#'
#' @examples
#' path <- getwd() # directory with outputs of the kuenm_ceval function
#' name <- "evaluation_results"
#'
#' html_file <- html_eval(path = path, file.name = name)

html_eval <- function(path = getwd(), file.name) {

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
