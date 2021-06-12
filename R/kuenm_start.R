#' Creation of an R markdown file for recording all analyses
#'
#' @description
#' Generate an R markdown file that serves as a guide for
#' performing most of the analyses included in this package.
#'
#' @param file.name (character) is the name of the R markdown file that will be
#' produced in your working directory. Extension is not needed
#'
#' @return An R markdown file with instructions and code for performing all
#' analyses included in this package.
#'
#' @export
#'
#' @usage
#' kuenm_start(file.name)
#'
#' @rdname kuenm_start
#'
#' @examples
#' kuenm_start(file.name = tempfile())

kuenm_start <- function(file.name){

  if (missing(file.name)) {
    stop("Argument 'file.name' must be defined")
  }

  file.name <- paste0(file.name, ".Rmd")
  suppressMessages(
    file.copy(from = system.file("extdata", "Rmd_start.Rmd",
                                 package = "kuenm"), to = file.name)
  )

  file.edit(file.name)
}


#' @usage
#' kuenm_start_swd(file.name)
#'
#' @rdname kuenm_start

kuenm_start_swd <- function(file.name){

  if (missing(file.name)) {
    stop("Argument 'file.name' must be defined")
  }

  file.name <- paste0(file.name, ".Rmd")
  suppressMessages(
    file.copy(from = system.file("extdata", "Rmd_start_swd.Rmd",
                                 package = "kuenm"), to = file.name)
  )

  file.edit(file.name)
}
