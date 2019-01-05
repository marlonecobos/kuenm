#' An ecological niche model created with Maxent with raw output
#'
#' A RasterLayer containing an ecological niche model for the a tick (*Amblyomma americanum*)
#' that was created with all occurrences and raw output.
#'
#' @name sp_mod_joint
#'
#' @format A RasterLayer with 150 rows, 249 columns, and 37350 cells:
#' \describe{
#'   \item{Suitability}{suitability, in probability values.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' model <- raster::raster(list.files(system.file("extdata", package = "kuenm"),
#'                                    pattern = "sp_model_joint.tif", full.names = TRUE))
#'
#' raster::plot(model)
NULL
