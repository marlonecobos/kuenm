#' An ecological niche model created with Maxent
#'
#' A RasterLayer containing an ecological niche model for the a tick (*Amblyomma americanum*)
#' that was created as part of the candidate models during the calibration process.
#'
#' @format A RasterLtack with 150 rows, 249 columns, and 37350 cells:
#' \describe{
#'   \item{Suitability}{suitability, in probability values.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#' @export
#'
#' @examples
#' model <- raster::raster(list.files(system.file("extdata", package = "kuenm"),
#'                                    pattern = "sp_model.tif", full.names = TRUE))
#'
#' raster::plot(model)
"sp_mod"
