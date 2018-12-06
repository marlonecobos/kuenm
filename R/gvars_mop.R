#' Variables of the area to which the ecological niche model is transferred.
#'
#' A dataset containing predictor variables of the area of projection of the ecological
#' niche model. Variables represent four future bioclimatic variables (2050) of the
#' NCAR-CCSM4 general circulation model under the RCP 8.5 emission scenario.
#'
#' @format A RasterStack with 900 rows, 2160 columns, 1944000 cells, and 4 layers:
#' \describe{
#'   \item{Temperature}{temperature, in Celsius degrees times 10.}
#'   \item{Precipitation}{precipitation, in milimeters.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#' @export
#'
#' @examples
#' gvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                   pattern = "Gbio_", full.names = TRUE))
#'
#' raster::plot(gvars[[1]])
"gvars_mop"
