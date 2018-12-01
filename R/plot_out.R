#' Helper function for detecting values out of the environmental range of M
#'
#' @description plot.out detects which environmental values in an area of projection are
#' out of the range of environmental values in the area where ecological niche models are
#' calibrated. This function is designed to be used specifically in the \code{\link{kuenm_mop}} function.
#'
#' @param M1 a numeric matrix containing values of all environmental variables in the calibration area.
#' @param G1 a numeric matrix containing values of all environmental variables in the full area of interest.
#'
#' @return A vector of environmental values in a projection area that are outside the range of values
#' in the calibration area of an ecological niche model.
#'
#' @export

plot_out <- function (M1, G1) {
  if(class(M1) == "RasterBrick" | class(M1) == "RasterStack" | class(M1) == "raster"){
    M1 <- raster::values(M1)
  }

  if(class(G1) == "RasterBrick" | class(G1) == "RasterStack" | class(G1) == "raster"){
    G1 <- raster::values(G1)
  }

  d1 <- dim(M1)
  AllVec <- vector()

  for (i in 1:d1[2]) {
    MRange <- range(M1[, i])
    l1 <- which(G1[, i] < range(M1[, i], na.rm = T)[1] | G1[, i] > range(M1[, i], na.rm = T)[2])
    AllVec <- c(l1, AllVec)
  }

  AllVec <- unique(AllVec)

  return(AllVec)
}
