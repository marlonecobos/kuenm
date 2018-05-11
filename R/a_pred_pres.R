#' Calculation of areas predicted under each value of suitability
#'
#' @description a_pred_pres calculates proportionate areas predicted under each value of suitability.
#' This function is designed to be used specifically in the \code{\link{kuenm_proc}} function.
#'
#' @param inrast a RasterLayer of the ecological niche model to be evaluated, but standardized
#' previously.

a_pred_pres <- function(inrast) {
  ### Calculate proportionate area predicted under each suitability
  classpixels <- raster::freq(inrast)
  ### Remove the NAs from table
  if (is.na(classpixels[dim(classpixels)[1], 1]) == TRUE)
  {
    classpixels <- classpixels[-dim(classpixels)[1], ]
  }
  classpixels <- classpixels[order(nrow(classpixels):1), ]
  totpixelperclass <- cumsum(classpixels[, 2])
  percentpixels <- totpixelperclass / sum(classpixels[, 2])

  classpixels <- cbind(classpixels, totpixelperclass, percentpixels)
  classpixels <- classpixels[order(nrow(classpixels):1), ]
  return(classpixels)
}
