#' Preparing an XY coordinates table for calculation of the partial ROC
#'
#' @description g_xy_tab is a function designed to be used specifically in the
#' \code{\link{kuenm_proc}} function. It generates an XY coordinates table to
#' calculate partial areas under the ROC curve.
#'
#' @param classpixels table with information about proportionate areas predicted under each
#' value of suitability.
#' @param occurinclass table with the percentage of points within each class of suitability.

g_xy_tab <- function(classpixels, occurinclass) {
  xytable <- classpixels[, c(1, 4)]
  xytable <- cbind(xytable,rep(-1, nrow(xytable)))

  ## Set the previous value for 1-omission, i.e Y-axis as the value of last
  ## class id in Occurrence table. Last class id will always smallest
  ## area predicted presence.
  prevyval <- occurinclass[1, 4]
  for (i in nrow(classpixels):1) {
    curclassid <- xytable[i, 1]
    yval <- occurinclass[which(occurinclass[, 2] == curclassid), 4]

    if (length(yval) == 0 ) {
      xytable[i, 3] <- prevyval
    }else {
      xytable[i, 3] <- yval
      prevyval <- yval
    }
  }

  ## Add A dummy class id in the xytable with coordinate as 0,0
  xytable <- rbind(xytable, c(xytable[nrow(xytable), 1] + 1, 0, 0))
  xytable <- as.data.frame(xytable)
  names(xytable) <- c("ClassID", "XCoor", "YCoor")
  ### Now calculate the area using trapezoid method.
  return(xytable)
}
