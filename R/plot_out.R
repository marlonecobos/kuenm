#' Detection of environmental values ouside the calibration area of a model
#'
#' @description plot.out for calculating a mobility-oriented parity layer.
#' This function is designed to be used specifically in the \code{\link{kuenm_mop}} function.
#'
#' @param M1 a numeric matrix containing values of all environmental variables in the calibration area.
#' @param G1 a numeric matrix containing values of all environmental variables in the full area of interest.


plot_out <- function (M1, G1) {
  d1 <- dim(M1)
  AllVec <- matrix(0, 0, 0)

  for (i in 3:d1[2]) {
    MRange <- range(M1[, i])
    l1 <- which(G1[, i] < range(M1[, i])[1] | G1[,4] > range(M1[, 4])[2])
    AllVec <- c(l1, AllVec)
  }

  AllVec <- unique(AllVec)

  return(AllVec)
}
