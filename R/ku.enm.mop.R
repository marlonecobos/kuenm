#' Extrapolation risk analysis for single comparisons
#'
#' @description ku.enm.mop calculates a Mobility-Oriented Parity layer by
#' comparing environmental values between the M area and the area or
#' scenario to which an ecological niche model is transferred.
#'
#' @param M.stack a raster stack of variables representing the M area.
#' @param G.stack a raster stack of variables representing the G area, areas or scenarios to
#' which models are transferred.
#' @param percent (numeric) percetage of values sampled from M to calculate the MOP.
#' @param normalized (logical) if true values of similarity are presented from 0 to 1,
#' default = TRUE.
#'
#' @return A Mobility-Oriented Parity raster layer.
#'
#' @details The MOP is calculated following Owens et al., 2013
#' \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}.
#'
#' @examples
#' data()
#' mvars <- mvars_mop
#' gvars <- gvars_mop
#' perc <- 10
#' norm <- TRUE
#'
#' mop <- ku.enm.mop(M.stack = mvars, G.stack = gvars,
#'                   percent = perc, normalized = norm)

ku.enm.mop <- function(M.stack, G.stack, percent = 10, normalized = TRUE) {
  mPoints <- raster::rasterToPoints(M.stack)
  gPoints <- raster::rasterToPoints(G.stack)

  m1 <- mPoints[, -(1:2)]
  m2 <- gPoints[, -(1:2)]

  if(dim(m1)[2] != dim(m2)[2]) {
    stop("Stacks must have the same dimensions")
  }

  steps <- seq(1, dim(m2)[1], 1000)
  kkk <- c(steps,  dim(m2)[1] + 1)
  out_index <- PlotOut(m1, m2)
  long_k <- length(kkk)

  mop1 <- lapply(1:(length(kkk) - 1), function(x) {
    seq_rdist <- kkk[x]:(kkk[x + 1] - 1)
    eudist <- fields::rdist(m2[seq_rdist, ], m1)
    mean_quantile <- parallel::mclapply(1:dim(eudist)[1], function(y) {
      di <- eudist[y, ]
      qdi <- quantile(di, probs = percent / 100, na.rm = TRUE)
      ii <-  which(di <= qdi)
      return(mean(di[ii]))
    })

    avance <- (x / long_k) * 100
    cat("Computation progress: ", avance,"%" ,"\n")

    return(unlist(mean_quantile))
  })

  mop2 <- unlist(mop1)
  mop_all <- data.frame(gPoints[, 1:2], mop2)
  mop_max <- max(na.omit(mop2))
  mop_all[out_index, 3] <- NA
  sp::coordinates(mop_all) <- ~x + y
  sp::gridded(mop_all) <- TRUE
  mop_raster <- raster::raster(mop_all)

  if(normalized) {
    mop_raster <- 1 - (mop_raster / mop_max)
  }

  return(mop_raster)
}


PlotOut <- function (M1, G1) {
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
