#' Extrapolation risk analysis for single comparisons
#'
#' @description kuenm_mop calculates a mobility-oriented parity layer by
#' comparing environmental values between the calibration area and the area or
#' scenario to which an ecological niche model is transferred.
#'
#' @param M.stack a RasterStack of variables representing the calibration area.
#' @param G.stack a RasterStack of variables representing the full area of interest, and areas
#' or scenarios to which models are transferred.
#' @param percent (numeric) percent of values sampled from te calibration region to calculate the MOP.
#' @param normalized (logical) if true values of similarity are presented from 0 to 1,
#' default = TRUE.
#'
#' @return A mobility-oriented parity RasterLayer.
#'
#' @details The MOP is calculated following Owens et al.
#' (2013; \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}). This function is a modification
#' of the \code{\link[ENMGadgets]{MOP}} funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
#'
#' @examples
#' mvars <- mvars_mop
#' gvars <- gvars_mop
#' perc <- 10
#' norm <- TRUE
#'
#' mop <- kuenm_mop(M.stack = mvars, G.stack = gvars,
#'                   percent = perc, normalized = norm)

kuenm_mop <- function(M.stack, G.stack, percent = 10, normalized = TRUE) {
  mPoints <- raster::rasterToPoints(M.stack)
  gPoints <- raster::rasterToPoints(G.stack)

  m1 <- mPoints[, -(1:2)]
  m2 <- gPoints[, -(1:2)]

  if(dim(m1)[2] != dim(m2)[2]) {
    stop("Stacks must have the same dimensions")
  }

  steps <- seq(1, dim(m2)[1], 1000)
  kkk <- c(steps,  dim(m2)[1] + 1)
  out_index <- plot_out(m1, m2)
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
  # What about using -1 instead of NA?
  mop_all[out_index, 3] <- NA
  sp::coordinates(mop_all) <- ~x + y
  sp::gridded(mop_all) <- TRUE
  mop_raster <- raster::raster(mop_all)

  if(normalized) {
    mop_raster <- 1 - (mop_raster / mop_max)
  }

  return(mop_raster)
}
