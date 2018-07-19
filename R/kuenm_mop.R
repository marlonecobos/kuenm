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
#'
#' @return A mobility-oriented parity RasterLayer normalized where values of 0 represent
#' strict extrapolation.
#'
#' @details The MOP is calculated following Owens et al.
#' (2013; \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}). This function is a modification
#' of the \code{\link[ENMGadgets]{MOP}} funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
#'
#'
#' @examples
#' mvars <- mvars_mop
#' gvars <- gvars_mop
#' perc <- 10
#'
#' mop <- kuenm_mop(M.stack = mvars, G.stack = gvars,
#'                   percent = perc)

kuenm_mop <- function(M.stack, G.stack, percent = 10) {
  mPoints <- raster::rasterToPoints(M.stack)
  gPoints <- raster::rasterToPoints(G.stack)

  m1 <- mPoints[, -(1:2)]
  m2 <- gPoints[, -(1:2)]

  if(dim(m1)[2] != dim(m2)[2]) {
    stop("Stacks must have the same dimensions")
  }

  steps <- seq(1, dim(m2)[1], comp_each)
  kkk <- c(steps,  dim(m2)[1] + 1)
  out_index <- plot_out(m1, m2)
  long_k <- length(kkk)
  suppressPackageStartupMessages(library("doParallel"))

  cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  mop1 <- lapply(1:(length(kkk) - 1), function(x) {
    seq_rdist <- kkk[x]:(kkk[x + 1] - 1)
    eudist <- fields::rdist(m2[seq_rdist, ], m1)
    percent <- percent
    mean_quantile <- foreach::foreach(y = 1:dim(eudist)[1],
                                      .packages = c("kuenm")) %dopar% {
                                        mop_dist(eudist_matrix = eudist,
                                                 irow=y,
                                                 percent = percent)

                                      }
    avance <- (x / long_k) * 100
    cat("Computation progress: ", avance,"%" ,"\n")

    return(unlist(mean_quantile))
  })

  parallel::stopCluster(cl)

  mop2 <- unlist(mop1)
  mop_all <- data.frame(gPoints[, 1:2], mop2)
  mop_max <- max(na.omit(mop2))
  mop_max <- max(mop2)
  mop_all[out_index, 3] <- mop_max * 1.05
  sp::coordinates(mop_all) <- ~x + y
  sp::gridded(mop_all) <- TRUE
  mop_raster <- raster::raster(mop_all)

  mop_raster <- 1 - (mop_raster / mop_max)

  return(mop_raster)
}
