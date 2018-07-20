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
#' @param comp.each (numeric) compute distance matrix for a each fixed number of rows (default = 1000).
#' @param parallel (logical) if TRUE, calculations will be performed in parallel using the available
#' cores of the computer. This will demand more RAM and almost full use of the CPU; hence, its use
#' is more recommended in high-performance computers. Using this option will speed up the analyses.
#' Default = FALSE
#'
#' @return A mobility-oriented parity RasterLayer where values of 0 represent strict extrapolation,
#' which means complete dissimilarity of environments between the calibration (M) and projection area (G).
#'
#' @details The MOP is calculated following Owens et al.
#' (2013; \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}). This function is a modification
#' of the \code{\link[ENMGadgets]{MOP}} funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
#'
#' @examples
#' mvars <- mvars_mop
#' gvars <- gvars_mop
#' perc <- 10
#'
#' mop <- kuenm_mop(M.stack = mvars, G.stack = gvars, percent = perc)

kuenm_mop <- function(M.stack, G.stack, percent = 10, comp.each = 1000, parallel = FALSE) {
  mPoints <- raster::rasterToPoints(M.stack)
  m_nona <- na.omit(mPoints)
  m_naID <- attr(m_nona,"na.action")
  gPoints <- raster::rasterToPoints(G.stack)
  g_nona <- na.omit(gPoints)
  g_naID <- attr(g_nona,"na.action")

  m1 <- m_nona[, -(1:2)]
  m2 <- g_nona[, -(1:2)]

  if(dim(m1)[2] != dim(m2)[2]) {
    stop("Stacks must have the same dimensions.")
  }

  steps <- seq(1, dim(m2)[1], comp.each)
  kkk <- c(steps,  dim(m2)[1] + 1)
  out_index <- plot_out(m1, m2)
  long_k <- length(kkk)

  if (parallel == FALSE) {
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

    mop_vals <- unlist(mop1)

  }else {
    suppressPackageStartupMessages(library("future"))
    future::plan(multiprocess)

    mop_env <- new.env()

    pasos <- 1:(length(kkk) - 1)
    pasosChar <- paste0(pasos)

    for (paso in pasosChar) {
      x <- as.numeric(paso)
      mop_env[[paso]] %<-% {
        seq_rdist <- kkk[x]:(kkk[x + 1] - 1)
        eudist <- fields::rdist(m2[seq_rdist, ], m1)
        mop_dist <- lapply(1:dim(eudist)[1], function(y){
          di <- eudist[y, ]
          qdi <- quantile(di, probs = percent / 100,
                          na.rm = TRUE)
          ii <-  which(di <= qdi)
          pond_mean <- mean(di,na.rm = TRUE)
          return(pond_mean)
        })
        mop <-unlist(mop_dist)
        return(mop)
      }
      avance <- (x / long_k) * 100
      cat("Computation progress: ", avance,"%" ,"\n")
    }

    mop_list <- as.list(mop_env)
    mop_names <- sort(as.numeric(names(mop_list)))
    mop_names <- as.character(mop_names)
    mop_vals <- unlist(mop_list[mop_names])
  }

  if(!is.null(g_naID)){
    mop_all <- data.frame(gPoints[,1:2])
    mop_all$mop <- NA
    mop_all$mop[-g_naID] <- mop_vals

  }else{
    mop_all <- data.frame(gPoints[, 1:2], mop = mop_vals)
  }

  mop_max <- max(na.omit(mop_vals)) * 1.05
  mop_all[out_index, 3] <- mop_max
  suppressWarnings({
    sp::coordinates(mop_all) <- ~ x + y
    sp::gridded(mop_all) <- TRUE
    mop_raster <- raster::raster(mop_all)

    mop_raster <- 1 - (mop_raster / mop_max)
  })

  if (parallel == TRUE) {
    if(.Platform$OS.type != "unix"){
      future:::ClusterRegistry("stop")
    }
  }

  return(mop_raster)
}
