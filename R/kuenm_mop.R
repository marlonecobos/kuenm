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
#' @param comp.each (numeric) compute distance matrix for a each fixed number of rows (default = 2000).
#' @param parallel (logical) if TRUE, calculations will be performed in parallel using the available
#' cores of the computer. This will demand more RAM and almost full use of the CPU; hence, its use
#' is more recommended in high-performance computers. Using this option will speed up the analyses.
#' Default = FALSE.
#'
#' @return A mobility-oriented parity RasterLayer where values of 0 represent strict extrapolation,
#' which means complete dissimilarity of environments between the calibration (M) and projection area (G).
#'
#' @details The MOP is calculated following Owens et al.
#' (2013; \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}). This function is a modification
#' of the \code{\link[ENMGadgets]{MOP}} funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
#'
#' @importFrom future %<-%
#' @import future
#' @export
#'
#' @examples
#' mvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                   pattern = "Mbio_", full.names = TRUE))
#' gvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                   pattern = "Gbio_", full.names = TRUE))
#' perc <- 10
#'
#' mop <- kuenm_mop(M.stack = mvars, G.stack = gvars, percent = perc)
#'
#' raster::plot(mop)

kuenm_mop <- function(M.stack, G.stack, percent = 10, comp.each = 2000, parallel = FALSE) {

  suppressPackageStartupMessages(library(future))

  mop_raster <- G.stack[[1]]
  mValues <- raster::getValues(M.stack)
  m_noNA <- stats::na.omit(mValues)
  m_naIDs <- attr(m_noNA, "na.action")
  gValues <- raster::getValues(G.stack)
  g_noNA <- stats::na.omit(gValues)
  g_naIDs <- attr(g_noNA, "na.action")

  ids_raster <- 1:dim(gValues)[1]
  ids_raster <- ids_raster[- g_naIDs]
  m1 <- m_noNA
  m2 <- g_noNA

  if(dim(m1)[2] != dim(m2)[2]) {
    stop("Stacks must have the same dimensions.")
  }

  out_index <- plot_out(mValues, gValues)

  steps <- seq(1, dim(m2)[1], comp.each)
  kkk <- c(steps,  dim(m2)[1] + 1)
  long_k <- length(kkk)

  if (parallel == FALSE) {
    mop1 <- lapply(1:(length(kkk) - 1), function(x) {
      seq_rdist <- kkk[x]:(kkk[x + 1] - 1)
      eudist <- fields::rdist(m2[seq_rdist, ], m1)
      mean_quantile <- lapply(1:dim(eudist)[1], function(y) {
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
    future::plan(future::multiprocess)
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
          pond_mean <- mean(di, na.rm = TRUE)
          return(pond_mean)
        })
        mop <- unlist(mop_dist)
        return(mop)
      }
      avance <- (x / long_k) * 100
      cat("Computation progress: ", avance,"%" ,"\n")
    }

    mop_list <- as.list(mop_env)
    mop_names <- sort(as.numeric(as.character(names(mop_list))))
    mop_names <- as.character(mop_names)
    mop_vals <- unlist(mop_list[mop_names])

    future::plan(future::sequential)
  }

  mop_raster[ids_raster] <- mop_vals
  mop_max <- raster::cellStats(mop_raster, "max") * 1.05
  mop_raster[out_index] <- mop_max
  mop_raster <- 1 - (mop_raster / mop_max)
  return(mop_raster)
}
