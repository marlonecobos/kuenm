#' Extrapolation risk analysis for single comparisons
#'
#' @description kuenm_mop calculates a mobility-oriented parity layer by
#' comparing environmental values between the calibration area and the area or
#' scenario to which an ecological niche model is transferred.
#'
#' @param M.variables a RasterStack of variables or a matrix with variables as columns
#' representing the calibration area. If matrix, columns must contain only
#' information for the variables to be used.
#' @param G.stack a RasterStack of variables representing the full area of interest, and areas
#' or scenarios to which models are transferred.
#' @param percent (numeric) percent of values sampled from te calibration region to calculate the MOP.
#' @param comp.each (numeric) compute distance matrix for a each fixed number of rows (default = 2000).
#' @param parallel (logical) if TRUE, calculations will be performed in parallel using \code{n.cores}
#' of the computer. This will demand more RAM and almost full use of the CPU; hence, its use
#' is more recommended in high-performance computers. Using this option will speed up the analyses.
#' Default = FALSE.
#' @param n.cores (numeric) number of cores to be used in parallel processing.
#' Default = NULL, in which case all CPU cores on current host - 1 will be used.
#'
#' @return A mobility-oriented parity RasterLayer where values of 0 represent strict extrapolation,
#' which means complete dissimilarity of environments between the calibration (M) or the background,
#' and the projection area (G).
#'
#' @details The MOP is calculated following Owens et al.
#' (2013; \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}). This function is a modification
#' of the \code{\link[ENMGadgets]{MOP}} funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
#'
#' @usage
#' kuenm_mop(M.variables, G.stack, percent = 10, comp.each = 2000,
#'           parallel = FALSE, n.cores = NULL)
#'
#' @export
#'
#' @examples
#' mvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                   pattern = "Mbio_", full.names = TRUE))
#' gvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                   pattern = "Gbio_", full.names = TRUE))
#' names(mvars) <- gsub("M", "", names(mvars))
#' names(gvars) <- names(mvars)
#'
#' perc <- 5
#'
#' mop <- kuenm_mop(M.variables = mvars, G.stack = gvars, percent = perc)


kuenm_mop <- function(M.variables, G.stack, percent = 10, comp.each = 2000,
                      parallel = FALSE, n.cores = NULL) {

  suppressPackageStartupMessages(library(doSNOW))
  suppressPackageStartupMessages(library(Kendall))
  suppressPackageStartupMessages(library(foreach))

  if (class(M.variables)[1] %in% c("RasterStack", "RasterBrick", "matrix", "data.frame")) {
    if (class(M.variables)[1] %in% c("RasterStack", "RasterBrick")) {
      mValues <- raster::getValues(M.variables)
    }
    if (class(M.variables)[1] == "data.frame") {mValues <- as.matrix(M.variables)}
    if (class(M.variables)[1] == "matrix") {mValues <- M.variables}
  } else {
    stop("Argument 'M.variables' is not valid.")
  }

  mop_raster <- G.stack[[1]]
  gValues <- raster::getValues(G.stack)

  mnames <- colnames(mValues)
  gnames <- colnames(gValues)
  gnames <- gnames[order(match(gnames, mnames))]
  gValues <- gValues[, gnames]

  if (!identical(mnames, gnames)) {
    stop("Variables in M and G must be the same.")
  }

  g_noNA <- stats::na.omit(gValues)
  g_naIDs <- attr(g_noNA, "na.action")
  m_noNA <- stats::na.omit(mValues)
  #m_naIDs <- attr(m_noNA, "na.action")

  ids_raster <- 1:dim(gValues)[1]
  if (!is.null(g_naIDs)) {
    ids_raster <- ids_raster[-g_naIDs]
  }

  out_index <- plot_out(mValues, gValues)

  steps <- seq(1, dim(g_noNA)[1], comp.each)
  kkk <- c(steps,  dim(g_noNA)[1] + 1)
  long_k <- length(kkk)

  if (parallel == FALSE) {
    pb <- txtProgressBar(min = 1, max = (length(kkk) - 1), style = 3)

    mop1 <- lapply(1:(length(kkk) - 1), function(x) {
      Sys.sleep(0.1)
      setTxtProgressBar(pb, x)

      seq_rdist <- kkk[x]:(kkk[x + 1] - 1)
      eudist <- fields::rdist(g_noNA[seq_rdist, ], m_noNA)
      mean_quantile <- lapply(1:dim(eudist)[1], function(y) {
        di <- eudist[y, ]
        qdi <- quantile(di, probs = percent / 100, na.rm = TRUE)
        ii <-  which(di <= qdi)
        return(mean(di[ii]))
      })

      return(unlist(mean_quantile))
    })

    close(pb)
    mop_vals <- unlist(mop1)

  }else {
    if (is.null(n.cores)) {n.cores <- parallel::detectCores() - 1}
    cl <- makeSOCKcluster(n.cores)
    registerDoSNOW(cl)

    pb <- txtProgressBar(min = 1, max = (length(kkk) - 1), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    mop_vals <- foreach(i = 1:(length(kkk) - 1), .packages = "Kendall", .inorder = TRUE,
                        .options.snow = opts, .combine = "c") %dopar% {
                          seq_rdist <- kkk[i]:(kkk[i + 1] - 1)
                          eudist <- fields::rdist(g_noNA[seq_rdist, ], m_noNA)
                          mean_quantile <- lapply(1:dim(eudist)[1], function(y) {
                            di <- eudist[y, ]
                            qdi <- quantile(di, probs = percent / 100, na.rm = TRUE)
                            ii <-  which(di <= qdi)
                            return(mean(di[ii]))
                          })
                          return(unlist(mean_quantile))
                        }

    close(pb)
    stopCluster(cl)
  }

  mop_raster[ids_raster] <- mop_vals
  mop_max <- raster::cellStats(mop_raster, "max") * 1.05
  mop_raster[out_index] <- mop_max
  mop_raster <- 1 - (mop_raster / mop_max)
  return(mop_raster)
}
