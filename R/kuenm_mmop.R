#' Extrapolation risk analysis for multiple comparisons
#'
#' @description kuenm_mmop calculates mobility-oriented parity (MOP) layers by
#' comparing environmental values between the calibration area and multiple areas
#' or scenarios to which ecological niche models are transferred.
#'
#' @param G.var.dir (character) if project is TRUE, name of the directory containing
#' folders names as the sets to be compared (\code{sets.var}), which contain
#' subfolders with variables representing the scenarios of projection.
#' @param M.var.dir (character) name of the folder containing either folders with
#' ascii raster layers or csv files (SWD format) representing representing sets
#' of variables used to calibrate the models.
#' @param is.swd (logical) whether model calibration and final models were produced
#' using SWD format.
#' @param sets.var (character) value or vector with the name(s) of the sets of
#' variables from G.var.dir and M.var.dir that are going to be compared to create
#' the MOP(s).
#' @param out.mop (character) name of the folder to which MOP results will be
#' written.
#' @param percent (numeric) percentage of values sampled from the calibration
#' region to calculate the MOP.
#' @param comp.each (numeric) compute distance matrix for a each fixed number
#' of rows (default 2000).
#' @param parallel (logical) option to be passed to the \code{\link{kuenm_mop}}
#' function (for each independent MOP analyses). If TRUE, calculations will be
#' performed in parallel using \code{n.cores} of the computer. This will demand
#' more RAM and almost full use of the CPU; hence, its use is more recommended
#' in high-performance computers. Using this option will speed up the analyses.
#' Default = FALSE.
#' @param n.cores (numeric) number of cores to be used in parallel processing.
#' Default = NULL, in which case all CPU cores on current host - 1 will be used.
#'
#' @return
#' A folder containing one or multiple mobility-oriented parity raster layers
#' depending on how many projection areas or scenarios are considered. This
#' results will be organized by the different sets of variables chosen for
#' creating final models. Values of 0 in resultant raster layers represent strict
#' extrapolation.
#'
#' @details
#' This function can be used after selection of parameters that produce the best
#' models (when chosen sets of variables are known), or after producing final
#' models with the function \code{\link{kuenm_mod}}. In a MOP layer, areas of
#' strict extrapolation are excluded and other values represent how similar areas
#' or scenarios are to environmental conditions in the calibration area. MOP is
#' calculated following Owens et al. (2013; \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}).
#'
#' @usage
#' kuenm_mmop(G.var.dir, M.var.dir, is.swd, sets.var, out.mop, percent = 10,
#'            comp.each = 2000, parallel = FALSE)
#'
#' @export

kuenm_mmop <- function(G.var.dir, M.var.dir, is.swd, sets.var, out.mop, percent = 10,
                       comp.each = 2000, parallel = FALSE, n.cores = NULL) {
  if (missing(G.var.dir)) {
    stop("Argument G.var.dir is not defined.")
  }
  if (!dir.exists(G.var.dir)) {
    stop(paste(G.var.dir, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (length(list.dirs(G.var.dir, recursive = FALSE)) == 0) {
    stop(paste(G.var.dir, "does not contain any subdirectory with sets of projection variables;",
               "\neach subdirectory inside", G.var.dir, "must containg at least one subdirectory",
               "\nwith the projection variables"))
  }
  if (missing(M.var.dir)) {
    stop("Argument M.var.dir is not defined.")
  }
  if (!dir.exists(M.var.dir)) {
    stop(paste(M.var.dir, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (length(dir(M.var.dir)) == 0) {
    stop(paste(M.var.dir, "is empty. Check function's help for details."))
  }
  if (missing(is.swd)) {
    stop("Argument is.swd is not defined.")
  }

  #MOP directory
  dir.create(out.mop)

  #Calculating MOP for each comparison set by set
  for (h in 1:length(sets.var)) {

    if (is.swd == TRUE) {
      dirsm <- dir(M.var.dir, pattern = paste0("^", sets.var[h], ".csv$"),
                   full.names = TRUE)
      m_vars <- read.csv(dirsm)[, -(1:3)] #stacking the variables
    } else {
      dirsm <- dir(M.var.dir, pattern = paste0("^", sets.var[h], "$"),
                   full.names = TRUE)
      m_var <- list.files(dirsm, pattern = "asc", full.names = TRUE) #listing vars in M
      m_vars <- raster::stack(m_var) #stacking the variables
    }

    dirsg <- dir(G.var.dir, pattern = paste0("^", sets.var[h], "$"),
                 full.names = TRUE)
    dirsg_in <- dir(dirsg, full.names = TRUE)
    namesg <- dir(dirsg)

    dir.create(paste(out.mop, sets.var[h], sep = "/"))

    dirs_mop <- paste(paste(out.mop, sets.var[h], "MOP", sep = "/"),
                      paste(percent, "%", sep = ""), namesg, sep = "_")

    if(.Platform$OS.type != "unix") {
      pb <- winProgressBar(title = "Progress bar", min = 0, max = length(dirsg_in),
                           width = 300) #progress bar
    }

    for(i in 1:length(dirsg_in)) {
      Sys.sleep(0.1)
      if(.Platform$OS.type != "unix") {
        setWinProgressBar(pb, i, title = paste(round(i / length(dirsg_in) * 100, 2),
                                               paste("% of the process for", sets.var[h], "has finished")))
      }

      g_var <- list.files(dirsg_in[i], pattern = "asc",
                               full.names = TRUE) #listing var of different Gs
      g_vars <- raster::stack(g_var)

      #MOP calculation
      mop_res <- kuenm_mop(M.variables = m_vars, G.stack = g_vars, percent = percent,
                           comp.each = comp.each, parallel = parallel, n.cores = n.cores)

      #Writing results
      raster::writeRaster(mop_res, filename = paste(dirs_mop[i],".tif", sep = ""),
                          format = "GTiff")

      if(.Platform$OS.type == "unix") {
        cat("\n\t", paste(i, "of", length(dirsg_in), "MOPs", sep = " "), "\n")
      }
    }

    if(.Platform$OS.type != "unix") {
      suppressMessages(close(pb))
    }
    cat("\n", paste(h, "of", length(sets.var), "processes", sep = " "), "\n")
  }

  # writting desciption
  result_description(process = "kuenm_mmop", out.dir = out.mop)
}
