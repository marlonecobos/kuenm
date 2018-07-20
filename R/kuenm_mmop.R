#' Extrapolation risk analysis for multiple comparisons
#'
#' @description kuenm_mmop calculates mobility-oriented parity (MOP) layers by
#' comparing environmental values between the calibration area and multiple areas or
#' scenarios to which ecological niche models are transferred.
#'
#' @param G.var.dir (character) if project is TRUE, name of the folder containing folders in wHich variables of
#' projection scenarios are placed.
#' @param M.var.dir (character) name of the forlder containing folders in which calibration environmental
#' datasets are placed.
#' @param sets.var (character) value or vector with the name(s) of the sets of variables
#' from G.var.dir and M.var.dir that are going to be compared to create the MOP(s).
#' @param out.mop (character) name of the folder to which MOP results will be written.
#' @param percent (numeric) percetage of values sampled from the calibration region to calculate the MOP.
#' @param comp.each (numeric) compute distance matrix for a each fixed number of rows (default 1000).
#' @param parallel (logical) option to be passed to the \code{kuenm_mop} function (for each independent
#' MOP analyses). If TRUE, calculations will be performed in parallel using the available cores of the
#' computer. This will demand more RAM and almost full use of the CPU; hence, its use is more
#' recommended in high-performance computers. Using this option will speed up the analyses.
#' Default = FALSE
#'
#' @return A folder containing one or multiple mobility-oriented parity raster layers depending on
#' how many projection areas or scenarios are considered. This results will be organized by the
#' different sets of variables chosen for creating final models. Values of 0 in resultant rasters
#' represent strict extrapolation.
#'
#' @details This function can be used after selection of parameters that produce the best
#' models (when chosen sets of variables are known), or after producing final models
#' with the function \code{\link{kuenm_mod}}. In a MOP layer, areas of strict extrapolation
#' are excluded and other values represent how similar areas or scenarios are to
#' environmental conditions in the calibration area. MOP is calculated following Owens et al.
#' (2013; \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}).

kuenm_mmop <- function(G.var.dir, M.var.dir, sets.var, out.mop,
                       percent = 10, comp.each = 1000, parallel = FALSE) {

  #MOP directory
  dir.create(out.mop)

  #Calculating MOP for each comparison set by set
  for (h in 1:length(sets.var)) {

    dirsm <- dir(M.var.dir, pattern = sets.var[h], full.names = TRUE)

    dirsg <- dir(G.var.dir, pattern = sets.var[h], full.names = TRUE)
    dirsg_in <- dir(dirsg, full.names = TRUE)
    namesg <- dir(dirsg)

    dir.create(paste(out.mop, sets.var[h], sep = "/"))

    dirs_mop <- paste(paste(out.mop, sets.var[h], "MOP", sep = "/"),
                      paste(percent, "%", sep = ""), namesg, sep = "_")


    m_var <- list.files(dirsm, pattern = "asc", full.names = TRUE) #listing vars in M
    m_vars <- raster::stack(m_var) #stacking the variables

    if(.Platform$OS.type == "unix") {
      pb <- txtProgressBar(min = 0, max = length(dirsg_in), style = 3)
    } else {
      pb <- winProgressBar(title = "Progress bar", min = 0, max = length(dirsg_in), width = 300) #progress bar
    }

    for(i in 1:length(dirsg_in)) {
      g_var <- list.files(dirsg_in[i], pattern = "asc",
                               full.names = TRUE) #listing var of different Gs
      g_vars <- raster::stack(g_var)

      #MOP calculation
      mop_res <- kuenm_mop(M.stack = m_vars, G.stack = g_vars, percent = percent,
                           comp.each = comp.each, parallel = parallel)

      #Writing results
      raster::writeRaster(mop_res, filename = paste(dirs_mop[i],".tif", sep = ""), format = "GTiff")

      Sys.sleep(0.1)
      if(.Platform$OS.type == "unix") {
        setTxtProgressBar(pb, i)
      } else {
        setWinProgressBar(pb, i, title = paste(round(i / length(dirsg_in) * 100, 2),
                                               paste("% of the process for", sets.var[h], "has finished")))
      }
    }

    if(.Platform$OS.type != "unix") {
      suppressMessages(close(pb))
    }
    cat("\n", paste(h, "of", length(sets.var), "processes", sep = " "), "\n")
  }
}
