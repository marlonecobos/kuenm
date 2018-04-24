#' Extrapolation risk analysis for multiple comparisons
#'
#' @description ku.enm.mmop calculates Mobility-Oriented Parity layers by
#' comparing environmental values between the M area and multiple areas or
#' scenarios to which ecological niche models are transferred.
#'
#' @param dirG a raster stack of variables representing the M area.
#' @param dirM a raster stack of variables representing the G area, areas or scenarios to
#' which models are transferred.
#' @param sets.var (character) value or vector with the name(s) of the sets of variables
#' from dirG and dirM that are going to be compared to create the MOP(s).
#' @param out.mop (character) name of the folder in which MOP results will be written.
#' @param percent (numeric) percetage of values sampled from M to calculate the MOP.
#' @param normalized (logical) if true values of similarity are presented from 0 to 1,
#' default = TRUE.
#'
#' @return A folder containing one or multiple Mobility-Oriented Parity raster layers depending on
#' how many projection areas or scenarios are considered. This results will be organized by the
#' different sets of variables chosen for creating final models.
#'
#' @details This function can be used after the selection of parameters that produce the best
#' candidate models (when chosen sets of variables are known), or after producing final models
#' with the function \code{\link{ku.enm.mod}}. In a MOP layer, areas of extrict extrapolation
#' are excluded and other values represent how similar are other areas or scenarios to
#' environmental conditions in the M. MOP is calculated following Owens et al., 2013
#' \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}.

ku.enm.mmop <- function(dirG, dirM, sets.var, out.mop, percent = 10, normalized = TRUE) {

  #MOP directory
  dir.create(out.mop)

  #Calculating MOP for each comparison set by set
  for (h in 1:length(sets.var)) {

    dirsm <- dir(dirM, pattern = sets.var[h], full.names = TRUE)

    dirsg <- dir(dirG, pattern = sets.var[h], full.names = TRUE)
    dirsg_in <- dir(dirsg, full.names = TRUE)
    namesg <- dir(dirsg)

    dir.create(paste(out.mop, sets.var[h], sep = "/"))

    dirs_mop <- paste(paste(out.mop, sets.var[h], "MOP", sep = "/"),
                      paste(percent, "%", sep = ""), namesg, sep = "_")


    m_var <- list.files(dirsm, pattern = "asc", full.names = TRUE) #listing vars in M
    m_vars <- raster::stack(m_var) #stacking the variables

    pb <- winProgressBar(title = "Progress bar", min = 0, max = length(dirsg_in), width = 300) #progress bar

    for(i in 1:length(dirsg_in)) {
      Sys.sleep(0.1)
      setWinProgressBar(pb, i, title = paste(round(i / length(dirsg_in) * 100, 2),
                                             paste("% of the process for", sets.var[h], "has finished")))

      g_var <- list.files(dirsg_in[i], pattern = "asc",
                               full.names = TRUE) #listing var of different Gs
      g_vars <- raster::stack(g_var)

      #MOP calculation
      mop_res <- ku.enm.mop(M.stack = m_vars, G.stack = g_vars,
                            percent = percent, normalized = normalized)

      #Writing results
      raster::writeRaster(mop_res, filename = paste(dirs_mop[i],".tif", sep = ""), format = "GTiff")
    }
    suppressMessages(close(pb))
    cat("\n", paste(h, "of", length(sets.var), "processes", sep = " "), "\n")
  }
}
