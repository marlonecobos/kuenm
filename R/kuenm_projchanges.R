#' Detection of changes in models projected in time
#'
#' @description kuenm_projchanges performs map algebra operations to represent how and where
#' models projected in time change compared to the current one. If more than one climate model
#' (GCM) was used, it gives the degree of agreement among all the GCMs per emission scenario.
#'
#' @param occ (character) name of the csv file with all the occurrences used to create final models;
#' columns must be: species, longitude, latitude. In the case of the kuenm package, this must be the
#' name of the file that was used to create final models.
#' @param fmod.stats (character) the  name of the folder in which final models are (i.e., the output
#' folder after using the \code{\link{kuenm_modstats}}) function.
#' @param threshold (numeric) value from 0 to 100 that will be used as threshold, default = 5.
#' @param current (character) pattern to look for when defining which is the scenario of current
#' projection. If not defined coparisons will be performed between the calibration area and time
#' projections defined by the three following arguments.
#' @param time.periods (character or numeric) pattern to be searched when identifying models from
#' distinct time projections. If not defined it is assumed that one time period was considered.
#' @param emi.scenarios (character) pattern to be searched for identifying distinct emission
#' scenarios (e.g., RCP numbers). If not defined it is asumed that only one emission scenario
#' was used.
#' @param clim.models (character) names of that identify climatic models used for project ENMs.
#' If not defined it is assumed that only one climate model was used.
#' @param ext.type (character) vector of pattern(s) to be searched in the folders inside
#' \code{fmod.dir} that identify the extrapolation type(s) of model projections. This pattern(s)
#' need to be clearly distinguishable from the rest of the name of the model folder name. For instance,
#' capital letter can be used to separate this pattern from the rest of the folder name (e.g., "EC" will
#' be the patter that denotes extrapolation and clamping in the folder named "M_0.1_F_l_set1_EC").
#' @param out.dir (character) name of the output directory to be created in which subdirectories
#' containing the results of analyses of changes are. Default = "Projection_changes".
#'
#' @return Folders named Changes_("pattern" depending on the ext.type) containing raster layers
#' of the results, which include: changes in suitability, changes in suitable areas, and binary
#' raster layers of models for all scenarios. All results will be written inside \code{out.dir}.
#'
#' @details
#' If any of the potential sources of variation is equal to one (e.g., only one parameter, or
#' only one climate model), this source of variation will not be considered.
#'
#' Users must be specific when defining the patterns that the function will search for. This patterns
#' must be part of the model (raster layer) names so the function can locate each file without problems.
#' This function uses this system of work to avoid demand of the RAM while perfomring these analyses.
#'
#' @export
#'
#' @examples
#' # Models statistics should have been calculated before starting. This can be done using the
#' # kuenm_modstats function.
#'
#' # Arguments
#' occ <- "Sp_occ.csv"
#' fmod_stats <- "Final_Model_Stats"
#' thres <- 5
#' curr <- "current"
#' emi_scenarios <- c("RCP4.5", "RCP8.5")
#' c_mods <- c("GCM1", "GCM2")
#' ext_type <- c("E", "EC", "NE")
#' out_dir1 <- "Projection_Changes"
#'
#' kuenm_projchanges(occ = occ, fmod.stats = fmod_stats, threshold = thres,
#'                   current = curr, emi.scenarios = emi_scenarios,
#'                   clim.models = c_mods, ext.type = ext_type, out.dir = out_dir1)

kuenm_projchanges <- function(occ, fmod.stats, threshold = 5, current, time.periods, emi.scenarios,
                              clim.models, ext.type, out.dir = "Projection_changes") {

  # installing needed packages if required
  # pcakages <- c("raster", "rgdal")
  # req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
  # if (length(req_packages) > 0) {
  #  install.packages(req_packages, dependencies = TRUE)
  #}

  cat("Preparing data for starting analyses, please wait...\n")

  if (missing(occ)) {
    stop("Argument occ needs to be defined.")
  }
  if (missing(fmod.stats)) {
    stop("Argument fmod.dir needs to be defined.")
  }
  if (!dir.exists(fmod.stats)) {
    stop(paste(fmod.stats, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (length(list.dirs(fmod.stats, recursive = FALSE)) == 0) {
    stop(paste(fmod.stats, "does not contain any subdirectory with sets of projection variables;",
               "\neach subdirectory inside", fmod.stats, "must containg at least one subdirectory",
               "\nwith the projection variables"))
  }
  if (missing(current)) {
    cat("Argument current is not defined, no current projection will be assumed.\n")
  }
  if (missing(time.periods)) {
    cat("Argument time.periods is not defined, only one time projection will be assumed.\n")
  }
  if (missing(clim.models)) {
    cat("Argument clim.models is not defined, only one cimatic model will be assumed.\n")
  }
  if (missing(emi.scenarios)) {
    cat("Argument emi.scenarios is not defined, only one emission scenario will be assumed.\n")
  }
  if (missing(ext.type)) {
    stop("Argument ext.type needs to be provided. See fucntion's help for details.")
  }

  # Reading model names
  nstas <- list.files(fmod.stats, pattern = "med.tif$",
                  full.names = TRUE, recursive = TRUE)

  # Separating by extrapolation type
  ext_types <- list()
  var_folders <- vector()

  for (i in 1:length(ext.type)) {
    ecl <- paste(".*_", ext.type[i], "/.*", sep = "")
    ecla <- gregexpr(ecl, nstas)
    eclam <- regmatches(nstas, ecla)
    ext_types[[i]] <- unlist(eclam)
    var_folders[i] <- paste(out.dir, paste("Changes", ext.type[i], sep = "_"), sep = "/")
  }

  # Folder for all outputs
  dir.create(out.dir)

  # Occurrences
  occ <- read.csv(occ)[, 2:3]

  cat("\nStarting analyses, please wait...\n")

  for (i in 1:length(ext_types)) {
    # Folders per each ext_type
    dir.create(var_folders[i])

    if (missing(time.periods)) {
      time.periods <- ""
      timep <- 1
    }else {
      timep <- time.periods
    }

    # Current model
    if (missing(current)) {
      current <- "calibration"
    }
    cu <- paste(".*", current, ".*", sep = "")
    cur <- gregexpr(cu, ext_types[[i]])
    curr <- regmatches(ext_types[[i]], cur)
    curre <- unlist(curr)

    ca <- paste(".*calibration.*", sep = "")
    cal <- gregexpr(ca, ext_types[[i]])
    cali <- regmatches(ext_types[[i]], cal)
    calib <- unlist(cali)

    for (j in 1:length(time.periods)) {
      # Separating by times if exist
      tp <- paste(".*", time.periods[j], ".*", sep = "")
      tpe <- gregexpr(tp, ext_types[[i]])
      tper <- regmatches(ext_types[[i]], tpe)
      tperi <- unlist(tper)

      # Folders per each time
      dir.create(paste(var_folders[i],
                       paste("Period", timep[j], sep = "_"), sep = "/"))

      ## Separating by scenarios if exist
      if (missing(emi.scenarios)) {
        emi.scenarios <- ""
      }

      ## If exist and more than one, separate by emission scenarios
      for (k in 1:length(emi.scenarios)) {
        ### Separating by scenarios if exist
        es <- paste(".*", emi.scenarios[k], ".*", sep = "")
        esc <- gregexpr(es, tperi)
        esce <- regmatches(tperi, esc)
        escen <- unlist(esce)

        ### Folders per each time
        in_folder <- paste(var_folders[i], paste("Period", timep[j], sep = "_"),
                           paste("Scenario", emi.scenarios[k], sep = "_"), sep = "/")
        dir.create(in_folder)

        ### Diferences between current and other time climatic models, continuous models
        if (missing(clim.models)) {
          clim.models <- ""
        }

        comp_models <- model_changes(calibration.model = calib, current.model = curre,
                                     fclim.models = escen, result = "continuous")


        ### Writing files
        raster::writeRaster(comp_models, filename = paste(in_folder, "continuous_comparison.tif",
                                                  sep = "/"), format = "GTiff")

        ### Comparison of binary models between current and future
        comp_models <- model_changes(calibration.model = calib, current.model = curre,
                                     fclim.models = escen, result = "binary",
                                     threshold = threshold, occ = occ,
                                     clim.models = clim.models,
                                     out.dir = paste(in_folder, "Binary", sep = "/"))

        ### Writing files
        raster::writeRaster(comp_models, filename = paste(in_folder, "binary_comparison.tif",
                                                          sep = "/"), format = "GTiff")

        cat(paste("   ", k, "of", length(emi.scenarios), "emission scenarios\n"))
      }
      cat(paste(" ", j, "of", length(time.periods), "time periods\n"))
    }
    cat(paste(i, "of", length(ext_types), "complete processes\n"))
  }
  cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
}

#' Helper function to calculate model changes
#'
#' @param calibration.model (character) name of calibration model raster name.
#' @param current.model (character) name of current model raster name. It can be the same than
#' \code{calibration.model}.
#' @param fclim.models (character) vector of climatic model raster names.
#' @param result (character) type of result needed. options are "continuous" and "binary".
#' Default = "continuous".
#' @param threshold (numeric) if \code{result} = "binary", value from 0 to 100 that will be used
#' as threshold, default = 5.
#' @param occ if \code{result} = "binary", a numerical matrix containing coordinates of
#' the occurrence data used to create the final models; columns must be: longitude and latitude.
#' @param clim.models (character) names of that identify climatic models used for project ENMs.
#' If not defined it is assumed that only one climate model was used.
#' @param out.dir (character) name of the folder that will be created to save the binary models
#' if \code{result} = "binary".
#'
#' @export

model_changes <- function(calibration.model, current.model, fclim.models,
                          result = "continuous", threshold = 5, occ,
                          clim.models, out.dir) {

  # stack of current and time projections
  fclim.models <- sort(fclim.models) # just to organice in order
  cmodels <- raster::stack(c(current.model, fclim.models))

  if (result == "continuous" | result == "binary") {
    if (result == "continuous") {
      mods <- raster::getValues(cmodels)
      mod <- cmodels[[1]]

      # continuous model differences
      if (length(fclim.models) == 1) {
        mod[] <- apply(mods, 1, diff)
      }
      if (length(fclim.models) > 1) {
        comp <- cbind(mods[, 1], apply(mods[, 2:(length(fclim.models) + 1)], 1, median))
        mod[] <- apply(comp, 1, diff)
      }
    }

    if (result == "binary") {
      # Folder for binary outputs
      dir.create(paste(out.dir))

      # Threshold value calculation
      model <- raster::raster(calibration.model)
      o_suit <- raster::extract(model, occ)
      o_suit_sort <- sort(o_suit)
      thres <- o_suit_sort[ceiling(length(occ[, 1]) * threshold / 100) + 1]

      # Binarization
      ## Calibration area
      bin <- model
      raster::values(bin)[raster::values(bin) < thres] <- 0
      raster::values(bin)[raster::values(bin) >= thres] <- 1

      raster::writeRaster(bin, filename = paste(out.dir, "binary_calibration.tif",
                                                sep = "/"), format = "GTiff")

      ## Current with 0 and 1
      bin <- raster::raster(current.model)
      raster::values(bin)[raster::values(bin) < thres] <- 0
      raster::values(bin)[raster::values(bin) >= thres] <- 1

      raster::writeRaster(bin, filename = paste(out.dir, "binary_current.tif",
                                                sep = "/"), format = "GTiff")

      ### Current 0 and (number of clim models + 1)
      bin <- raster::raster(current.model)
      raster::values(bin)[raster::values(bin) < thres] <- 0
      raster::values(bin)[raster::values(bin) >= thres] <- (length(fclim.models) + 1)

      ### Climate models
      bins <- raster::stack(fclim.models)
      raster::values(bins)[raster::values(bins) < thres] <- 0
      raster::values(bins)[raster::values(bins) >= thres] <- 1

      raster::writeRaster(bins, filename = paste(out.dir, "binary.tif", sep = "/"),
                          format = "GTiff", bylayer = TRUE, suffix = sort(clim.models))

      ### Comparison
      all <- raster::stack(bin, bins)
      mods <- raster::getValues(all)
      mod <- all[[1]]

      mod[] <- apply(mods, 1, sum)
    }
  }else {
    stop("result option selected is not valid. Check the function's help.")
  }

  return(mod)
}
