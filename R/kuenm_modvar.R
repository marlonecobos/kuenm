#' Prediction variance coming from distinct sources in ENMs
#'
#' @description kuenm_modvar calculates the variance in model predictions
#' distinguishing the source from which this is coming. In this version
#' potential sources of variation are: replicates, parameterizations, general
#' circulation models (GCMs), and emission scenarios. The last two considered
#' only when projections in time are performed.
#'
#' @param sp.name (character) name of the species. This name must be the one
#' that appears as part of the raster file of each model replicate. If results
#' are from Maxent, this is the name that is in the first column of the csv
#' containing species occurrence data (species) but excluding spaces.
#' @param fmod.dir (character) the  name of the folder in which final models
#' are (i.e., the output folder after using the \code{\link{kuenm_mod}}) function.
#' @param is.swd (logical) whether model calibration and final models were
#' produced using SWD format.
#' @param replicated (logical) whether or not final models were created performing
#' replicates.
#' @param format (character) format of model raster files. Options are: "asc"
#' or "tif"; default = "asc".
#' @param project (logical) if TRUE, assumes that models were projected to other
#' scenarios. These scenarios can be current (projections in space), and/or
#' future or past (projections in time).
#' @param current (character) pattern to look for when defining which is the
#' scenario of current projection. If not defined variance maps will be produced
#' for the area of calibration, and if #' any of \code{time.periods},
#' \code{clim.models}, or \code{emi.scenarios} exist, variance maps will be
#' produced for these layers as well.
#' @param time.periods (character or numeric) pattern to be searched when
#' identifying models from distinct time projections. If not defined it is
#' assumed that one time period was considered.
#' @param emi.scenarios (character) pattern to be searched for identifying
#' distinct emission scenarios (e.g., RCP numbers). If not defined, it is
#' assumed that only one emission scenario was used. Therefore, this source of
#' variation will not be considered.
#' @param clim.models (character) names that identify climatic models used to
#' project ENMs. If not defined it is assumed that only one climate model was
#' used. Therefore, this source of variation will not be considered.
#' @param ext.type (character) valid if \code{project} = TRUE, vector of
#' pattern(s) to be searched in the folders inside \code{fmod.dir} that identify
#' the extrapolation type(s) of model projections. This pattern(s) need to be
#' clearly distinguishable from the rest of the name of the model folder name.
#' For instance, capital letter can be used to separate this pattern from the
#' rest of the folder name (e.g., "EC" will be the patter that denotes
#' extrapolation and clamping in the folder named "M_0.1_F_l_set1_EC").
#' @param split.length (numeric) limit number of models to be processed at the
#' time. Bigger numbers would demand more from the RAM. Default = 100.
#' @param out.dir (character) name of the output directory to be created in
#' which subdirectories containing raster layers of model variance will be
#' written. Default = "Variation_from_sources".
#'
#' @return
#' Folders named Variation or Variation_("pattern" depending on the ext.type)
#' containing subdirectories named according to where/when models were projected.
#' Inside this folder, raster layers of variance coming from distinct sources.
#' All results will be written inside \code{out.dir}.
#'
#' @details
#' If any of the potential sources of variation is equal to one (e.g., only one
#' parameter, or only one climate model), this source of variation will not be
#' considered.
#'
#' Users must be specific when defining the patterns that the function will
#' search for. This patterns must be part of the model (raster layer) names so
#' the function can locate each file without problems. This function uses this
#' system of work to avoid high demands of the RAM while perfomring these
#' analyses.
#'
#' @usage
#' kuenm_modvar(sp.name, fmod.dir, is.swd, replicated, format = "asc", project,
#'              current, time.periods, emi.scenarios, clim.models, ext.type,
#'              split.length = 100, out.dir = "Variation_from_sources")
#'
#' @export
#'
#' @examples
#' # Models should be ready before starting these analyses, for an example of
#' # how to create them see https://github.com/marlonecobos/kuenm
#'
#' # Arguments
#' sp_name <- "sp1"
#' fmod_dir <- "Final_Models"
#' is_swd <- FALSE
#' rep <- TRUE
#' format <- "asc"
#' project <- TRUE
#' curr <- "current"
#' emi_scenarios <- c("RCP4.5", "RCP8.5")
#' c_mods <- c("GCM1", "GCM2")
#' ext_type <- c("E", "EC", "NE")
#' split <- 100
#' out_dir2 <- "Variation_from_sources"
#'
#' kuenm_modvar(sp.name = sp_name, fmod.dir = fmod_dir, is.swd = is_swd,
#'              replicated = rep, format = format, project = project,
#'              current = curr, emi.scenarios = emi_scenarios,
#'              clim.models = c_mods, ext.type = ext_type, split.length = split,
#'              out.dir = out_dir2)

kuenm_modvar <- function(sp.name, fmod.dir, is.swd, replicated, format = "asc",
                         project, current, time.periods, emi.scenarios,
                         clim.models, ext.type, split.length = 100,
                         out.dir = "Variation_from_sources") {

  cat("Preparing data for starting analyses, please wait...\n")

  if (missing(sp.name)) {
    stop("Argument sp.name needs to be defined.")
  }
  if (missing(fmod.dir)) {
    stop("Argument fmod.dir needs to be defined.")
  }
  if (!dir.exists(fmod.dir)) {
    stop(paste(fmod.dir, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (length(list.dirs(fmod.dir, recursive = FALSE)) == 0) {
    stop(paste(fmod.dir, "does not contain any subdirectory with sets of projection variables;",
               "\neach subdirectory inside", fmod.dir, "must containg at least one subdirectory",
               "\nwith the projection variables"))
  }
  if (missing(replicated)) {
    stop("Argument replicated needs to be provided. See fucntion's help for details.")
  }
  if (missing(project)) {
    stop("Argument project needs to be provided. See fucntion's help for details.")
  }
  if (project == TRUE) {
    if (missing(current)) {
      cat("Argument current is not defined, no current projection will be assumed.")
    }
    if (missing(time.periods)) {
      cat("Argument time.periods is not defined, an only time period will be assumed.")
    }
    if (missing(emi.scenarios)) {
      cat("Argument emi.scenarios is not defined, an only emission scenario will be assumed.")
    }
    if (missing(clim.models)) {
      cat("Argument clim.models is not defined, an only cimatic model will be assumed.")
    }
    if (missing(ext.type)) {
      stop("Argument ext.type needs to be provided. See fucntion's help for details.")
    }
  }
  if (missing(is.swd)) {
    stop("Argument is.swd is not defined.")
  }

  # All asc files
  a <- list.files(fmod.dir, pattern = paste0(".", format, "$"),
                  full.names = TRUE, recursive = TRUE)

  # Species name
  # Not statistics
  if (replicated == TRUE) {
    sn <- paste0(".*", sp.name, "_\\d.*")
  } else {
    sn <- paste0(".*", sp.name, ".*")
  }
  stn <- gregexpr(sn, a)
  stan <- regmatches(a, stn)
  statn <- unlist(stan)

  # No clamp no mess
  nc <- ".*clamping.*"
  ncl <- gregexpr(nc, statn)
  ncla <- regmatches(statn, ncl)
  nclam <- unlist(ncla)

  nm <- ".*novel.*"
  nme <- gregexpr(nm, statn)
  nmes <- regmatches(statn, nme)
  nmess <- unlist(nmes)

  # Final
  nstas <- statn[!statn %in% nclam & !statn %in% nmess]

  # Replicates
  if (replicated == TRUE) {
    nr <- paste0(sp.name, "_\\d")
    nre <- gregexpr(nr, nstas)
    nrep <- regmatches(nstas, nre)
    nrepl <- unlist(nrep)
    nrepli <- as.numeric(gsub(paste0(sp.name, "_"), "", unique(nrepl)))
  } else {
    nrepli <- 0
  }

  # Parameters
  p <- dir(fmod.dir)

  pa <- ".*_"
  para <- gregexpr(pa, p)
  param <- regmatches(p, para)
  parame <- unique(unlist(param))
  paramet <- gsub("_$", "", parame)

  # Folder to save all results
  dir.create(out.dir)

  #####
  # No projection
  if (project == FALSE) {
    ## Folder to save partial results
    dir.create(paste(out.dir, "Variation", sep = "/"))
    in_folder <- paste(out.dir, "Variation/Calibration_area", sep = "/")
    dir.create(in_folder)

    ## Calibration area
    nump <- ifelse(replicated == TRUE, "_\\d", "")
    if (is.swd == TRUE) {
      ca <- paste0(".*", sp.name, nump, ".asc")
      b <- list.files(fmod.dir, pattern = ca, full.names = T, recursive = T)
      cal <- gregexpr(ca, b)
      cali <- regmatches(b, cal)
    } else {
      ca <- paste0(".*", sp.name, nump, ".", format)
      cal <- gregexpr(ca, nstas)
      cali <- regmatches(nstas, cal)
    }

    calib <- unlist(cali)

    ### Var from replicates
    if (length(nrepli) > 1) {
      cat("\n   Variation from replicates:\n")
      cal_var <- var_models(model.names = calib, is.swd = is.swd, format = format,
                            sp.name = sp.name, source.codes = nrepli,
                            source = "replicates", split.length = split.length)

      if (is.swd == TRUE) {
        write.csv(cal_var, paste0(in_folder, "/Cal_var_replicates.csv"),
                  row.names = FALSE)
      } else {
        raster::writeRaster(cal_var,
                            filename = paste0(in_folder, "/Cal_var_replicates.tif"),
                            format = "GTiff")
      }
    }

    ### Var from parameters
    if (length(paramet) > 1) {
      cat("\n   Variation from parameters:\n")
      cal_var <- var_models(model.names = calib, is.swd = is.swd, format = format,
                            sp.name = sp.name, source.codes = paramet,
                            source = "parameters", split.length = split.length)

      if (is.swd == TRUE) {
        write.csv(cal_var, paste0(in_folder, "/Cal_var_parameters.csv"),
                  row.names = FALSE)
      } else {
        raster::writeRaster(cal_var,
                            filename = paste0(in_folder, "/Cal_var_parameters.tif"),
                            format = "GTiff")
      }
    }

  } else {
    # Splitting data according to extrapolation type
    ext_types <- list()
    var_folders <- vector()

    for (i in 1:length(ext.type)) {
      ecl <- paste0(".*_", ext.type[i], "/.*")
      ecla <- gregexpr(ecl, nstas)
      eclam <- regmatches(nstas, ecla)
      ext_types[[i]] <- unlist(eclam)
      var_folders[i] <- paste(out.dir, paste0("Variation_", ext.type[i]), sep = "/")
    }

    for (i in 1:length(ext_types)) {
      ## Folder to save partial results
      dir.create(var_folders[i])

      ## Calibration area
      in_folder <- paste(var_folders[i], "Calibration_area", sep = "/")
      dir.create(in_folder)

      nump <- ifelse(replicated == TRUE, "_\\d", "")
      if (is.swd == TRUE) {
        ca <- paste0(".*", sp.name, nump, ".csv")
        b <- list.files(fmod.dir, pattern = ca, full.names = T, recursive = T)
        cal <- gregexpr(gsub(".*", ".*_E/", ca, fixed = T), b)
        cali <- regmatches(b, cal)
      } else {
        ca <- paste0(".*", sp.name, nump, paste0(".", format))
        cal <- gregexpr(ca, ext_types[[i]])
        cali <- regmatches(ext_types[[i]], cal)
      }

      calib <- unlist(cali)

      ### Var from replicates
      if (length(nrepli) > 1) {
        cat("\n   Calibration area, variation from replicates:\n")
        cal_var <- var_models(model.names = calib, is.swd = is.swd, format = format,
                              sp.name = sp.name, source.codes = nrepli,
                              source = "replicates", split.length = split.length)

        if (is.swd == TRUE) {
          write.csv(cal_var, paste0(in_folder, "/Cal_var_replicates.csv"),
                    row.names = FALSE)
        } else {
          raster::writeRaster(cal_var,
                              filename = paste0(in_folder, "/Cal_var_replicates.tif"),
                              format = "GTiff")
        }
      }

      ### Var from parameters
      if (length(paramet) > 1) {
        cat("\n   Calibration area, variation from parameters:\n")
        cal_var <- var_models(model.names = calib, is.swd = is.swd, format = format,
                              sp.name = sp.name, source.codes = paramet,
                              source = "parameters", split.length = split.length)

        if (is.swd == TRUE) {
          write.csv(cal_var, paste0(in_folder, "/Cal_var_parameters.csv"),
                    row.names = FALSE)
        } else {
          raster::writeRaster(cal_var,
                              filename = paste0(in_folder, "/Cal_var_parameters.tif"),
                              format = "GTiff")
        }
      }

      # Current projections
      if (!missing(current)) {
        ## Current projection
        in_folder <- paste(var_folders[i], "Current_proj_area", sep = "/")
        dir.create(in_folder)

        cur <- paste(".*", current, format, sep = ".")
        curr <- gregexpr(cur, ext_types[[i]])
        curre <- regmatches(ext_types[[i]], curr)
        currente <- unlist(curre)

        ### Var from replicates
        if (length(nrepli) > 1) {
          cat("\n   Projection area (Current), variation from replicates:\n")
          cal_var <- var_models(model.names = currente, is.swd = F, format = format,
                                sp.name = sp.name, source.codes = nrepli,
                                source = "replicates", split.length = split.length)

          raster::writeRaster(cal_var,
                              filename = paste0(in_folder, "/Cur_proj_var_replicates.tif"),
                              format = "GTiff")
        }

        ### Var from parameters
        if (length(paramet) > 1) {
          cat("\n   Projection area (Current), variation from parameters:\n")
          cal_var <- var_models(model.names = currente, is.swd = F, format = format,
                                sp.name = sp.name, source.codes = paramet,
                                source = "parameters", split.length = split.length)

          raster::writeRaster(cal_var,
                              filename = paste0(in_folder, "/Cur_proj_var_parameters.tif"),
                              format = "GTiff")
        }
      }

      # Projections in time
      if (!missing(time.periods) | !missing(emi.scenarios) | !missing(clim.models)) {
        ## Models for other time periods
        if (!missing(current)) {
          timep <- ext_types[[i]][!ext_types[[i]] %in% calib & !ext_types[[i]] %in% currente]
        } else {
          timep <- ext_types[[i]][!ext_types[[i]] %in% calib]
        }

        if (missing(time.periods)) {
          timeper <- list(timep)
          timpe <- 1
        } else {
          timeper <- list()
          for (j in 1:length(time.periods)) {
            tp <- paste(paste0(".*", sp.name, "_\\d\\D.*", time.periods[j]),
                        "*", paste0(format, "$"), sep = ".")
            tip <- gregexpr(tp, timep)
            timp <- regmatches(timep, tip)
            timeper[[j]] <- unlist(timp)
          }
        }

        for (j in 1:length(timeper)) {
          in_folder <- paste(var_folders[i], paste0("Time", timpe[j] ,"proj_area"), sep = "/")
          dir.create(in_folder)

          ### Var from replicates
          if (length(nrepli) > 1) {
            cat("\n   Projection area", paste0("(Time ", j, "),"),
                "variation from replicates:\n")
            cal_var <- var_models(model.names = timeper[[j]], is.swd = F, format = format,
                                  sp.name = sp.name, source.codes = nrepli,
                                  source = "replicates", split.length = split.length)

            fname <- paste0(in_folder, paste0("/Time_", timpe[j],
                                              "_var_replicates.tif"))
            raster::writeRaster(cal_var, filename = fname, format = "GTiff")
          }

          ### Var from parameters
          if (length(paramet) > 1) {
            cat("\n   Projection area", paste0("(Time ", j, "),"),
                "variation from parameters:\n")
            cal_var <- var_models(model.names = timeper[[j]], is.swd = F, format = format,
                                  sp.name = sp.name, source.codes = paramet,
                                  source = "parameters", split.length = split.length)

            fname <- paste0(in_folder, paste0("/Time_", timpe[j],
                                              "_var_parameters.tif"))
            raster::writeRaster(cal_var, filename = fname, format = "GTiff")
          }

          ### Var from climatic models
          if (!missing(clim.models)) {
            if (length(clim.models) > 1) {
              cat("\n   Projection area", paste0("(Time ", j, "),"),
                  "variation from climatic models:\n")
              cal_var <- var_models(model.names = timeper[[j]], is.swd = F, format = format,
                                    sp.name = sp.name, source.codes = clim.models,
                                    source = "clim_models", split.length = split.length)

              fname <- paste0(in_folder, paste0("/Time_", timpe[j],
                                                "_var_clim_models.tif"))
              raster::writeRaster(cal_var, filename = fname, format = "GTiff")
            }
          }

          ### Var from emission scenarios
          if (!missing(emi.scenarios)) {
            if (length(emi.scenarios) > 1) {
              cat("\n   Projection area", paste0("(Time ", j, "),"),
                  "variation from emission scenarios:\n")
              cal_var <- var_models(model.names = timeper[[j]], is.swd = F, format = format,
                                    sp.name = sp.name, source.codes = emi.scenarios,
                                    source = "emi_scenarios", split.length = split.length)

              fname <- paste0(in_folder, paste0("/Time_", timpe[j],
                                                "_var_emi_scenarios.tif"))
              raster::writeRaster(cal_var, filename = fname, format = "GTiff")
            }
          }
          cat(paste(" ", j, "of", length(timpe), "time periods\n"))
        }
      }
      cat(paste(i, "of", length(ext_types), "complete processes\n"))
    }
  }

  # writing desciption
  result_description(process = "kuenm_modvar", out.dir = out.dir)

  cat(paste("\nCheck your working directory:", getwd(), "\n", sep = "\t"))
}

#' Helper function to calculate raster layers of model variance
#'
#' @param model.names (character) vector of model names.
#' @param is.swd (logical) whether model calibration and final models were
#' produced using SWD format.
#' @param format (character) format of model raster files. Options are: "asc" or
#' "tif"; default = "asc".
#' @param sp.name (character) species names. This name must be the one that
#' appears as part of the raster file of each model replicate.
#' @param source.codes (character or numeric) vector of names or numbers that
#' will be part of the pattern that will be searched.
#' @param source (character) source of variation to be evaluated. Options are:
#' "replicates", "parameters", "clim_models", and "emi_scenarios".
#' @param split.length (numeric) limit number of models to be processed at the
#' time. Bigger numbers would demand more from the RAM. Default = 100.
#'
#' @usage
#' var_models(model.names, is.swd, format = "asc", sp.name, source.codes,
#'            source, split.length = 100)
#'
#' @export
#'
#' @return
#' RasterLayer or data.frame depending on \code{is.swd} argument.

var_models <- function(model.names, is.swd, format = "asc", sp.name, source.codes,
                       source, split.length = 100) {
  if (missing(model.names)) {
    stop("Argument model.names is not defined.")
  }
  if (missing(sp.name)) {
    stop("Argument sp.name is not defined.")
  }
  if (missing(source.codes)) {
    stop("Argument source.codes is not defined.")
  }
  if (missing(source)) {
    stop("Argument source is not defined.")
  }
  if (missing(is.swd)) {
    stop("Argument is.swd is not defined.")
  }

  format <- ifelse(is.swd == TRUE, "csv", format)

  means <- list()
  for (i in 1:length(source.codes)) {
    if (source == "replicates") {
      pattern <- paste0(".*", sp.name, "_", source.codes[i], ".*", format)
    }
    if (source %in% c("parameters", "clim_models", "emi_scenarios")) {
      pattern <- paste0(".*", source.codes[i], ".*", format)
    }

    repl <- gregexpr(pattern, model.names)
    repli <- regmatches(model.names, repl)
    replic <- unlist(repli)
    lre <- length(replic)

    if (is.swd == TRUE) {
      if (lre == 1) {
        means[[i]] <- read.csv(replic)
      } else {
        mods <- lapply(1:lre, function(x) {
          tk <- ifelse(x == 1, 1, 3)
          read.csv(replic[x])[, tk:3]
        })
        means[[i]] <- apply(do.call(cbind, mods), 1, mean)
      }
    } else {
      if (lre <= split.length) {
        mod <- raster::stack(replic)
        mods <- raster::getValues(mod)
        mod <- mod[[1]]

        if (lre == 1) {
          means[[i]] <- c(mods)
        } else {
          means[[i]] <- suppressWarnings(apply(mods, 1, mean))
        }
        mods <- 1
      } else {
        replic <- split(replic, ceiling(seq_along(replic) / (split.length / 2)))

        mean <- list()
        for (j in 1:length(replic)) {
          mod <- raster::stack(replic[[j]])
          mods <- raster::getValues(mod)
          mod <- mod[[1]]

          mean[[j]] <- suppressWarnings(apply(mods, 1, mean))

          mods <- 1
        }

        mean <- suppressWarnings(do.call(cbind, mean))
        means[[i]] <- suppressWarnings(apply(mean, 1, mean))
      }
    }
    cat("     Process", i, "of", length(source.codes), "\n")
  }

  means <- suppressWarnings(do.call(cbind, means))

  if (is.swd == TRUE) {
    mod <- data.frame(means[, 1:2], variance = apply(means[, -(1:2)], 1, var))
  } else {
    mod[] <- suppressWarnings(apply(means, 1, var))
  }
  return(mod)
}
