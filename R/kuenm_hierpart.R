#' Hierarchical partition of the variance coming from distinct sources in ENMs
#'
#' @description kuenm_hierpart performs a hierarchical partitioning analysis of
#' the variance coming from distinct sources in ENMs. In this version potential
#' sources of variation are: replicates, parameterizations, general circulation
#' models (GCMs), and emission scenarios. The last two are considered only when
#' projections in time are performed. At least two of these sources of variation
#' must be present in results.
#'
#' @param sp.name (character) name of the species. This name must be the one
#' that appears as part of the raster file of each model replicate. If results
#' are from Maxent, this is the name that is in the first column of the csv
#' containing species occurrence data (species) but spaces replaced by "_".
#' @param fmod.dir (character) name of the folder where all models are (e.g.,
#' the output folder after using the \code{\link{kuenm_mod}}) function.
#' @param format (character) format of model raster files. Options are "ascii",
#' "GTiff", and "EHdr" = bil. Default = "ascii".
#' @param replicated (logical) whether or not models were created with
#' replicates.
#' @param project (logical) if TRUE, it is assumed that models were projected
#' to other scenarios. These scenarios can be current (projections in space),
#' and/or future or past (projections in time).
#' @param current (character) pattern to look for when defining which is the
#' raster file representing current projections. If NULL, results will be
#' produced for the area of calibration, and if any of \code{time.periods},
#' \code{clim.models}, or \code{emi.scenarios} is defined, results will be
#' be produced for these variance sources as well.
#' @param time.periods (character) pattern to be searched to identify model
#' projections to distinct time periods. If NULL, the default, it is assumed
#' that only one time period was considered.
#' @param emi.scenarios (character) pattern to be searched to identify
#' distinct emission scenarios (e.g., "recp45"). If NULL, the default, it is
#' assumed that only one emission scenario was used. Therefore, this source of
#' variation will not be considered.
#' @param clim.models (character) names that identify climatic models used for
#' project ENMs. If NULL, the default, it is assumed that only one climate model
#' was used. Therefore, this source of variation will not be considered.
#' @param ext.type (character) pattern(s) to be searched in the folders inside
#' \code{fmod.dir} that identify the extrapolation type(s) used in model
#' projections. This pattern(s) needs to be clearly distinguishable from the
#' other parts of the name of the folder name containing the model. For instance,
#' "EC" will be the patter that denotes extrapolation and clamping in the folder
#' named "M_0.1_F_l_set1_EC".
#' @param iterations (numeric) number of iterations to be performed in the
#' hierarchical partitioning analysis. Default = 100.
#' @param sample.size (numeric) number of pixels to be sampled per each model.
#' Default = 100. Increasing this number is recommended when the number of
#' models and the computer features allow it.
#' @param keep.tables (logical) if TRUE, tables that are written in
#' \code{out.dir} for each iteration of the hierarchical partitioning analyses
#' are kept. Default = FALSE.
#' @param factors.col a vector of colors for the bars to be plotted; if not
#' defined, a gray color palette is used.
#' @param out.dir (character) name of the output directory to be created
#' where results of the hierarchical partitioning analysis will be written.
#' @param verbose (logical) whether to print messages; default = TRUE.
#'
#' @return
#' The function returns a data.frame containing the summary of total effects of
#' factors on variance contained in the models. A plot of these effects is also
#' returned.
#'
#' Other results are written in \code{out.dir}. Folders named Variation or
#' HP_results_(EC, NE, and/or E, depending on \code{ext.type}) containing
#' csv files with the results of the hierarchical partitioning analyses an a
#' plot summarizing the total effects of the sources of variation on the
#' variance in the models.
#'
#' @details
#' If the length of any of the potential sources of variation is equal to one
#' (e.g., only one parameter, or only one climate model), this source of
#' variation will not be considered.
#'
#' Users must be specific when defining the patterns that the function will
#' search for. These patterns must be part of the raster file names of the
#' models so the function can locate each file without problems.
#'
#' Error whiskers in resulting plots represent the 95% Confidence Interval of
#' the mean. This interval is calculated using a bootstrap approach.
#'
#' @usage
#' kuenm_hierpart(sp.name, fmod.dir, format = "ascii", replicated, project,
#'                current = NULL, time.periods = NULL, emi.scenarios = NULL,
#'                clim.models = NULL, ext.type, iterations = 100,
#'                sample.size = 1000, keep.tables = FALSE,
#'                factors.col, out.dir)
#'
#' @export
#'
#' @importFrom graphics arrows
#' @importFrom stats sd
#' @importFrom grDevices jpeg dev.off
#'
#' @examples
#' # Models should be ready before starting these analyses, for an example of
#' how to create them see https://github.com/marlonecobos/kuenm
#'
#' # Here an example of how to use the function once the models are ready
#' \dontrun{
#' ## Arguments
#' sp_name <- "sp1"
#' fmod_dir <- "Final_Models"
#' rep <- TRUE
#' format <- "ascii"
#' project <- TRUE
#' curr <- "current"
#' emi_scenarios <- c("RCP4.5", "RCP8.5")
#' c_mods <- c("GCM1", "GCM2")
#' ext_type <- c("E", "EC")
#' iter <- 100
#' s_size <- 1000
#' out_dir3 <- "Hierarchical_partitioning"
#'
#' # Running the function
#' kuenm_hierpart(sp.name = sp_name, fmod.dir = fmod_dir, format = format,
#'                replicated = rep, project = project, current = curr,
#'                emi.scenarios = emi_scenarios, clim.models = c_mods,
#'                ext.type = ext_type, iterations = iter,
#'                sample.size = s_size, out.dir = out_dir3)
#' }

kuenm_hierpart <- function(sp.name, fmod.dir, format = "ascii", replicated,
                           project, current = NULL, time.periods = NULL,
                           emi.scenarios = NULL, clim.models = NULL, ext.type,
                           iterations = 100, sample.size = 1000,
                           keep.tables = FALSE, factors.col = NULL,
                           out.dir, verbose = TRUE) {

  if (missing(sp.name)) {
    stop("Argument 'sp.name' needs to be defined")
  }
  if (missing(fmod.dir)) {
    stop("Argument 'fmod.dir' needs to be defined")
  }
  if (!format %in% c("ascii", "GTiff", "EHdr")) {
    stop("Argument 'format' is not valid")
  }
  if (!dir.exists(fmod.dir)) {
    stop(fmod.dir, " does not exist in the working directory, check folder name")
  }
  if (length(list.dirs(fmod.dir, recursive = FALSE)) == 0) {
    stop(fmod.dir, " does not contain any subdirectories, see function's help")
  }
  if (missing(replicated)) {
    stop("Argument 'replicated' needs to be provided. See fucntion's help for details")
  }
  if (missing(project)) {
    stop("Argument 'project' needs to be provided. See fucntion's help for details")
  }
  if (project == TRUE) {
    if (missing(current)) {
      if (verbose == TRUE) {
        message("Argument 'current' is not defined, no current projection will be assumed")
      }
    }
    if (missing(time.periods)) {
      if (verbose == TRUE) {
        message("Argument 'time.periods' is not defined, an only time period will be assumed")
      }
    }
    if (missing(emi.scenarios)) {
      if (verbose == TRUE) {
        message("Argument 'emi.scenarios' is not defined, an only emission scenario will be assumed")
      }
    }
    if (missing(clim.models)) {
      if (verbose == TRUE) {
        message("Argument 'clim.models' is not defined, an only cimatic model will be assumed")
      }
    }
    if (missing(ext.type)) {
      stop("Argument 'ext.type' needs to be defined. See fucntion's help for details")
    }
  }
  if (missing(factors.col)) {
    if (verbose == TRUE) {
      message("Argument 'factors.col' is not defined, a grey color palette will be used for ploting bars")
    }
  }

  if (verbose == TRUE) {
    message("Preparing data to start analyses, please wait...")
  }

  # All asc files
  a <- list.files(fmod.dir, pattern = paste0(rformat_type(format), "$"),
                  full.names = TRUE, recursive = TRUE)

  # Species name
  # Not statistics
  sn <- paste0(".*", sp.name, "_\\d.*")
  stn <- gregexpr(sn, a)
  stan <- regmatches(a, stn)
  statn <- unlist(stan)

  # No clamp no mess
  nc <- paste0(".*clamping.*")
  ncl <- gregexpr(nc, statn)
  ncla <- regmatches(statn, ncl)
  nclam <- unlist(ncla)

  nm <- paste0(".*novel.*")
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

  # par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  #####
  # No projection
  if (project == FALSE) {
    ## Folder to save partial results
    in_folder <- paste0(out.dir, "/HP_results")
    dir.create(in_folder)

    ## Calibration area
    if (replicated == TRUE) {
      ca <- paste0(".*", sp.name, "_\\d", paste0(rformat_type(format)))
    } else {
      ca <- paste0(".*", sp.name, paste0(rformat_type(format)))
    }
    cal <- gregexpr(ca, nstas)
    cali <- regmatches(nstas, cal)
    calib <- unlist(cali)

    if (length(nrepli) > 1 & length(paramet) > 1) {
      if (verbose == TRUE) {
        message("  Data preparation:")
      }
      tab_folder <- paste0(in_folder, "/hierpart_tables")
      dir.create(tab_folder)

      hierpart_tables(model.names = calib, sp.name = sp.name, format = format,
                      replicate.numbers = nrepli, parameters = paramet,
                      iterations = iterations, sample.size = sample.size,
                      out.dir = tab_folder)

      if (verbose == TRUE) {
        message("  Hierarchical partitioning analyses:")
      }
      hres_folder <- paste0(in_folder, "/hierpart_results")
      dir.create(hres_folder)
      hp_res <- hierpart_analyses(tables.folder = tab_folder,
                                  out.dir = in_folder, kept = keep.tables)

      ## Mean and confidence limits
      hp_mean <- apply(hp_res, 2, mean)
      names(hp_mean) <- colnames(hp_res)

      hp_se <- apply(hp_res, 2, function(x) {
        sd(apply(matrix(sample(x, size = length(x) * 1000, replace = TRUE),
                        1000, length(x)), 1, mean))
      })
      hp_lcu <- hp_mean + (2 * hp_se)
      hp_lcl <- hp_mean - (2 * hp_se)

      ## Figure
      if (missing(factors.col)) {
        factors.col <- "#E6E6E6"
      }

      plot_hierpart(hp_mean, hp_lcu, hp_lcl, factors.col, "Calibration area",
                    stacked = FALSE)

      ## Saving the figure
      jpeg(paste0(out.dir, "/Hier_par_results.jpg"), width = 166,
           height = 166, units = "mm", res = 600) #image to be saved

      plot_hierpart(hp_mean, hp_lcu, hp_lcl, factors.col, "Calibration area",
                    stacked = FALSE)

      invisible(dev.off())
    } else {
      stop("None or only one source of variation is being considered.\n",
           "These analyses require more than one source to be analyzed.")
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
      var_folders[i] <- paste0(out.dir, paste0("/HP_results_", ext.type[i]))
    }

    for (i in 1:length(ext_types)) {
      ## Folder to save partial results
      dir.create(var_folders[i])

      ## Calibration area
      ca <- paste0(".*", sp.name, "_\\d", paste0(".", rformat_type(format)))
      cal <- gregexpr(ca, ext_types[[i]])
      cali <- regmatches(ext_types[[i]], cal)
      calib <- unlist(cali)

      if (length(nrepli) > 1 & length(paramet) > 1) {
        if (verbose == TRUE) {
          message("  Calibration area, data preparation:")
        }
        tab_folder <- paste0(var_folders[i], "/Cal_hierpart_tables")
        dir.create(tab_folder)

        hierpart_tables(model.names = calib, sp.name = sp.name,
                        format = format, replicate.numbers = nrepli,
                        parameters = paramet, iterations = iterations,
                        sample.size = sample.size, out.dir = tab_folder)
        if (verbose == TRUE) {
          message("  Calibration area, hierarchical partitioning analyses:")
        }
        hres_folder <- paste0(var_folders[i], "/Cal_hierpart_results")
        dir.create(hres_folder)
        hp_res <- hierpart_analyses(tables.folder = tab_folder,
                                    out.dir = hres_folder, kept = keep.tables)

        ## Mean and confidence limits
        cal_hp_mean <- apply(hp_res, 2, mean)
        names(cal_hp_mean) <- colnames(hp_res)

        cal_hp_se <- apply(hp_res, 2, function(x) {
          sd(apply(matrix(sample(x, size = length(x) * 1000, replace = TRUE),
                          1000, length(x)), 1, mean))
        })
        cal_hp_lcu <- cal_hp_mean + (2 * cal_hp_se)
        cal_hp_lcl <- cal_hp_mean - (2 * cal_hp_se)
      } else {
        if (verbose == TRUE) {
          message("Calibration area:\n",
                  "None or only one source of variation is being considered.\n",
                  "These analyses require more than one source to be analyzed.")
        }
      }

      # Current projections
      if (!missing(current)) {
        ## Current projection
        cur <- paste(".*", current, rformat_type(format), sep = ".")
        curr <- gregexpr(cur, ext_types[[i]])
        curre <- regmatches(ext_types[[i]], curr)
        currente <- unlist(curre)

        if (length(nrepli) > 1 & length(paramet) > 1) {
          if (verbose == TRUE) {
            message("  Projection area (Current), data preparation:")
          }
          tab_folder <- paste0(var_folders[i], "/Curr_hierpart_tables")
          dir.create(tab_folder)

          hierpart_tables(model.names = currente, sp.name = sp.name,
                          format = format, replicate.numbers = nrepli,
                          parameters = paramet, iterations = iterations,
                          sample.size = sample.size, out.dir = tab_folder)
          if (verbose == TRUE) {
            message("  Projection area (Current), hierarchical partitioning analyses:")
          }
          hres_folder <- paste0(var_folders[i], "/Curr_hierpart_results")
          dir.create(hres_folder)
          hp_res <- hierpart_analyses(tables.folder = tab_folder,
                                      out.dir = hres_folder, kept = keep.tables)

          ## Mean and confidence limits
          cur_hp_mean <- apply(hp_res, 2, mean)
          names(cur_hp_mean) <- colnames(hp_res)

          cur_hp_se <- apply(hp_res, 2, function(x) {
            sd(apply(matrix(sample(x, size = length(x) * 1000, replace = TRUE),
                            1000, length(x)), 1, mean))
          })
          cur_hp_lcu <- cur_hp_mean + (2 * cur_hp_se)
          cur_hp_lcl <- cur_hp_mean - (2 * cur_hp_se)
        } else {
          if (verbose == TRUE) {
            message("Projection area (Current):\n",
                    "None or only one source of variation is being considered.\n",
                    "These analyses require more than one source to be analyzed.")
          }
        }
      }

      # Projections in time
      if (!missing(time.periods) | !missing(emi.scenarios) | !missing(clim.models)) {
        ## Models for other time periods
        if (!missing(current)) {
          timep <- ext_types[[i]][!ext_types[[i]] %in% calib &
                                    !ext_types[[i]] %in% currente]
        } else {
          timep <- ext_types[[i]][!ext_types[[i]] %in% calib]
        }

        if (missing(time.periods)) {
          timeper <- list(timep)
          timeps <- 1
        } else {
          timeps <- time.periods
          timeper <- list()
          for (j in 1:length(time.periods)) {
            tp <- paste(paste0(".*", sp.name, "_\\d\\D.*", time.periods[j]),
                        "*", paste0(rformat_type(format), "$"), sep = ".")
            tip <- gregexpr(tp, timep)
            timp <- regmatches(timep, tip)
            timeper[[j]] <- unlist(timp)
          }
        }

        time_hp_mean <- list()
        time_hp_lcu <- list()
        time_hp_lcl <- list()

        for (j in 1:length(timeper)) {
          if (verbose == TRUE) {
            message("  Projection area ", paste0("(Time ", timeps[j], ")"),
                    " data preparation:")
          }
          tab_folder <- paste0(var_folders[i], paste0("/Time_", timeps[j],
                                                      "_hierpart_tables"))
          dir.create(tab_folder)

          if (missing(nrepli) | length(nrepli) == 1) {
            nrepli <- ""
          }
          if (missing(paramet) | length(paramet) == 1) {
            paramet <- ""
          }
          if (missing(clim.models) | length(clim.models) == 1) {
            clim.models <- ""
          }
          if (missing(emi.scenarios) | length(emi.scenarios) == 1) {
            emi.scenarios <- ""
          }

          if (sum(c(length(nrepli) > 1, length(paramet) > 1, length(clim.models) > 1,
                    length(emi.scenarios) > 1)) >= 2) {

            if (length(nrepli) > 1 & length(paramet) > 1 & length(clim.models) > 1 &
                length(emi.scenarios) > 1) {
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name,
                              format = format, replicate.numbers = nrepli,
                              parameters = paramet, clim.models = clim.models,
                              emi.scenarios = emi.scenarios,
                              iterations = iterations, sample.size = sample.size,
                              out.dir = tab_folder)
            }

            if (length(nrepli) > 1 & length(paramet) > 1 & length(clim.models) > 1 &
                length(emi.scenarios) == 1) {
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name,
                              format = format, replicate.numbers = nrepli,
                              parameters = paramet, clim.models = clim.models,
                              iterations = iterations, sample.size = sample.size,
                              out.dir = tab_folder)
            }

            if (length(nrepli) > 1 & length(paramet) > 1 & length(clim.models) == 1 &
                length(emi.scenarios) > 1) {
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name,
                              format = format, replicate.numbers = nrepli,
                              parameters = paramet, emi.scenarios = emi.scenarios,
                              iterations = iterations, sample.size = sample.size,
                              out.dir = tab_folder)
            }

            if (length(nrepli) > 1 & length(paramet) == 1 & length(clim.models) > 1 &
                length(emi.scenarios) > 1) {
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name,
                              format = format, replicate.numbers = nrepli,
                              clim.models = clim.models, emi.scenarios = emi.scenarios,
                              iterations = iterations, sample.size = sample.size,
                              out.dir = tab_folder)
            }

            if (length(nrepli) == 1 & length(paramet) > 1 & length(clim.models) > 1 &
                length(emi.scenarios) > 1) {
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name,
                              format = format, parameters = paramet,
                              clim.models = clim.models, emi.scenarios = emi.scenarios,
                              iterations = iterations, sample.size = sample.size,
                              out.dir = tab_folder)
            }

            if (verbose == TRUE) {
              message("  Projection area ", paste0("(Time ", timeps[j], ")"),
                      " hierarchical partitioning analyses:")
            }
            hres_folder <- paste0(var_folders[i], paste0("/Time_", timeps[j],
                                                         "_hierpart_results"))
            dir.create(hres_folder)
            hp_res <- hierpart_analyses(tables.folder = tab_folder,
                                        out.dir = hres_folder, kept = keep.tables)

            ## Mean and confidence limits
            time_hp_mean[[j]] <- apply(hp_res, 2, mean)
            names(time_hp_mean[[j]]) <- colnames(hp_res)

            time_hp_se <- apply(hp_res, 2, function(x) {
              sd(apply(matrix(sample(x, size = length(x) * 1000, replace = TRUE),
                              1000, length(x)), 1, mean))
            })
            time_hp_lcu[[j]] <- time_hp_mean[[j]] + (2 * time_hp_se)
            time_hp_lcl[[j]] <- time_hp_mean[[j]] - (2 * time_hp_se)
          } else {
            if (verbose == TRUE) {
              message("Projection area ", paste0("(Time ", j, "):\n"),
                      "None or only one source of variation is being considered.\n",
                      "These analyses require more than one source to be analyzed.")
            }
          }
          if (verbose == TRUE) {
            message(" ", j, " of ", length(timeps), " time periods")
          }
        }
      }

      if (!exists("time_hp_mean") & !exists("cur_hp_mean") &
          !exists("cal_hp_mean")) {
        stop("None of the areas where models were projected to has enough sources of variation")
      }

      if (exists("time_hp_mean")) {
        ## Preparing data for figures
        time_hp_mean <- do.call(rbind, time_hp_mean)
        time_hp_lcu <- do.call(rbind, time_hp_lcu)
        time_hp_lcl <- do.call(rbind, time_hp_lcl)

        if (!exists("cur_hp_mean") & exists("cal_hp_mean")) {
          cal_hp_mean <- c(cal_hp_mean, rep(0, (dim(time_hp_mean)[2] - 2)))
          cal_hp_lcu <- c(cal_hp_lcu, rep(0, (dim(time_hp_mean)[2] - 2)))
          cal_hp_lcl <- c(cal_hp_lcl, rep(0, (dim(time_hp_mean)[2] - 2)))

          hp_mean <- rbind(cal_hp_mean, time_hp_mean)
          hp_lcu <- rbind(cal_hp_lcu, time_hp_lcu)
          hp_lcl <- rbind(cal_hp_lcl, time_hp_lcl)
          colnames(hp_mean) <- colnames(time_hp_mean)
          colnames(hp_lcu) <- colnames(time_hp_lcu)
          colnames(hp_lcl) <- colnames(time_hp_lcl)

          if (missing(time.periods)) {
            time_names <- paste("Projection in time", timeps)
            area_names <- c("Calibration area", time_names)
          } else {
            area_names <- c("Calibration area", "Projection in time")
          }
        }
        if (exists("cur_hp_mean") & exists("cal_hp_mean")) {
          cal_hp_mean <- c(cal_hp_mean, rep(0, (dim(time_hp_mean)[2] - 2)))
          cal_hp_lcu <- c(cal_hp_lcu, rep(0, (dim(time_hp_mean)[2] - 2)))
          cal_hp_lcl <- c(cal_hp_lcl, rep(0, (dim(time_hp_mean)[2] - 2)))

          cur_hp_mean <- c(cur_hp_mean, rep(0, (dim(time_hp_mean)[2] - 2)))
          cur_hp_lcu <- c(cur_hp_lcu, rep(0, (dim(time_hp_mean)[2] - 2)))
          cur_hp_lcl <- c(cur_hp_lcl, rep(0, (dim(time_hp_mean)[2] - 2)))

          hp_mean <- rbind(cal_hp_mean, cur_hp_mean, time_hp_mean)
          hp_lcu <- rbind(cal_hp_lcu, cur_hp_lcu, time_hp_lcu)
          hp_lcl <- rbind(cal_hp_lcl, cur_hp_lcl, time_hp_lcl)
          colnames(hp_mean) <- colnames(time_hp_mean)
          colnames(hp_lcu) <- colnames(time_hp_lcu)
          colnames(hp_lcl) <- colnames(time_hp_lcl)

          if (missing(time.periods)) {
            time_names <- paste("Projection in time", timeps)
            area_names <- c("Calibration area", "Current projection", time_names)
          } else {
            area_names <- c("Calibration area", "Current projection",
                            "Projection in time")
          }

        }

        if (exists("cur_hp_mean") & !exists("cal_hp_mean")) {
          cur_hp_mean <- c(cur_hp_mean, rep(0, (dim(time_hp_mean)[2] - 2)))
          cur_hp_lcu <- c(cur_hp_lcu, rep(0, (dim(time_hp_mean)[2] - 2)))
          cur_hp_lcl <- c(cur_hp_lcl, rep(0, (dim(time_hp_mean)[2] - 2)))

          hp_mean <- rbind(cur_hp_mean, time_hp_mean)
          hp_lcu <- rbind(cur_hp_lcu, time_hp_lcu)
          hp_lcl <- rbind(cur_hp_lcl, time_hp_lcl)
          colnames(hp_mean) <- colnames(time_hp_mean)
          colnames(hp_lcu) <- colnames(time_hp_lcu)
          colnames(hp_lcl) <- colnames(time_hp_lcl)

          if (missing(time.periods)) {
            time_names <- paste("Projection in time", timeps)
            area_names <- c("Current projection", time_names)
          } else {
            area_names <- c("Current projection", "Projection in time")
          }
        }

        if (!exists("cur_hp_mean") & !exists("cal_hp_mean")) {
          hp_mean <- time_hp_mean
          hp_lcu <-  time_hp_lcu
          hp_lcl <- time_hp_lcl
          colnames(hp_mean) <- colnames(time_hp_mean)
          colnames(hp_lcu) <- colnames(time_hp_lcu)
          colnames(hp_lcl) <- colnames(time_hp_lcl)

          if (missing(time.periods)) {
            time_names <- paste("Projection in time", timeps)
            area_names <- time_names
          } else {
            area_names <- "Projection in time"
          }
        }
      } else {
        if (!exists("cur_hp_mean") & exists("cal_hp_mean")) {
          hp_mean <- cal_hp_mean
          hp_lcu <- cal_hp_lcu
          hp_lcl <- cal_hp_lcl

          area_names <- c("Calibration area")
        }
        if (exists("cur_hp_mean") & exists("cal_hp_mean")) {
          hp_mean <- rbind(cal_hp_mean, cur_hp_mean)
          hp_lcu <- rbind(cal_hp_lcu, cur_hp_lcu)
          hp_lcl <- rbind(cal_hp_lcl, cur_hp_lcl)

          area_names <- c("Calibration area", "Current projection")
        }
        if (exists("cur_hp_mean") & !exists("cal_hp_mean")) {
          hp_mean <- cur_hp_mean
          hp_lcu <- cur_hp_lcu
          hp_lcl <- cur_hp_lcl

          area_names <- c("Current projection")
        }
        if (!exists("cur_hp_mean") & !exists("cal_hp_mean")) {
          stop(paste("None of the areas where the models were projected to has at least two sources of",
                     "\nvariation with contribution of more than one class."))
        }
      }

      ## Figure
      if (missing(factors.col)) {
        factors.col <- sort(gray.colors(dim(hp_mean)[1] + 1),
                            decreasing = TRUE)[1:dim(hp_mean)[1]]
      }

      suppressWarnings({
        plot_hierpart(hp_mean, hp_lcu, hp_lcl, factors.col, area_names,
                      stacked = TRUE)

        ## Saving the figure
        jpeg(paste0(var_folders[i], "/Hier_par_results.jpg"), width = 166,
             height = 166, units = "mm", res = 600) #image to be saved

        plot_hierpart(hp_mean, hp_lcu, hp_lcl, factors.col, area_names,
                      stacked = TRUE)

        invisible(dev.off())
      })

      if (!missing(current) | exists("cal_hp_mean")) {
        if (verbose == TRUE) {
          message("Warnings derive from current models having no effects for all sources")
        }
      }
      if (verbose == TRUE) {
        message(i, " of ", length(ext_types), " complete processes")
      }
    }
  }

  # writing description
  result_description(process = "kuenm_hierpart", out.dir = out.dir)

  if (verbose == TRUE) {
    message("Check your working directory:  ", getwd())
  }
}



# Helper to plot hierarchical partition results
plot_hierpart <- function(hp.mean, hp.ucl, hp.lcl, factors.col, area.names,
                          stacked = FALSE, length = 0.1) {
  if (missing(hp.mean)) {
    stop("Argument 'hp.mean' needs to be defined")
  }
  if (missing(hp.ucl)) {
    stop("Argument 'hp.ucl' needs to be defined")
  }
  if (missing(hp.lcl)) {
    stop("Argument 'hp.lcl' needs to be defined")
  }
  if (missing(factors.col)) {
    stop("Argument 'factors.col' needs to be defined")
  }
  if (missing(area.names)) {
    stop("Argument 'area.names' needs to be defined")
  }

  par(mar = c(4.5, 4.5, 0.5, 0.5), cex = 1.2)

  if (stacked == TRUE) {
    barx <- barplot(hp.mean, ylim = c(0, 100), border = "gray25",
                    main = NULL, xlab = "", ylab = "", las = 1,
                    beside = TRUE, col = factors.col)
  } else {
    barx <- barplot(hp.mean, ylim = c(0, 100), border = "gray25", space = 1.6,
                    main = NULL, xlab = "", ylab = "", las = 1,
                    col = factors.col)
  }
  arrows(barx, hp.ucl, barx, hp.lcl, angle = 90, code = 3, length = length)

  title(xlab = "Sources of variation", ylab = "Total effects (%)",
        cex.lab = 1.2)
  legend("topright", legend = area.names, pch = 22, col = "gray25",
         pt.bg = factors.col, bty = "n", cex = 0.9, inset = 0.02,
         title = "Means and 95% CI for")

  box(bty = "l")
}
