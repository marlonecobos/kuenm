#' Hierarchical partition of the variance comming from distinct sources in ENMs
#'
#' @description kuenm_hierpart performs a hierarchical partitioning of the variance comming
#' from distinct sources in ENMs. In this version potential sources of variation are:
#' replicates, parameterizations, general circulation models (GCMs), and emission scenarios.
#' The last two considered only when projections in time are performed.
#'
#' @param sp.name (character) name of the species. This name must be the one that appears as part
#' of the raster file of each model repliate. If results are from Maxent, this is the name that
#' is in the first column of the csv containing species occurrence data (species) but ecxcluding spaces.
#' @param fmod.dir (character) the  name of the folder in which final models are (i.e., the output
#' folder after using the \code{\link{kuenm_mod}}) function.
#' @param replicated (logical) whether or not final models were created performing replicates.
#' @param format (character) format of model raster files. Options are: "asc" or "tif"; default = "asc".
#' @param project (logical) if TRUE, assumes that models were projected to other scenarios.
#' These scenarios can be current (projections in space), and/or future or past (projections in time).
#' @param current (character) pattern to look for when defining which is the scenario of current
#' projection. If not defined variance maps will be produced for the area of calibration, and if
#' any of \code{time.periods}, \code{clim.models}, or \code{emi.scenarios} exist, variace maps will
#' be produced for these layers as well.
#' @param time.periods (character or numeric) pattern to be searched when identifying models from
#' distinct time projections. If not defined it is assumed that one time period was considered.
#' @param emi.scenarios (character) pattern to be searched for identifying distinct emission
#' scenarios (e.g., RCP numbers). If not defined, it is
#' asumed that only one emission scenario was used. Therefore, this source of variation will not be
#' considered.
#' @param clim.models (character) names that identify climatic models used for project ENMs.
#' If not defined it is assumed that only one climate model was used. Therefore, this source of
#' variation will not be considered.
#' @param ext.type (character) vector of pattern(s) to be searched in the folders inside
#' \code{fmod.dir} that identify the extrapolation type(s) of model projections. This pattern(s)
#' need to be clearly distinguishable from the rest of the name of the model folder name. For instance,
#' capital letter can be used to separate this pattern from the rest of the folder name (e.g., "EC" will
#' be the patter that denotes extrapolation and clamping in the folder named "M_0.1_F_l_set1_EC").
#' @param iterations (numeric) number of iterations to be performed in the hierarchical partitioning
#' analysis. Default = 100.
#' @param sample.size (numeric) number of pixels to be sampled per each model. Default = 1000. Increasing
#' this number is recommended if each source of variation is not too numerous (i.e., few replicates, parameter
#' settings, climate models, or emission scenarios, are used when creating the models). Increasing this
#' number implies that more RAM will be needed.
#' @param keep.tables (logical) if TRUE, tables that are written in \code{out.dir} for each iteration
#' of the hirearchical partitioning analyses are kept. Default = FALSE.
#' @param factors.col a vector of colors for the bars; if not defined, a grey color palette is used.
#' @param out.dir (character) name of the output directory to be created in which subdirectories
#' (according to the \code{ext.type}) containing results of the hierarchical partitioning of the variance
#' of models will be written. Default = "Hierarchical_partitioning".
#'
#' @return Folders named Variation or HP_results_(EC, NE, and/or E, depending on the ext.type) conatining
#' csv files with the results of the hierarchichal partitioning analyses an a plot sumarizing the total
#' effects of the sourcen of variation on the total variance of the models. The plot will return to the
#' plot window as well. All results will be written inside \code{out.dir}.
#'
#' @details If any of the potential sources of variation is equal to one (e.g., only one parameter, or
#' only one climate model), this source of variation will not be considered.
#'
#' Error whiskers in resultant plots represent the 95% Confidence Interval of the mean. This interval is
#' calculated using a bootstrap approach.
#'
#' Users must be specific when defining the patterns that the function will search for. This patterns
#' must be part of the model (raster layer) names so the function can locate each file without problems.
#' This function uses this system of work to avoid high demands of the RAM while perfomring these analyses.
#'
#' @export
#'
#' @examples
#' # Models should be ready before starting these analyses, for an example of how to create them see
#' # https://github.com/marlonecobos/kuenm
#'
#' # Arguments
#' sp_name <- "sp1"
#' fmod_dir <- "Final_Models"
#' rep <- TRUE
#' format <- "asc"
#' project <- TRUE
#' curr <- "current"
#' emi_scenarios <- c("RCP4.5", "RCP8.5")
#' c_mods <- c("GCM1", "GCM2")
#' ext_type <- c("E", "EC", "NE")
#' iter <- 100
#' s_size <- 1000
#' out_dir3 <- "Hierarchical_partitioning"
#'
#' kuenm_hierpart(sp.name = sp_name, fmod.dir = fmod_dir, replicated = rep, format = format,
#'                project = project, current = curr, emi.scenarios = emi_scenarios,
#'                clim.models = c_mods, ext.type = ext_type, iterations = iter,
#'                sample.size = s_size, out.dir = out_dir3)

kuenm_hierpart <- function(sp.name, fmod.dir, replicated, format = "asc", project, current, time.periods, emi.scenarios,
                           clim.models, ext.type, iterations = 100, sample.size = 1000, keep.tables = FALSE,
                           factors.col, out.dir = "Hierarchical_partitioning") {

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
      cat("Argument current is not defined, no current projection will be assumed.\n")
    }
    if (missing(time.periods)) {
      cat("Argument time.periods is not defined, an only time period will be assumed.\n")
    }
    if (missing(emi.scenarios)) {
      cat("Argument emi.scenarios is not defined, an only emission scenario will be assumed.\n")
    }
    if (missing(clim.models)) {
      cat("Argument clim.models is not defined, an only cimatic model will be assumed.\n")
    }
    if (missing(ext.type)) {
      stop("Argument ext.type needs to be provided. See fucntion's help for details.\n")
    }
  }
  if (missing(factors.col)) {
    cat("Argument factors.col is not defined, a grey color palette will be used for ploting bars.")
  }

  # All asc files
  a <- list.files(fmod.dir, pattern = paste(".", format, "$", sep = ""),
                  full.names = TRUE, recursive = TRUE)

  # Species name
  # Not statistics
  sn <- paste(".*", sp.name, "_\\d.*", sep = "")
  stn <- gregexpr(sn, a)
  stan <- regmatches(a, stn)
  statn <- unlist(stan)

  # No clamp no mess
  nc <- paste(".*clamping.*", sep = "")
  ncl <- gregexpr(nc, statn)
  ncla <- regmatches(statn, ncl)
  nclam <- unlist(ncla)

  nm <- paste(".*novel.*", sep = "")
  nme <- gregexpr(nm, statn)
  nmes <- regmatches(statn, nme)
  nmess <- unlist(nmes)

  # Final
  nstas <- statn[!statn %in% nclam & !statn %in% nmess]

  # Replicates
  if (replicated == TRUE) {
    nr <- paste(sp.name, "_\\d", sep = "")
    nre <- gregexpr(nr, nstas)
    nrep <- regmatches(nstas, nre)
    nrepl <- unlist(nrep)
    nrepli <- as.numeric(gsub(paste(sp.name, "_", sep = ""), "", unique(nrepl))) #0:(length(unique(nrepl)) - 1)
  }else {
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

  # Funtion to graphic error whiskers in bars
  error_whiskers <- function(x, upper, lower, length = 0.1,...){
    if(length(x) !=length(lower) | length(lower) != length(upper))
      stop("vectors must have the same length")
    arrows(x, upper, x, lower, angle = 90, code = 3, length = length, ...)
  }

  #####
  # No projection
  if (project == FALSE) {
    ## Folder to save partial results
    in_folder <- paste(out.dir, "HP_results", sep = "/")
    dir.create(in_folder)

    ## Calibration area
    if (replicated == TRUE) {
      ca <- paste(".*", sp.name, "_\\d", paste(".", format, sep = ""), sep = "")
    }else {
      ca <- paste(".*", sp.name, paste(".", format, sep = ""), sep = "")
    }
    cal <- gregexpr(ca, nstas)
    cali <- regmatches(nstas, cal)
    calib <- unlist(cali)

    if (length(nrepli) > 1 & length(paramet) > 1) {
      cat("\n   Hierarchical partitioning table preparation:\n")
      tab_folder <- paste(in_folder, "hierpart_tables", sep = "/")
      dir.create(tab_folder)

      hierpart_tables(model.names = calib, sp.name = sp.name, replicate.numbers = nrepli, format = format,
                      parameters = paramet, iterations = iterations, sample.size = sample.size,
                      out.dir = tab_folder)

      cat("\n   Hierarchical partitioning analyses:\n")
      hres_folder <- paste(in_folder, "hierpart_results", sep = "/")
      dir.create(hres_folder)
      hp_res <- hierpart_analyses(tables.folder = tab_folder, out.dir = in_folder, kept = keep.tables)

      ## Mean and confidence limits
      hp_mean <- apply(hp_res, 2, mean)
      names(hp_mean) <- colnames(hp_res)

      hp_se <- apply(hp_res, 2, function(x){
        sd(apply(matrix(sample(x, size = length(x) * 1000, replace = TRUE), 1000, length(x)),
                 1, mean))
      })
      hp_lcu <- hp_mean + (2 * hp_se)
      hp_lcl <- hp_mean - (2 * hp_se)

      ## Figure
      if (missing(factors.col)) {
        factors.col <- "#E6E6E6"
      }

      par(mar = c(4.5,4.5,0.5,0.5), cex = 1.2)
      barx <- barplot(hp_mean, ylim = c(0, 100), border = "gray25", space = 1.6,
                      main = NULL, xlab = "", ylab = "", las = 1, col = factors.col)

      title(xlab = "Sources of variation", ylab = "Total effects (%)", cex.lab = 1.2)
      legend("topright", legend = "Means and 95% CIs for  ", bty = "n", cex = 0.9, inset = -0.01)
      legend("topright", legend = "Calibration area", pch = 22, col = "gray25", pt.bg = factors.col,
             bty = "n", cex = 0.9, inset = 0.03)

      error_whiskers(x = barx, upper = hp_lcu, lower = hp_lcl)
      box(bty = "l")

      ## Saving the figure
      jpeg(paste(out.dir, "Hier_par_results.jpg", sep = "/"), width = 166,
           height = 166, units = "mm", res = 600) #image to be saved

      par(mar = c(4.5,4.5,0.5,0.5), cex = 1.2)
      barx <- barplot(hp_mean, ylim = c(0, 100), border = "gray25", space = 1.6,
                      main = NULL, xlab = "", ylab = "", las = 1, col = factors.col)

      title(xlab = "Sources of variation", ylab = "Total effects (%)", cex.lab = 1.2)
      legend("topright", legend = "Means and 95% CIs for  ", bty = "n", cex = 0.9, inset = -0.01)
      legend("topright", legend = "Calibration area", pch = 22, col = "gray25", pt.bg = factors.col,
             bty = "n", cex = 0.9, inset = 0.03)

      error_whiskers(x = barx, upper = hp_lcu, lower = hp_lcl)
      box(bty = "l")

      invisible(dev.off())
    }else {
      stop(paste("None or only one source of variation is being considered.\n",
                 "These analyses require more than one source to be analyzed.", sep = ""))
    }
  }else {
    # Splitting data according to extrapolation type
    ext_types <- list()
    var_folders <- vector()

    for (i in 1:length(ext.type)) {
      ecl <- paste(".*_", ext.type[i], "/.*", sep = "")
      ecla <- gregexpr(ecl, nstas)
      eclam <- regmatches(nstas, ecla)
      ext_types[[i]] <- unlist(eclam)
      var_folders[i] <- paste(out.dir, paste("HP_results", ext.type[i], sep = "_"), sep = "/")
    }

    for (i in 1:length(ext_types)) {
      ## Folder to save partial results
      dir.create(var_folders[i])

      ## Calibration area
      ca <- paste(".*", sp.name, "_\\d", paste(".", format, sep = ""), sep = "")
      cal <- gregexpr(ca, ext_types[[i]])
      cali <- regmatches(ext_types[[i]], cal)
      calib <- unlist(cali)

      if (length(nrepli) > 1 & length(paramet) > 1) {
        cat("\n   Calibration area, hierarchical partitioning table preparation:\n")
        tab_folder <- paste(var_folders[i], "Cal_hierpart_tables", sep = "/")
        dir.create(tab_folder)

        hierpart_tables(model.names = calib, sp.name = sp.name, replicate.numbers = nrepli, format = format,
                        parameters = paramet, iterations = iterations, sample.size = sample.size,
                        out.dir = tab_folder)

        cat("\n   Calibration area, hierarchical partitioning analyses:\n")
        hres_folder <- paste(var_folders[i], "Cal_hierpart_results", sep = "/")
        dir.create(hres_folder)
        hp_res <- hierpart_analyses(tables.folder = tab_folder, out.dir = hres_folder, kept = keep.tables)

        ## Mean and confidence limits
        cal_hp_mean <- apply(hp_res, 2, mean)
        names(cal_hp_mean) <- colnames(hp_res)

        cal_hp_se <- apply(hp_res, 2, function(x){
          sd(apply(matrix(sample(x, size = length(x) * 1000, replace = TRUE), 1000, length(x)),
                   1, mean))
        })
        cal_hp_lcu <- cal_hp_mean + (2 * cal_hp_se)
        cal_hp_lcl <- cal_hp_mean - (2 * cal_hp_se)
      }else {
        cat(paste("\nCalibration area:\n",
                  "None or only one source of variation is being considered.\n",
                  "These analyses require more than one source to be analyzed.\n", sep = ""))
      }

      # Current projections
      if (!missing(current)) {
        ## Current projection
        cur <- paste(".*", current, format, sep = ".")
        curr <- gregexpr(cur, ext_types[[i]])
        curre <- regmatches(ext_types[[i]], curr)
        currente <- unlist(curre)

        if (length(nrepli) > 1 & length(paramet) > 1) {
          cat("\n   Projection area (Current), hierarchical partitioning table preparation:\n")
          tab_folder <- paste(var_folders[i], "Curr_hierpart_tables", sep = "/")
          dir.create(tab_folder)

          hierpart_tables(model.names = currente, sp.name = sp.name, replicate.numbers = nrepli, format = format,
                          parameters = paramet, iterations = iterations, sample.size = sample.size,
                          out.dir = tab_folder)

          cat("\n   Projection area (Current), hierarchical partitioning analyses:\n")
          hres_folder <- paste(var_folders[i], "Curr_hierpart_results", sep = "/")
          dir.create(hres_folder)
          hp_res <- hierpart_analyses(tables.folder = tab_folder, out.dir = hres_folder, kept = keep.tables)

          ## Mean and confidence limits
          cur_hp_mean <- apply(hp_res, 2, mean)
          names(cur_hp_mean) <- colnames(hp_res)

          cur_hp_se <- apply(hp_res, 2, function(x){
            sd(apply(matrix(sample(x, size = length(x) * 1000, replace = TRUE), 1000, length(x)),
                     1, mean))
          })
          cur_hp_lcu <- cur_hp_mean + (2 * cur_hp_se)
          cur_hp_lcl <- cur_hp_mean - (2 * cur_hp_se)
        }else {
          cat(paste("\nProjection area (Current):\n",
                    "None or only one source of variation is being considered.\n",
                    "These analyses require more than one source to be analyzed.\n", sep = ""))
        }
      }

      # Projections in time
      if (!missing(time.periods) | !missing(emi.scenarios) | !missing(clim.models)) {
        ## Models for other time periods
        if (!missing(current)) {
          timep <- ext_types[[i]][!ext_types[[i]] %in% calib & !ext_types[[i]] %in% currente]
        }else {
          timep <- ext_types[[i]][!ext_types[[i]] %in% calib]
        }

        if (missing(time.periods)) {
          timeper <- list(timep)
          timeps <- 1
        }else {
          timeps <- time.periods
          timeper <- list()
          for (j in 1:length(time.periods)) {
            tp <- paste(paste(".*", sp.name, "_\\d\\D.*", time.periods[j], sep = ""),
                        "*", paste(format, "$", sep = ""), sep = ".")
            tip <- gregexpr(tp, timep)
            timp <- regmatches(timep, tip)
            timeper[[j]] <- unlist(timp)
          }
        }

        time_hp_mean <- list()
        time_hp_lcu <- list()
        time_hp_lcl <- list()

        for (j in 1:length(timeper)) {
          cat("\n   Projection area", paste("(Time ", timeps[j], ")", sep = ""), "hierarchical partitioning table preparation:\n")
          tab_folder <- paste(var_folders[i], paste("Time", timeps[j] ,"hierpart_tables", sep = "_"), sep = "/")
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
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name, replicate.numbers = nrepli, format = format,
                              parameters = paramet, clim.models = clim.models, emi.scenarios = emi.scenarios,
                              iterations = iterations, sample.size = sample.size, out.dir = tab_folder)
            }

            if (length(nrepli) > 1 & length(paramet) > 1 & length(clim.models) > 1 &
                length(emi.scenarios) == 1) {
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name, replicate.numbers = nrepli, format = format,
                              parameters = paramet, clim.models = clim.models, iterations = iterations,
                              sample.size = sample.size, out.dir = tab_folder)
            }

            if (length(nrepli) > 1 & length(paramet) > 1 & length(clim.models) == 1 &
                length(emi.scenarios) > 1) {
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name, replicate.numbers = nrepli, format = format,
                              parameters = paramet, emi.scenarios = emi.scenarios,
                              iterations = iterations, sample.size = sample.size, out.dir = tab_folder)
            }

            if (length(nrepli) > 1 & length(paramet) == 1 & length(clim.models) > 1 &
                length(emi.scenarios) > 1) {
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name, replicate.numbers = nrepli, format = format,
                              clim.models = clim.models, emi.scenarios = emi.scenarios,
                              iterations = iterations, sample.size = sample.size, out.dir = tab_folder)
            }

            if (length(nrepli) == 1 & length(paramet) > 1 & length(clim.models) > 1 &
                length(emi.scenarios) > 1) {
              hierpart_tables(model.names = timeper[[j]], sp.name = sp.name, parameters = paramet, format = format,
                              clim.models = clim.models, emi.scenarios = emi.scenarios,
                              iterations = iterations, sample.size = sample.size, out.dir = tab_folder)
            }

            cat("\n   Projection area", paste("(Time ", timeps[j], ")", sep = ""), "hierarchical partitioning analyses:\n")
            hres_folder <- paste(var_folders[i], paste("Time", timeps[j] ,"hierpart_results", sep = "_"), sep = "/")
            dir.create(hres_folder)
            hp_res <- hierpart_analyses(tables.folder = tab_folder, out.dir = hres_folder, kept = keep.tables)

            ## Mean and confidence limits
            time_hp_mean[[j]] <- apply(hp_res, 2, mean)
            names(time_hp_mean[[j]]) <- colnames(hp_res)

            time_hp_se <- apply(hp_res, 2, function(x){
              sd(apply(matrix(sample(x, size = length(x) * 1000, replace = TRUE), 1000, length(x)),
                       1, mean))
            })
            time_hp_lcu[[j]] <- time_hp_mean[[j]] + (2 * time_hp_se)
            time_hp_lcl[[j]] <- time_hp_mean[[j]] - (2 * time_hp_se)
          }else {
            cat(paste("\nProjection area", paste("(Time ", j, "):\n,", sep = ""),
                      "None or only one source of variation is being considered.\n",
                      "These analyses require more than one source to be analyzed.\n", sep = ""))
          }

          cat(paste(" ", j, "of", length(timeps), "time periods\n"))
        }
      }

      if (!exists("time_hp_mean") & !exists("cur_hp_mean") & !exists("cal_hp_mean")) {
        stop(paste("None of the areas where the models were projected to has at least two sources of",
                   "\nvariation with contribution of more than one class."))
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
          }else {
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
          }else {
            area_names <- c("Calibration area", "Current projection", "Projection in time")
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
          }else {
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
          }else {
            area_names <- "Projection in time"
          }
        }
      }else {
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
        factors.col <- sort(gray.colors(dim(hp_mean)[1] + 1), decreasing = TRUE)[1:dim(hp_mean)[1]]
      }

      suppressWarnings({
        par(mar = c(4.5,4.5,0.5,0.5), cex = 1.2)
        barx <- barplot(hp_mean, ylim = c(0, 100), border = "gray25",
                        main = NULL, xlab = "", ylab = "", las = 1,
                        beside = TRUE, col = factors.col)

        title(xlab = "Sources of variation", ylab = "Total effects (%)", cex.lab = 1.2)
        legend("topright", legend = "Means and 95% CIs for  ", bty = "n", cex = 0.9,
               inset = -0.01)
        legend("topright", legend = area_names, pch = 22, col = "gray25", pt.bg = factors.col,
               bty = "n", cex = 0.9, inset = 0.03)

        box(bty = "l")
        error_whiskers(x = barx, upper = hp_lcu, lower = hp_lcl)

        ## Saving the figure

        jpeg(paste(var_folders[i], "Hier_par_results.jpg", sep = "/"), width = 166,
             height = 166, units = "mm", res = 600) #image to be saved

        par(mar = c(4.5,4.5,0.5,0.5), cex = 1.2)
        barx <- barplot(hp_mean, ylim = c(0, 100), border = "gray25",
                        main = NULL, xlab = "", ylab = "", las = 1,
                        beside = TRUE, col = factors.col)

        title(xlab = "Sources of variation", ylab = "Total effects (%)", cex.lab = 1.2)
        legend("topright", legend = "Means and 95% CIs for  ", bty = "n", cex = 0.9,
               inset = -0.01)
        legend("topright", legend = area_names, pch = 22, col = "gray25", pt.bg = factors.col,
               bty = "n", cex = 0.9, inset = 0.03)

        box(bty = "l")
        error_whiskers(x = barx, upper = hp_lcu, lower = hp_lcl)

        invisible(dev.off())
      })

      if (!missing(current) | exists("cal_hp_mean")) {
        cat("\nWarning messages are produced because current models have cero total \neffects for RCPs and GCMs.\n")
      }

      cat(paste(i, "of", length(ext_types), "complete processes\n"))
    }
  }
  cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
}
