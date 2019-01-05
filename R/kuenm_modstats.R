#' Calculation of descriptive statistics of models
#'
#' @description kuenm_modstats calculates raster layers with some descriptive statistics of all
#' model replicates across multiple parameter settings. All of this, discriminating among models
#' transfered to distinct projection areas (scenarios).
#'
#' @param sp.name (character) name of the species. This name must be the one that appears as part
#' of the raster file of each model. If results are from Maxent, this is the name that
#' is in the first column of the csv containing species occurrence data (species) but changing spaces
#' (if there is any) by underscore.
#' @param fmod.dir (character) the  name of the folder in which final models are (e.g., the output
#' folder after using the \code{\link{kuenm_mod}}) function. It is important to have only the folders
#' containing the models in this directory. It can be only one folder or multiple folders containing models
#' for the same species, created with distinct parameter settings. If models were projected, and the
#' distinct types of extrapolation were used, the name of the folders contained in this directory should
#' include a pattern describing the type of extrapolation used (e.g., "EC" for extrapolation and
#' clamping in Maxent).
#' @param format (character) format of model raster files. Options are: "asc" or "tif"; default = "asc".
#' @param project (logical) if TRUE, it is assumed that models were projected to other scenarios.
#' @param statistics (character) vector of descriptive statistics to be calculated. Options include
#' med = median, mean, min = minimum, max = maximum, range, sd = stadard deviation, and se = standard error.
#' By default c("med", "min", "max", "range") are calculated, unless a character vector with the desired
#' statistics is provided.
#' @param replicated (logical) whether or not final models were created performing replicates.
#' @param proj.scenarios (character) valid if \code{project} = TRUE, vector of pattern(s) that identify
#' each projection area (scenario) to which models were projected.
#' @param ext.type (character) valid if \code{project} = TRUE, vector of pattern(s) to be searched in the
#' folders inside \code{fmod.dir} that identify the extrapolation type(s) of model projections. This pattern(s)
#' need to be clearly distinguishable from the rest of the name of the model folder name. For instance,
#' capital letter can be used to separate this pattern from the rest of the folder name (e.g., "EC" will be
#' the patter that denotes extrapolation and clamping in the folder named "M_0.1_F_l_set1_EC").
#' @param out.dir (character) name of the output directory to be created in which subdirectories
#' containing raster layers of model statistics will be written. Default = "Final_Model_Stats".
#'
#' @return Folders named Statistics or Statistics_("pattern" depending on the ext.type)
#' with all the raster layers of the descriptive statistics for models in \code{fmod.dir}.
#' Folders will be written inside \code{out.dir}.
#'
#' @details
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
#' format <- "asc"
#' project <- TRUE
#' stats <- c("med", "range")
#' rep <- TRUE
#' scenarios <- c("current", "GCM1_RCP4.5", "GCM1_RCP8.5", "GCM2_RCP4.5", "GCM2_RCP8.5")
#' ext_type <- c("E", "EC", "NE") # you can select only one type of extrapolation if needed
#' out_dir <- "Final_Model_Stats"
#'
#' kuenm_modstats(sp.name = sp_name, fmod.dir = fmod_dir, format = format, project = project,
#'                statistics = stats, replicated = rep, proj.scenarios = scenarios,
#'                ext.type = ext_type, out.dir = out_dir)

kuenm_modstats <- function(sp.name, fmod.dir, format = "asc", project, statistics, replicated,
                           proj.scenarios, ext.type, out.dir = "Final_Model_Stats") {

  # installing needed packages if required
  # pcakages <- c("raster", "rgdal")
  # req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
  # if (length(req_packages) > 0) {
  #  install.packages(req_packages, dependencies = TRUE)
  # }

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
    if (missing(proj.scenarios)) {
      stop("Argument proj.scenarios is not defined.")
    }
    if (missing(ext.type)) {
      stop("Argument ext.type needs to be provided. See fucntion's help for details.")
    }
  }

  # Folders depending on extrapolation settings
  if (project == FALSE) {
    parameters <- list(list.dirs(fmod.dir, recursive = FALSE))

    # Patterns to be found
    if (replicated == FALSE) {
      scenarios <- "calibration"
      m_names <- paste(sp.name, paste(".", format, "$", sep = ""), sep = "")
    }else {
      scenarios <- "calibration"
      m_names <- paste(sp.name, "_\\d", paste(".", format, "$", sep = ""), sep = "")
    }

    # Folders to save statistics
    res_folders <- paste(out.dir, "Statistics", sep = "/")

    # Name for results
    res_names <- "calibration"

  }else {
    parameters <- list()
    res_folders <- vector()

    for (i in 1:length(ext.type)) {
      parameters[[i]] <- dir(fmod.dir, pattern = ext.type[i], full.names = TRUE)
      res_folders[i] <- paste(out.dir, paste("Statistics", ext.type[i], sep = "_"), sep = "/")
    }

    # Patterns to be found
    if (replicated == FALSE) {
      scenarios <- c(sp.name, proj.scenarios)
      m_names <- paste(scenarios, paste(".", format, "$", sep = ""), sep = "")
    }else {
      sp.name1 <- paste(sp.name, "_\\d", sep = "")
      scenarios <- c(sp.name1, proj.scenarios)
      m_names <- paste(scenarios, paste(".", format, "$", sep = ""), sep = "")
    }

    # Name for results
    res_names <- c("calibration", proj.scenarios)
  }

  # Folder to save all results
  dir.create(out.dir)

  for (i in 1:length(res_folders)) {
    dir.create(res_folders[i])

    for (j in 1:length(scenarios)) {
      mod <- list()

      for (k in 1:length(parameters[[i]])) {
        mod[[k]] <- list.files(path = parameters[[i]][k], pattern = m_names[j], full.names = TRUE)
      }

      mod <- raster::stack(unlist(mod))
      n <- dim(mod)[3]
      mods <- raster::getValues(mod)
      mod <- mod[[1]]

      if (missing(statistics)) {
        statistics <- c("med", "min", "max", "range")
      }

      for (l in 1:length(statistics)) {
        if (statistics[l] == "med") {
          mod[] <- apply(mods, 1, median)
        }
        if (statistics[l] == "min") {
          mod[] <- apply(mods, 1, min)
        }
        if (statistics[l] == "max") {
          mod[] <- apply(mods, 1, max)
        }
        if (statistics[l] == "range") {
          mod[] <- apply(mods, 1, function(x) {max(x) - min(x)})
        }
        if (statistics[l] == "mean") {
          mod[] <- apply(mods, 1, mean)
        }
        if (statistics[l] == "sd") {
          mod[] <- apply(mods, 1, sd)
        }
        if (statistics[l] == "se") {
          mod[] <- apply(mods, 1, function(x) {sd(x) / sqrt(n)})
        }

        raster::writeRaster(mod, filename = paste(res_folders[i],
                                                  paste(res_names[j],
                                                        paste(statistics[l],
                                                              ".tif", sep = ""),
                                                        sep = "_"),
                                                  sep = "/"), format = "GTiff")

        cat(paste("     ", l, "of", length(statistics), "statistics\n"))
      }
      cat(paste("   ", j, "of", length(scenarios), "scenarios\n"))
    }
    cat(paste(i, "of", length(res_folders), "complete processes\n"))
  }
  cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
}
