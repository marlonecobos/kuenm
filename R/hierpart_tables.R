#' Prepare tables for hierarchical partitioning analysis
#'
#' @description Helper function to prepare tables for hierarchical partitioning
#' analysis of the variance coming from distinct sources in ecological niche
#' models.
#'
#' @param model.names (character) vector of file names for model predictions
#' located in one or more directories.
#' @param sp.name (character) species name as in the csv file used for creating
#' models.
#' @param is.swd (logical) whether models were produced using SWD format.
#' @param format (character) format of model raster files. Options are "ascii",
#' "GTiff", and "EHdr" = bil. Default = "ascii".
#' @param replicate.numbers (numeric) vector of replicate numbers.
#' @param parameters (character) vector of parameter settings names.
#' @param clim.models (character) vector of climatic model names.
#' @param emi.scenarios (character) vector of emission scenario names.
#' @param iterations (numeric) number of iterations to be performed;
#' default = 100.
#' @param sample.size (numeric) number of pixels to be sampled per each model.
#' Default = 1000.
#' @param set.seed (numeric) initial seed to be set before running analysis.
#' @param out.dir (character) name of the output directory where matrices will
#' be written.
#' @param verbose (logical) whether to print messages; default = TRUE.
#'
#' @return
#' All the tables needed to performing hierarchical partitioning analyses of the
#' variance coming from distinct sources in ecological niche models. All results
#' are written in \code{out.dir}.
#'
#' @details
#' At least two of the following sources of variation must exist to perform
#' the analysis: \code{replicate.numbers}, \code{parameters}, \code{clim.models},
#' and \code{emi.scenarios}.
#'
#' @usage
#' hierpart_tables(model.names, sp.name, format = "ascii", replicate.numbers,
#'                 parameters, clim.models, emi.scenarios, iterations = 100,
#'                 sample.size = 1000, set.seed = 1, out.dir, verbose = TRUE)
#'
#' @importFrom raster stack
#' @importFrom stats na.omit
#' @export

hierpart_tables <- function(model.names, sp.name, is.swd, format = "ascii",
                            replicate.numbers, parameters, clim.models,
                            emi.scenarios, iterations = 100, sample.size = 1000,
                            set.seed = 1, out.dir, verbose = TRUE) {

  # initial tests
  if (missing(model.names)) {
    stop("Argument 'model.names' must be defined, see function's help")
  }
  if (missing(sp.name)) {
    stop("Argument 'sp.name' must be defined, see function's help")
  }
  if (missing(is.swd)) {
    stop("Argument 'is.swd' must be defined, see function's help")
  }
  if (!format %in% c("ascii", "GTiff", "EHdr")) {
    stop("Argument 'format' is not valid")
  }
  if (missing(replicate.numbers) & missing(parameters) &
      missing(clim.models) & missing(emi.scenarios)) {
    stop("At least two of the following sources of variation must exist:\n",
         "replicates, parameter settings, climate models, or emision scenarios")
  }
  if (sum(c(!missing(replicate.numbers), !missing(parameters),
            !missing(clim.models), !missing(emi.scenarios))) <= 1) {
    stop("At least two of the following sources of variation must exist:\n",
         "replicates, parameter settings, climate models, or emision scenarios")
  }
  if (missing(out.dir)) {
    stop("Argument 'out.dir' must be defined, see function's help")
  }

  # running
  ## format of rasters
  format <- rformat_type(format)

  ## seeds
  set.seed(set.seed)
  seeds <- sample(1:10000, size = iterations)

  ## loop for analysis
  for (h in 1:iterations) {
    filenam <- paste0(out.dir, "/", paste0("table_hp_", h, ".csv"))

    if (!missing(replicate.numbers) & !missing(parameters)) {
      if (missing(clim.models) & missing(emi.scenarios)) {

        cat("Value,Replicates,Parameters\n", file = filenam)

        for (i in 1:length(parameters)) {
          if (is.swd == FALSE) {
            pattern <- paste0(".*", parameters[i], paste0(".*", format, "$"))
          } else {
            pattern <- paste0(".*", parameters[i], ".*.csv$")
          }

          vals <- get_rvals(model.names, pattern, is.swd, sample.size, seeds[h])

          cat(paste(vals, rep(replicate.numbers, each = sample.size),
                    paste0(i, "\n"), sep = ","),
              append = TRUE, file = filenam)
        }
      }

      if (!missing(clim.models) & missing(emi.scenarios)) {
        cat("Value,Replicates,Parameters,GCMs\n", file = filenam)
        source <- clim.models

        for (i in 1:length(source)) {
          pattern <- paste0(".*", source[i], paste0(".*", format, "$"))

          gepl <- gregexpr(pattern, model.names)
          gepli <- regmatches(model.names, gepl)
          geplic <- unlist(gepli)

          for (j in 1:length(parameters)) {
            pattern <- paste0(".*", parameters[j], paste0(".*", format, "$"))

            vals <- get_rvals(geplic, pattern, FALSE, sample.size, seeds[h])

            cat(paste(vals, rep(replicate.numbers, each = sample.size),
                      j, paste0(i, "\n"), sep = ","),
                append = TRUE, file = filenam)
          }
        }
      }

      if (missing(clim.models) & !missing(emi.scenarios)) {
        cat("Value,Replicates,Parameters,RCPs\n", file = filenam)
        source <- emi.scenarios

        for (i in 1:length(source)) {
          pattern <- paste0(".*", source[i], paste0(".*", format, "$"))

          gepl <- gregexpr(pattern, model.names)
          gepli <- regmatches(model.names, gepl)
          geplic <- unlist(gepli)

          for (j in 1:length(parameters)) {
            pattern <- paste0(".*", parameters[j], paste0(".*", format, "$"))

            vals <- get_rvals(geplic, pattern, FALSE, sample.size, seeds[h])

            cat(paste(vals, rep(replicate.numbers, each = sample.size),
                      j, paste0(i, "\n"), sep = ","),
                append = TRUE, file = filenam)
          }
        }
      }

      if (!missing(clim.models) & !missing(emi.scenarios)) {
        cat("Value,Replicates,Parameters,GCMs,RCPs\n", file = filenam)

        for (i in 1:length(emi.scenarios)) {
          pattern <- paste0(".*", emi.scenarios[i], paste0(".*", format, "$"))

          eepl <- gregexpr(pattern, model.names)
          eepli <- regmatches(model.names, eepl)
          eeplic <- unlist(eepli)

          for (j in 1:length(clim.models)) {
            pattern <- paste0(".*", clim.models[j], paste0(".*", format, "$"))

            gepl <- gregexpr(pattern, eeplic)
            gepli <- regmatches(eeplic, gepl)
            geplic <- unlist(gepli)

            for (k in 1:length(parameters)) {
              pattern <- paste0(".*", parameters[k], paste0(".*", format, "$"))

              vals <- get_rvals(geplic, pattern, FALSE, sample.size, seeds[h])

              cat(paste(vals, rep(replicate.numbers, each = sample.size),
                        k, j, paste0(i, "\n"), sep = ","),
                  append = TRUE, file = filenam)
            }
          }
        }
      }
    }

    if (!missing(replicate.numbers) & missing(parameters)) {
      if (!missing(clim.models) & missing(emi.scenarios)) {
        cat("Value,Replicates,GCMs\n", file = filenam)
        source <- clim.models

        for (i in 1:length(source)) {
          pattern <- paste0(".*", source[i], paste0(".*", format, "$"))

          vals <- get_rvals(model.names, pattern, FALSE, sample.size, seeds[h])

          cat(paste(vals, rep(replicate.numbers, each = sample.size),
                    paste0(i, "\n"), sep = ","),
              append = TRUE, file = filenam)
        }
      }

      if (missing(clim.models) & !missing(emi.scenarios)) {
        cat("Value,Replicates,RCPs\n", file = filenam)
        source <- emi.scenarios

        for (i in 1:length(source)) {
          pattern <- paste0(".*", source[i], paste0(".*", format, "$"))

          vals <- get_rvals(model.names, pattern, FALSE, sample.size, seeds[h])

          cat(paste(vals, rep(replicate.numbers, each = sample.size),
                    paste0(i, "\n"), sep = ","),
              append = TRUE, file = filenam)
        }
      }

      if (!missing(clim.models) & !missing(emi.scenarios)) {
        cat("Value,Replicates,GCMs,RCPs\n", file = filenam)

        for (i in 1:length(emi.scenarios)) {
          pattern <- paste0(".*", emi.scenarios[i], paste0(".*", format, "$"))

          eepl <- gregexpr(pattern, model.names)
          eepli <- regmatches(model.names, eepl)
          eeplic <- unlist(eepli)

          for (j in 1:length(clim.models)) {
            pattern <- paste0(".*", clim.models[j], paste0(".*", format, "$"))

            vals <- get_rvals(eeplic, pattern, FALSE, sample.size, seeds[h])

            cat(paste(vals, rep(replicate.numbers, each = sample.size),
                      j, paste0(i, "\n"), sep = ","),
                append = TRUE, file = filenam)
          }
        }
      }
    }

    if (missing(replicate.numbers) & !missing(parameters)) {
      if (!missing(clim.models) & missing(emi.scenarios)) {
        cat("Value,Parameters,GCMs\n", file = filenam)
        source <- clim.models

        for (i in 1:length(source)) {
          pattern <- paste0(".*", source[i], paste0(".*", format, "$"))

          vals <- get_rvals(model.names, pattern, FALSE, sample.size, seeds[h])

          cat(paste(vals, rep(1:length(parameters), each = sample.size),
                    paste0(i, "\n"), sep = ","),
              append = TRUE, file = filenam)
        }
      }

      if (missing(clim.models) & !missing(emi.scenarios)) {
        cat("Value,Parameters,RCPs\n", file = filenam)
        source <- emi.scenarios

        for (i in 1:length(source)) {
          pattern <- paste0(".*", source[i], paste0(".*", format, "$"))

          vals <- get_rvals(model.names, pattern, FALSE, sample.size, seeds[h])

          cat(paste(vals, rep(1:length(parameters), each = sample.size),
                    paste0(i, "\n"), sep = ","),
              append = TRUE, file = filenam)
        }
      }

      if (!missing(clim.models) & !missing(emi.scenarios)) {
        cat("Value,Parameters,GCMs,RCPs\n", file = filenam)

        for (i in 1:length(emi.scenarios)) {
          pattern <- paste0(".*", emi.scenarios[i], paste0(".*", format, "$"))

          eepl <- gregexpr(pattern, model.names)
          eepli <- regmatches(model.names, eepl)
          eeplic <- unlist(eepli)

          for (j in 1:length(clim.models)) {
            pattern <- paste0(".*", clim.models[j], paste0(".*", format, "$"))

            vals <- get_rvals(eeplic, pattern, FALSE, sample.size, seeds[h])

            cat(paste(vals, rep(1:length(parameters), each = sample.size),
                      j, paste0(i, "\n"), sep = ","),
                append = TRUE, file = filenam)
          }
        }
      }
    }

    if (verbose == TRUE) {
      message("   ", h, " of ", iterations, " tables prepared")
    }
  }
}



# Helper to get values from raster files based on a pattern
get_rvals <- function(model.names, pattern, is.swd, sample.size = 1000,
                      set.seed = 1) {
  gepl <- gregexpr(pattern, model.names)
  gepli <- regmatches(model.names, gepl)
  geplic <- unlist(gepli)

  if (is.swd == FALSE) {
    mods <- raster::stack(geplic)
    mods <- na.omit(mods[])
  } else {
    mods <- lapply(geplic, function(x) {read.csv(x)[, 3]})
    mods <- do.call(cbind, mods)
    mods <- na.omit(mods)
  }

  set.seed(set.seed)
  if (nrow(mods) > sample.size) {
    mods <- mods[sample(1:nrow(mods), sample.size), ]
  }
  return(c(mods))
}
