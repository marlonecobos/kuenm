#' Prepare tables for hierarchical partitioning
#'
#' @description Helper function to prepare tables for hierarchical partitioning analyses
#' of the varianze comming from distinct sources in ecological niche models.
#'
#' @param model.names (character) vector of model names.
#' @param sp.name (character) species names as it is in the csv file used for creating
#' the model.
#' @param format (character) format of model raster files. Options are: "asc" or "tif"; default = "asc".
#' @param replicate.numbers (numeric) vector of replicate numbers.
#' @param parameters (character) vector of parameter settings names.
#' @param clim.models (character) vector of climatic model names.
#' @param emi.scenarios (character) vector of emission scenario names.
#' @param iterations (numeric) number of iterations to be performed; default = 100.
#' @param sample.size (numeric) number of pixels to be sampled per each model. Default = 100.
#' @param out.dir (character) name of the output directory where matrices will be written.
#'
#' @return All the tables needed to performing hierarchical partitioning analyses of the
#' varianze in ecological niche models.
#'
#' @details At least two of the following sources of variation must exist to perform
#' the analysis: \code{replicate.numbers}, \code{parameters}, \code{clim.models},
#' or \code{emi.scenarios}.
#'
#' @export

hierpart_tables <- function(model.names, sp.name, replicate.numbers, format = "asc", parameters, clim.models,
                            emi.scenarios, iterations = 100, sample.size = 100, out.dir) {

  # installing needed packages if required
  # pcakages <- c("raster", "rgdal")
  # req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
  # if (length(req_packages) > 0) {
  #  install.packages(req_packages, dependencies = TRUE)
  # }

  if (missing(replicate.numbers) & missing(parameters) &
      missing(clim.models) & missing(emi.scenarios)) {
    stop(paste("At least two of the following sources of variation must exist:\n",
               "replicates, parameter settings, climate models, or emision scenarios.", sep = ""))
  }
  if (sum(c(!missing(replicate.numbers), !missing(parameters),
            !missing(clim.models), !missing(emi.scenarios))) <= 1) {
    stop(paste("At least two of the following sources of variation must exist:\n",
               "replicates, parameter settings, climate models, or emision scenarios.", sep = ""))
  }

  for (h in 1:iterations) {
    if (!missing(replicate.numbers) & !missing(parameters)) {
      if (missing(clim.models) & missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))
        cat("Value,Replicates,Parameters\n")

        for (i in 1:length(parameters)) {
          pattern <- paste(".*", parameters[i], paste(".*", format, sep = ""), sep = "")

          repl <- gregexpr(pattern, model.names)
          repli <- regmatches(model.names, repl)
          replic <- unlist(repli)

          mod <- raster::stack(replic)
          mods <- na.omit(raster::getValues(mod))

          set.seed(seeds[h])
          val <- mods[sample(1:nrow(mods), sample.size), ]
          vals <- c(val)

          cat(paste(vals, rep(replicate.numbers, each = sample.size),
                    paste(i, "\n", sep = ""), sep = ","))
        }
        sink()
      }

      if (!missing(clim.models) & missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))

        cat("Value,Replicates,Parameters,GCMs\n")
        source <- clim.models

        for (i in 1:length(source)) {
          pattern <- paste(".*", source[i], paste(".*", format, sep = ""), sep = "")

          gepl <- gregexpr(pattern, model.names)
          gepli <- regmatches(model.names, gepl)
          geplic <- unlist(gepli)

          for (j in 1:length(parameters)) {
            pattern <- paste(".*", parameters[j], paste(".*", format, sep = ""), sep = "")

            repl <- gregexpr(pattern, geplic)
            repli <- regmatches(geplic, repl)
            replic <- unlist(repli)

            mod <- raster::stack(replic)
            mods <- na.omit(raster::getValues(mod))

            set.seed(seeds[h])
            val <- mods[sample(1:nrow(mods), sample.size), ]
            vals <- c(val)

            cat(paste(vals, rep(replicate.numbers, each = sample.size),
                      j, paste(i, "\n", sep = ""), sep = ","))
          }
        }
        sink()
      }

      if (missing(clim.models) & !missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))

        cat("Value,Replicates,Parameters,RCPs\n")
        source <- emi.scenarios

        for (i in 1:length(source)) {
          pattern <- paste(".*", source[i], paste(".*", format, sep = ""), sep = "")

          gepl <- gregexpr(pattern, model.names)
          gepli <- regmatches(model.names, gepl)
          geplic <- unlist(gepli)

          for (j in 1:length(parameters)) {
            pattern <- paste(".*", parameters[j], paste(".*", format, sep = ""), sep = "")

            repl <- gregexpr(pattern, geplic)
            repli <- regmatches(geplic, repl)
            replic <- unlist(repli)

            mod <- raster::stack(replic)
            mods <- na.omit(raster::getValues(mod))

            set.seed(seeds[h])
            val <- mods[sample(1:nrow(mods), sample.size), ]
            vals <- c(val)

            cat(paste(vals, rep(replicate.numbers, each = sample.size),
                      j, paste(i, "\n", sep = ""), sep = ","))
          }
        }
        sink()
      }

      if (!missing(clim.models) & !missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))

        cat("Value,Replicates,Parameters,GCMs,RCPs\n")

        for (i in 1:length(emi.scenarios)) {
          pattern <- paste(".*", emi.scenarios[i], paste(".*", format, sep = ""), sep = "")

          eepl <- gregexpr(pattern, model.names)
          eepli <- regmatches(model.names, eepl)
          eeplic <- unlist(eepli)

          for (j in 1:length(clim.models)) {
            pattern <- paste(".*", clim.models[j], paste(".*", format, sep = ""), sep = "")

            gepl <- gregexpr(pattern, eeplic)
            gepli <- regmatches(eeplic, gepl)
            geplic <- unlist(gepli)

            for (k in 1:length(parameters)) {
              pattern <- paste(".*", parameters[k], paste(".*", format, sep = ""), sep = "")

              repl <- gregexpr(pattern, geplic)
              repli <- regmatches(geplic, repl)
              replic <- unlist(repli)

              mod <- raster::stack(replic)
              mods <- na.omit(raster::getValues(mod))

              set.seed(seeds[h])
              val <- mods[sample(1:nrow(mods), sample.size), ]
              vals <- c(val)

              cat(paste(vals, rep(replicate.numbers, each = sample.size),
                        k, j, paste(i, "\n", sep = ""), sep = ","))
            }
          }
        }
        sink()
      }
    }

    if (!missing(replicate.numbers) & missing(parameters)) {
      if (!missing(clim.models) & missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))

        cat("Value,Replicates,GCMs\n")
        source <- clim.models

        for (i in 1:length(source)) {
          pattern <- paste(".*", source[i], paste(".*", format, sep = ""), sep = "")

          gepl <- gregexpr(pattern, model.names)
          gepli <- regmatches(model.names, gepl)
          geplic <- unlist(gepli)

          mod <- raster::stack(geplic)
          mods <- na.omit(raster::getValues(mod))

          set.seed(seeds[h])
          val <- mods[sample(1:nrow(mods), sample.size), ]
          vals <- c(val)

          cat(paste(vals, rep(replicate.numbers, each = sample.size),
                    paste(i, "\n", sep = ""), sep = ","))
        }
        sink()
      }

      if (missing(clim.models) & !missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))

        cat("Value,Replicates,RCPs\n")
        source <- emi.scenarios

        for (i in 1:length(source)) {
          pattern <- paste(".*", source[i], paste(".*", format, sep = ""), sep = "")

          gepl <- gregexpr(pattern, model.names)
          gepli <- regmatches(model.names, gepl)
          geplic <- unlist(gepli)

          mod <- raster::stack(geplic)
          mods <- na.omit(raster::getValues(mod))

          set.seed(seeds[h])
          val <- mods[sample(1:nrow(mods), sample.size), ]
          vals <- c(val)

          cat(paste(vals, rep(replicate.numbers, each = sample.size),
                    paste(i, "\n", sep = ""), sep = ","))
        }
        sink()
      }

      if (!missing(clim.models) & !missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))

        cat("Value,Replicates,GCMs,RCPs\n")

        for (i in 1:length(emi.scenarios)) {
          pattern <- paste(".*", emi.scenarios[i], paste(".*", format, sep = ""), sep = "")

          eepl <- gregexpr(pattern, model.names)
          eepli <- regmatches(model.names, eepl)
          eeplic <- unlist(eepli)

          for (j in 1:length(clim.models)) {
            pattern <- paste(".*", clim.models[j], paste(".*", format, sep = ""), sep = "")

            gepl <- gregexpr(pattern, eeplic)
            gepli <- regmatches(eeplic, gepl)
            geplic <- unlist(gepli)

            mod <- raster::stack(geplic)
            mods <- na.omit(raster::getValues(mod))

            set.seed(seeds[h])
            val <- mods[sample(1:nrow(mods), sample.size), ]
            vals <- c(val)

            cat(paste(vals, rep(replicate.numbers, each = sample.size),
                      j, paste(i, "\n", sep = ""), sep = ","))
          }
        }
        sink()
      }
    }

    if (missing(replicate.numbers) & !missing(parameters)) {
      if (!missing(clim.models) & missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))

        cat("Value,Parameters,GCMs\n")
        source <- clim.models

        for (i in 1:length(source)) {
          pattern <- paste(".*", source[i], paste(".*", format, sep = ""), sep = "")

          gepl <- gregexpr(pattern, model.names)
          gepli <- regmatches(model.names, gepl)
          geplic <- unlist(gepli)

          mod <- raster::stack(geplic)
          mods <- na.omit(raster::getValues(mod))

          set.seed(seeds[h])
          val <- mods[sample(1:nrow(mods), sample.size), ]
          vals <- c(val)

         cat(paste(vals, rep(1:length(parameters), each = sample.size),
                    paste(i, "\n", sep = ""), sep = ","))
        }
        sink()
      }

      if (missing(clim.models) & !missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))

        cat("Value,Parameters,RCPs\n")
        source <- emi.scenarios

        for (i in 1:length(source)) {
          pattern <- paste(".*", source[i], paste(".*", format, sep = ""), sep = "")

          gepl <- gregexpr(pattern, model.names)
          gepli <- regmatches(model.names, gepl)
          geplic <- unlist(gepli)

          mod <- raster::stack(geplic)
          mods <- na.omit(raster::getValues(mod))

          set.seed(seeds[h])
          val <- mods[sample(1:nrow(mods), sample.size), ]
          vals <- c(val)

          cat(paste(vals, rep(1:length(parameters), each = sample.size),
                    paste(i, "\n", sep = ""), sep = ","))
        }
        sink()
      }

      if (!missing(clim.models) & !missing(emi.scenarios)) {
        set.seed(1)
        seeds <- sample(1:1000, size = iterations)
        sink(paste(out.dir, paste("table_hp_", h, ".csv", sep = ""), sep = "/"))

        cat("Value,Parameters,GCMs,RCPs\n")

        for (i in 1:length(emi.scenarios)) {
          pattern <- paste(".*", emi.scenarios[i], paste(".*", format, sep = ""), sep = "")

          eepl <- gregexpr(pattern, model.names)
          eepli <- regmatches(model.names, eepl)
          eeplic <- unlist(eepli)

          for (j in 1:length(clim.models)) {
            pattern <- paste(".*", clim.models[j], paste(".*", format, sep = ""), sep = "")

            gepl <- gregexpr(pattern, eeplic)
            gepli <- regmatches(eeplic, gepl)
            geplic <- unlist(gepli)

            mod <- raster::stack(geplic)
            mods <- na.omit(raster::getValues(mod))

            set.seed(seeds[h])
            val <- mods[sample(1:nrow(mods), sample.size), ]
            vals <- c(val)

            cat(paste(vals, rep(1:length(parameters), each = sample.size),
                      j, paste(i, "\n", sep = ""), sep = ","))
          }
        }
        sink()
      }
    }

    cat(paste("   ", h, "of", iterations, "tables prepared\n"))
  }
}
