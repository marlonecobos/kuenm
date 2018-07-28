#' Evaluation of final Maxent models with independent data
#'
#' @description kuenm_feval evaluates final Maxent models in terms of statistical
#' significance (partial ROC) and omission rates with a user-defined threshold (E).
#'
#' @param path (character) directory in which folders containig final models are being created or
#' were created.
#' @param occ.joint (character) the  csv file with training and testing occurrences combined,
#' or the file containing occurrences used to create final models; columns must be: species,
#' longitude, latitude.
#' @param occ.ind (character) the name of the csv file with independent occurrences for model
#' evaluation; these occurrences were not used when creating final models; columns as in occ.joint.
#' @param replicates (logical) whether or not final models were created performing replicates.
#' @param out.eval (character) name of the folder where evaluation results will be written.
#' @param threshold (numeric) the percentage of omission error allowed (E), default = 5.
#' @param rand.percent (numeric) the percentage of data to be used for the bootstrapping process
#' when calculating partial ROCs; default = 50.
#' @param iterations (numeric) the number of times that the bootstrap is going to be repeated;
#' default = 500.
#'
#' @return A list with two dataframes containing results from the evaluation process, and
#' a folder, in the working directory, containing a csv file with the final models evaluation
#' results.
#'
#' @details This function is used after or during the creation of final models.

kuenm_feval <- function(path, occ.joint, occ.ind, replicates, out.eval, threshold = 5,
                        rand.percent = 50, iterations = 500) {

  #Checking potential issues
  if (missing(path)) {
    stop(paste("Argument path is not defined, this is necessary for reading the",
               "\nfinal models created with the kuenm_mod function."))
  }
  if (!dir.exists(path)) {
    stop(paste(path, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (!file.exists(occ.joint)) {
    stop(paste(occ.joint, "does not exist in the working directory, check file name",
               "\nor extension, example: species_joint.csv"))
  }
  if (!file.exists(occ.ind)) {
    stop(paste(occ.ind, "does not exist in the working directory, check file name",
               "\nor extension, example: species_ind.csv"))
  }
  if (missing(out.eval)) {
    stop(paste("Argument out.eval is not defined, this is necessary for creating a folder",
               "\nwith the outputs of this function."))
  }
  if (missing(replicates)) {
    stop(paste("Logical argument replicates is not defined, this is necessary for",
               "\nselecting the layer that will be evaluated; it can be TRUE or FALSE."))
  }

  #####
  #Data
  ##Source of initial information for model evaluation
  fmodels <- dir(path)

  ###Model(s) for evaluation
  fm <- strsplit(fmodels, split = "_")

  fms <- vector()
  for (i in 1:length(fm)) {
    fms[i] <- paste(fm[[i]][-length(fm[[i]])], collapse = "_")
  }

  u_fmodels <- vector()
  for (i in 1:length(unique(fms))) {
    u_fmodels[i] <- dir(path, pattern = unique(fms)[i], full.names = TRUE)[1]
  }

  ###Names of the models to be evaluated
  mod_nam <- unique(fms)

  ##Joint set and independent occurrences
  occ <- read.csv(occ.joint) #read joint occurrences
  occ <- occ[, -1] #erase species name column

  occ1 <- read.csv(occ.ind) #read test occurrences
  occ1 <- occ1[, -1] #erase species name column

  ###ascii model name
  sp <- as.character(read.csv(occ.joint)[1, 1])
  sp <- gsub(" ", "_", sp)

  #####
  #pROCs and omission rates calculation
  cat("\nPartial ROCs and omission rates calculation, please wait...\n")

  proc_res <- list() #empty list of pROC values
  om_rates <- vector() #empty vector of omision rates

  if(.Platform$OS.type == "unix") {
    pb <- txtProgressBar(min = 0, max = length(u_fmodels), style = 3)
  } else {
    pb <- winProgressBar(title = "Progress bar", min = 0, max = length(u_fmodels),
                         width = 300) #progress bar
  }
  options(list(show.error.messages = FALSE, suppressWarnings = TRUE))

  for(i in 1:length(u_fmodels)) {
    Sys.sleep(0.1)
    if(.Platform$OS.type == "unix") {
      setTxtProgressBar(pb, i)
    } else {
      setWinProgressBar(pb, i, title = paste(round(i / length(u_fmodels) * 100, 2),
                                             "% of the evaluation process has finished"))
    }

    #Models to be evaluated
    if(replicates == TRUE) {
      mods1 <- list.files(u_fmodels[i], pattern = paste(sp, "median.asc", sep = "_"),
                          full.names = TRUE) #ascii models
      mod1 <- try(raster::raster(mods1), silent = TRUE)

      #partialROC calculation
      proc <- try(kuenm_proc(occ.test = occ1, model = mod1, threshold = threshold,
                             rand.percent = rand.percent, iterations = iterations),
                  silent = TRUE)
      proc_class <- class(proc)
      while (proc_class == "try-error") {
        mods1 <- list.files(u_fmodels[i], pattern = paste(sp, "median.asc", sep = "_"),
                            full.names = TRUE) #ascii models
        mod1 <- try(raster::raster(mods1), silent = TRUE)
        proc <- try(kuenm_proc(occ.test = occ1, model = mod1, threshold = threshold,
                               rand.percent = rand.percent, iterations = iterations),
                    silent = TRUE)
        proc_class <- class(proc)
        if(proc_class == "list") {
          break()
        }
      }
    }else {
      mods1 <- list.files(u_fmodels[i], pattern = paste(sp, ".asc", sep = "_"),
                          full.names = TRUE) #ascii models
      mod1 <- try(raster::raster(mods1), silent = TRUE)

      #partialROC calculation
      proc <- try(kuenm_proc(occ.test = occ1, model = mod1, threshold = threshold,
                             rand.percent = rand.percent, iterations = iterations),
                  silent = TRUE)
      proc_class <- class(proc)
      while (proc_class == "try-error") {
        mods1 <- list.files(u_fmodels[i], pattern = paste(sp, ".asc", sep = "_"),
                            full.names = TRUE) #ascii models
        mod1 <- try(raster::raster(mods1), silent = TRUE)
        proc <- try(kuenm_proc(occ.test = occ1, model = mod1, threshold = threshold,
                               rand.percent = rand.percent, iterations = iterations),
                    silent = TRUE)
        proc_class <- class(proc)
        if(proc_class == "list") {
          break()
        }
      }
    }

    #pROCs table
    proc_res[[i]] <- proc[[1]]

    #Omission rates calculation
    om_rates[i] <- kuenm_omrat(model = mod1, threshold = threshold,
                               occ.tra = occ, occ.test = occ1)

  }
  if(.Platform$OS.type != "unix") {
    suppressMessages(close(pb))
  }
  n.mod <- i

  ##Creating final tables
  ###From pROC analyses
  proc_res1 <- do.call(rbind, proc_res) #joining tables of the pROC results
  proc_res_m <- data.frame(mod_nam, proc_res1) #adding a new column with the number of AUC ratios interations < 1

  #####
  #Joining the results
  ku_enm_eval <- data.frame(proc_res_m, om_rates)
  colnames(ku_enm_eval) <- c("Model", "Mean_AUC_ratio", "Partial_ROC",#changing column names in the final table
                             paste("Omission_rate_at_", threshold, "%", sep = ""))

  #####
  #Statistics of the process
  ##Counting
  ku_enm_pROC <- ku_enm_eval[ku_enm_eval[, 3] <= threshold / 100, ]

  ku_enm_OR <- ku_enm_eval[ku_enm_eval[, 4] <= threshold / 100, ]

  ku_enm_pROC_OR <- ku_enm_pROC[ku_enm_pROC[, 4] <= threshold / 100, ]

  ##Preparing the table
  r_names <- c("All final models", "Statistically significant models",
               "Models meeting omission rate criteria", "Models meeting pROC and omission rate critera")
  statis <- c(length(ku_enm_eval[, 1]),
              length(ku_enm_pROC[, 3]),
              length(ku_enm_OR[, 4]),
              length(ku_enm_pROC_OR[, 4]))

  ku_enm_stats <- data.frame(r_names, statis)
  colnames(ku_enm_stats) <- c("Criteria", "Number of models")

  #####
  #Writing the results
  ##csv files
  cat("\nWriting kuenm_feval results...\n")
  dir.create(out.eval)

  name <- paste(out.eval, "fm_evaluation_results.csv", sep = "/")

  write.csv(ku_enm_eval, file = name, eol = "\n", na = "NA", row.names = FALSE)

  ##Retuning objects
  ###Dataframes in a list
  list_res <- list(ku_enm_stats, ku_enm_eval)
  names(list_res) <- c("Summary", "Evaluated models")

  #####
  #Finalizing the function
  cat("\nProcess finished\n")
  cat(paste("A folder containing results of the evaluation of", n.mod,
            "\nfinal models has been written\n", sep = " "))

  cat(paste("\nThe folder", out.eval, "contains a\n", sep = " "))
  cat("csv file with the final models' evaluation results.\n")
  cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
  cat("\n\nHere a summary of the evaluation process.\n\n")

  options(list(show.error.messages = TRUE, suppressWarnings = FALSE))

  print(ku_enm_stats)

  return(list_res)
}
