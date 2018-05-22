#' Evaluation of candidate Maxent models during calibration
#'
#' @description kuenm_ceval evaluates candidate models in terms of statistical
#' significance (partial ROC), prediction ability (omission rates), and model complexity (AICc).
#' After evaluation, this function selects the best models based on user-defined criteria.
#'
#' @param path (character) directory in which folders containig calibration models are being created
#' or were created.
#' @param occ.joint (character) the name of csv file with training and testing occurrences combined;
#' columns must be: species, longitude, latitude.
#' @param occ.tra (character) the name of the csv file with the training occurrences;
#' columns as in occ.joint.
#' @param occ.test (character) the name of the csv file with the evaluation occurrences;
#' columns as in occ.joint.
#' @param batch (character) name of the batch file with the code to create all candidate models
#' for calibration.
#' @param out.eval (character) name of the folder where evaluation results will be written.
#' @param threshold (numeric) the percentage of training data omission error allowed (E); default = 5.
#' @param rand.percent (numeric) the percentage of data to be used for the bootstraping process
#' when calculating partial ROCs; default = 50.
#' @param iterations (numeric) the number of times that the bootstrap is going to be repeated;
#' default = 500.
#' @param kept (logical) if FALSE, all candidate models will be erased after evaluation, default = TRUE.
#' @param selection (character) model selection criterion, can be "OR_AICc", "AICc", or "OR";
#' OR = omission rates. Default = "OR_AICc", which means that among models that are statistically significant
#' and that present omission rates below the threshold, those with delta AICc up to 2 will be
#' selected. If "AICc" criterion is chosen, all significant models with delta AICc up to 2 will be selected
#' If "OR" is chosen, the 10 first significant models with the lowest omission rates will be selected.
#'
#' @return  A list with three dataframes containing information about the calibration process. In addition,
#' a folder, in the working directory, containing a csv file with information about models meeting
#' the user-defined selection criterion, another csv file with a summary of the evaluation and selection
#' process, an extra csv file containing all the statistics of model performance (pROC, AICc, and
#' omission rates) for all candidate models, a png scatterplot of all models based on the AICc values and
#' omission rates, and an HTML file sumarizing all the information produced after evaluation for helping with
#' further interpretation.
#'
#' @details This function is used after or during the creation of Maxent candidate models for
#' calibration.

kuenm_ceval <- function(path, occ.joint, occ.tra, occ.test, batch, out.eval, threshold = 5,
                        rand.percent = 50, iterations = 500, kept = TRUE, selection = "OR_AICc") {
  #Checking potential issues
  if (!dir.exists(path)) {
    stop(paste(path, "does not exist in the working directory, check folder name",
               "\nor its existense."))
  }
  if (!file.exists(occ.joint)) {
    stop(paste(occ.joint, "does not exist in the working directory, check file name",
               "\nor extension, example: species_joint.csv"))
  }
  if (!file.exists(occ.tra)) {
    stop(paste(occ.tra, "does not exist in the working directory, check file name",
               "\nor extension, example: species_train.csv"))
  }
  if (!file.exists(occ.test)) {
    stop(paste(occ.test, "does not exist in the working directory, check file name",
               "\nor extension, example: species_test.csv"))
  }
  if (missing(out.eval)) {
    stop(paste("Argument out.eval is not defined, it is necessary for creating a folder",
               "\nwith the outputs of this function."))
  }
  if (missing(batch)) {
    warning(paste("Argument batch is not defined, the default name candidate_models",
                  "\nwill be used."))
    batch <- "candidate_models"
  }


  #####
  #Data
  ##Source of initial information for model evaluation order
  bat <- readLines(paste(batch, ".bat", sep = "")) #reading the batch file written to create the calibration models

  ###Recognizing the folders names and separating them for different procedures
  fol <- gregexpr("outputd.\\S*", bat)
  fold <- regmatches(bat, fol)
  folde <- unlist(fold)
  extract <- paste("outputdirectory=", path, "\\", sep = "")
  folder <- gsub(extract, "", folde, fixed = T) #names of all the calibration models folders

  folder_a <- gregexpr("\\S*all", folder)
  folder_al <- regmatches(folder, folder_a)
  folder_all <- unlist(folder_al) #folders with the models for calculating AICcs

  folder_c <- gregexpr("\\S*cal", folder)
  folder_ca <- regmatches(folder, folder_c)
  folder_cal <- unlist(folder_ca) #folder with the models for calculating pROCs and omission rates

  ##Models
  ###For AICc
  dir_names <- as.vector(paste(getwd(), "/", path, "/", folder_all, sep = "")) #vector of the subdirectories with the models

  ###For pROC and omission rates
  dir_names1 <- as.vector(paste(getwd(), "/", path, "/", folder_cal, sep = "")) #vector of the subdirectories with the models

  ###Names of the models to be evaluated
  mod_nam <- as.vector(gsub("_all", "", folder_all, fixed = TRUE)) #names of the models (taken from the folders names)

  ##Complete set and calibration and evaluation occurrences
  oc <- read.csv(occ.joint) #read all occurrences
  oc <- oc[, -1] #erase species name column

  occ <- read.csv(occ.tra) #read calibration occurrences
  occ <- occ[, -1] #erase species name column

  occ1 <- read.csv(occ.test) #read test occurrences
  occ1 <- occ1[, -1] #erase species name column

  #####
  #pROCs, omission rates, and AICcs calculation
  cat("\nPartial ROCs, omission rates, and AICcs calculation, please wait...\n")

  aiccs <- list() #empty list of AICc results
  proc_res <- list() #empty list of pROC values
  om_rates <- vector() #empty vector of omision rates

  if(.Platform$OS.type == "unix") {
    pb <- txtProgressBar(min = 0, max = length(dir_names), style = 3) #progress bar
  } else {
    pb <- winProgressBar(title = "Progress bar", min = 0, max = length(dir_names),
                         width = 300) #progress bar
  }

  for(i in 1:length(dir_names)) {
    Sys.sleep(0.1)
    if(.Platform$OS.type == "unix") {
      setTxtProgressBar(pb, i)
    } else {
      setWinProgressBar(pb, i, title = paste(round(i / length(dir_names) * 100, 2),
                                             "% of the evaluation process has finished"))
    }


    #AICc calculation
    suppressWarnings(while (!file.exists(as.vector(list.files(dir_names[i], pattern = ".lambdas",
                                                              full.names = TRUE)))) {
      Sys.sleep(1)
    })

    lbds <- as.vector(list.files(dir_names[i], pattern = ".lambdas",
                                 full.names = TRUE)) #lambdas file
    lambdas <- readLines(lbds)
    par_num <- n.par(lambdas) #getting the number of parameters for each model

    suppressWarnings(while (!file.exists(list.files(dir_names[i], pattern = "asc",
                                                    full.names = TRUE))) {
      Sys.sleep(1)
    })

    mods <- list.files(dir_names[i], pattern = "asc", full.names = TRUE) #name of ascii model
    mod <- raster::raster(mods) #reading each ascii model created with the complete set of occurrences
    aiccs[[i]] <- suppressWarnings(ENMeval::calc.aicc(nparam = par_num, occ = oc,
                                             predictive.maps = mod)) #calculating AICc for each model

    #pROCs calculation
    suppressWarnings(while (!file.exists(list.files(dir_names1[i], pattern = "asc",
                                                    full.names = TRUE))) {
      Sys.sleep(1)
    })

    mods1 <- list.files(dir_names1[i], pattern = "asc", full.names = TRUE) #ascii models
    mod1 <- raster::raster(mods1) #reading each ascii model created with the calibration occurrences

    proc <- kuenm_proc(occ.test = occ1, model = mod1, threshold = threshold,
                        rand.percent = rand.percent, iterations = iterations) #Partial ROC analyses for each model
    proc_res[[i]] <- proc[[1]]

    #Omission rates calculation
    om_rates[i] <- kuenm_omrat(model = mod1, threshold = threshold,
                                occ.tra = occ, occ.test = occ1)

    #Erasing calibration models after evaluating them if kept = FALSE
    if(kept == FALSE) {
      unlink(dir_names[i], recursive = T)
      unlink(dir_names1[i], recursive = T)
    }
  }
  if(.Platform$OS.type != "unix") {
    suppressMessages(close(pb))
  }
  n.mod <- i

  ##Erasing main folder of candidate models if kept = FALSE
  if(kept == FALSE) {
    unlink(path, recursive = T)
    cat("\nAll candidate models were deleted\n")
  }else{
    cat("\nAll candidate models were kept\n")
  }

  ##Creating the final tables
  ###From AICc analyses
  aiccs <- do.call(rbind, aiccs) #joining tables
  for (i in 1:length(aiccs[, 1])) {
    aiccs[i, 2] <- (aiccs[i, 1] - min(aiccs[, 1], na.rm = TRUE))
    aiccs[i, 3] <- (exp(-0.5 * aiccs[i,2])) / (sum(exp(-0.5 * aiccs[,2]), na.rm = TRUE))
  }

  ###From pROC analyses
  proc_res1 <- do.call(rbind, proc_res) #joining tables of the pROC results
  proc_res_m <- data.frame(mod_nam, proc_res1) #adding a new column with the number of AUC ratios interations < 1

  #####
  #Joining the results
  ku_enm_eval <- data.frame(proc_res_m, om_rates, aiccs)
  colnames(ku_enm_eval) <- c("Model", "Mean_AUC_ratio", "Partial_ROC",#changing column names in the final table
                             paste("Omission_rate_at_", threshold, "%", sep = ""), "AICc",
                             "delta_AICc", "W_AICc", "num_parameters")


  #####
  #Choosing the best models
  cat("\nSelecting the best candidate models...\n")

  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR") {
    if(selection == "OR_AICc") {
      ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
      ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] <= threshold / 100, ]
      if(length(ku_enm_best[, 4]) != 0) {
        for (i in 1:length(ku_enm_best[,1])) {
          ku_enm_best[i, 6] <- (ku_enm_best[i, 5] - min(ku_enm_best[, 5], na.rm = TRUE))
          ku_enm_best[i, 7] <- (exp(-0.5 * ku_enm_best[i, 6])) / (sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE))
        }
        ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]

      }else {
        cat(paste("\nNone of the significant candidate models met the omission rate criterion,",
                  "\nmodels with the smallest omission rate and lowest AICc will be presented\n"))

        ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] == min(ku_enm_bes[, 4]), ]
        for (i in 1:length(ku_enm_best[, 1])) {
          ku_enm_best[i, 6] <- (ku_enm_best[i, 5] - min(ku_enm_best[, 5], na.rm = TRUE))
          ku_enm_best[i, 7] <- (exp(-0.5 * ku_enm_best[i, 6])) / (sum(exp(-0.5 * ku_enm_best[i, 6]), na.rm = TRUE))
        }
        ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
      }
    }

    if(selection == "AICc") {
      ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
      ku_enm_best <- ku_enm_bes[ku_enm_bes[, 6] <= 2, ]
      if(length(ku_enm_best[, 6] != 0)) {
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
      }else {
        cat(paste("\nNone of the significant candidate models met the AICc criterion,",
                  "\ndelta AICc will be recalculated for significant models\n"))

        for (i in 1:length(ku_enm_best[, 6])) {
          ku_enm_best[i, 6] <- (ku_enm_best[i, 5] - min(ku_enm_best[, 5], na.rm = TRUE))
          ku_enm_best[i, 7] <- (exp(-0.5 * ku_enm_best[i, 6])) / (sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE))
        }
        ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
      }
    }

    if(selection == "OR") {
      ku_enm_b <- ku_enm_eval[!is.na(ku_enm_eval[, 3]), ]
      ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
      ku_enm_bes1 <- ku_enm_b[ku_enm_b[, 3] <= 0.05, ]
      ku_enm_best <- ku_enm_bes1[ku_enm_bes1[, 4] <= threshold / 100, ]
      if(length(ku_enm_best[, 4]) != 0) {
        if(length(ku_enm_best[, 4]) > 10) {
          ku_enm_best <- ku_enm_best[order(ku_enm_best[, 4]), ][1:10, ]
        }else {
          ku_enm_best <- ku_enm_best[order(ku_enm_best[, 4]), ]
        }
      }else {
        cat(paste("\nNone of the significant candidate models met the omission rate criterion,",
                  "\nmodels with the smallest omission rate will be presented\n"))

        ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] == min(ku_enm_bes[, 4]), ][1:10, ]
      }
    }
  }else {
  cat("\nNo valid model selection criterion has been defined,\n
        no file containing best models was created.\n
        Select your best models from the complete list.\n")
  }

  #####
  #Statistics of the process
  ##Counting
  ku_enm_sign <- ku_enm_eval[!is.na(ku_enm_eval[, 3]), ]
  ku_enm_sign <- ku_enm_sign[ku_enm_sign[, 3] <= 0.05, ]

  ku_enm_or <- ku_enm_eval[ku_enm_eval[, 4] <= 0.05, ]

  ku_enm_AICc <- ku_enm_eval[!is.na(ku_enm_eval[, 6]), ]
  ku_enm_AICc <- ku_enm_AICc[ku_enm_AICc[, 6] <= 2, ]

  ku_enm_best_OR <- ku_enm_sign[ku_enm_sign[, 4] <= threshold / 100, ]

  ku_enm_best_AICc <- ku_enm_bes[ku_enm_bes[, 6] <= 2, ]

  ku_enm_best_OR_AICc <- ku_enm_bes[ku_enm_bes[, 4] <= threshold / 100, ]
  if(length(ku_enm_best_OR_AICc[, 4]) != 0) {
    for (i in 1:length(ku_enm_best_OR_AICc[, 1])) {
      ku_enm_best_OR_AICc[i, 6] <- (ku_enm_best_OR_AICc[i, 5] - min(ku_enm_best_OR_AICc[, 5],
                                                                    na.rm = TRUE))
      ku_enm_best_OR_AICc[i, 7] <- (exp(-0.5 * ku_enm_best_OR_AICc[i, 6])) / (sum(exp(-0.5 * ku_enm_best_OR_AICc[, 6]),
                                                                                  na.rm = TRUE))
    }
    ku_enm_best_OR_AICc <- ku_enm_best_OR_AICc[ku_enm_best_OR_AICc[, 6] <= 2, ]
  }

  ##Preparing the table
  r_names <- c("All candidate models", "Statistically significant models", "Models meeting omission rate criteria",
               "Models meeting AICc criteria", "Statistically significant models meeting omission rate criteria",
               "Statistically significant models meeting AICc critera",
               "Statistically significant models meeting omission rate and AICc criteria")
  statis <- c(length(ku_enm_eval[, 1]),
              length(ku_enm_sign[, 3]),
              length(ku_enm_or[, 4]),
              length(ku_enm_AICc[, 6]),
              length(ku_enm_best_OR[, 4]),
              length(ku_enm_best_AICc[, 6]),
              length(ku_enm_best_OR_AICc[, 2]))

  ku_enm_stats <- data.frame(r_names, statis)
  colnames(ku_enm_stats) <- c("Criteria", "Number of models")

  #####
  #Writing the results
  ##csv files
  cat("\nWriting kuenm_ceval results...\n")
  dir.create(out.eval)

  name <- paste(out.eval, "calibration_results.csv", sep = "/")
  name0 <- paste(out.eval, "calibration_stats.csv", sep = "/")
  name1 <- paste(out.eval, "best_candidate_models_OR_AICc.csv", sep = "/")
  name2 <- paste(out.eval, "best_candidate_models_AICc.csv", sep = "/")
  name3 <- paste(out.eval, "best_candidate_models_OR.csv", sep = "/")


  write.csv(ku_enm_eval, file = name, eol = "\n", na = "NA", row.names = FALSE)
  write.csv(ku_enm_stats, file = name0, eol = "\n", na = "NA", row.names = FALSE)

  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR"){
    if(selection == "OR_AICc"){
      write.csv(ku_enm_best, file = name1, eol = "\n", na = "NA", row.names = FALSE)
    }
    if(selection == "AICc"){
      write.csv(ku_enm_best, file = name2, eol = "\n", na = "NA", row.names = FALSE)
    }
    if(selection == "OR"){
      write.csv(ku_enm_best, file = name3, eol = "\n", na = "NA", row.names = FALSE)
    }
  }

  ##Plot
  png(paste(out.eval, "calibration_figure.png", sep = "/"), width = 80, height = 80,
      units = "mm", res = 600)

  par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.58)
  plot(na.omit(ku_enm_eval)[, 4]~log(na.omit(ku_enm_eval)[, 5]),
       xlab = "Natural logarithm of AICc", ylab = paste("Omission rates at",
                                                        paste(threshold, "%", sep = ""),
                                                        "threshold value", sep = " "),
       las = 1, col = "grey35")

  points(na.omit(ku_enm_eval[!ku_enm_eval[, 1] %in% ku_enm_bes[, 1], ])[, 4]~log(na.omit(ku_enm_eval[!ku_enm_eval[, 1] %in% ku_enm_bes[, 1], ])[, 5]),
         col = "red1", pch = 19, cex = 1.1)

  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR") {
    if(selection == "OR_AICc") {
      points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
             col = "dodgerblue1", pch = 17, cex = 1.4)
      legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
             pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("dodgerblue1", "red1", "gray35"), bty = "n",
             inset = c(0.01, 0))
    }
    if(selection == "AICc") {
      points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
             col = "darkorchid1", pch = 17, cex = 1.4)
      legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
             pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("darkorchid1", "red1", "gray35"), bty = "n",
             inset = c(0.01, 0))
    }
    if(selection == "OR") {
      points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
             col = "orange2", pch = 17, cex = 1.4)
      legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
             pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("orange2", "red1", "gray35"), bty = "n",
             inset = c(0.01, 0))
    }else {
      cat("All selected models had NAs as AICc values, imposible to plot them.")
      legend("bottomright", legend = c("Non significant models", "All candidate models"),
             pt.cex = c(1.1, 1), pch = c(19, 1), col = c("red1", "gray35"), bty = "n",
             inset = c(0.01, 0))
    }
  }
  dev.off()

  ##html file
  ###Writing the html file
  html_eval(path = out.eval, file.name = "calibration_results")

  ##Retuning objects
  ###Dataframes in a list
  list_res <- list(ku_enm_stats, ku_enm_best, ku_enm_eval)
  names(list_res) <- c("Summary", "Best models",
                       "All models")

  ###Plot
  ku_enm_plot <- {
    par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.85)
    plot(na.omit(ku_enm_eval)[, 4]~log(na.omit(ku_enm_eval)[, 5]),
         xlab = "Natural logarithm of AICc", ylab = paste("Omission rates at",
                                                          paste(threshold, "%", sep = ""),
                                                          "threshold value", sep = " "),
         las = 1, col = "grey35")

    points(na.omit(ku_enm_eval[!ku_enm_eval[, 1] %in% ku_enm_bes[, 1], ])[, 4]~log(na.omit(ku_enm_eval[!ku_enm_eval[, 1] %in% ku_enm_bes[, 1], ])[, 5]),
           col = "red1", pch = 19, cex = 1.1)

    if(selection == "OR_AICc" | selection == "AICc" | selection == "OR") {
      if(selection == "OR_AICc") {
        points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
               col = "dodgerblue1", pch = 17, cex = 1.4)
        legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
               pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("dodgerblue1", "red1", "gray35"), bty = "n",
               inset = c(0.01, 0))
      }
      if(selection == "AICc") {
        points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
               col = "darkorchid1", pch = 17, cex = 1.4)
        legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
               pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("darkorchid1", "red1", "gray35"), bty = "n",
               inset = c(0.01, 0))
      }
      if(selection == "OR") {
        points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
               col = "orange2", pch = 17, cex = 1.4)
        legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
               pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("orange2", "red1", "gray35"), bty = "n",
               inset = c(0.01, 0))
      }else {
        cat("All selected models had NAs as AICc values, imposible to plot them.")
        legend("bottomright", legend = c("Non significant models", "All candidate models"),
               pt.cex = c(1.1, 1), pch = c(19, 1), col = c("red1", "gray35"), bty = "n",
               inset = c(0.01, 0))
      }
    }
  }

  #####
  #Finalizing the function
  cat("\nProcess finished\n")
  cat(paste("A folder containing results of the calibration of", n.mod,
            "\ncandidate models has been written\n", sep = " "))

  cat(paste("\nThe folder", out.eval, "contains:\n", sep = " "))
  cat("   -A html file and its dependencies that sum all the results, check\n")
  cat(paste("    ", "calibration_results.html\n", sep = ""))

  cat("   -Two csv files with all models' calibration results and stats,\n")

  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR") {
    if(selection == "OR_AICc"){
      cat("    and an aditional csv file containing the best models selected by OR and AICc.\n")
    }
    if(selection == "AICc") {
      cat("    and an aditional csv file containing the best models selected by AICc.\n")
    }
    if(selection == "OR") {
      cat("    and an aditional csv file containing the best models selected by OR.\n")
    }
  }

  cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
  return(list_res)
}

n.par <- function(x) {
  lambdas <- x[1:(length(x) - 4)]
  countNonZeroParams <- function(x) {
    if (strsplit(x, split = ", ")[[1]][2] != "0.0")
      1
  }
  no.params <- sum(unlist(sapply(lambdas, countNonZeroParams)))
  return(no.params)
}
