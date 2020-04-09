#' AICc calculation of Maxent SWD predictions
#'
#' @description aicc calculates the Akaike information criterion corrected for
#' small sample sizes (AICc) to single or multiple models produced with Maxent.
#'
#' @param occ numerical matrix containing coordinates of the occurrences used to
#' create the ecological niche models to be evaluated; columns must be: longitude
#' and latitude.
#' @param prediction matrix of longitude and latidue coordinates, and Maxent Raw
#' predictions obtained using the SWD format. Prediction coordinates must include
#' the ones in \code{occ}
#' @param npar (numeric) number of parameters of the model. Use function \code{\link{n.par}}
#' to obtain number of parameters for each model from lambdas file.
#'
#' @return A data.frame with results from AICc analyses.
#'
#' @export
#'
#' @details
#' This function is a modification of calc.aicc from the ENMeval package to allow
#' calculations using numerical values instead raster predictions.

aicc <- function(occ, prediction, npar) {
  AIC.valid <- npar < nrow(occ)
  if (nrow(prediction) == 0) {
    res <- data.frame(cbind(AICc = NA, delta_AICc = NA,
                            weight_AICc = NA, parameters = npar))
    warning("Cannot calculate AICc when prediction has 0 rows.")
  } else {
    vals <- prediction[paste(prediction[, 1], prediction[, 2]) %in%
                         paste(occ[, 1], occ[, 2]), 3]
    vals <- na.omit(vals)
    probsum <- sum(prediction[, 3], na.rm = TRUE)
    LL <- colSums(log(t(t(vals)/probsum)), na.rm = TRUE)
    AICc <- ((2 * npar) - (2 * LL)) + (2 * npar * (npar + 1) /
                                                 (nrow(occ) - npar - 1))
    AICc[AIC.valid == FALSE] <- NA
    AICc[is.infinite(AICc)] <- NA
    if (sum(is.na(AICc)) == length(AICc)) {
      warning("AICc not valid: too many parameters, or likelihood = Inf... returning NA.")
      res <- data.frame(cbind(AICc, delta_AICc = NA, weight_AICc = NA,
                              parameters = npar))
    } else {
      delta_AICc <- AICc - min(AICc, na.rm = TRUE)
      weight_AICc <- exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc), na.rm = TRUE)
      res <- data.frame(AICc, delta_AICc, weight_AICc, parameters = npar)
      rownames(res) <- NULL
    }
  }
  rownames(res) <- NULL
  return(res)
}


#' Omission rates calculation for Maxent SWD predictions
#'
#' @description or calculates omission rates of numerical projections of ecological
#' niche models based on one or multiple user-specified thresholds.
#'
#' @param prediction matrix of longitude and latidue coordinates, and Maxent
#' prediction obtained using the SWD format. Prediction coordinates must include
#' the ones in \code{occ.tra}, and \code{occ.test}.
#' @param occ.tra numerical matrix containing coordinates of the occurrence data
#' used to create the prediction to be evaluated; columns must be: longitude and
#' latitude.
#' @param occ.test numerical matrix containing coordinates of the occurrences
#' used to test the prediction to be evaluated; columns must be: longitude and
#' latitude.
#' @param threshold (numeric) vector of value(s) from 0 to 100 that will be used
#' as thresholds, default = 5.
#'
#' @return A named numeric value or numeric vector with the result(s).
#'
#' @export

or <- function(prediction, occ.tra, occ.test, threshold = 5) {
  if (min(prediction, na.rm = T) == max(prediction, na.rm = T)) {
    warning("Model imput has no variability, omission rate = NA.")
    om_rate <- NA
  } else {
    vals <- prediction[paste(prediction[, 1], prediction[, 2]) %in%
                         paste(occ.tra[, 1], occ.tra[, 2]), 3]
    tvals <- prediction[paste(prediction[, 1], prediction[, 2]) %in%
                          paste(occ.test[, 1], occ.test[, 2]), 3]

    vals <- na.omit(vals); tvals <- na.omit(tvals)

    om_rate <- vector("numeric", length = length(threshold))
    for (i in 1:length(threshold)) {
      val <- ceiling(nrow(occ.tra) * threshold[i] / 100) + 1
      omi_val_suit <- sort(vals)[val]
      om_rate[i] <- length(tvals[tvals < omi_val_suit]) / length(tvals)
    }
    names(om_rate) <- paste("om_rate_", threshold, "%", sep = "")
  }
  return(om_rate)
}


#' Partial ROC, omission rates, and AICc calculations in concert (helper)
#'
#' @description proc_or_aicc performs a series of step by step processes that
#' help to read files from directores, extract necessary data, and evaluate
#' Maxent predictions based on partial ROC, omission rates, and AICc values.
#'
#' @param occ.joint (character) the name of csv file with training and testing
#' occurrences combined; columns must be: species, longitude, and latitude.
#' @param occ.tra (character) the name of the csv file with the training
#' occurrences; columns as in \code{occ.joint}.
#' @param occ.test (character) the name of the csv file with the evaluation
#' occurrences; columns as in \code{occ.joint}.
#' @param raw.folders (character) vector of names of directories containing
#' models created with all occurrences and raw outputs.
#' @param log.folders (character) vector of names of directories containing
#' models created with training occurrences and logistic outputs.
#' @param threshold (numeric) the percentage of training data omission error
#' allowed (E); default = 5.
#' @param rand.percent (numeric) the percentage of data to be used for the
#' bootstraping process when calculating partial ROCs; default = 50.
#' @param iterations (numeric) the number of times that the bootstrap is going
#' to be repeated; default = 500.
#' @param kept (logical) if FALSE, all candidate models will be erased after
#' evaluation, default = TRUE.
#'
#' @return
#' A data.frame with the results of partial ROC, omission rates, and AICc metrics
#' for all candidate models.
#'
#' @export
#'
#' @usage
#' proc_or_aicc(occ.joint, occ.tra, occ.test, raw.folders, log.folders,
#'              threshold = 5, rand.percent = 50, iterations = 500, kept = TRUE)
#'
#' @export

proc_or_aicc <- function(occ.joint, occ.tra, occ.test,
                         raw.folders, log.folders, threshold = 5,
                         rand.percent = 50, iterations = 500, kept = TRUE) {
  #pROCs, omission rates, and AICcs calculation
  message("Evaluation using partial ROC, omission rates, and AICc")

  # Slash
  if(.Platform$OS.type == "unix") {sl <- "/"; dl <- "/"} else {sl <- "\\"; dl <- "\\\\"}

  # model names
  model_names <- gsub(paste0(".*", dl), "", gsub("_all$", "", raw.folders))

  # occurrences
  oc <- read.csv(occ.joint)
  spn <- gsub(" ", "_", as.character(oc[1, 1]))
  oc <- oc[, -1]
  occ <- read.csv(occ.tra)[, -1]
  occ1 <- read.csv(occ.test)[, -1]

  longitude <- colnames(oc)[1]
  latitude <- colnames(oc)[2]

  aics <- list()
  proc_res <- list()
  om_rates <- numeric()
  nm <- length(raw.folders)

  if(.Platform$OS.type == "unix") {
    pb <- txtProgressBar(min = 0, max = nm, style = 3)
  } else {
    pb <- winProgressBar(title = "Progress bar", min = 0, max = nm, width = 300)
  }

  for(i in 1:nm) {
    Sys.sleep(0.1)
    if(.Platform$OS.type == "unix") {
      setTxtProgressBar(pb, i)
    } else {
      setWinProgressBar(pb, i, title = paste(round(i / nm * 100, 2),
                                             "% of the evaluation has finished"))
    }

    #AICc calculation
    lbds <- paste0(raw.folders[i], sl, spn, ".lambdas")
    waiting <- wait_written_done(lbds)
    lambdas <- readLines(lbds)

    par_num <- n.par(lambdas)

    mods <- paste0(raw.folders[i], sl, spn, ".csv")
    waiting <- wait_written_done(mods)
    mod <- read.csv(mods)

    aic <- aicc(oc, mod, par_num)
    aics[[i]] <- aic

    #pROCs and omission rates calculation
    mods1 <- paste0(log.folders[i], sl, spn, ".csv")
    waiting <- wait_written_done(mods1)
    mod1 <- read.csv(mods1)

    tval <- mod1[paste(mod1[, 1], mod1[, 2]) %in% paste(occ1[, 1], occ1[, 2]), 3]
    proc <- kuenm_proc(tval, mod1[, 3], threshold, rand.percent, iterations)

    proc_res[[i]] <- proc[[1]]

    om_rates[i] <- or(mod1, occ, occ1, threshold)

    #Erasing calibration models after evaluating them if kept = FALSE
    if(kept == FALSE) {
      unlink(raw.folders[i], recursive = T)
      unlink(log.folders[i], recursive = T)
    }
  }

  if(.Platform$OS.type != "unix") {suppressMessages(close(pb))}

  # From AICc analyses few calculations
  aiccs <- do.call(rbind, aics)

  aiccs[, 2] <- aiccs[, 1] - min(aiccs[, 1], na.rm = TRUE)
  aiccs[, 3] <- exp(-0.5 * aiccs[, 2]) / sum(exp(-0.5 * aiccs[, 2]), na.rm = TRUE)

  # From pROC analyses
  proc_res_m <- data.frame(model_names, do.call(rbind, proc_res))[, 1:3]

  # Joining the results
  ku_enm_eval <- data.frame(proc_res_m, om_rates, aiccs)
  colnames(ku_enm_eval) <- c("Model", "Mean_AUC_ratio", "pval_pROC",
                             paste0("Omission_rate_at_", threshold, "%"), "AICc",
                             "delta_AICc", "W_AICc", "N_parameters")

  return(ku_enm_eval)
}


#' Helper to summarize all results from model calibration exercises
#'
#' @param proc.or.aicc.results data.frame with results from evaluation of all
#' candidate models. Generally the output of \code{\link{proc_or_aicc}}.
#' @param selection (character) model selection criterion, can be "OR_AICc",
#' "AICc", or "OR"; OR = omission rates. Default = "OR_AICc", which means that
#' among models that are statistically significant and that present omission
#' rates below the threshold, those with delta AICc up to 2 will be selected.
#' See details for other selection criteria.
#'
#' @details
#' Other selecton criteria are described below: If "AICc" criterion is chosen,
#' all significant models with delta AICc up to 2 will be selected If "OR" is
#' chosen, the 10 first significant models with the lowest omission rates will
#' be selected.
#'
#' @return
#' A list with all results that need to be written to produce the evaluation report.
#'
#' @export

summary_calibration <- function(proc.or.aicc.results, selection = "OR_AICc") {

  ku_enm_eval <- proc.or.aicc.results
  threshold <- gsub("Omission_rate_at_", "", colnames(ku_enm_eval)[4])
  threshold <- as.numeric(gsub("%", "", threshold))

  # Choosing the best models
  if(selection == "OR_AICc") {
    ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
    ku_enm_best <- na.omit(ku_enm_bes[which(ku_enm_bes[, 4] <= threshold / 100), ])
    if(length(ku_enm_best[, 4]) != 0) {
      ku_enm_best[, 6] <- ku_enm_best[, 5] - min(ku_enm_best[, 5], na.rm = TRUE)
      ku_enm_best[, 7] <- exp(-0.5 * ku_enm_best[, 6]) /
        sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE)
      ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
      ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]

    }else {
      message("None of the significant candidate models met the omission rate criterion,",
              "\nmodels with the lowest omission rate and lowest AICc will be presented")

      ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] == min(ku_enm_bes[, 4]), ]
      ku_enm_best[, 6] <- ku_enm_best[, 5] - min(ku_enm_best[, 5], na.rm = TRUE)
      ku_enm_best[, 7] <- exp(-0.5 * ku_enm_best[, 6]) /
        sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE)
      ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
      ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
    }
  }

  if(selection == "AICc") {
    ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
    ku_enm_best <- ku_enm_bes[ku_enm_bes[, 6] <= 2, ]
    if(length(ku_enm_best[, 6]) != 0) {
      ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
    }else {
      message("None of the significant candidate models met the AICc criterion,",
              "\ndelta AICc will be recalculated for significant models")

      ku_enm_best[, 6] <- ku_enm_best[, 5] - min(ku_enm_best[, 5], na.rm = TRUE)
      ku_enm_best[, 7] <- exp(-0.5 * ku_enm_best[, 6]) /
        sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE)
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
      message("None of the significant candidate models met the omission rate criterion,",
              "\nmodels with the lowest omission rate will be presented")

      ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] == min(ku_enm_bes[, 4]), ][1:10, ]
    }
  }

  #####
  #Statistics of the process
  ##Counting
  ku_enm_sign <- ku_enm_eval[!is.na(ku_enm_eval[, 3]), ]
  ku_enm_sign <- ku_enm_sign[ku_enm_sign[, 3] <= 0.05, ]

  ku_enm_or <- ku_enm_eval[ku_enm_eval[, 4] <= threshold / 100, ]

  ku_enm_AICc <- ku_enm_eval[!is.na(ku_enm_eval[, 6]), ]
  ku_enm_AICc <- ku_enm_AICc[ku_enm_AICc[, 6] <= 2, ]

  ku_enm_best_OR <- ku_enm_sign[ku_enm_sign[, 4] <= threshold / 100, ]

  ku_enm_best_AICc <- ku_enm_bes[ku_enm_bes[, 6] <= 2, ]

  ku_enm_best_OR_AICc <- ku_enm_bes[ku_enm_bes[, 4] <= threshold / 100, ]
  if(length(ku_enm_best_OR_AICc[, 4]) != 0) {
    ku_enm_best_OR_AICc[, 6] <- ku_enm_best_OR_AICc[, 5] -
      min(ku_enm_best_OR_AICc[, 5], na.rm = TRUE)
    ku_enm_best_OR_AICc[, 7] <- exp(-0.5 * ku_enm_best_OR_AICc[, 6]) /
      sum(exp(-0.5 * ku_enm_best_OR_AICc[, 6]), na.rm = TRUE)
    ku_enm_best_OR_AICc <- ku_enm_best_OR_AICc[ku_enm_best_OR_AICc[, 6] <= 2, ]
  }

  # Preparing the table
  r_names <- c("All candidate models", "Statistically significant models",
               "Models meeting omission rate criteria",
               "Models meeting AICc criteria",
               "Statistically significant models meeting omission rate criteria",
               "Statistically significant models meeting AICc criteria",
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



  # returning results
  results <- list(calibration_statistics = ku_enm_stats,
                  selected_models = ku_enm_best,
                  calibration_results = ku_enm_eval,
                  threshold = threshold, significant_models = ku_enm_sign)
  return(results)
}
