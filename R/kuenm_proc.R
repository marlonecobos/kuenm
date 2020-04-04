#' Partial ROC calculation for ecological niche models
#'
#' @description kuenm_proc applies partial ROC tests to model predictions.
#'
#' @param model RasterLayer or numeric vector of ecological niche model
#' predictions to be evaluated. If RasterLayer, layer of predicted suitability.
#' If numeric vector, predicted suitability values.
#' @param occ.test matrix, data.frame, or numeric vector containing coordinates
#' of occurrences to test model predictions to be evaluated. If matrix or
#' data.frame, columns must include longitude and latitude in that order.
#' If numeric, values of suitability in such occurrences. If a matrix or a
#' data.frame is provided, \code{model} must be a RasterLayer.
#' @param threshold (numeric) value from 0 to 100 to represent the percentage of
#' potential error (E) that the data could have due to any source of uncertainty.
#' Default = 5.
#' @param iterations (numeric) number of bootstrap iterations to be performed;
#' default = 500.
#' @param rand.percent (numeric) percentage of testing data to be used in each
#' bootstrapped process for calculating the partial ROC. Default = 50.
#' @param parallel (logical) if TRUE, calculations will be performed in parallel
#' using the available cores of the computer. This will demand more RAM and
#' almost full use of the CPU; hence, its use is recommended in high-performance
#' computers. Using this option will speed up the analyses only if \code{model}
#' is a large RasterLayer or if \code{iterations} are more than 5000.
#' Default = FALSE.
#'
#' @return A list with the summary of the results and a data.frame containing
#' the AUC values and AUC ratios calculated for all iterations.
#'
#' @usage
#' kuenm_proc(occ.test, model, threshold = 5, rand.percent = 50,
#'            iterations = 500, parallel = FALSE)
#'
#' @details Partial ROC is calculated following Peterson et al. (2008;
#' \url{http://dx.doi.org/10.1016/j.ecolmodel.2007.11.008}).
#'
#' @importFrom purrr map_df
#' @useDynLib kuenm
#' @export
#'
#' @examples
#' occ <- read.csv(list.files(system.file("extdata", package = "kuenm"),
#'                            pattern = "sp_test.csv", full.names = TRUE))
#' model <- raster::raster(list.files(system.file("extdata", package = "kuenm"),
#'                                    pattern = "sp_model.tif", full.names = TRUE))
#' thres <- 5
#' rand_perc <- 50
#' iterac <- 500
#'
#' p_roc <- kuenm_proc(occ.test = occ, model = model, threshold = thres,
#'                    rand.percent = rand_perc, iterations = iterac)

kuenm_proc <- function(occ.test, model, threshold = 5, rand.percent = 50,
                       iterations = 500, parallel = FALSE) {

  # -----------
  # detecting potential errors, other potential problems tested in code
  if (missing(model)) {
    stop("Argument 'model' is necessary to perform the analysis.")
  }
  if (missing(occ.test)) {
    stop("Argument 'occ.test' is necessary to perform the analysis.")
  }
  c_pred <- class(model)[1]
  if (!c_pred %in% c("RasterLayer", "numeric")) {
    stop("'model' must be of class RasterLayer or numeric.")
  }
  c_tdat <- class(occ.test)[1]
  if (!c_tdat %in% c("matrix", "data.frame", "numeric")) {
    stop("'occ.test' must be of class matrix, data.frame, or numeric.")
  }
  if (c_pred == "numeric" & c_tdat != "numeric") {
    stop("'occ.test' must be of class numeric if model is a numeric vector.")
  }


  # -----------
  # package and function needed
  suppressPackageStartupMessages(library(dplyr))

  calc_aucDF <- function(big_classpixels, fractional_area, test_data, n_data,
                         n_samp, error_sens) {
    rowsID <- sample(x = n_data, size = n_samp, replace = TRUE)
    test_data1 <- test_data[rowsID]
    omssion_matrix <- big_classpixels > test_data1
    sensibility <- 1 - colSums(omssion_matrix) / n_samp
    xyTable <- data.frame(fractional_area, sensibility)
    less_ID <- which(xyTable$sensibility <= error_sens)
    xyTable <- xyTable[-less_ID, ]
    xyTable <- xyTable[order(xyTable$fractional_area, decreasing = F), ]
    auc_pmodel <- kuenm:::trap_roc(xyTable$fractional_area, xyTable$sensibility)
    auc_prand <- kuenm:::trap_roc(xyTable$fractional_area, xyTable$fractional_area)
    auc_ratio <- auc_pmodel / auc_prand
    auc_table <- data.frame(auc_pmodel, auc_prand, auc_ratio = auc_ratio )
    return(auc_table)
  }

  # -----------
  # preparing data
  if (c_pred == "RasterLayer") {model <- raster::setMinMax(model)}
  min_pred <- ifelse(c_pred == "numeric", min(model, na.rm = TRUE),
                     model@data@min)
  max_pred <- ifelse(c_pred == "numeric", max(model, na.rm = TRUE),
                     model@data@max)

  model <- round((model / max_pred) * 1000)
  if (c_pred == "RasterLayer") {
    if (c_tdat != "numeric") {
      test_data <- na.omit(raster::extract(model,
                                           occ.test[, 1:2]))
    } else {
      test_data <- round((occ.test / max_pred) * 1000)
    }
    classpixels <- data.frame(raster::freq(model, useNA = "no"))
  } else {
    test_data <- round((occ.test / max_pred) * 1000)
    vals <- na.omit(unique(model))
    classpixels <- data.frame(value = vals, count = c(table(model)),
                              row.names = 1:length(vals))
  }

  # -----------
  # analysis
  if(min_pred == max_pred){
    warning("\nmodel has no variability, pROC will return NA.\n")

    p_roc <- rep(NA, 2)
    names(p_roc) <- c(paste0("Mean_AUC_ratio_at_", threshold, "%"), "pval_pROC")

    auc_ratios <- rep(NA, 3)
    names(auc_ratios) <- c("Model_partial_AUC", "Random_curve_partial_AUC",
                           "AUC_ratio")

    p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)

    return(p_roc_res)
  }else {
    classpixels <- classpixels %>%
      dplyr::mutate_(value = ~rev(value),
                     count = ~rev(count),
                     totpixperclass = ~cumsum(count),
                     percentpixels = ~totpixperclass/sum(count)) %>%
      dplyr::arrange(value)


    error_sens <- 1 - (threshold / 100)
    prediction_errors <- classpixels[, "value"]
    fractional_area <- classpixels[, "percentpixels"]
    n_data <- length(test_data)
    n_samp <- ceiling((rand.percent / 100) * n_data)

    big_classpixels <- matrix(rep(prediction_errors, each = n_samp),
                              ncol = length(prediction_errors))

    if(parallel){
      future::plan(future::multiprocess)
      roc_env <- new.env()
      n_cores <- future::availableCores()
      niter_big <- floor(iterations/n_cores)
      n_runs <- rep(niter_big,n_cores)
      sum_n_runs <- sum(n_runs)
      n_runs[1] <- n_runs[1] + (iterations - sum_n_runs)

      for(i in 1:length(n_runs)){
        x <- as.character(i)
        roc_env[[x]] %<-% {
          x1 <- 1:n_runs[i]
          auc_matrix1 <- x1 %>%
            purrr::map_df(~calc_aucDF(big_classpixels,
                                      fractional_area,
                                      test_data,n_data,n_samp,
                                      error_sens))
        }
      }
      partial_AUC <- as.list(roc_env)
      rm(roc_env)
      partial_AUC <- do.call(rbind.data.frame,partial_AUC)
      rownames(partial_AUC) <- NULL
      future::plan(future::sequential)

    }else{
      partial_AUC <- 1:iterations %>%
        purrr::map_df(~calc_aucDF(big_classpixels, fractional_area, test_data,
                                  n_data, n_samp, error_sens))
    }

    naID <- !is.na(partial_AUC$auc_ratio)
    nona_valproc <- partial_AUC$auc_ratio[naID]
    mauc <- mean(nona_valproc)
    proc <- sum(nona_valproc <= 1) / length(nona_valproc)

    p_roc <- c(mauc, proc)
    names(p_roc) <- c(paste0("Mean_AUC_ratio_at_", threshold, "%"), "pval_pROC")

    auc_ratios <- partial_AUC
    names(auc_ratios) <- c("Model_partial_AUC", "Random_curve_partial_AUC",
                           "AUC_ratio")

    p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)

    return(p_roc_res)
  }
}
