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
#' @param parallel (logical) argument deprecated. Default = NULL.
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
#' data("sp_test", package = "kuenm")
#' model <- raster::raster(system.file("extdata/sp_model.tif",
#'                                        package = "kuenm"))
#' thres <- 5
#' rand_perc <- 50
#' iterac <- 500
#'
#' p_roc <- kuenm_proc(occ.test = sp_test, model = model, threshold = thres,
#'                    rand.percent = rand_perc, iterations = iterac)

kuenm_proc <- function(occ.test, model, threshold = 5, rand.percent = 50,
                       iterations = 500, parallel = NULL) {

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
  if (!is.null(parallel)) {
    message("'parallel' is a deprecated argument.")
  }


  # -----------
  # package needed
  suppressPackageStartupMessages(library(dplyr))

  # -----------
  # preparing data
  if (c_pred == "RasterLayer") {
    model <- raster::setMinMax(model)
  }

  min_pred <- ifelse(c_pred == "numeric", min(model, na.rm = TRUE),
                     model@data@min)
  max_pred <- ifelse(c_pred == "numeric", max(model, na.rm = TRUE),
                     model@data@max)

  if (c_pred == "RasterLayer") {
    if (c_tdat != "numeric") {
      test_data <- na.omit(raster::extract(model,
                                           occ.test[, 1:2]))
    } else {
      test_data <- na.omit(occ.test)
    }
  } else {
    test_data <- na.omit(occ.test)
  }

  vals <- na.omit(model[])

  # ndec <- dec_places_proc(min_pred, min(test_data))
  # fix_dec <- as.numeric(paste0("1e+", ndec))
  #
  # test_data <- test_data * fix_dec
  # vals <- vals * fix_dec
  #
  # minmin <- min(c(vals, test_data))
  #
  # test_data <- round(test_data / minmin)
  # vals <- round(vals / minmin)

  nvals <- length(vals)
  vals <- c(vals, test_data)
  vals <- as.numeric(cut(vals, 500))
  test_data <- vals[(nvals + 1):length(vals)]
  vals <- vals[1:nvals]

  classpixels <- as.data.frame(table(vals), stringsAsFactors = FALSE)
  colnames(classpixels) <- c("value", "count")

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

  } else {
    classpixels <- classpixels %>%
      dplyr::mutate(value = rev(value),
                    count = rev(count),
                    totpixperclass = cumsum(count),
                    percentpixels = totpixperclass/sum(count)) %>%
      dplyr::arrange(value)

    error_sens <- 1 - (threshold / 100)
    prediction_errors <- classpixels[, "value"]
    fractional_area <- classpixels[, "percentpixels"]
    n_data <- length(test_data)
    n_samp <- ceiling((rand.percent / 100) * n_data)

    big_classpixels <- matrix(rep(prediction_errors, each = n_samp),
                              ncol = length(prediction_errors))

    st <- Sys.time()
    partial_AUC <- 1:iterations %>%
      purrr::map_df(~calc_aucDF(big_classpixels, fractional_area, test_data,
                                n_data, n_samp, error_sens))
    ed <- Sys.time()
    print(ed - st)

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
  }

  return(p_roc_res)
}


dec_places_proc <- function(model_vals, test_vals) {

  x <- c(model_vals, test_vals)
  x <- x[abs(x - round(x)) > .Machine$double.eps^0.5]

  x <- do.call(rbind, strsplit(sub('0+$', '', as.character(x)), ".",
                              fixed = TRUE))[, 2]

  return(max(nchar(x), na.rm = TRUE))
}


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
