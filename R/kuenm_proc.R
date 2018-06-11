#' Partial ROC calculation of single models
#'
#' @description kuenm_proc applies partial ROC tests to single models.
#'
#' @param occ.test a numerical matrix containing coordinates of the occurrences used to test
#' the ecological niche model to be evaluated; columns must be: longitude and latitude.
#' @param model a RasterLayer of the ecological niche model to be evaluated.
#' @param threshold (numeric) value from 0 to 100 that will be used as threshold (E);
#' default = 5.
#' @param rand.percent (numeric) value from 0 to 100 representing the percent of testing data
#' to be used for performing the bootstrap process for calculating the partial ROC;
#' default = 50.
#' @param iterations (numeric) number of bootstrap iterations to be performed;
#' default = 1000.
#'
#' @return A list containing a named vector with the final partial ROC results and
#' a matrix containing the AUC values and AUC ratios calculated for each iteration.
#'
#' @details Partial ROC is calculated following Peterson et al.
#' (2008; \url{https://doi.org/10.1016/j.ecolmodel.2007.11.008}). This function is a modification
#' of the \code{\link[ENMGadgets]{PartialROC}} funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
#'
#' @examples
#' occ <- sp_test
#' model <- sp_mod
#' thres <- 5
#' rand_perc <- 50
#' iterac <- 100
#'
#' p_roc <- kuenm_proc(occ.test = occ, model = model, threshold = thres,
#'                    rand.percent = rand_perc, iterations = iterac)

kuenm_proc <- function(occ.test, model, threshold = 5, rand.percent = 50,
                        iterations = 1000) {

  suppressMessages(library(dplyr))

  if(raster::cellStats(model,"min") == raster::cellStats(model,"max")) {
    warning("\nModel with no variability, pROC will return NA.\n")

    p_roc <- rep(NA, 2)
    names(p_roc) <- c(paste("Mean_AUC_ratio_at_", threshold, "%", sep = ""), "Partial_ROC")

    auc_ratios <- rep(NA, 4)
    names(auc_ratios) <- c("Iteration", paste("AUC_at_", 100 - threshold, "%", sep = ""),
                           paste("AUC_at_", threshold, "%", sep = ""), "AUC_ratio")

    p_roc_res <- list(p_roc, auc_ratios)

    return(p_roc_res)
  }else {
    omissionval <- (100 - threshold) / 100

    inrastlog <- model

    ## Currently fixing the number of classes to 100. But later flexibility should be given in the parameter.
    inrast <- round((inrastlog / raster::cellStats(inrastlog, max)) * 1000)

    ## This function should be called only once outside the loop. This function generates values for x-axis.
    ## As x-axis is not going to change.
    classpixels <- a_pred_pres(inrast)

    occur <- occ.test
    extrast <- raster::extract(inrast, occur)

    ## Remove all the occurrences in the class NA. As these points are not used in the calibration.
    occurtbl <- cbind(occur, extrast)
    occurtbl <- occurtbl[which(is.na(occurtbl[, 3]) == FALSE), ]

    pointid <- seq(1:nrow(occurtbl))
    occurtbl <- cbind(pointid, occurtbl)
    names(occurtbl) <- c("PointID", "Longitude", "Latitude", "ClassID")

    ## Partial ROC iterations
    output_auc <- parallel::mclapply(1:(iterations),
                                     function(x) auc_comp(x, occurtbl,
                                                          rand.percent,
                                                          omissionval,
                                                          classpixels))
    auc_ratios <- data.frame(t(sapply(output_auc, c)))


    mauc <- apply(auc_ratios, 2, mean)[4] #mean of AUC ratios interations
    proc <- sum(auc_ratios[, 4] <= 1) / length(auc_ratios[, 4]) #proportion of AUC ratios <= 1
    p_roc <- c(mauc, proc)
    names(p_roc) <- c(paste("Mean_AUC_ratio_at_", threshold, "%", sep = ""), "Partial_ROC")

    p_roc_res <- list(p_roc, auc_ratios)

    return(p_roc_res)
  }
}
