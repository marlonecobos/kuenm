#' Partial ROC calculation of single models
#'
#' @description ku.enm.proc applies partial ROC tests to single models.
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
#' p_roc <- ku.enm.proc(occ.test = occ, model = model, threshold = thres,
#'                    rand.percent = rand_perc, iterations = iterac)

ku.enm.proc <- function(occ.test, model, threshold = 5, rand.percent = 50,
                        iterations = 1000) {

  if(min(na.omit(raster::getValues(model))) == max(na.omit(raster::getValues(model)))) {
    warning("\nModel with no variability, pROC will return NA.\n")

    p_roc <- rep(NA, 2)
    names(p_roc) <- c(paste("Mean_AUC_ratio_at_", threshold, "%", sep = ""), "Partial_ROC")

    auc_ratios <- rep(NA, 4)
    names(auc_ratios) <- c("Iteration", paste("AUC_value"),
                           paste("AUC_at_", threshold, "%", sep = ""), "AUC_ratio")

    p_roc_res <- list(p_roc, auc_ratios)

    return(p_roc_res)
  }else {
    inrastlog <- model

    ## Currently fixing the number of classes to 100. But later flexibility should be given in the parameter.
    inrast <- round((inrastlog / raster::cellStats(inrastlog, max)) * 1000)

    ## This function should be called only once outside the loop. This function generates values for x-axis.
    ## As x-axis is not going to change.
    classpixels <- a.pred.pres(inrast)

    occur <- occ.test
    extrast <- raster::extract(inrast, occur)

    ## Remove all the occurrences in the class NA. As these points are not used in the calibration.
    occurtbl <- cbind(occur, extrast)
    occurtbl <- occurtbl[which(is.na(occurtbl[, 3]) == FALSE), ]

    pointid <- seq(1:nrow(occurtbl))
    occurtbl <- cbind(pointid, occurtbl)
    names(occurtbl) <- c("PointID", "Longitude", "Latitude", "ClassID")

    ## Use option cl.cores to choose an appropriate cluster size.
    auc_ratio <- lapply(X = 1:iterations, FUN = function(x) {
      ll <- sample(nrow(occurtbl), round(rand.percent / 100 * nrow(occurtbl)), replace = TRUE)
      occurtbl1 <- occurtbl[ll, ]
      ## Generate the % points within each class in this table. Write SQL, using sqldf package
      occurinclass <- sqldf::sqldf("Select count(*), ClassID from occurtbl1 group by ClassID order by ClassID desc")
      occurinclass <- cbind(occurinclass, cumsum(occurinclass[, 1]),
                            cumsum(occurinclass[, 1]) / nrow(occurtbl1))
      names(occurinclass) <- c("OccuCount", "ClassID", "OccuSumBelow", "Percent")

      #### Raster file will contain all the classes in ClassID column, while
      #### occurrences table may not have all the classes.
      xytable <- g.xy.tab(classpixels, occurinclass)
      arearow <- calc.auc(xytable, threshold / 100, x)
      names(arearow) <- c("Iteration", paste("AUC_value"),
                          paste("AUC_at_", threshold, "%", sep = ""), "AUC_ratio")
      return(arearow)
    })

    auc_ratios <- as.data.frame(do.call(rbind, auc_ratio)) #converting each list of AUC ratios interations in a table
    mauc <- apply(auc_ratios, 2, mean)[4] #mean of AUC ratios interations
    proc <- sum(auc_ratios[, 4] <= 1) / length(auc_ratios[, 4]) #proportion of AUC ratios <= 1
    p_roc <- c(mauc, proc)
    names(p_roc) <- c(paste("Mean_AUC_ratio_at_", threshold, "%", sep = ""), "Partial_ROC")

    p_roc_res <- list(p_roc, auc_ratios)

    return(p_roc_res)
  }
}
