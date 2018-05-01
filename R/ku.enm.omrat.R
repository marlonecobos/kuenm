#' Omission rates calculation for single models
#'
#' @description ku.enm.omrat calculates omission rates of geographic projections
#' of ecological niche models based on one or multiple user-specified thresholds.
#'
#' @param model a RasterLayer of the model to be evaluated.
#' @param threshold (numeric vector) value(s) from 0 to 100 that will be used as thresholds,
#' default = 5.
#' @param occ.tra a numerical matrix containing coordinates of the occurrence data used to create
#' the ecological niche model to be evaluated; columns must be: longitude and latitude.
#' @param occ.test a numerical matrix containing coordinates of the occurrences used to test
#' the ecological niche model to be evaluated; columns must be: longitude and latitude.
#'
#' @return A named numeric value or numeric vector with the result(s).
#'
#' @examples
#' # single threshold
#' model <- sp_mod
#' thres <- 5
#' octr <- sp_train
#' octe <- sp_test
#'
#' om_rate <- ku.enm.omrat(model, threshold = thres,
#'                         occ.tra = octr, occ.test = octe)
#'
#' # multiple thresholds
#' thres1 <- c(5, 10, 20)
#'
#' om_rate <- ku.enm.omrat(model, threshold = thres1,
#'                         occ.tra = octr, occ.test = octe)

ku.enm.omrat <- function(model, threshold = 5, occ.tra, occ.test) {

  if(min(na.omit(raster::getValues(model))) == max(na.omit(raster::getValues(model)))) {
    warning("\nModel imput has no variability, omission rate will return NA.\n")

    om_rate <- NA
  }else {
    om_rate <- vector("numeric", length = length(threshold))

    for (i in 1:length(threshold)) {
      suit_val_cal <- na.omit(raster::extract(model, occ.tra))
      suit_val_eval <- na.omit(raster::extract(model, occ.test))
      val <- ceiling(length(occ.tra[, 1]) * threshold[i] / 100) + 1
      omi_val_suit <- sort(suit_val_cal)[val]
      om_rate[i] <- as.numeric(length(suit_val_eval[suit_val_eval < omi_val_suit]) / length(suit_val_eval))
    }

    names(om_rate) <- paste("om_rate_", threshold, "%", sep = "")
    return(om_rate)
  }
}
