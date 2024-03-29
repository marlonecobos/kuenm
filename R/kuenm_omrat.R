#' Omission rates calculation for single models
#'
#' @description kuenm_omrat calculates omission rates of geographic projections
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
#' @export
#'
#' @examples
#' # single threshold
#' model <- raster::raster(system.file("extdata/sp_model.tif",
#'                                        package = "kuenm"))
#' thres <- 5
#' data("sp_train", package = "kuenm")
#' data("sp_test", package = "kuenm")
#'
#' om_rate <- kuenm_omrat(model, threshold = thres,
#'                         occ.tra = sp_train, occ.test = sp_test)
#'
#' # multiple thresholds
#' thres1 <- c(5, 10, 20)
#'
#' om_rate <- kuenm_omrat(model, threshold = thres1,
#'                         occ.tra = sp_train, occ.test = sp_test)

kuenm_omrat <- function(model, threshold = 5, occ.tra, occ.test) {

  if (missing(model)) {
    stop("Argument model is not defined.")
  }
  if (missing(occ.tra)) {
    stop("Argument occ.tra is not defined.")
  }
  if (missing(occ.test)) {
    stop("Argument occ.test is not defined.")
  }

  ran_mod <- range(na.omit(model[]))

  if(ran_mod[1] == ran_mod[2]) {
    warning("\nModel imput has no variability, omission rate will return NA.\n")

    om_rate <- rep(NA, length(threshold))
  }else {

    suit_val_cal <- na.omit(raster::extract(model, occ.tra))
    suit_val_eval <- na.omit(raster::extract(model, occ.test))

    om_rate <- vector("numeric", length = length(threshold))

    for (i in 1:length(threshold)) {
      val <- ceiling(length(occ.tra[, 1]) * threshold[i] / 100) + 1
      omi_val_suit <- sort(suit_val_cal)[val]
      om_rate[i] <- as.numeric(length(suit_val_eval[suit_val_eval < omi_val_suit]) / length(suit_val_eval))
    }
  }

  names(om_rate) <- paste("om_rate_", threshold, "%", sep = "")
  return(om_rate)
}
