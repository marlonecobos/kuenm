#' AICc calculation for Maxent models
#'
#' @description kuenm_aicc calculates the Akaike information criterion corrected
#' for small sample sizes (AICc) for single or multiple models produced with
#' Maxent.
#'
#' @param occ matrix or data.frame with coordinates of the occurrences used to
#' create the model (raster) to be evaluated; columns must be: longitude and
#' latitude.
#' @param model a RasterLayer or RasterStack of model projections created using
#' Maxent with raw outputs.
#' @param npar (numeric) number of parameters for \code{model}. Length must
#' correspond with number of layers in \code{model}. See function
#' \code{\link{n_par}}.
#'
#' @return
#' A data.frame containing values of AICc, delta AICc, weight of AICc, and
#' number of parameters. The number of rows of the data.frame corresponds to
#' the number of models evaluated.
#'
#' @details
#' Calculations are done following
#' [Warren and Seifert (2011)](https://doi.org/10.1890/10-1171.1).
#'
#' @usage
#' kuenm_aicc(occ, model, npar)
#'
#' @export
#'
#' @seealso
#' \code{\link{aicc}} for results obtained using the SWD format.
#'
#' @examples
#' data("sp_joint", package = "kuenm")
#' model <- raster::raster(system.file("extdata/sp_model_joint.tif",
#'                                              package = "kuenm"))
#'
#' lbds <- readLines(system.file("extdata/lambdas_model_joint.lambdas",
#'                               package = "kuenm"))
#' npar <- n_par(lbds) # counting number of parameters
#'
#' aicc <- kuenm_aicc(occ = sp_joint, model = model, npar = npar)

kuenm_aicc <- function (occ, model, npar) {
  if (missing(occ)) {
    stop("Argument 'occ' must be defined, see function's help.")
  }
  if (missing(model)) {
    stop("Argument 'model' must be defined, see function's help.")
  }
  if (!class(model)[1] %in% c("RasterLayer", "RasterStack", "RasterBrick")) {
    stop("'model' must be a RasterLayer or RasterStack object. See function's help.")
  }
  if (missing(npar)) {
    stop("Argument 'npar' must be defined, see function's help.")
  }
  if (dim(model)[3] != length(npar)) {
    stop("Number of 'models' to evaluate must correspond with length of 'npar'.")
  }

  AIC.valid <- npar < nrow(occ)
  if (dim(model)[3] == 0) {
    res <- data.frame(cbind(AICc = NA, delta.AICc = NA,
                            w.AIC = NA, parameters = npar))
    warning("Cannot calculate AICc when model = FALSE. Returning NA's.")
  } else {
    vals <- raster::extract(model, occ)
    probsum <- sum(raster::values(model), na.rm = TRUE)
    LL <- colSums(log(t(t(vals)/probsum)), na.rm = TRUE)
    AICc <- (2 * npar - 2 * LL) + (2 * (npar) * (npar + 1)/(nrow(occ) - npar - 1))
    AICc[AIC.valid == FALSE] <- NA
    AICc[is.infinite(AICc)] <- NA
    if (sum(is.na(AICc)) == length(AICc)) {
      warning("AICc not valid: too many parameters, or likelihood = Inf. Returning NA.")
      res <- data.frame(cbind(AICc, delta.AICc = NA, w.AIC = NA,
                              parameters = npar))
    } else {
      delta.AICc <- (AICc - min(AICc, na.rm = TRUE))
      w.AIC <- (exp(-0.5 * delta.AICc))/(sum(exp(-0.5 *
                                                   delta.AICc), na.rm = TRUE))
      res <- data.frame(AICc, delta.AICc, w.AIC, parameters = npar)
      rownames(res) <- NULL
    }
  }
  rownames(res) <- NULL
  return(res)
}

