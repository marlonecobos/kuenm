#' A set of occurrence records for ecological niche models
#'
#' A data.frame containing occurrence records of a tick (*Amblyomma americanum*)
#' across North America. The data combines records for training and testing.
#'
#' @name sp_joint
#'
#' @format A data frame with 178 rows and 2 columns.
#' \describe{
#'   \item{Longitude}{longitude, in decimal degrees.}
#'   \item{Latitude}{latitude, in decimal degrees.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' data("sp_joint", package = "kuenm")
#'
#' head(sp_joint)
NULL


#' A set of occurrence records to test candidate ecological niche models
#'
#' A data.frame containing occurrence records of a tick (*Amblyomma americanum*)
#' in North America, used to test candidate models during calibration.
#'
#' @name sp_test
#'
#' @format A data frame with 89 rows and 2 columns.
#' \describe{
#'   \item{Longitude}{longitude, in decimal degrees.}
#'   \item{Latitude}{latitude, in decimal degrees.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' data("sp_test", package = "kuenm")
#'
#' head(sp_test)
NULL



#' A set of occurrence records for training candidate ecological niche models
#'
#' A data.frame containing occurrence records of a tick (*Amblyomma americanum*)
#' across North America, used to train candidate models during calibration.
#'
#' @name sp_train
#'
#' @format A data frame with 89 rows and 2 columns.
#' \describe{
#'   \item{Longitude}{longitude, in decimal degrees.}
#'   \item{Latitude}{latitude, in decimal degrees.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' data("sp_train", package = "kuenm")
#'
#' head(sp_train)
NULL




#' A lambdas file resulted from a modeling process in Maxent
#'
#' A lambdas file resulted from a model created in Maxent with raw output for
#' *Amblyomma americanum* in North America. This file is used to calculate number
#' of parameters in the model, which is needed while calculating AICc values.
#'
#' @name sp_lambdas
#'
#' @format A lambdas file.
#' \describe{
#'   \item{parameters}{number of parameters in the Maxent model.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' lbds <- readLines(system.file("extdata/lambdas_model_joint.lambdas",
#'                               package = "kuenm"))
#'
#' head(lbds)
NULL




#' Raster variables masked to the area where a model is calibrated
#'
#' A RasterStack of predictor variables masked to the calibration area where
#' a model is calibrated. Variables represent four current bioclimatic variables
#' downloaded from the WorldClim database (\url{http://www.worldclim.org/}).
#'
#' @name mvars
#'
#' @format A RasterStack with 150 rows, 249 columns, 37350 cells, and 4 layers:
#' \describe{
#'   \item{Temperature}{temperature, in Celsius degrees times 10.}
#'   \item{Precipitation}{precipitation, in milimeters.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' mvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                   pattern = "Mbio_", full.names = TRUE))
#'
#' summary(mvars)
NULL





#' Variables masked to the area where a model is transferred
#'
#' A RasterStack containing predictor variables masked to the area where a model
#' is projected. Variables represent four future bioclimatic variables (2050) of
#' the NCAR-CCSM4 general circulation model under the RCP 8.5 emission scenario.
#'
#' @name gvars
#'
#' @format A RasterStack with 900 rows, 2160 columns, 1944000 cells, and 4 layers:
#' \describe{
#'   \item{Temperature}{temperature, in Celsius degrees times 10.}
#'   \item{Precipitation}{precipitation, in milimeters.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' gvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                   pattern = "Gbio_", full.names = TRUE))
#'
#' summary(gvars)
NULL





#' A raster output of an ecological niche model created with Maxent (logistic)
#'
#' A RasterLayer containing an ecological niche model for the tick
#' (*Amblyomma americanum*) that was created as part of the candidate models
#' during a calibration process.
#'
#' @name sp_model
#'
#' @format A RasterLayer with 150 rows, 249 columns, and 37350 cells:
#' \describe{
#'   \item{Suitability}{suitability values.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' sp_model <- raster::raster(system.file("extdata/sp_model.tif",
#'                                        package = "kuenm"))
#'
#' summary(sp_model)
NULL






#' A raster output of an ecological niche model created with Maxent (raw)
#'
#' A RasterLayer containing an ecological niche model for the a tick
#' (*Amblyomma americanum*) that was created with all occurrences.
#'
#' @name sp_mod_joint
#'
#' @format A RasterLayer with 150 rows, 249 columns, and 37350 cells:
#' \describe{
#'   \item{Suitability}{suitability values.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' sp_model_joint <- raster::raster(system.file("extdata/sp_model_joint.tif",
#'                                              package = "kuenm"))
#'
#' summary(sp_model_joint)
NULL
