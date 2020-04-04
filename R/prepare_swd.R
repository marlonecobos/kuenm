#' Prepare data for SWD maxent calibration processes
#'
#' @description prepare_swd helps to create csv files containing occurrence
#' records (all, train, and test records) and background coordinates, together
#' with values of predictor variables, that later can be used to run model
#' calibration in Maxent using the SWD format.
#'
#' @param occ data.frame containing occurrence records of the species of interest.
#' Mandatory columns are: species, longitude, and latitude. Other columns will
#' be ignored.
#' @param species (character) name of column containing species name.
#' @param longitude (character) name of column containing longitude values.
#' @param latitude (character) name of column containing latitude values.
#' @param data.split.method (character) name of the method to split training and
#' testing records. Default and only option for now = "random".
#' @param train.proportion (numeric) proportion of records to be used for training
#' models. Default = 0.5
#' @param raster.layers RasterStack of predictor variables masked to the area
#' where the model will be calibrated.
#' @param sample.size (numeric) number of points to represent the background for
#' the model. Default = 10000
#' @param save (logical) whether or not to write csv files containing all, train,
#' and test occurrences, as well as the background. All files will contain
#' additional columns with the values of the variables for each coordinate.
#' Default = FALSE
#' @param name.occ (character) name to be used for files with occurrence records.
#' Only one name is needed, a sufix will be added to represent all (_join),
#' _train, and _test records (e.g., "occurrences").
#' @param name.back name for the csv file containing background coordinates
#' (e.g., "background").
#' @param set.seed seed to be used when sampling background and splitting records.
#' Default = 1
#'
#' @usage
#' prepare_swd(occ, species, longitude, latitude,
#'             data.split.method = "random", train.proportion = 0.5,
#'             raster.layers, sample.size = 10000, save = FALSE,
#'             name.occ, name.back, set.seed = 1)
#' @export
#'
#' @examples
#' # data
#' occ <- read.csv(list.files(system.file("extdata", package = "kuenm"),
#'                            pattern = "sp_joint.csv", full.names = TRUE))
#' occ <- data.frame(Species = "A_americanum", occ)
#'
#' mvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                   pattern = "Mbio_", full.names = TRUE))
#'
#' # preparing swd
#' prep <- prepare_swd(occ, species = "Species", longitude = "Longitude",
#'                     latitude = "Latitude", raster.layers = mvars,
#'                     sample.size = 5000)

prepare_swd <- function(occ, species, longitude, latitude,
                        data.split.method = "random", train.proportion = 0.5,
                        raster.layers, sample.size = 10000, save = FALSE,
                        name.occ, name.back, set.seed = 1) {
  xy <- occ[, c(longitude, latitude)]
  xyval <- raster::extract(raster.layers, xy, cellnumbers = TRUE)
  xyras <- raster::xyFromCell(raster.layers, xyval[, 1])
  occ <- data.frame(occ[, species], xyras, xyval[, -1])
  colnames(occ)[1:3] <- c(species, longitude, latitude)

  back <- raster::rasterToPoints(raster.layers)
  set.seed(set.seed)
  back <- back[sample(nrow(back), sample.size), ]
  back <- data.frame(background = "background", back)
  names(back)[1:3] <- c("background", longitude, latitude)

  octi <- which(!paste(occ[, longitude], occ[, latitude]) %in%
                  paste(back[, longitude], back[, latitude]))
  octid <- occ[octi, ]
  names(octid)[1:3] <- c("background", longitude, latitude)
  octid$background <- "background"

  back <- rbind(octid, back)

  occ <- kuenm_occsplit(occ, train.proportion, data.split.method, save, name.occ)
  if (save == TRUE) {write.csv(back, file = name.back, row.names = FALSE)}

  occ$background <- back
  return(occ)
}
