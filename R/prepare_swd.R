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
#' @param var.sets (character or list) if character the only option is "all_comb",
#' which will prepare the background to obtain all potential combinations of
#' variables considering the ones in \code{raster.layers}. The minimum number of
#' variables per set is defied by \code{min.number}. If list, a list
#' of character vectors with the names of the variables per each set. Names of
#' variables in sets must match names of layers in \code{raster.layers}.
#' The default (NULL) produces only one set of variables for the background.
#' @param min.number (numeric) minimum number of variables per set when option
#' "all_comb" is used in \code{var.sets}. Default = 2.
#' @param save (logical) whether or not to write csv files containing all, train,
#' and test occurrences, as well as the background. All files will contain
#' additional columns with the values of the variables for each coordinate.
#' Default = FALSE.
#' @param name.occ (character) name to be used for files with occurrence records.
#' Only one name is needed, a sufix will be added to represent all (_join),
#' _train, and _test records (e.g., "occurrences").
#' @param back.folder name for the csv file containing background coordinates
#' (e.g., "background").
#' @param set.seed seed to be used when sampling background and splitting records.
#' Default = 1
#'
#' @usage
#' prepare_swd(occ, species, longitude, latitude, data.split.method = "random",
#'             train.proportion = 0.5, raster.layers, sample.size = 10000,
#'             var.sets = NULL, min.number = 2, save = FALSE, name.occ,
#'             back.folder, set.seed = 1)
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
#' # preparing swd data one set of variables
#' prep <- prepare_swd(occ, species = "Species", longitude = "Longitude",
#'                     latitude = "Latitude", raster.layers = mvars,
#'                     sample.size = 5000)
#'
#' # various sets of variables
#' preps <- prepare_swd(occ, species = "Species", longitude = "Longitude",
#'                      latitude = "Latitude", raster.layers = mvars,
#'                      var.sets = "all_comb", min.number = 3,
#'                      sample.size = 5000)

prepare_swd <- function(occ, species, longitude, latitude,
                        data.split.method = "random", train.proportion = 0.5,
                        raster.layers, sample.size = 10000, var.sets = NULL,
                        min.number = 2, save = FALSE, name.occ, back.folder,
                        set.seed = 1) {
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
  bnames <- c("background", longitude, latitude)
  names(octid)[1:3] <- bnames
  octid$background <- "background"

  back <- rbind(octid, back)

  occ <- kuenm_occsplit(occ, train.proportion, data.split.method, save, name.occ)

  if (save == TRUE) {dir.create(back.folder)}
  if (!is.null(var.sets)) {
    if (class(var.sets)[1] %in% c("character", "list")) {
      if (class(var.sets)[1] == "character") {
        if (var.sets == "all_comb") {
          if (min.number == 1) {
            message("Minimum number of variables in background sets is 1, do not use product features.")
          }
          var_names <- colnames(back)[-(1:3)]
          var.sets <- all_var_comb(var_names, min.number)
        } else {
          warning("Argument 'var.sets' is not valid returning one set of background variables.")
        }
      } else {
        ls <- sapply(var.sets, length)
        if (any(ls == 1)) {
          message("Minimum number of variables in background sets is 1, do not use product features.")
        }
        names(var.sets) <- paste0("Set_", 1:length(var.sets))
      }
    } else {
      warning("Argument 'var.sets' is not valid returning one set of background variables.")
    }
    if (save == TRUE) {
      nambs <- names(var.sets)
      sv <- sapply(nambs, function(x) {
        nms <- c(bnames, var.sets[[x]])
        write.csv(back[, nms], file = paste0(back.folder, "/", x, ".csv"),
                  row.names = FALSE)
      })
    }
  } else {
    var.sets <- list(Set_1 = colnames(back)[-(1:3)])
    if (save == TRUE) {
      write.csv(back, file = paste0(back.folder, "/Set_1.csv"), row.names = FALSE)
    }
  }


  occ$background <- back
  occ$sets <- var.sets

  return(occ)
}

#' Helper to create all variable combinations
#' @param var.names (character) vector of variable names
#' @param min.number (numeric) minimum number of variables per set.
#' @export
#' @return A list of vectors containing variable names per set.
all_var_comb <- function(var.names, min.number = 2) {
  var_comb <- lapply(min.number:length(var.names), function(x) {
    comb <- combn(var.names, m = x)
    comb_vs <- lapply(1:dim(comb)[2], function(y) {comb[, y]})
  })

  var_combs <- do.call(c, var_comb)
  names(var_combs) <- paste0("Set_", 1:length(var_combs))
  return(var_combs)
}
