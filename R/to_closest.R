#' Move occurrences to closest pixel with environmental data
#'
#' @description to_closest helps in changing the longitude and latitude values
#' of occurrences with no environmental data, so they move to the closest pixel
#' of a raster layer that contains relevant information. This process prevents
#' NAs in future analyses.
#'
#' @param data data.frame or matrix of occurrence records. Columns must include
#' longitude and latitude. Other columns are optional and wont be changed.
#' @param longitude (character) name of the column with longitude data.
#' @param latitude (character) name of the column with latitude data.
#' @param raster.layer RasterLayer to be used as a reference.
#'
#' @return data.frame or matrix with the corrected coordinates.
#'
#' @export
#'
#' @examples
#' data <- read.csv(list.files(system.file("extdata", package = "kuenm"),
#'                             pattern = "sp_test.csv", full.names = TRUE))
#'
#' var <- raster::raster(list.files(system.file("extdata", package = "kuenm"),
#'                                  pattern = "Mbio_", full.names = TRUE)[1])
#'
#' raster::plot(var)
#'
#' out <- data.frame(as.character(data$species[1]), rbind(c(-103, 27), c(-90, 26.5),
#'                                                        c(-109, 40), c(-70, 41)))
#' colnames(out) <- colnames(data)
#'
#' data <- rbind.data.frame(data, out)
#'
#' points(data[, 2:3])
#'
#' data1 <- to_closest(data, longitude = "Longitude", latitude = "Latitude",
#'                     raster.layer = var)
#'
#' points(data1[, 2:3], col = "red")

to_closest <- function(data, longitude, latitude, raster.layer) {
  # detecting potential errors
  if (missing(data)) {
    stop("Argument data is necessary to perform the analysis")
  }
  if (missing(longitude)) {
    stop("Argument longitude is not defined.")
  }
  if (missing(latitude)) {
    stop("Argument latitude is not defined.")
  }
  if (missing(raster.layer)) {
    stop("Argument raster_layer is not defined.")
  }

  # preparing data
  xy <- data[, c(longitude, latitude)]
  vals <- raster::extract(raster.layer, xy)

  tomove <- which(is.na(vals))
  xyout <- xy[tomove, ]

  xyras <- raster::rasterToPoints(raster.layer)[, 1:2]

  dists <- raster::pointDistance(xyout, xyras, lonlat = TRUE)

  # running process
  cat("\nMoving occurrences to closest pixels:\n")
  for (i in 1:nrow(xyout)) {
    xyin <- xyras[dists[i, ] == min(dists[i, ])[1], ]
    if (class(xyin)[1] == "matrix") {
      xyin <- xyin[1, ]
    }
    data[tomove[i], longitude] <- xyin[1]
    data[tomove[i], latitude] <- xyin[2]
    cat("\tOccurrence", i, "of", nrow(xyout), "moved\n")
  }

  return(data)
}
