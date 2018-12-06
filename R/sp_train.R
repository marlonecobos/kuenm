#' A set of occurrence records for training candidate ecological niche models
#'
#' A dataset containing occurrence recods of a tick (*Amblyomma americanum*) across
#' North America, used to train candidate ecological niche models during calibration.
#'
#' @format A data frame with 89 rows and 2 columns.
#' \describe{
#'   \item{Longitude}{longitude, in decimal degrees.}
#'   \item{Latitude}{latitude, in decimal degrees.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#' @export
#'
#' @examples
#' occtr <- read.csv(list.files(system.file("extdata", package = "kuenm"),
#'                             pattern = "sp_train.csv", full.names = TRUE))
#'
#' head(occtr)
"sp_train"
