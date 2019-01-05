#' A set of occurrence records for ecological niche models
#'
#' A dataset containing occurrence recods of a tick (*Amblyomma americanum*) across
#' North America. The dataset combines records for training and testing.
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
#' occjt <- read.csv(list.files(system.file("extdata", package = "kuenm"),
#'                             pattern = "sp_joint.csv", full.names = TRUE))
#'
#' head(occjt)
NULL
