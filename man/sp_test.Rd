% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_documentation.R
\name{sp_test}
\alias{sp_test}
\title{A set of occurrence records to test candidate ecological niche models}
\format{
A data frame with 89 rows and 2 columns.
\describe{
\item{Longitude}{longitude, in decimal degrees.}
\item{Latitude}{latitude, in decimal degrees.}
}
}
\source{
\url{https://kuscholarworks.ku.edu/handle/1808/26376}
}
\description{
A data.frame containing occurrence records of a tick (\emph{Amblyomma americanum})
in North America, used to test candidate models during calibration.
}
\examples{
data("sp_test", package = "kuenm")

head(sp_test)
}
