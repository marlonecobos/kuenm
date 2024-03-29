% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SWD_analysis_helpers.R
\name{aicc}
\alias{aicc}
\title{AICc calculation of Maxent SWD predictions}
\usage{
aicc(occ, prediction, npar)
}
\arguments{
\item{occ}{matrix or data.frame with coordinates of the occurrences used to
create the model (raster) to be evaluated; columns must be: longitude and
latitude.}

\item{prediction}{matrix or data.frame of longitude and latitude coordinates,
and Maxent Raw predictions obtained using the SWD format in Maxent.
Coordinates in this prediction must include the ones in \code{occ}}

\item{npar}{(numeric) number of parameters of the model. Use function
\code{\link{n_par}} to obtain number of parameters in the model from
the lambdas file.}
}
\value{
A data.frame containing values of AICc, delta AICc, weight of AICc, and
number of parameters. The number of rows of the data.frame corresponds to
the number of models evaluated.
}
\description{
aicc calculates the Akaike information criterion corrected for
small sample sizes (AICc) for predictions produced with Maxent.
}
\details{
Calculations are done following
\href{https://doi.org/10.1890/10-1171.1}{Warren and Seifert (2011)}.
}
