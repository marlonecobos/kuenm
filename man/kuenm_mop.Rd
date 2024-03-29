% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kuenm_mop.R
\name{kuenm_mop}
\alias{kuenm_mop}
\title{Extrapolation risk analysis for single comparisons}
\usage{
kuenm_mop(M.variables, G.stack, percent = 10, comp.each = 2000,
          parallel = FALSE, n.cores = NULL)
}
\arguments{
\item{M.variables}{a RasterStack of variables or a matrix with variables as columns
representing the calibration area. If matrix, columns must contain only
information for the variables to be used.}

\item{G.stack}{a RasterStack of variables representing the full area of interest, and areas
or scenarios to which models are transferred.}

\item{percent}{(numeric) percent of values sampled from te calibration region to calculate the MOP.}

\item{comp.each}{(numeric) compute distance matrix for a each fixed number of rows (default = 2000).}

\item{parallel}{(logical) if TRUE, calculations will be performed in parallel using \code{n.cores}
of the computer. This will demand more RAM and almost full use of the CPU; hence, its use
is more recommended in high-performance computers. Using this option will speed up the analyses.
Default = FALSE.}

\item{n.cores}{(numeric) number of cores to be used in parallel processing.
Default = NULL, in which case all CPU cores on current host - 1 will be used.}
}
\value{
A mobility-oriented parity RasterLayer where values of 0 represent strict extrapolation,
which means complete dissimilarity of environments between the calibration (M) or the background,
and the projection area (G).
}
\description{
kuenm_mop calculates a mobility-oriented parity layer by
comparing environmental values between the calibration area and the area or
scenario to which an ecological niche model is transferred.
}
\details{
The MOP is calculated following Owens et al.
(2013; \url{https://doi.org/10.1016/j.ecolmodel.2013.04.011}). This function is a modification
of the \code{\link[ENMGadgets]{MOP}} funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
}
\examples{
mvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
                                  pattern = "Mbio_", full.names = TRUE))
gvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
                                  pattern = "Gbio_", full.names = TRUE))
names(mvars) <- gsub("M", "", names(mvars))
names(gvars) <- names(mvars)

perc <- 5

mop <- kuenm_mop(M.variables = mvars, G.stack = gvars, percent = perc)
}
