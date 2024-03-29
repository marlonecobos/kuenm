% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kuenm_rpca.R
\name{kuenm_rpca}
\alias{kuenm_rpca}
\title{Principal componens for raster layers and projections}
\usage{
kuenm_rpca(variables, in.format, var.scale = TRUE, write.result = TRUE,
           out.format = "GTiff", out.dir = "PCA_results", project = FALSE,
           proj.vars, n.pcs)
}
\arguments{
\item{variables}{(character or RasterStack) if character, name of the folder where raster layers are located.
If RasterStack, stack of raster layers to be used in principal component analyses.}

\item{in.format}{(character) valid only if \code{variables} is character. Format of variables in the directory.
Options are "ascii", "GTiff", and "EHdr" = bil.}

\item{var.scale}{(logical) wheter or not to scale variables before performing principal component
analyses. Default = TRUE.}

\item{write.result}{(logical) whether or not to write PCA results and raster layers (PCs) in \code{out.dir}.}

\item{out.format}{(character) if \code{write.result} = TRUE, format of variables to be written in distinct
sets inside \code{out.dir}. Options are "ascii", "GTiff", and "EHdr" = bil. Default = "GTiff".}

\item{out.dir}{(character) valid if \code{write.result} = TRUE. Name of the folder to be created to save the
results of the analyses. Default = "PCA_results".}

\item{project}{(logical) whether or not to project the species niche to other scenario(s).
If TRUE, argument \code{proj.variables} needs to be defined. Default = FALSE.}

\item{proj.vars}{(character or RasterStack) if character, name of the folder where subfolders with environmental
variables of scenarios for projections are (useful if multiple projections are needed). If RasterStack, object
containing stacked variables of only one projection scenario. Variables must correspond with variables in \code{vars.folder}
(i.e., their names must correspond but they should represent conditions in other scenario).}

\item{n.pcs}{(numeric) number of principal components to be returned as rasters. By default all principal
components are returned as RasterLayers.}
}
\value{
A list containing PCA loadings and PCA summary as matrices, as well as one or multiple (if projected) RasterStacks
of principal components.

If \code{write.result} = TRUE, all results are written in \code{out.dir}.
}
\description{
kuenm_rpca performs a principal component analysis with a set of variables and
produces raster layers of them. If needed the pricipal components are projected to other
scenarios.
}
\details{
If \code{var.scale} = TRUE, variables are centered to cero and scaled using \code{\link[base]{scale}}.
}
\examples{
# Data
variab <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
                                   pattern = "Mbio_", full.names = TRUE))
names(variab) <- paste0("bio_", c(1, 12, 15, 17))

proj_var <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
                                     pattern = "Gbio_", full.names = TRUE))
names(proj_var) <- paste0("bio_", c(1, 12, 15, 17))

# Example with no projection
npcs <- 3

rpca <- kuenm_rpca(variables = variab, var.scale = TRUE, write.result = FALSE,
                   n.pcs = npcs)

# Example with projection
project <- TRUE

rpca1 <- kuenm_rpca(variables = variab, var.scale = TRUE, write.result = FALSE,
                    project = project, proj.vars = proj_var, n.pcs = npcs)
}
