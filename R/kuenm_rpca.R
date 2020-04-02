#' Principal componens for raster layers and projections
#'
#' @description kuenm_rpca performs a principal component analysis with a set of variables and
#' produces raster layers of them. If needed the pricipal components are projected to other
#' scenarios.
#'
#' @param variables (character or RasterStack) if character, name of the folder where raster layers are located.
#' If RasterStack, stack of raster layers to be used in principal component analyses.
#' @param in.format (character) valid only if \code{variables} is character. Format of variables in the directory.
#' Options are "ascii", "GTiff", and "EHdr" = bil.
#' @param var.scale (logical) wheter or not to scale variables before performing principal component
#' analyses. Default = TRUE.
#' @param write.result (logical) whether or not to write PCA results and raster layers (PCs) in \code{out.dir}.
#' @param out.dir (character) valid if \code{write.result} = TRUE. Name of the folder to be created to save the
#' results of the analyses. Default = "PCA_results".
#' @param out.format (character) if \code{write.result} = TRUE, format of variables to be written in distinct
#' sets inside \code{out.dir}. Options are "ascii", "GTiff", and "EHdr" = bil. Default = "GTiff".
#' @param project (logical) whether or not to project the species niche to other scenario(s).
#' If TRUE, argument \code{proj.variables} needs to be defined. Default = FALSE.
#' @param proj.vars (character or RasterStack) if character, name of the folder where subfolders with environmental
#' variables of scenarios for projections are (useful if multiple projections are needed). If RasterStack, object
#' containing stacked variables of only one projection scenario. Variables must correspond with variables in \code{vars.folder}
#' (i.e., their names must correspond but they should represent conditions in other scenario).
#' @param n.pcs (numeric) number of principal components to be returned as rasters. By default all principal
#' components are returned as RasterLayers.
#'
#' @return
#' A list containing PCA loadings and PCA summary as matrices, as well as one or multiple (if projected) RasterStacks
#' of principal components.
#'
#' If \code{write.result} = TRUE, all results are written in \code{out.dir}.
#'
#' @details
#' If \code{var.scale} = TRUE, variables are centered to cero and scaled using \code{\link[base]{scale}}.
#'
#' @usage
#' kuenm_rpca(variables, in.format, var.scale = TRUE, write.result = TRUE,
#'            out.format = "GTiff", out.dir = "PCA_results", project = FALSE,
#'            proj.vars, n.pcs)
#'
#' @export
#'
#' @examples
#' # Data
#' variab <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                    pattern = "Mbio_", full.names = TRUE))
#' names(variab) <- paste0("bio_", c(1, 12, 15, 17))
#'
#' proj_var <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                      pattern = "Gbio_", full.names = TRUE))
#' names(proj_var) <- paste0("bio_", c(1, 12, 15, 17))
#'
#' # Example with no projection
#' npcs <- 3
#'
#' rpca <- kuenm_rpca(variables = variab, var.scale = TRUE, write.result = FALSE, n.pcs = npcs)
#'
#' # Example with projection
#' project <- TRUE
#'
#' rpca1 <- kuenm_rpca(variables = variab, var.scale = TRUE, write.result = FALSE, project = project,
#'                     proj.vars = proj_var, n.pcs = npcs)


kuenm_rpca <- function(variables, in.format, var.scale = TRUE, write.result = TRUE, out.format = "GTiff",
                       out.dir = "PCA_results", project = FALSE, proj.vars, n.pcs) {

  # testing potential errors
  if (missing(variables)) {
    stop("Argument variables must be defined. See functions help.")
  }
  if (project == TRUE) {
    if (missing(proj.vars)) {
      stop("If projections are needed, argument proj.vars must be defined. See functions help.")
    }
  }

  # formatting
  if (class(variables)[1] == "character") {
    if (missing(in.format)) {
      stop("Argument variables is a character, in.format needs to be defined.")
    }
    if (in.format == "ascii") {
      patt <- ".asc$"
    }
    if (in.format == "GTiff") {
      patt <- ".tif$"
    }
    if (in.format == "EHdr") {
      patt <- ".bil$"
    }
  }

  if (!missing(write.result)) {
    if (out.format == "ascii") {
      patt1 <- ".asc"
    }
    if (out.format == "GTiff") {
      patt1 <- ".tif"
    }
    if (out.format == "EHdr") {
      patt1 <- ".bil"
    }
  }

  # reading variables
  if (class(variables)[1] == "character") {
    var <- list.files(variables, pattern = patt, full.names = TRUE)
    variables <- raster::stack(var)
  }

  var_points <- na.omit(raster::values(variables))

  # pca analyses
  if (var.scale == TRUE) {
    pca <- prcomp(var_points, center = TRUE, scale = TRUE)
  } else {
    pca <- prcomp(var_points, center = TRUE, scale = FALSE)
  }

  scores <- pca$x

  if (missing(n.pcs)) {
    n.pcs <- length(var)
  }

  pcras <- list()

  if (write.result == TRUE) {
    cat("\nWriting raster PCs in Output folder, please wait...\n")
    dir.create(out.dir)

    pca_fol <- paste(out.dir, "Initial", sep = "/")
    dir.create(pca_fol)
  }

  for (i in 1:n.pcs) {
    pcra <- variables[[1]]
    pcra[!is.na(raster::values(pcra))] <- scores[, i]

    if (write.result == TRUE) {
      filenam <- paste(pca_fol, "/pc_", i, patt1, sep = "")
      raster::writeRaster(pcra, filenam, format = out.format)
    }

    pcras[[i]] <- pcra
  }

  pcras <- do.call(raster::stack, pcras)
  names(pcras) <- paste0("pc_", 1:dim(pcras)[3])

  StdDev <- pca$sdev
  VarExp <- pca$sdev^2/sum(pca$sdev^2)
  CumVar <- cumsum(VarExp)
  SumPCAMat <- rbind(StdDev, VarExp, CumVar)
  colnames(SumPCAMat) <- paste("PC", seq(1, length(StdDev)), sep = "")
  row.names(SumPCAMat) <- c("Standard deviation", "Proportion of Variance",
                            "Cumulative Proportion")

  if (write.result == TRUE) {
    sink(paste(paste(pca_fol, "pca_results.txt", sep = "/")))
    cat("Principal component analysis results\n")
    cat("\nPCA loadings\n")
    print(pca$rotation)

    cat("\n\nPCA summary\n")
    print(SumPCAMat)
    sink()
  }

  # pca results to be returned
  loadings <- pca$rotation
  respca <- SumPCAMat

  # projecting PCs
  if (project == TRUE) {
    ppcrass <- list()

    if (write.result == TRUE) {
      cat("\nProjecting and writing projected raster PCs in Output folder, please wait...\n")
    } else {
      cat("\nProjecting raster PCs\n")
    }

    if (class(proj.vars)[1] == "character") {
      proj_dirs <- list.dirs(proj.vars, recursive = FALSE)
      proj_names <- list.dirs(proj.vars, recursive = FALSE, full.names = FALSE)
      fol_names <- paste(out.dir, proj_names, sep = "/")
    }
    if (class(proj.vars)[1] %in% c("RasterStack", "RasterBrick")) {
      proj_dirs <- "projection"
      proj_names <- "Projected_PCs"
      fol_names <- paste(out.dir, proj_names, sep = "/")
    }


    for (h in 1:length(proj_dirs)) {
      if (class(proj.vars)[1] == "character") {
        pvar <- list.files(proj_dirs[h], pattern = patt, full.names = TRUE)
        p_stack <- raster::stack(pvar)
      }
      if (class(proj.vars)[1] %in% c("RasterStack", "RasterBrick")) {
        p_stack <- proj.vars
      }
      if (write.result == TRUE) {
        dir.create(fol_names[h])
      }

      ppcras <- list()

      p_stackp <- na.omit(raster::values(p_stack))
      colnames(p_stackp) <- names(pca[[4]])
      p_pcs <- predict(pca, newdata = p_stackp)

      for (i in 1:n.pcs) {
        pcra <- p_stack[[1]]
        pcra[!is.na(raster::values(pcra))] <- p_pcs[, i]

        if (write.result == TRUE) {
          filenam <- paste(fol_names[h], "/pc_", i, patt1, sep = "")
          raster::writeRaster(pcra, filenam, format = out.format)
        }

        ppcras[[i]] <- pcra
      }

      ppcrass[[h]] <- do.call(raster::stack, ppcras)
      names(ppcrass[[h]]) <- paste0("pc_", 1:dim(ppcrass[[h]])[3])
    }

    names(ppcrass) <- paste("PCRasters", proj_names, sep = "_")
  }

  if (project == TRUE) {
    results <- c(list(loadings, respca, pcras), ppcrass)
    names(results)[1:3] <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
  }else {
    results <- list(loadings, respca, pcras)
    names(results) <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
  }

  if (write.result == TRUE) {
    cat("\nRaster PCA finished. Check your output directory", paste(getwd(), out.dir, sep = "/"), "\n")
  }

  return(results)
}
