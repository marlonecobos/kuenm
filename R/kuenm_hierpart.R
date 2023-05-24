#' Hierarchical partition of the variance coming from distinct sources in ENMs
#'
#' @description kuenm_hierpart has been deprecated for the moment. A future,
#' more complete version will be available soon.
#'
#' @param sp.name (character) name of the species. This name must be the one
#' that appears as part of the raster file of each model replicate. If results
#' are from Maxent, this is the name that is in the first column of the csv
#' containing species occurrence data (species) but spaces replaced by "_".
#' @param fmod.dir (character) name of the folder where all models are (e.g.,
#' the output folder after using the \code{\link{kuenm_mod}}) function.
#' @param is.swd (logical) whether model calibration and final models were
#' produced using SWD format.
#' @param format (character) format of model raster files. Options are "ascii",
#' "GTiff", and "EHdr" = bil. Default = "ascii".
#' @param replicated (logical) whether or not models were created with
#' replicates.
#' @param project (logical) if TRUE, it is assumed that models were projected
#' to other scenarios. These scenarios can be current (projections in space),
#' and/or future or past (projections in time).
#' @param current (character) pattern to look for when defining which is the
#' raster file representing current projections. If NULL, results will be
#' produced for the area of calibration, and if any of \code{time.periods},
#' \code{clim.models}, or \code{emi.scenarios} is defined, results will be
#' be produced for these variance sources as well.
#' @param time.periods (character) pattern to be searched to identify model
#' projections to distinct time periods. If NULL, the default, it is assumed
#' that only one time period was considered.
#' @param emi.scenarios (character) pattern to be searched to identify
#' distinct emission scenarios (e.g., "recp45"). If NULL, the default, it is
#' assumed that only one emission scenario was used. Therefore, this source of
#' variation will not be considered.
#' @param clim.models (character) names that identify climatic models used for
#' project ENMs. If NULL, the default, it is assumed that only one climate model
#' was used. Therefore, this source of variation will not be considered.
#' @param ext.type (character) pattern(s) to be searched in the folders inside
#' \code{fmod.dir} that identify the extrapolation type(s) used in model
#' projections. This pattern(s) needs to be clearly distinguishable from the
#' other parts of the name of the folder name containing the model. For instance,
#' "EC" will be the patter that denotes extrapolation and clamping in the folder
#' named "M_0.1_F_l_set1_EC".
#' @param iterations (numeric) number of iterations to be performed in the
#' hierarchical partitioning analysis. Default = 100.
#' @param sample.size (numeric) number of pixels to be sampled per each model.
#' Default = 1000. Increasing this number is recommended when the number of
#' models and the computer features allow it.
#' @param set.seed (numeric) initial seed to be set before running analysis.
#' @param keep.tables (logical) if TRUE, tables that are written in
#' \code{out.dir} for each iteration of the hierarchical partitioning analyses
#' are kept. Default = FALSE.
#' @param factors.col a vector of colors for the bars to be plotted; if not
#' defined, a gray color palette is used.
#' @param out.dir (character) name of the output directory to be created
#' where results of the hierarchical partitioning analysis will be written.
#' @param verbose (logical) whether to print messages; default = TRUE.
#'
#' @return
#' The function returns a list containing the summary of total effects of
#' factors on variance contained in the models (mean and confidence intervals
#' of total effects). A plot of these values is also returned.
#'
#' Other results are written in \code{out.dir}. Folders named Variation or
#' HP_results_(EC, NE, and/or E, depending on \code{ext.type}) containing
#' csv files with the results of the hierarchical partitioning analyses an a
#' plot summarizing the total effects of the sources of variation on the
#' variance in the models.
#'
#' @details
#' If the length of any of the potential sources of variation is equal to one
#' (e.g., only one parameter, or only one climate model), this source of
#' variation will not be considered.
#'
#' Users must be specific when defining the patterns that the function will
#' search for. These patterns must be part of the raster file names of the
#' models so the function can locate each file without problems.
#'
#' Error whiskers in resulting plots represent the 95% Confidence Interval of
#' the mean. This interval is calculated using a bootstrap approach.
#'
#' @usage
#' kuenm_hierpart(sp.name, fmod.dir, is.swd, format = "ascii", replicated, project,
#'                current = NULL, time.periods = NULL, emi.scenarios = NULL,
#'                clim.models = NULL, ext.type, iterations = 100,
#'                sample.size = 1000, set.seed = 1, keep.tables = FALSE,
#'                factors.col = NULL, out.dir, verbose = TRUE)
#'
#' @export


kuenm_hierpart <- function(sp.name, fmod.dir, is.swd, format = "ascii", replicated,
                           project, current = NULL, time.periods = NULL,
                           emi.scenarios = NULL, clim.models = NULL, ext.type,
                           iterations = 100, sample.size = 1000, set.seed = 1,
                           keep.tables = FALSE, factors.col = NULL,
                           out.dir, verbose = TRUE) {

  stop("The function 'kuenm_hierpart' has been excluded from 'kuenm' for the moment.")
}

