#' Evaluation of variable contribution to Maxent final models
#'
#' @description model_var_contrib helps to explore variable contribution of
#' Maxent models created as final models with the functions \code{\link{kuenm_mod}}
#' or \code{\link{kuenm_mod_swd}}. Variable contribution is measured based on
#' metrics of contribution percentage, permutation importance, and, if existent,
#' a jackknife analysis.
#'
#' @param fmod.dir (character) the  name of the folder in which final models are
#' (e.g., the output folder after using the \code{\link{kuenm_mod}}) function.
#' It is important to have only the folders containing the models in this
#' directory. It can be only one folder or multiple subfolders containing models
#' for the same species, created with distinct parameter settings. If models were
#' projected, and the distinct types of extrapolation were used, the name of the
#' folders contained in this directory should include a pattern describing the
#' type of extrapolation used (e.g., "EC" for extrapolation and clamping in
#' Maxent).
#' @param model_name (character) pattern to be searched when finding the model of
#' interest. This pattern does not include the pattern of \code{ext.type}. By
#' default, NULL, all models are considered.
#' @param project (logical) if TRUE, it is assumed that models were projected to
#' other scenarios (this must be always true if models were produced in SWD
#' format).
#' @param ext.type (character) vector of pattern(s) to be searched in the
#' folders inside \code{fmod.dir} that identify the extrapolation type(s) of
#' model projections of interest (e.g., "E", "EC", "NE", or a vector of more
#' than one of them). Ignored if \code{project} = FALSE.
#'
#' @return
#' A list with results of variable contribution, permutation importance, and
#' jackknife results. If multiple models are evaluated, a nested list with results
#' for all models is returned.
#'
#' @details
#' When models are created with replicates, the values returned correspond to the
#' average of such replicates.
#'
#' @usage
#' model_var_contrib(fmod.dir, model_name = NULL, project, ext.type)
#'
#' @export


model_var_contrib <- function(fmod.dir, model_name = NULL, project, ext.type) {
  # tests
  if (missing(fmod.dir)) {
    stop("Argument fmod.dir needs to be defined.")
  }
  if (!dir.exists(fmod.dir)) {
    stop(paste(fmod.dir, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (missing(project)) {
    stop("Argument project needs to be defined.")
  } else {
    if (project == TRUE) {
      if (missing(ext.type)) {
        stop("Argument ext.type needs to be defined.")
      }
    }
  }

  # Folders
  if (project == FALSE) {
    if (!is.null(model_name)) {
      parameters <- list(dir(fmod.dir, pattern = model_name[1],
                             full.names = TRUE))
      models <- list(dir(fmod.dir, pattern = model_name[1]))
    } else {
      parameters <- list(list.dirs(fmod.dir, recursive = FALSE))
      models <- list.dirs(fmod.dir, full.names = FALSE, recursive = FALSE)
    }
  } else {
    if (!is.null(model_name)) {
      parameters <- lapply(ext.type, function(i) {
        dir(fmod.dir, pattern = paste0(model_name, "_", i, "$"),
            full.names = TRUE, recursive = FALSE)
      })
      models <- lapply(ext.type, function(i) {
        dir(fmod.dir, pattern = paste0(model_name, "_", i, "$"),
            full.names = FALSE, recursive = FALSE)
      })
    } else {
      parameters <- lapply(ext.type, function(i) {
        dir(fmod.dir, pattern = paste0("_", i, "$"),
            full.names = TRUE, recursive = FALSE)
      })
      models <- lapply(ext.type, function(i) {
        dir(fmod.dir, pattern = paste0("_", i, "$"),
            full.names = FALSE, recursive = FALSE)
      })
    }
  }


  # preparing results
  var_cont_res <- lapply(1:length(parameters), function(x) {
    re <- lapply(1:length(parameters[[x]]), function(y) {
      allres <- read.csv(paste0(parameters[[x]][[y]], "/maxentResults.csv"))
      cols <- colnames(allres)
      nro <- nrow(allres)

      ## relevant columns
      colcont <- grep("contribution", cols)
      colperm <- grep("permutation.importance", cols)
      colgain <- grep("Regularized.training.gain", cols)
      colwith <- grep("gain.with.only", cols)
      colwout <- grep("gain.without", cols)

      ## relevant values
      varnamesmx <- gsub(".contribution", "", cols[colcont])
      contrib <- data.frame(Variable = varnamesmx,
                            Contribution = unlist(allres[nro, colcont]))
      rownames(contrib) <- NULL

      permimp <- data.frame(Variable = varnamesmx,
                            Permutation_importance = unlist(allres[nro, colperm]))
      rownames(permimp) <- NULL

      if (length(colwith) > 0) {
        jackkni <- data.frame(Variable = varnamesmx,
                              Training_gain_with = unlist(allres[nro, colwith]),
                              Training_gain_without = unlist(allres[nro, colwout]))
        rownames(jackkni) <- NULL

        rtg <- allres[nro, colgain]
      } else {
        jackkni <- NULL
        rtg <- NULL
      }

      list(Contribution = contrib, Permutation_importance = permimp,
           Jackknife_results = list(Regularized_training_gain_model = c(rtg),
                                    Training_gain_with_without = jackkni))
    })

    names(re) <- models[[x]]
    re
  })

  if (project == FALSE) {
    var_cont_res <- var_cont_res[[1]]
  } else {
    names(var_cont_res) <- ext.type
  }

  # results
  return(var_cont_res)
}
