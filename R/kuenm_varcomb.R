#' All potential combinations of a group of variables
#'
#' @description kuenm_varcomb creates multiple sets of variables by grouping them in all their potential combinations.
#'
#' @param var.dir (character) the name of the folder containing variables that will be combined.
#' @param out.dir (character) the name of the folder in which subfolders with distinct combinations of
#' variables will be written.
#' @param min.number (integer) the minimum number of variables per combination. This number must be > 1.
#' Default = 2.
#' @param in.format (character) format of variables in \code{var.dir}. Options are "ascii", "GTiff", and "EHdr" = bil.
#' Default = "ascii".
#' @param out.format (character) format of variables to be written in distinct sets inside \code{out.dir}.
#' Options are "ascii", "GTiff", and "EHdr" = bil. Default = "ascii".
#'
#' @return A list containing vectors of all the potential combinations of variables. In addition, a folder
#' named \code{out.dir} with subfolders in which distinct combinations of variables produced are written.
#'
#' @details
#' Sest of variables are written in the working directory and not retained as RasterStacks to avoid
#' problems related to RAM limitations.
#'
#' Time of processing will be reduced considerably if \code{in.format} and \code{out.format} coincide
#' because files will be copied and not loaded and written.
#'
#' @usage
#' kuenm_varcomb(var.dir, out.dir, min.number = 2, in.format = "ascii",
#'               out.format = "ascii")
#'
#' @export
#'
#' @examples
#' # This example depends on data stored in your directory
#' var_dir <- "Variables" # your directory with variables to be combined
#' out_dir <- "M_variables" # output directory to be created
#' min_n <- 2
#' in_format <- "ascii"
#' out_format <- "GTiff"
#'
#' comb <- kuenm_varcomb(var.dir = var_dir, out.dir = out_dir, min.number = min_n,
#'                       in.format = in_format, out.format = out_format)


kuenm_varcomb <- function(var.dir, out.dir, min.number = 2, in.format = "ascii",
                          out.format = "ascii") {

  # Setting things up
  if (min.number < 2) {
    stop("min.number must be an integer > 1.")
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

  # List variable names
  variables <- list.files(path = var.dir, pattern = patt)

  if (length(variables) == 0) {
    stop(paste("No variables with format", in.format, "were found in the directory", var.dir))
  }

  if (min.number > length(variables)) {
    stop("min.number must be < the total number of variables.")
  }

  # Generating all combinations of variable names
  var_combinations <- all_var_comb(variables, min.number)

  # Output directory (General)
  dir.create(out.dir)

  # Subdirectories for variable sets
  sub_paths <- paste(out.dir, paste("Set", 1:length(var_combinations), sep = "_"), sep = "/")

  # Telling users how many sets they will create
  cat("\nA total of", length(sub_paths), "sets of variables resulted from combinations of",
      length(variables), "variables will be written.\n")

  # Copying or writing variable sin new sets
  if (in.format == out.format) {
    # Copying variables
    if(.Platform$OS.type == "unix") {
      pb <- txtProgressBar(min = 0, max = length(sub_paths), style = 3) #progress bar
    } else {
      pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sub_paths), width = 300) #progress bar
    }

    for (k in 1:length(sub_paths)) {
      Sys.sleep(0.1)
      if(.Platform$OS.type == "unix") {
        setTxtProgressBar(pb, i)
      } else {
        setWinProgressBar(pb, k, title = paste(round(k / length(sub_paths) * 100, 2),
                                               paste("% of the process has finished")))
      }

      dir.create(sub_paths[k])
      vars_comb <- paste(var.dir, var_combinations[[k]], sep = "/")

      vars_set <- paste(sub_paths[k], var_combinations[[k]], sep = "/")

      file.copy(from = vars_comb, to = vars_set)
    }

    if(.Platform$OS.type != "unix") {
      suppressMessages(close(pb))
    }

  } else {
    # Formats
    if (out.format == "ascii") {
      patt1 <- ".asc"
    }
    if (out.format == "GTiff") {
      patt1 <- ".tif"
    }
    if (out.format == "EHdr") {
      patt1 <- ".bil"
    }

    # change format names
    var_combinations <- lapply(var_combinations, function(x) {gsub(patt, patt1, x)})

    # Preparing folders, variable combinations, and writing results
    vars_all <- raster::stack(paste(var.dir, variables, sep = "/"))

    if(.Platform$OS.type == "unix") {
      pb <- txtProgressBar(min = 0, max = length(sub_paths), style = 3) #progress bar
    } else {
      pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sub_paths), width = 300) #progress bar
    }

    for (k in 1:length(sub_paths)) {
      Sys.sleep(0.1)
      if(.Platform$OS.type == "unix") {
        setTxtProgressBar(pb, i)
      } else {
        setWinProgressBar(pb, k, title = paste(round(k / length(sub_paths) * 100, 2),
                                               paste("% of the process has finished")))
      }

      dir.create(sub_paths[k])
      vars_set <- vars_all[[gsub(patt1, "", var_combinations[[k]])]]

      for (l in 1:dim(vars_set)[3]) {
        raster::writeRaster(vars_set[[l]], filename = paste(sub_paths[k], var_combinations[[k]][l], sep = "/"),
                            format = out.format)
      }
    }

    if(.Platform$OS.type != "unix") {
      suppressMessages(close(pb))
    }
  }

  return(var_combinations)
}


