#' Creation of Maxent models with selected parameters in SWD format
#'
#' @description kuenm_mod_swd creates and executes a batch file (bash for Unix)
#' for generating Maxent models using parameters previously selected with the
#' \code{\link{kuenm_cal_swd}} function.
#'
#' @param occ.joint (character) the name of csv file with training and testing
#' occurrences combined; columns must be: species, longitude, latitude, and two
#' or more columns representing distinct variables. See details in
#' \code{\link{kuenm_cal_swd}}.
#' @param back.dir (character) the name of the folder containing one or more csv
#' files representing one or more sets of predictor variables for a background.
#' Columns in background files must be: background, longitude, latitude, and two
#' or more columns representing distinct variables. See details in
#' \code{\link{kuenm_cal_swd}}.
#' @param out.eval (character) name of the folder where evaluation results
#' (from calibration) were written.
#' @param batch (character) name for the batch file (bash for Unix) with the code
#' to create final Maxent models.
#' @param rep.n (numeric) number of model replicates, default = 10.
#' @param rep.type (character) the replicate type; it can be: "Crossvalidate",
#' "Bootstrap", or "Subsample".
#' @param jackknife (logical) if TRUE, a jackknife process is performed while
#' runing Maxent models, default = FALSE.
#' @param max.memory (numeric) maximum memory (in megabytes) to be used by maxent
#' while creating the models. Default = 1000.
#' @param out.format (character) the model output format; it can be: "raw",
#' "logistic", "cloglog", or "cumulative".
#' @param project (logical) if TRUE, models will be projected to scenarios
#' in G.var.dir, default = FALSE.
#' @param G.var.dir (character) if project is TRUE, name of the folder containing
#' subfolders named as in Sets of \code{back.dir}, in which other subfolders where
#' variables of projection scenarios are placed.
#' @param ext.type (character) if project is TRUE, is the extrapolation type of
#' projections; can be: "all", "ext_clam", "ext", and "no_ext", default = "all".
#' ext = free extrapolation, ext_clam = extrapolation and clamping,
#' no_ext = no extrapolation, and all = all three options listed above.
#' @param write.mess (logical) if TRUE, grids of MESS analysis results will be
#' written, default = FALSE.
#' @param write.clamp (logical) if TRUE, a grid of the spatial distribution of
#' clamping will be written, default = FALSE.
#' @param maxent.path (character) the complete path were maxent.jar file is in
#' the computer.
#' @param args (character) additional arguments that can be passed to Maxent.
#' See the Maxent help for more information on how to write these arguments,
#' default = NULL. Note that some arguments cannot be changed here because they
#' are part of the parameters of the function already (e.g., "writemess").
#' See details for other options.
#' @param out.dir (character) name of the output directory to be created and
#' in which all model subdirectories will be created.
#' @param wait (logical) if TRUE, R will wait until all the Maxent models are
#' created. If FALSE the process of model creation will be performed separately
#' and R could be used at the same time. Default = FALSE
#' @param run (logical) if TRUE, the batch runs after its creation; if FALSE,
#' it will only be created and its running would be manual, default = TRUE.
#'
#' @return A folder named as out.dir with all the subfolders to save Maxent final
#' model results when running the .bat file (.sh for Unix). A batch file
#' (bash for Unix) for creating all the final Maxent models with their projections
#' if \code{project} = TRUE.
#'
#' @details
#' Same requirements regarding Java and maxent than in \code{\link{kuenm_cal_swd}}.
#'
#' The way to include further arguments is as follows:
#' args = "biasfile=COMPLETE_PATH\\bias.asc biastype=3" in windows,
#' or args = "biasfile=COMPLETE_PATH/bias.asc biastype=3" in Unix based systems.
#' If the path contains spaces the way to write it is:
#' args = "biasfile=\"COMPLETE PATH\\bias.asc\" biastype=3" in windows, or
#' args = "biasfile=\"COMPLETE PATH/bias.asc\" biastype=3" in Unix based systems.
#'
#' Other options that can be included in args are all "Flags" from the following
#' list:
#'
#' Flag | Abbrv | Type | Default | Meaning
#' - maximumbackground | MB | integer | 10000 | If the number of background points / grid cells is larger than this number, then this number of cells is chosen randomly for background points.
#' - togglelayertype | t | string | | Toggle continuous/categorical for environmental layers whose names begin with this prefix (default: all continuous).
#' - biasfile | | file | | Sampling is assumed to be biased according to the sampling distribution given in this grid file. Values in this file must not be zero or negative. MaxEnt will factor out the bias. We recomend to create this file as a kernell density of geographic points representing all localities were samplings of similar organisms have been performed (multiply this layer by 1000 and round it to reduce number of decimals). IMPORTANT: A biasfile must be included with its entire path, as indicated above above.
#' - biastype | | integer | | If biasfile is defined, this integer needs to be definef depending on the type of bias added. If the bias file is prepared as recomended, biastype=3.
#' - writebackgroundpredictions | | boolean | FALSE | Write .csv file with predictions at background points.
#' - maximumiterations | m | integer | 500 | Stop training after this many iterations of the optimization algorithm.
#' - convergencethreshold | c | double | 0.00001 | Stop training when the drop in log loss per iteration drops below this number.
#' - threads | | integer | 1 | Number of processor threads to use. Matching this number to the number of cores on your computer speeds up some operations, especially variable jackknifing.
#' - logfile | | string | maxent.log | File name to be used for writing debugging information about a run in output directory.
#' - cache | | boolean | TRUE | Make a .mxe cached version of ascii files, for faster access.
#' - defaultprevalence | | double | 0.5 | Default prevalence of the species: probability of presence at ordinary occurrence points. See Elith et al., Diversity and Distributions, 2011 for details.
#'
#' Other more advanced arguments are (use these ones only if you understand them completely):
#' - lq2lqptthreshold | | integer | 80 | Number of samples at which product and threshold features start being used.
#' - l2lqthreshold | | integer | 10 | Number of samples at which quadratic features start being used.
#' - hingethreshold | | integer | 15 | Number of samples at which hinge features start being used.
#' - beta_threshold | | double | -1 | Regularization parameter to be applied to all threshold features; negative value enables automatic setting.
#' - beta_categorical | | double | -1 | Regularization parameter to be applied to all categorical features; negative value enables automatic setting.
#' - beta_lqp | | double | -1 | Regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting.
#' - beta_hinge | | double | -1 | Regularization parameter to be applied to all hinge features; negative value enables automatic setting.
#'
#' @usage
#' kuenm_mod_swd(occ.joint, back.dir, out.eval, batch, rep.n = 10,
#'               rep.type = "Bootstrap", jackknife = FALSE,
#'               max.memory = 1000, out.format = "logistic",
#'               project = FALSE, G.var.dir, ext.type = "all",
#'               write.mess = FALSE, write.clamp = FALSE, maxent.path,
#'               args = NULL, out.dir, wait = FALSE, run = TRUE)
#'
#' @export

kuenm_mod_swd <- function(occ.joint, back.dir, out.eval, batch, rep.n = 10,
                          rep.type = "Bootstrap", jackknife = FALSE,
                          max.memory = 1000, out.format = "logistic",
                          project = FALSE, G.var.dir, ext.type = "all",
                          write.mess = FALSE, write.clamp = FALSE, maxent.path,
                          args = NULL, out.dir, wait = FALSE, run = TRUE) {

  #Checking potential issues
  if (!file.exists(occ.joint)) {
    stop(paste(occ.joint, "does not exist in the working directory, check file name",
               "\nor extension, example: species_joint.csv"))
  }
  if (missing(back.dir)) {
    stop("Argument back.dir is not defined.")
  }
  if (!dir.exists(back.dir)) {
    stop(paste(back.dir, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (length(list.files(back.dir, pattern = ".csv$")) == 0) {
    stop(paste(back.dir, "does not contain any csv file."))
  }
  if (missing(out.eval)) {
    stop(paste("Argument out.eval is not defined, it is necessary for reading selected",
               "\nsets of parameters."))
  }
  if (project == TRUE) {
    if (missing(G.var.dir)) {
      stop("Argument G.var.dir is not defined.")
    }
    if (!dir.exists(G.var.dir)) {
      stop(paste(G.var.dir, "does not exist in the working directory, check folder name",
                 "\nor its existence."))
    }
    if (length(list.dirs(G.var.dir, recursive = FALSE)) == 0) {
      stop(paste(G.var.dir, "does not contain any subdirectory with sets of projection variables;",
                 "\neach subdirectory inside", G.var.dir, "must containg at least one subdirectory",
                 "\nwith the projection variables"))
    }
  }
  if (missing(maxent.path)) {
    stop("Argument 'maxent.path' is not defined, it is necessary for executing Maxent.")
  }


  # Slash
  if(.Platform$OS.type == "unix") {sl <- "/"; dl <- "/"} else {sl <- "\\"; dl <- "\\\\"}

  # Data
  ## Data from selected models table
  sett <- read.csv(paste0(out.eval, "/selected_models.csv"), stringsAsFactors = FALSE)
  setts <- strsplit(sett[, 1], split = "_")

  ### RM and FC
  rm <- sapply(setts, function(x) {x[2]})
  f.clas <- sapply(setts, function(x) {x[4]})

  ### Background
  var.di <- sapply(setts, function(x) {paste(x[5:length(x)], collapse = "_")})
  var.dir <- paste0("\"", paste(gsub("/", dl, paste(getwd(), back.dir, sep = sl)),
                                paste0(var.di, ".csv"), sep = sl), "\"")

  # output directory
  dir.create(out.dir)
  out.dir <- gsub("/", dl, paste(getwd(), out.dir, sep = sl))

  # Defining maximum ram to be used (50% of free memory)
  ram <- paste0("-mx", max.memory, "m")

  #####
  #Maxent settings
  ##Environmental calibration variable sets
  env <- paste0("environmentallayers=", var.dir)

  ##Species occurrences
  samp <- paste0("samplesfile=",
                 gsub("/", dl, paste0("\"", paste(getwd(), occ.joint, sep = sl), "\"")))

  ##Feature classes combinations
  fea <- feature_classes(f.clas)

  ## Projection (G) variables folders and subfolders, extrapolation types, and writting clamp and MESS
  if(project == TRUE) {
    G.dir <- paste(gsub("/", dl, paste(getwd(), G.var.dir, sep = sl)),
                   var.di, sep = sl)
    G.dirs <- sapply(1:length(G.dir), function(x) {
      dires <- sapply(dir(G.dir[x]), function(y) {
        paste0("\"", paste(gsub("/", dl, paste(getwd(), G.var.dir, sep = sl)),
                           var.di[x], y, sep = sl), "\"")
      })
      paste0("projectionlayers=", paste(dires, collapse = ","))
    })

    mid.com <- ext_type(ext.type)[[1]]; ext.nam <- ext_type(ext.type)[[2]]

    w.mess <- ifelse(write.mess == FALSE, "writeclampgrid=false", "writeclampgrid=true")
    w.clamp <- ifelse(write.clamp == FALSE, "writemess=false", "writemess=true")

    mid.com <- paste(mid.com, w.mess, w.clamp, "responsecurves=true")
  }else {
    ext.nam <- NULL
    mid.com <- "extrapolate=false doclamp=false writeclampgrid=false writemess=false responsecurves=true"
  }

  ## Jackknife
  jack <- ifelse(jackknife == TRUE, "jackknife=true", "jackknife=false")

  ## Output format
  out <- paste0("outputformat=", out.format)

  ## Number of replicates
  rep <-  paste0("replicates=", rep.n)

  ##Replicate type
  rept <- paste0("replicatetype=", rep.type)

  # Fixed commands
  ## Intitial command
  in.comm <- paste("java", ram, paste("-jar", gsub("/", dl, paste0("\"", paste(maxent.path, "maxent.jar", sep = sl), "\""))))

  ## Autofeature
  a.fea <- "autofeature=false"

  ## Other maxent settings
  fin.com <- "warnings=false visible=false autorun\n"

  #####
  # Final part
  ## Creating subdirectories
  rtimes <- ifelse(is.null(ext.nam), 1, length(ext.nam))
  par <- rep(sett[, 1], each = rtimes)
  subdir <- paste0(par, ext.nam)
  subfol <- paste0("outputdirectory=", paste0("\"", out.dir, sl, subdir, "\""))
  di <- sapply(subdir, function(x) {dir.create(paste0(out.dir, sl, x))})

  ## Repeating things as needed
  env <- rep(env, each = rtimes)
  reg.m <- rep(paste0("betamultiplier=", rm), each = rtimes)
  fea <- rep(fea, each = rtimes)
  G.dirs <- rep(G.dirs, each = rtimes)

  ## Writing java code
  if(project == TRUE) {
    jmx <- paste(in.comm, env, samp, G.dirs, subfol, reg.m, a.fea, fea, rep, rept,
                 jack, out, mid.com, args, fin.com)
  } else {
    message("Argument 'project' is FALSE, not geogrpahic projections will be produced.\n")
    jmx <- paste(in.comm, env, samp, subfol, reg.m, a.fea, fea, rep, rept,
                 jack, out, mid.com, args, fin.com)
  }

  if(.Platform$OS.type == "unix") {
    cat(c("#! /bin/csh\n", jmx), file = paste0(batch, ".sh"))
  } else {
    cat(jmx, file = paste0(batch, ".bat"))
  }

  # Running models
  message("If asked, RUN as administrator.")
  run_maxent(batch, maxent.path, wait = wait)

  # Finishing
  if (wait == FALSE) {
    message(paste0("A total of ", length(jmx), " final models will be created."))
  } else {
    message(paste0("A total of ", length(jmx), " final models were created."))
  }
  message(paste0("Check your working directory!!!\t", getwd()))
}
