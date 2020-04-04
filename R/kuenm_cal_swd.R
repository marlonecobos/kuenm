#' Complete Maxent model calibration in SWD format
#'
#' @description kuenm_cal_swd performs the whole process of model calibration
#' (i.e., candidate model creation and evaluation) using Maxent in SWD format.
#' Models are created with multiple parameter combinations, including distinct
#' regularization multiplier values, various feature classes, and different sets
#' of environmental variables represented by csv files that contain the background.
#' Evaluation is done in terms of statistical significance (partial ROC),
#' prediction ability (omission rates), and model complexity (AICc). After
#' evaluation, this function selects the best models based on user-defined
#' criteria.
#'
#' @param occ.joint (character) the name of csv file with training and testing
#' occurrences combined; columns must be: species, longitude, latitude, and two
#' or more columns representing distinct variables. See details.
#' @param occ.tra (character) the name of the csv file with the training
#' occurrences; columns as in occ.joint.
#' @param occ.test (character) the name of the csv file with the evaluation
#' occurrences; columns as in occ.joint.
#' @param back.dir (character) the name of the folder containing one or more csv
#' files representing one or more sets of predictor variables for a background.
#' Columns in background files must be: background, longitude, latitude, and two
#' or more columns representing distinct variables. See details.
#' @param batch (character) name of the batch file (bash for Unix) with the
#' code to create all candidate models for calibration.
#' @param out.dir.models (character) name of the folder that will contain all
#' calibration model subfolders.
#' @param reg.mult (numeric vector) regularization multiplier(s) to be evaluated.
#' @param f.clas (character) feature clases can be selected from five different
#' combination sets or manually. Combination sets are: "all", "basic", "no.t.h",
#' "no.h", and "no.t". Default = "all". basic = "l", "lq", "lqp", "lqpt", "lqpth".
#' Combinations "no.t.h", "no.h", and "no.t", exclude t and/or h. See details for
#' all the available potential combinations of feature classes.
#' @param max.memory (numeric) maximum memory (in megabytes) to be used by maxent
#' while creating the models. Default = 1000.
#' @param args (character) additional arguments that can be passed to Maxent.
#' See the Maxent help for more information on how to write these arguments,
#' default = NULL. Note that some arguments cannot be changed here because they
#' are part of the parameters of the function already (e.g., "betamultiplier" or
#' "plots"). See details for other options.
#' @param maxent.path (character) the path were maxent.jar file is in your computer.
#' @param selection (character) model selection criterion, can be "OR_AICc",
#' "AICc", or "OR"; OR = omission rates. Default = "OR_AICc", which means that
#' among models that are statistically significant and that present omission
#' rates below the threshold, those with delta AICc up to 2 will be selected.
#' See details for other selection criteria.
#' @param threshold (numeric) the percentage of training data omission error
#' allowed (E); default = 5.
#' @param rand.percent (numeric) the percentage of data to be used for the
#' bootstraping process when calculating partial ROCs; default = 50.
#' @param iterations (numeric) the number of times that the bootstrap is going
#' to be repeated; default = 500.
#' @param kept (logical) if FALSE, all candidate models will be erased after
#' evaluation, default = TRUE.
#' @param out.dir.eval (character) name of the folder where evaluation results
#' will be written.
#'
#' @details
#' Java needs to be installed in the computer and maxent.jar needs to be in a
#' known place in the computer. Java can be obtained from https://java.com/es/download/manual.jsp.
#' Users of Linux and Mac need the entire Java Development Kit available in
#' http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html.
#' Maxent can be downloaded from https://biodiversityinformatics.amnh.org/open_source/maxent/
#'
#' Below all potential combinations of feature classes are shown. Manual selection
#' can be done by creating a vector of one or more of the combinations of this
#' list. l = linear, q = quadratic, p = product, t = threshold, and h = hinge.
#' "l", "q", "p", "t", "h", "lq", "lp", "lt", "lh", "qp", "qt", "qh", "pt", "ph",
#' "th", "lqp", "lqt", "lqh", "lpt", "lph", "qpt", "qph", "qth", "pth", "lqpt",
#' "lqph", "lqth", "lpth", and "lqpth".
#'
#' Other selecton criteria are described below: If "AICc" criterion is chosen,
#' all significant models with delta AICc up to 2 will be selected If "OR" is
#' chosen, the 10 first significant models with the lowest omission rates will
#' be selected.
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
#' @return
#' A folder named \code{out.dir.models} with all the subfolders to save Maxent
#' results when running the .bat file (.sh for Unix). A .bat file (.sh for Unix)
#' containing the java codes to run the calibration models, it will run auotmatically
#' or on some computers a dialog box will ask if running is allowed.
#'
#' A list with three dataframes containing results from the calibration process
#' and a scatterplot of all models based on the AICc values and omission rates.
#' In addition, a folder, in the working directory, containing a csv file with i
#' nformation about models meeting the user-defined selection criterion, another
#' csv file with a summary of the evaluation and selection process, an extra csv
#' file containing all the statistics of model performance (pROC, AICc, and
#' omission rates) for all candidate models, a png scatterplot of all models
#' based on the AICc values and rates, and an HTML file sumarizing all the i
#' nformation produced after evaluation for helping with further interpretation.
#'
#' @usage
#' kuenm_cal_swd(occ.joint, occ.tra, occ.test, back.dir, batch,
#'               out.dir.models, reg.mult, f.clas = "all",
#'               max.memory = 1000, args = NULL, maxent.path,
#'               selection = "OR_AICc", threshold = 5,
#'               rand.percent = 50, iterations = 500,
#'               kept = TRUE, out.dir.eval)
#'
#' @export

kuenm_cal_swd <- function(occ.joint, occ.tra, occ.test, back.dir, batch,
                          out.dir.models, reg.mult, f.clas = "all",
                          max.memory = 1000, args = NULL, maxent.path,
                          selection = "OR_AICc", threshold = 5,
                          rand.percent = 50, iterations = 500,
                          kept = TRUE, out.dir.eval) {

  # Slash
  if(.Platform$OS.type == "unix") {sl <- "/"; dl <- "/"} else {sl <- "\\"; dl <- "\\\\"}

  #####
  # candidate models

  # Data
  ## Environmental variables sets
  m <- dir(back.dir)
  ms <- paste(gsub("/", dl, paste(getwd(), back.dir, sep = sl)), sl, m, sep = "")
  env <- paste("environmentallayers=", paste("\"", ms, "\"", sep = ""), sep = "")
  m <- gsub(".csv$", "", m)

  ## Species occurrences
  oc <- occ.joint
  samp <- paste("samplesfile=", gsub("/", dl, paste("\"", paste(getwd(), oc, sep = sl),
                                                    "\"", sep = "")), sep = "")
  occ <- occ.tra
  samp1 <- paste("samplesfile=", gsub("/", dl, paste("\"", paste(getwd(), occ, sep = sl),
                                                     "\"", sep = "")), sep = "")

  # Maxent settings
  ## Feature classes combinations
  fea <- feature_classes(f.clas)

  # output directories
  dir.create(out.dir.models)
  out.dir <- gsub("/", dl, paste(getwd(), out.dir.models, sep = sl))

  # Getting ram to be used
  ram <- paste("-mx", max.memory, "m", sep = "")

  # Fixed commands
  ## Intitial command
  in.comm <- paste("java", ram, paste("-jar", gsub("/", dl, paste0("\"", paste(maxent.path, "maxent.jar", sep = sl), "\""))))

  ## Autofeature
  a.fea <- "autofeature=false"

  ## Other maxent settings
  fin.com <- "extrapolate=false doclamp=false replicates=1 replicatetype=Crossvalidate responsecurves=false jackknife=false plots=false pictures=false outputformat=raw warnings=false visible=false redoifexists autorun\n"
  fin.com1 <- "extrapolate=false doclamp=false replicates=1 replicatetype=Crossvalidate responsecurves=false jackknife=false plots=false pictures=false outputformat=logistic warnings=false visible=false redoifexists autorun\n"

  # Final set of calibration models
  ## preparin final arguments
  treg <- length(reg.mult); tfea <- length(fea); tenv <- length(env)
  total_comb <- treg * tfea * tenv
  repfea <- total_comb / tfea

  reg.mult <- rep(reg.mult, each = total_comb / treg)
  fea <- rep(rep(fea, repfea / tenv), each = tenv)
  m <- rep(m, total_comb / tenv)

  reg.m <- paste0("betamultiplier=", reg.mult)
  env <- rep(env, total_comb / tenv)

  ## creating subdirectories
  subdir <- paste("M", reg.mult, "F", names(fea), m, "all", sep = "_")
  subfol <- paste0("outputdirectory=", paste0("\"", out.dir, sl, subdir, "\""))
  di <- sapply(subdir, function(x) {dir.create(paste0(out.dir, sl, x))})

  subdir1 <- paste("M", reg.mult, "F", names(fea), m, "train", sep = "_")
  subfol1 <- paste0("outputdirectory=", paste0("\"", out.dir, sl, subdir1, "\""))
  di <- sapply(subdir1, function(x) {dir.create(paste0(out.dir, sl, x))})

  ## writing java code
  allc <- paste(in.comm, env, samp, subfol, reg.m, a.fea, fea, args, fin.com)
  trac <- paste(in.comm, env, samp1, subfol1, reg.m, a.fea, fea, args, fin.com1)
  jmx <- unlist(lapply(1:length(allc), function(x) {c(allc[x], trac[x])}))
  if(.Platform$OS.type == "unix") {
    cat(c("#! /bin/csh\n", jmx), file = paste0(batch, ".sh"))
  } else {
    cat(jmx, file = paste0(batch, ".bat"))
  }

  # running models
  message("If asked, RUN as administrator")
  run_maxent(batch, maxent.path)

  # candidate model messages
  message(paste0("\nA total of ", total_comb, " candidate models will be created"))


  #####
  # evaluation
  message("\nStarting evaluation process")

  # data
  ## For AICc
  raw_folders <- paste0(out.dir, sl, subdir)

  ## For pROC and omission rates
  log_folders <- paste0(out.dir, sl, subdir1)

  # evaluation process
  poa <- proc_or_aicc(occ.joint, occ.tra, occ.test,
                      raw_folders, log_folders, threshold, rand.percent,
                      iterations, kept)

  ## Erasing main folder of candidate models if kept = FALSE
  if(kept == FALSE) {
    unlink(out.dir.models, recursive = T)
    message("All candidate models were deleted")
  }

  # summary of results and model selection
  list_res <- summary_calibration(poa, selection)

  # writing results
  ## csv files
  message("\nWriting calibration results")
  dir.create(out.dir.eval)

  name <- paste0(out.dir.eval, "/calibration_results.csv")
  name0 <- paste0(out.dir.eval, "/calibration_stats.csv")
  name1 <- paste0(out.dir.eval, "/selected_models.csv")

  write.csv(list_res[[3]], file = name, row.names = FALSE)
  write.csv(list_res[[1]], file = name0, row.names = FALSE)
  write.csv(list_res[[2]], file = name1, row.names = FALSE)

  ## plot
  png(paste0(out.dir.eval, "/calibration_figure.png"), width = 80, height = 80,
      units = "mm", res = 600)
  par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.58)
  plot_proc_aicc(list_res)
  dev.off()

  ## writing the html file
  html_calibration(path = out.dir.eval, file.name = "calibration_results")

  # finishing evaluation
  message("\nProcess finished")
  message("A folder containing results of model calibration for ", total_comb,
          "\ncandidate models has been written")

  message("\nThe folder ", out.dir.eval, " contains:")
  message("\t-A html file and its dependencies that summarize all the results")
  message("\t-Two csv files with results and statistics from models calibration")
  if(selection == "OR_AICc"){
    message("\t-And an aditional csv file containing the models selected by OR and AICc\n")
  }
  if(selection == "AICc") {
    message("\t-And  an aditional csv file containing the models selected by AICc\n")
  }
  if(selection == "OR") {
    message("\t-And  an aditional csv file containing the models selected by OR\n")
  }

  message(paste0("Check your working directory!!!\t", getwd()))

  # returning results
  return(list_res[1:3])
}
