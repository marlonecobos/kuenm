#' Evaluation and plot of variable contribution to Maxent models
#'
#' @description Exploration of variable contribution to Maxent models based on
#' metrics of contribution percentage, permutation importance, and a jackknife
#' analysis.
#'
#' @param occ a data.frame with occurrence records. If \code{M_variables} is of
#' class "RasterStack", columns must be (in that order): Species, Longitude,
#' Latitude. If \code{M_variables} is of class "data.frame" (SWD), columns
#' containing the values of the environmental variables in \code{M_variables}
#' must be added.
#' @param M_variables RasterStack or data.frame. If RasterStack, object containing
#' the variables to be used for modeling. If data.frame (SWD), columns must be:
#' "background", "longitude", "latitude", and a column for each of the variables
#' of interest. If data.frame, \code{occ} must be a data.frame.
#' @param maxent.path (character) the path were maxent.jar file is in your computer.
#' @param reg.mult (numeric vector) regularization multiplier(s) to be evaluated.
#' @param f.clas (character) feature classes can be selected from five different
#' combination sets or manually. Combination sets are: "all", "basic", "no.t.h",
#' "no.h", and "no.t". Default = "all". basic = "l", "lq", "lqp", "lqpt", "lqpth".
#' Combinations "no.t.h", "no.h", and "no.t", exclude t and/or h. See details for
#' all the available potential combinations of feature classes.
#' @param max.memory (numeric) maximum memory (in megabytes) to be used by Maxent
#' while creating the models. Default = 1000.
#' @param args (character) additional arguments that can be passed to Maxent.
#' See the Maxent help for more information on how to write these arguments,
#' default = NULL. Note that some arguments cannot be changed here because they
#' are part of the parameters of the function already.
#' See details for other options.
#' @param sample.size (numeric) number of points to represent the background for
#' the model. Default = 10000
#' @param plot (logical) whether to produce a default plot of results.
#' @param verbose (logical) whether to print messages; default = TRUE.
#'
#' @return
#' A list with results of variable contribution, permutation importance, and
#' jackknife results. If \code{plot} = TRUE results are plotted as horizontal
#' bars using \code{\link{plot_contribution}}.
#'
#' @details
#' All potential combinations of feature classes
#' (l = linear, q = quadratic, p = product, t = threshold, and h = hinge) are:
#' "l", "q", "p", "t", "h", "lq", "lp", "lt", "lh", "qp", "qt", "qh", "pt", "ph",
#' "th", "lqp", "lqt", "lqh", "lpt", "lph", "lth", "qpt", "qph", "qth", "pth",
#' "lqpt", "lqph", "lqth", "lpth", "qpth", and "lqpth".
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
#' explore_var_contrib(occ, M_variables, maxent.path, reg.mult = 1,
#'                     f.clas = NULL, max.memory = 1000, args = NULL,
#'                     plot = TRUE, verbose = TRUE)
#'
#' @export
#'
#' @importFrom utils write.csv read.csv
#'
#' @examples
#' # data
#' data("sp_joint", package = "kuenm")
#' occ <- data.frame(Species = "A_americanum", sp_joint)
#'
#' mvars <- raster::stack(list.files(system.file("extdata", package = "kuenm"),
#'                                   pattern = "Mbio_", full.names = TRUE))
#' # analysis
#' var_cont <- explore_var_contrib(occ = occ, M_variables = mvars,
#'                                 maxent.path = "C:/Maxent/3.4.1", plot = FALSE)
#'
#' # plot
#' plot_contribution(var_cont)

explore_var_contrib <- function(occ, M_variables, maxent.path, reg.mult = 1,
                                f.clas = NULL, max.memory = 1000, args = NULL,
                                sample.size = 10000, plot = TRUE, verbose = TRUE) {
  # initial tests
  if (missing(occ)) {
    stop("Argument 'occ' must be defined, see function's help")
  }
  if (missing(M_variables)) {
    stop("Argument 'M_variables' must be defined, see function's help")
  }
  if (missing(maxent.path)) {
    stop("Argument 'maxent.path' must be defined, see function's help")
  }
  if (!class(M_variables)[1] %in% c("RasterLayer", "RasterStack", "RasterBrick", "data.frame")) {
    stop("'M_variables' must be of class RasterLayer, RasterStack, or data.frame")
  }
  if (!class(occ)[1] %in% c("data.frame")) {
    stop("'occ' must be a data.frame")
  }

  # preparing data
  out.dir <- file.path(tempdir(), "jack_maxent")
  if (dir.exists(out.dir)) {unlink(out.dir, recursive = TRUE)}
  dir.create(out.dir)

  ## occs
  occs <- paste0(out.dir, "/occ_joint.csv")

  ## variables
  dirM <- paste0(out.dir, "/background/Set_1.csv")
  spc <- colnames(occ)[1]
  loc <- colnames(occ)[2]
  lac <- colnames(occ)[3]

  if (class(M_variables)[1] == "RasterStack") {
    pr <- prepare_swd(occ = occ, species = spc, longitude = loc, latitude = lac,
                      raster.layers = M_variables, sample.size = sample.size,
                      save = TRUE, name.occ = paste0(out.dir, "/occ"),
                      back.folder = paste0(out.dir, "/background"))
  } else {
    nameso <- sort(colnames(occ)[-(1:3)])
    namesm <- sort(colnames(M_variables)[-(1:3)])

    if (!all.equal(nameso, namesm)) {
      stop("Variables in 'occ' and 'M_variables' must be the same.")
    }

    write.csv(occ, occs, row.names = FALSE)
    dir.create(paste0(out.dir, "/background"))
    write.csv(M_variables, dirM, row.names = FALSE)
  }

  # Slash
  if(.Platform$OS.type == "unix") {sl <- "/"; dl <- "/"} else {sl <- "\\"; dl <- "\\\\"}

  # Data
  ## Environmental variables sets
  ms <- gsub("/", dl, dirM)
  env <- paste("environmentallayers=", paste("\"", ms, "\"", sep = ""), sep = "")

  ## Species occurrences
  samp <- paste("samplesfile=", gsub("/", dl, paste0("\"", occs, "\"")), sep = "")

  # Maxent settings
  ## Feature classes combinations
  if (!is.null(f.clas)) {
    fea <- feature_classes(f.clas)
    a.fea <- "autofeature=false"
  } else {
    a.fea <- "autofeature=true"
  }

  # output directory
  odir <- paste0(out.dir, "/model")
  dir.create(odir)
  subfol <- paste0("outputdirectory=", gsub("/", dl, paste0("\"", odir, "\"")))

  # Getting ram to be used
  ram <- paste("-mx", max.memory, "m", sep = "")

  # Fixed commands
  ## Intitial command
  in.comm <- paste("java", ram, paste("-jar", gsub("/", dl, paste0("\"", paste(maxent.path, "maxent.jar", sep = sl), "\""))))

  ## Other maxent settings
  fin.com <- "extrapolate=false doclamp=false replicates=1 replicatetype=Crossvalidate responsecurves=false jackknife=true plots=false pictures=false outputformat=raw warnings=false visible=false autorun\n"

  # Final set of calibration models
  ## preparing final arguments
  reg.m <- paste0("betamultiplier=", reg.mult)

  ## writing java code
  jmx <- ifelse(!is.null(f.clas),
                paste(in.comm, env, samp, subfol, reg.m, a.fea, fea, args, fin.com),
                paste(in.comm, env, samp, subfol, reg.m, a.fea, args, fin.com))

  batch <- paste0(out.dir, "/mx_jackniffe")

  if(.Platform$OS.type == "unix") {
    cat(c("#! /bin/csh\n", jmx), file = paste0(batch, ".sh"))
  } else {
    cat(jmx, file = paste0(batch, ".bat"))
  }

  # running models
  if (verbose == TRUE) {
    message("If asked, RUN as administrator")
  }
  run_maxent(batch, maxent.path, add_path = FALSE, wait = TRUE)


  # preparing results
  allres <- read.csv(paste0(odir, "/maxentResults.csv"))
  cols <- colnames(allres)

  ## relevant columns
  colcont <- grep("contribution", cols)
  colperm <- grep("permutation.importance", cols)
  colgain <- grep("Regularized.training.gain", cols)
  colwith <- grep("gain.with.only", cols)
  colwout <- grep("gain.without", cols)

  ## relevant values
  varnamesmx <- gsub(".contribution", "", cols[colcont])
  contrib <- data.frame(Variable = varnamesmx,
                        Contribution = unlist(allres[, colcont]))
  rownames(contrib) <- NULL

  permimp <- data.frame(Variable = varnamesmx,
                        Permutation_importance = unlist(allres[, colperm]))
  rownames(permimp) <- NULL

  jackkni <- data.frame(Variable = varnamesmx,
                        Training_gain_with = unlist(allres[, colwith]),
                        Training_gain_without = unlist(allres[, colwout]))
  rownames(jackkni) <- NULL

  rtg <- allres[, colgain]

  results <- list(Contribution = contrib, Permutation_importance = permimp,
                  Jackknife_results = list(Regularized_training_gain_model = c(rtg),
                                           Training_gain_with_without = jackkni))

  ## plotting
  if (plot == TRUE) {
    plot_contribution(results)
  }

  # results
  return(results)
}


#' Helper to plot variable contribution to single models
#' @param contribution_list a list of results obtained with
#' \code{\link{explore_var_contrib}}, or an element of a list obtained with
#' \code{\link{model_var_contrib}}.
#' @param col.cont color of contribution bars; default = "gray25".
#' @param col.imp color of importance bars; default = "gray25".
#' @param col.with color of regularized training gain bars when using individual
#' variables; default = "gray15".
#' @param col.without color of regularized training gain bars when using all
#' variables except the one in question; default = "gray65".
#' @param col.all color of line representing regularized training gain of the
#' model using all variables; default = "gray1".
#' @export
#' @importFrom graphics par barplot abline title box plot.new legend text layout
#' @rdname explore_var_contrib
#' @usage
#' plot_contribution(contribution_list, col.cont = "gray25",
#'                   col.imp = "gray25", col.with = "gray15",
#'                   col.without = "gray65", col.all = "gray1")


plot_contribution <- function(contribution_list, col.cont = "gray25",
                              col.imp = "gray25", col.with = "gray15",
                              col.without = "gray65", col.all = "gray1") {
  if (missing(contribution_list)) {
    stop("Argument 'contribution_list' must be defined")
  }

  cont <- matrix(contribution_list$Contribution$Contribution, nrow = 1)
  perm <- matrix(contribution_list$Permutation_importance$Permutation_importance,
                 nrow = 1)
  tgain <- contribution_list$Jackknife_results$Regularized_training_gain_model
  vgain <- rbind(contribution_list$Jackknife_results$Training_gain_with_without[, 2],
                 contribution_list$Jackknife_results$Training_gain_with_without[, 3])

  # par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # label space and place settings
  vars <- contribution_list$Contribution$Variable
  maxnc <- max(sapply(vars, nchar), na.rm = T)
  fc <- maxnc * 0.3

  yl <- 1:ncol(cont) / ncol(cont)
  yl <- yl - (min(yl) / 2)

  # plot
  layout(matrix(1:5, 1), widths = c(fc, 10, 10, 10, 3.5))
  par(mar = c(3.5, 0, 3, 0))
  plot.new()
  text(0.5, yl, vars)

  par(mar = c(3.5, 1, 3, 0))
  barplot(cont, las = 1, col = col.cont, horiz = T, border = NA,
          main = "Contribution")
  box(bty = "l")

  barplot(perm, las = 1, col = col.imp, horiz = T, border = NA,
          main = "Permutation importance")
  box(bty = "l")

  barplot(vgain, las = 1, col = c(col.with, col.without), horiz = T,
          beside = TRUE, xlim = c(0, c(tgain)), border = NA,
          main = "Jackkniffe results")
  abline(v = tgain, col = col.all, lwd = 1.5, lty = 2)
  title(xlab = "Regularized training gain", line = 2.2)
  box(bty = "l")

  par(mar = c(3.5, 0.15, 0, 0))
  plot.new()
  legend("bottom", legend = c("all", "with", "without"), lty = c(2, NA, NA),
         lwd = c(1.5, NA, NA), pch = c(NA, 22, 22), bty = "n", pt.lwd = 2,
         pt.bg = c(NA, col.with, col.without),
         col = c(col.all, col.with, col.without), x.intersp = 0.5)
}
