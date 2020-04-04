#' Creation of Maxent models with selected parameters
#'
#' @description kuenm_mod creates and executes a batch file (bash for Unix) for generating Maxent models using
#' parameters previously selected with the \code{\link{kuenm_ceval}} function.
#'
#' @param occ.joint (character) the  csv file with all the occurrences; columns must be: species, longitude, latitude.
#' @param M.var.dir (character) name of the forlder containing folders in which calibration environmental
#' datasets are placed.
#' @param out.eval (character) name of the folder where evaluation results (from calibration) were written.
#' @param batch (character) name of the batch file (bash for Unix) with the code to create final Maxent models.
#' @param rep.n (numeric) number of model replicates, default = 10.
#' @param rep.type (character) the replicate type; it can be: "Crossvalidate", "Bootstrap", or "Subsample".
#' @param jackknife (logical) if TRUE, a jackknife process is performed while runing Maxent models, default = FALSE.
#' @param out.dir (character) name of the output directory to be created and in which all model subdirectories
#' will be created.
#' @param max.memory (numeric) maximum memory (in megabytes) to be used by maxent while creating the models. Default = 1000.
#' @param out.format (character) the model output format; it can be: "raw", "logistic", "cloglog", or "cumulative".
#' @param project (logical) if TRUE, your models will be projected to scenarios in G.var.dir, default = FALSE.
#' @param G.var.dir (character) if project is TRUE, name of the folder containing folders in which variables of
#' projection scenarios are placed.
#' @param ext.type (character) if project is TRUE, is the extrapolation type of projections; can be: "all", "ext_clam",
#' "ext", and "no_ext", default = "all". ext = free extrapolation, ext_clam = extrapolation and clamping,
#' no_ext = no extrapolation, and all = all three of the options listed above.
#' @param write.mess (logical) if TRUE, grids of MESS analysis results will be written, default = FALSE.
#' @param write.clamp (logical) if TRUE, a grid of the spatial distribution of clamping will be written, default = FALSE.
#' @param maxent.path (character) the path were the maxent.jar file is in your computer.
#' @param args (character) additional arguments that can be passed to Maxent. See the Maxent help for more information
#' on how to write these arguments, default = NULL. Note that some arguments cannot be changed here because they are
#' part of the parameters of the function already (e.g., "writemess"). See details for other options.
#' @param wait (logical) if TRUE R will wait until all the Maxent models are created. If FALSE the process of
#' model creation will be performed separately and R could be used at the same time. This may be useful for evaluating
#' candidate models parallelly. Default = TRUE.
#' @param run (logical) if TRUE, the batch runs after its creation; if FALSE, it will only be created and its running
#' would be manual, default = TRUE.
#'
#' @return A folder named as out.dir with all the subfolders to save Maxent final model results when running the .bat file
#' (.sh for Unix). A batch file (bash for Unix) for creating all the final Maxent models with their projections if project = TRUE.
#'
#' @details Same requirements regarding Java and maxent than in \code{\link{kuenm_cal}}.
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
#' kuenm_mod(occ.joint, M.var.dir, out.eval, batch, rep.n = 10, rep.type = "Bootstrap",
#'           jackknife = FALSE, out.dir, max.memory = 1000, out.format = "logistic",
#'           project = FALSE, G.var.dir, ext.type = "all", write.mess = FALSE,
#'           write.clamp = FALSE, maxent.path, args = NULL, wait = TRUE, run = TRUE)
#'
#' @export
#'
#' @examples
#' # To run this function model evaluation and selection using the kuenm_ceval function should have been used before.
#' # The evaluation function generates one of the imputs needed.
#'
#' # Variables with information to be used as arguments.
#' occ_joint <- "aame_joint.csv"
#' M_var_dir <- "M_variables"
#' out_eval <- "Calibration_results"
#' batch_fin <- "Final_models"
#' mod_dir <- "Final_Models"
#' rep_n <- 10
#' rep_type <- "Bootstrap"
#' jackknife <- FALSE
#' G_var_dir <- "G_variables"
#' out_format <- "logistic"
#' project <- TRUE
#' ext_type <- "all"
#' write_mess <- FALSE
#' write_clamp <- FALSE
#' wait1 <- FALSE
#' run1 <- TRUE
#' args <- NULL
#'
#' kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
#'           rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
#'           G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp,
#'           maxent.path = maxent_path, args = args, wait = wait1, run = run1)

kuenm_mod <- function(occ.joint, M.var.dir, out.eval, batch, rep.n = 10, rep.type = "Bootstrap",
                      jackknife = FALSE, out.dir, max.memory = 1000, out.format = "logistic",
                      project = FALSE, G.var.dir, ext.type = "all", write.mess = FALSE,
                      write.clamp = FALSE, maxent.path, args = NULL, wait = TRUE, run = TRUE) {

  #Checking potential issues
  if (!file.exists(occ.joint)) {
    stop(paste(occ.joint, "does not exist in the working directory, check file name",
               "\nor extension, example: species_joint.csv"))
  }
  if (missing(M.var.dir)) {
    stop("Argument M.var.dir is not defined.")
  }
  if (!dir.exists(M.var.dir)) {
    stop(paste(M.var.dir, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (length(list.dirs(M.var.dir, recursive = FALSE)) == 0) {
    stop(paste(M.var.dir, "does not contain any subdirectory with environmental variables,",
               "\neach set of variables must be in a subdirectory inside",
               paste(M.var.dir, ".", sep = "")))
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
    stop(paste("Argument maxent.path is not defined, it is necessary for executing",
               "\nthe Maxent software."))
  }


  #####
  if(.Platform$OS.type == "unix") {
    sl <- "/"
    dl <- "/"
  } else {
    sl <- "\\"
    dl <- "\\\\"
  }

  #Data
  ##Data from best models table
  best <- list.files(path = out.eval, pattern = "best", full.names = TRUE)
  sett <- read.csv(best)
  sett1 <- as.character(sett[, 1])
  setts <- strsplit(sett1, split = "_")

  ###Regularization multipliers
  rm <- vector()
  for (i in 1:length(setts)) {
    rm[i] <- setts[[i]][2]
  }

  ###Feature classes
  f.clas <- vector()
  for (i in 1:length(setts)) {
    f.clas[i] <- setts[[i]][4]
  }

  ###Calibration (M) variables
  var.di <- vector()
  for (i in 1:length(setts)) {
    var.di[i] <- paste(setts[[i]][5:length(setts[[i]])], collapse = "_")
  }
  var.dir <- paste("\"", paste(gsub("/", dl, paste(getwd(), M.var.dir, sep = sl)), var.di, sep = sl),
                   "\"", sep = "")

  #output directory
  dir.create(out.dir)
  out.dir <- gsub("/", dl, paste(getwd(), out.dir, sep = sl))

  #Defining maximum ram to be used (50% of free memory)
  ram <- paste("-mx", max.memory, "m", sep = "")

  #####
  #Maxent settings
  ##Environmental calibration variables sets
  env <- vector()
  for (i in 1:length(var.dir)) {
    env[i] <- paste("environmentallayers=", var.dir[i], sep = "")
  }

  ##Species occurrences
  samp <- paste("samplesfile=", gsub("/", dl,
                                     paste("\"", paste(getwd(), occ.joint, sep = sl),
                                     "\"", sep = "")), sep = "")

  ##Feature classes combinations
  fea <- c("linear=true quadratic=false product=false threshold=false hinge=false",
           "linear=false quadratic=true product=false threshold=false hinge=false",
           "linear=false quadratic=false product=true threshold=false hinge=false",
           "linear=false quadratic=false product=false threshold=true hinge=false",
           "linear=false quadratic=false product=false threshold=false hinge=true",
           "linear=true quadratic=true product=false threshold=false hinge=false",
           "linear=true quadratic=false product=true threshold=false hinge=false",
           "linear=true quadratic=false product=false threshold=true hinge=false",
           "linear=true quadratic=false product=false threshold=false hinge=true",
           "linear=false quadratic=true product=true threshold=false hinge=false",
           "linear=false quadratic=true product=false threshold=true hinge=false",
           "linear=false quadratic=true product=false threshold=false hinge=true",
           "linear=false quadratic=false product=true threshold=true hinge=false",
           "linear=false quadratic=false product=true threshold=false hinge=true",
           "linear=false quadratic=false product=false threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=false hinge=false",
           "linear=true quadratic=true product=false threshold=true hinge=false",
           "linear=true quadratic=true product=false threshold=false hinge=true",
           "linear=true quadratic=false product=true threshold=true hinge=false",
           "linear=true quadratic=false product=true threshold=false hinge=true",
           "linear=false quadratic=true product=true threshold=true hinge=false",
           "linear=false quadratic=true product=true threshold=false hinge=true",
           "linear=false quadratic=true product=false threshold=true hinge=true",
           "linear=false quadratic=false product=true threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=true hinge=false",
           "linear=true quadratic=true product=true threshold=false hinge=true",
           "linear=true quadratic=true product=false threshold=true hinge=true",
           "linear=true quadratic=false product=true threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=true hinge=true")

  names(fea) <- c("l", "q", "p", "t", "h", "lq", "lp", "lt", "lh", "qp", "qt", "qh",
                  "pt", "ph", "th", "lqp", "lqt", "lqh", "lpt", "lph", "qpt", "qph",
                  "qth", "pth", "lqpt", "lqph", "lqth", "lpth", "lqpth")

  ###Selected feature classes using data from the best models table
  fea <- fea[f.clas]

  ##Projection (G) variables folders and subfolders, extrapolation types, and writting clamp and MESS
  if(project == TRUE) {
    G.dir <- paste(gsub("/", dl, paste(getwd(), G.var.dir, sep = sl)), var.di, sep = sl)
    G.dirs <- vector()
    for (i in 1:length(G.dir)) {
      dirs <- dir(G.dir[i])
      dires <- vector()
      for (j in 1:length(dirs)) {
        dires[j] <- paste("\"", paste(gsub("/", dl, paste(getwd(), G.var.dir, sep = sl)),
                          var.di[i], dirs[j], sep = sl), "\"", sep = "")
      }
      G.dirs[i] <- paste("projectionlayers=", paste(dires, collapse = ","), sep = "")
    }

    if(ext.type == "ext_clam") {
      mid.com <- "extrapolate=true doclamp=true responsecurves=true"
      ext.nam <- "_EC"
    }
    if(ext.type == "ext") {
      mid.com <- "extrapolate=true doclamp=false responsecurves=true"
      ext.nam <- "_E"
    }
    if(ext.type == "no_ext") {
      mid.com <- "extrapolate=false doclamp=false responsecurves=true"
      ext.nam <- "_NE"
    }
    if(ext.type == "all") {
      mid.com <- c("extrapolate=true doclamp=true responsecurves=true",
                   "extrapolate=true doclamp=false responsecurves=true",
                   "extrapolate=false doclamp=false responsecurves=true")
      ext.nam <- c("_EC", "_E", "_NE")
    }

    if(write.mess == FALSE) {
      w.mess <- "writeclampgrid=false"
    }
    if(write.clamp == FALSE) {
      w.clamp <- "writemess=false"
    }
  }else {
    mid.com <- "extrapolate=false doclamp=false writeclampgrid=false writemess=false responsecurves=true"
  }

  ##Jackknife
  if(jackknife == TRUE) {
    jack <- "jackknife=true"
  }else {
    jack <- "jackknife=false"
  }

  ##Output format
  out <- paste("outputformat=", out.format, sep = "")

  ##Number of replicates
  rep <-  paste("replicates=", rep.n, sep = "")

  ##Replicate type
  rept <- paste("replicatetype=", rep.type, sep = "")

  #Fixed commands
  ##Intitial command
  in.comm <- paste("java", ram,
                   paste("-jar",
                         gsub("/", dl,
                              paste("\"", paste(maxent.path, "maxent.jar", sep = sl), "\"", sep = ""))),
                   sep = " ")

  ##Autofeature
  a.fea <- "autofeature=false"

  ##Other maxent settings
  fin.com <- "warnings=false visible=false redoifexists autorun\n"

  #####
  #Final code
  if (exists("args") == TRUE) {
    if(project == TRUE) {
      if(write.clamp == FALSE | write.mess == FALSE) {
        if(write.clamp == FALSE & write.mess == TRUE) {
          if(.Platform$OS.type == "unix") {
            cat("\nCreating directories and maxent batch file, please wait...\n")
            sink(paste(batch, ".sh", sep = "")) # beginning file preparation
            cat("#! /bin/csh\n")
          } else {
            pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
            sink(paste(batch, ".bat", sep = "")) # beginning file preparation
          }

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            if(.Platform$OS.type == "unix") {

            } else {
              setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
            }

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i], ext.nam[j],
                                                        "\"", sep = ""), sep = "")
              dir.create(paste(out.dir, sl, sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.clamp, args, fin.com, sep = " "))
            }
          }
        }
        if(write.mess == FALSE & write.clamp == TRUE) {
          if(.Platform$OS.type == "unix") {
            cat("\nCreating directories and maxent batch file, please wait...\n")
            sink(paste(batch, ".sh", sep = "")) # beginning file preparation
            cat("#! /bin/csh\n")
          } else {
            pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
            sink(paste(batch, ".bat", sep = "")) # beginning file preparation
          }

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            if(.Platform$OS.type == "unix") {

            } else {
              setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
            }

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i], ext.nam[j],
                                                        "\"", sep = ""), sep = "")
              dir.create(paste(out.dir, sl, sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.mess, args, fin.com, sep = " "))
            }
          }
        }
        if(write.clamp == FALSE & write.mess == FALSE) {
          if(.Platform$OS.type == "unix") {
            cat("\nCreating directories and maxent batch file, please wait...\n")
            sink(paste(batch, ".sh", sep = "")) # beginning file preparation
            cat("#! /bin/csh\n")
          } else {
            pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
            sink(paste(batch, ".bat", sep = "")) # beginning file preparation
          }

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            if(.Platform$OS.type == "unix") {

            } else {
              setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
            }

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i], ext.nam[j],
                                                        "\"", sep = ""), sep = "")
              dir.create(paste(out.dir, sl, sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.clamp, w.mess, args, fin.com, sep = " "))
            }
          }
        }
      }else {
        if(.Platform$OS.type == "unix") {
          cat("\nCreating directories and maxent batch file, please wait...\n")
          sink(paste(batch, ".sh", sep = "")) # beginning file preparation
          cat("#! /bin/csh\n")
        } else {
          pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
          sink(paste(batch, ".bat", sep = "")) # beginning file preparation
        }

        for (i in 1:length(sett1)) {
          Sys.sleep(0.1)
          if(.Platform$OS.type == "unix") {

          } else {
            setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
          }

          for (j in 1:length(ext.nam)) {
            subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i], ext.nam[j],
                                                      "\"", sep = ""), sep = "")
            dir.create(paste(out.dir, sl, sett1[i], ext.nam[j], sep = ""))

            reg.m <- paste("betamultiplier=", rm[i], sep = "")
            cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], args, fin.com, sep = " "))
          }
        }
      }

      sink()
      if(.Platform$OS.type != "unix") {
        suppressMessages(close(pb))
      }

      cat("\nIf asked, allow runing as administrator.")
      if(run == TRUE){
        if(.Platform$OS.type == "unix") {
          batfile_path <- file.path(getwd(), paste(batch, ".sh", sep = "")) # bat file
          r_wd <- getwd() # real working directory
          setwd(maxent.path) # change temporally the working directory

          system(paste("bash", batfile_path), wait = wait)
        } else {
          batfile_path <- file.path(getwd(), paste(batch, ".bat", sep = "")) # bat file
          r_wd <- getwd() # real working directory
          setwd(maxent.path) # change temporally the working directory

          system2(batfile_path, wait = wait, invisible = FALSE)
        }
        setwd(r_wd) # return actual working directory
      }

      cat("\nProcess finished\n")
      if(.Platform$OS.type == "unix") {
        cat(paste("A maxent shell script for creating", length(sett1) * length(ext.nam), "final models and their projections has been written", sep = " "))
      } else {
        cat(paste("A maxent batch file for creating", length(sett1) * length(ext.nam), "final models and their projections has been written", sep = " "))
      }
      cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
    }else {
      if(.Platform$OS.type == "unix") {
        cat("\nCreating directories and maxent batch file, please wait...\n")
        sink(paste(batch, ".sh", sep = "")) # beginning file preparation
        cat("#! /bin/csh\n")
      } else {
        pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
        sink(paste(batch, ".bat", sep = "")) # beginning file preparation
      }

      for (i in 1:length(sett1)) {
        Sys.sleep(0.1)
        if(.Platform$OS.type == "unix") {

        } else {
          setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
        }

        subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i],
                                                  "\"", sep = ""), sep = "")
        dir.create(paste(out.dir, sl, sett1[i], sep = ""))

        reg.m <- paste("betamultiplier=", rm[i], sep = "")
        cat(paste(in.comm, env[i], samp, subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com, args, fin.com, sep = " "))
      }

      sink()
      if(.Platform$OS.type != "unix") {
        suppressMessages(close(pb))
      }

      cat("\nIf asked, allow runing as administrator.")
      if(run == TRUE){
        if(.Platform$OS.type == "unix") {
          batfile_path <- file.path(getwd(), paste(batch, ".sh", sep = "")) # bat file
          r_wd <- getwd() # real working directory
          setwd(maxent.path) # change temporally the working directory

          system(paste("bash", batfile_path), wait = wait)
        } else {
          batfile_path <- file.path(getwd(), paste(batch, ".bat", sep = "")) # bat file
          r_wd <- getwd() # real working directory
          setwd(maxent.path) # change temporally the working directory

          system2(batfile_path, wait = wait, invisible = FALSE)
        }
        setwd(r_wd) # return actual working directory
      }

      cat("\nProcess finished\n")
      if(.Platform$OS.type == "unix") {
        cat(paste("A maxent shell script for creating", length(sett1), "final models without projections has been written", sep = " "))
      } else {
        cat(paste("A maxent batch file for creating", length(sett1), "final models without projections has been written", sep = " "))
      }
      cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
    }
  }else {
    if(project == TRUE) {
      if(write.clamp == FALSE | write.mess == FALSE) {
        if(write.clamp == FALSE & write.mess == TRUE) {
          if(.Platform$OS.type == "unix") {
            cat("\nCreating directories and maxent batch file, please wait...\n")
            sink(paste(batch, ".sh", sep = "")) # beginning file preparation
            cat("#! /bin/csh\n")
          } else {
            pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
            sink(paste(batch, ".bat", sep = "")) # beginning file preparation
          }

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            if(.Platform$OS.type == "unix") {

            } else {
              setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
            }

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i], ext.nam[j],
                                                        "\"", sep = ""), sep = "")
              dir.create(paste(out.dir, sl, sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.clamp, fin.com, sep = " "))
            }
          }
        }
        if(write.mess == FALSE & write.clamp == TRUE) {
          if(.Platform$OS.type == "unix") {
            cat("\nCreating directories and maxent batch file, please wait...\n")
            sink(paste(batch, ".sh", sep = "")) # beginning file preparation
            cat("#! /bin/csh\n")
          } else {
            pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
            sink(paste(batch, ".bat", sep = "")) # beginning file preparation
          }

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            if(.Platform$OS.type == "unix") {

            } else {
              setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
            }

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i], ext.nam[j],
                                                        "\"", sep = ""), sep = "")
              dir.create(paste(out.dir, sl, sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.mess, fin.com, sep = " "))
            }
          }
        }
        if(write.clamp == FALSE & write.mess == FALSE) {
          if(.Platform$OS.type == "unix") {
            cat("\nCreating directories and maxent batch file, please wait...\n")
            sink(paste(batch, ".sh", sep = "")) # beginning file preparation
            cat("#! /bin/csh\n")
          } else {
            pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
            sink(paste(batch, ".bat", sep = "")) # beginning file preparation
          }

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            if(.Platform$OS.type == "unix") {

            } else {
              setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
            }

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i], ext.nam[j],
                                                        "\"", sep = ""), sep = "")
              dir.create(paste(out.dir, sl, sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.clamp, w.mess, fin.com, sep = " "))
            }
          }
        }
      }else {
        if(.Platform$OS.type == "unix") {
          cat("\nCreating directories and maxent batch file, please wait...\n")
          sink(paste(batch, ".sh", sep = "")) # beginning file preparation
          cat("#! /bin/csh\n")
        } else {
          pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
          sink(paste(batch, ".bat", sep = "")) # beginning file preparation
        }

        for (i in 1:length(sett1)) {
          Sys.sleep(0.1)
          if(.Platform$OS.type == "unix") {

          } else {
            setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
          }

          for (j in 1:length(ext.nam)) {
            subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i], ext.nam[j],
                                                      "\"", sep = ""), sep = "")
            dir.create(paste(out.dir, sl, sett1[i], ext.nam[j], sep = ""))

            reg.m <- paste("betamultiplier=", rm[i], sep = "")
            cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], fin.com, sep = " "))
          }
        }
      }

      sink()
      if(.Platform$OS.type != "unix") {
        suppressMessages(close(pb))
      }

      cat("\nIf asked, allow runing as administrator.")
      if(run == TRUE){
        if(.Platform$OS.type == "unix") {
          batfile_path <- file.path(getwd(), paste(batch, ".sh", sep = "")) # bat file
          r_wd <- getwd() # real working directory
          setwd(maxent.path) # change temporally the working directory

          system(paste("bash", batfile_path), wait = wait)
        } else {
          batfile_path <- file.path(getwd(), paste(batch, ".bat", sep = "")) # bat file
          r_wd <- getwd() # real working directory
          setwd(maxent.path) # change temporally the working directory

          system2(batfile_path, wait = wait, invisible = FALSE)
        }
        setwd(r_wd) # return actual working directory
      }

      cat("\nProcess finished\n")
      if(.Platform$OS.type == "unix") {
        cat(paste("A maxent shell script for creating", length(sett1) * length(ext.nam), "final models and their projections has been written", sep = " "))
      } else {
        cat(paste("A maxent batch file for creating", length(sett1) * length(ext.nam), "final models and their projections has been written", sep = " "))
      }
      cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
    }else {
      if(.Platform$OS.type == "unix") {
        cat("\nCreating directories and maxent batch file, please wait...\n")
        sink(paste(batch, ".sh", sep = "")) # beginning file preparation
        cat("#! /bin/csh\n")
      } else {
        pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
        sink(paste(batch, ".bat", sep = "")) # beginning file preparation
      }

      for (i in 1:length(sett1)) {
        Sys.sleep(0.1)
        if(.Platform$OS.type == "unix") {

        } else {
          setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))
        }

        subfol <- paste("outputdirectory=", paste("\"", out.dir, sl, sett1[i],
                                                  "\"", sep = ""), sep = "")
        dir.create(paste(out.dir, sl, sett1[i], sep = ""))

        reg.m <- paste("betamultiplier=", rm[i], sep = "")
        cat(paste(in.comm, env[i], samp, subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com, fin.com, sep = " "))
      }

      sink()
      if(.Platform$OS.type != "unix") {
        suppressMessages(close(pb))
      }

      cat("\nIf asked, allow runing as administrator.")
      if(run == TRUE){
        if(.Platform$OS.type == "unix") {
          batfile_path <- file.path(getwd(), paste(batch, ".sh", sep = "")) # bat file
          r_wd <- getwd() # real working directory
          setwd(maxent.path) # change temporally the working directory

          system(paste("bash", batfile_path), wait = wait)
        } else {
          batfile_path <- file.path(getwd(), paste(batch, ".bat", sep = "")) # bat file
          r_wd <- getwd() # real working directory
          setwd(maxent.path) # change temporally the working directory

          system2(batfile_path, wait = wait, invisible = FALSE)
        }
        setwd(r_wd) # return actual working directory
      }

      cat("\nProcess finished\n")
      if(.Platform$OS.type == "unix") {
        cat(paste("A maxent shell script for creating", length(sett1), "final models without projections has been written", sep = " "))
      } else {
        cat(paste("A maxent batch file for creating", length(sett1), "final models without projections has been written", sep = " "))
      }
      cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
    }
  }
}
