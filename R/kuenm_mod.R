#' Creation of Maxent models with selected parameters
#'
#' @description kuenm_mod creates and executes a batch file for generating Maxent models using
#' parameters previously selected with the \code{\link{kuenm_ceval}} function.
#'
#' @param occ.joint (character) the  csv file with all the occurrences; columns must be: species, longitude, latitude.
#' @param M.var.dir (character) name of the forlder containing folders in which calibration environmental
#' datasets are placed.
#' @param out.eval (character) name of the folder where evaluation results were written.
#' @param batch (character) name of the batch file with the code to create final Maxent models.
#' @param rep.n (numeric) number of model replicates.
#' @param rep.type (character) the replicate type; it can be: "Crossvalidate", "Bootstrap", or "Subsample".
#' @param jackknife (logical) if TRUE, a jackknife process is performed while runing Maxent models, default = FALSE.
#' @param out.dir (character) name of the output directory to be created and in which all model subdirectories
#' will be created.
#' @param out.format (character) the model output format; it can be: "raw", "logistic", "cloglog", or "cumulative".
#' @param project (logical) if TRUE, your models will be projected to the scenarios in G.var.dir, default = FALSE.
#' @param G.var.dir (character) if project is TRUE, name of the folder containing folders in wHich variables of
#' projection scenarios are placed.
#' @param ext.type (character) if project is TRUE, is the extrapolation type of projections; can be: "all", "ext_clam",
#' "ext", and "no_ext", default = "all". ext = free extrapolation, ext_clam = extrapolation and clamping,
#' no_ext = no extrapolation, and all = all three of the options listed above.
#' @param write.mess (logical) if TRUE, grids of MESS analysis results will be written, default = FALSE.
#' @param write.clamp (logical) if TRUE, a grid of the spatial distribution of clamping will be written, default = FALSE.
#' @param run (logical) if TRUE, the batch runs after its creation; if FALSE, it will only be created and its running
#' would be manual, default = TRUE.
#' @param args (character) additional arguments that can be passed to Maxent. See the Maxent help for more information
#' on how to write these arguments, default = NULL.
#'
#' @return A folder named as out.dir with all the subfolders to save Maxent final model results when running the .bat file.
#' A batch file for creating all the final Maxent models with their projections if project = TRUE.
#'
#' @details Same requirements as in \code{\link{kuenm_cal}}

kuenm_mod <- function(occ.joint, M.var.dir, out.eval, batch, rep.n = 10, rep.type = "Bootstrap",
                       jackknife = FALSE, out.dir, out.format = "logistic", project = FALSE, G.var.dir,
                       ext.type = "all", write.mess = FALSE, write.clamp = FALSE, run = TRUE, args = NULL) {

  #####
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
  var.dir <- paste(M.var.dir, var.di, sep = "\\")

  #output directory
  dir.create(out.dir)

  #Defining maximum ram to be used (50% of free memory)
  ram <- paste("-mx", (round((get_free_ram()/1000)*0.5)), "m", sep = "")

  #####
  #Maxent settings
  ##Environmental calibration variables sets
  env <- vector()
  for (i in 1:length(var.dir)) {
    env[i] <- paste("environmentallayers=", var.dir[i], sep = "")
  }

  ##Species occurrences
  samp <- paste("samplesfile=", occ.joint, sep = "")

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
    G.dir <- paste(G.var.dir, var.di, sep = "/")
    G.dirs <- vector()
    for (i in 1:length(G.dir)) {
      dirs <- dir(G.dir[i])
      dires <- vector()
      for (j in 1:length(dirs)) {
        dires[j] <- paste(G.var.dir, var.di[i], dirs[j], sep = "\\")
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
  in.comm <- paste("java", ram, "-jar maxent.jar", sep = " ")

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
          pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
          sink(paste(batch, ".bat", sep = ""))

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
              dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.clamp, args, fin.com, sep = " "))
            }
          }
        }
        if(write.mess == FALSE & write.clamp == TRUE) {
          pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
          sink(paste(batch, ".bat", sep = ""))

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
              dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.mess, args, fin.com, sep = " "))
            }
          }
        }
        if(write.clamp == FALSE & write.mess == FALSE) {
          pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
          sink(paste(batch, ".bat", sep = ""))

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
              dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.clamp, w.mess, args, fin.com, sep = " "))
            }
          }
        }
      }else {
        pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
        sink(paste(batch, ".bat", sep = ""))

        for (i in 1:length(sett1)) {
          Sys.sleep(0.1)
          setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

          for (j in 1:length(ext.nam)) {
            subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
            dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))

            reg.m <- paste("betamultiplier=", rm[i], sep = "")
            cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], args, fin.com, sep = " "))
          }
        }
      }

      sink()
      suppressMessages(close(pb))

      cat("\nIf asked, allow runing as administrator.")
      shell.exec(file.path(getwd(), paste(batch, ".bat", sep = "")))

      cat("\nProcess finished\n")
      cat(paste("A maxent batch file for creating", length(sett1) * length(ext.nam), "final models and their projections has been written", sep = " "))
      cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
    }else {
      pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
      sink(paste(batch, ".bat", sep = ""))

      for (i in 1:length(sett1)) {
        Sys.sleep(0.1)
        setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

        subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], sep = "")
        dir.create(paste(out.dir, "/", sett1[i], sep = ""))

        reg.m <- paste("betamultiplier=", rm[i], sep = "")
        cat(paste(in.comm, env[i], samp, subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com, args, fin.com, sep = " "))
      }

      sink()
      suppressMessages(close(pb))

      cat("\nIf asked, allow runing as administrator.")
      shell.exec(file.path(getwd(), paste(batch, ".bat", sep = "")))

      cat("\nProcess finished\n")
      cat(paste("A maxent batch file for creating", length(sett1), "final models without projections has been written", sep = " "))
      cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
    }
  }else {
    if(project == TRUE) {
      if(write.clamp == FALSE | write.mess == FALSE) {
        if(write.clamp == FALSE & write.mess == TRUE) {
          pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
          sink(paste(batch, ".bat", sep = ""))

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
              dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.clamp, fin.com, sep = " "))
            }
          }
        }
        if(write.mess == FALSE & write.clamp == TRUE) {
          pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
          sink(paste(batch, ".bat", sep = ""))

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
              dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.mess, fin.com, sep = " "))
            }
          }
        }
        if(write.clamp == FALSE & write.mess == FALSE) {
          pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
          sink(paste(batch, ".bat", sep = ""))

          for (i in 1:length(sett1)) {
            Sys.sleep(0.1)
            setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

            for (j in 1:length(ext.nam)) {
              subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
              dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))

              reg.m <- paste("betamultiplier=", rm[i], sep = "")
              cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], w.clamp, w.mess, fin.com, sep = " "))
            }
          }
        }
      }else {
        pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
        sink(paste(batch, ".bat", sep = ""))

        for (i in 1:length(sett1)) {
          Sys.sleep(0.1)
          setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

          for (j in 1:length(ext.nam)) {
            subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
            dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))

            reg.m <- paste("betamultiplier=", rm[i], sep = "")
            cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com[j], fin.com, sep = " "))
          }
        }
      }

      sink()
      suppressMessages(close(pb))

      cat("\nIf asked, allow runing as administrator.")
      shell.exec(file.path(getwd(), paste(batch, ".bat", sep = "")))

      cat("\nProcess finished\n")
      cat(paste("A maxent batch file for creating", length(sett1) * length(ext.nam), "final models and their projections has been written", sep = " "))
      cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
    }else {
      pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
      sink(paste(batch, ".bat", sep = ""))

      for (i in 1:length(sett1)) {
        Sys.sleep(0.1)
        setWinProgressBar(pb, i, title = paste(round(i / length(sett1) * 100, 0), "% finished"))

        subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], sep = "")
        dir.create(paste(out.dir, "/", sett1[i], sep = ""))

        reg.m <- paste("betamultiplier=", rm[i], sep = "")
        cat(paste(in.comm, env[i], samp, subfol, reg.m, a.fea, fea[i], rep, rept, jack, out, mid.com, fin.com, sep = " "))
      }

      sink()
      suppressMessages(close(pb))

      cat("\nIf asked, allow runing as administrator.")
      shell.exec(file.path(getwd(), paste(batch, ".bat", sep = "")))

      cat("\nProcess finished\n")
      cat(paste("A maxent batch file for creating", length(sett1), "final models without projections has been written", sep = " "))
      cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
    }
  }
}
