#' Helper function to select feature classes
#' @param f.clas (character) feature classes can be selected from five different
#' combination sets or manually. Combination sets are: "all", "basic", "no.t.h",
#' "no.h", and "no.t". Default = "all". basic = "l", "lq", "lqp", "lqpt", "lqpth".
#' Combinations "no.t.h", "no.h", and "no.t", exclude t and/or h. See details for
#' all the available potential combinations of feature classes.
#' @details
#' Below all potential combinations of feature classes are shown. Manual selection
#' can be done by creating a vector of one or more of the combinations of this
#' list. l = linear, q = quadratic, p = product, t = threshold, and h = hinge.
#' "l", "q", "p", "t", "h", "lq", "lp", "lt", "lh", "qp", "qt", "qh", "pt", "ph",
#' "th", "lqp", "lqt", "lqh", "lpt", "lph", "lth", "qpt", "qph", "qth", "pth",
#' "lqpt", "lqph", "lqth", "lpth", "qpth", and "lqpth".
#' @return character containing java code for defining feature classes in Maxent
#' candidate models.
#' @export

feature_classes <- function(f.clas = "all") {
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
           "linear=true quadratic=false product=false threshold=true hinge=true",
           "linear=false quadratic=true product=true threshold=true hinge=false",
           "linear=false quadratic=true product=true threshold=false hinge=true",
           "linear=false quadratic=true product=false threshold=true hinge=true",
           "linear=false quadratic=false product=true threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=true hinge=false",
           "linear=true quadratic=true product=true threshold=false hinge=true",
           "linear=true quadratic=true product=false threshold=true hinge=true",
           "linear=true quadratic=false product=true threshold=true hinge=true",
           "linear=false quadratic=true product=true threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=true hinge=true")

  names(fea) <- c("l", "q", "p", "t", "h", "lq", "lp", "lt", "lh", "qp", "qt", "qh",
                  "pt", "ph", "th", "lqp", "lqt", "lqh", "lpt", "lph", "lth", "qpt",
                  "qph", "qth", "pth", "lqpt", "lqph", "lqth", "lpth", "qpth", "lqpth")

  if(any(f.clas %in% c("all", "basic", "no.t.h", "no.h", "no.t"))) {
    if(f.clas == "all"){fea <- fea}
    if(f.clas == "basic"){fea <- fea[c(1, 6, 16, 26, 31)]}
    if(f.clas == "no.t.h"){fea <- fea[c(1:3, 6:7, 10, 16)]}
    if(f.clas == "no.h"){fea <- fea[c(1:4, 6:8, 10:11, 13, 16:17, 19, 22, 26)]}
    if(f.clas == "no.t"){fea <- fea[c(1:3, 5:7, 9:10, 12, 14, 16, 18, 20, 23, 27)]}
  }else{
    if (any(f.clas %in% names(fea))) {
      fea <- fea[f.clas]
    } else {
      stop("Argument 'f.clas' is not valid.")
    }
  }
  return(fea)
}


#' Helper function to run maxent.jar from R
#' @param batch (character) name of the batch file (bash for Unix) with the code
#' to create all candidate models.
#' @param maxent.path (character) the path were the maxent.jar file is in your
#' computer.
#' @param add_path (logical) whether to add full path to \code{batch}.
#' @param wait (logical) whether R waits until the running is done or not.
#' @export
run_maxent <- function(batch, maxent.path, add_path = TRUE, wait = FALSE) {

  if(.Platform$OS.type == "unix") {
    if (add_path == TRUE) {
      batfile_path <- file.path(getwd(), paste0(batch, ".sh"))
    } else {
      batfile_path <- file.path(paste0(batch, ".sh"))
    }

    r_wd <- getwd()
    setwd(maxent.path)

    system(paste("bash", batfile_path), wait = wait)

  } else {
    if (add_path == TRUE) {
      batfile_path <- file.path(getwd(), paste0(batch, ".bat"))
    } else {
      batfile_path <- file.path(paste0(batch, ".bat"))
    }

    r_wd <- getwd() # real working directory
    setwd(maxent.path) # change temporally the working directory

    system2(batfile_path, wait = wait, invisible = FALSE)
  }
  setwd(r_wd)
}


#' Helper function to wait until a file writing is done
#' @param file (character) name of the file of interest.
#' @export
wait_written_done <- function(file) {
  while (!file.exists(file)) {
    if (file.exists(file)) {break()}
    Sys.sleep(1)
  }
  stime <- Sys.time()
  Sys.sleep(0.1)
  fi <- file.info(file)
  di <- as.numeric(fi$mtime - stime)

  while (di >= 0) {
    stime <- Sys.time()
    Sys.sleep(0.1)
    fi <- file.info(file)
    di <- as.numeric(fi$mtime - stime)

    if (di < 0) {break()}
    Sys.sleep(1)
  }
  return(di < 0)
}


#' Helper function to plot omission rate and AICc results
#' @param summary.calibration data.frame containing the summary of all metrics
#' calculated for all models during calibration.
#' @export
plot_proc_aicc <- function(summary.calibration) {
  plot(na.omit(summary.calibration[[3]])[, 4] ~
         log(na.omit(summary.calibration[[3]])[, 5]),
       xlab = "Natural logarithm of AICc", las = 1, col = "#a6cee3",
       ylab = paste0("Omission rates at ", summary.calibration[[4]], "% threshold value"))

  nonsig <- summary.calibration[[3]][!summary.calibration[[3]][, 1] %in%
                                       summary.calibration[[5]][, 1], ]

  points(na.omit(nonsig)[, 4] ~ log(na.omit(nonsig)[, 5]),
         col = "#b2df8a", pch = 19, cex = 1.1)

  points(na.omit(summary.calibration[[2]])[, 4] ~
           log(na.omit(summary.calibration[[2]])[, 5]), col = "#1f78b4",
         pch = 17, cex = 1.4)

  legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
         pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), box.col = "white",
         col = c("#1f78b4", "#b2df8a", "#a6cee3"), bg = "white") #, inset = c(0.01, 0)
  box()
}


#' Helper function to calculate the AICc values (number of parameters).
#'
#' @param x An object derived from reading the lambdas file created for Maxent.
#' Use \code{\link[base]{readLines}} function to read the file.
#'
#' @export

n_par <- function(x) {
  lambdas <- x[1:(length(x) - 4)]
  countNonZeroParams <- function(x) {
    if (strsplit(x, split = ", ")[[1]][2] != "0.0")
      1
  }
  no.params <- sum(unlist(sapply(lambdas, countNonZeroParams)))
  return(no.params)
}



#' Helper function to select extrapolation options
#' @param ext.type (character) extrapolation type to be used. Options are:
#' "all", "ext_clam", "ext", and "no_ext". Default = "all".
#' @export
ext_type <- function(ext.type = "all") {
  if(ext.type == "ext_clam") {
    mid.com <- "extrapolate=true doclamp=true"
    ext.nam <- "_EC"
  }
  if(ext.type == "ext") {
    mid.com <- "extrapolate=true doclamp=false"
    ext.nam <- "_E"
  }
  if(ext.type == "no_ext") {
    mid.com <- "extrapolate=false doclamp=false"
    ext.nam <- "_NE"
  }
  if(ext.type == "all") {
    mid.com <- c("extrapolate=true doclamp=true",
                 "extrapolate=true doclamp=false",
                 "extrapolate=false doclamp=false")
    ext.nam <- c("_EC", "_E", "_NE")
  }
  return(list(code = mid.com, name = ext.nam))
}


# finds raster format type according to format name
rformat_type <- function(format) {
  if (missing(format)) {stop("Argument 'format' needs to be defined")}
  if (format == "GTiff") {format1 <- ".tif"}
  if (format == "EHdr") {format1 <- ".bil"}
  if (format == "ascii") {format1 <- ".asc"}
  return(format1)
}
