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
    if(f.clas == "basic"){fea <- fea[c(1, 6, 16, 25, 29)]}
    if(f.clas == "no.t.h"){fea <- fea[c(1:3, 6:7, 10, 16)]}
    if(f.clas == "no.h"){fea <- fea[c(1:4, 6:8, 10:11, 13, 16:17, 19, 21, 25)]}
    if(f.clas == "no.t"){fea <- fea[c(1:3, 5:7, 9:10, 12, 14, 16, 18, 20, 22, 26)]}
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
#' @param wait (logical) whether R waits until the runing is done or not.
#' @export
run_maxent <- function(batch, maxent.path, wait = FALSE) {
  if(.Platform$OS.type == "unix") {
    batfile_path <- file.path(getwd(), paste0(batch, ".sh"))
    r_wd <- getwd()
    setwd(maxent.path)

    system(paste("bash", batfile_path), wait = wait)

  } else {
    batfile_path <- file.path(getwd(), paste0(batch, ".bat"))
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


#' Helper function to create an HTML file summarizing results
#' @param path directory where the HTML file will be written; default = current.
#' @param file.name (character) name of the HTML file.
#' @return HTML file in \code{path}
#' @export
html_calibration <- function (path = getwd(), file.name) {
  rmdfile <- paste(path, paste0(file.name, ".Rmd"), sep = "/")
  cat("---\ntitle: \"ku_enm: calibration results\"\noutput:\n  html_document:\n      toc: true\n      toc_depth: 4\n---\n\n```{r setup, include=FALSE}\nknitr::opts_chunk$set(echo = TRUE)\n```\n\n<br>\n\n### Brief description of the model calibration and selection process\n\n```{r, echo=FALSE}\nst4 <- read.csv(\"calibration_results.csv\")\nsett <- as.character(st4[,1])\nsetts <- strsplit(sett, split = \"_\")\nrm <- vector()\nfor (i in 1:length(setts)) {\nrm[i] <- setts[[i]][2]\n}\nf.clas <- vector()\nfor (i in 1:length(setts)) {\nf.clas[i] <- setts[[i]][4]\n}\nvar.di <- vector()\nfor (i in 1:length(setts)) {\nvar.di[i] <- paste(setts[[i]][5:length(setts[[i]])], collapse = \"_\")\n}\nrm1 <- paste(unique(rm), collapse = \", \")\nf.clas1 <- paste(unique(f.clas), collapse = \", \")\nvar.di1 <- paste(unique(var.di), collapse = \", \")\npar <- rbind(rm1, f.clas1, var.di1)\n```\n\nThis is the final report of the ku_enm_ceval function implemented in the ku_enm R package.\n\nIn all, `r length(st4[,1])` candidate models, with parameters reflecting all combinations of `r length(unique(rm))` regularization multiplier settings, `r length(unique(f.clas))` feature class combinations, and `r length(unique(var.di))` distinct sets of environmental variables, have been evaluated. Model peformance was evaluated based on statistical significance (Partial_ROC), omission rates (OR), and the Akaike information criterion corrected for small sample sizes (AICc).\n\n```{r par, echo=FALSE}\ncolnames(par) <- \"Parameters\"\nrow.names(par) <- c(\"Regularization multipliers\", \"Feature classes\", \"Sets of predictors\")\nknitr::kable(par, digits=c(0,0), row.names = TRUE, caption = \"Table 1. Parameters of the candidate models.\")\n```\n\n<br>\n\nThe results presented below can be found in the calibration output folder if desired for further analyses.\n\n<br>\n<br>\n\n### Model calibration statistics\n\nIn the following table is information about how many models met the four selection criteria that this function uses.\n\n```{r, echo=FALSE}\nst <- read.csv(\"calibration_stats.csv\")\ncolnames(st) <- c(\"Criteria\",\t\"Number_of_models\")\nknitr::kable(st, digits=c(0,0), caption = \"Table 2. General statistics of models that met distinct criteria.\")\n```\n\n<br>\n<br>\n\n### Models selected according to user-defined criteria\n\nThe following table contains the models selected according to the user's pre-defined criteria.\n\nNote that if the selection argument was \"OR_AICc\", delta AICc values were recalculated only among models meeting the omission rate criterion (*E*).\n\n```{r, echo=FALSE}\nst1 <- read.csv(\"selected_models.csv\")\ncolnames(st1) <- c(\"Model\",\t\"Mean_AUC_ratio\",\t\"Partial_ROC\", gsub(\"[.]\", \"%\", colnames(st1)[4]), \"AICc\",\t\"delta_AICc\",\t\"W_AICc\",\t\"num_parameters\")\nknitr::kable(st1, digits=c(0,3,3,3,3,3,3,0), caption = \"Table 3. Performance statistics for the models selected based on the user's pre-defined critera.\")\n```\n\n<br>\n<br>\n\n### Model performance plot\n\nThe figure below shows the position of the selected models in the distribution of all candidate models in terms of statistical significance, omission rates, and AICc values.\n\n![Figure 1. Distribution of all models, non-statistically significant models, and selected models in terms of omission rates and AICc values.](calibration_figure.png){width=60%}\n\n<br>\n<br>\n\n### Performance statistics for all models\n\nFollowing are the performance statistics for all candidate models (a sample if more than 500 models). See file calibration_results.csv for the complete list.\n\n```{r, echo=FALSE}\nst4 <- read.csv(\"calibration_results.csv\")\nif (dim(st4)[1] > 500) {\n   st4 <- st4[1:500, ]\n}\ncolnames(st4) <-  c(\"Model\",\t\"Mean_AUC_ratio\",\t\"Partial_ROC\", gsub(\"[.]\", \"%\", colnames(st4)[4]), \"AICc\",\t\"delta_AICc\",\t\"W_AICc\",\t\"num_parameters\")\nknitr::kable(st4, digits=c(0,3,3,3,3,3,3,0), caption = \"Table 4. Performance statistics for all candidate models.\")\n```",
      file = rmdfile)
  rmarkdown::render(rmdfile, "html_document", quiet = TRUE)
  unlink(rmdfile)
}


#' Helper function to select extrapolation options
#' @param ext.type (character) extrapolation type to be used. Options are:
#' "all", "ext_clam", "ext", and "no_ext", default = "all".
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
