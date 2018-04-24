#' Creation of an html file with results from model calibration
#'
#' @description html.eval creates an html file that summarizes all outputs from
#' the model calibration and selection processes.
#'
#' @param path directory in wich the html file will be written, default = current directory.
#' @param file.name (character) name of the html file.
#'
#' @return An html file summarizing all results from the model calibration and selection process.
#'
#' @details This function is used along with the \code{\link{ku.enm.ceval}} function.
#'
#' @examples
#' path <- getwd()
#' name <- "evaluation_results"
#'
#' html_file <- html.eval(path = path, file.name = name)

html.eval <- function(path = getwd(), file.name) {

  sink(paste(path, paste(file.name, ".Rmd", sep = ""), sep = "/"))
  cat(
"---
title: \"ku.enm: calibration results\"
output:
  html_document:
      toc: true
      toc_depth: 4
---

\```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
\```

<br>

### Brief description of the model calibration and selection process

\```{r, echo=FALSE}
st4 <- read.csv(\"calibration_results.csv\")
sett <- as.character(st4[,1])
setts <- strsplit(sett, split = \"_\")
rm <- vector()
for (i in 1:length(setts)) {
rm[i] <- setts[[i]][2]
}
f.clas <- vector()
for (i in 1:length(setts)) {
f.clas[i] <- setts[[i]][4]
}
var.di <- vector()
for (i in 1:length(setts)) {
var.di[i] <- paste(setts[[i]][5:length(setts[[i]])], collapse = \"_\")
}
rm1 <- paste(unique(rm), collapse = \", \")
f.clas1 <- paste(unique(f.clas), collapse = \", \")
var.di1 <- paste(unique(var.di), collapse = \", \")
par <- rbind(rm1, f.clas1, var.di1)
\```

This is the final report of the ku.enm.ceval function implemented in the ku.enm R package.

A total of \`r length(st4[,1])\` candidate models with parameters reflecting all combinations of \`r length(unique(rm))\` regularization multiplier settings, \`r length(unique(f.clas))\` feature classes combinations, and \`r length(unique(var.di))\` distinct sets of environmental variables, have been evaluated. Model peformance was evaluated based on statistical significance (Partial_ROC), omission rates, and the Akaike information criterion corrected for small sample sizes (AICc).

\```{r par, echo=FALSE}
colnames(par) <- \"Parameters\"
row.names(par) <- c(\"Regularization multipliers\", \"Feature classes\", \"Sets of predictors\")
knitr::kable(par, digits=c(0,0), row.names = TRUE, caption = \"Table 1. Parameters of the candidate models.\")
\```

<br>
All the results presented below can be found in the calibration output folder for further analyses.

<br>
<br>

### Model calibration statistics

In the following table is information about how many models met the four criteria of selection that this function uses.

\```{r, echo=FALSE}
st <- read.csv(\"calibration_stats.csv\")
colnames(st) <- c(\"Criteria\",	\"Number_of_models\")
knitr::kable(st, digits=c(0,0), caption = \"Table 2. General statistics of models that met distinct criteria.\")
\```

<br>
<br>

### Best models according to user-defined criteria

The following table contains the best models selected by your pre-defined criteria.

Notice that if the selection criterion was \"OR_AICc\", models below the omission rate and among them those with lower AICc values, the delta AICc values were recalculated only among models below the omission rate (*E*).

\```{r, echo=FALSE}
best <- list.files(pattern = \"best\")
st1 <- read.csv(best)
colnames(st1) <- c(\"Model\",	\"Mean_AUC_ratio\",	\"Partial_ROC\", gsub(\"[.]\", \"%\", colnames(st1)[4]), \"AICc\",	\"delta_AICc\",	\"W_AICc\",	\"num_parameters\")
knitr::kable(st1, digits=c(0,3,3,3,3,3,3,0), caption = \"Table 3. Performance statistics of the best models selected based on the pre-defined critera.\")
\```

<br>
<br>

### Model performance plot

The figure below shows the position of your selected models in the distribution of all candidate models in terms of statistical significance, omission rates, and AICc values.

![Figure 1. Distribution of all, non-statistically significant, and selected models according to their omission rates and AICc values.](calibration_figure.png){width=60%}

<br>
<br>

### Performance statistics for all models

Following you will find the performance statistics for all candidate models.

\```{r, echo=FALSE}
st4 <- read.csv(\"calibration_results.csv\")
colnames(st4) <-  c(\"Model\",	\"Mean_AUC_ratio\",	\"Partial_ROC\", gsub(\"[.]\", \"%\", colnames(st4)[4]), \"AICc\",	\"delta_AICc\",	\"W_AICc\",	\"num_parameters\")
knitr::kable(st4, digits=c(0,3,3,3,3,3,3,0), caption = \"Table 4. Performance statistics for all of the candidate models.\")
\```"
      )
  sink()
  rmarkdown::render(paste(path, paste(file.name, ".Rmd", sep = ""), sep = "/"), "html_document", quiet = TRUE)
  unlink(paste(path, paste(file.name, ".Rmd", sep = ""), sep = "/"))

}
