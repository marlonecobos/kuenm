#' Creation of an R Markdown file for recording all analyses
#'
#' @description ku.enm.start generates an R Markdown file that serves as a guide for
#' performing all the analyses that this package includes.
#'
#' @param file.name (character) is the name of the R Markdown file that will be
#' produced in your working directory.
#'
#' @return An R Markdown file with instructions and code for performing all
#' analyses included in this package.
#'
#' @examples
#' ku.enm.start("complete_process")
#' ku.enm.start("bird_modeling_process")

ku.enm.start <- function(file.name){

  sink(paste(file.name, ".Rmd", sep = ""))
  cat(
"---
title: \"ku.enm: modeling process\"
output:
  html_document:
      toc: yes
      toc_depth: 4
---

This R Markdown file is created in the working directory, and is designed to make more reproducible the processes of model calibration and final model creation.

Information on using this R Markdown file:

- Try executing code chunks by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.
- Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

A brief tutorial for using functions of the ku.enm R package can be found in https://github.com/manubio13/ku.enm/blob/master/Tutorial/ku.enm_vignette.md. Additionally, function help can be checked to change arguments according to specific needs.

### Candidate models

Candidate models are a large set of candidate models created to respond to the need to test multiple parameter combinations: for example, distinct regularization multiplier values, various feature classes, and different sets of environmental variables. The following code calls the help page of the function ku.enm.cal.

\```{r, eval=FALSE}
help(ku.enm.cal)
\```

The next chunk of code is for preparing the arguments for using the function following the modularity principle. These variables can be changed according to each case.

\```{r, eval=FALSE}
occ_joint <- \"Sp_joint.csv\"
occ_tra <- \"Sp_train.csv\"
M_var_dir <- \"M_variables\"
batch_cal <- \"Candidate_models\"
out_dir <- \"Candidate_Models\"
reg_mult <- c(seq(0.1,1,0.1),seq(2,6,1),8,10)
f_clas <- \"all\"
run <- TRUE
\```

The following is the code for using the function.

\```{r, eval=FALSE}
ku.enm.cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
           out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, run = run)
\```


### Evaluation and selection of best models

Evaluation is a crucial step in model calibration. This step centers on selecting candidate models and their associated parameters to identify the very best models for the purposes of the study. The ku.enm.eval function evaluates candidate models based on three distinct criteria: statistical significance (based on partial ROC), prediction ability (we use omission rates, but other metrics, such as overall correct classification rate, can also be used), and complexity (here evaluated using AICc). The following code chunk calls the function help window.

\```{r, eval=FALSE}
help(ku.enm.ceval)
\```

Below arguments for this functions will be defined.

\```{r, eval=FALSE}
occ_test <- \"Sp_test.csv\"
out_eval <- \"Calibration_results\"
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- \"OR_AICc\"
# Notice, some of the variables used here as arguments were already created for the previous function
\```

And this code allows evaluating candidate models that were created previously, selecting those with best performance based on the three criteria.

\```{r, eval=FALSE}
ku.enm.ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
            out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
            kept = kept, selection = selection)
\```

### Final model creation

After selecting parametrizations that produce the best models, the next step is to create the final models, and if needed transfer them to other environmental data sets (e.g., to other time periods or other geographic regions). The function help is called via this code:

\```{r, eval=FALSE}
help(ku.enm.mod)
\```

For preparing the arguments for this function use the following chunk of code.

\```{r, eval=FALSE}
batch_fin <- \"Final_models\"
mod_dir <- \"Final_Models\"
rep_n <- 10
rep_type <- \"Bootstrap\"
jackknife <- FALSE
G_var_dir <- \"G_variables\"
out_format <- \"logistic\"
project <- TRUE
ext_type <- \"all\"
write_mess <- FALSE
write_clamp <- FALSE
run1 <- TRUE
# Again, some of the variables used here as arguments were already created for the previous functions
\```

And the ku.enm.mod function has the following syntax:

\```{r, eval=FALSE}
ku.enm.mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
           rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
           G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, run = run1)
\```


### Final model evaluation

Final model can be evaluated after being created, for this step independent data is needed (data that has not being used in the calibration process and usually comming from different sources). The function help is called via this code:

\```{r, eval=FALSE}
help(ku.enm.feval)
\```

The next chunk of code is for preparing the arguments for using the function.

\```{r, eval=FALSE}
occ_ind <- \"Sp_ind.csv\"
replicates <- TRUE
out_feval <- \"Final_Models_evaluation\"
# Most of the variables used here as arguments were already created for the previous functions
\```

The following is the code for using the function.

\```{r, eval=FALSE}
ku.enm.feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
             out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
             iterations = iterations)
\```


### MOP analysis

If transfers were performed at the moment of creating final models, the MOP analysis will healp to identify areas of strict extrapolation and levels of similarity between the M and the area or scenario of projection. The code below will help to see the function's description:

\```{r, eval=FALSE}
help(ku.enm.mmop)
\```

Below arguments for this functions will be defined.

\```{r, eval=FALSE}
sets_var <- c(\"Set3\")
out_mop <- \"MOP_results\"
percent <- 10
normalized <- TRUE
# Some of the variables used here as arguments were already created for the previous functions
\```

And the ku.enm.mmop function has the following syntax:

\```{r, eval=FALSE}
ku.enm.mmop(dirG = G_var_dir, dirM = M_var_dir, sets.var = sets_var, out.mop = out_mop,
            percent = percent, normalized = normalized)
\```"
  )
  sink()
  file.edit(paste(file.name, ".Rmd", sep = ""))
}
