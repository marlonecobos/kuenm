*A. ammericanum*: modeling process
================

-   [Candidate models](#candidate-models)
-   [Evaluation and selection of best models](#evaluation-and-selection-of-best-models)
-   [Final model creation](#final-model-creation)
-   [Final model evaluation](#final-model-evaluation)
-   [MOP analysis](#mop-analysis)

The R markdown file is created in the working directory, and is designed to make the processes of model calibration and final model creation more reproducible.

Information on using this R Markdown file:

-   Try executing code chunks by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.
-   Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

A brief tutorial for using functions of the ku.enm R package can be found in the <a href="https://github.com/manubio13/ku.enm" target="_blank">package vignette</a>. Additionally, function help can be checked to change arguments according to specific needs.

<br>

### Candidate models

Candidate models are a large set of candidate models created to respond to the need to test broad suites of parameter combinations, such as, distinct regularization multiplier values, various feature classes, and different sets of environmental variables. The following code calls the help page of the function ku.enm.cal.

``` r
help(ku.enm.cal)
```

<br>

The next chunk of code is for preparing the arguments for using the function following the modularity principle. These variables can be changed according to each case.

``` r
occ_joint <- "aame_joint.csv"
occ_tra <- "aame_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "all"
run <- TRUE
```

<br>

The following is the code for using the function.

``` r
ku.enm.cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
           out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, run = run)
```

<br>

### Evaluation and selection of best models

Evaluation is a crucial step in model calibration. This step centers on selecting candidate models and their associated parameters to identify the best models for the purposes of the study. The ku.enm.eval function evaluates candidate models based on three distinct criteria: statistical significance (based on partial ROC analyses), prediction ability (we use omission rates, but other metrics, such as overall correct classification rate, can also be used), and model complexity (here evaluated using AICc). The following code chunk calls the function help window.

``` r
help(ku.enm.ceval)
```

<br>

Below, arguments for this functions will be defined.

``` r
occ_test <- "aame_test.csv"
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 100
kept <- TRUE
selection <- "OR_AICc"
# Notice, some of the variables used here as arguments were already created for the previous function
```

<br>

This code also allows evaluating candidate models that were created previously, selecting those with best performance based on the three criteria.

``` r
cal_eval <- ku.enm.ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                         out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                         kept = kept, selection = selection)
```

<br>

### Final model creation

After selecting parametrizations that produce best models, the next step is to create the final models, and if needed transfer them to other environmental data sets (e.g., to other time periods or other geographic regions). The function help is called via this code:

``` r
help(ku.enm.mod)
```

<br>

For preparing the arguments for this function use the following chunk of code.

``` r
batch_fin <- "Final_models"
mod_dir <- "Final_Models"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
G_var_dir <- "G_variables"
out_format <- "logistic"
project <- TRUE
ext_type <- "all"
write_mess <- FALSE
write_clamp <- FALSE
run1 <- TRUE
args <- "outputgrids=false"
# Again, some of the variables used here as arguments were already created for the previous functions
```

<br>

The ku.enm.mod function has the following syntax:

``` r
ku.enm.mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
           rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
           G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, run = run1, args = args)
```

<br>

### Final model evaluation

Final model can be evaluated after being created; for this step, independent data are needed (data not used in the calibration process, ideally coming from different sources). The function help is called via this code:

``` r
help(ku.enm.feval)
```

<br>

For preparing the arguments for this function, use the following chunk code.

``` r
occ_ind <- "aame_ind.csv"
replicates <- TRUE
out_feval <- "Final_Models_evaluation"
# Again, some of the variables used here as arguments were already created for the previous functions
```

<br>

The following is the code for using the function.

``` r
fin_eval <- ku.enm.feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                         out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                         iterations = iterations)
```

<br>

### MOP analysis

If transfers were performed when creating final models, the MOP analysis will help to identify areas of strict extrapolation and levels of similarity between the calibration area and the region or scenario of projection. The code below will help to see the function's documentation:

``` r
help(ku.enm.mmop)
```

<br>

Below, arguments for this functions will be defined.

``` r
sets_var <- "Set3"
out_mop <- "MOP_results"
percent <- 10
normalized <- TRUE
# Again, some of the variables used here as arguments were already created for the previous functions
```

<br>

The ku.enm.mmop function has the following syntax:

``` r
ku.enm.mmop(dirG = G_var_dir, dirM = M_var_dir, sets.var = sets_var, out.mop = out_mop,
            percent = percent, normalized = normalized)
```
