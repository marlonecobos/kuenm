kuenm vignette
================
Marlon E. Cobos and A. Townsend Peterson.
April 25, 2018

-   [Introduction](#introduction)
-   [Getting started](#getting-started)
    -   [Directory structure and necessary data](#directory-structure-and-necessary-data)
    -   [Installing the package](#installing-the-package)
    -   [Downloading the example data](#downloading-the-example-data)
-   [Doing analyses for a single species project](#doing-analyses-for-a-single-species-project)
    -   [Workflow recording](#workflow-recording)
    -   [Calibration of models](#calibration-of-models)
        -   [Creation of candidate models](#creation-of-candidate-models)
        -   [Evaluation and selection of best models](#evaluation-and-selection-of-best-models)
    -   [Final model creation](#final-model-creation)
    -   [Final model evaluation](#final-model-evaluation)
    -   [Extrapolation risk analysis](#extrapolation-risk-analysis)
    -   [References](#references)

<br>

Introduction
------------

**kuenm** is an R package designed to make the process of model calibration and final model creation easier and more reproducible, and at the same time more robust. The aim of this package is to design suites of candidate models to create diverse calibrations of Maxent models to allow selection of optimal parameterizations for each study. Other objectives of this program are to make the task of creating final models and their transfers easier, as well to permit assessing extrapolation risks when model transfers are needed.

This document is a brief tutorial for using the functions of the **kuenm** R package. An example of a disease vector (a tick) is used in this tutorial to make it more clear and understandable. Functions help can be checked while performing the processes.

<br>

Getting started
---------------

### Directory structure and necessary data

Since this package was designed to perform complex analyses while avoiding excessive demands on the computer (especially related to RAM memory used for R), it needs certain data and to be organized carefully in the working directory. Following this structure (Figure 1) will allow working with one or more species in a project, and avoid potential problems during the analyses.

Before starting the analyses, the user must make sure that the working directory (the folder with information for an individual species) has the following components:

-   A folder containing the distinct sets of environmental variables (i.e., M\_variables in Figure 1) to be used (more than one is highly recommended, but not mandatory). These variables must represent the area over which models are calibrated.
-   A csv file containing training and testing occurrence data together (preferably after cleaning and thinning original data to avoid problems like wrong records and spatial auto-correlation). This data set consists of three fields: species name, longitude, and latitude. See Sp\_joint.csv, in Figure 1.
-   A csv file containing occurrence data for training models. This file and the next file generally represent exclusive subsets of the full set of records. Occurrences can be subsetted in multiple ways (Muscarella et al. 2014), but some degree of independence of training and testing data is desired. See Sp\_train.csv in figure 1.
-   A csv file containing species occurrence data for testing models as part of the calibration process (i.e., Sp\_test.csv in figure 1).
-   If available, a csv file containing a completely independent subset of occurrence data—external to training and testing data—for a final, formal model evaluation. This dataset (i.e., for final model evaluation) is given as Sp\_ind.csv, in Figure 1.

A crucial requirement is to have the maxent.jar application in any user-defined directory (we encourage you to maintain it in a fixed directory). This software is available in the <a href="https://biodiversityinformatics.amnh.org/open_source/maxent/" target="_blank">Maxent repository</a>. Another important requirement for using Maxent and therefore the kuenm package is to have the Java Development Kit installed in the computer. The Java Development Kit is available in <a href="http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html" target="_blank">this repository</a>.

<br>

![Figure 1. Directory structure and data for starting processing, as well as directory structure when the processes finish using the kuenm R package. Background colors represent data necessary before starting the analyses (blue) and data generated after the following steps: running the start function (yellow), creating candidate models (lighter green), evaluating candidate models (purple), preparing projection layers (light orange), generating final models and its projections (light gray), evaluating final models with independent data (brown), and analyzing extrapolation risks in projection areas or scenarios (darker green).](Structure.png)

<br>

### Installing the package

The **kuenm** R package is in a GitHub repository and can be installed and/or loaded using the following code (make sure to have Internet connection).

``` r
# Installing and loading packages
if(!require(devtools)){
    install.packages("devtools")
}

if(!require(kuenm)){
    devtools::install_github("marlonecobos/kuenm")
}

library(kuenm)
```

<br>

### Downloading the example data

Data used as an example for testing this package correspond to the turkey tick *Amblyomma americanum*, a vector of various diseases, including human monocytotropic ehrlichiosis, canine and human granulocytic ehrlichiosis, tularemia, and southern tick-associated rash illness. This species is distributed broadly in North America and a complete analysis of the risk of its invasion in other regions is being developed by Raghavan et al. (in review).

These data are already structured as needed for doing analysis with this package, and can be downloaded (from <a href="http://doi.org/10.17161/1808.26376" target="_blank">kuenm example data</a>) and extracted using the code below.

``` r
# Change "YOUR/DIRECTORY" by your actual directory.
download.file(url = "https://kuscholarworks.ku.edu/bitstream/handle/1808/26376/ku.enm_example_data.zip?sequence=3&isAllowed=y", 
              destfile = "YOUR/DIRECTORY/ku.enm_example_data.zip", mode = "wb",
              quiet = FALSE) # donwload the zipped example folder in documents

unzip(zipfile = "YOUR/DIRECTORY/ku.enm_example_data.zip",
      exdir = "YOUR/DIRECTORY") # unzip the example folder in documents

unlink("YOUR/DIRECTORY/ku.enm_example_data.zip") # erase zip file

setwd("YOUR/DIRECTORY/ku.enm_example_data/A_americanum") # set the working directory

dir() # check what is in your working directory

# If you have your own data and they are organized as in the first part of Figure 1, change 
# your directory and follow the instructions below.
```

Your working directory will be structured similarly to that presented for Species\_1 in the left part of Figure 1.

<br>

Doing analyses for a single species project
-------------------------------------------

### Workflow recording

Once the working directory and data are ready, the function *kuenm\_start* (.Rmd) will allow generating an *R Markdown* file as a guide to performing all the analyses that this package includes (Figure 1, yellow area). By recording all the code chunks used in the process, this file also helps to make analyses more reproducible. This file will be written in the working directory. The usage of this function is optional, but it is recomended if recording individual workflows per each species is desired.

``` r
help(kuenm_start)
```

``` r
# Preparing variables to be used in arguments
file_name <- "aame_enm_process"
```

``` r
kuenm_start(file.name = file_name)
```

<br>

### Calibration of models

Notice that, from this point, the following procedures will be performed in the *R Markdown* file previously created, but only if the *kuenm\_start* function was used.

<br>

#### Creation of candidate models

The function *kuenm\_cal* creates and executes a batch file for generating Maxent candidate models that will be written in subdirectories, named as the parameterizations selected, inside the output directory (Figure 1, light green area). Calibration models will be created with multiple combinations of regularization multipliers, feature classes, and sets of environmental predictors. For each combination, this function creates one Maxent model with the complete set of occurrences and another with training occurrences only. On some computers, the user will be asked if ruinning the batch file is allowed before the modeling process starts in Maxent.

Maxent will run in command-line interface (**do not close the application**) and its graphic interface will not show up, to avoid interfering with activities other than the modeling process.

``` r
help(kuenm_cal)
```

``` r
# Variables with information to be used as arguments. Change "YOUR/DIRECTORY" by your actual directory.
occ_joint <- "aame_joint.csv"
occ_tra <- "aame_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "all"
background <- 10000
maxent_path <- "YOUR/DIRECTORY/ku.enm_example_data/A_americanum"
wait <- TRUE
invisible <- TRUE
run <- TRUE
```

``` r
kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, background = background,
          maxent.path = maxent_path, wait = wait, invisible = invisible, run = run)
```

<br>

#### Evaluation and selection of best models

The function *kuenm\_ceval* evaluates model performance based on statistical significance (partial ROC), omission rate (*E* = a user-selected proportion of occurrence data that may present meaninful errors; see Peterson et al. (2008)), and model complexity (AICc), and selects best models based on distinct, user-set criteria (see selection in function help). Partial ROC and omission rates are evaluated based on models created with training occurrences, whereas AICc values are calculated for models created with the full set of occurrences (Warren and Seifert 2011). Outputs are stored in a folder that will contain a .csv file with the statistics of models meeting each of the evaluation criteria, another with only the models selected based on the user-specified criteria, a third with performance metrics for all candidate models, a plot of model performance, and an HTML file reporting all the results of the model evaluation and selection process designed to guide further interpretations (Figure 1, purple area).

``` r
help(kuenm_ceval)
```

``` r
occ_test <- "aame_test.csv"
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 100
kept <- TRUE
selection <- "OR_AICc"
# Note, some of the variables used here as arguments were already created for previous function
```

``` r
cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection)
```

<br>

### Final model creation

After selecting parameterizations producing the best models, the next step is that of creating final models and, if needed, transferring them to other areas or scenarios. The *kuenm\_mod* function takes the .csv file with the best models from the model selection process, and writes and executes a batch file for creating final models with the selected parameterizations. Models and projections are stored in subdirectories inside an output folder; these subdirectories will be named as with the candidate models. By allowing projections (i.e., project = TRUE) and defining the folder holding the data for transfers (i.e., folder name in G.var.dir argument), this function automatically performs those transfers.

Maxent will run in command-line interface, as it did when creating the calibration models (**again, do not close the application**). However, the process of creating final models may take considerably more time, especially when transferring to other regions or scenarios.

``` r
help(kuenm_mod)
```

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
wait1 <- TRUE
invisible1 <- TRUE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
             # "outputgrids=false" which avoids writing grids of replicated models and only writes the 
             # summary of them (e.g., average, median, etc.) when rep.n > 1
             # note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions
```

``` r
kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, invisible = invisible1, run = run1)
```

<br>

### Final model evaluation

Final models should be evaluated using independent occurrence data (i.e., data that have not been used in the calibration process that usually come from different sources). The *kuenm\_feval* function evaluates final models based on statistical significance (partial ROC) and omission rate (*E*). This function will return a folder containing a .csv file with the results of the evaluation (see Figure 1, brown color).

``` r
help(kuenm_feval)
```

``` r
occ_ind <- "aame_ind.csv"
replicates <- TRUE
out_feval <- "Final_Models_evaluation"
# Most of the variables used here as arguments were already created for previous functions
```

``` r
fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations)
```

<br>

### Extrapolation risk analysis

If transfers were performed when creating final models, risks of extrapolation can be assessed using the *kuenm\_mmop* function. This function calculates mobility-oriented parity (MOP) layers (Owens et al. 2013) by comparing environmental values between the calibration area and one or multiple regions or scenarios to which ecological niche models were transferred. The layers produced with this function help to visualize were strict extrapolation risks exist, and different similarity levels between the projection regions or scenarios and the calibration area. Results from this analysis will be written for each set of variables inside an specific data (Figure 1, dark green areas).

``` r
help(kuenm_mmop)
```

``` r
sets_var <- "Set3" # here a vector of various sets can be used
out_mop <- "MOP_results"
percent <- 10
normalized <- TRUE
# Two of the variables used here as arguments were already created for previous functions
```

``` r
kuenm_mmop(G.var.dir = G_var_dir, M.var.dir = M_var_dir, sets.var = sets_var, out.mop = out_mop,
           percent = percent, normalized = normalized)
```

<br>

### References

Muscarella, R. et al. 2014. ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. - Methods in Ecology and Evolution 5: 1198–1205.

Owens, H. L. et al. 2013. Constraints on interpretation of ecological niche models by limited environmental ranges on calibration areas. - Ecological Modelling 263: 10–18.

Peterson, A. T. et al. 2008. Rethinking receiver operating characteristic analysis applications in ecological niche modeling. - Ecological Modelling 213: 63–72.

Warren, D. L. and Seifert, S. N. 2011. Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria. - Ecological Applications 21: 335–342.
