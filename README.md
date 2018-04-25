ku.enm vignette
================
Marlon E. Cobos and A. Townsend Peterson
March 15, 2018

-   [Introduction](#introduction)
-   [Getting started](#getting-started)
    -   [Directory structure and necessary data](#directory-structure-and-necessary-data)
    -   [Installing the package](#installing-the-package)
    -   [Downloading the example data](#downloading-the-example-data)
-   [Doing the analyses](#doing-the-analyses)
    -   [Workflow recording](#workflow-recording)
    -   [Calibration of models](#calibration-of-models)
        -   [Candidate models creation](#candidate-models-creation)
        -   [Evaluation, selection of best models](#evaluation-selection-of-best-models)
    -   [Final model creation](#final-model-creation)
    -   [Final model evaluation](#final-model-evaluation)
    -   [Extrapolation risk analysis](#extrapolation-risk-analysis)

<br>

Introduction
------------

**ku.enm** is an R package designed to make the process of model calibration and final model creation easier and more reproducible, and at the same time more robust. The aim of this package is to desing suites of candidate models to create calibrations of Maxent models that will select optimal parametrizations for each study. Other objectives of this program are to make the task of creating final models and their transfers easier, as well as presenting tools for assessing extrapolation risks when model transfers are needed.

This document is a brief tutorial for using the functions of the **ku.enm** R package. An example of a disease vector is used in this tutorial to make it more clear and understandable. Functions help can be checked while performing the processes.

<br>

Getting started
---------------

### Directory structure and necessary data

Since this package was designed to perform complex analyses while avoiding excessive demands on the computer (especially related to RAM memory used for R), it needs certain data and organization in the working directory. Following this structure (Fig.1) will allow working with one or more species in a project and avoid potential problems during the analyses.

Before starting the analyses, make sure your working directory has the following components:

-   A folder containing the distinct sets of environmental variables (i.e., M\_variables in Figure 1) that are going to be used (more than one is highly recommended, but not mandatory). These variables must represent the area over which your models are going to be calibrated.
-   The maxent.jar application, available from the <a href="https://biodiversityinformatics.amnh.org/open_source/maxent/" target="_blank">Maxent repository</a>.
-   A csv file containing training and testing occurrence data together (preferably after cleaning and thinning original data to avoid problems like wrong records or spatial auto-correlation). This data set consists of three fields: species name, longitud, and latitude. See Sp\_joint.csv, in figure 1.
-   A csv file containing occurrence data for training models. This and the next files generally represent exclusive subsets of the complete set of records. Occurrences can be subsetted in multiple ways, but independence of training and testing data is desired. See Sp\_train.csv in figure 1.
-   A csv file containing species occurrence data for testing models as part of the calibration process (i.e., Sp\_test.csv in figure 1).
-   A csv file containing a completely independent subset of occurrence data--external to training and testing data--for a final, formal model evaluation. This data (i.e., for final model evaluation) is given as Sp\_ind.csv, in figure 1.

Another important requirement for using Maxent and therefore the ku.enm package is to have Java installed in the computer. Java can be downloaded from the <a href="https://java.com/es/download/" target="_blank">Java download page</a>.

<br>

<img src="Structure.png" alt="Figure 1. Directory structure and data for starting and when finished the processes using the ku.enm R package. Background colors represent necessary data before starting the analyses (blue) and data generated after: using the start function (yellow), creating candidate models (lighter green), evaluating candidate models (purple), preparing projection layers (light orange), generating final models and its projections (light grey), evaluating final models with independent data (brown), and analyzing extrapolation risks in projection areas or scenarios (darker green)." width="4114" />
<p class="caption">
Figure 1. Directory structure and data for starting and when finished the processes using the ku.enm R package. Background colors represent necessary data before starting the analyses (blue) and data generated after: using the start function (yellow), creating candidate models (lighter green), evaluating candidate models (purple), preparing projection layers (light orange), generating final models and its projections (light grey), evaluating final models with independent data (brown), and analyzing extrapolation risks in projection areas or scenarios (darker green).
</p>

<br>

### Installing the package

The **ku.enm** R package is in a GitHub repository and can be installed and/or loaded using the following code (make sure to have internet connection).

``` r
if(!require(ku.enm)){
    devtools::install_github("manubio13/ku.enm")
    library(ku.enm)
}
```

<br>

### Downloading the example data

Data used as an example for testing this package correspond to the turkey tick *Amblyomma americanum*, a vector of various diseases, including human monocytotropic ehrlichiosis, canine and human granulocytic ehrlichiosis, tularemia, and southern tick-associated rash illness. This species is distributed in North America and a complete analysis of the risk of its invasion in other regions is being developed by Raghavan et al. (in prep.).

These data are already structured as needed for doing analysis with this package, and can be downloaded (from <a href="https://biodiversityinformatics.amnh.org/open_source/maxent/" target="_blank">ku.enm example data</a>) and extracted using the code below.

``` r
download.file(url = "url", destfile = "C:/Users/###YOUR_USER###/Documents/ku.enm_example_data.zip", 
              quiet = FALSE) #donwload the zipped example folder

unzip(zipfile = "C:/Users/###YOUR_USER###/Documents/ku.enm_example_data.zip",
      exdir = "C:/Users/###YOUR_USER###/Documents") #unzip the example folder

unlink("C:/Users/###YOUR_USER###/Documents/ku.enm_example_data.zip") #erase zip file

setwd("C:/Users/###YOUR_USER###/Documents/ku.enm_example_data/A_americanum") #set the working directory

dir() #check what is in your working directory

# If you have your own data and they are organized as in the first part of Figure 1, change 
# your directory and follow the instructions below.
```

Your working directory will be structured similarly to that presented in Figure 1.

<br>

Doing the analyses
------------------

### Workflow recording

Once the working diretory and data are ready, the function *ku.enm.start* (.Rmd) will allow generating an *R Markdown* file as a guide to perform all the analyses that this package includes (Fig 1. yellow area). By recording all the code chunks to be used in the process, this file also helps to make analyses more reproducible. This file will be written in the working directory. The usage of this function is optional, but recomended to treat each species as a single case.

``` r
help("ku.enm.start")
```

``` r
#Preparing variables to be used in arguments
file_name = "aame_enm_process"
```

``` r
ku.enm.start(file.name = file_name)
```

<br>

### Calibration of models

Notice that, from this point, the following procedures will be performed in the *R Markdown* file previously created, but only if the *ku.enm.start* function was used.

#### Candidate models creation

The function *ku.enm.cal* creates and executes a batch file for generating Maxent candidate models that will be written in sub-directories, named as the parameterizations selected, inside the output directory (Fig. 2, light green area). Calibration models will be created with multiple combinations of regularization multipliers, feature classes, and sets of environmental predictors. For each combination, this function creates one Maxent model with the complete set of occurrences and another with training occurrences only. On some computers, the user will be asked if ruining the batch file is allowed before starting the modeling process in Maxent.

Maxent will run in command-line interface (**do not close the application**) and its graphic interface will not show up, to avoid interfering with activities other than the modeling process.

``` r
help("ku.enm.cal")
```

``` r
#Variables with information to be used as arguments
occ_joint <- "aame_joint.csv"
occ_tra <- "aame_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "all"
run <- TRUE
```

``` r
ku.enm.cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
           out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, run = run)
```

<br>

#### Evaluation, selection of best models

The function *ku.enm.ceval* evaluates model performance based on statistical significance (partial ROC), omission rate (E = user-selected%), and complexity (AICc), and selects best models based on distinct criteria (see selection in function help). Partial ROC and omission rates are evaluated on models created with training occurrences, and AICc values are calculated on models created with the full set of occurrences. The outputs are stored in a folder that will contain a csv file with the statistics of models meeting various criteria, another with only the models selected based on the user's criteria, a third with performance metrics for all candidate models, a plot of the model performance based on the metrics, and an HTML file reporting all the results of the model evaluation and selection process designed to guide further interpretations (Fig. 2, purple area).

``` r
help("ku.enm.ceval")
```

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

``` r
ku.enm.ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
             out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
             kept = kept, selection = selection)
```

<br>

### Final model creation

After selecting parameterizations producing the best models, the next step is that of creating the final models and, if needed, transferring them to other areas or scenarios. The *ku.enm.mod* function takes the csv file with the best models from the model selection process, and writes and executes a batch file for creating final models with the selected parameterizations. Models and projections are stored in subdirectories inside an output folder and these subdirectories will be named as with the candidate models. By allowing projections (i.e., project = TRUE) and defining the folder holding the data for transfers (i.e., folder name in G.var.dir argument), this function automatically performs those transfers.

Maxent will run in command-line interface as it did when creating the calibration models (**again, do not close the application**). However, the process of creating final models may take considerably more time, especially when transfering.

``` r
help("ku.enm.mod")
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
run1 <- TRUE
# Again, some of the variables used here as arguments were already created for the previous functions
```

``` r
ku.enm.mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
           rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
           G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, run = run1)
```

<br>

### Final model evaluation

Final models should be evaluated using independent occurrence data (i.e., data that has not being used during the calibration process that usually come from different sources). The *ku.enm.feval* function evaluates final models based on statistical significance (partial ROC) and omission rate (E = user-selected%). This function will return a folder containing a csv file with the results of the evaluation (see Figure 1, brown color).

``` r
help("ku.enm.feval")
```

``` r
occ_ind <- "aame_ind.csv"
replicates <- TRUE
out_feval <- "Final_Models_evaluation"
# Most of the variables used here as arguments were already created for the previous functions
```

``` r
ku.enm.feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
             out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
             iterations = iterations)
```

<br>

### Extrapolation risk analysis

If transfers were performed at the moment of creating final models, risks of extrapolation can be assessed using the *ku.enm.mmop* function. This function calculates Mobility-Oriented Parity (MOP) layers by comparing environmental values between the M area and one or multiple areas or scenarios to which ecological niche models were transferred. The layers produced with this function help to visualize were there is strict extrapolation risks and different similarity levels between the projection areas or scenarios and the M (calibration area). Results from this analysis will be written for each set of variables inside an specific data (Figure 1, dark green areas).

``` r
help("ku.enm.mmop")
```

``` r
sets_var <- "Set3"
out_mop <- "MOP_results"
percent <- 10
normalized <- TRUE
# Two of the variables used here as arguments were already created for the previous functions
```

``` r
ku.enm.mmop(dirG = G_var_dir, dirM = M_var_dir, sets.var = sets_var, out.mop = out_mop,
            percent = percent, normalized = normalized)
```

<br>
