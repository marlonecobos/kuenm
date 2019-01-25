## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Fig.1, echo=FALSE, message=FALSE, warning=FALSE, fig.width=4, fig.cap="Figure 1. Directory structure and data for starting using the R functions for the assessment of variation in ecological niche model outputs. This figure represents models that were created using two distinct parameters settings and were projected to a bigger area in the current time and to future scenarios. Future conditions are represented by two climate models (GCM) in two emission scenarios (RCP). E, EC, and NE, represent three distinct options of extrapolation, free extrapolation, extrapolation and clamping, and no extrapolation, respectively."----

knitr::include_graphics("Structure_variation.png")

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  help(kuenm_modstats)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  sp_name <- "sp1"
#  fmod_dir <- "Final_Models"
#  format <- "asc"
#  project <- TRUE
#  stats <- c("med", "range")
#  rep <- TRUE
#  scenarios <- c("current", "GCM1_RCP4.5", "GCM1_RCP8.5", "GCM2_RCP4.5", "GCM2_RCP8.5")
#  ext_type <- c("E", "EC", "NE") # the type of extrapolation can be selected according to user requirements
#  out_dir <- "Final_Model_Stats"
#  
#  # argument "time.periods" is not included in the example but it can be used when models
#  # are projected to more than one time period, other than current.

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  kuenm_modstats(sp.name = sp_name, fmod.dir = fmod_dir, format = format, project = project,
#                 statistics = stats, replicated = rep, proj.scenarios = scenarios,
#                 ext.type = ext_type, out.dir = out_dir)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  help(kuenm_projchanges)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  # other arguments were defined before
#  occ <- "Sp_occ.csv"
#  fmod_stats <- "Final_Model_Stats"
#  thres <- 5
#  curr <- "current"
#  emi_scenarios <- c("RCP4.5", "RCP8.5")
#  c_mods <- c("GCM1", "GCM2")
#  ext_type <- c("E", "EC", "NE")
#  out_dir1 <- "Projection_Changes"

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  kuenm_projchanges(occ = occ, fmod.stats = fmod_stats, threshold = thres, current = curr,
#                    emi.scenarios = emi_scenarios, clim.models = c_mods, ext.type = ext_type,
#                    out.dir = out_dir1)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  help(kuenm_modvar)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  sp_name <- "sp1"
#  fmod_dir <- "Final_Models"
#  rep <- TRUE
#  format <- "asc"
#  project <- TRUE
#  curr <- "current"
#  emi_scenarios <- c("RCP4.5", "RCP8.5")
#  c_mods <- c("GCM1", "GCM2")
#  ext_type <- c("E", "EC", "NE")
#  split <- 100
#  out_dir2 <- "Variation_from_sources"

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  kuenm_modvar(sp.name = sp_name, fmod.dir = fmod_dir, replicated = rep, format = format,
#               project = project, current = curr, emi.scenarios = emi_scenarios,
#               clim.models = c_mods, ext.type = ext_type, split.length = split, out.dir = out_dir2)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  help(kuenm_hierpart)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  sp_name <- "sp1"
#  fmod_dir <- "Final_Models"
#  rep <- TRUE
#  format <- "asc"
#  project <- TRUE
#  curr <- "current"
#  emi_scenarios <- c("RCP4.5", "RCP8.5")
#  c_mods <- c("GCM1", "GCM2")
#  ext_type <- c("E", "EC", "NE")
#  iter <- 100
#  s_size <- 1000
#  out_dir3 <- "Hierarchical_partitioning"
#  # argument "factors_col" is not defined here, but if default colors (grey scale) need to be changed,
#  # you can use this argument.

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  kuenm_hierpart(sp.name = sp_name, fmod.dir = fmod_dir, replicated = rep, format = format,
#                 project = project, current = curr, emi.scenarios = emi_scenarios,
#                 clim.models = c_mods, ext.type = ext_type, iterations = iter,
#                 sample.size = s_size, out.dir = out_dir3)

