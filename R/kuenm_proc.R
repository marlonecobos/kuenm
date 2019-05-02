#' Partial ROC calculation of single models
#'
#' @description kuenm_proc applies partial ROC tests to single models.
#'
#' @param occ.test a numerical matrix containing coordinates of the occurrences used to test
#' the ecological niche model to be evaluated; columns must be: longitude and latitude.
#' @param model a RasterLayer of the ecological niche model to be evaluated.
#' @param threshold (numeric) value from 0 to 100 that will be used as threshold (E);
#' default = 5.
#' @param rand.percent (numeric) value from 0 to 100 representing the percent of testing data
#' to be used for performing the bootstrap process for calculating the partial ROC;
#' default = 50.
#' @param iterations (numeric) number of bootstrap iterations to be performed;
#' default = 500.
#' @param parallel (logical) if TRUE, calculations will be performed in parallel using the available
#' cores of the computer. This will demand more RAM and almost full use of the CPU; hence, its use
#' is more recommended in high-performance computers. Using this option will speed up the analyses
#' only if \code{model} is a large RasterLayer or if \code{iterations} are more than 5000.
#' Default = FALSE.
#'
#' @return A data.frame containing the AUC values and AUC ratios calculated for each iteration.
#'
#' @details Partial ROC is calculated following Peterson et al.
#' (2008; \url{http://dx.doi.org/10.1016/j.ecolmodel.2007.11.008}). This function is a modification
#' of the \code{\link[ntbox]{pROC}} funcion, available at \url{https://github.com/luismurao/ntbox}.
#'
#' @importFrom purrr map_df
#' @useDynLib kuenm
#' @export
#'
#' @examples
#' occ <- read.csv(list.files(system.file("extdata", package = "kuenm"),
#'                            pattern = "sp_test.csv", full.names = TRUE))
#' model <- raster::raster(list.files(system.file("extdata", package = "kuenm"),
#'                                    pattern = "sp_model.tif", full.names = TRUE))
#' thres <- 5
#' rand_perc <- 50
#' iterac <- 500
#'
#' p_roc <- kuenm_proc(occ.test = occ, model = model, threshold = thres,
#'                    rand.percent = rand_perc, iterations = iterac)

kuenm_proc <- function(occ.test, model, threshold = 5, rand.percent = 50,
                       iterations = 500, parallel = FALSE){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(future))

  if(model@data@min == model@data@max){
    warning("\nModel with no variability, pROC will return NA.\n")

    p_roc <- rep(NA, 2)
    names(p_roc) <- c(paste("Mean_AUC_ratio_at_", threshold, "%", sep = ""), "pval_pROC")

    auc_ratios <- rep(NA, 3)
    names(auc_ratios) <- c("Model_partial_AUC", "Random_curve_partial_AUC", "AUC_ratio")

    p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)

    return(p_roc_res)
  }else {
    model <- round((model/raster::cellStats(model, max)) * 1000)
    test_value <- na.omit(raster::extract(model, occ.test))
    #test_value <- unique(test_value)
    classpixels <- data.frame(raster::freq(model))
    classpixels <- data.frame(stats::na.omit(classpixels))

    classpixels <- classpixels  %>% dplyr::mutate_(value= ~rev(value),
                                                   count= ~rev(count),
                                                   totpixperclass = ~cumsum(count),
                                                   percentpixels= ~ totpixperclass/sum(count)) %>%
      dplyr::arrange(value)


    error_sens <- 1-(threshold/100)
    models_thresholds <- classpixels[,"value"]
    fractional_area <- classpixels[,"percentpixels"]
    n_data <- length(test_value)
    n_samp <- ceiling((rand.percent/100)*(n_data))

    big_classpixels <- matrix(rep(models_thresholds,each=n_samp),
                              ncol=length(models_thresholds))


    calc_aucDF <- function(big_classpixels,fractional_area,
                           test_value,n_data,n_samp,error_sens){

      rowsID <- sample(x = n_data,
                       size = n_samp,
                       replace=TRUE)
      test_value1 <- test_value[rowsID]
      omssion_matrix <-   big_classpixels >  test_value1
      sensibility <- 1 - colSums(omssion_matrix)/n_samp
      xyTable <- data.frame(fractional_area,sensibility)
      less_ID <- which(xyTable$sensibility<=error_sens)
      xyTable <- xyTable[-less_ID,]

      xyTable <- xyTable[order(xyTable$fractional_area,
                               decreasing = F),]

      auc_pmodel <- trap_roc(xyTable$fractional_area,
                             xyTable$sensibility)

      auc_prand <- trap_roc(xyTable$fractional_area,
                            xyTable$fractional_area)
      auc_ratio <- auc_pmodel/auc_prand

      auc_table <- data.frame(auc_pmodel,
                              auc_prand,
                              auc_ratio =auc_ratio )
      return(auc_table)

    }


    if(parallel){

      future::plan(future::multiprocess)
      roc_env <- new.env()
      n_cores <- future::availableCores()
      niter_big <- floor(iterations/n_cores)
      n_runs <- rep(niter_big,n_cores)
      sum_n_runs <- sum(n_runs)
      n_runs[1] <- n_runs[1] + (iterations - sum_n_runs)

      for(i in 1:length(n_runs)){
        x <- as.character(i)
        roc_env[[x]] %<-% {
          x1 <- 1:n_runs[i]
          auc_matrix1 <- x1 %>%
            purrr::map_df(~calc_aucDF(big_classpixels,
                                      fractional_area,
                                      test_value,n_data,n_samp,
                                      error_sens))
        }
      }
      partial_AUC <- as.list(roc_env)
      rm(roc_env)
      partial_AUC <- do.call(rbind.data.frame,partial_AUC)
      rownames(partial_AUC) <- NULL
      future::plan(future::sequential)

    }else{

      partial_AUC <- 1:iterations %>%
        purrr::map_df(~calc_aucDF(big_classpixels,
                                  fractional_area,
                                  test_value,n_data,n_samp,
                                  error_sens))

    }

    mauc <- mean(partial_AUC$auc_ratio, na.rm = TRUE)
    #proc <- sum(partial_AUC$auc_ratio <= 1) / length(partial_AUC$auc_ratio)
    proc <- (sum(partial_AUC$auc_ratio[!is.na(partial_AUC$auc_ratio)] <= 1) + sum(is.na(partial_AUC$auc_ratio)))/
      length(partial_AUC$auc_ratio)
    p_roc <- c(mauc, proc)
    names(p_roc) <- c(paste("Mean_AUC_ratio_at_", threshold, "%", sep = ""), "pval_pROC")

    auc_ratios <- partial_AUC
    names(auc_ratios) <- c("Model_partial_AUC", "Random_curve_partial_AUC", "AUC_ratio")

    p_roc_res <- list(pROC_summary = p_roc, pROC_results = auc_ratios)
    return(p_roc_res)
  }
}
