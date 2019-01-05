#' A lambdas file resulted from a modeling process in Maxent
#'
#' A lambdas file resulted from a model created in Maxent with raw output for a tick
#' (*Amblyomma americanum*) in North America. This file is used to calculate number
#' of parameters in the model, which is needed while calculating AICc values.
#'
#' @name sp_lambdas
#'
#' @format A lambdas file.
#' \describe{
#'   \item{parameters}{number of parameters in the Maxent model.}
#' }
#'
#' @source \url{https://kuscholarworks.ku.edu/handle/1808/26376}
#'
#' @examples
#' lbds <- readLines(list.files(system.file("extdata", package = "kuenm"), # lambdas file
#'                                           pattern = "lambdas_model_joint.lambdas", full.names = TRUE))
#'
#' head(lbds)
NULL
