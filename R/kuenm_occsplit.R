#' Split occurrence files in training and testing data
#'
#' @description kuenm_occsplit splits occurrences contained in a data.frame to obtain training
#' and testing data based on distinct methods for calibrating models.
#'
#' @param occ data.frame of occurrence records containing at least species,
#' longitude, and latitude columns.
#' @param train.proportion (numeric) proportion (from 0 to 1) of data to be used as training
#' occurrences. The remaining data will be used for testing.
#' @param method (character) method for selecting training and testing occurrences. Current
#' option is "random".
#' @param save (logical) whether or not to save the results in the working
#' directory. Default = FALSE.
#' @param name (character) common name for csv files to be written. A suffix will be added
#' depending on if the data is the complete set, training set, or testing set of occurrences.
#'
#' @return
#' List with all, training, and testing occurrences. Three csv files will be written in the
#' working directory according to the name defined in \code{name} plus the suffix _joint
#' for all records, _train for the training set, and _test for the testing set.
#'
#' @usage
#' kuenm_occsplit(occ, train.proportion = 0.5, method = "random",
#'                save = FALSE, name = "occ")
#'
#' @export
#'
#' @examples
#' # arguments
#' data("sp_joint", package = "kuenm")
#'
#' occs <- data.frame(Species = "A_americanum", sp_joint)
#' train_prop <- 0.5
#' method = "random"
#'
#' # running
#' data_split <- kuenm_occsplit(occ = occs, train.proportion = train_prop,
#'                              method = method)

kuenm_occsplit <- function(occ, train.proportion = 0.5, method = "random",
                           save = FALSE, name = "occ") {

  if (missing(occ)) {stop("Argument 'occ' needs to be defined.")}

  occ <- na.omit(occ)

  if (method == "random") {
    files <- occ_randsplit(occ, train.proportion = train.proportion)
  }

  if (save == TRUE) {
    names <- paste0(name, c("_joint", "_train", "_test"), ".csv")
    wrt <- sapply(1:length(files), function(x){
      write.csv(files[[x]], file = names[x], row.names = FALSE)
    })
  }

  return(files)
}


#' Split occurrences randomly in training and testing data
#'
#' @description occ_randsplit splits a set of occurrences to obtain training and testing
#' data randomly.
#'
#' @param occ matrix or data.frame with the occurrences to be split. Columns may vary but
#' species, longitude, and latitue are recommended.
#' @param train.proportion (numeric) proportion (from 0 to 1) of data to be used as training
#' occurrences. The remaining data will be used for testing.
#'
#' @return
#' List with all occurrences (joint), training occurrences (train), and testing (test)
#' occurrences.
#'
#' @usage
#' occ_randsplit(occ, train.proportion = 0.5)
#'
#' @export
#'
#' @examples
#' # arguments
#' occs <- read.csv(list.files(system.file("extdata", package = "kuenm"),
#'                             pattern = "sp_test.csv", full.names = TRUE))
#' occs <- data.frame(Species = "Species_1", occs)
#' train_prop <- 0.5
#'
#' # running
#' occ_rsplit <- occ_randsplit(occ = occs, train.proportion = train_prop)

occ_randsplit <- function(occ, train.proportion = 0.5) {
  ndata <- nrow(occ)
  ids <- sample(ndata, size = round(train.proportion * ndata))
  data <- list(joint = occ, train = occ[ids, ], test = occ[-ids, ])

  return(data)
}
