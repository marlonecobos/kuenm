#' Split occurrence files in training and testing data
#'
#' @description kuenm_occsplit splits occurrences contained in a csv file to obtain training
#' and testing data based on distinct methods for calibrating models.
#'
#' @param occ_file (character) name of the csv file with all the occurrences; columns must be:
#' species, longitude, latitude.
#' @train.proportion (numeric) proportion (from 0 to 1) of data to be used as training
#' occurrences. The remaining data will be used for testing.
#' @param method (character) method for selecting training and testing occurrences. Current
#' option is "random".
#' @param name (character) common name for csv files to be written. A suffix will be added
#' depending on if the data is the complete set, training set, or testing set of occurrences.
#'
#' @return
#' List with all, training, and testing occurrences. Three csv files will be written in the
#' working directory according to the name defined in \code{name} plus the suffix _joint
#' for all records, _train for the training set, and _test for the testing set.
#'
#' @export
#'
#' @examples
#' # arguments
#' occs <- "occurrences.csv"
#' train_prop <- 0.5
#' method = "random"
#' name <- "occ"
#'
#' # running
#' data_split <- kuenm_occsplit(occ.file = occs, train.proportion = train_prop,
#'                              method = method, name = name)

kuenm_occsplit <- function(occ.file, train.proportion = 0.5, method = "random", name = "occ") {

  if (!file.exists(occ.file)) {
    stop(paste(occ.file, "does not exist in the working directory, check file name",
               "\nor extension, example: species_all.csv"))
  }

  occ <- na.omit(read.csv(occ.file))

  if (method == "random") {
    files <- occ_randsplit(occ, train.proportion = train.proportion)
  }

  names <- paste0(name, c("_joint", "_train", "_test"))
  wrt <- sapply(1:length(files), function(x){
    write.csv(files[[x]], file = names[x], row.names = FALSE)
    })

  return(files)
}


#' Split occurrences randomly in training and testing data
#'
#' @description occ_randsplit splits a set of occurrences to obtain training and testing
#' data randomly.
#'
#' @param occ matrix or data.frame with the occurrences to be split. Columns may vary but
#' species, longitude, and latitue are recommended.
#' @train.proportion (numeric) proportion (from 0 to 1) of data to be used as training
#' occurrences. The remaining data will be used for testing.
#'
#' @return
#' List with all occurrences (joint), training occurrences (train), and testing (test)
#' occurrences.
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
