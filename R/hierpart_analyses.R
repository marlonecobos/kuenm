#' Performing hierarchical partitioning analyses
#'
#' @description Helper function to perform hierarchical partitioning analyses
#' of the varianze comming from distinct sources in ecological niche models.
#'
#' @param tables.folder (character) name of the folder where the tables created with the
#' function {\link{hierpart_tables}} are.
#' @param out.dir (character) name of the output directory where results will be written.
#'
#' @return All the results from performing hierarchical partitioning analyses of the
#' varianze in ecological niche models.
#'
#' @export

hierpart_analyses <- function(tables.folder, out.dir, kept = FALSE) {

  # installing needed packages if required
  #install.packages("hier.part", dependencies = TRUE)

  fls <- list.files(tables.folder, pattern = "*.csv", full.names = TRUE) #Files to be used

  r1_hpar_ran <- list()
  r2_hpar_ran <- list()
  r3_hpar_ran <- list()

  for (i in 1:length(fls)) {
    ra <- read.csv(fls[i])
    ran <- as.data.frame(ra)
    ran[, 2:dim(ran)[2]] <- lapply(ran[, 2:dim(ran)[2]], factor)
    ran[, 1] <- (ran[, 1] + 1) # adding 1 to all the Suitability to use gamma family

    # Hierarchical Partitioning analyses
    hpar_ran <- hier.part::hier.part(ran[, 1], ran[, 2:dim(ran)[2]], fam = "Gamma",
                                     gof = "logLik", barplot = F)
    r1_hpar_ran[[i]] <- hpar_ran$gfs
    r2_hpar_ran[[i]] <- data.frame(source = rownames(hpar_ran$IJ), hpar_ran$IJ,
                                   row.names = c())
    r3_hpar_ran[[i]] <- hpar_ran$IJ$Total * 100 / (sum(hpar_ran$IJ$Total))

    cat(paste("   ", i, "of", length(fls), "analyses finished\n"))
  }

  # Create the final tables
  r1_hpar_ran <- do.call(rbind, r1_hpar_ran)
  r2_hpar_ran <- do.call(rbind, r2_hpar_ran)
  r2_hpar_ran <- data.frame(rep(paste("Iteration", 1:length(fls)), each = dim(ra)[2] - 1),
                            r2_hpar_ran) # names of the rows
  colnames(r2_hpar_ran) <- c("", "Source", "Independent", "Joint", "Total")

  r3_hpar_ran <- do.call(rbind, r3_hpar_ran)
  colnames(r3_hpar_ran) <- rownames(hpar_ran$IJ)

  r4_hpar_ran <- list()
  for (i in 1:length(rownames(hpar_ran$IJ))) {
    r4_hpar_ran[[i]] <- apply(r2_hpar_ran[r2_hpar_ran$Source == rownames(hpar_ran$IJ)[i],
                                          3:dim(r2_hpar_ran)[2]], 2, mean)
  }

  r4_hpar_ran <- cbind(rownames(hpar_ran$IJ), do.call(rbind, r4_hpar_ran))
  colnames(r4_hpar_ran) <- c("Source", "Independent", "Joint", "Total")

  row.names(r1_hpar_ran) <- paste("Iteration", 1:length(fls)) # names of the rows
  row.names(r3_hpar_ran) <- paste("Iteration", 1:length(fls)) # names of the rows

  # Write the ressults
  write.csv(r1_hpar_ran, paste(out.dir, "hierpart_Goodness_fit.csv", sep = "/"),
            row.names = TRUE)
  write.csv(r2_hpar_ran, paste(out.dir, "hierpart_Raw_effects.csv", sep = "/"),
            row.names = FALSE)
  write.csv(r3_hpar_ran, paste(out.dir, "hierpart_Total_effects_percent.csv", sep = "/"),
            row.names = TRUE)
  write.csv(r4_hpar_ran, paste(out.dir, "hierpart_Mean_effects.csv", sep = "/"),
            row.names = FALSE)

  if (kept == FALSE) {
    unlink(tables.folder, recursive = TRUE, force = TRUE)
  }

  return(r3_hpar_ran)
}
