#' Helper function to compute partial AUC, AUC ratio.
#' @param iterationno number of iteration to compute partial AUC values
#' @param occurtbl Validation data. Must have 3 columns SpName, Longitude, Latitude.
#' @param rand.percent Occurrence points to be sampled randomly from the test data for bootstrapping.
#' @param omissionval 1-E.
#' @param classpixels Pixel classes.

auc_comp <- function(iterationno,occurtbl,rand.percent,omissionval,classpixels) {
  ClassID <- NULL
  n <- NULL
  OccuSumBelow <- NULL

  if (iterationno > 0) {
    ll <- sample(nrow(occurtbl),
                 round(rand.percent / 100 * nrow(occurtbl)),
                 replace = TRUE)
    occurtbl1 <- occurtbl[ll, ]
  }else {
    occurtbl1 <- occurtbl
  }

  occurinclass <- occurtbl1 %>% dplyr::group_by(ClassID) %>%
    dplyr::count() %>% dplyr::arrange(dplyr::desc(ClassID))
  occurinclass <- occurinclass %>%
    dplyr::ungroup() %>% dplyr::mutate(OccuSumBelow= cumsum(n)) %>%
    dplyr::mutate(Percent = OccuSumBelow / nrow(occurtbl1))
  occurinclass <- as.data.frame(occurinclass)

  names(occurinclass) <- c("ClassID", "OccuCount",
                           "OccuSumBelow", "Percent")

  xytable <- generatexytableb(classpixels,occurinclass)
  arearow <- calculateauc(xytable, omissionval, iterationno)

  names(arearow) <- c("Iteration", paste("AUC_at_", omissionval, "%", sep = ""),
                      paste("AUC_at_", 100 - omissionval, "%", sep = ""), "AUC_ratio")
  return(arearow)
}


#' Helper function to compute the area (number of pixels) that a certain threshold has.
#' @param inrast A raster class object of a continuous model output.

a_pred_pres <- function(inrast) {
  ### Now calculate proportionate area predicted under each suitability
  classpixels <- raster::freq(inrast)
  ### Remove the NA pixels from the table.
  if (is.na(classpixels[dim(classpixels)[1], 1]) == TRUE) {
    classpixels <- classpixels[-dim(classpixels)[1], ]
  }

  classpixels <- classpixels[order(nrow(classpixels):1), ]
  totpixperclass <- cumsum(classpixels[, 2])
  percentpixels <- totpixperclass / sum(classpixels[, 2])

  classpixels <- cbind(classpixels, totpixperclass, percentpixels)
  classpixels <- classpixels[order(nrow(classpixels):1), ]
  return(classpixels)
}


#' Helper function to compute the area (number of pixels) that a certain threshold has.
#' @param classpixels Pixel threshold class .
#' @param occurinclass Ocurrence points that lies in a certain class.


generatexytableb <- function(classpixels, occurinclass) {
  xytable <- data.frame(classpixels[, c(1, 4)])
  names(xytable) <- c("ClassID", "PercentPixels")
  xytable <- suppressMessages(dplyr::full_join(xytable, occurinclass))

  xytable$Percent[1] <- 1
  xytable$Percent[is.na(xytable$Percent)] <- occurinclass[1, "Percent"]
  xytable <- xytable[, c("ClassID","PercentPixels", "Percent")]
  xytable <- rbind(xytable, c(nrow(xytable) + 1, 0, 0))
  names(xytable) <- c("ClassID", "XCoor", "YCoor")
  return(xytable)
}

#' Helper function to compute AUC (partialAUC, AUC at Random, AUC ratio) values
#' @param xytable A table with the output of the function generatexytableb
#' @param omissionval Omission value.
#' @param iterationno Number of boostrap interation.

calculateauc <- function(xytable, omissionval, iterationno) {
  ## if omissionval is 0, then calculate the complete area under the curve. Otherwise calculate only partial area
  if (omissionval > 0) {
    partialxytable <- xytable[which(xytable[, 3] >= omissionval), ]
    ### Here calculate the X, Y coordinate for the parallel line to x-axis depending upon the omissionval
    ### Get the classid which is bigger than the last row of the xytable and get the xcor and ycor for that class
    ### So that slope of the line is calculated and then intersection point between line parallel to x-axis and passing through
    ### ommissionval on Y-axis is calculated.
    prevxcor <- xytable[which(xytable[, 1] == partialxytable[nrow(partialxytable), 1]) + 1, 2]
    prevycor <- xytable[which(xytable[, 1] == partialxytable[nrow(partialxytable), 1]) + 1, 3]
    xcor1 <- partialxytable[nrow(partialxytable), 2]
    ycor1 <- partialxytable[nrow(partialxytable), 3]
    ## Calculate the point of intersection of line parallel to x-asiz and this line. Use the equation of line
    ## in point-slope form y1 = m(x1-x2)+y2
    slope <- (ycor1 - prevycor) / (xcor1 - prevxcor)
    ycor0 <- omissionval
    xcor0 <- (ycor0 - prevycor + (slope * prevxcor)) / slope
    ### Add this coordinate in the partialxytable with classid greater than highest class id in the table.
    ### Actually class-id is not that important now, only the place where we add this xcor0 and ycor0 is important.
    ### add this as last row in the table
    partialxytable <- rbind(partialxytable,
                            c(partialxytable[nrow(partialxytable), 1] + 1,
                              xcor0, ycor0))
  }else {
    partialxytable <- xytable
  } ### if omissionval > 0

  ## Now calculate the area under the curve on this table.
  xcor1 <- partialxytable[nrow(partialxytable), 2]
  ycor1 <- partialxytable[nrow(partialxytable), 3]
  aucvalue <- 0
  aucvalueatrandom <- 0
  for (i in (nrow(partialxytable) - 1):1) {
    xcor2 <- partialxytable[i,2]
    ycor2 <- partialxytable[i,3]

    # This is calculating the AUCArea for 2 point trapezoid.
    traparea <- (ycor1 * (abs(xcor2 - xcor1))) +
      (abs(ycor2 - ycor1) * abs(xcor2 - xcor1)) / 2
    aucvalue <- aucvalue + traparea
    # now caluclate the area below 0.5 line.
    # Find the slope of line which goes to the point
    # Equation of line parallel to Y-axis is X=k and equation of line at 0.5 is y = x
    trapareaatrandom <- (xcor1 * (abs(xcor2 - xcor1))) +
      (abs(xcor2 - xcor1) * abs(xcor2 - xcor1)) / 2
    aucvalueatrandom <- aucvalueatrandom + trapareaatrandom
    xcor1 <- xcor2
    ycor1 <- ycor2

  }

  newrow <- c(iterationno, aucvalue,
              aucvalueatrandom,
              aucvalue / aucvalueatrandom)

  return(newrow)

}
