#' Calculation of partial areas under the ROC curve
#'
#' @description calc_auc calculates the area under the curve using a threshold (E) value.
#' This function is designed to be used specifically in the \code{\link{kuenm_proc}} function.
#'
#' @param xytable a XY coordinates table prepared with the \code{\link{g_xy_tab}} function.
#' @param omissionval (numeric) value from 0 to 100 that will be used as threshold (E).
#' @param iterationno (numeric) number of bootstrap iterations to be performed.

calc_auc <- function(xytable, omissionval, iterationno) {
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
    ## in point-slope form y1 <- m(x1-x2)+y2
    slope <- (ycor1 - prevycor) / (xcor1 - prevxcor)
    ycor0 <- omissionval
    xcor0 <- (ycor0 - prevycor + (slope * prevxcor)) / slope
    ### Add this coordinate in the partialxytable with classid greater than highest class id in the table.
    ### Actually class-id is not that important now, only the place where we add this xcor0 and ycor0 is important.
    ### add this as last row in the table
    partialxytable <- rbind(partialxytable, c(partialxytable[nrow(partialxytable), 1] + 1, xcor0, ycor0))
  }else {
    partialxytable <- xytable
  } ### if omissionval > 0

  ## Now calculate the area under the curve on this table.
  xcor1 <- partialxytable[nrow(partialxytable), 2]
  ycor1 <- partialxytable[nrow(partialxytable), 3]
  aucvalue <- 0
  aucvalueatrandom <- 0

  for (i in (nrow(partialxytable) - 1):1) {
    xcor2 <- partialxytable[i, 2]
    ycor2 <- partialxytable[i, 3]

    # This is calculating the AUCArea for 2 point trapezoid.
    traparea <- (ycor1 * (abs(xcor2 - xcor1))) + (abs(ycor2 - ycor1) * abs(xcor2 - xcor1)) / 2
    aucvalue <- aucvalue + traparea
    # now caluclate the area below 0.5 line.
    # Find the slope of line which goes to the point
    # Equation of line parallel to Y-axis is X=k and equation of line at 0.5 is y = x
    trapareaatrandom <- (xcor1 * (abs(xcor2 - xcor1))) + (abs(xcor2 - xcor1) * abs(xcor2 - xcor1)) / 2
    aucvalueatrandom <- aucvalueatrandom + trapareaatrandom
    xcor1 <- xcor2
    ycor1 <- ycor2
  }
  newrow <- c(iterationno, aucvalue, aucvalueatrandom, aucvalue / aucvalueatrandom)

  return(newrow)
}
