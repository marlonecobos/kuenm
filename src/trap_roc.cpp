#include <Rcpp.h>
using namespace Rcpp;


// Trapezoidal integration for Partial ROC curve
// @param x is the fractional area predicted presence
// @param sensibility sensibility at threshold (fract_area value).
// @return Returns the area under the Partial ROC curve.
// [[Rcpp::export]]

double trap_roc(NumericVector x, NumericVector y) {
  int x_s = x.size();
  int y_s = y.size();
  if(x_s != y_s)
    ::Rf_error("x  and y must have the same length x_size (%d) != y_size (%d)", x_s, y_s);
  double auc = 0;

  for(int i = 1; i < x_s; ++i) {
    auc += 0.5*(y[i-1] + y[i])*(x[i]-x[i-1]);
  }
  return auc;
}
