#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

NumericMatrix conditionalExpectionVec_OnlyY_Cpp(NumericMatrix mat_spline,
                                                NumericMatrix mat_x,
                                                NumericVector vec_y,
                                                NumericVector tauVec,
                                                NumericVector betaVec,
                                                NumericVector sigmaScalar)
{
  int i, j, k;
  int n_non_missing = mat_spline.nrow(), n_missing = mat_x.nrow();
  int n_spline = mat_spline.ncol(), n_coef = mat_x.ncol();
  double tempSum, pi_i;
  
  if (n_missing < 1 || n_non_missing < 1)
    return 1;
  
  NumericMatrix unscaled_condtional_expectation(n_non_missing, n_missing);
  
  for (i=0; i<n_non_missing; i++)
  {
    for (j=0; j<n_missing; j++)
    {
      tempSum = 0;
      for (k=0; k<n_spline; k++)
      {
        tempSum += mat_spline(i,k)*tauVec(k);
      }
      pi_i = std::exp(tempSum)/(1+std::exp(tempSum));
      
      tempSum = 0;
      for (k=0; k<n_coef; k++)
      {
        tempSum += mat_x(j, k)*betaVec(k);
      }
  
      unscaled_condtional_expectation(i,j) = (1-pi_i)*R::dnorm(vec_y(i), tempSum, sigmaScalar(0), 0);
    }
  }
  
  NumericMatrix newMat = transpose(unscaled_condtional_expectation);
  NumericVector rowsums1 = rowSums(newMat);
  
  for (i=0; i<n_missing;i++)
  {
    for (j=0; j<n_non_missing; j++)
    {
      newMat(i,j) = newMat(i,j)/rowsums1(i);
    }
  }
  
  return newMat;
}
