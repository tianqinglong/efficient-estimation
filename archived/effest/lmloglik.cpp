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

NumericVector betaLogLik(NumericVector betaNew,
                         NumericVector sigmaNew,
                         NumericMatrix mat_non_missing_X,
                         NumericMatrix mat_X,
                         NumericVector yVec,
                         NumericMatrix condProb)
{
  int i, j;
  int ncoef = sigmaNew.length(), n_non_missing = mat_non_missing_X.nrow(), n_missing = mat_X.nrow();
  double logLik1 = 0, mutemp, weighttemp;
  
  for (i=0; i<n_non_missing; i++)
  {
    mutemp = 0;
    for (j=0; j<ncoef; j++)
    {
      mutemp += mat_non_missing_X(i,j)*betaNew(j);
    }
    logLik1 += R::dnorm(yVec(i), mutemp, sigmaNew(0), 1);
  }
  
  for (i=0; i<n_missing; i++)
  {
    mutemp = 0;
    for (j=0; j<ncoef; j++)
    {
      mutemp += mat_X(i,j)*betaNew(j);
    }
    
    for (j=0; j<n_non_missing; j++)
    {
      weighttemp = condProb(i,j);
      logLik1 += R::dnorm(yVec(j), mutemp, sigmaNew(0), 1)*weighttemp;
    }
  }
  
  NumericVector negloglik(1);
  negloglik(0) = -logLik1;
  
  return negloglik;
}
