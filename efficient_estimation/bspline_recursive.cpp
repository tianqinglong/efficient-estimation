#include <Rcpp.h>
using namespace Rcpp;

// Use dynamic programming to evaluate the b-spline functions at y\in[0, 1)
// Caution: this function is ill-defined at y=1.
// Do NOT evaluate at y=1!

// [[Rcpp::export]]
NumericVector BSplinesForNewY(NumericVector y, NumericVector q,
                              NumericVector bn, NumericVector delta) {
  int i, j;
  double tLower, tUpper, z = y(0), w1, w2;
  NumericMatrix dict(bn(0)+2*q(0)-1, q(0));
  
  for(i=0; i<bn(0)+2*q(0)-1; i++)
  {
    tLower = delta(i);
    tUpper = delta(i+1);
    
    dict(i,0) = (z >= tLower && z < tUpper) ? 1 : 0;
  }
  
  for(i=1; i<q(0); i++)
    for(j=0; j<bn(0)+2*q(0)-1-i; j++)
    {
      w1 = (delta(j+i) == delta(j)) ? 0 : (z-delta(j))/(delta(j+i)-delta(j));
      w2 = (delta(j+i+1) == delta(j+1)) ? 0 : (delta(j+i+1)-z)/(delta(j+i+1)-delta(j+1));
      dict(j,i) = w1*dict(j,i-1)+w2*dict(j+1, i-1);
    }
  
  NumericMatrix output = dict(Range(0, bn(0)+q(0)-1), Range(q(0)-1, q(0)-1));
  return output;
}
