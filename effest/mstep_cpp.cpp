#include <RcppGSL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

extern "C" {
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include "fdf.h"
#include "fdf.c"
}

// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]

Rcpp::NumericVector mStepSecondPart(Rcpp::List pars)
{
  int nMiss = pars[0], nObs = pars[1], nSpline = pars[2];
  
  Rcpp::NumericVector splineMatix = pars[3];
  Rcpp::NumericVector condProbMat = pars[4];
  
  double *splineMat = splineMatix.begin();
  double *condProb = condProbMat.begin();
  
  struct params extraParam = {nMiss, nObs, nSpline, splineMat, condProb};
  
  size_t iter = 0;
  int i, status;
  
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  
  gsl_vector *x;
  gsl_multimin_function_fdf my_func;
  
  my_func.n = nSpline;
  my_func.f = &obj_f;
  my_func.df = &gr_f;
  my_func.fdf = &f_df;
  my_func.params = &extraParam;
  
  x = gsl_vector_alloc(nSpline);
  for(i=0; i<nSpline; i++)
  {
    gsl_vector_set(x, i, 0);
  }
  
  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc (T, nSpline);
  
  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.005, 1e-5);
  
  do
  {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status)
      break;
    
    status = gsl_multimin_test_gradient (s->gradient, 1e-3);
    
  }
  while (status == GSL_CONTINUE && iter < 100);
  
  Rcpp::NumericVector output (nSpline+1);
  for (i = 0; i<nSpline; i++)
  {
    output[i] = gsl_vector_get(s->x, i);
  }
  
  output[nSpline] = obj_f(s->x, &extraParam);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
  
  return output;
}
