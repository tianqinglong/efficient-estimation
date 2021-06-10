struct params
{
	int nMiss;
	int nObs;
	int nSpline;
	double *splineMat;
	double *condProb;
};

double inv_logit(double x);
double obj_f(const gsl_vector *v, void *extraParams);
void gr_f(const gsl_vector *v, void *extraParams, gsl_vector *df);
void f_df(const gsl_vector *v, void *params, double *f, gsl_vector *df);