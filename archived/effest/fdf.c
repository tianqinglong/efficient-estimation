double inv_logit(double x)
{
	double out;
	out = exp(x)/(1+exp(x));

	return out;
}

double
obj_f (const gsl_vector *v, void *extraParams)
{
	struct params *paramList = (struct params *) extraParams;

	int nMiss =  paramList->nMiss;
	int nObs = paramList->nObs;
	int nSpline = paramList->nSpline;
	double* splineMat = paramList->splineMat; // length = nObs*nSpline; by row
	double* condProb = paramList->condProb; // length = nMiss*nObs; by row

	double loglik = 0;
	double tempSum, tempWeight, templik, negloglik;

	int i, j, k;

	for (i=0; i<nObs; i++)
	{
		tempSum = 0;
		for (j=0; j<nSpline; j++)
		{
			tempSum += splineMat[i*nSpline+j]*gsl_vector_get(v, j);
		}
		loglik += log(inv_logit(tempSum));
	}

	for (i=0; i<nMiss; i++)
	{
		for (j=0; j<nObs; j++)
		{
			tempWeight = condProb[i*nObs+j];
			tempSum = 0;
			for (k=0; k<nSpline; k++)
			{
				tempSum += splineMat[j*nSpline+k]*gsl_vector_get(v, k);
			}
			templik = log(1-inv_logit(tempSum));
			loglik += tempWeight*templik;
		}
	}

	negloglik = -loglik;

	return negloglik;
}

void
gr_f (const gsl_vector *v, void *extraParams, gsl_vector *df)
{
	struct params *paramList = (struct params *) extraParams;

	int nMiss =  paramList->nMiss;
	int nObs = paramList->nObs;
	int nSpline = paramList->nSpline;
	double* splineMat = paramList->splineMat; // length = nObs*nSpline; by row
	double* condProb = paramList->condProb; // length = nMiss*nObs; by row

	int i, j, k, m;
	double grTemp, tempSum, tempWeight;

	for (i=0; i<nSpline; i++)
	{
		grTemp = 0;
		for (j=0; j<nObs; j++)
		{
			tempSum=0;
			for (k=0; k<nSpline;k++)
			{
				tempSum += splineMat[j*nSpline+k]*gsl_vector_get(v, k);
			}
			grTemp += splineMat[j*nSpline+i]/(1+exp(tempSum));
		}

		for (j=0; j<nMiss; j++)
		{
			for (k=0; k<nObs; k++)
			{
				tempSum=0;
				for (m=0; m<nSpline; m++)
				{
					tempSum += splineMat[k*nSpline+m]*gsl_vector_get(v, m);
				}
				tempWeight = condProb[j*nObs+k];
				grTemp -= tempWeight*inv_logit(tempSum)*splineMat[k*nSpline+i];
			}
		}
		gsl_vector_set(df, i, -grTemp);
	}
}

void
f_df (const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
	*f = obj_f(v, params);
	gr_f(v, params, df);
}
