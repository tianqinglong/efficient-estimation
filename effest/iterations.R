#-------------------------------
# This function is to do the EM iterations
# Author: Qinglong Tian
# Date: June 9, 2021
#-------------------------------

iterationEM <- function(datSpline, max_iter = 500, tol = 1e-4)
{
	options(digits=4)
	# Split the data into missing and non-missing parts
	datNonMissing <- datSpline[datSpline$Obs == 1,]
	datMissing <- datSpline[datSpline$Obs == 0,]

	# Extract some information

	## Non-missing y values
	yVec <- datNonMissing[,1]

	## B-spline matrix
	obsIndex <- which(colnames(datNonMissing)=="Obs")
	endIndex <- ncol(datNonMissing)
	mat_Spline <- datNonMissing[,(obsIndex+1):endIndex]

	## Data matrices
	mat_non_missing_X <- cbind(Intercept=1, datNonMissing[, c(2, obsIndex)])
	mat_X <- cbind(Intercept=1, datMissing[, c(2, obsIndex-1)])

	# Initial values

	## beta
	initLM <- lm(Y~., data = datNonMissing[,1:(obsIndex-1)])
	betaVec <- initLM$coefficients

	## sigma 
	sigmaSCL <- sigma(initLM)

	## tau
	nSpline <- endIndex-obsIndex
	tauVec <- numeric(nSpline)

	# EM algorithm
	iter <- 0
	success <- 0
	while (iter < max_iter & !success)
	{
		iter <- iter+1

		## E-step

		condProb <- conditionalExpectionVec_OnlyY_Cpp(as.matrix(mat_Spline), as.matrix(mat_X), yVec, tauVec, betaVec, sigmaSCL)

		## M-step

		### Tau part
		
		paramList <- list(nMiss = nrow(datMissing), nObs = nrow(datNonMissing), nSpline = ncol(mat_Spline),
			splineMatrix = c(t(mat_Spline)), condProbMat = c(t(condProb)))
		tauNew <- mStepSecondPart(paramList)

		### beta/sigma part
		betaSigmaNew <- betaTargetFuncMax_LM_light(condProb, mat_non_missing_X, mat_Spline, mat_X, yVec)

		betaNew <- betaSigmaNew[1:length(betaVec)]
		sigmaNew <- betaSigmaNew[length(betaSigmaNew)]

		## Check convergence
		diff <- (sigmaNew-sigmaSCL)^2 + sum((betaNew-betaVec)^2) + sum((tauNew-tauVec)^2)

		print(paste("Iteration ", iter,"; Diff=", diff, "; Sigma=", format(round(sigmaNew, 4), nsmall = 4), "; Beta=", format(round(max(betaNew), 4), nsmall = 4), "; Tau=", format(round(max(tauNew), 2), nsmall = 2), sep = ""))

		if (diff < tol)
		{
			success <- 1
		}

		betaVec <- betaNew
		tauVec <- tauNew
		sigmaSCL <- sigmaNew
	}


	if (success)
	{
		print("Success!")
		return(list(beta = betaVec, tau = tauVec, sigma = sigmaSCL))
	}
	else
	{
		warning("Failed to converge!")
	}
}