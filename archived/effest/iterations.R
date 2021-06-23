#-------------------------------
# This function is to do the EM iterations
# Author: Qinglong Tian
# Date: June 9, 2021
#-------------------------------

# Stopping criteria: log-likelihood does not change much

iterationEM <- function(datSpline, max_iter = 500, tol = 1e-4)
{
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
	mat_non_missing_X <- cbind(Intercept=1, datNonMissing[, 2:(obsIndex-1)])
	mat_X <- cbind(Intercept=1, datMissing[, 2:(obsIndex-1)])

	xName <- colnames(datNonMissing)[2:(obsIndex-1)]
	colnames(mat_non_missing_X)[2:(obsIndex-1)] <- xName
	colnames(mat_X)[[2:(obsIndex-1)]] <- xName

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
		outNew <- mStepSecondPart(paramList)
		tauNew <- outNew[1:(length(outNew)-1)]

		loglik2 <- outNew[length(outNew)]

		### beta/sigma part
		betaSigmaNew <- betaTargetFuncMax_LM_light(condProb, mat_non_missing_X, mat_Spline, mat_X, yVec)

		betaNew <- betaSigmaNew[1:length(betaVec)]
		sigmaNew <- betaSigmaNew[length(betaSigmaNew)]

		# loglik1 <- betaTargetFunc_OnlyY(c(betaNew, sigmaNew), betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

		loglik1 <- betaLogLik(betaNew, sigmaNew, mat_non_missing_X, mat_X, yVec, condProb)

		loglik1_true <- betaLogLik(c(3, -3.5), 1, mat_non_missing_X, mat_X, yVec, condProb)

		## Check convergence
		newlogLik <- loglik1+loglik2
		newlogLik_true_beta <- loglik1_true+loglik2
		if (iter == 1)
		{
			loglikSum <- newlogLik
			diff <- 1
		}
		else
		{

			diff <- abs(loglikSum-newlogLik)/abs(loglikSum)
		}

		print(paste("Iteration ", iter,": Change=", format(round(100*diff, 4), nsmall = 4),
			"%; Sigma=", format(round(sigmaNew, 4), nsmall = 4),
			"; Max_Beta=", format(round(max(betaNew), 4), nsmall = 4),
			"; Max_Tau=", format(round(max(tauNew), 2), nsmall = 2),
			"; Log-lik=", newlogLik,
			"; Log-lik-2=", loglik2, sep = ""))

		if (diff < tol)
		{
			success <- 1
		}

		betaVec <- betaNew
		tauVec <- tauNew
		sigmaSCL <- sigmaNew
		loglikSum <- newlogLik
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

# Stopping criteria: parameters not change much

iterationEM_2 <- function(datSpline, max_iter = 500, tol = 1e-4)
{
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
	mat_non_missing_X <- cbind(Intercept=1, datNonMissing[, 2:(obsIndex-1)])
	mat_X <- cbind(Intercept=1, datMissing[, 2:(obsIndex-1)])

	xName <- colnames(datNonMissing)[2:(obsIndex-1)]
	colnames(mat_non_missing_X)[2:(obsIndex-1)] <- xName
	colnames(mat_X)[[2:(obsIndex-1)]] <- xName

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
		outNew <- mStepSecondPart(paramList)
		tauNew <- outNew[1:(length(outNew)-1)]

		loglik2 <- outNew[length(outNew)]

		### beta/sigma part
		betaSigmaNew <- betaTargetFuncMax_LM_light(condProb, mat_non_missing_X, mat_Spline, mat_X, yVec)

		betaNew <- betaSigmaNew[1:length(betaVec)]
		sigmaNew <- betaSigmaNew[length(betaSigmaNew)]

		# loglik1 <- betaTargetFunc_OnlyY(c(betaNew, sigmaNew), betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

		# loglik1 <- betaLogLik(betaNew, sigmaNew, mat_non_missing_X, mat_X, yVec, condProb)

		## Check convergence
		# newlogLik <- loglik1+loglik2

		max_abs_diff <- max(abs(sigmaNew-sigmaSCL), abs(betaNew-betaVec), abs(tauVec-tauNew))
		diff <- max_abs_diff^2

		print(paste("Iteration ", iter,": Change=", diff, "; Sigma=", format(round(sigmaNew, 4), nsmall = 4), "; Max_Beta=", format(round(max(betaNew), 4), nsmall = 4), "; Max_Tau=", format(round(max(tauNew), 2), nsmall = 2), sep = ""))

		if (diff < tol)
		{
			success <- 1
		}

		betaVec <- betaNew
		tauVec <- tauNew
		sigmaSCL <- sigmaNew
		# loglikSum <- newlogLik
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