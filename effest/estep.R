inv_logit <- function(x)
{
  return(exp(x)/(1+exp(x)))
}
#-------------------------------
# Function to compute the conditional expection of the indicator function (E-step)
# Author: Qinglong Tian
# Date: 5 June, 2021
#-------------------------------

conditionalExpectionVec_OnlyY <- function(mat_spline, mat_x, vec_y, tauVec, betaVec, sigmaScalar)
# mat_spline: is the matrix for the b-spline (Y Only, no X involved), columns is the spline, row is the observation
# mat_x: is the model matrix (only missing), the first column should be all 1;
# vec_y: is the non-missing response; 
# tauVec: is the estimate of tau from the last iteration
# betaVec: is the estiamte of beta from the last iteration
# sigmaScalar: is the estimate of sigma from the last iteration
{
  n_non_missing <- nrow(mat_spline)
  n_missing <- nrow(mat_x)
  
  if (n_missing < 1 || n_non_missing <1)
  {
    stop("Data input not valid!")
  }
  
  # Each row is for "k", each column is for the missing obs "i"
  unscaled_condtional_expectation <- matrix(nrow = n_non_missing, ncol = n_missing)
  
  # Check dimensions
  if (ncol(mat_x) != length(betaVec))
  {
    stop("Model dimension is wrong!")
  }
  
  if (ncol(mat_spline) != length(tauVec))
  {
    stop("Dimension of the spline does not match!")
  }
  
  # Fill the matrix
  for (k in 1:n_non_missing)
    for (i in 1:n_missing)
    {
      cur_spline <- mat_spline[k,]
      yk <- vec_y[k]
      
      pi_k <- inv_logit(sum(cur_spline*tauVec))
      
      cur_xVec <- mat_x[i,]
      mu_temp <- sum(cur_xVec*betaVec)
      sigma_temp <- sigmaScalar
      likModel <- dnorm(yk, mean = mu_temp, sd = sigma_temp)
      
      unscaled_condtional_expectation[k, i] <- (1-pi_k)*likModel
    }
  
  unscaled_condtional_expectation <- t(unscaled_condtional_expectation)
  sumAcrossK <- rowSums(unscaled_condtional_expectation)
  
  for (i in 1:n_missing)
  {
    unscaled_condtional_expectation[i,] <- unscaled_condtional_expectation[i,]/sumAcrossK[i]
  }
  
  return(unscaled_condtional_expectation)
}

#-------------------------------
# The target functions (that need to be maximized)
# Author: Qinglong Tian
# Date: 5 June, 2021
#-------------------------------

betaTargetFunc_OnlyY <- function(pars, betaVec, sigmaVec, tauVec, X_obs_mat,
                                 mat_spline, X_missing_mat, Y_non_missing)
# Compute the target function of beta (as well as sigma)
{
  betaNew <- pars[-length(pars)]
  sigmaNew <- pars[length(pars)]
  
  n_non_missing <- nrow(X_obs_mat)
  n_missing <- nrow(X_missing_mat)
  
  logLik1 <- 0
  
  for (i in 1:n_non_missing)
  {
    mu_temp <- sum(X_obs_mat[i,]*betaNew)
    logLik1 <- logLik1 + dnorm(Y_non_missing[i], mu_temp, sigmaNew, log = TRUE)
  }
  
  conditionalExpectionVec_OnlyY(mat_spline, X_missing_mat, Y_non_missing, tauVec,
                                betaVec, sigmaVec) -> condProbMat
  for (i in 1:n_missing)
  {
    mu_temp <- sum(X_missing_mat[i,]*betaNew)

    logLikModel <- dnorm(Y_non_missing, mu_temp, sigmaNew, log = TRUE)
    weightVec <- condProbMat[i,]
    
    weightedSum <- sum(weightVec*logLikModel)
    
    logLik1 <- logLik1+weightedSum
  }
  
  return(logLik1)
}

etaTargetFunc_OnlyY <- function(pars, betaVec, sigmaVec, tauVec, X_obs_mat,
                                mat_spline, X_missing_mat, Y_non_missing)
{
  tauNew <- pars
  
  n_non_missing <- nrow(X_obs_mat)
  n_missing <- nrow(X_missing_mat)
  
  logLik1 <- 0
  
  for (i in 1:n_non_missing)
  {
    splineValX <- mat_spline[i,]
    logLik1 <- log(inv_logit(sum(splineValX*tauNew)))+logLik1
  }
  
  conditionalExpectionVec_OnlyY(mat_spline, X_missing_mat, Y_non_missing, tauVec,
                                betaVec, sigmaVec) -> condProbMat
  
  for (i in 1:n_missing)
  {
    for (k in 1:n_non_missing)
    {
      weight <- condProbMat[i, k]
      splineValX <- mat_spline[k,]
      logLik_temp <- log(1-inv_logit(sum(splineValX*tauNew)))
      logLik1 <- logLik1 + logLik_temp*weight
    }
  }
  
  return(logLik1)
}
