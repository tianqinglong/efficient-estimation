#-------------------------------
# The objective functions (that need to be maximized)
# The first part still needs investigation but the second part is done.
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
  
  return(-logLik1)
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
  
  return(-logLik1)
}

#-------------------------------
# Use existing algorithms to derive the beta values
# Comments: beta estimates match the optim function, but sigma does not. Need to investigate
#-------------------------------

# betaTargetFuncMax_LM <- function(betaVec, sigmaVec, tauVec, X_obs_mat,
#                                  mat_spline, X_missing_mat, Y_non_missing)
# {
#   n_non_missing <- nrow(X_obs_mat)
#   n_missing <- nrow(X_missing_mat)
  
#   Y_new <- c(Y_non_missing, rep(Y_non_missing, n_missing))
#   X_new <- as.matrix(X_obs_mat)
#   for (i in 1:n_missing)
#   {
#     X_inst <- X_missing_mat[i,]
#     X_add <- matrix(rep(X_inst, n_non_missing), ncol = length(X_inst), byrow = TRUE)
#     X_new <- rbind(X_new, X_add)
#   }
  
#   X_new <- matrix(unlist(X_new), ncol = 3)
  
#   WeightVec <- rep(1, n_non_missing)
#   conditionalExpectionVec_OnlyY(mat_spline, X_missing_mat, Y_non_missing, tauVec,
#                                 betaVec, sigmaVec) -> condProbMat
#   WeightVec <- c(WeightVec, c(t(condProbMat)))
  
#   dat <- cbind(Y_new, X_new)
#   colnames(dat) <- c("Y", "Intercept", "U", "Z")
#   dat <- as.data.frame(dat)
  
#   glm(Y~U+Z, data = dat, weights = WeightVec, family=gaussian(link = "identity")) -> lmFit
  
#   # betaNew <- lmFit$coefficients
#   # sigmaNew <- sqrt(sum(WeightVec*lmFit$residuals^2)/sum(WeightVec))
  
#   return(lmFit)
# }

# This version does the E step outside (as input)
betaTargetFuncMax_LM_light <- function(condProbMat, X_obs_mat, mat_spline, X_missing_mat, Y_non_missing)
{
  n_non_missing <- nrow(X_obs_mat)
  n_missing <- nrow(X_missing_mat)
  
  Y_new <- c(Y_non_missing, rep(Y_non_missing, n_missing))
  Y_new <- as.matrix(Y_new, ncol = 1)
  X_new <- as.matrix(X_obs_mat)
  for (i in 1:n_missing)
  {
    X_inst <- X_missing_mat[i,]
    X_add <- matrix(rep(X_inst, n_non_missing), ncol = length(X_inst), byrow = TRUE)
    X_new <- rbind(X_new, X_add)
  }
  
  WeightVec <- rep(1, n_non_missing)
  WeightVec <- c(WeightVec, c(t(condProbMat)))
  
  # dat <- cbind(Y_new, X_new)

  # colnames(dat) <- c("Y", colnames(X_obs_mat))
  # dat <- as.data.frame(dat)

  # glm(Y~0+., data = dat, weights = WeightVec, family=gaussian(link = "identity")) -> lmFit

  wts <- sqrt(WeightVec)
  lmFit <- .Call(stats:::C_Cdqrls, X_new*wts, Y_new*wts, 1e-4, FALSE)

  betaNew <- lmFit$coefficients
  sigmaNew <- sqrt(sum(WeightVec*lmFit$residuals^2)/sum(WeightVec))
  
  return(c(betaNew, sigmaNew))
}

#-------------------------------
# The gradient function of the second part of the EM algorithm
#-------------------------------

etaTargetFunc_OnlyY_gr <- function(params, betaVec, sigmaVec, tauVec, X_obs_mat,
                                   mat_spline, X_missing_mat, Y_non_missing)
# params: the tau vector
{
  tauNew <- params
  
  n_non_missing <- nrow(X_obs_mat)
  n_missing <- nrow(X_missing_mat)
  
  sn <- ncol(mat_spline)
  grVec <- numeric(sn)
  
  condProb <- conditionalExpectionVec_OnlyY(mat_spline, X_missing_mat, Y_non_missing,
                                            tauVec, betaVec, sigmaVec)
  
  for (j in 1:sn)
  {
    gr_j <- 0
    for (i in 1:n_non_missing)
    {
      Z_i <- sum(mat_spline[i,]*tauNew)
      gr_j <- gr_j+mat_spline[i,j]/(1+exp(Z_i))
    }
    
    for (i in 1:n_missing)
      for (k in 1:n_non_missing)
      {
        Z_k <- sum(mat_spline[k,]*tauNew)
        w_ki <- condProb[i,k]
        gr_j <- gr_j-w_ki*inv_logit(Z_k)*mat_spline[k,j]
      }
    grVec[j] <- gr_j
  }
  return(-grVec)
}

