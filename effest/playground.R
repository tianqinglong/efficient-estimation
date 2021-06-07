library(tidyverse)
library(lbfgs)
#-------------------------------
# Test the functions
#-------------------------------

SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))

#-------------------------------
# Test
#-------------------------------

num_of_sieve <- 4

dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))
datWBspline <- AddBsplineColumn(dat, splinesForY, num_of_sieve, 2)

#-------------------------------
# Test
#-------------------------------

datWBspline %>% filter(Obs == 1) -> datNonMissing
datWBspline %>% filter(Obs == 0) -> datMissing

yVec <- datNonMissing[,1]
mat_Spline <- datNonMissing[,-c(1, 2, 3, 4)]
mat_non_missing_X <- cbind(Intercept=1, datNonMissing[, c(2, 3)])
mat_X <- cbind(Intercept=1, datMissing[, c(2, 3)])

betaVec <- c(0.9, -0.5, 0.9)
tauVec <- rnorm(num_of_sieve)
sigmaSCL <- sd(yVec)

conditionalExpectionVec_OnlyY(mat_Spline, mat_X, yVec, tauVec, betaVec, sigmaSCL) -> condProb

#-------------------------------
# Test Target Functions
#-------------------------------

newPars <- c(0.8, -0.4, 1, 0.5)
betaTargetFunc_OnlyY(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

newPars <- rnorm(num_of_sieve)
etaTargetFunc_OnlyY(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)
etaTargetFunc_OnlyY_gr(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

# Test the first part of the EM

newPars <- c(0.8, -0.4, 1, 0.5)
optim(newPars, betaTargetFunc_OnlyY, betaVec = betaVec,
      sigmaVec = sigmaSCL, tauVec = tauVec, X_obs_mat = mat_non_missing_X,
      mat_spline = mat_Spline, X_missing_mat = mat_X, Y_non_missing = yVec) -> op1

betaTargetFuncMax_LM(betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec) -> glmfit

# Test the second part of the EM

newPars <- rnorm(num_of_sieve)
optim(newPars, etaTargetFunc_OnlyY, betaVec = betaVec,
      sigmaVec = sigmaSCL, tauVec = tauVec, X_obs_mat = mat_non_missing_X,
      mat_spline = mat_Spline, X_missing_mat = mat_X, Y_non_missing = yVec) -> op2

lbfgs(etaTargetFunc_OnlyY, etaTargetFunc_OnlyY_gr,
      betaVec = betaVec, sigmaVec = sigmaSCL, tauVec = tauVec, X_obs_mat = mat_non_missing_X,
      mat_spline = mat_Spline, X_missing_mat = mat_X, Y_non_missing = yVec,
      vars = rnorm(num_of_sieve))
