library(tidyverse)
#-------------------------------
# Test the functions
#-------------------------------

SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))

#-------------------------------
# Test
#-------------------------------

dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))
datWBspline <- AddBsplineColumn(dat, splinesForY, 6, 2)

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
tauVec <- rnorm(6)
sigmaSCL <- sd(yVec)

conditionalExpectionVec_OnlyY(mat_Spline, mat_X, yVec, tauVec, betaVec, sigmaSCL) -> condProb

#-------------------------------
# Test Target Functions
#-------------------------------

newPars <- c(0.8, -0.4, 1, 0.5)
betaTargetFunc_OnlyY(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

newPars <- c(1, 0.5, 1, 1, 1, 1)
etaTargetFunc_OnlyY(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

newPars <- c(0.8, -0.4, 1, 0.5)
optim(newPars, betaTargetFunc_OnlyY, betaVec = betaVec,
      sigmaVec = sigmaSCL, tauVec = tauVec, X_obs_mat = mat_non_missing_X,
      mat_spline = mat_Spline, X_missing_mat = mat_X, Y_non_missing = yVec) -> op1

betaTargetFuncMax_LM(betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)
