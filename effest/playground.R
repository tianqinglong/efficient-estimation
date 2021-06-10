library(tidyverse)
library(lbfgs)
#-------------------------------
# Test simulating the data
#-------------------------------

# dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))

#-------------------------------
# Test adding b-splines
#-------------------------------

num_of_sieve <- 4

dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))
datWBspline <- AddBsplineColumn(dat, splinesForY, num_of_sieve, 2)

#-------------------------------
# Test the conditional expectation function
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
# Test Objectives Functions
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
      vars = rnorm(num_of_sieve)) -> op_lbfgs

# Try to solve the inconsistency
# 
# b0 <- 2
# b1 <- 3
# sd <- 0.7
# 
# nn <- 12
# 
# xx <- rnorm(nn)
# yy <- b0+b1*xx+sd*rnorm(nn)
# 
# ws <- abs(rnorm(nn))
# 
# lmFit <- lm(yy~xx, weights = ws)
# glmFit <- glm(yy~xx, family = gaussian(link = "identity"), weights = ws)

# res <- lmFit$residuals
# sum(res^2*ws)/sum(ws) %>% sqrt

# sigma(lmFit)

# Weighted log-likelihood

# objFunc1 <- function(pars, xx, yy, ws)
# {
#    b0 <- pars[1]
#    b1 <- pars[2]
#    sd <- pars[3]
#    
#    muVec <- b0+b1*xx
#    loglikVec <- ws*dnorm(yy, muVec, sd, log = T)
#    
#    return(-sum(loglikVec))
# }
# 
# pars <- c(1.5, 2.5, 0.05)
# optim(pars, objFunc1, xx = xx, yy = yy, ws = ws) -> op0
# 
# sqrt(sum(ws*lmFit$residuals^2)/sum(ws))

library(tidyverse)
source("simulate_data.R")
source("bspline.R")
source("estep.R")
source("mstep.R")
Rcpp::sourceCpp("mstep_cpp.cpp")

num_of_sieve <- 4

dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))
datWBspline <- AddBsplineColumn(dat, splinesForY, num_of_sieve, 2)

#-------------------------------
# Test the conditional expectation function
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

newPars <- rnorm(num_of_sieve)
etaTargetFunc_OnlyY(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)
etaTargetFunc_OnlyY_gr(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

paramPackage <- list(nMiss = nrow(datMissing),
                     nObs = nrow(datNonMissing),
                     nSpline = ncol(mat_Spline),
                     splineMatrix = c(t(mat_Spline)),
                     condProbMat = c(t(condProb)))
mStepSecondPart(paramPackage)
