# library(tidyverse)
# library(lbfgs)
# source("simulate_data.R")
# source("bspline.R")
# source("estep.R")
# source("mstep.R")
# Rcpp::sourceCpp("mstep_cpp.cpp")
# Rcpp::sourceCpp("estep_cpp.cpp")
# #-------------------------------
# # Test simulating the data
# #-------------------------------

# # dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))

# #-------------------------------
# # Test adding b-splines
# #-------------------------------

# num_of_sieve <- 4

# dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))
# datWBspline <- AddBsplineColumn(dat, splinesForY, num_of_sieve, 2)

# #-------------------------------
# # Test the conditional expectation function
# #-------------------------------

# datWBspline %>% filter(Obs == 1) -> datNonMissing
# datWBspline %>% filter(Obs == 0) -> datMissing

# yVec <- datNonMissing[,1]
# mat_Spline <- datNonMissing[,-c(1, 2, 3, 4)]
# mat_non_missing_X <- cbind(Intercept=1, datNonMissing[, c(2, 3)])
# mat_X <- cbind(Intercept=1, datMissing[, c(2, 3)])

# betaVec <- c(0.9, -0.5, 0.9)
# tauVec <- rnorm(num_of_sieve)
# sigmaSCL <- sd(yVec)

# mat_Spline <- as.matrix(mat_Spline)
# mat_X <- as.matrix(mat_X)

# conditionalExpectionVec_OnlyY(mat_Spline, mat_X, yVec, tauVec, betaVec, sigmaSCL) -> condProb
# conditionalExpectionVec_OnlyY_Cpp(mat_Spline, mat_X, yVec, tauVec, betaVec, sigmaSCL) -> condProb_cpp

# microbenchmark::microbenchmark(conditionalExpectionVec_OnlyY(mat_Spline, mat_X, yVec, tauVec, betaVec, sigmaSCL),
#                                conditionalExpectionVec_OnlyY_Cpp(mat_Spline, mat_X, yVec, tauVec, betaVec, sigmaSCL))
# #-------------------------------
# # Test Objectives Functions
# #-------------------------------

# newPars <- c(0.8, -0.4, 1, 0.5)
# betaTargetFunc_OnlyY(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

# newPars <- rnorm(num_of_sieve)
# etaTargetFunc_OnlyY(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)
# etaTargetFunc_OnlyY_gr(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

# # Test the first part of the EM

# newPars <- c(0.8, -0.4, 1, 0.5)
# optim(newPars, betaTargetFunc_OnlyY, betaVec = betaVec,
#       sigmaVec = sigmaSCL, tauVec = tauVec, X_obs_mat = mat_non_missing_X,
#       mat_spline = mat_Spline, X_missing_mat = mat_X, Y_non_missing = yVec) -> op1

# betaTargetFuncMax_LM(betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec) -> glmfit

# # Test the second part of the EM

# newPars <- rnorm(num_of_sieve)
# optim(newPars, etaTargetFunc_OnlyY, betaVec = betaVec,
#       sigmaVec = sigmaSCL, tauVec = tauVec, X_obs_mat = mat_non_missing_X,
#       mat_spline = mat_Spline, X_missing_mat = mat_X, Y_non_missing = yVec) -> op2

# lbfgs(etaTargetFunc_OnlyY, etaTargetFunc_OnlyY_gr,
#       betaVec = betaVec, sigmaVec = sigmaSCL, tauVec = tauVec, X_obs_mat = mat_non_missing_X,
#       mat_spline = mat_Spline, X_missing_mat = mat_X, Y_non_missing = yVec,
#       vars = rnorm(num_of_sieve)) -> op_lbfgs

# # Try to solve the inconsistency
# # 
# # b0 <- 2
# # b1 <- 3
# # sd <- 0.7
# # 
# # nn <- 12
# # 
# # xx <- rnorm(nn)
# # yy <- b0+b1*xx+sd*rnorm(nn)
# # 
# # ws <- abs(rnorm(nn))
# # 
# # lmFit <- lm(yy~xx, weights = ws)
# # glmFit <- glm(yy~xx, family = gaussian(link = "identity"), weights = ws)

# # res <- lmFit$residuals
# # sum(res^2*ws)/sum(ws) %>% sqrt

# # sigma(lmFit)

# # Weighted log-likelihood

# # objFunc1 <- function(pars, xx, yy, ws)
# # {
# #    b0 <- pars[1]
# #    b1 <- pars[2]
# #    sd <- pars[3]
# #    
# #    muVec <- b0+b1*xx
# #    loglikVec <- ws*dnorm(yy, muVec, sd, log = T)
# #    
# #    return(-sum(loglikVec))
# # }
# # 
# # pars <- c(1.5, 2.5, 0.05)
# # optim(pars, objFunc1, xx = xx, yy = yy, ws = ws) -> op0
# # 
# # sqrt(sum(ws*lmFit$residuals^2)/sum(ws))

# library(tidyverse)
# library(lbfgs)
# source("simulate_data.R")
# source("bspline.R")
# source("estep.R")
# source("mstep.R")
# Rcpp::sourceCpp("mstep_cpp.cpp")
# Rcpp::sourceCpp("estep_cpp.cpp")

# num_of_sieve <- 4

# dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))
# datWBspline <- AddBsplineColumn(dat, splinesForY, num_of_sieve, 2)

# #-------------------------------
# # Test the conditional expectation function
# #-------------------------------

# datWBspline %>% filter(Obs == 1) -> datNonMissing
# datWBspline %>% filter(Obs == 0) -> datMissing

# yVec <- datNonMissing[,1]
# mat_Spline <- datNonMissing[,-c(1, 2, 3, 4)]
# mat_non_missing_X <- cbind(Intercept=1, datNonMissing[, c(2, 3)])
# mat_X <- cbind(Intercept=1, datMissing[, c(2, 3)])

# betaVec <- c(0.9, -0.5, 0.9)
# tauVec <- rnorm(num_of_sieve)
# sigmaSCL <- sd(yVec)

# conditionalExpectionVec_OnlyY_Cpp(as.matrix(mat_Spline),
#                                   as.matrix(mat_X),
#                                   yVec, tauVec, betaVec, sigmaSCL) -> condProb

# newPars <- rnorm(num_of_sieve)
# etaTargetFunc_OnlyY(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)
# etaTargetFunc_OnlyY_gr(newPars, betaVec, sigmaSCL, tauVec, mat_non_missing_X, mat_Spline, mat_X, yVec)

# paramPackage <- list(nMiss = nrow(datMissing),
#                      nObs = nrow(datNonMissing),
#                      nSpline = ncol(mat_Spline),
#                      splineMatrix = c(t(mat_Spline)),
#                      condProbMat = c(t(condProb)))
# mStepSecondPart(paramPackage)

# lbfgs(etaTargetFunc_OnlyY, etaTargetFunc_OnlyY_gr,
#       betaVec = betaVec, sigmaVec = sigmaSCL, tauVec = tauVec, X_obs_mat = mat_non_missing_X,
#       mat_spline = mat_Spline, X_missing_mat = mat_X, Y_non_missing = yVec,
#       vars = rnorm(num_of_sieve)) -> op_lbfgs
# op_lbfgs

# ### Test the EM algorithm

# library(tidyverse)

# source("simulate_data.R")
# source("bspline.R")
# source("estep.R")
# source("mstep.R")
# source("iterations.R")
# Rcpp::sourceCpp("mstep_cpp.cpp")
# Rcpp::sourceCpp("estep_cpp.cpp")

# num_of_sieve <- 4

# dat <- SimulateData(30, linear.model.additive, c(3, -2, 4), missing.single.y, c(1, -1, 1, 0.5))
# datWBspline <- AddBsplineColumn(dat, splinesForY, num_of_sieve, 2)

#       datNonMissing <- datWBspline[datSpline$Obs == 1,]
#       datMissing <- datWBspline[datSpline$Obs == 0,]

#       # Extract some information

#       ## Non-missing y values
#       yVec <- datNonMissing[,1]

#       ## B-spline matrix
#       obsIndex <- which(colnames(datNonMissing)=="Obs")
#       endIndex <- ncol(datNonMissing)
#       mat_Spline <- datNonMissing[,(obsIndex+1):endIndex]

#       ## Data matrices
#       mat_non_missing_X <- cbind(Intercept=1, datNonMissing[, c(2, obsIndex)])
#       mat_X <- cbind(Intercept=1, datMissing[, c(2, obsIndex-1)])

#       # Initial values

#       ## beta
#       initLM <- lm(Y~., data = datNonMissing[,1:(obsIndex-1)])
#       betaVec <- initLM$coefficients

# source("iterations.R")
# iterationEM(datWBspline)

#-------------------------------
# Rewrite the data simulating functions
#-------------------------------

source("effest/simulate_data.R")
source("effest/bspline.R")
source("effest/estep.R")
source("effest/mstep.R")
source("effest/iterations.R")
Rcpp::sourceCpp("effest/mstep_cpp.cpp")
Rcpp::sourceCpp("effest/estep_cpp.cpp")
Rcpp::sourceCpp("effest/lmloglik.cpp")

# missParam <- c(-1, 1, -0.5, 3) # parameters for function yTransFunc1()
missParam <- c(-1, 1, 0.8) # parameters for function yTransFunc2()


n <- 500
xx <- rnorm(n, 0, 0.5)
xx <- as.matrix(xx, ncol = 1)
xMat <- cbind(1, xx)

betaVal <- c(3, -3.5)

yy <- ySimulatorLM(xMat, betaVal, sd=1)
obsVec <- numeric(n)

FUN <- yTransFunc3

for (i in 1:length(obsVec))
{
      tempProb <- PiY(yy[i], FUN, missParam)
      obsVec[i] <- ifelse(runif(1) < tempProb, 1, 0)
}

num_of_sieve <- 3
dat <- data.frame(Y = yy, Z = xx, Obs = obsVec)
datSpline <- AddBsplineColumn(dat, splinesForY, num_of_sieve, 2)
table(obsVec)
#-------------------
# datNonMissing <- datSpline[datSpline$Obs == 1,]
# datMissing <- datSpline[datSpline$Obs == 0,]

# # Extract some information

# ## Non-missing y values
# yVec <- datNonMissing[,1]

# ## B-spline matrix
# obsIndex <- which(colnames(datNonMissing)=="Obs")
# endIndex <- ncol(datNonMissing)
# mat_Spline <- datNonMissing[,(obsIndex+1):endIndex]

# ## Data matrices
# mat_non_missing_X <- cbind(Intercept=1, datNonMissing[, 2:(obsIndex-1)])
# mat_X <- cbind(Intercept=1, datMissing[, 2:(obsIndex-1)])
#-------------------

t0 <- Sys.time()
iterationEM(datSpline)

lm(yy~0+xMat)

dat_no_missing <- dat[dat$Obs == 1,]

lm(Y~Z, data = dat_no_missing) -> fittemp
summary(fittemp)

iterationEM_2(datSpline)
Sys.time()-t0

table(obsVec)
