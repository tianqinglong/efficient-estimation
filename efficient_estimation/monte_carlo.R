#-----------------------
# Simulate N samples and compute the empirical covariance matrix
# Both covariates, X_1 and X_2, are centered at 0
# Use the intercept to control the proportion of missing
# Missing mechanism is to use the sigmoid function, could add some variation by adding sin(2*pi*y)
#-----------------------
library(fastGHQuad)
library(tidyverse)
source("simulate_data.R")
Rcpp::sourceCpp("bspline_recursive.cpp")
source("add_splines.R")
source("estep.R")
source("mstep.R")
source("variance.R")
#-----------------------
# Auxiliary functions
#-----------------------
monteCarloInDataFrame <- function(monteCarloResults)
{
  numBeta <- length(monteCarloResults[[1]]$Beta_EM)
  
  BetaEM <- t(sapply(monteCarloResults, function(x) {x$Beta_EM}))
  colnames(BetaEM) <- paste("BetaEM_", 1:numBeta, sep = "")
  
  BetaNM <- t(sapply(monteCarloResults, function(x) {x$Beta_NM}))
  colnames(BetaNM) <- paste("BetaNM_", 1:numBeta, sep = "")
  
  BetaOR <- t(sapply(monteCarloResults, function(x) {x$Beta_OR}))
  colnames(BetaOR) <- paste("BetaOR_", 1:numBeta, sep = "")
  
  BetaTotal <- cbind(BetaEM, BetaNM, BetaOR)
  
  SigmaEM <- sapply(monteCarloResults, function(x) {x$Sigma_EM})
  SigmaNM <- sapply(monteCarloResults, function(x) {x$Sigma_NM})
  SigmaOR <- sapply(monteCarloResults, function(x) {x$Sigma_OR})
  SigmaTotal <- cbind(SigmaEM, SigmaNM, SigmaOR)
  
  colnames(SigmaTotal) <- c("Sigma_EM", "Sigma_NM", "Sigma_OR")
  
  Convergence <- sapply(monteCarloResults, function(x) {x$Convergence})
  NumMissing <- sapply(monteCarloResults, function(x) {x$Num_Miss})
  MissingProportion <- sapply(monteCarloResults, function(x) {x$Prop_Miss})
  
  output <- cbind(BetaTotal, SigmaTotal, Convergence, NumMissing, MissingProportion)
  
  return(output)
}

#-----------------------
# Single trial
#-----------------------

# Hyper-parameters
n <- 300

ratio <- 3

pr_non_missing <- 0.8 # does not mean proportion of non-missing rate, but can be used to control the missing rate
uni_radius_1 <- 1 # control how spread-out X1 is
uni_radius_2 <- 1 # control how spread-out X2 is
# std <- abs(uni_radius_1)/ratio # standard deviation (sigma) of the linear data model
std <- 1
bn <- 2 # interior knots
q <- 2 # order of basis-spline
gHNodes <- 8 # Gauss-Hermite nodes
max_iter <- 500
tol <- 1e-4

# Data model coefficients
# coef_intercept <- log(pr_non_missing/(1-pr_non_missing))
coef_intercept <- 1
coef_x1 <- uni_radius_1
coef_x2 <- uni_radius_2
# coef1 <- c(coef_intercept, coef_x1, coef_x2)
coef1 <- c(coef_intercept, coef_x1)

# Simulate complete data
# X <- matrix(runif(2*n, min = -1, max = 1), ncol = 2)
X <- matrix(rnorm(n), ncol = 1)
Y <- simuY(cbind(1, X), coef1, std)
Z <- X
U <- NULL # \pi(Y)

# Missing model
yy <- log(Y/(1-Y))
# YU <- cbind(yy, sin(2*pi*yy))
YU <- cbind(1, yy)
coef2 <- c(-1, 1)
Obs <- simuMiss(YU, coef2)

# Check overlap between 0/1 and proportion of missing
dat <- cbind(Y, Obs, Z, U)
# colnames(dat) <- c("Y", "Obs", "X1", "X2")
colnames(dat) <- c("Y", "Obs", "X1")

dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = yy))+geom_density(aes(fill = OBS), alpha = 0.5)

# EM algorithm
# df_MNAR <- list(data = dat, Z_indices = c(3, 4), U_indices = NULL)
df_MNAR <- list(data = dat, Z_indices = 3, U_indices = NULL)
emEstimate <- main(df_MNAR, 2*coef1, 2*std, runif(bn+q), bn, q, gHNodes, max_iter, tol)

# Non-missing data
yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]
yLM <- log(yObs/(1-yObs))
obsLM <- lm(yLM~xObs)

(beta_non_missing <- obsLM$coefficients)
(sigma_non_missing <- sigma(obsLM))

# EM results
emEstimate$Beta
emEstimate$Sigma

# Oracle
oracleLM <- lm(yy~X)
(beta_oracle <- oracleLM$coefficients)
(sigma_oracle <- sigma(oracleLM))

# True values
coef1
std

# Data summary
print(paste("There are", table(Obs)[1],"out of", n,"missing observations,",
            paste("the proportion of missing is",
                  paste(100*table(Obs)[1]/n,"%.", sep=""))))

#-----------------------
# Monte Carlo Using the same setting as the single trial
# mclapply not working under Windows (use mclapply)
#-----------------------
# 
# B <- 30
# 
# df_MNAR_list <- lapply(1:B, function(x)
# {
#   X <- matrix(rnorm(n), ncol = 1)
#   Y <- simuY(cbind(1, X), coef1, std)
#   Z <- X
#   U <- NULL # \pi(Y)
# 
#   # Missing model
#   yy <- log(Y/(1-Y))
#   # YU <- cbind(yy, sin(2*pi*yy))
#   YU <- cbind(1, yy)
#   coef2 <- c(-1, 1)
#   Obs <- simuMiss(YU, coef2)
# 
#   # Check overlap between 0/1 and proportion of missing
#   dat <- cbind(Y, Obs, Z, U)
#   # colnames(dat) <- c("Y", "Obs", "X1", "X2")
#   colnames(dat) <- c("Y", "Obs", "X1")
# 
#   df_MNAR <- list(data = dat, Z_indices = 3, U_indices = NULL)
# 
#   return(df_MNAR)
# }
# )
# 
# monteCarloResults <- lapply(df_MNAR_list, function(x)
# {
#   df_MNAR <- x
#   emEstimate <- main(df_MNAR, 2*coef1, 2*std, runif(bn+q), bn, q, gHNodes, max_iter, tol)
# 
#   varEst <- ProfileCov(df_MNAR, sqrt(min(abs(vcov(obsLM)))),
#                        emEstimate$Beta, emEstimate$Sigma, emEstimate$Tau,
#                        bn, q, gHNodes, bn+q, 1, max_iter, 1e-4)
#   
#   betaCIs <- matrix(ncol = 4)
#   sd_b0 <- sqrt(varEst[1,1])
#   sd_b1 <- sqrt(varEst[2,2])
#   betaCIs[1] <- emEstimate$Beta[1]-1.96*sd_b0
#   betaCIs[2] <- emEstimate$Beta[1]+1.96*sd_b0
#   
#   betaCIs[3] <- emEstimate$Beta[2]-1.96*sd_b1
#   betaCIs[4] <- emEstimate$Beta[2]+1.96*sd_b1
#   
#   dat <- df_MNAR$data
#   Obs <- dat[,"Obs"]
#   Y <- dat[,"Y"]
#   X <- dat[,-c(1,2)]
# 
#   yObs <- Y[which(Obs == 1)]
# # xObs <- X[which(Obs == 1),]
#   xObs <- as.matrix(X)[which(Obs == 1),]
# 
#   yLM <- log(yObs/(1-yObs))
#   obsLM <- lm(yLM~xObs)
# 
#   beta_non_missing <- obsLM$coefficients
#   sigma_non_missing <- sigma(obsLM)
# 
#   yy <- log(Y/(1-Y))
#   oracleLM <- lm(yy~X)
#   beta_oracle <- oracleLM$coefficients
#   sigma_oracle <- sigma(oracleLM)
# 
#   return(list(Beta_EM = emEstimate$Beta, Sigma_EM = emEstimate$Sigma, Convergence = emEstimate$Success,
#               Beta_NM = beta_non_missing, Sigma_NM = sigma_non_missing,
#               Beta_OR = beta_oracle, Sigma_OR = sigma_oracle,
#               Num_Miss = table(Obs)[1], Prop_Miss = table(Obs)[1]/length(Obs),
#               CI = betaCIs))
# }
# )
# 
# MCResults <- monteCarloInDataFrame(monteCarloResults)
# ThetaMat <- MCResults[,c("BetaEM_1", "BetaEM_2", "Sigma_EM")]
# 
# CIMat <- sapply(monteCarloResults, function(x) {x$CI})
# CIMat <- CIMat[,complete.cases(t(CIMat))]
# Count0 <- 0
# Count1 <- 0
# for(i in 1:ncol(CIMat))
# {
#   if (CIMat[1,i] <= 1 & 1 <= CIMat[2,i])
#   {
#     Count0 <- Count0+1
#   }
#   if (CIMat[3,i] <= 1 & 1 <= CIMat[4,i])
#   {
#     Count1 <- Count1+1
#   }
# }
# cp <- list(Count0/ncol(CIMat), Count1/ncol(CIMat))

# Empirical covariance matrix
cov(ThetaMat)

# Covariance matrix using profile likelihood
(ProfileCov(df_MNAR, sqrt(min(diag(vcov(obsLM)))), emEstimate$Beta, emEstimate$Sigma, emEstimate$Tau, bn, q, gHNodes, bn+q, 1, max_iter, 1e-4) -> varEst)
