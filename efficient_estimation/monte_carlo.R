#-----------------------
# Simulate N samples and compute the empirical covariance matrix
# Both covariates, X_1 and X_2, are centered at 0
# Use the intercept to control the proportion of missing
# Missing mechanism is to use the sigmoid function, could add some variation by adding sin(2*pi*y)
#-----------------------
library(fastGHQuad)
library(tidyverse)
library(parallel)
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
gHNodes <- 9 # Gauss-Hermite nodes
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
X <- matrix(truncnorm::rtruncnorm(n, b = 0.7), ncol = 1)
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

#dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
#  ggplot(aes(x = yy))+geom_density(aes(fill = OBS), alpha = 0.5)

dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = X))+geom_density(aes(fill = OBS), alpha = 0.5)

# EM algorithm
# df_MNAR <- list(data = dat, Z_indices = c(3, 4), U_indices = NULL)
df_MNAR <- list(data = dat, Z_indices = 3, U_indices = NULL)
emEstimate <- main(df_MNAR, 2*coef1, 2*std, runif(bn+q, min = -0.4, max = 0.4), bn, q, gHNodes, max_iter, tol)

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

# Covariance matrix using profile likelihood
(varEst <- ProfileCov(df_MNAR, 1/sqrt(n), emEstimate$Beta, emEstimate$Sigma, emEstimate$Tau, bn, q, gHNodes, bn+q, 1, max_iter, 1e-4))

#-----------------------
# Dump probalematic samples
#-----------------------

total <- 1
cnt <- 1
B <- 300
n <- 300

beta_0 <- 1
beta_1 <- 1
std <- 1
bn <- 2
q <- 2
gHNodes <- 10 # Gauss-Hermite nodes
max_iter <- 500
tol <- 1e-4

coef1 <- c(beta_0, beta_1)

df_MNAR_list <- NULL
while (total <= B) {
  X <- matrix(rnorm(n), ncol = 1)
  Y <- simuY(cbind(1, X), coef1, std)
  Z <- X
  U <- NULL
  
  yy <- log(Y/(1-Y))
  YU <- cbind(1, yy)
  coef2 <- c(-1, 1)
  Obs <- simuMiss(YU, coef2)
  
  dat <- cbind(Y, Obs, Z, U)
  colnames(dat) <- c("Y", "Obs", "X1")
  
  df_MNAR <- list(data = dat, Z_indices = 3, U_indices = NULL)
  skip_to_next <- FALSE
  tryCatch(expr = {
    emEstimate <- main(df_MNAR, coef1+.1, std+.1, runif(bn+q), bn, q, gHNodes, max_iter, tol)
    varEst <- ProfileCov(df_MNAR, 1/sqrt(n), emEstimate$Beta, emEstimate$Sigma, emEstimate$Tau, bn, q, gHNodes, bn+q, 1, max_iter, tol)
    },
    warning=function(w){
      skip_to_next <<- TRUE
      },
    finally = {
      if(any(diag(varEst)<0))
      {
        skip_to_next <<- TRUE
      }
    }
    )
  cnt <- cnt+1
  if(skip_to_next)
  {
    next
  }
  df_MNAR_list[[total]] <- list(df = df_MNAR, EM = emEstimate, Var = diag(varEst))
  total <- total+1
}

# saveRDS(df_MNAR_list, file = "../results.rds")

count0 <- 0
count1 <- 0
count2 <- 0
for (i in 1:length(df_MNAR_list))
{
  res <- df_MNAR_list[[i]]
  # beta0
  beta0 <- res$EM$Beta[1]
  sd_beta0 <- sqrt(res$Var[1])
  count0 <- count0+( (beta_0 >= beta0-1.96*sd_beta0) & (beta_0 <= beta0+1.96*sd_beta0) )
  
  # beta1
  beta1 <- res$EM$Beta[2]
  sd_beta1 <- sqrt(res$Var[2])
  count1 <- count1+( (beta_1 >= beta1-1.96*sd_beta1) & (beta_1 <= beta1+1.96*sd_beta1) )
  
  # sigma
  sigma0 <- res$EM$Sigma
  sd_sigma0 <- sqrt(res$Var[3])
  count2 <- count2+( (std >= sigma0-1.96*sd_sigma0) & (std <= sigma0+1.96*sd_sigma0) )
}
count0/length(df_MNAR_list)
count1/length(df_MNAR_list)
count2/length(df_MNAR_list)

#-----------------------
# Monte Carlo Using the same setting as the single trial
# mclapply not working under Windows (use mclapply)
#-----------------------
# 
# B <- 40
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
#   varEst <- ProfileCov(df_MNAR, 1/2/sqrt(n),
#                        emEstimate$Beta, emEstimate$Sigma, emEstimate$Tau,
#                        bn, q, gHNodes, bn+q, 1, max_iter, tol)
# 
#   betaCIs <- matrix(ncol = 6)
#   sdVec <- sqrt(diag(varEst))
#     
#   betaCIs[1] <- emEstimate$Beta[1]-1.96*sdVec[1]
#   betaCIs[2] <- emEstimate$Beta[1]+1.96*sdVec[1]
# 
#   betaCIs[3] <- emEstimate$Beta[2]-1.96*sdVec[2]
#   betaCIs[4] <- emEstimate$Beta[2]+1.96*sdVec[2]
# 
#   betaCIs[5] <- max(emEstimate$Sigma-1.96*sdVec[3], 0)
#   betaCIs[6] <- emEstimate$Sigma+1.96*sdVec[3]
#   
#   return(list(Beta_EM = emEstimate$Beta, Sigma_EM = emEstimate$Sigma, Convergence = emEstimate$Success,
#               Beta_NM = beta_non_missing, Sigma_NM = sigma_non_missing,
#               Beta_OR = beta_oracle, Sigma_OR = sigma_oracle,
#               Num_Miss = table(Obs)[1], Prop_Miss = table(Obs)[1]/length(Obs),
#               CIs = betaCIs)
#          )
# }
# )
# 
# MCResults <- monteCarloInDataFrame(monteCarloResults)
# ThetaMat <- MCResults[,c("BetaEM_1", "BetaEM_2", "Sigma_EM")]
# 
# CIMat <- sapply(monteCarloResults, function(x) {x$CI})
# CIMat <- CIMat[,complete.cases(t(CIMat[1:4,]))][1:4,]
# Count0 <- 0
# Count1 <- 0
# for(i in 1:ncol(CIMat))
# {
#   if (CIMat[1,i] <= coef_intercept & coef_intercept <= CIMat[2,i])
#   {
#     Count0 <- Count0+1
#   }
#   if (CIMat[3,i] <= coef_x1 & coef_x1 <= CIMat[4,i])
#   {
#     Count1 <- Count1+1
#   }
# }
# cp <- list(Count0/ncol(CIMat), Count1/ncol(CIMat))
# 
# 
# # Empirical covariance matrix
# cov(ThetaMat)
# 
# # Covariance matrix using profile likelihood
# varEst
