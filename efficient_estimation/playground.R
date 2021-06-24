library(fastGHQuad)
library(tidyverse)
source("simulate_data.R")
Rcpp::sourceCpp("bspline_recursive.cpp")
source("add_splines.R")
source("estep.R")
source("mstep.R")

# # Test data generating functions
# 
# n <- 100
# X <- matrix(rnorm(2*n), ncol = 2)
# Z <- X[,1]
# U <- X[,2]
# 
# coef1 <- c(1, -1, 1)
# Y <- simuY(cbind(1, X), coef1, 1)
# 
# YU <- cbind(1, Y, U)
# coef2 <- c(1, -2, 1)
# Obs <- simuMiss(YU, coef2)
# 
# dat <- cbind(Y, Obs, Z, U)
# colnames(dat) <- c("Y", "Obs", "Z", "U")
# 
# df_MNAR <- list(data = dat, Z_indices = 3, U_indices = 4)
# 
# AppendSplines(df_MNAR, 3,3) -> datList
# 
# # Test "YSplinePrep"
# 
# xVec <- c(1, 1.5, 2)
# rules <- gaussHermiteData(10)
# q <- 3
# bn <- 3
# delta <- seq(0, 1, length.out = bn+2)
# delta <- c(rep(0, q-1), delta, rep(1, q+1))
# betaOld <- c(1, 1, 1)
# sigmaOld <- 1
# 
# YSplinePrep(xVec, rules, delta, bn, q, betaOld, sigmaOld)
# 
# # Test YUSplinePrep
# 
# x <- datList$bs_u1[1,]
# xBS <- as.matrix(x)
# yBS <- datList$bs_y[1,]
# 
# YUSplinePrep(yBS, xBS)
# 
# # Test PseudoObservation
# 
# yVec <- runif(10)
# w <- rules$w
# ubsMat <- YUSplinePrep(yBS, xBS)
# tauOld <- rnorm(length(ubsMat))
# 
# PseudoObservation(yVec, w, 2, 2.4, xBS,  delta, bn, q, tauOld)
# 
# # Test Make Dataset
# 
# df_MNAR <- list(data = dat, Z_indices = 3, U_indices = 4)
# AppendSplines(df_MNAR, bn = 3, q = 3) -> datList
# df1 <- MakeFullDataSetObs(datList)
# df2 <- MakeFullDataSetMissing(datList, rules, betaOld, sigmaOld, tauOld)

# Test MStep function

# datList <- AppendSplines(df_MNAR, 3, 3)
# rules <- fastGHQuad::gaussHermiteData(10)
# 
# df1 <- MakeFullDataSetObs(datList)
# df2 <- MakeFullDataSetMissing(datList,rules,c(0.8, -0.5, 0.8), 0.7, rnorm(6))
# 
# nu <- 0
# nz <- 2
# nsieve <- 6
# 
# MStep(df1, df2, nu, nz, nsieve)

## Generate data
n <- 2000
X <- matrix(runif(2*n, min = -1, max = 1), ncol = 2)
Z <- X
U <- NULL

## Data model
coef1 <- c(1, -2, 2)
std <- 1
Y <- simuY(cbind(1, X), coef1, std) # (0, 1)
yy <- log(Y/(1-Y)) # (-inf, inf)

## Missing model
YU <- cbind(yy, sin(2*pi*yy))
coef2 <- c(1, 1)
Obs <- simuMiss(YU, coef2) # 1: observed; 0: missing

table(Obs)

yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]

# Only use non-missing data
yLM <- log(yObs/(1-yObs))

obsLM <- lm(yLM~xObs)
beta_init <- obsLM$coefficients
sigma_init <- sigma(obsLM)

# Oracle
oracleLM <- lm(yy~X)
beta_oracle <- oracleLM$coefficients
sigma_oracle <- sigma(oracleLM)

dat_df <- as.data.frame(dat)
dat_df %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>%
  ggplot(aes(x = yy))+geom_density(aes(fill = OBS), alpha=0.5)

# EM
dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", "X1", "X2")
df_MNAR <- list(data = dat, Z_indices = c(3, 4), U_indices = NULL)

bn <- 3
q <- 3

EM_estimate <- main(df_MNAR, beta_init*3, sigma_init*3, rep(0, bn+q),
                    bn, q, gaussHermiteNodes = 8, tol = 1e-4)

# Proportion of missing
table(Obs)

# EM algorithm
EM_estimate[1:2]

# Only non-missing data
beta_init
sigma_init

# Oracle (use all data)
beta_oracle
sigma_oracle

# True
coef1
std

ProfileEM(df_MNAR, EM_estimate$Beta, EM_estimate$Sigma, rep(0, bn+q),
          bn, q, gaussHermiteNodes = 8, tol = 1e-4)
