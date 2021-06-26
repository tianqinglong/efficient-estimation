#-----------------------
# Try to make \pi(Y, U) work
# Author: Qinglong Tian
# Date: June 24, 2021
#-----------------------
library(fastGHQuad)
library(tidyverse)
source("simulate_data.R")
Rcpp::sourceCpp("bspline_recursive.cpp")
source("add_splines.R")
source("estep.R")
source("mstep.R")
source("estep_additive.R")
source("mstep_additive.R")
#-----------------------

# Hyper-parameters
n <- 20000

ratio <- 3

pr_non_missing <- .9 # does not mean proportion of non-missing rate, but can be used to control the missing rate
uni_radius_1 <- 2 # control how spread-out X1 is
uni_radius_2 <- 1.5 # control how spread-out X2 is
std <- (abs(uni_radius_1)+abs(uni_radius_2))/ratio # standard deviation (sigma) of the linear data model
bn <- 6 # interior knots
q <- 4 # order of basis-spline
gHNodes <- 8 # Gauss-Hermite nodes
max_iter <- 200
tol <- 1e-4

# Data model coefficients
coef_intercept <- log(pr_non_missing/(1-pr_non_missing))
coef_x1 <- uni_radius_1
coef_x2 <- uni_radius_2
coef1 <- c(coef_intercept, coef_x1, coef_x2)

# Simulate complete data
X <- matrix(runif(2*n, min = -1, max = 1), ncol = 2)
Y <- simuY(cbind(1, X), coef1, std)
Z <- X[,1] # Z first, U second
U <- X[,2] # \pi(Y,U)

# Missing model

hyper <- log(U/2+1/2)-log(1-U/2)
yy <- log(Y/(1-Y))
YU <- cbind(yy, hyper)
coef2 <- c(2, 2)
Obs <- simuMiss(YU, coef2)

# Check overlap between 0/1 and proportion of missing
dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", "X1", "X2")

dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = yy))+geom_density(aes(fill = OBS), alpha = 0.5)

# Data summary
print(paste("There are", table(Obs)[1],"out of", n,"missing observations,",
            paste("the proportion of missing is",
                  paste(100*table(Obs)[1]/n,"%.", sep=""))))

# Run EM algorithm
df_MNAR <- list(data = dat, Z_indices = 3, U_indices = 4)
# emEstimate <- main(df_MNAR, coef1+0.1, std+0.1, runif((bn+q)^2, min = -0.5, max = .5), bn, q, gHNodes, max_iter, tol)
emEstimate_add <- main_add(df_MNAR, coef1+0.1, std+0.1, runif((bn+q)*2, min = -0.5, max = .5), bn, q, gHNodes, max_iter, tol)


# Non-missing data
yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]
yLM <- log(yObs/(1-yObs))
obsLM <- lm(yLM~xObs)

(beta_non_missing <- obsLM$coefficients)
(sigma_non_missing <- sigma(obsLM))

# Oracle
oracleLM <- lm(yy~X)
(beta_oracle <- oracleLM$coefficients)
(sigma_oracle <- sigma(oracleLM))

# EM
emEstimate_add$Beta
emEstimate_add$Sigma

# True
coef1
std