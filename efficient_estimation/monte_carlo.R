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
source("parallel.R")
#-----------------------
# Auxiliary functions
#-----------------------
monteCarloInDataFrame <- function(monteCarloResults)
{
  BetaEM <- t(sapply(monteCarloResults, function(x) {x$Beta_EM}))
  colnames(BetaEM) <- c("Beta0_EM", "Beta1_EM", "Beta2_EM")
  BetaNM <- t(sapply(monteCarloResults, function(x) {x$Beta_NM}))
  colnames(BetaNM) <- c("Beta0_NM", "Beta1_NM", "Beta2_NM")
  BetaOR <- t(sapply(monteCarloResults, function(x) {x$Beta_OR}))
  colnames(BetaOR) <- c("Beta0_OR", "Beta1_OR", "Beta2_OR")
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
n <- 100
pr_non_missing <- 0.8 # does not mean 10% missing rate, but can be use to control the missing rate
uni_radius_1 <- 2 # control how spread-out X1 is
uni_radius_2 <- 2 # control how spread-out X2 is
std <- 1 # standard deviation (sigma) of the linear data model
bn <- 3 # interior knots
q <- 3 # order of basis-spline
gHNodes <- 9 # Gauss-Hermite nodes
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
Z <- X
U <- NULL # \pi(Y)

# Missing model
yy <- log(Y/(1-Y))
YU <- cbind(yy, sin(2*pi*yy))
coef2 <- c(1, 1)
Obs <- simuMiss(YU, coef2)


# Check overlap between 0/1 and proportion of missing
dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", "X1", "X2")

dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = yy))+geom_density(aes(fill = OBS), alpha = 0.5)

# EM algorithm
df_MNAR <- list(data = dat, Z_indices = c(3, 4), U_indices = NULL)
emEstimate <- main(df_MNAR, 2*coef1, 2*std, runif(bn+q), bn, q, gHNodes, max_iter, tol)

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

B <- 20

df_MNAR_list <- mclapply(1:B, function(x)
{
  X <- matrix(runif(2*n, min = -1, max = 1), ncol = 2)
  Y <- simuY(cbind(1, X), coef1, std)
  Z <- X
  U <- NULL # \pi(Y)
  
  yy <- log(Y/(1-Y))
  YU <- cbind(yy, sin(2*pi*yy))
  coef2 <- c(1, 1)
  Obs <- simuMiss(YU, coef2)
  
  dat <- cbind(Y, Obs, Z, U)
  colnames(dat) <- c("Y", "Obs", "X1", "X2")
  
  df_MNAR <- list(data = dat, Z_indices = c(3, 4), U_indices = NULL)
  
  return(df_MNAR)
}
)

monteCarloResults <- lapply(df_MNAR_list, function(x)
{
  df_MNAR <- x
  emEstimate <- main(df_MNAR, 2*coef1, 2*std, runif(bn+q), bn, q, gHNodes, max_iter, tol)
  
  dat <- df_MNAR$data
  Obs <- dat[,"Obs"]
  Y <- dat[,"Y"]
  X <- dat[,-c(1,2)]
  
  yObs <- Y[which(Obs == 1)]
  xObs <- X[which(Obs == 1),]
  
  yLM <- log(yObs/(1-yObs))
  obsLM <- lm(yLM~xObs)
  
  beta_non_missing <- obsLM$coefficients
  sigma_non_missing <- sigma(obsLM)
  
  yy <- log(Y/(1-Y))
  oracleLM <- lm(yy~X)
  beta_oracle <- oracleLM$coefficients
  sigma_oracle <- sigma(oracleLM)
  
  return(list(Beta_EM = emEstimate$Beta, Sigma_EM = emEstimate$Sigma, Convergence = emEstimate$Success,
              Beta_NM = beta_non_missing, Sigma_NM = sigma_non_missing,
              Beta_OR = beta_oracle, Sigma_OR = sigma_oracle,
              Num_Miss = table(Obs)[1], Prop_Miss = table(Obs)[1]/length(Obs)))
}
)

monteCarloInDataFrame(monteCarloResults)
