library(fastGHQuad)
library(tidyverse)
library(parallel)
source("../simulate_data.R")
Rcpp::sourceCpp("../bspline_recursive.cpp")
source("../add_splines.R")
source("../estep.R")
source("../mstep.R")
source("../loglikelihood.R")

# Simulation 1

ghn <- 8
max_iter <- 350
tol <- 1e-5

n <- 1000
mu <- 0.5
sigX <- 0.5

beta0 <- 0.25
beta1 <- -0.5
sigY <- 1

coef1 <- c(beta0, beta1)

X <- matrix(rnorm(n, mu, sigX), ncol = 1)
Y <- simuY(cbind(1, X), coef1, sigY)
Z <- X
U <- NULL # \pi(Y)

yy <- log(Y/(1-Y))
YU <- cbind(1, yy)
coef2 <- c(1, 1)
Obs <- simuMiss(YU, coef2)
dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", "X1")

# Missing at random
yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]
yLM <- log(yObs/(1-yObs))
lmMAR <- lm(yLM~xObs)
hn <- min(sqrt(diag(vcov(lmMAR))))

# EM
df_MNAR <- list(data = dat, Z_indices = 3, U_indices = NULL)

bn <- 4
q <- 2
nsieves <- bn+q
emEstimate <- main(df_MNAR, coef1+0.1, sigY, runif(bn+q, min = -1, max = 1), bn, q, ghn, max_iter, tol, fixed_sigma = T)
emEstimate
lmMAR

# Simulation 2
n <- 1000
Mu <- c(0, 0, 0)
Sigma <- matrix(nrow = 3, ncol = 3)
for (i in 1:3)
  for (j in 1:3)
  {
    Sigma[i,j] <- 0.5^(abs(i-j))
  }

X <- MASS::mvrnorm(n, Mu, Sigma)
sigY <- 1
coef1 <- c(0, 0.1, -0.2, -0.3)
Y <- simuY(cbind(1, X), coef1, sigY)

Z <- X
U <- NULL # \pi(Y)

yy <- log(Y/(1-Y))
YU <- cbind(1, yy)
coef2 <- c(1, 1)
Obs <- simuMiss(YU, coef2)
dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", "X1", "X2", "X3")

# Missing at random
yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]
yLM <- log(yObs/(1-yObs))
lmMAR <- lm(yLM~xObs)
hn <- min(sqrt(diag(vcov(lmMAR))))

df_MNAR <- list(data = dat, Z_indices = c(3,4,5), U_indices = NULL)
bn <- 4
q <- 2
nsieves <- bn+q
emEstimate <- main(df_MNAR, coef1+0.1, sigY, runif(bn+q, min = -1, max = 1), bn, q, ghn, max_iter, tol, fixed_sigma = T)
emEstimate
