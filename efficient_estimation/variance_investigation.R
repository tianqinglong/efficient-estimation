library(fastGHQuad)
library(tidyverse)
library(parallel)
source("simulate_data.R")
Rcpp::sourceCpp("bspline_recursive.cpp")
source("add_splines.R")
source("estep.R")
source("mstep.R")
source("variance.R")

n <- 300
coef1 <- c(1,1)
std <- 1
ghn <- 8
bn <- 2
q <- 2
max_iter <- 300
tol <- 1e-7

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
colnames(dat) <- c("Y", "Obs", "X1")
dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = X))+geom_density(aes(fill = OBS), alpha = 0.5)

df_MNAR <- list(data = dat, Z_indices = 3, U_indices = NULL)
emEstimate <- main(df_MNAR, 2*coef1, 2*std, runif(bn+q, min = -1, max = 1), bn, q, ghn, max_iter, tol)

#
nsieves <- bn+q
beta <- emEstimate$Beta
sd <- emEstimate$Sigma
tau <- emEstimate$Tau
computeLogLikelihood(df_MNAR, beta, sd, tau, ghn, nsieves, bn, q)

em1 <- main(df_MNAR, beta+c(0.01, 0.01), sd, emEstimate$Tau, bn, q, ghn, max_iter, tol, T)
computeLogLikelihood(df_MNAR, beta+c(0.01, 0.01), sd, em1$Tau, ghn, nsieves, bn, q)

em2 <- main(df_MNAR, beta-c(0.01, 0.01), sd, emEstimate$Tau, bn, q, ghn, max_iter, tol, T)
computeLogLikelihood(df_MNAR, beta-c(0.01, 0.01), sd, em2$Tau, ghn, nsieves, bn, q)

#
ProfileCov(df_MNAR, 1/sqrt(n), beta, sd, tau, bn, q, ghn, bn+q, 1, max_iter, tol)
