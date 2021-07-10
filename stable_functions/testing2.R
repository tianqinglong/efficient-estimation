library(fastGHQuad)
library(tidyverse)
source("simulate_data.R")
Rcpp::sourceCpp("bspline_recursive.cpp")
source("add_splines.R")
source("estep.R")
source("mstep.R")
source("loglikelihood.R")
source("analysis.R")
####################################

n <- 1000
coef1 <- c(1, 1, 1, 1)
sd <- 1
ghn <- 8
bn <- 2
q <- 2
max_iter <- 500
tol <- 1e-5

X1 <- matrix(truncnorm::rtruncnorm(n, b = .5), ncol = 1)
X2 <- matrix(truncnorm::rtruncnorm(2*n, a = -1, b = 1), ncol = 2)
X <- cbind(X1, X2)
Y <- simuY(cbind(1, X), coef1, sd)
Z <- X1
U <- X2

n_covarites <- ncol(X)
yy <- log(Y/(1-Y))
YU <- cbind(1, yy, U)
coef2 <- c(-.5, .75, .75, .75)
Obs <- simuMiss(YU, coef2, use_logit = T)
table(Obs)

dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", "X1", "X2", "X3")
dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = yy))+geom_density(aes(fill = OBS), alpha = 0.5)
df_MNAR <- list(data = dat, Z_indices = 3, U_indices = c(4,5))

yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]
yLM <- log(yObs/(1-yObs))
lmMAR <- lm(yLM~xObs)
hn <- min(sqrt(diag(vcov(lmMAR))))

# # Multiplicative model
# nsieves <- (bn+q)^2
# (emEstimate0 <- main(df_MNAR, coef1+0.5, sd+0.5, runif(nsieves, min = -1, max = 1), bn, q, ghn, max_iter, tol))

# Additive model
nsieves_add <- (bn+q-1)*(ncol(U)+1)+1
emEstimate <- main_additive(df_MNAR, coef1+0.5, sd+0.5, runif(nsieves_add, min = -1, max = 1), bn, q, ghn, max_iter, tol)
emEstimate$Beta
lmMAR
(lmBD <- lm(yy~X))

## Compute covariance matrix
beta_mle <- emEstimate$Beta
sd_mle <- emEstimate$Sigma
tau_mle <- emEstimate$Tau

temp <- ProfileCov_additive(df_MNAR, min(hn, 1/sqrt(n)), beta_mle, sd_mle, tau_mle, bn, q, ghn, nsieves_add, n_covarites, max_iter, tol)
temp
emEstimate
lmMAR