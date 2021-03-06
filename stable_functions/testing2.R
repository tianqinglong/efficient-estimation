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
coef1 <- c(1, 2, 2)
sd <- 0.5

X1 <- NULL
X2 <- matrix(truncnorm::rtruncnorm(2*n), ncol = 2)
X <- cbind(X1, X2)
Y <- simuY(cbind(1, X), coef1, sd)
Z <- NULL
U <- X

n_covarites <- ncol(X)
yy <- log(Y/(1-Y))
YU <- cbind(1, yy, U[,1]^2, U[,2]^2)
coef2 <- c(1, 1, 1, -1)
Obs <- simuMiss(YU, coef2, use_logit =F)
table(Obs)

dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", paste("X", 1:(ncol(dat)-2), sep = ""))
dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = yy))+geom_density(aes(fill = OBS), alpha = 0.5)
df_MNAR <- list(data = dat, Z_indices = NULL, U_indices = c(3, 4))

yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]
yLM <- log(yObs/(1-yObs))
lmMAR <- lm(yLM~xObs)
hn <- min(sqrt(diag(vcov(lmMAR))))

# # Multiplicative model
# nsieves <- (bn+q)^2
# (emEstimate0 <- main(df_MNAR, coef1+0.5, sd+0.5, runif(nsieves, min = -1, max = 1), bn, q, ghn, max_iter, tol))

ghn <- 7
bn <- 2
q <- 3
max_iter <- 500
tol <- 1e-4

# Additive model
nsieves_add <- (bn+q-1)*(ncol(U)+1)+1
emEstimate <- main_additive(df_MNAR, coef1+0.1, sd+0.1, rep(0, nsieves_add), bn, q, ghn, max_iter, tol)
emEstimate
lmMAR
(lmBD <- lm(yy~X))
hn <- min(sqrt(diag(vcov(lmBD))))

## Compute covariance matrix
beta_mle <- emEstimate$Beta
sd_mle <- emEstimate$Sigma
tau_mle <- emEstimate$Tau

temp <- ProfileCov_additive(df_MNAR, min(hn, 1/sqrt(n)), beta_mle, sd_mle, tau_mle, bn, q, ghn, nsieves_add, n_covarites, max_iter, tol)
cbind(beta_mle-qnorm(0.975)*sqrt(diag(temp))[1:length(beta_mle)],
      beta_mle+qnorm(0.975)*sqrt(diag(temp))[1:length(beta_mle)])
confint(lmMAR)
confint(lmBD)

####################################
coef1 <- c(1, 2, 2)
sd <- 0.5

## My setting 3
rout <- readRDS("no_git/rout_tian3_new.rds")
analysis(rout, c(coef1, sd))

rout <- readRDS("no_git/rout_tian3_new_small_sd.rds")
analysis(rout, c(coef1, sd))

rout <- readRDS("no_git/rout_tian3_resample.rds")
analysis(rout, c(coef1, sd))
