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
coef1 <- c(2, 2)
sd <- 1
ghn <- 8
bn <- 2
q <- 2
max_iter <- 500
tol <- 1e-5

X <- matrix(truncnorm::rtruncnorm(n), ncol = 1)
Y <- simuY(cbind(1, X), coef1, sd)
Z <- NULL
U <- X

yy <- log(Y/(1-Y))
YU <- cbind(1, yy, U)
coef2 <- c(-1, 2, -2)
Obs <- simuMiss(YU, coef2, use_logit = T)
table(Obs)

dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", "X1")
dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = X))+geom_density(aes(fill = OBS), alpha = 0.5)

yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]
yLM <- log(yObs/(1-yObs))
lmMAR <- lm(yLM~xObs)
hn <- min(sqrt(diag(vcov(lmMAR))))

# EM algorithm

df_MNAR <- list(data = dat, Z_indices = NULL, U_indices = 3)
nsieves_add <- (bn+q-1)*(ncol(U)+1)+1
(emEstimate <- main_additive(df_MNAR, coef1+0.1, sd, rep(0, times = nsieves_add), bn, q, ghn, max_iter, tol, fixed_sigma = T))
# emEstimate <- main_additive(df_MNAR, coef1+0.1, sd, rep(0, times = nsieves_add), bn, q, ghn, max_iter, tol)
emEstimate$Beta -> beta_mle
emEstimate$Sigma -> sigma_mle
emEstimate$Tau -> tau_mle

n_covariate <- ifelse(is.null(ncol(Z)), 0, ncol(Z))+ifelse(is.null(ncol(U)), 0, ncol(U))
ProfileCov_additive_sigma_fixed(df_MNAR, hn, beta_mle, sigma_mle, tau_mle,
                                bn, q, ghn, nsieves_add, n_covariate, max_iter, tol) -> cov_mat
# 
# ProfileCov_additive(df_MNAR, hn, beta_mle, sigma_mle, tau_mle,
#                     bn, q, ghn, nsieves_add, n_covariate, max_iter, tol) -> cov_mat

Beta <- beta_mle
Upper <- beta_mle+1.96*sqrt(diag(cov_mat))
Lower <- beta_mle-1.96*sqrt(diag(cov_mat))
t(rbind(Beta, Lower, Upper))
# MAR

yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]
yLM <- log(yObs/(1-yObs))
lmMAR <- lm(yLM~xObs)
confint(lmMAR)
