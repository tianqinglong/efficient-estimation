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

## Setting 1

n <- 1000
coef1 <- c(1,1)
sd <- 1
ghn <- 8
bn <- 2
q <- 2
max_iter <- 500
tol <- 1e-5


X <- matrix(truncnorm::rtruncnorm(n), ncol = 1)
Y <- simuY(cbind(1, X), coef1, sd)
Z <- X
U <- NULL # \pi(Y)
nsieves <- bn+q

yy <- log(Y/(1-Y))
YU <- cbind(1, yy, yy^2)
coef2 <- c(1, -1, .5)
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

df_MNAR <- list(data = dat, Z_indices = 3, U_indices = NULL)

opt1 <- optim(c(0.5, 0.5, 0.5), pseudo_loglikelihood_tian1, df_MNAR = df_MNAR, ghxw = gaussHermiteData(8), hessian = T)

(emEstimate <- main(df_MNAR, coef1+0.5, sd+0.5, runif(bn+q, min = -1, max = 1), bn, q, ghn, max_iter, tol))
# Compute covariance matrix
beta_mle <- emEstimate$Beta
sd_mle <- emEstimate$Sigma
tau_mle <- emEstimate$Tau

temp <- ProfileCov(df_MNAR, min(hn, 1/sqrt(n)), beta_mle, sd_mle, tau_mle, bn, q, ghn, nsieves, 1, max_iter, tol)
temp
cbind(beta_mle-qnorm(0.975)*sqrt(diag(temp))[1:length(beta_mle)],
      beta_mle+qnorm(0.975)*sqrt(diag(temp))[1:length(beta_mle)])

opt1$par-qnorm(0.975)*sqrt(diag(solve(opt1$hessian)))
opt1$par+qnorm(0.975)*sqrt(diag(solve(opt1$hessian)))

## Setting 2

n <- 1000
coef1 <- c(-1, 4, -2)
sd <- 1

Z <- matrix(rnorm(n), ncol = 1)
U <- matrix(rnorm(n, 1-Z, 1), ncol = 1)
X <- cbind(Z, U)
Y <- simuY(cbind(1, X), coef1, sd)

n_covarites <- ncol(X)
yy <- log(Y/(1-Y))
YU <- cbind(1, yy, U)
coef2 <- c(4, -3, 1)
Obs <- simuMiss(YU, coef2, use_logit = T)
table(Obs)

dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", paste("X", 1:(ncol(dat)-2), sep = ""))
dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = U))+geom_density(aes(fill = OBS), alpha = 0.5)
df_MNAR <- list(data = dat, Z_indices = 3, U_indices = 4)

# opt2 <- optim(c(coef1, sd), pseudo_loglikelihood_tian2, df_MNAR = df_MNAR, ghxw = gaussHermiteData(20), hessian = T)

# opt2 <- optim(c(coef1, sd), pseudo_loglikelihood_tian2, df_MNAR = bad_df, ghxw = gaussHermiteData(20), hessian = T)

# EM Method

yObs <- Y[which(Obs == 1)]
xObs <- X[which(Obs == 1),]
yLM <- log(yObs/(1-yObs))
lmMAR <- lm(yLM~xObs)
hn <- min(sqrt(diag(vcov(lmMAR))))

ghn <- 10
bn <- 3
q <- 2
max_iter <- 500
tol <- 1e-4

# Additive model
nsieves_add <- (bn+q-1)*(ncol(U)+1)+1
emEstimate <- main_additive(df_MNAR, coef1+0.1, sd+0.1, rep(0, times = nsieves_add), bn, q, ghn, max_iter, tol)
emEstimate

beta_mle <- emEstimate$Beta
sd_mle <- emEstimate$Sigma
tau_mle <- emEstimate$Tau

temp <- ProfileCov_additive(df_MNAR, min(hn, 1/sqrt(n))*2, beta_mle, sd_mle, tau_mle, bn, q, ghn, nsieves_add, n_covarites, max_iter, tol)
cbind(beta_mle-qnorm(0.975)*sqrt(diag(temp))[1:length(beta_mle)],
      beta_mle+qnorm(0.975)*sqrt(diag(temp))[1:length(beta_mle)])

opt2$par-qnorm(0.975)*sqrt(diag(solve(opt2$hessian)))
opt2$par+qnorm(0.975)*sqrt(diag(solve(opt2$hessian)))

## Setting 4

n <- 1000
coef1 <- c(-1, 4, -1, -1)
sd <- 1

ghn <- 8
bn <- 4
q <- 2
max_iter <- 500
tol <- 1e-5

Z <- matrix(rnorm(n), ncol = 1)
U <- matrix(rnorm(2*n, 1-Z, 1), ncol = 2)
X <- cbind(Z, U)
Y <- simuY(cbind(1, X), coef1, sd)

n_covarites <- ncol(X)
yy <- log(Y/(1-Y))
YU <- cbind(1, yy, U)
coef2 <- c(4, -3, 1, 1)
Obs <- simuMiss(YU, coef2, use_logit =T)
table(Obs)

dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", paste("X", 1:(ncol(dat)-2), sep = ""))
dat %>% as.data.frame %>% mutate(OBS = as.factor(Obs), yy = log(Y/(1-Y))) %>% 
  ggplot(aes(x = yy))+geom_density(aes(fill = OBS), alpha = 0.5)
df_MNAR <- list(data = dat, Z_indices = 3, U_indices = c(4, 5))

opt4 <- optim(c(-0.5, 3, -0.5, -0.5, 0.5), pseudo_loglikelihood_tian4, df_MNAR = df_MNAR, ghxw = gaussHermiteData(8), hessian = T)
opt4$par-qnorm(0.975)*sqrt(diag(solve(opt4$hessian)))
opt4$par+qnorm(0.975)*sqrt(diag(solve(opt4$hessian)))

nsieves_add <- (bn+q-1)*(ncol(U)+1)+1
emEstimate <- main_additive(df_MNAR, coef1+0.1, sd+0.1, rep(0, nsieves_add), bn, q, ghn, max_iter, tol)
emEstimate

beta_mle <- emEstimate$Beta
sd_mle <- emEstimate$Sigma
tau_mle <- emEstimate$Tau

temp1 <- ProfileCov_additive(df_MNAR, min(hn, 1/sqrt(n)), beta_mle, sd_mle, tau_mle, bn, q, ghn, nsieves_add, ncovariate, max_iter, tol)
cbind(beta_mle-qnorm(0.975)*sqrt(diag(temp1))[1:length(beta_mle)],
      beta_mle+qnorm(0.975)*sqrt(diag(temp1))[1:length(beta_mle)])

coef1 <- c(-1, 4, -1, -1)
sd <- 1

rout <- readRDS("no_git/rout_tian4_final.rds")
analysis(rout, c(coef1, sd))

## Draw probability

coef1 <- c(1,1)
sd <- 1
rout <- readRDS("no_git/rout_tian1_final.rds")

q <- 2
tauMat <- lapply(rout, function(x) {x$EM$Tau})
df_List <- lapply(rout, function(x) {x$Data})

max_of_min_Y <- -Inf
min_of_max_Y <- Inf
y_all <- NULL
for (i in 1:length(tauMat))
{
  yVec <- rout[[i]]$Data$data[,"Y"]
  max_of_min_Y <- max(max_of_min_Y, min(yVec))
  min_of_max_Y <- min(min_of_max_Y, max(yVec))
  y_all <- c(y_all, yVec)
}

yPseudoVec <- seq(max_of_min_Y, min_of_max_Y, length.out = 1000)
outmat <- NULL

for (i in 1:length(tauMat))
{
  emEstimate <- rout[[i]]$EM
  tauVec <- tauMat[[i]]
  if (length(tauVec) == 4)
  {
    bn <- 2
  }
  else
  {
    bn <- 3
  }
  
  df_MNAR <- df_List[[i]]
  dfspline <- AppendSplines(df_MNAR, bn, q)
  delta <- dfspline$delta
  yOriginal <- log(1/(1-yPseudoVec))
  
  yPseudoSpine <- matrix(nrow = bn+q, ncol = length(yPseudoVec))
  for (i in 1:length(yPseudoVec))
  {
    yPseudoSpine[,i] <- BSplinesForNewY(yPseudoVec[i], q, bn, delta)
  }
  yPseudoSpine <- t(yPseudoSpine)
  nonMissingProb <- yPseudoSpine %*% emEstimate$Tau
  nonMissingProb <- 1/(1+exp(-nonMissingProb))
  outmat <- cbind(outmat, nonMissingProb)
}

ci <- t(apply(outmat, MARGIN = 1, FUN = function(x) {quantile(x, probs = c(0.025, 0.975, 0.5))}))
nonMissingProb_true <- 1-yOriginal+0.5*yOriginal^2
nonMissingProb_true <- 1/(1+exp(-nonMissingProb_true))

plot(x = yOriginal, y = ci[,1], type = "l",
     xlab = "Y", ylab = expression(paste(pi,"(Y)",sep = "")), lty = 2, col = 2, ylim = c(0.5, 1), lwd = 2)
lines(x = yOriginal, y = nonMissingProb_true, lty = 1, col = 1, lwd = 2)
lines(x = yOriginal, y = ci[,2], lty = 2, col = 2, lwd = 2)
lines(x = yOriginal, y = ci[,3], lty = 4, col = 4, lwd = 2.2)
dsty <- density(log(y_all/(1-y_all)))
legend(x = c(3.4, 4.5), y = c(0.5, .65),c("True",  "Median", "95% Band"), lty = c(1, 4,2), col = c(1,4,2), lwd = c(1,1.25,1))

# Draw proportion of missing

dfspline <- AppendSplines(df_MNAR, bn, q)
delta <- dfspline$delta
yPseudoVec <- seq(range(Y)[1], range(Y)[2], length.out = 1000)

yOriginal <- log(1/(1-yPseudoVec))

yPseudoSpine <- matrix(nrow = bn+q, ncol = length(yPseudoVec))
for (i in 1:length(yPseudoVec))
{
  yPseudoSpine[,i] <- BSplinesForNewY(yPseudoVec[i], q, bn, delta)
}
yPseudoSpine <- t(yPseudoSpine)
nonMissingProb <- yPseudoSpine %*% emEstimate$Tau
nonMissingProb <- 1/(1+exp(-nonMissingProb))
plot(x = yOriginal, y = nonMissingProb_true, type="l",
     xlab = "Y", ylab = "Missing Probability", lty = 1, col = 1)
nonMissingProb_true <- 1-yOriginal+0.5*yOriginal^2
nonMissingProb_true <- 1/(1+exp(-nonMissingProb_true))
lines(x = yOriginal, y = nonMissingProb, col =2, lty = 3)
legend("topright",c("Estimated", "True"), lty = c(1, 3), col = c(1,2))
