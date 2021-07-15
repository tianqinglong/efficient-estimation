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


