library(fastGHQuad)
library(tidyverse)
source("simulate_data.R")
Rcpp::sourceCpp("bspline_recursive.cpp")
source("add_splines.R")
source("estep.R")
source("mstep.R")
source("loglikelihood.R")
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
