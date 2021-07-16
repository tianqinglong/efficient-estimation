#-----------------------
# Functions to implement the pseudolikelihood method
# Author: Qinglong Tian
# Date: July 15, 2021
#-----------------------

## First consider the situation that Z is independent of U
organize_data <- function(df_MNAR)
{
  zIndex <- df_MNAR$Z_indices
  uIndex <- df_MNAR$U_indices
  
  dat <- df_MNAR$data
  dat_non_missing <- dat[dat[,"Obs"] == 1,]
  
  # Outputs
  yVec <- dat_non_missing[,"Y"]
  yVec <- log(yVec/(1-yVec))
  
  if (is.null(zIndex))
  {
    zMat <- NULL
  }
  else
  {
    zMat <- dat_non_missing[,zIndex]
  }
  
  if (is.null(uIndex))
  {
    uMat <- NULL
  }
  else
  {
    uMat <- dat_non_missing[,uIndex]
  }
  
  return(list(yVec = yVec,
              zMat = zMat,
              uMat = uMat))
}

### normal assumption: setting 1
pseudo_loglikelihood_tian1 <- function(theta, df_MNAR, ghxw)
{
  # parameters
  beta <- theta[-length(theta)]
  sigma <- theta[length(theta)]
  
  # extract data
  yzu <- organize_data(df_MNAR)
  yVec <- yzu$yVec
  XObs <- cbind(yzu$zMat, yzu$uMat)
  
  data <- df_MNAR$data
  X <- data[,3]
  
  # Fit X - compute log-likelihood part2
  f0 <- lm(X~1)
  mu_x <- f0$coefficients
  sigma_x <- sigma(f0)
  
  xRaw <- ghxw$x
  wVec <- ghxw$w
  
  xReal <- matrix(sqrt(2)*xRaw*sigma_x+mu_x, ncol = 1)
  muReal <- cbind(1, xReal) %*% matrix(beta, ncol = 1)
  
  log_int <- numeric(length(yVec))
  for (i in 1:length(yVec)) {
    log_int[i] <- log(sum(dnorm(yVec[i], mean = muReal, sd = sigma)*wVec)/sqrt(pi))
  }
  loglik2 <- -sum(log_int)
  
  # Compute log-likelihood part1
  muRealObs <- cbind(1, XObs) %*% matrix(beta, ncol = 1)
  loglik1 <- sum(dnorm(yVec, muRealObs, sigma, log = T))
  
  return(-(loglik1+loglik2))
}

### normal assumption: setting 2
pseudo_loglikelihood_tian2 <- function(theta, df_MNAR, ghxw)
{
  # parameters
  beta <- theta[-length(theta)]
  sigma <- theta[length(theta)]
  
  # extract data
  yzu <- organize_data(df_MNAR)
  yVec <- yzu$yVec
  XObs <- cbind(yzu$zMat, yzu$uMat)
  
  data <- df_MNAR$data
  X <- data[,c(3,4)]
  
  # Fit X - compute log-likelihood part2
  f0 <- lm(X[,1]~X[,2])
  zCoef <- f0$coefficients
  z_sigma <- sigma(f0)
  
  xRaw <- ghxw$x
  wVec <- ghxw$w
  
  log_int <- numeric(length(yVec))
  z_real_mu <- cbind(1, XObs[,2]) %*% matrix(zCoef, ncol=1)
  for (i in 1:length(yVec)) {
    xReal <- matrix(sqrt(2)*xRaw*z_sigma+z_real_mu[i], ncol = 1)
    muReal <- cbind(1, xReal, rep(XObs[i,2], length(xReal))) %*% matrix(beta, ncol = 1)
    log_int[i] <- log(sum(dnorm(yVec[i], mean = muReal, sd = sigma)*wVec)/sqrt(pi))
  }
  loglik2 <- -sum(log_int)
  
  # Compute log-likelihood part1
  muRealObs <- cbind(1, XObs) %*% matrix(beta, ncol = 1)
  loglik1 <- sum(dnorm(yVec, muRealObs, sigma, log = T))
  
  return(-(loglik1+loglik2))
}

### normal assumption: setting 4
pseudo_loglikelihood_tian4 <- function(theta, df_MNAR, ghxw)
{
  # parameters
  beta <- theta[-length(theta)]
  sigma <- theta[length(theta)]
  
  # extract data
  yzu <- organize_data(df_MNAR)
  yVec <- yzu$yVec
  XObs <- cbind(yzu$zMat, yzu$uMat)
  
  data <- df_MNAR$data
  X <- data[,c(3,4,5)]
  
  # Fit X - compute log-likelihood part2
  f0 <- lm(X[,1]~X[,2]+X[,3])
  zCoef <- f0$coefficients
  z_sigma <- sigma(f0)
  
  xRaw <- ghxw$x
  wVec <- ghxw$w
  
  log_int <- numeric(length(yVec))
  z_real_mu <- cbind(1, XObs[,c(2,3)]) %*% matrix(zCoef, ncol=1)
  for (i in 1:length(yVec)) {
    xReal <- matrix(sqrt(2)*xRaw*z_sigma+z_real_mu[i], ncol = 1)
    muReal <- cbind(1, xReal, matrix(rep(XObs[i,c(2,3)], length(xReal)), ncol = 2, byrow = T)) %*% matrix(beta, ncol = 1)
    log_int[i] <- log(sum(dnorm(yVec[i], mean = muReal, sd = sigma)*wVec)/sqrt(pi))
  }
  loglik2 <- -sum(log_int)
  
  # Compute log-likelihood part1
  muRealObs <- cbind(1, XObs) %*% matrix(beta, ncol = 1)
  loglik1 <- sum(dnorm(yVec, muRealObs, sigma, log = T))
  
  return(-(loglik1+loglik2))
}
