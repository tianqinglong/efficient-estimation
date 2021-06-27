#-----------------------
# Use profile likelihood to obtain the covariance matrix for \beta and \sigma
# Author: Qinglong Tian
# Date: June 23, 2021
#-----------------------

#-----------------------
# EM algorithm to maximize the tau
# A way to do is not to update beta and sigma, which is simple to implement but not efficient.
# Could also focus on the second part of the log-likelihood.
#-----------------------

LogLikFixedBetaSigma <- function(tauVec, matObsFull, matMissFull, nsieve, fixedBeta, fixedSigma, n_covariate)
{
  
  n_obs <- nrow(matObsFull)
  n_miss <- nrow(matMissFull)
  obsVec <- c(rep(1, n_obs), rep(0, n_miss))
  
  Y2 <- c(rep(1, nrow(matObsFull)), rep(0, nrow(matMissFull)))
  
  matAllFull <- as.data.frame(rbind(matObsFull, matMissFull))
  Y1 <- log(matAllFull[,"Y"])-log(1-matAllFull[,"Y"])
  matAllFull <- cbind(matAllFull, Y2, Y1)
  
  matAllComp <- matAllFull[complete.cases(matAllFull),]
  
  # Linear part
  dat_model <- as.matrix(cbind(1, matAllComp[, 3:(2+n_covariate)]))
  mu_vec <- dat_model %*% matrix(fixedBeta, ncol = 1)
  logLikLinear <- sum(dnorm(matAllComp$Y1, mu_vec, fixedSigma, log = T)*matAllComp$Weight)
  
  # Logistic part
  bspline_model <- as.matrix(matAllComp[,(3+n_covariate):(2+n_covariate+nsieve)])
  logit_pi <- bspline_model %*% tauVec 
  for_non_missing <- logit_pi-log(1+exp(logit_pi))
  for_missing <- -log(1+exp(logit_pi))
  
  output <- sum(obsVec*for_non_missing+(1-obsVec)*for_missing)+logLikLinear
  
  return(-output)
}

ProfileEMOneStep <- function(tau_init, matObsFull, matMissFull, nsieve, fixedBeta, fixedSigma, n_covariate)
{
  opt <- optim(tau_init, LogLikFixedBetaSigma, matObsFull = matObsFull,
               matMissFull = matMissFull, nsieve = nsieve,
               fixedBeta = fixedBeta, fixedSigma = fixedSigma, n_covariate = n_covariate)
  return(opt)
}

ProfileEM <- function(df_MNAR, beta_init, sigma_init, tau_init,
                      bn, q, gaussHermiteNodes, nsieve, n_covariate,
                      max_iter = 200, tol = 1e-4, evaluate_only = F)
{
  rules <- gaussHermiteData(gaussHermiteNodes)
  datList <- AppendSplines(df_MNAR, bn, q)
  nu <- length(df_MNAR$U_indices)
  nz <- length(df_MNAR$Z_indices)
  
  iter <- 1
  SUCCESS <- 0
  df1 <- MakeFullDataSetObs(datList)
  beta_fixed <- beta_init
  sigma_fixed<- sigma_init
  tau_old <- tau_init
  
  if (evaluate_only)
  {
    df2 <- MakeFullDataSetMissing(datList, rules, beta_fixed, sigma_fixed, tau_old)
    out <- LogLikFixedBetaSigma(tauVec = tau_init, matObsFull = df1, matMissFull = df2,
                                nsieve = nsieve, fixedBeta = beta_init, fixedSigma = sigma_init,
                                n_covariate = n_covariate)
    return(-out)
  }
  
  while (iter <= max_iter & !SUCCESS)
  {
    df2 <- MakeFullDataSetMissing(datList, rules, beta_fixed, sigma_fixed, tau_old)
    opt <- ProfileEMOneStep(tau_old, df1, df2, nsieve, beta_fixed, sigma_fixed, n_covariate)

    tau_new <- opt$par
    dis <- sum((tau_new-tau_old)^2)
    
    if (dis < tol)
    {
      SUCCESS <- 1
    }
    
    # if (iter %% 5  == 0 | iter == 1)
    # {
    #   print(paste(paste("iter ", iter, ":", sep = ""), "distance=", round(dis, digits = 6)))
    # }
    
    tau_old <- tau_new
    iter <- iter+1
  }
  
  # if (SUCCESS)
  # {
  #   print("EM algorithm converged!")
  # }
  # else
  # {
  #   print("EM algorithm failed to converge!")
  # }
  
  # n_covariate <- nu+nz
  # 
  # matAllFull <- as.data.frame(rbind(df1, df2))
  # Y1 <- log(matAllFull[,"Y"])-log(1-matAllFull[,"Y"])
  # matAllFull <- cbind(matAllFull, Y1)
  # matAllFullSub <- matAllFull[!(matAllFull$Y1 %in% c(NA, NaN, Inf, -Inf)) & ! is.nan(matAllFull$Weight) ,]
  # 
  # dat_model <- as.matrix(cbind(1, matAllFullSub[, 3:(2+n_covariate)]))
  # mu_vec <- dat_model %*% matrix(beta_fixed, ncol = 1)
  # 
  # LogLik1 <- sum( dnorm(matAllFullSub$Y1, mu_vec, sigma_fixed, log = T) * matAllFullSub$Weight )
  
  return(-opt$value)
  
#  return(list(LogLik = as.numeric(newList$LogLik+LogLik1)))
}

#-----------------------
# Compute the log-likelihood
#-----------------------


#-----------------------
# Compute Covariance Matrix Component
#-----------------------

ProfileCov <- function(df_MNAR, hn, beta_mle, sigma_mle, tau_init,
                       bn, q, gHNodes, n_sieve, n_covariate, max_iter, tol)
{
  numPara <- length(beta_mle)+1
  covMat <- matrix(nrow = numPara, ncol = numPara)
  # hn <- h_coef/sqrt(n)
  diagMat <- diag(numPara)
  theta_mle <- c(beta_mle, sigma_mle)
  pf4 <- ProfileEM(df_MNAR, beta_mle, sigma_mle, tau_init,
                   bn, q, gHNodes, n_sieve, n_covariate, max_iter, tol, evaluate_only = T)
  pf1Vec <- numeric(numPara)
  for (i in 1:numPara)
  {
    theta_1 <- theta_mle+(diagMat[i,])*hn
    pf1Vec[i] <- ProfileEM(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                           bn, q, gHNodes, n_sieve, n_covariate, max_iter, tol)
  }
  
  pf2Mat <- numeric(choose(numPara,2))
  combn_mat <- combn(1:numPara, 2)
  
  for (i in 1:length(pf2Mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    theta_1 <- theta_mle+(diagMat[k,])*hn+(diagMat[l,])*hn
    pf2Mat[i] <- ProfileEM(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                           bn, q, gHNodes, n_sieve, n_covariate, max_iter, tol)
  }
  
  pf2Diag <- numeric(numPara)
  for(i in 1:numPara)
  {
    theta_1 <- theta_mle+(diagMat[i,])*2*hn
    pf2Diag[i] <- ProfileEM(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                            bn, q, gHNodes, n_sieve, n_covariate, max_iter, tol)
  }
  
  for (i in 1:length(pf2Mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    
    covMat[k, l] <- (pf2Mat[i]-pf1Vec[k]-pf1Vec[l]+pf4)
    covMat[l, k] <- covMat[k, l]
  }
  
  for (i in 1:numPara)
  {
    covMat[i,i] <- (pf2Diag[i]-2*pf1Vec[i]+pf4)
  }
  
  covMat <- covMat/hn^2
  return(-solve(covMat))
}
