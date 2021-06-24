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

ProfileEMOneStep <- function(matObsFull, matMissFull, nsieve)
{
  matAllFull <- as.data.frame(rbind(matObsFull, matMissFull))
  Y2 <- c(rep(1, nrow(matObsFull)), rep(0, nrow(matMissFull)))
  matAllFull <- cbind(matAllFull, Y2)
  
  matAllFullComp <- matAllFull[complete.cases(matAllFull),]
  fmla2 <- as.formula(paste("Y2 ~ 0+", paste(paste("bs", 1:nsieve, sep = ""), collapse = "+")))
  
  f2 <- glm(formula = fmla2, weights = Weight, family = binomial(link = "logit"), data = matAllFullComp)
  tauNew <- f2$coefficients
  
  return(list(Tau = tauNew, LogLik = logLik(f2)))
}

ProfileEM <- function(df_MNAR, beta_init, sigma_init, tau_init,
                      bn = 3, q = 3, gaussHermiteNodes = 10,
                      max_iter = 200, tol = 1e-4)
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
  
  while (iter <= max_iter & !SUCCESS)
  {
    df2 <- MakeFullDataSetMissing(datList, rules, beta_fixed, sigma_fixed, tau_old)
    newList <- ProfileEMOneStep(df1, df2, bn+q)
    
    tau_new <- newList$Tau
    
    dis <- max((tau_new-tau_old)^2)
    
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
  
  n_covariate <- nu+nz
  
  matAllFull <- as.data.frame(rbind(df1, df2))
  Y1 <- log(matAllFull[,"Y"])-log(1-matAllFull[,"Y"])
  matAllFull <- cbind(matAllFull, Y1)
  matAllFullSub <- matAllFull[!(matAllFull$Y1 %in% c(NA, NaN, Inf, -Inf)) & ! is.nan(matAllFull$Weight) ,]
  
  dat_model <- as.matrix(cbind(1, matAllFullSub[, 3:(2+n_covariate)]))
  mu_vec <- dat_model %*% matrix(beta_fixed, ncol = 1)
  
  LogLik1 <- sum( dnorm(matAllFullSub$Y1, mu_vec, sigma_fixed, log = T) * matAllFullSub$Weight )
  
  return(list(LogLik = as.numeric(newList$LogLik+LogLik1)))
}

#-----------------------
# Compute Covariance Matrix Component
#-----------------------

ProfileCov <- function(df_MNAR, n, h_coef, beta_mle, sigma_mle, tau_init,
                       bn, q, gHNodes, max_iter, tol)
{
  numPara <- length(beta_mle)+1
  covMat <- matrix(nrow = numPara, ncol = numPara)
  hn <- h_coef/sqrt(n)
  diagMat <- diag(numPara)
  theta_mle <- c(beta_mle, sigma_mle)
  for (k in 1:numPara)
    for (l in 1:numPara)
    {
      e_k <- diagMat[k,]
      e_l <- diagMat[l,]
      
      theta_1 <- theta_mle+(e_k+e_l)*hn
      theta_2 <- theta_mle+e_k*hn
      theta_3 <- theta_mle+e_l*hn
      
      pf1 <- ProfileEM(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                       bn, q, gHNodes, max_iter, tol)$LogLik
      pf2 <- ProfileEM(df_MNAR, theta_2[1:(numPara-1)], theta_2[numPara], tau_init,
                       bn, q, gHNodes, max_iter, tol)$LogLik
      pf3 <- ProfileEM(df_MNAR, theta_3[1:(numPara-1)], theta_3[numPara], tau_init,
                       bn, q, gHNodes, max_iter, tol)$LogLik
      pf4 <- ProfileEM(df_MNAR, beta_mle, sigma_mle, tau_init,
                       bn, q, gHNodes, max_iter, tol)$LogLik
      covMat[k,l] <- (pf1-pf2-pf3+pf4)/hn^2
    }
  return(solve(-covMat))
}
