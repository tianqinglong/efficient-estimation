computeLogLikelihood <- function(df_MNAR, beta, sd, tau, ghnodes, nsieves, bn, q)
{
  dataList <- AppendSplines(df_MNAR, bn, q)
  nz <- length(dataList$Z_indices)
  nu <- length(dataList$U_indices)
  rules <- gaussHermiteData(ghnodes)
  
  dataObsYUSpline <- MakeFullDataSetObs(dataList)
  dataMissYUSpline <- MakeFullDataSetMissing(dataList, rules, beta, sd, tau)
  datYUSpline <- rbind(dataObsYUSpline, dataMissYUSpline)
  
  # datCompIndex <- complete.cases(datYUSpline)
  datYUSpline <- as.data.frame(datYUSpline)
  
  fake01 <- c(rep(1, nrow(dataObsYUSpline)),
              rep(0, nrow(dataMissYUSpline)))
  # fake01 <- fake01[datCompIndex]
  
  datYUSpline %>% mutate(Y1 = log(Y/(1-Y)), Y2 = fake01) -> datYUSpline
  
  # Summing up log-likelihood
  if(nz == 0)
  {
    ZCol <- NULL
  }
  else
  {
    ZCol <- paste("Z", 1:nz, sep = "")
  }
  
  if(nu == 0)
  {
    UCol <- NULL
  }
  else
  {
    UCol <- paste("U", 1:nu, sep = "")
  }
  
  covariateCol <- c(ZCol, UCol)
  
  rr <- datYUSpline[, "Y2"]
  yy <- datYUSpline[, "Y1"]
  xx <- datYUSpline %>% select(all_of(covariateCol))
  bb <- datYUSpline[, paste("bs", 1:nsieves, sep = "")]
  mu <- as.matrix(cbind(1,xx)) %*% as.matrix(beta, ncol=1)
  ww <- datYUSpline[,"Weight"]
  tauB <- as.matrix(bb) %*% as.matrix(tau, ncol=1)
  
  # logliksum <- rr*ww*dnorm(yy, mu, sd, log = T)+(1-rr)*ww*dnorm(yy, mu, sd, log = T)+
  #   rr*ww*(tauB-log(1+exp(tauB)))-(1-rr)*ww*log(1+exp(tauB))
  
  logliksum1 <- sum(rr*dnorm(yy, mu, sd, log = T)+rr*tauB-rr*log(1+exp(tauB)), na.rm = T)
  
  logliksum2 <- 0
  nRealMiss <- sum(!dataList$data[,"Obs"])
  datMiss <- datYUSpline %>% filter(Y2 == 0)
  w <- rules$w
  for (i in 1:nRealMiss)
  {
    startIndex <- 1+(i-1)*ghnodes
    endIndex <- i*ghnodes
    
    bsMat <- datMiss[startIndex:endIndex, paste("bs", 1:nsieves, sep = "")]
    tauB <- as.matrix(bsMat) %*% matrix(tau, ncol = 1)
    sigmoidTauB <- 1/(1+exp(-tauB))
    integral <- sum(sigmoidTauB*w)/sqrt(pi)
    
    logliksum2 <- logliksum2+log(1-integral)
  }
  
  return(logliksum1+logliksum2)
}

computeLogLikelihood_additive <- function(df_MNAR, beta, sd, tau, ghnodes, nsieves, bn, q)
{
  dataList <- AppendSplines(df_MNAR, bn, q)
  nz <- length(dataList$Z_indices)
  nu <- length(dataList$U_indices)
  rules <- gaussHermiteData(ghnodes)
  
  dataObsYUSpline <- MakeFullDataSetObs_additive(dataList)
  dataMissYUSpline <- MakeFullDataSetMissing_additive(dataList, rules, beta, sd, tau)
  datYUSpline <- rbind(dataObsYUSpline, dataMissYUSpline)
  
  # datCompIndex <- complete.cases(datYUSpline)
  datYUSpline <- as.data.frame(datYUSpline)
  
  fake01 <- c(rep(1, nrow(dataObsYUSpline)),
              rep(0, nrow(dataMissYUSpline)))
  # fake01 <- fake01[datCompIndex]
  
  datYUSpline %>% mutate(Y1 = log(Y/(1-Y)), Y2 = fake01) -> datYUSpline
  
  # Summing up log-likelihood
  if(nz == 0)
  {
    ZCol <- numeric(0)
  }
  else
  {
    ZCol <- paste("Z", 1:nz, sep = "")
  }
  
  if(nu == 0)
  {
    UCol <- numeric(0)
  }
  else
  {
    UCol <- paste("U", 1:nu, sep = "")
  }
  
  covariateCol <- c(ZCol, UCol)
  
  rr <- datYUSpline[, "Y2"]
  yy <- datYUSpline[, "Y1"]
  xx <- datYUSpline %>% select(all_of(covariateCol))
  bb <- datYUSpline[, paste("bs", 1:nsieves, sep = "")]
  mu <- as.matrix(cbind(1,xx)) %*% as.matrix(beta, ncol=1)
  ww <- datYUSpline[,"Weight"]
  tauB <- as.matrix(bb) %*% as.matrix(tau, ncol=1)
  
  # logliksum <- rr*ww*dnorm(yy, mu, sd, log = T)+(1-rr)*ww*dnorm(yy, mu, sd, log = T)+
  #   rr*ww*(tauB-log(1+exp(tauB)))-(1-rr)*ww*log(1+exp(tauB))
  
  logliksum1 <- sum(rr*dnorm(yy, mu, sd, log = T)+rr*tauB-rr*log(1+exp(tauB)), na.rm = T)
  
  logliksum2 <- 0
  nRealMiss <- sum(!dataList$data[,"Obs"])
  datMiss <- datYUSpline %>% filter(Y2 == 0)
  w <- rules$w
  for (i in 1:nRealMiss)
  {
    startIndex <- 1+(i-1)*ghnodes
    endIndex <- i*ghnodes
    
    bsMat <- datMiss[startIndex:endIndex, paste("bs", 1:nsieves, sep = "")]
    tauB <- as.matrix(bsMat) %*% matrix(tau, ncol = 1)
    sigmoidTauB <- 1/(1+exp(-tauB))
    integral <- sum(sigmoidTauB*w)/sqrt(pi)
    
    logliksum2 <- logliksum2+log(1-integral)
  }
  
  return(logliksum1+logliksum2)
}

ProfileCov <- function(df_MNAR, hn, beta_mle, sigma_mle, tau_mle,
                        bn, q, gHNodes, n_sieve, n_covariate, max_iter, tol)
{
  numPara <- length(beta_mle)+1
  covMat <- matrix(nrow = numPara, ncol = numPara)
  diagMat <- diag(numPara)
  theta_mle <- c(beta_mle, sigma_mle)
  tau_init <- tau_mle
  combn_mat <- cbind(combn(1:numPara, 2), matrix(rep(1:numPara, each = 2), nrow = 2))
  
  plusplus <- numeric(ncol(combn_mat))
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    theta_1 <- theta_mle+(diagMat[k,]+diagMat[l,])*hn/2
    theta_1[numPara] <- max(theta_1[numPara], hn/2)
    emTemp <- main(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                   bn, q, gHNodes, max_iter, tol, T)
    plusplus[i] <- computeLogLikelihood(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], emTemp$Tau,
                                        gHNodes, n_sieve, bn, q)
  }
  
  minusplus <- numeric(ncol(combn_mat))
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    theta_1 <- theta_mle+(-diagMat[k,]+diagMat[l,])*hn/2
    theta_1[numPara] <- max(theta_1[numPara], hn/2)
    emTemp <- main(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                   bn, q, gHNodes, max_iter, tol, T)
    minusplus[i] <- computeLogLikelihood(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], emTemp$Tau,
                                         gHNodes, n_sieve, bn, q)
  }
  
  plusminus <- numeric(ncol(combn_mat))
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    theta_1 <- theta_mle+(diagMat[k,]-diagMat[l,])*hn/2
    theta_1[numPara] <- max(theta_1[numPara], hn/2)
    emTemp <- main(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                   bn, q, gHNodes, max_iter, tol, T)
    plusminus[i] <- computeLogLikelihood(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], emTemp$Tau,
                                         gHNodes, n_sieve, bn, q)
  }
  
  minusminus <- numeric(ncol(combn_mat))
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    theta_1 <- theta_mle+(-diagMat[k,]-diagMat[l,])*hn/2
    theta_1[length(theta_1)] <- max(theta_1[length(theta_1)], 0.01)
    emTemp <- main(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                   bn, q, gHNodes, max_iter, tol, T)
    minusminus[i] <- computeLogLikelihood(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], emTemp$Tau,
                                         gHNodes, n_sieve, bn, q)
  }
  
  diffs <- plusplus-minusplus-plusminus+minusminus
  
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    
    val <- diffs[i]
    covMat[k,l] <- val
    covMat[l,k] <- val
  }
  
  cov_mat <- cov_mat/hn^2
  tryCatch(
    expr = {
      cov_mat <- solve(-covMat)
    },
    error = function(e)
    {
      cov_mat <- solve(-covMat+diag(hn, nrow = numPara))
    }
  )
  
  return(cov_mat)
}

ProfileCov_additive <- function(df_MNAR, hn, beta_mle, sigma_mle, tau_mle,
                                bn, q, gHNodes, n_sieve, n_covariate, max_iter, tol)
{
  numPara <- length(beta_mle)+1
  covMat <- matrix(nrow = numPara, ncol = numPara)
  diagMat <- diag(numPara)
  theta_mle <- c(beta_mle, sigma_mle)
  tau_init <- tau_mle
  combn_mat <- cbind(combn(1:numPara, 2), matrix(rep(1:numPara, each = 2), nrow = 2))
  
  plusplus <- numeric(ncol(combn_mat))
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    theta_1 <- theta_mle+(diagMat[k,]+diagMat[l,])*hn/2
    theta_1[numPara] <- max(theta_1[numPara], hn/2)
    emTemp <- main_additive(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                   bn, q, gHNodes, max_iter, tol, T)
    plusplus[i] <- computeLogLikelihood_additive(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], emTemp$Tau,
                                        gHNodes, n_sieve, bn, q)
  }
  
  minusplus <- numeric(ncol(combn_mat))
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    theta_1 <- theta_mle+(-diagMat[k,]+diagMat[l,])*hn/2
    theta_1[numPara] <- max(theta_1[numPara], hn/2)
    emTemp <- main_additive(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                   bn, q, gHNodes, max_iter, tol, T)
    minusplus[i] <- computeLogLikelihood_additive(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], emTemp$Tau,
                                         gHNodes, n_sieve, bn, q)
  }
  
  plusminus <- numeric(ncol(combn_mat))
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    theta_1 <- theta_mle+(diagMat[k,]-diagMat[l,])*hn/2
    theta_1[numPara] <- max(theta_1[numPara], hn/2)
    emTemp <- main_additive(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                   bn, q, gHNodes, max_iter, tol, T)
    plusminus[i] <- computeLogLikelihood_additive(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], emTemp$Tau,
                                         gHNodes, n_sieve, bn, q)
  }
  
  minusminus <- numeric(ncol(combn_mat))
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    theta_1 <- theta_mle+(-diagMat[k,]-diagMat[l,])*hn/2
    theta_1[length(theta_1)] <- max(theta_1[length(theta_1)], 0.01)
    emTemp <- main_additive(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], tau_init,
                   bn, q, gHNodes, max_iter, tol, T)
    minusminus[i] <- computeLogLikelihood_additive(df_MNAR, theta_1[1:(numPara-1)], theta_1[numPara], emTemp$Tau,
                                         gHNodes, n_sieve, bn, q)
  }
  
  diffs <- plusplus-minusplus-plusminus+minusminus
  
  for(i in 1:ncol(combn_mat))
  {
    chosen <- combn_mat[,i]
    k <- chosen[1]
    l <- chosen[2]
    
    val <- diffs[i]
    covMat[k,l] <- val
    covMat[l,k] <- val
  }
  
  covMat <- covMat/hn^2
  
  tryCatch(
    expr = {
      cov_mat <- solve(-covMat)
    },
    error = function(e)
    {
      cov_mat <- solve(-covMat+diag(hn, nrow = numPara))
    }
  )
  
  return(cov_mat)
}
