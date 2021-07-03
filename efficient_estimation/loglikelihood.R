computeLogLikelihood <- function(df_MNAR, beta, sd, tau, ghnodes, nsieves, bn, q)
{
  
  dataList <- AppendSplines(df_MNAR, bn, q)
  nz <- length(dataList$Z_indices)
  nu <- length(dataList$U_indices)
  rules <- gaussHermiteData(ghnodes)
  
  dataObsYUSpline <- MakeFullDataSetObs(dataList)
  dataMissYUSpline <- MakeFullDataSetMissing(dataList, rules, beta, sd, tau)
  datYUSpline <- rbind(dataObsYUSpline, dataMissYUSpline)
  
  datCompIndex <- complete.cases(datYUSpline)
  datYUSpline <- as.data.frame(datYUSpline[datCompIndex,])
  
  fake01 <- c(rep(1, nrow(dataObsYUSpline)),
              rep(0, nrow(dataMissYUSpline)))
  fake01 <- fake01[datCompIndex]
  
  datYUSpline %>% mutate(Y1 = log(Y/(1-Y)), Y2 = fake01) -> datYUSpline
  
  # Summing up log-likelihood
  logLikSum <- 0
  for (i in 1:nrow(datYUSpline))
  {
    r_i <- datYUSpline[i, "Y2"]
    y_i <- datYUSpline[i, "Y1"]
    
    ZCol <- ifelse(nz == 0, numeric(0), paste("Z", 1:nz, sep = ""))
    UCol <- ifelse(nu == 0, numeric(0), paste("U", 1:nu, sep = ""))
    
    if (is.na(ZCol))
    {
      covariateCol <- UCol
    }
    else if (is.na(UCol))
    {
      covariateCol <- ZCol
    }
    else
    {
      covariateCol <- c(ZCol, UCol)
    }
    
    x_i <- datYUSpline[i, covariateCol]
    b_i <- datYUSpline[i, paste("bs", 1:nsieves, sep = "")]
    mu_i <- sum(c(1, x_i)*beta)
    w_i <- datYUSpline[i, "Weight"]
    tauB_i <- sum(tau*b_i)
    
    logLikSum <- logLikSum+r_i*dnorm(y_i, mu_i, sd, log = T)+
      (1-r_i)*w_i*dnorm(y_i, mu_i, sd, log = T)+
      r_i*(tauB_i-log(1+exp(tauB_i)))-
      (1-r_i)*w_i*log(1+exp(tauB_i))
  }
  return(logLikSum)
}