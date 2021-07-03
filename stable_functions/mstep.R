#-----------------------
# Reorganize the dataset
#-----------------------

MakeFullDataSetObs <- function(dataList)
# This function is to make a data frame for the M-step for the non-missing data
{
  data <- dataList$data
  dataObs <- data[data[,"Obs"] == 1,]
  
  # For the observed data
  uIndices <- dataList$U_indices
  bn <- dataList$bn
  q <- dataList$q
  
  nObs <- nrow(dataObs)
  nu <- length(uIndices)
  bs_u_list <- paste("bs_u", 1:nu, sep = "")
  
  bsMatObsDat <- matrix(nrow = nObs, ncol = (bn+q)^(1+nu))
  
  for (i in 1:nObs)
  {
    if(nu > 0)
    {
      ubsMat <- matrix(ncol = nu, nrow = bn+q)
      for (j in 1:nu)
      {
        bs_u <- dataList[[bs_u_list[j]]]
        ubsMat[,j] <- bs_u[i,]
      }
    }
    else
    {
      ubsMat <- NULL
    }
    YUSplinePrep(dataList[["bs_y"]][i,], ubsMat) -> bsMatObsDat[i,]
  }
  colnames(bsMatObsDat) <- paste("bs", 1:((bn+q)^(1+nu)), sep = "")
  
  zIndices <- dataList$Z_indices
  if(length(zIndices) > 0)
  {
    ZMat <- as.matrix(dataObs[, zIndices])
    colnames(ZMat) <- paste("Z", 1:length(zIndices), sep = "")
  }
  else
  {
    ZMat <- NULL
  }
  
  if(nu > 0)
  {
    UMat <- as.matrix(dataObs[, uIndices])
    colnames(UMat) <- paste("U", 1:nu, sep = "")
  }
  else
  {
    UMat <- NULL
  }
  
  Y <- dataObs[,"Y"]
  Weight <- 1
  glmMatObsDat <- cbind(Y, Weight, ZMat, UMat, bsMatObsDat)
  
  return(glmMatObsDat)
}

MakeFullDataSetMissing <- function(dataList, rules, betaOld, sigmaOld, tauOld)
{
  data <- dataList$data
  dataMiss <- data[data[,"Obs"] == 0,]
  nMiss <- nrow(dataMiss)
  nObs <- nrow(data)-nMiss
  
  bn <- dataList$bn
  q <- dataList$q
  uIndices <- dataList$U_indices
  zIndices <- dataList$Z_indices
  
  delta <- dataList$delta
  
  nu <- length(uIndices)
  bs_u_list <- paste("bs_u", 1:nu, sep = "")
  
  # For missing data
  x <- rules$x
  w <- rules$w
  
  bsMatMissDat <- matrix(ncol = (bn+q)^(nu+1)+nu+length(zIndices)+2,
                         nrow = nMiss*length(x))
  for (i in 1:nMiss)
  {
    covariates <- c(1, dataMiss[i, c(zIndices, uIndices)])
    mu <- sum(betaOld*covariates)
    yVec <- sqrt(2)*sigmaOld*x+mu
    
    if (is.null(zIndices))
    {
      Z <- NULL
    }
    else
    {
      Z <- dataMiss[i, zIndices]
    }
    
    if (is.null(uIndices))
    {
      U <- NULL
    }
    else
    {
      U <- dataMiss[i, uIndices]
    }
    
    if(nu > 0)
    {
      ubsMat <- matrix(ncol = nu, nrow = bn+q)
      for (j in 1:nu)
      {
        bs_u <- dataList[[bs_u_list[j]]]
        ubsMat[,j] <- bs_u[i+nObs,]
      }
    }
    else
    {
      ubsMat <- NULL
    }
    bsMatMissDat[(1+(i-1)*length(x)):(i*length(x)),] <- PseudoObservation(yVec, w, Z, U, ubsMat,  delta, bn, q, tauOld)
  }
  
  if(is.null(zIndices))
  {
    zNames <- NULL
  }
  else
  {
    zNames <- paste("Z",1:length(zIndices), sep = "")
  }

  if(is.null(uIndices))
  {
    uNames <- NULL
  }
  else
  {
    uNames <- paste("U",1:length(uIndices), sep = "")
  }
  
  colnames(bsMatMissDat) <- c("Y", "Weight", zNames, uNames, paste("bs", 1:(bn+q)^(nu+1), sep = ""))
  return(bsMatMissDat)
}

#-----------------------
# Use glm to do the estimation
#-----------------------

MStep_1Step <- function(matObsFull, matMissFull, nu, nz, nsieve, fixed_beta_sigma = FALSE, beta_old, sigma_old)
{
  # Beta and Sigma Part
  if (nu > 0)
  {
    allU <- paste("U", 1:nu, sep = "")
  }
  else
  {
    allU <- NULL
  }
  
  if (nz > 0)
  {
    allZ <- paste("Z", 1:nz, sep = "")
  }
  else
  {
    allZ <- NULL
  }
    
  matAllFull <- as.data.frame(rbind(matObsFull, matMissFull))
  Y1 <- log(matAllFull[,"Y"])-log(1-matAllFull[,"Y"])
  Y2 <- c(rep(1, nrow(matObsFull)), rep(0, nrow(matMissFull)))
  
  matAllFull %>% mutate(Y1 = log(Y/(1-Y)), Y2=Y2) -> matAllFull
  comIndex <- complete.cases(matAllFull)
  matAllFull <- matAllFull[comIndex,]
  
  if(fixed_beta_sigma)
  {
    betaNew <- beta_old
    sigmaNew <- sigma_old
  }
  else
  {
    fmla1 <- as.formula(paste("Y1 ~ ", paste(c(allZ, allU), collapse = "+")))
    # matAllFullSub <- matAllFull[!(matAllFull$Y1 %in% c(NA, NaN, Inf, -Inf)) & ! is.nan(matAllFull$Weight) ,]
    f1 <- lm(formula = fmla1, weights = Weight, data = matAllFull)
    betaNew <- f1$coefficients
    sigmaNew <- sqrt(sum(matAllFull$Weight*f1$residuals^2, na.rm = T)/sum(matAllFull$Weight, na.rm = T))
  }
  
  # Tau Part
  fmla2 <- as.formula(paste("Y2 ~ 0+", paste(paste("bs", 1:nsieve, sep = ""), collapse = "+")))
  f2 <- glm(formula = fmla2, weights = Weight, family = binomial(link = "logit"), data = matAllFull)
  tauNew <- f2$coefficients
  
  return(list(Beta = betaNew, Sigma = sigmaNew, Tau = tauNew))
}

#-----------------------
# EM Algorithm Iteration Function
#-----------------------
# df_MNAR: is the list without splines. Only with Y, Obs, X and info about U and Z
main <- function(df_MNAR, beta_init, sigma_init, tau_init,
                 bn = 3, q = 3, gaussHermiteNodes = 10,
                 max_iter = 200, tol = 1e-4, fixed_beta_sigma=F)
{
  rules <- gaussHermiteData(gaussHermiteNodes)
  datList <- AppendSplines(df_MNAR, bn, q)
  nu <- length(df_MNAR$U_indices)
  nz <- length(df_MNAR$Z_indices)
  nsieve <- (bn+q)^(nu+1)
  
  iter <- 1
  SUCCESS <- 0
  df1 <- MakeFullDataSetObs(datList)
  beta_old <- beta_init
  sigma_old <- sigma_init
  tau_old <- tau_init
  
  while(iter <= max_iter & !SUCCESS)
  {
    df2 <- MakeFullDataSetMissing(datList, rules, beta_old, sigma_old, tau_old)
    newList <- MStep_1Step(df1, df2, nu, nz, nsieve, fixed_beta_sigma, beta_old, sigma_old)
    
    beta_new <- newList$Beta
    sigma_new <- newList$Sigma
    tau_new <- newList$Tau
    
    dis <- sum((beta_new-beta_old)^2)+(sigma_new-sigma_old)^2+sum((tau_old-tau_new)^2)
    
    if (dis < tol)
    {
      SUCCESS <- 1
    }
    
    if (iter %% 5 == 0 | iter == 1)
    {
      print(paste(paste("iter ", iter, ":", sep = ""), "distance=", round(dis, digits = 6)))
    }
    
    beta_old <- beta_new
    sigma_old <- sigma_new
    tau_old <- tau_new
    
    iter <- iter+1
  }
  
  if (SUCCESS)
  {
    print("EM algorithm converges!")
  }
  else
  {
    print("Failed to converge!")
  }
  
  newList[["Success"]] <- SUCCESS
  return(newList)
}
