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
  
  colnames(bsMatMissDat) <- c("Y", "Weights", zNames, uNames, paste("bs", 1:(bn+q)^(nu+1), sep = ""))
  return(bsMatMissDat)
}
