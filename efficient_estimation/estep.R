#-----------------------
# Use Gauss-Hermite to Expand the Expectation
# Author: Qinglong Tian
# Date: June 18, 2021
#-----------------------
# library(fastGHQuad)
# gaussHermiteData(10)
#-----------------------

# Input: data list, number of nodes, old beta (sigma), new beta (sigma), old tau
# Output: the quadrature weight times the weight inside integrand, for each x_i (missing obs)

YSplinePrep <- function(y, delta, bn, q)
# y: nodes from the Gauss-Hermite rules, without inv-logit transformation, with normal transformation
# xVec: contains intercept
{
  y_inv_logit <- exp(y)/(exp(y)+1)
  
  ybs <- matrix(nrow = length(y_inv_logit), ncol = bn+q)
  for (i in 1:length(y_inv_logit))
  {
    BSplinesForNewY(y_inv_logit[i], q, bn, delta) -> ybs[i,]
  }
  
  colnames(ybs) <- paste("bs", 1:(bn+q), sep = "")
  return(ybs)
}

YUSplinePrep <- function(ybs1, ubsMat)
# ubsMat: row is b-spline, col is different u
# Combine the Y spline and U spline
{
  if (is.null(ubsMat))
  {
    return(ybs1)
  }

  yuBs <- numeric(nrow(ubsMat)^(ncol(ubsMat)+1))
  combs <- expand.grid(rep(list(1:nrow(ubsMat)), ncol(ubsMat)+1))
  
  for (i in 1:length(yuBs))
  {
    xProd <- 1
    for (j in 1:ncol(ubsMat))
    {
      xProd <- xProd*ubsMat[combs[i,j+1],j]
    }
    
    yuBs[i] <- ybs1[combs[i,1]]*xProd
  }
  
  return(yuBs)
}

PseudoObservation <- function(y, w, Z, U, ubsMat, delta, bn, q, tauOld)
# Contains three components: y, and weights. Can be further used as new observations
# y: nodes from Gauss-Hermite on scale (-\inf,\inf), after transformation.
{
  YSplinePrep(y, delta, bn, q) -> yBSMat
  de <- 0
  realWeight <- numeric(nrow(yBSMat))
  B_yu <- matrix(nrow = length(y), ncol = (bn+q)^(ifelse(is.null(ubsMat), 0, ncol(ubsMat))+1))
  for (i in 1:nrow(yBSMat))
  {
    ybs1 <- yBSMat[i,]
    yuBs <- YUSplinePrep(ybs1, ubsMat)
    B_yu[i,] <- yuBs
    
    exp_tauB_inv <- 1-1/(exp(-sum(yuBs*tauOld))+1)
    realWeight[i] <- exp_tauB_inv*w[i]/sqrt(pi)
    
    de <- de+w[i]*exp_tauB_inv/sqrt(pi)
  }
  Weight <- realWeight/de
  # if(sum(is.nan(Weight))>0)
  # {
  #   browser()
  # }
  Y <- 1/(1+exp(-y))
  colnames(B_yu) <- paste("bs", 1:ncol(B_yu), sep = "")
  
  if (is.null(Z))
  {
    ZMat <- NULL
  }
  else
  {
    ZMat <- matrix(rep(Z, length(y)), ncol = length(Z), byrow = T)
    colnames(ZMat) <- paste("Z", 1:length(Z), sep = "")
  }
  
  if (is.null(U))
  {
    UMat <- NULL
  }
  else
  {
    UMat <- matrix(rep(U, length(y)), ncol = length(U), byrow = T)
    colnames(UMat) <- paste("U", 1:length(U), sep = "")
  }
  
  return(cbind(Y, Weight, ZMat, UMat, B_yu))
}
