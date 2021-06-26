PseudoObs_additive <- function(y, w, Z, U, ubsMat, delta, bn, q, tauOld)
{
  yBSMat <- YSplinePrep(y, delta, bn, q)
  de <- 0
  realWeight <- numeric(nrow(yBSMat))
  B_yu <- matrix(nrow = length(y), ncol = (bn+q)*(ifelse(is.null(ubsMat), 0, ncol(ubsMat))+1))
  for (i in 1:nrow(yBSMat))
  {
    ybs1 <- yBSMat[i,]
    yuBs <- c(ybs1, c(ubsMat))
    B_yu[i,] <- yuBs
    
    exp_tauB_inv <- 1-1/(exp(-sum(yuBs*tauOld))+1)
    realWeight[i] <- exp_tauB_inv*w[i]/sqrt(pi)
    de <- de+w[i]*exp_tauB_inv/sqrt(pi)
  }
  Weight <- realWeight/de
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
