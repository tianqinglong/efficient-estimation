#-----------------------
# Function to do the analysis of the output
# Author: Qinglong Tian
# Date: July 10, 2021
#-----------------------

NegLogLikLinearRegression <- function(theta, data)
{
  beta <- theta[-length(theta)]
  sigma <- theta[length(theta)]
  dat <- data[,-2]
  
  yDat <- dat[,1]
  xDat <- cbind(1, dat[,-1])
  xMu <- as.matrix(xDat) %*% matrix(beta, ncol = 1)
  
  -sum(dnorm(log(yDat/(1-yDat)), mean = xMu, sd = sigma, log = T))
}

analysis <- function(rout, true_theta)
{
  ntotal <- length(rout)
  rout_save <- rout
  rout[sapply(rout, function(x) {is.character(x$Var)})] <- NULL
  nclean <- length(rout)
  
  # Proportion of successful trials
  success_rate <- nclean/ntotal
  
  ## EM
  
  # Coverage probability
  sdVec <- t(sapply(rout, function(x) {sqrt(diag(x$Var))}))
  coefVec <- t(sapply(rout, function(x) {c(x$EM$Beta, x$EM$Sigma)}))
  lowerConf <- coefVec-1.96*sdVec
  upperConf <- coefVec+1.96*sdVec
  
  trueVal <- matrix(rep(true_theta, times = length(rout)), byrow = T, ncol = length(true_theta))
  cpVec <- colMeans((trueVal <= upperConf) & (trueVal >= lowerConf))
  num_covariate <- length(cpVec)-2
  names(cpVec) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
  
  # Empirical Bias
  ebVec <- colMeans(coefVec-trueVal)
  names(ebVec) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
  
  # Mean of SE
  seVec <- colMeans(sdVec)
  names(seVec) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
  
  # Empirical standard deviation
  esdVec <- apply(coefVec, MARGIN = 2, sd)
  names(esdVec) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
  
  emOut <- rbind(cpVec, ebVec, seVec, esdVec)
  rownames(emOut) <- c("EM_CP", "EM_Bias", "EM_SE", "EM_SD")
  
  ## MAR & Complete
  
  lapply(rout, function(x){
    dat <- as.data.frame(x$Data$data)
    fmla <- formula(paste("log(Y/(1-Y)) ~ ", paste("X", 1:num_covariate, sep = "", collapse = "+"), sep = ""))
    # Complete
    compLM <- lm(fmla, data = dat)
    
    compopt <- optim(c(compLM$coefficients, sigma(compLM)), NegLogLikLinearRegression,
                     data = dat, hessian = T)
    compCoef <- compopt$par
    names(compCoef) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    compSE <- sqrt(diag(solve(compopt$hessian)))
    names(compSE) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    # Missing at random
    datMiss <- dat[dat$Obs == 1,]
    marLM <- lm(fmla, data = datMiss)
    
    maropt <- optim(c(marLM$coefficients, sigma(marLM)), NegLogLikLinearRegression,
                     data = datMiss, hessian = T)
    sigma_se_mar <- sqrt(diag(compopt$hessian))[length(maropt$coefficients)]
    
    marCoef <- maropt$par
    
    names(marCoef) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    marSE <- sqrt(diag(solve(maropt$hessian)))
    names(marSE) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    return(list(
      BD_coef = compCoef,
      BD_se = compSE,
      MAR_coef = marCoef,
      MAR_se = marSE
    ))
  }
  ) -> mar_comp
  
  # Complete
  BD_coef_mat <- t(sapply(mar_comp, function(x){x$BD_coef}))
  BD_se_mat <- t(sapply(mar_comp, function(x){x$BD_se}))
  BD_lowerConf <- BD_coef_mat-1.96*BD_se_mat
  BD_upperConf <- BD_coef_mat+1.96*BD_se_mat
  BD_cpVec <- (trueVal <= BD_upperConf) & (trueVal >= BD_lowerConf)
  
  BD_CP <- colMeans(BD_cpVec)
  BD_Bias <- colMeans(BD_coef_mat-trueVal)
  BD_SE <- colMeans(BD_se_mat)
  BD_SD <- apply(BD_coef_mat, MARGIN = 2, sd)
  
  bdOut <- rbind(BD_CP, BD_Bias, BD_SE, BD_SD)
  
  # MAR
  MAR_coef_mat <- t(sapply(mar_comp, function(x){x$MAR_coef}))
  MAR_se_mat <- t(sapply(mar_comp, function(x){x$MAR_se}))
  MAR_lowerConf <- MAR_coef_mat-1.96*MAR_se_mat
  MAR_upperConf <- MAR_coef_mat+1.96*MAR_se_mat
  MAR_cpVec <- (trueVal <= MAR_upperConf) & (trueVal >= MAR_lowerConf)
  
  MAR_CP <- colMeans(MAR_cpVec)
  MAR_Bias <- colMeans(MAR_coef_mat-trueVal)
  MAR_SE <- colMeans(MAR_se_mat)
  MAR_SD <- apply(MAR_coef_mat, MARGIN = 2, sd)
  
  marOut <- rbind(MAR_CP, MAR_Bias, MAR_SE, MAR_SD)
  
  ## Pseudo-likelihood
  
  lapply(rout, function(x)
  {
    pseudo <- x$Pseudo
    coefPseudo <- pseudo$par
    sePseudo <- sqrt(diag(solve(pseudo$hessian)))
    lowerPseudo <- coefPseudo-1.96*sePseudo
    upperPseudo <- coefPseudo+1.96*sePseudo
    coverPseudo <- ( (true_theta <= upperPseudo) & (true_theta >= lowerPseudo))
    
    return(list(
      coef = coefPseudo,
      se = sePseudo,
      lower = lowerPseudo,
      upper = upperPseudo,
      cover = coverPseudo
    ))
  }
  ) -> pseudo_out
  
  Pseudo_CP <- rowMeans(sapply(pseudo_out, function(x) {
    x$cover
  })
  )
  
  Pseudo_Bias <- rowMeans(sapply(pseudo_out, function(x) {
    x$coef-true_theta
  })
  )
  
  Pseudo_SE <- rowMeans(sapply(pseudo_out, function(x) {
    x$se
  })
  )
  
  Pseudo_SD <- apply(sapply(pseudo_out, function(x) {x$coef}), MARGIN = 1, sd)
  
  pseudoOut <- rbind(Pseudo_CP, Pseudo_Bias, Pseudo_SE, Pseudo_SD)
  
  ## Diagnosis
  rout_save[sapply(rout_save, function(x) {!is.character(x$Var)})] <- NULL
  if(length(rout_save) == 0)
  {
    Diagnosis1 <- "No Failure."
  }
  else
  {
    Diagnosis1 <- table(sapply(rout_save, function(x) {x$Var}))
  }
  
  Diagnosis2 <- c(nclean, ntotal, success_rate)
  names(Diagnosis2) <- c("Successfull Trials", "Total Trials", "Success Rate")
  
  return(list(emOut, pseudoOut, bdOut, marOut, Diagnosis1, Diagnosis2))
}
