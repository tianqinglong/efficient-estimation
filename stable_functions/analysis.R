#-----------------------
# Function to do the analysis of the output
# Author: Qinglong Tian
# Date: July 10, 2021
#-----------------------
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
    
    sigma_lower_comp <- sqrt(sum(compLM$residuals^2)/qchisq(0.975, df = compLM$df.residual))
    sigma_upper_comp <- sqrt(sum(compLM$residuals^2)/qchisq(0.025, df = compLM$df.residual))
    
    sigma_se_comp <- (sigma_upper_comp-sigma_lower_comp)/2/qnorm(0.975)
    
    compCoef <- c(compLM$coefficients, sigma(compLM))
    names(compCoef) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    compSE <- c(sqrt(diag(vcov(compLM))), sigma_se_comp)
    names(compSE) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    # Missing at random
    datMiss <- dat[dat$Obs == 1,]
    marLM <- lm(fmla, data = datMiss)
    
    sigma_lower_mar <- sqrt(sum(marLM$residuals^2)/qchisq(0.975, df = marLM$df.residual))
    sigma_upper_mar <- sqrt(sum(marLM$residuals^2)/qchisq(0.025, df = marLM$df.residual))
    
    sigma_se_mar <- (sigma_upper_mar-sigma_lower_mar)/2/qnorm(0.975)
    
    marCoef <- c(marLM$coefficients, sigma(marLM))
    names(marCoef) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    marSE <- c(sqrt(diag(vcov(marLM))), sigma_se_mar)
    names(marSE) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    return(list(
      BD_coef = compCoef,
      BD_se = compSE,
      BD_sigma_lower = sigma_lower_comp,
      BD_sigma_upper = sigma_upper_comp,
      MAR_coef = marCoef,
      MAR_se = marSE,
      MAR_sigma_lower = sigma_lower_mar,
      MAR_sigma_upper = sigma_upper_mar
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
  
  BD_sigma_cp <- t(sapply(mar_comp, function(x)
    {
    (x$BD_coef[length(x$BD_coef)] >= x$BD_sigma_lower) & (x$BD_coef[length(x$BD_coef)] <= x$BD_sigma_upper)
    }
    )
  )
  BD_CP[length(BD_CP)] <- mean(BD_sigma_cp)
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
  
  MAR_sigma_cp <- t(sapply(mar_comp, function(x)
    {
    (x$MAR_coef[length(x$MAR_coef)] >= x$MAR_sigma_lower) & (x$MAR_coef[length(x$MAR_coef)] <= x$BD_sigma_upper)
    }
    )
  )
  MAR_CP[length(MAR_CP)] <- mean(MAR_sigma_cp)
  marOut <- rbind(MAR_CP, MAR_Bias, MAR_SE, MAR_SD)
  
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
  
  return(list(emOut, bdOut, marOut, Diagnosis1, Diagnosis2))
}
