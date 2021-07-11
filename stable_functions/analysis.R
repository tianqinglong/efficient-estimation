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
  lowerConf <- coefVec-2*sdVec
  upperConf <- coefVec+2*sdVec
  
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
    
    compCoef <- c(compLM$coefficients, sigma(compLM))
    names(compCoef) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    compSE <- c(sqrt(diag(vcov(compLM))), NA)
    names(compSE) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    # Missing at random
    datMiss <- dat[dat$Obs == 1,]
    marLM <- lm(fmla, data = datMiss)
    
    marCoef <- c(marLM$coefficients, sigma(marLM))
    names(marCoef) <- c("Intercept", paste("X", 1:num_covariate, sep = ""), "Sigma")
    
    marSE <- c(sqrt(diag(vcov(marLM))), NA)
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
