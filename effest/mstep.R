#-------------------------------
# Use optimization to derive the beta values
# Comments: If we can use lm(), no need to use the target function in e-step.
# Comments: weighted likelihood estimation is different from weighted least square!
# Comments: need to change, the current version is wrong!
#-------------------------------

betaTargetFuncMax_LM <- function(betaVec, sigmaVec, tauVec, X_obs_mat,
                                 mat_spline, X_missing_mat, Y_non_missing)
{
  n_non_missing <- nrow(X_obs_mat)
  n_missing <- nrow(X_missing_mat)
  
  Y_new <- c(Y_non_missing, rep(Y_non_missing, n_missing))
  X_new <- as.matrix(X_obs_mat)
  for (i in 1:n_missing)
  {
    X_inst <- X_missing_mat[i,]
    X_add <- matrix(rep(X_inst, n_non_missing), ncol = length(X_inst), byrow = TRUE)
    X_new <- rbind(X_new, X_add)
  }
  
  X_new <- matrix(unlist(X_new), ncol = 3)
  
  WeightVec <- rep(1, n_non_missing)
  conditionalExpectionVec_OnlyY(mat_spline, X_missing_mat, Y_non_missing, tauVec,
                                betaVec, sigmaVec) -> condProbMat
  WeightVec <- c(WeightVec, c(t(condProbMat)))
  
  dat <- cbind(Y_new, X_new)
  colnames(dat) <- c("Y", "Intercept", "U", "Z")
  dat <- as.data.frame(dat)
  
  lm(Y~U+Z, data = dat, weights = WeightVec) -> lmFit
  return(lmFit)
}
