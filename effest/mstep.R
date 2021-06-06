#-------------------------------
# Use optimization to derive the beta values
# Comments: beta estimates match the optim function, but sigma does not. Need to investigate
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
  
  glm(Y~U+Z, data = dat, weights = WeightVec, family=gaussian(link = "identity")) -> lmFit
  return(lmFit)
}
