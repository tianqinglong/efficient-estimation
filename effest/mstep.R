#-------------------------------
# Function to compute the conditional expection of the indicator function
# Author: Qinglong Tian
# Date: 5 June, 2021
#-------------------------------

inv_logit <- function(x)
{
  return(exp(x)/(1+exp(x)))
}

conditionalExpectionVec <- function(mat_spline, mat_x, vec_y, tauVec, betaVec, sigmaScalar)
# mat_spline: is the matrix for the b-spline, columns is the spline, row is the observation
# mat_x: is the model matrix, the first column should be all 1;
# vec_y: is the non-missing response; 
# tauVec: is the estimate of tau from the last iteration
# betaVec: is the estiamte of beta from the last iteration
# sigmaScalar: is the estimate of sigma from the last iteration
{
  n_non_missing <- nrow(mat_spline)
  # n_non_missing <- nrow(mat_x)
  unscaled_condtional_expectation <- numeric(n_non_missing)
  
  if (ncol(mat_x) != length(betaVec))
  {
    stop("Model dimension is wrong!")
  }
  
  if (ncol(mat_spline) != length(tauVec))
  {
    stop("Dimension of the spline does not match!")
  }
  
  for (i in 1:n_non_missing)
  {
    cur_spline <- mat_spline[i,]
    cur_xVec <- mat_x[i,]
    cur_y <- vec_y[i]
    
    prob_non_missing <- inv_logit(sum(cur_spline*tauVec))
    
    mu_temp <- sum(cur_xVec*betaVec)
    sigma_temp <- sigmaScalar
    likModel <- dnorm(cur_y, mean = mu_temp, sd = sigma_temp)
    
    unscaled_condtional_expectation[i] <- (1-prob_non_missing)*likModel
  }
  
  return(unscaled_condtional_expectation/(sum(unscaled_condtional_expectation)))
}
