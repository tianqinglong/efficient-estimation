#-----------------------
# A general function to simulate data with no missing
# Author: Qinglong Tian
# Date: June 17, 2021
#-----------------------

simuY_LM <- function(X, coef, sd=1)
{
  coef <- matrix(coef, ncol = 1)
  rnorm(nrow(X), X%*%coef, sd)
}

simuY <- function(X, coef, sd)
{
  y <- simuY_LM(X, coef, sd)
  return(1/(1+exp(-y)))
}

simuMiss <- function(YU, coef, use_logit = T)
{
  coef <- matrix(coef, ncol = 1)
  if(use_logit)
  {
    prob <- 1/(1+exp(-YU %*% coef))
  }
  else
  {
    prob <- pnorm(YU %*% coef)
  }
  return(ifelse(runif(nrow(YU))<prob, 1, 0))
}
