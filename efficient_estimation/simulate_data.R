#-----------------------
# A general function to simulate data with no missing
# Author: Qinglong Tian
# Date: June 17, 2021
#-----------------------

simuY_LM <- function(X, coef, sd)
{
  coef <- matrix(coef, ncol = 1)
  return(X%*%coef+rnorm(nrow(X), 0, sd))
}

simuY <- function(X, coef, sd)
{
  y <- simuY_LM(X, coef, sd)
  return(exp(y)/(1+exp(y)))
}

simuMiss <- function(YU, coef)
{
  coef <- matrix(coef, ncol = 1)
  expLinear <- exp(YU%*%coef)
  prob <- expLinear/(1+expLinear)
  return(ifelse(runif(nrow(YU))<prob, 1, 0))
}
