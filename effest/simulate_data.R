#-------------------------------
# Function: Missing Mechanism
# Authors: Qinglong Tian
#-------------------------------

# linear.model.additive <- function(x1, x2, params, sigma = 1)
# {
#   X <- cbind(1, x1, x2)
#   params <- matrix(params, ncol=1)
#   return(X %*% params+sigma*rnorm(nrow(X)))
# }

# linear.model.interaction <- function(x1, x2, params)
# # x1 and x2 are vectors
# {
#   x12 <- x1*x2
#   X <- cbind(1, x1, x2, x12)
#   params <- matrix(params, ncol=1)
#   return(X %*% params)
# }

# missing.single.y <- function(y, params)
# {
#   yVec <- c(1, 3*sin(y), y, y^2)
#   return(sum(yVec*params))
# }

# MissingMechanism <- function(y, U, linear.model, params, link="probit")
# # 
# {
#   n <- length(y)
#   linear.sum <- linear.model(y, U, params)
#   if (link == "probit")
#   {
#     non.missing.prob <- pnorm(linear.sum)
#   }
#   else if (link == "logit")
#   {
#     non.missing.prob <- 1/(1+exp(-linear.sum))
#   }
  
#   non.missing.ind <- rbinom(n, size = 1, prob = non.missing.prob)
#   return(non.missing.ind)
# }

#-------------------------------
# Function: simulate fully observed data with missing indicator
#-------------------------------

SimulateData <- function(n, response.model, params.response,
                         missing.model, params.missing,
                         dim_U=1, dim_Z=1)
{
  U <- matrix(rnorm(n*dim_U), nrow = n)
  Z <- matrix(rnorm(n*dim_Z), nrow = n)
  Y <- response.model(U, Z, params.response)
  
  non.missing.ind <- MissingMechanism(Y, U, missing.model, params.missing)
  
  output <- data.frame(Y=Y, U=U, Z=Z, Obs=non.missing.ind)
  return(output)
}

#-------------------------------
# Rewrite the data generating functions
# Author: Qinglong Tian
# Date: June 10, 2021
#-------------------------------

yTransFunc1 <- function(y)
{
  return(c(1, y, y^2, sin(y)))
}

yTransFunc2 <- function(y)
{
  return(c(1, y))
}

yTransFunc3 <- function(y)
{
  return(c(1, y, y^2))
}

PiY <- function(y, FUN, PARAS, probitLink = T)
{
  yVec <- FUN(y)
  if (length(PARAS) != length(yVec))
  {
    stop("Dimensions of the missing mechanism function did not match!")
  }
  sumY <- sum(yVec*PARAS)

  if (probitLink)
  {
    nonMissProb <- pnorm(sumY)
  }
  else
  {
    nonMissProb <- exp(sumY)/(1+exp(sumY))
  }

  return(nonMissProb)
}

ySimulatorLM <- function(X_mat, betaCoef, sd)
{
  betaCoef <- matrix(betaCoef, ncol = 1)
  return(X_mat%*%betaCoef+rnorm(nrow(X_mat), 0, sd))
}
