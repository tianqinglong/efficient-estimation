#-------------------------------
# Function: Missing Mechanism
# Authors: Qinglong Tian
#-------------------------------

linear.model.additive <- function(x1, x2, params, sigma = 1)
{
  X <- cbind(1, x1, x2)
  params <- matrix(params, ncol=1)
  return(X %*% params+sigma*rnorm(nrow(X)))
}

linear.model.interaction <- function(x1, x2, params)
# x1 and x2 are vectors
{
  x12 <- x1*x2
  X <- cbind(1, x1, x2, x12)
  params <- matrix(params, ncol=1)
  return(X %*% params)
}

missing.single.y <- function(y, params)
{
  yVec <- c(1, 3*sin(y), y, y^2)
  return(sum(yVec*params))
}

MissingMechanism <- function(y, U, linear.model, params, link="probit")
# 
{
  n <- length(y)
  linear.sum <- linear.model(y, U, params)
  if (link == "probit")
  {
    non.missing.prob <- pnorm(linear.sum)
  }
  else if (link == "logit")
  {
    non.missing.prob <- 1/(1+exp(-linear.sum))
  }
  
  non.missing.ind <- rbinom(n, size = 1, prob = non.missing.prob)
  return(non.missing.ind)
}

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
