#-------------------------------
# Functions use to compute the B-splines
# Authors: Qinglong Tian
# Y, U are vector (scalar covariate)
#-------------------------------

library(splines)
# dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))

#-------------------------------
# B-Spline Function for \pi(Y)
#-------------------------------

splinesForY1 <- function(simY, N_SIEVE)
{
  n <- length(simY)
  Bspline_Y <- matrix(nrow = n, ncol = N_SIEVE)
  
  cut_y <- cut(simY,
               breaks = quantile(simY, probs = seq(0, 1, 1/N_SIEVE)),
               include.lowest = TRUE)
  
  for (i in 1:N_SIEVE)
  {
    Bspline_Y[,i] <- as.numeric(cut_y == names(table(cut_y))[i])
  }
  colnames(Bspline_Y) <- paste("bs", 1:N_SIEVE, sep="")
  
  return(Bspline_Y)
}

splinesForY <- function(simY, simU, N_SIEVE, degree)
{
  if (degree == 1)
  {
    return(splinesForY1(simY, N_SIEVE))
  }
  
  n <- length(simY)
  bsDegree <- degree-1
  bs1 <- bs(simY, df = N_SIEVE, degree = bsDegree,
            Boundary.knots = range(simY), intercept = TRUE)
  colnames(bs1) <- paste("bs", 1:ncol(bs1), sep = "")
  return(bs1)
}

#-------------------------------
# B-Spline Function for \pi(Y, U)
#-------------------------------

splinesForYU1 <- function(simY, simU, N_SIEVE)
# simU: original U after removing obs that missing Y
{
  n <- length(simY)
  # A quick check
  if (n != length(simU))
  {
    stop("Y and U do not match!")
  }
  
  Bspline_Y <- matrix(nrow = n, ncol = N_SIEVE)
  Bspline_U <- matrix(nrow = n, ncol = N_SIEVE)
  
  cut_y <- cut(simY, breaks = quantile(simY, probs = seq(0, 1, 1/N_SIEVE)),
               include.lowest = TRUE)
  cut_u <- cut(simU, breaks = quantile(simU, probs = seq(0, 1, 1/N_SIEVE)),
               include.lowest = TRUE)
  
  for (i in 1:N_SIEVE) {
    Bspline_Y[,i] = as.numeric(cut_y == names(table(cut_y))[i])
    Bspline_U[,i] = as.numeric(cut_u == names(table(cut_u))[i])
  }
  
  Bspline_YU <- matrix(nrow = n, ncol=N_SIEVE^2)
  
  for (i in 1:ncol(Bspline_Y))
    for (j in 1:ncol(Bspline_U))
    {
      idx <- i*N_SIEVE+j-N_SIEVE
      Bspline_YU[,idx] <- Bspline_Y[,i]*Bspline_U[,j]
    }
  colnames(Bspline_YU) <- paste("bs", 1:ncol(Bspline_YU), sep="")
  return(Bspline_YU)
}

splinesForYU <- function(simY, simU, N_SIEVE, degree)
{
  if (degree == 1)
  {
    return(splinesForYU1(simY, simU, N_SIEVE))
  }
  
  n <- length(simY)
  bsDegree <- degree-1
  bsY <- bs(simY, df = N_SIEVE, degree = bsDegree, Boundary.knots = range(simY), intercept = TRUE)
  bsU <- bs(simU, df = N_SIEVE, degree = bsDegree, Boundary.knots = range(simU), intercept = TRUE)
  
  Bspline_YU <- matrix(ncol = N_SIEVE^2, nrow = n)
  for (i in 1:ncol(bsY))
    for (j in 1:ncol(bsU))
    {
      idx <- i*N_SIEVE+j-N_SIEVE
      Bspline_YU[, idx] <- bsY[,i]*bsU[,j]
    }
  
  colnames(Bspline_YU) <- paste("bs", 1:ncol(Bspline_YU), sep = "")
  return(Bspline_YU)
}

#-------------------------------
# Add B-spline column to the data set
#-------------------------------

AddBsplineColumn <- function(dat, FUN, N_SIEVE, degree)
{
  nTot <- length(dat$Obs)
  numMissing <- sum(dat$Obs == 0)
  numNonMissing <- nTot-numMissing
  
  dat <- dat[order(dat$Obs),]
  datMissing <- dat[1:numMissing,]
  datNonMissing <- dat[(numMissing+1):nTot,]
  
  bspline_non_missing <- FUN(datNonMissing$Y, datNonMissing$U, N_SIEVE, degree)
  bspline_ncol <- ncol(bspline_non_missing)
  bspline_missing <- matrix(nrow = numMissing, ncol = bspline_ncol)
  
  bspline_mat <- rbind(bspline_missing, bspline_non_missing)
  datOut <- cbind(dat, bspline_mat)
  
  return(datOut)
}
