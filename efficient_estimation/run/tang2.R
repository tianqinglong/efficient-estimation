library(fastGHQuad)
library(tidyverse)
library(parallel)
library(truncnorm)
source("../simulate_data.R")
Rcpp::sourceCpp("../bspline_recursive.cpp")
source("../add_splines.R")
source("../estep.R")
source("../mstep.R")
source("../variance.R")

total <- 1
cnt <- 1
B <- 1000
n <- 300

beta_0 <- 1
beta_1 <- 1
std <- 1
bn <- 2
q <- 2
gHNodes <- 8 # Gauss-Hermite nodes
max_iter <- 500
tol <- 1e-4

coef1 <- c(beta_0, beta_1)

df_MNAR_list <- NULL
while (total <= B) {
  X <- matrix(rtruncnorm(n, b = 0.7), ncol = 1)
  Y <- simuY(cbind(1, X), coef1, std)
  Z <- X
  U <- NULL
  
  yy <- log(Y/(1-Y))
  YU <- cbind(1, yy)
  coef2 <- c(-1, 1)
  Obs <- simuMiss(YU, coef2)
  
  dat <- cbind(Y, Obs, Z, U)
  colnames(dat) <- c("Y", "Obs", "X1")
  
  df_MNAR <- list(data = dat, Z_indices = 3, U_indices = NULL)
  skip_to_next <- FALSE
  tryCatch(expr = {
    emEstimate <- main(df_MNAR, coef1+.1, std+.1, runif(bn+q), bn, q, gHNodes, max_iter, tol)
    varEst <- ProfileCov(df_MNAR, 1/sqrt(n), emEstimate$Beta, emEstimate$Sigma, emEstimate$Tau, bn, q, gHNodes, bn+q, 1, max_iter, tol)
  },
  warning=function(w){
    skip_to_next <<- TRUE
  },
  finally = {
    if(any(diag(varEst)<0))
    {
      skip_to_next <<- TRUE
    }
  }
  )
  cnt <- cnt+1
  if(skip_to_next)
  {
    next
  }
  df_MNAR_list[[total]] <- list(df = df_MNAR, EM = emEstimate, Var = diag(varEst))
  total <- total+1
}

# saveRDS(df_MNAR_list, file = "tang2.rds")

#-----------------
# Analysis
#-----------------

df_MNAR_list <- readRDS(file = "tang2_d5_1.rds")

count0 <- 0
count1 <- 0
count2 <- 0

beta_0 <- 1
beta_1 <- 1
std <- 1

for (i in 1:length(df_MNAR_list))
{
  res <- df_MNAR_list[[i]]
  # beta0
  beta0 <- res$EM$Beta[1]
  sd_beta0 <- sqrt(res$Var[1])
  count0 <- count0+( (beta_0 >= beta0-1.96*sd_beta0) & (beta_0 <= beta0+1.96*sd_beta0) )
  
  # beta1
  beta1 <- res$EM$Beta[2]
  sd_beta1 <- sqrt(res$Var[2])
  count1 <- count1+( (beta_1 >= beta1-1.96*sd_beta1) & (beta_1 <= beta1+1.96*sd_beta1) )
  
  # sigma
  sigma0 <- res$EM$Sigma
  sd_sigma0 <- sqrt(res$Var[3])
  count2 <- count2+( (std >= max(0,sigma0-1.96*sd_sigma0)) & (std <= sigma0+1.96*sd_sigma0) )
}

count0/length(df_MNAR_list)
count1/length(df_MNAR_list)
count2/length(df_MNAR_list)

betaVal <- sapply(df_MNAR_list, function(x)
{
  x$EM$Beta
})

betaVal <- t(betaVal)
colMeans(betaVal)

