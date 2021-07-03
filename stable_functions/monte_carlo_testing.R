library(fastGHQuad)
library(tidyverse)
library(parallel)
source("../simulate_data.R")
Rcpp::sourceCpp("../bspline_recursive.cpp")
source("../add_splines.R")
source("../estep.R")
source("../mstep.R")
source("../loglikelihood.R")

B <- 3000
n <- 300
coef1 <- c(1,1)
sd <- 1
ghn <- 8
bn <- 2
q <- 2
max_iter <- 300
tol <- 1e-5

out <- mclapply(1:B, function(x)
{
  # Simulate data
  X <- matrix(truncnorm::rtruncnorm(n, b = 0.7), ncol = 1)
  Y <- simuY(cbind(1, X), coef1, sd)
  Z <- X
  U <- NULL # \pi(Y)
  nsieves <- bn+q
  yy <- log(Y/(1-Y))
  YU <- cbind(1, yy)
  coef2 <- c(-1, 1)
  Obs <- simuMiss(YU, coef2)
  dat <- cbind(Y, Obs, Z, U)
  colnames(dat) <- c("Y", "Obs", "X1")
  
  # Missing at random
  yObs <- Y[which(Obs == 1)]
  xObs <- X[which(Obs == 1),]
  yLM <- log(yObs/(1-yObs))
  lmMAR <- lm(yLM~xObs)
  hn <- min(sqrt(diag(vcov(lmMAR))))
  
  # EM
  df_MNAR <- list(data = dat, Z_indices = 3, U_indices = NULL)
  
  stop_this_trial <- FALSE
  tryCatch(
    expr = {
      emEstimate <- main(df_MNAR, coef1+0.5, sd+0.2, runif(bn+q, min = -1, max = 1), bn, q, ghn, max_iter, tol)
    },
    warning = function(w){
      stop_this_trial <<- TRUE
    },
    error = function(e){
      stop_this_trial <<- TRUE
    }
  )
  if(stop_this_trial)
  {
    return(list(Data = df_MNAR, EM = "EM abnormal", Var = "EM abnormal!"))
  }
  
  # Profile likelihood
  beta_mle <- emEstimate$Beta
  sd_mle <- emEstimate$Sigma
  tau_mle <- emEstimate$Tau
  
  change_hn <- FALSE
  tryCatch(
    expr = {
      covMat <- ProfileCov(df_MNAR, min(hn, 1/sqrt(n)), beta_mle, sd_mle, tau_mle, bn, q, ghn, nsieves, 1, max_iter, tol)
    },
    error = function(e){
      change_hn <<- TRUE
    }
  )
  
  if( change_hn == FALSE & (any(is.nan(covMat)) | any(is.na(covMat))) )
  {
    change_hn <- TRUE
  }
  
  if( change_hn == FALSE & any(diag(covMat)<0) )
  {
    change_hn <- TRUE
  }
  
  if( change_hn )
  {
    continue_hn2 <- FALSE
    hn1 <- min(hn, 1/sqrt(n))/5
    tryCatch(
      expr = {
        covMat <- ProfileCov(df_MNAR, hn1, beta_mle, sd_mle, tau_mle, bn, q, ghn, nsieves, 1, max_iter, tol)
      },
      error = function(e){
        continue_hn2 <<- TRUE
      }
    )
    
    if (exists("covMat"))
    {
      if (any(is.nan(covMat)) | any(is.na(covMat)))
      {
        continue_hn2 <- TRUE
      }
      else if(any(diag(covMat)<0))
      {
        continue_hn2 <- TRUE
      }
    }
    
    if(continue_hn2)
    {
      hn1 <- min(hn, 1/sqrt(n))/10
      tryCatch(
        expr = {
          covMat <- ProfileCov(df_MNAR, hn1, beta_mle, sd_mle, tau_mle, bn, q, ghn, nsieves, 1, max_iter, tol)
        },
        error = function(e){
          stop_this_trial <<- TRUE
        }
      )
    }
    
    if(stop_this_trial)
    {
      return(list(Data = df_MNAR, EM = emEstimate, Var = "Error profiling!"))
    }
  }

  if (any(is.nan(covMat)) | any(is.na(covMat)))
  {
    return(list(Data = df_MNAR, EM = emEstimate, Var = "NA/NaN in covariance matrix!"))
  }
  else if(any(diag(covMat)<0))
  {
    return(list(Data = df_MNAR, EM = emEstimate, Var = "Hessian is not PD!"))
  }

  return(list(Data = df_MNAR, EM=emEstimate, Var = covMat))
},
mc.cores = 16
)

saveRDS(out, file = "rout1.rds")
