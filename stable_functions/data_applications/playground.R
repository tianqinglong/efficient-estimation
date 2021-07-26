library(readr)
library(tidyverse)
library(fastGHQuad)
Rcpp::sourceCpp("../bspline_recursive.cpp")
source("../add_splines.R")
source("../estep.R")
source("../mstep.R")
source("../loglikelihood.R")

MIMIC <- read_csv("../../data/MIMICchenchi.csv",
                  col_types = cols(gender = col_factor(levels = c())))

# Select the rearrange data

data <- MIMIC %>%
  select(calcium_mean, redbloodcell_mean, magnesium_mean, sofa, albumin_mean) %>%
  mutate(Obs = !is.na(albumin_mean))
data <- data[, c(5, 6, 1, 2, 3, 4)]

# Use the simulation routines

original_names <- names(data)
colnames(data) <- c("Y", "Obs", "X1", "X2", "X3", "X4")
data2 <- data %>% mutate(Y = 1/(1+exp(-Y)))

data_MAR <- data2 %>% filter(Obs == 1)
data2 <- data2 %>% as.matrix

fitMar <- lm(Y~X1+X2+X3+X4, data = data_MAR)

beta_mar <- fitMar$coefficients
sigma_mar <- sigma(fitMar)
hn <- min(sqrt(diag(vcov(fitMar))))
numObs <- nrow(data)
hn <- min(hn, 1/sqrt(numObs))

df_MNAR <- list(data = data2, Z_indices = 3, U_indices = c(4, 5, 6))

num_u <- 3
num_z <- 1
bn <-2
q <- 2
nsieves_add <- (bn+q-1)*(num_u+1)+1
n_covarites <- (num_z+num_u)
ghn <- 8
max_iter <- 500
tol <- 1e-4
emEstimate <- main_additive(df_MNAR, beta_mar, sigma_mar, rep(0, times = nsieves_add), bn, q, ghn, max_iter, tol)
beta_mle <- emEstimate$Beta
sd_mle <- emEstimate$Sigma
tau_mle <- emEstimate$Tau

temp <- ProfileCov_additive(df_MNAR, hn*2, beta_mle, sd_mle, tau_mle, bn, q, ghn, nsieves_add, n_covarites, max_iter, tol)
names(beta_mle) <- c("(Intercept)", original_names[-c(1,2)])
cbind(beta_mle-qnorm(0.975)*sqrt(diag(temp))[1:length(beta_mle)],
      beta_mle+qnorm(0.975)*sqrt(diag(temp))[1:length(beta_mle)]) -> ci
ci
