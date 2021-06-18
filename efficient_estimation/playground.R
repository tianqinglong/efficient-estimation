source("simulate_data.R")
Rcpp::sourceCpp("bspline_recursive.cpp")
source("add_splines.R")

n <- 100
X <- matrix(rnorm(2*n), ncol = 2)
Z <- X[,1]
U <- X[,2]

coef1 <- c(1, -1, 1)
Y <- simuY(cbind(1, X), coef1, 1)

YU <- cbind(1, Y, U)
coef2 <- c(1, -2, 1)
Obs <- simuMiss(YU, coef2)

dat <- cbind(Y, Obs, Z, U)
colnames(dat) <- c("Y", "Obs", "Z", "U")

df_MNAR <- list(data = dat, Z_indices = 3, U_indices = 4)

AppendSplines(df_MNAR, 3,3)
