source("simulate_data.R")
Rcpp::sourceCpp("bspline_recursive.cpp")
source("add_splines.R")

# Test data generating functions

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

AppendSplines(df_MNAR, 3,3) -> datList

# Test "YSplinePrep"

xVec <- c(1, 1.5, 2)
rules <- gaussHermiteData(10)
q <- 3
bn <- 3
delta <- seq(0, 1, length.out = bn+2)
delta <- c(rep(0, q-1), delta, rep(1, q+1))
betaOld <- c(1, 1, 1)
sigmaOld <- 1

YSplinePrep(xVec, rules, delta, bn, q, betaOld, sigmaOld)

# Test YUSplinePrep

x <- datList$bs_u1[1,]
xBS <- as.matrix(x)
yBS <- datList$bs_y[1,]

YUSplinePrep(yBS, xBS)

# Test PseudoObservation

yVec <- runif(10)
w <- rules$w
ubsMat <- YUSplinePrep(yBS, xBS)
tauOld <- rnorm(length(ubsMat))

PseudoObservation(yVec, w, 2, 2.4, xBS,  delta, bn, q, tauOld)

# Test Make Dataset
MakeFullDataSetMissing(datList, rules, betaOld, sigmaOld, tauOld)

#
df_MNAR <- list(data = dat, Z_indices = c(3, 4), U_indices = NULL)
AppendSplines(df_MNAR, bn = 3, q = 3) -> datList
