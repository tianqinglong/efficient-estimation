library(TwoPhaseReg)
n = 2000
n2 = 600
true_beta = 0.3
true_gamma = 0.4
true_eta = 0.5
seed = 12345
r = 0.3
N_SIEVE = 8

#### Sieve with histogram bases
set.seed(12345)
U2 = runif(n)
simX = runif(n)
simZ = r*simX+U2
simY = true_beta*simX+true_gamma*simZ+rnorm(n)
order.simY = order(simY)
phase2.id = c(order.simY[1:(n2/2)], order.simY[(n-(n2/2)+1):n])
Bspline_Z = matrix(NA, nrow=n, ncol=N_SIEVE)
cut_z = cut(simZ, breaks=quantile(simZ, probs=seq(0, 1, 1/N_SIEVE)), include.lowest = TRUE)
for (i in 1:N_SIEVE) {
    Bspline_Z[,i] = as.numeric(cut_z == names(table(cut_z))[i])
}
colnames(Bspline_Z) = paste("bs", 1:N_SIEVE, sep="")
dat = data.frame(Y=simY, X=simX, Z=simZ, Bspline_Z)
dat[-phase2.id,"X"] = NA

res = smle(Y="Y", X="X", Z="Z", Bspline_Z=colnames(Bspline_Z), data=dat, verbose = T)
res

#### Sieve with linear bases
library(splines)
set.seed(12345)
U1 = runif(n)
U2 = runif(n)
simX = runif(n)
simZ_1 = r*simX+U1
simZ_2 = r*simX+U2
simY = true_beta*simX+true_gamma*simZ_1+true_eta*simZ_2+rnorm(n)
order.simY = order(simY)
phase2.id = c(order.simY[1:(n2/2)], order.simY[(n-(n2/2)+1):n])
bs1 = bs(simZ_1, df=N_SIEVE, degree=1, Boundary.knots=range(simZ_1), intercept=TRUE)
bs2 = bs(simZ_2, df=N_SIEVE, degree=1, Boundary.knots=range(simZ_2), intercept=TRUE)
Bspline_Z = matrix(NA, ncol=N_SIEVE^2, nrow=n)
for (i in 1:ncol(bs1)) {
    for (j in 1:ncol(bs2)) {
        idx = i*N_SIEVE+j-N_SIEVE
        Bspline_Z[,idx] = bs1[,i]*bs2[,j]
    }
}
colnames(Bspline_Z) = paste("bs", 1:ncol(Bspline_Z), sep="")
dat = data.frame(Y=simY, X=simX, Z1=simZ_1, Z2=simZ_2, Bspline_Z)
dat[-phase2.id,"X"] = NA

res = smle(Y="Y", X="X", Z=c("Z1", "Z2"), Bspline_Z=colnames(Bspline_Z), data=dat)
res
