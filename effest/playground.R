library(tidyverse)
#-------------------------------
# Test the functions
#-------------------------------

SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))

#-------------------------------
# Test
#-------------------------------

dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))
datWBspline <- AddBsplineColumn(dat, splinesForY, 6, 2)

#-------------------------------
# Test
#-------------------------------

datWBspline %>% filter(Obs == 1) -> datNonMissing
datWBspline %>% filter(Obs == 0) -> datMissing

yVec <- datNonMissing[,1]
mat_Spline <- datNonMissing[,-c(1, 2, 3, 4)]
mat_X <- cbind(Intercept=1, datMissing[, c(2, 3)])

betaVec <- c(0.9, -0.5, 0.9)
tauVec <- rnorm(6)
sigmaSCL <- sd(yVec)

conditionalExpectionVec_OnlyY(mat_Spline, mat_X, yVec, tauVec, betaVec, sigmaSCL) -> condProb
