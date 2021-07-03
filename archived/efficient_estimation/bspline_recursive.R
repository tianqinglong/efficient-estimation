BSplineForNewY_R <- function(y, q, bn, delta)
{
  dict <- matrix(nrow = bn+2*q-1, ncol = q)
  
  for (j in 1:(bn+2*q-1))
  {
    tLower <- delta[j]
    tUpper <- delta[j+1]
    dict[j,1] <- ifelse(y > tLower & y <= tUpper, 1, 0)
  }
  
  for (i in 2:q)
    for (j in 1:(bn+2*q-i))
    {
      # if (i == 2 & j == 3)
      # {
      #   browser()
      # }
      w1 <- ifelse(delta[j+i-1]==delta[j], 0, (y-delta[j])/(delta[j+i-1]-delta[j]))
      w2 <- ifelse(delta[j+i]==delta[j+1], 0, (delta[j+i]-y)/(delta[j+i]-delta[j+1]))
      dict[j,i] <- w1*dict[j,i-1]+w2*dict[j+1,i-1]
    }
  
  return(dict[1:(bn+q),q])
}
# 
# q <- 4
# bn <- 5
# delta <- seq(0, 1, length.out = bn+2)
# delta <- c(rep(0, q-1), delta, rep(1, q+1))
# 
# sampledPoint <- 1000
# yVec <- seq(0, 1-1e-6, length.out = sampledPoint)
# bspline <- matrix(nrow = sampledPoint, ncol = bn+q)
# 
# for (i in 1:sampledPoint)
# {
#   yVal <- yVec[i]
#   bsVal <- BSplinesForNewY(yVal, q, bn, delta)
#   for (j in 1:(bn+q))
#   {
#     bspline[i, j] <- bsVal[j]
#   }
# }
# 
# bspline <- cbind(yVec, bspline)
# bspline <- as.data.frame(bspline)
# 
# reshape2::melt(bspline, id.vars = 'yVec', variable.names='bs') -> bspline_plot
# 
# library(ggplot2)
# ggplot(bspline_plot, aes(x = yVec, y = value, col = variable))+geom_line(aes(linetype = variable))
# 
# yVal <- 0.4
# microbenchmark::microbenchmark(BSplineForNewY_R(yVal, q, bn, delta), BSplinesForNewY(yVal, q, bn, delta))
