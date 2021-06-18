#-----------------------
# Append basis-splines to the data list object
# Author: Qinglong Tian
# Date: June 17, 2021
# Note: Input is a data frame containing covariate U
#-----------------------

AppendSplines <- function(dat_list, bn, q)
{
  # Sort the dat
  dat <- dat_list$data
  dat[order(dat[,"Obs"]),] -> dat
  dat_list$data <- dat
  
  # B-splines for Y
  yObs <- dat[dat[,"Obs"] == 1, "Y"]
  bs_y <- matrix(nrow = length(yObs), ncol = bn+q)
  
  delta <- quantile(yObs, probs = seq(0, 1, length.out = (bn+2))[2:(bn+1)])
  delta <- c(rep(0, q), delta, rep(1, q+2))
  for (i in 1:length(yObs))
  {
    bs_y[i,] <- BSplinesForNewY(yObs[i], q, bn, delta)
  }
  colnames(bs_y) <- paste("bs", 1:(bn+q), sep="")
  dat_list$bs_y <- bs_y
  
  # B-splines for U (Use all U, not only the U with non-missing Y)
  posi_u <- dat_list$U_indices
  ## If there is no u
  if (posi_u < 0)
  {
    return(dat_list)
  }
  U <- as.matrix(dat[,posi_u])
  dim_U <- length(posi_u)
  
  for (i in 1:dim_U)
  {
    nam <- paste("bs_u", i, sep = "")
    temp <- splines::bs(U[,i], df = bn+q, degree = q-1, Boundary.knots = range(U[,i]))
    colnames(temp) <- paste("bs", 1:(bn+q), sep="")
    assign(nam, temp)
  }
  
  bs_u_names <- paste("bs_u", 1:dim_U, sep = "")

  for (i in 1:dim_U)
  {
    dat_list[[bs_u_names[i]]] <- get(bs_u_names[i])
  }
  return(dat_list)
}
