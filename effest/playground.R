
#-------------------------------
# Test the functions
#-------------------------------

SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))

#-------------------------------
# Test
#-------------------------------

dat <- SimulateData(50, linear.model.additive, c(1, -1, 1), linear.model.interaction, c(1, -1, 1, 0.5))
AddBsplineColumn(dat, splinesForYU, 4, 2)

