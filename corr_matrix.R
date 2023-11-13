corr_matrix <- function(x1, x2, rhosq) {
  "
  x1 is the first set of points in the correlation matrix
  x2 is the second set of points in the correlation matrix
  N1 is the number of points in x1
  N2 is the number of points in x2
  rhosq are the range parameters
  "
  N1 = dim(x1)[1]
  N2 = dim(x2)[1]
  
  rhosq = as.numeric(rhosq)
  mat52 = function(d) (1 + sqrt(5) * d + (5/3) * d^2) * exp(-sqrt(5) * d)  # Matern 5/2 correlation function
  
  x1_rep = array(x1, dim = c(N1, dim(x1)[2],N2))
  x2_rep = array(x2, dim = c(N2, dim(x2)[2],N1))
  x1_diff = x1_rep - aperm(x2_rep, c(3,2,1))
  d = aperm(sqrt(aperm(x1_diff^2, c(2,1,3)) / rhosq), c(2,1,3))
  M = apply(mat52(d), MARGIN = c(1,3), prod)
  
  return(t(M))
}
