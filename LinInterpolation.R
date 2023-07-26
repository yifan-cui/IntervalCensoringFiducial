# This code corresponds to Lines 21-25 in Algorithm 1 (Linear interpolations via quadratic programming)
library('pracma')

linear_interpolation = function(grid, fid, lambda = 1) {

  ngrid = length(grid)
  fid.lower = fid[1:ngrid]   # obtain u^L
  fid.upper = fid[(ngrid + 1):(2 * ngrid)] # obtain u^U
  
  dgrid = 1 / ((ngrid * diff(grid)) ^ 2)
  
  
  off_diag_ind = 2:ngrid     # off diagonal index
  quad_prod_mat = diag(c(lambda, dgrid) + c(dgrid, lambda))
  quad_prod_mat[cbind(off_diag_ind, off_diag_ind - 1)] = -dgrid
  quad_prod_mat[cbind(off_diag_ind - 1, off_diag_ind)] = -dgrid
  quad_prod_mat = matrix(quad_prod_mat, ngrid, ngrid) # symmetric matrix representing the quadratic term
  
  w0 = rbeta(1, 0.5, 0.5) * (fid.upper[1] - fid.lower[1]) + fid.lower[1] # randomize u_0 
  wEnd = rbeta(1, 0.5, 0.5) * (fid.upper[ngrid] - fid.lower[ngrid]) + fid.lower[ngrid] # randomize u_{m+1} 
  quad_prod_lin = lambda * c(-w0, rep(0, max(0, ngrid - 2)), -wEnd) # linear term
  
  inter.grid = quadprog(quad_prod_mat, quad_prod_lin, lb = fid.lower, ub =
                          fid.upper)$xmin   # Solve quadratic programming
  return(inter.grid)
}