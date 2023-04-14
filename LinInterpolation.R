library('pracma')

linear_interpolation = function(grid, fid, lambda = 1) {
  ## Linear interpolations via quadratic programming
  ngrid = length(grid)
  fid.lower = fid[1:ngrid]
  fid.upper = fid[(ngrid + 1):(2 * ngrid)]
  
  dgrid = 1 / ((ngrid * diff(grid)) ^ 2)
  
  
  off_diag_ind = 2:ngrid
  quad_prod_mat = diag(c(lambda, dgrid) + c(dgrid, lambda))
  quad_prod_mat[cbind(off_diag_ind, off_diag_ind - 1)] = -dgrid
  quad_prod_mat[cbind(off_diag_ind - 1, off_diag_ind)] = -dgrid
  quad_prod_mat = matrix(quad_prod_mat, ngrid, ngrid)
  
  w0 = rbeta(1, 0.5, 0.5) * (fid.upper[1] - fid.lower[1]) + fid.lower[1]
  wEnd = rbeta(1, 0.5, 0.5) * (fid.upper[ngrid] - fid.lower[ngrid]) + fid.lower[ngrid]
  quad_prod_lin = lambda * c(-w0, rep(0, max(0, ngrid - 2)), -wEnd)
  
  inter.grid = quadprog(quad_prod_mat, quad_prod_lin, lb = fid.lower, ub =
                          fid.upper)$xmin
  return(inter.grid)
}