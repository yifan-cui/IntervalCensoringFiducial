# Fiducial algorithm for Interval-Censoring Data
algorithm.R includes the implementation of fiducial MCMC estimators for interval-censoring data; LinInterpolation.R includes the implementation of linear interpolations via quadratic programming.

The main algorithm algorithm.R inputs:
n (# of samples) 
l (a vector representing the left end) 
r (a vector representing the right end)
grid (t_grid in the paper)
testgrid (a grid of values one would like to test on)
and outputs:
FiducialLower1 (lower fiducial samples) 
FiducialUpper1 (upper fiducial samples)
point_li (fiducial point estimator)

To reproduce tables in the simulations of the paper, run 

To reproduce results in the real data analysis of the paper, run 
