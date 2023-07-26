# Scenario 1 point estimation
source('LinInterpolation.R')
source('rdunif.R')
library(interval)
#########
n = 50 # change n to other sample sizes
nrep = 1000
mfid = 1000
mburn = 100
alpha = 0.05

ngrid = 101
grid.low = 0 #min(l)
grid.high = 5 #max(r)
grid = seq(grid.low,grid.high,length.out=ngrid)
testgrid = c(log(4/3),log(2),log(4))
truth = pexp(testgrid)
ntestgrid = length(testgrid)

mse_mix = matrix(0,ntestgrid,1)
mse_li = matrix(0,ntestgrid,1)
mse_mle1 = matrix(0,ntestgrid,1)
mse_mle2 = matrix(0,ntestgrid,1)
mse_mle3 = matrix(0,ntestgrid,1)

for (k in 1:nrep){
  
set.seed(2020*k)
t <- rexp(n, 1)
l <- rep(NA,n)
r <- rep(NA,n)
for (i in 1:n){
c <- rexp(1,1)
if(length(c[c < t[i]])>0){
l[i] <- max(c[c < t[i]])
} else l[i] <- 0
if(length(c[c > t[i]])>0) {
r[i] <- min(c[c > t[i]])
} else r[i] <- Inf
}



############## generate Fid on a pre-sprecified grid ############
nsample <- n
mid=(l+r)/2
vu <- sort(runif(nsample))
vu[order(mid)] <- vu

fid.u <- matrix(NA,mfid+mburn,nsample)
fid.lower <- matrix(NA,mfid+mburn,ngrid)
fid.upper <- matrix(NA,mfid+mburn,ngrid)
fid.lower.test <- matrix(NA,mfid+mburn,ntestgrid)
fid.upper.test <- matrix(NA,mfid+mburn,ntestgrid)

for (j in 1:(mfid+mburn)){
for (i in 1:nsample){
  vu.pre <- vu[-i]
  l.pre <- l[-i]
  r.pre <- r[-i]
  index1 <- which(l[i]>=r.pre)
  if(length(index1)>0){
  u.lower <- max(vu.pre[index1]) }
  else u.lower <- 0  ##### l 
  index2 <- which(r[i]<=l.pre)
  if(length(index2)>0){
  u.upper <- min(vu.pre[index2])  
  } else u.upper <- 1 ##### r
  vu[i] <- runif(1,min=u.lower, max=u.upper)
} 
  
temp <- runif(nsample)   
vu <- sort(temp)[rank(vu)]  
fid.u[j,] <- vu

fiducial.sample=sapply(grid,function(s){
  index1 <- which(s>=r)
  if(length(index1)>0){
    lower <- max(vu[index1]) }
  else lower <- 0  ##### l 
  index2 <- which(s<l)
  if(length(index2)>0){
    upper <- min(vu[index2])  
  } else upper <- 1 ##### r
  return(rbind(lower,upper))
})

fid.lower[j,]=fiducial.sample[1,]
fid.upper[j,]=fiducial.sample[2,]


fiducial.sample.test=sapply(testgrid,function(s){
  index1 <- which(s>=r)
  if(length(index1)>0){
    lower <- max(vu[index1]) }
  else lower <- 0  ##### l 
  index2 <- which(s<l)
  if(length(index2)>0){
    upper <- min(vu[index2])  
  } else upper <- 1 ##### r
  return(rbind(lower,upper))
})

fid.lower.test[j,]=fiducial.sample.test[1,]
fid.upper.test[j,]=fiducial.sample.test[2,]

}

############# fiducial estimator ##############
FiducialLower1=t(fid.lower.test[(mburn+1):(mburn+mfid),])
FiducialUpper1=t(fid.upper.test[(mburn+1):(mburn+mfid),])

linefidm=apply(cbind(FiducialLower1,FiducialUpper1),1,sort)
point_mix=linefidm[2*mfid*0.5,]           

# linefidsl=apply(FiducialLower1,1,sort)
# linefidsu=apply(FiducialUpper1,1,sort)
# lowervar_ds=linefidsl[mfid*(alpha/2),]
# uppervar_ds=linefidsu[mfid*(1-alpha/2),]

fid.grid=cbind(fid.lower,fid.upper)[(mburn+1):(mburn+mfid),]
linefid1=t(apply(fid.grid,1,linear_interpolation,grid=grid))
FiducialMidLine1= apply(linefid1,1,function(x){approx(grid,x,testgrid, method = 'linear')$y})
linefidi=apply(FiducialMidLine1,1,sort)
point_li=linefidi[mfid*0.5,]


fit <- icfit(l,r)
point_mle1 = 1 - getsurv(times = testgrid, fit, nonUMLE.method = "interpolation")[[1]]$S
point_mle2 = 1 - getsurv(times = testgrid, fit, nonUMLE.method = "left")[[1]]$S
point_mle3 = 1 - getsurv(times = testgrid, fit, nonUMLE.method = "right")[[1]]$S


############# evaluate mse ############

mse_mix = mse_mix + (point_mix-truth)^2
mse_li = mse_li + (point_li-truth)^2
mse_mle1 = mse_mle1 + (point_mle1-truth)^2
mse_mle2 = mse_mle2 + (point_mle2-truth)^2
mse_mle3 = mse_mle3 + (point_mle3-truth)^2


}


mse_mix = mse_mix/nrep
mse_li = mse_li/nrep
mse_mle1 = mse_mle1/nrep
mse_mle2 = mse_mle2/nrep
mse_mle3 = mse_mle3/nrep


