source('LinInterpolation.R')
library(interval)
#########
# input n l r
# input ngrid grid.low grid.high
# input testgrid
mfid = 1000
mburn = 100
alpha = 0.05

grid = seq(grid.low,grid.high,length.out=ngrid)
ntestgrid = length(testgrid)

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

############# fiducial CIs ##############
FiducialLower1=t(fid.lower.test[(mburn+1):(mburn+mfid),])
FiducialUpper1=t(fid.upper.test[(mburn+1):(mburn+mfid),])

linefidm=apply(cbind(FiducialLower1,FiducialUpper1),1,sort)
point_mix=linefidm[2*mfid*0.5,]           

fid.grid=cbind(fid.lower,fid.upper)[(mburn+1):(mburn+mfid),]
linefid1=t(apply(fid.grid,1,linear_interpolation,grid=grid))
FiducialMidLine1= apply(linefid1,1,function(x){approx(grid,x,testgrid, method = 'linear')$y})
linefidi=apply(FiducialMidLine1,1,sort)
point_li=linefidi[mfid*0.5,]



