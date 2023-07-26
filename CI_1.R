# Scenario 1 CI
source('LinInterpolation.R')
source('rdunif.R')
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

low_mix = matrix(0,ntestgrid,1)
upp_mix = matrix(0,ntestgrid,1)
length_mix = matrix(0,ntestgrid,1)
low_ds = matrix(0,ntestgrid,1)
upp_ds = matrix(0,ntestgrid,1)
length_ds = matrix(0,ntestgrid,1)
low_li = matrix(0,ntestgrid,1)
upp_li = matrix(0,ntestgrid,1)
length_li = matrix(0,ntestgrid,1)

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

############# fiducial CIs ##############

FiducialLower1=t(fid.lower.test[(mburn+1):(mburn+mfid),])
FiducialUpper1=t(fid.upper.test[(mburn+1):(mburn+mfid),])

linefidm=apply(cbind(FiducialLower1,FiducialUpper1),1,sort)
lowervar_mix=linefidm[2*mfid*(alpha/2),]           
uppervar_mix=linefidm[2*mfid*(1-alpha/2),]

linefidsl=apply(FiducialLower1,1,sort)
linefidsu=apply(FiducialUpper1,1,sort)
lowervar_ds=linefidsl[mfid*(alpha/2),]
uppervar_ds=linefidsu[mfid*(1-alpha/2),]

fid.grid=cbind(fid.lower,fid.upper)[(mburn+1):(mburn+mfid),]
linefid1=t(apply(fid.grid,1,linear_interpolation,grid=grid))
FiducialMidLine1= apply(linefid1,1,function(x){approx(grid,x,testgrid, method = 'linear')$y})
linefidi=apply(FiducialMidLine1,1,sort)
lowervar_li=linefidi[mfid*(alpha/2),]
uppervar_li=linefidi[mfid*(1-alpha/2),]


############# evaluate coverage and length ############
templ_mix =(truth < lowervar_mix)
tempu_mix =(truth > uppervar_mix)
low_mix = low_mix + templ_mix
upp_mix = upp_mix + tempu_mix
templength_mix = uppervar_mix - lowervar_mix
length_mix = length_mix + templength_mix

templ_ds = (truth < lowervar_ds)
tempu_ds = (truth > uppervar_ds)
low_ds = low_ds + templ_ds
upp_ds = upp_ds + tempu_ds
templength_ds = uppervar_ds - lowervar_ds
length_ds = length_ds + templength_ds

templ_li = (truth < lowervar_li)
tempu_li = (truth > uppervar_li)
low_li = low_li + templ_li
upp_li = upp_li + tempu_li
templength_li = uppervar_li - lowervar_li
length_li = length_li + templength_li

}


ave_length_mix=length_mix/nrep
ave_length_ds=length_ds/nrep
ave_length_li=length_li/nrep

