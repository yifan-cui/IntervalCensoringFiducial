# current status data of "A unified nonparametric fiducial approach to interval-censored data, Cui et al. (2023)"
library('curstatCI')
source('LinInterpolation.R')
dt <- rubella

#########
n = 230
mfid = 1000
mburn = 100
alpha = 0.05

ngrid = 101 # 325  163
grid.low = 0 #min(l)
grid.high = 81 #max(r)
grid = seq(grid.low,grid.high,length.out=ngrid)

set.seed(2020)

ll <- ifelse(dt$freq1==1, 0, dt$t)
rr <- ifelse(dt$freq1==1, dt$t, Inf)
dt.temp <- dt[which(dt$freq2==2),]
lll <- ifelse(dt.temp$freq1==2, 0, dt.temp$t)
rrr <- ifelse(dt.temp$freq1==2, dt.temp$t, Inf)
l <- c(ll,lll)
r <- c(rr,rrr)
#l <- ll
#r <- rr

############## generate Fid on a pre-sprecified grid ############
nsample <- n
mid=(l+r)/2
vu <- sort(runif(nsample))
vu[order(mid)] <- vu

fid.u <- matrix(NA,mfid+mburn,nsample)
fid.lower <- matrix(NA,mfid+mburn,ngrid)
fid.upper <- matrix(NA,mfid+mburn,ngrid)

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
  
}

############# fiducial CIs ##############

FiducialLower1=t(fid.lower[(mburn+1):(mburn+mfid),])
FiducialUpper1=t(fid.upper[(mburn+1):(mburn+mfid),])

linefidm=apply(cbind(FiducialLower1,FiducialUpper1),1,sort)
lowervar_mix.1=linefidm[2*mfid*(alpha/2),]           
uppervar_mix.1=linefidm[2*mfid*(1-alpha/2),] #alternative, not used

linefidsl=apply(FiducialLower1,1,sort)
linefidsu=apply(FiducialUpper1,1,sort)
lowervar_ds.1=linefidsl[mfid*(alpha/2),]
uppervar_ds.1=linefidsu[mfid*(1-alpha/2),] #alternative, not used

fid.grid=cbind(fid.lower,fid.upper)[(mburn+1):(mburn+mfid),]
linefid1=t(apply(fid.grid,1,linear_interpolation,grid=grid))
FiducialMidLine1= t(linefid1)
linefidi=apply(FiducialMidLine1,1,sort)
lowervar_li.1=linefidi[mfid*(alpha/2),]
var_li.1=linefidi[mfid*0.5,]
uppervar_li.1=linefidi[mfid*(1-alpha/2),]


save.image(file = "rubella_real.RData")

mycol1 <- rgb(0, 0, 255, max = 255, alpha = 20, names = "blue5y0")

pdf("rubella.pdf") 

plot(grid,lowervar_li.1,type = 'l',lty = 'dashed',lwd=2,col = 'blue',ylim=c(0, 1),xlab="t", ylab="F(t)",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(grid,uppervar_li.1,type = 'l',lty = 'dashed',lwd=2,col = 'blue')
lines(grid,var_li.1,type = 'l',lwd=2, col = 'blue')
polygon(c(grid, rev(grid)), c(lowervar_li.1, rev(uppervar_li.1)),
        col=mycol1, border = NA)

dev.off() 

