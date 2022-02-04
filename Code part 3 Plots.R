# load in the required functions
setwd(dir = choose.dir())
source(file="Code part 1 Functions.R")


###################################################
### contour plots for the 3 introductory models ###
###################################################

# another function for plotting the density
n=400 # grid points in 1 direction
densityf=function(point,basefunc,alpha,mu,B,tpars){
  
  d=length(point)
  
  B=matrix(B,nrow=d)
  
  V=point%*%B
  VV=mu%*%B
  
  dens=det(B)
  
  for(i in 1:d){
    if(basefunc[i]=="laplace"){
      dens=dens*twopiecelaplacedensity(V[i]-VV[i],alpha[i])
    } else if(basefunc[i]=="normal"){
      dens=dens*twopiecenormaldensity(V[i]-VV[i],alpha[i])
    } else if(basefunc[i]=="logistic"){
      dens=dens*twopiecelogisticdensity(V[i]-VV[i],alpha[i])
    } else {
      dens=dens*twopiecestudentdensity(V[i]-VV[i],alpha[i],nu=tpars[i])
    }
  }
  
  return(dens)
}

## model 1
dims=2
A=matrix(c(12,-5,4,8),nrow=dims,ncol=dims) # transpose of A is used to generate samples according to the definition
location=c(20,20)
alpha=c(0.25,0.65)
basefunc=c("normal","logistic")

xtemp=seq(-40,140,length.out=n) # grid in 1 direction
ytemp=seq(-60,100,length.out=n)
xx=as.matrix(expand.grid(xtemp,ytemp)) # square grid

fmat=apply(xx,1,densityf,basefunc=basefunc,alpha=alpha,mu=location,B=solve(A),tpars=0)
fmat=matrix(fmat,nrow=n)

# plot
x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
contour(xtemp,ytemp,fmat,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2)



## model 2
dims=2
A=matrix(c(7,0,-6,3),nrow=dims,ncol=dims) # transpose of A is used to generate samples according to the definition
location=c(20,20)
alpha=c(0.25,0.65)
basefunc=c("normal","logistic")

xtemp=seq(-20,100,length.out=n) # grid in 1 direction
ytemp=seq(-50,50,length.out=n)
xx=as.matrix(expand.grid(xtemp,ytemp)) # square grid

fmat=apply(xx,1,densityf,basefunc=basefunc,alpha=alpha,mu=location,B=solve(A),tpars=0)
fmat=matrix(fmat,nrow=n)
x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
contour(xtemp,ytemp,fmat,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2)


## model 3
dims=2
A=matrix(c(12,0,0,8),nrow=dims,ncol=dims) # transpose of A is used to generate samples according to the definition
location=c(20,20)
alpha=c(0.25,0.65)
basefunc=c("normal","t")
df=c(NA,5)

xtemp=seq(-40,140,length.out=n) # grid in 1 direction
ytemp=seq(-60,100,length.out=n)
xx=as.matrix(expand.grid(xtemp,ytemp)) # square grid

fmat=apply(xx,1,densityf,basefunc=basefunc,alpha=alpha,mu=location,B=solve(A),tpars=df)
fmat=matrix(fmat,nrow=n)
x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
contour(xtemp,ytemp,fmat,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2)



###################################################
### Plots and tables for the simulation studies ###
###################################################


### simulation model 1 ###
##########################
dims=2
A=matrix(c(4,-3,1,4),nrow=dims,ncol=dims) # transpose of A is used to generate samples according to the definition
mu=c(0,0)
alpha=c(0.35,0.7)
basefunc=c("normal","t")
df=c(NA,6)

xtemp=seq(-20,40,length.out=n) # grid in 1 direction
ytemp=seq(-40,20,length.out=n)
xx=as.matrix(expand.grid(xtemp,ytemp)) # square grid

fmat=apply(xx,1,densityf,basefunc=basefunc,alpha=alpha,mu=mu,B=solve(A),tpars=df)
fmat=matrix(fmat,nrow=n)
x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
contour(xtemp,ytemp,fmat,xlab=expression(X[1]),ylab=expression(X[2]),cex.lab=1.5,cex.axis=1.5,lwd=2,cex=2)


reps=400
dims=2
A=matrix(c(4,-3,1,4),nrow=dims,ncol=dims) # transpose of A is used to generate samples according to the definition
mu=c(0,0)
alpha=c(0.35,0.7)
basefunc=c("normal","t")
tpars=c(6)

load("~/M4corrected.Rdata")
M1=Model4_Corrected[[1]]
M2=Model4_Corrected[[2]]
M4=Model4_Corrected[[3]]
M8=Model4_Corrected[[4]]

x11()
boxplot(cbind(M1[,1],M2[,1],M4[,1],M8[,1]),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = alpha[1],col=1)
x11()
boxplot(cbind(M1[,4],M2[,4],M4[,4],M8[,4]),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = mu[2],col=1)
x11()
boxplot(cbind(M1$A12,M2$A12,M4$A12,M8$A12),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = A[1,2],col=1)

FI=FIM(alpha = alpha,mu = mu,A = A,basefunc = basefunc,tpars = tpars)
n=800
vars=diag(solve(FI)/n)

x11()
hist(M8[,1],freq=F,ylim=c(0,16),main="",xlab = "",breaks=seq(0.25,0.45,length.out=13),cex.axis=1.5,cex.lab=1.5,col=NULL)
lines(x=seq(0,1,length.out=1000),dnorm(seq(0,1,length.out=1000),mean = alpha[1],sd = sqrt(vars[1])),lwd=2)

x11()
hist(M8[,4],freq=F,ylim=c(0,0.5),main="",xlab = "",breaks=seq(-3,3,length.out=13),cex.axis=1.5,cex.lab=1.5,col=NULL)
lines(x=seq(-5,5,length.out=1000),dnorm(seq(-5,5,length.out=1000),mean = mu[2],sd = sqrt(vars[4])),lwd=2)

x11()
hist(M8$A12,freq=F,ylim=c(0,2),main="",xlab = "",breaks=seq(0.1,1.9,length.out=13),cex.axis=1.5,cex.lab=1.5,col=NULL)
lines(x=seq(-0,2,length.out=1000),dnorm(seq(-0,2,length.out=1000),mean = A[1,2],sd = sqrt(vars[7])),lwd=2)

### Tables simulation setting 1 (similar one can obtain those for other settings)
M41=Model4_Corrected[[1]]
M42=Model4_Corrected[[2]]
M44=Model4_Corrected[[3]]
M48=Model4_Corrected[[4]]

M41[,9]=round(M41[,9],digits=0)
M42[,9]=round(M42[,9],digits=0)
M44[,9]=round(M44[,9],digits=0)
M48[,9]=round(M48[,9],digits=0)

FIM4=FI(alpha=c(0.35,0.7),mu = c(0,0),A = matrix(c(4,-3,1,4),nrow=2,ncol=2), basefunc = c("normal","t"),tpars = 6 )
varcov4=solve(FIM4)

# boxplot per parameter
pars4=c(0.35,0.7,0,0,4,-3,1,4,6)

a=intToUtf8(945) # character for alpha
m=intToUtf8(956) # character for mu
# tau=intToUtf8(964)
names=c(paste0(a,1:2),paste0(m,1:2),paste0("A",apply(expand.grid(1:2, 1:2), 1, paste, collapse="")),"df")

for(i in 1:9){
  x11()
  boxplot(M41[,i],M42[,i],M44[,i],M48[,i],names=c("n=100","n=200","n=400","n=800"))
  abline(h = pars4[i],col=2,lwd=1)
}

# bias:
bias4=matrix(NA,nrow=4,ncol=9)
colnames(bias4)=names
bias4[1,]=apply(M41[,1:9],2,mean)-pars4
bias4[2,]=apply(M42[,1:9],2,mean)-pars4
bias4[3,]=apply(M44[,1:9],2,mean)-pars4
bias4[4,]=apply(M48[,1:9],2,mean)-pars4

# relative error of variances per parameter
thvars4=diag(varcov4)
relvars4=matrix(NA,nrow=4,ncol=9)
relvars4[1,]=abs(apply(M41[,1:9],2,var)-thvars4/100)/(thvars4/100)
relvars4[2,]=abs(apply(M42[,1:9],2,var)-thvars4/200)/(thvars4/200)
relvars4[3,]=abs(apply(M44[,1:9],2,var)-thvars4/400)/(thvars4/400)
relvars4[4,]=abs(apply(M48[,1:9],2,var)-thvars4/800)/(thvars4/800)

# absolute error of variances per parameter
absvars4=matrix(NA,nrow=4,ncol=9)
absvars4[1,]=(apply(M41[,1:9],2,var)-thvars4/100)
absvars4[2,]=(apply(M42[,1:9],2,var)-thvars4/200)
absvars4[3,]=(apply(M44[,1:9],2,var)-thvars4/400)
absvars4[4,]=(apply(M48[,1:9],2,var)-thvars4/800)

# relative difference of covariance matrix elements
reldifcov41=abs((var(M41[,1:9])-varcov4/100)/(varcov4/100))
reldifcov42=abs((var(M42[,1:9])-varcov4/200)/(varcov4/200))
reldifcov44=abs((var(M44[,1:9])-varcov4/400)/(varcov4/400))
reldifcov48=abs((var(M48[,1:9])-varcov4/800)/(varcov4/800))









### Simulation setting 2 ###
############################
reps=400
dims=6
A=matrix(
  c(10,0 ,5 ,0 ,1 ,0 ,
    0 ,10,1 ,0 ,-4,2 ,
    -5,-1,10,0 ,6 ,0 ,
    0 ,0 ,0 ,10,0 ,-2,
    -1,4 ,-6,0 ,10,0 ,
    0 ,-2,0 ,2 ,0 ,10),nrow=dims,ncol=dims,byrow=T) 
mu=1:dims
alpha=seq(0.2,0.4,length.out = dims)
basefunc=rep("laplace",dims)
tpars=NULL


load("~/D6S100sorted.Rdata")
M1=M
load("~/D6S200sorted.Rdata")
M2=M
load("~/D6S400sorted.Rdata")
M4=M
load("~/D6S800sorted.Rdata")
M8=M

x11()
boxplot(cbind(M1[,3],M2[,3],M4[,3],M8[,3]),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = alpha[3],col=1)
x11()
boxplot(cbind(M1[,5],M2[,5],M4[,5],M8[,5]),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = alpha[5],col=1)
x11()
boxplot(cbind(M1[,9],M2[,9],M4[,9],M8[,9]),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = mu[3],col=1)
x11()
boxplot(cbind(M1[,11],M2[,11],M4[,11],M8[,11]),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = mu[5],col=1)
x11()
boxplot(cbind(M1$A33,M2$A33,M4$A33,M8$A33),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = A[3,3],col=1)
x11()
boxplot(cbind(M1$A35,M2$A35,M4$A35,M8$A35),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = A[3,5],col=1)
x11()
boxplot(cbind(M1$A53,M2$A53,M4$A53,M8$A53),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = A[5,3],col=1)
x11()
boxplot(cbind(M1$A55,M2$A55,M4$A55,M8$A55),names=c("n=100","n=200","n=400","n=800"),cex.axis=1.5)
abline(h = A[5,5],col=1)





### introductory model ###
##########################
reps=400
n= 800
dims=2
A=matrix(c(12,-5,4,8),nrow=dims,ncol=dims) # transpose of A is used to generate samples according to the definition
mu=c(20,20)
alpha=c(0.25,0.65)
basefunc=c("normal","logistic")
tpars=c(NA,2)

FI=FIM(alpha = alpha,mu = mu,A = A,basefunc = basefunc,tpars = NULL)
vars=diag(solve(FI)/n)

load("~/M1_SS800.Rdata")
M=dataset

x11()
hist(M[,2],freq=F,ylim=c(0,20),main="",xlab = "",breaks=seq(0.57,0.72,length.out = 15),cex.axis=1.5,cex.lab=1.5,col=NULL)
lines(x=seq(-0,1,length.out=1000),dnorm(seq(-0,1,length.out=1000),mean = alpha[2],sd = sqrt(vars[2])),lwd=2)

x11()
hist(M$A22,freq=F,ylim=c(0,1.5),main="",xlab = "",breaks=seq(6.7,9,length.out = 15),cex.axis=1.5,cex.lab=1.5,col=NULL)
lines(x=seq(6.4,9.4,length.out=1000),dnorm(seq(6.4,9.4,length.out=1000),mean = A[2,2],sd = sqrt(vars[8])),lwd=2)



xtemp=seq(-40,140,length.out=n) # grid in 1 direction
ytemp=seq(-60,100,length.out=n)
xx=as.matrix(expand.grid(xtemp,ytemp)) # square grid

fmat=apply(xx,1,densityf,basefunc=basefunc,alpha=alpha,mu=mu,B=solve(A),tpars=tpars)
fmat=matrix(fmat,nrow=n)
x11()
contour(xtemp,ytemp,fmat,xlab=expression(X[1]),ylab=expression(X[2]))



##############################################################
### Plots for the impact of sample size and dimensionality ###
##############################################################

### impact sample size and dimension ###
########################################
load("~/D2S10000sorted.Rdata")
D2=M
load("~/D4S10000sorted.Rdata")
D4=M
load("~/D6S10000sorted.Rdata")
D6=M
load("~/D8S10000sorted.Rdata")
D8=M
load("~/D10S10000sorted.Rdata")
D10=M

x11()
boxplot(cbind(D2[,10],D4[,26],D6[,50],D8[,82],D10[,122]),names=c("d=2","d=4","d=6","d=8","d=10"),cex.axis=1.5,cex.lab=1.5,ylab="Time (in s)")

load("~/D6S2000sorted.Rdata")
N2=M
load("~/D6S4000sorted.Rdata")
N4=M
load("~/D6S6000sorted.Rdata")
N6=M
load("~/D6S8000sorted.Rdata")
N8=M
load("~/D6S10000sorted.Rdata")
N10=M

x11()
boxplot(cbind(N2[,50],N4[,50],N6[,50],N8[,50],N10[,50]),names=c("n=2000","n=4000","n=6000","n=8000","n=10000"),cex.axis=1.5,cex.lab=1.5,ylab="Time (in s)")





############################################
### Plots for the real data applications ###
############################################

### contourplots AIS data ###
#############################
n=1000
ais=DAAG::ais
aisdata=ais[,c(6,9)] 


xtemp=seq(15,35,length.out=n) # grid in 1 direction
ytemp=seq(30,110,length.out=n)
xx=as.matrix(expand.grid(xtemp,ytemp)) # square grid

fmat=apply(xx,1,densityf,basefunc=c("normal","normal"),alpha=c(0.2178,0.3020),mu=c(20.0532,54.6272),B=solve(matrix(c(0.7490,0.6493,0.9029,5.2251),nrow=2,ncol=2)),tpars=NULL)
fmat=matrix(fmat,nrow=n)

x11()
plot(aisdata,cex.axis=1.5,cex.lab=1.5,pch=1.5)
contour(xtemp,ytemp,fmat,xlab=expression(X[1]),ylab=expression(X[2]),add=T,lwd=1.5,nlevels = 7)

fmat=apply(xx,1,densityf,basefunc=c("logistic","normal"),alpha=c(0.2246, 0.3005),mu=c(20.1003, 54.5137),B=solve(matrix(c(0.4262, 0.6305, 0.5402, 5.1894),nrow=2,ncol=2)),tpars=NULL)
fmat=matrix(fmat,nrow=n)

x11()
plot(aisdata,cex.axis=1.5,cex.lab=1.5,pch=1.5)
contour(xtemp,ytemp,fmat,xlab=expression(X[1]),ylab=expression(X[2]),add=T,lwd=1.5,nlevels = 7)


fmat=dmsn(x=xx,xi = c(20.1355,61.7612),Omega = matrix(c(16.116,35.3676,35.3676,179.6722),nrow = 2),alpha = c(5.5153,-2.3022))
fmat=matrix(fmat,nrow=n)

x11()
plot(aisdata,cex.axis=1.5,cex.lab=1.5,pch=1.5)
contour(xtemp,ytemp,fmat,xlab=expression(X[1]),ylab=expression(X[2]),add=T,lwd=1.5,nlevels = 7)

fmat=dmst(x=xx,xi = c(20.1979,61.9651),Omega = matrix(c(14.8864,32.6333,32.6333,171.7735),nrow = 2),alpha = c(5.2424,-2.2349),nu = 51.0020)
fmat=matrix(fmat,nrow=n)

x11()
plot(aisdata,cex.axis=1.5,cex.lab=1.5,pch=1.5)
contour(xtemp,ytemp,fmat,xlab=expression(X[1]),ylab=expression(X[2]),add=T,lwd=1.5,nlevels = 7)



### DD plots AIS data ###
#########################
ais=ais
aisdata=ais[,c(6,9)]

Z=as.matrix(aisdata)
dims=length(Z[1,])
n=length(Z[,1])


sampsize=10000
alpha=c(0.21776, 0.30200)
mu=c(20.05317, 54.62719)
A=matrix(c(0.74897, 0.64934, 0.90294, 5.22512),nrow=2)
tpars=c(7.301664)
basefunc=c("normal","normal")
seed=12457

sample=Xsample(A=A,location=mu,basefunc=basefunc,alpha=alpha,sampsize=sampsize,dims=dims,tpars=tpars,seed=seed+23)
X=sample[[2]]

DinDist=hdepth(x=X,z=Z)$depthZ
DinSample=hdepth(x=Z,z=Z)$depthZ
x11()
plot(DinDist,DinSample,xlab="Theoretical depth",ylab="Depth in sample",cex.lab=1.5,cex.axis=1.5)
abline(a=0,b=1,lwd=2)



sampsize=10000
alpha=c(0.2262012, 0.3001808)
mu=c(20.10881, 54.49455)
A=matrix(c(0.6667076, 0.6274283, 0.8514435, 5.182682),nrow=2)
tpars=c(7.301664)
basefunc=c("t","normal")
seed=12457

sample=Xsample(A=A,location=mu,basefunc=basefunc,alpha=alpha,sampsize=sampsize,dims=dims,tpars=tpars,seed=seed+23)
X=sample[[2]]

DinDist=hdepth(x=X,z=Z)$depthZ
DinSample=hdepth(x=Z,z=Z)$depthZ
x11()
plot(DinDist,DinSample,xlab="Theoretical depth",ylab="Depth in sample",cex.lab=1.5,cex.axis=1.5)
abline(a=0,b=1,lwd=2)



sampsize=10000
azfit=msn.mle(x=rep(1,length(Z[,1])),y=Z)
X=rmsn(n=sampsize,dp=azfit$dp)

DinDist=hdepth(x=X,z=Z)$depthZ
DinSample=hdepth(x=Z,z=Z)$depthZ
x11()
plot(DinDist,DinSample,xlab="Theoretical depth",ylab="Depth in sample",cex.lab=1.5,cex.axis=1.5)
abline(a=0,b=1,lwd=2)


sampsize=10000
azfitt=mst.mple(x=rep(1,length(Z[,1])),y=Z)
X=rmst(n=sampsize,dp=azfitt$dp)

DinDist=hdepth(x=X,z=Z)$depthZ
DinSample=hdepth(x=Z,z=Z)$depthZ
x11()
plot(DinDist,DinSample,xlab="Theoretical depth",ylab="Depth in sample",cex.lab=1.5,cex.axis=1.5)
abline(a=0,b=1,lwd=2)




### DD plots Pokemon data ###
#############################

# the data
Pokemon <- read.csv("~/Pokemon.csv", header=FALSE, comment.char="#")
Z=as.matrix(Pokemon[,6:11])

# the fits
load("~/fittedmodels_pokemon.Rdata")

pokemonfits=modelfit
logl=pokemonfits$logl

# 75 corresponds to the best fitting model
A=matrix(pokemonfits$A[75,],nrow=6)
alpha=pokemonfits$alpha[75,]
mu=pokemonfits$mu[75,]
basefunc=pokemonfits$basefunc[75,]
tpars=pokemonfits$df[75,]

# generate a sample to approximate the theoretical distribution
sampsize=5000
dims=6
seed=12457

sample=Xsample(A=A,location=mu,basefunc=basefunc,alpha=alpha,sampsize=sampsize,dims=dims,tpars=tpars,seed=seed+23)
X=sample[[2]]

DinDist=hdepth(x=X,z=Z)$depthZ
DinSample=hdepth(x=Z,z=Z)$depthZ
x11()
plot(DinDist,DinSample,xlab="Theoretical depth",ylab="Depth in sample",cex.axis=1.5,cex.lab=1.5)
abline(a=0,b=1,lwd=2)



# generate a sample to approximate the theoretical distribution
sampsize=10000
azfitt=mst.mple(x=rep(1,length(Z[,1])),y=Z)
X=rmst(n=sampsize,dp=azfitt$dp)

DinDist=hdepth(x=X,z=Z)$depthZ
DinSample=hdepth(x=Z,z=Z)$depthZ
x11()
plot(DinDist,DinSample,xlab="Theoretical depth",ylab="Depth in sample",cex.axis=1.5,cex.lab=1.5)
abline(a=0,b=1,lwd=2)


### heatmaps Pokemon data ###
#############################

Pokemon <- read.csv("~/Pokemon.csv", header=FALSE, comment.char="#")
Z=as.matrix(Pokemon[,6:11])

# cov pokemon data
Z=as.matrix(Pokemon[,6:11])
cormatPD=cor(Z)

# estimated covariance linear combination
load("~/fittedmodels_pokemon.Rdata")

pokemonfits=modelfit
logl=pokemonfits$logl

# 75 is best fitting model
A=matrix(pokemonfits$A[75,],nrow=6)
alpha=pokemonfits$alpha[75,]
mu=pokemonfits$mu[75,]
basefunc=pokemonfits$basefunc[75,]
tpars=pokemonfits$df[75,]

d=6
df=tpars

# mu1 for the contruction of diag(var(Z_i))
mu1=rep(NA,d)
for(i in 1:d){
  mu1[i]=switch(basefunc[i],normal=sqrt(2/pi),
                laplace=1,
                logistic=2*log(2),
                t=(sqrt(df[i])*gamma((df[i]-1)/2))/(sqrt(pi)*gamma(df[i]/2))   )
}

# mu2 for teh contruction of diag(var(Z_i))
mu2=rep(NA,d)
for(i in 1:d){
  mu2[i]=switch(basefunc[i],normal=1,
                laplace=2,
                logistic=pi^2/3,
                t=df[i]/(df[i]-2)   )
}

# diagonal matrix with variances of Z_i
varZ=diag( ( (1-2*alpha)^2*(mu2-mu1^2)+alpha*(1-alpha)*mu2 )/( alpha^2*(1-alpha)^2 ) )


# covariance matrix of the data according to the fitted model
cormatLC=cov2cor(t(A)%*%varZ%*%A)

# estimated covariance mst
azfitt=mst.mple(x=rep(1,length(Z[,1])),y=Z)
S=rmst(n = 100000,dp = azfitt$dp)
cormatMST=cor(S)

# plots of difference in correlation
x11()
image.plot(cormatLC-cormatPD,col = gray.colors(n = 100,start = 1,end = 0),xaxt="n",yaxt="n",breaks = seq(-0.01,0.25,length.out = 101))
axis(side = 1,at=seq(0,1,length.out=6),labels = c("HP","Atk","Def","Sp. Atk","Sp. Def","Spd"),cex.axis=1.2)
axis(side = 2,at=seq(0,1,length.out=6),labels = c("HP","Atk","Def","Sp. Atk","Sp. Def","Spd"),cex.axis=1.2)

x11()
image.plot(cormatMST-cormatPD,col = gray.colors(n = 100,start = 1,end = 0),xaxt="n",yaxt="n",breaks = seq(-0.01,0.25,length.out = 101))
axis(side = 1,at=seq(0,1,length.out=6),labels = c("HP","Atk","Def","Sp. Atk","Sp. Def","Spd"),cex.axis=1.2)
axis(side = 2,at=seq(0,1,length.out=6),labels = c("HP","Atk","Def","Sp. Atk","Sp. Def","Spd"),cex.axis=1.2)

x11()
heatmap(cormatPD/cormatLC,Rowv=NA,Colv = NA,revC = T,labRow = c("HP","Atk","Def","Sp. Atk","Sp. Def","Spd"),labCol = c("HP","Atk","Def","Sp. Atk","Sp. Def","Spd"))
x11()
heatmap(cormatPD/cormatMST,Rowv=NA,Colv = NA,revC = T,labRow = c("HP","Atk","Def","Sp. Atk","Sp. Def","Spd"),labCol = c("HP","Atk","Def","Sp. Atk","Sp. Def","Spd"))
