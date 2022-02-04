##############################################################################
### Simulations + dimensionality test simulation and impact of sample size ###
##############################################################################

# Uncomment the desired model.
# The output is written in a pre specified folder, this can be changed below.
# Please do so when running the simulation. The written output is later also 
# read in again, please also change the correct directory there when needed.

# this is the directory where everything is saved. 
setwd(dir = choose.dir())

# If the file "Code part 1 Functions.R" is not saved in the previously set directory, add the correct path.
source("Code part 1 Functions.R")




### Models ###
##############

## Simulation setting 1
reps=400                      # repetitions
sampsize=800                  # sample size
dims=2                        # dimensions
A=matrix(c(4,-3,1,4),nrow=dims,ncol=dims) # transpose of A is used to generate samples according to the definition
location=c(0,0)               # mu
alpha=c(0.35,0.7)             # alpha
basefunc=c("normal","t")      # type of base function
tpars=c(6)                    # degrees of freedom


# # Simulation setting 2
# reps=400
# dims=6
# sampsize=400 #100, 200, 400, 800
# A=matrix(
#   c(10,0 ,5 ,0 ,1 ,0 ,
#     0 ,10,1 ,0 ,-4,2 ,
#     -5,-1,10,0 ,6 ,0 ,
#     0 ,0 ,0 ,10,0 ,-2,
#     -1,4 ,-6,0 ,10,0 ,
#     0 ,-2,0 ,2 ,0 ,10),nrow=dims,ncol=dims,byrow=T) # of 20 op hoofddiagonaal
# location=1:dims
# alpha=seq(0.2,0.4,length.out = dims)
# basefunc=rep("laplace",dims)
# tpars=NULL





## model for testing the impact of the sample size:

# # model dimensions 6
# reps=400
# dims=6
# sampsize=400 #200, 400, 800
# A=matrix(
#   c(10,0 ,5 ,0 ,1 ,0 ,
#     0 ,10,1 ,0 ,-4,2 ,
#     -5,-1,10,0 ,6 ,0 ,
#     0 ,0 ,0 ,10,0 ,-2,
#     -1,4 ,-6,0 ,10,0 ,
#     0 ,-2,0 ,2 ,0 ,10),nrow=dims,ncol=dims,byrow=T) # of 20 op hoofddiagonaal
# location=1:dims
# alpha=seq(0.2,0.4,length.out = dims)
# basefunc=rep("laplace",dims)
# tpars=NULL

# # for the fitter
# seed=126
# numstarts=50
# maxiter=1000
# ftol=10^-3




# ## model for testing the impact of the dimensionality on the computing time:
# parameters and other requirements

# ## model dimensions 2
# reps=100
# sampsize=10000
# dims=2
# A=matrix(
#   c(20,0,
#     0 ,20),nrow=dims,ncol=dims,byrow=T)
# location=1:dims
# alpha=seq(0.2,0.4,length.out = dims)
# basefunc=rep("laplace",dims)
# tpars=NULL

# ## model dimensions 4
# reps=100
# sampsize=10000
# dims=4
# A=matrix(
#   c(20,0 ,5 ,0 ,
#     0 ,20,1 ,0 ,
#     -5,-1,20,3 ,
#     0 ,0 ,-3,20),nrow=dims,ncol=dims,byrow=T)
# location=1:dims
# alpha=seq(0.2,0.4,length.out = dims)
# basefunc=rep("laplace",dims)
# tpars=NULL



# ## model dimensions 6
# reps=100
# sampsize=10000
# dims=6
# A=matrix(
#   c(20,0 ,5 ,0 ,1 ,0 ,
#     0 ,20,1 ,0 ,-4,2 ,
#     -5,-1,20,3 ,5 ,0 ,
#     0 ,0 ,-3,20,0 ,-2,
#     -1 ,4,-5,0 ,20,0 ,
#     0 ,-2,0 ,2 ,0 ,20),nrow=dims,ncol=dims,byrow=T)
# location=1:dims
# alpha=seq(0.2,0.4,length.out = dims)
# basefunc=rep("laplace",dims)
# tpars=NULL


# ## model dimensions 8
# reps=100
# sampsize=10000
# dims=8
# A=matrix(
#   c(20,0 ,5 ,0 ,1 ,0 ,2 ,0 ,
#     0 ,20,1 ,0 ,-4,2 ,0 ,1 ,
#     -5,-1,20,3 ,5 ,0 ,4 ,1 ,
#     0 ,0 ,-3,20,0 ,-2,-3,2 ,
#     -1 ,4,-5,0 ,20,0 ,0 ,0 ,
#     0 ,-2,0 ,2 ,0 ,20,1 ,-3,
#     -2,0 ,-4,3 ,0 ,-1,20,0 ,
#     0 ,-1,-1,-2,0 ,3 ,0 ,20),nrow=dims,ncol=dims,byrow=T)
# location=1:dims
# alpha=seq(0.2,0.4,length.out = dims)
# basefunc=rep("laplace",dims)
# tpars=NULL


# ## model dimensions 10
# reps=100
# sampsize=10000
# dims=10
# A=matrix(
#   c(20,0 ,5 ,0 ,1 ,0 ,2 ,0 ,-3,0 ,
#     0 ,20,1 ,0 ,-4,2 ,0 ,1 ,-2,3 ,
#     -5,-1,20,3 ,5 ,0 ,4 ,1 ,0 ,-4,
#     0 ,0 ,-3,20,0 ,-2,-3,2 ,-1,0 ,
#     -1 ,4,-5,0 ,20,0 ,0 ,0 ,0 ,2 ,
#     0 ,-2,0 ,2 ,0 ,20,1 ,-3,2 ,0 ,
#     -2,0 ,-4,3 ,0 ,-1,20,0 ,0 ,5 ,
#     0 ,-1,-1,-2,0 ,3 ,0 ,20,5 ,-3,
#     3 ,2 ,0 ,1 ,0 ,-2,0 ,-5,20,0 ,
#     0 ,-3,4 ,0 ,-2,0 ,-5,3 ,0 ,20),nrow=dims,ncol=dims,byrow=T)
# location=1:dims
# alpha=seq(0.2,0.4,length.out = dims)
# basefunc=rep("laplace",dims)
# tpars=NULL

# for the fitter
seed=126
numstarts=50
maxiter=35000 # 75000 solely for 8D and 10D, otherwise 35000
ftol=10^-9
tol=10^-6


### main loop ###
#################
cores=detectCores()
cl <- makeCluster(4)
registerDoParallel(cl)

# create folder to store output in
dir.create(paste0("~/D",dims,"S",sampsize))

result=foreach(i=1:reps,.packages=c('combinat','mrfDepth','nloptr',"optimx"),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 # required files
                 source(file="Code part 1 Functions.R")

                 # generating a sample
                 sample=Xsample(A=A,location=location,basefunc=basefunc,alpha=alpha,sampsize=sampsize,dims=dims,tpars=tpars,seed=seed-i+1)
                 X=sample[[2]]
                 
                 # tracking runtime
                 start.time <- Sys.time()
                 
                 # fitting the model to the data
                 fit=MVQBAFD_fit(data = X,basefunc = basefunc,seed = seed+i,maxiter = maxiter,tol = tol,numstarts=numstarts)
                 
                 
                 fits=c(as.matrix(fit$`fitted parameters`,nrow=1),fit$`log likelihood fit`)
                 
                 # tracking runtime
                 end.time <- Sys.time()
                 time.taken = difftime(end.time,start.time,units = "secs")
                 
                 # writing the output to the specified directory
                 output=c(fits,time.taken)
                 filename=paste0("~/D",dims,"S",sampsize,"/run",i,".csv")
                 write.csv(output,file=filename)
               }
stopCluster(cl)
x11()



### Bundling the output ###
###########################
dataset=data.frame()
ind=c()
for (i in 1:reps){
  try({
    temp_data = t(read.csv(paste0("~/D",dims,"S",sampsize,"/run",i,".csv"), row.names=1))
    
    dataset <- rbind(dataset, temp_data) #for each iteration, bind the new data to the building dataset
    ind=c(ind,i)
  },silent=T)
}

# cleaning up de dataframe
rownames(dataset)=ind
if(is.element("t",basefunc)){
  indt=which(basefunc=="t")
} else {
  indt=NULL
}
a=intToUtf8(945) # character for alpha
m=intToUtf8(956) # character for mu
if(is.null(indt)){
  colnames(dataset)=c(paste0(a,1:dims),paste0(m,1:dims),paste0("A",apply(expand.grid(1:dims, 1:dims), 1, paste, collapse="")),"log likelihood","time taken (s)")
} else {
  colnames(dataset)=c(paste0(a,1:dims),paste0(m,1:dims),paste0("A",apply(expand.grid(1:dims, 1:dims), 1, paste, collapse="")),paste0("df",indt),"log likelihood","time taken (s)")
}
save(dataset,file=paste0("~/D",dims,"S",sampsize,".Rdata"))



### load the datafile ###
#########################
load("~/D6S10000.Rdata")

# change to correct sample size and dimension (for naming purposes)
d=6
sampsize=10000
M=dataset


### sorting parameter estimates because of row-swapping ###
###########################################################

# sorting function specifically for the considered models
# based on the large diagonal elements

sortingdimtest=function(fits,d){
  # fits is a nxd matrix containing the fitted parameters
  # d is the dimensionality of the problem
  
  n=length(fits[,1])
  undecided=c()
  out=fits
  
  for(i in 1:n){
    # create temporary storage for alpha and A as these can be swapped
    af=fits[i,1:d]
    Af=matrix(fits[i,(2*d+1):(2*d+d^2)],nrow=d,ncol=d)
    
    # largest element in each column, this element is to placed on the diagonal
    inds=apply(abs(Af),2,which.max)
    
    
    if(length(unique(inds))==d){
      # if each row has a unique maximal element in the columns, we put these on the diagonal
      Af[1:d,]=Af[inds,]
      af=af[inds]
      
      # negative diagonal elements are not allowed and have to be transformed
      negdiag=which(diag(Af)<0)
      Af[negdiag,]=-Af[negdiag,]
      af[negdiag]=1-af[negdiag]
      
      # # for the t-t case
      # indt=length(which(basefunc=="t"))
      # df=fits[,c(9,10)]
      # df=df[inds]
      # out[,c(9,10)]=df
      
      # output the sorted parameters
      out[i,1:d]=af
      out[i,(2*d+1):(2*d+d^2)]=Af
    } else {
      # if a row has multiple column-wise maxima, a manual sorting has to be applied
      undecided=c(undecided,i)
    }
  }
  return(list("fits"=out,"manual"=undecided))
}
sorted=sortingdimtest(fits=as.matrix(M[,1:(d^2+2*d+2)]),d = d)
M[,1:(d^2+2*d)]=sorted$fits

for(i in (1:100)[-sorted$manual]){
  print(matrix(M[i,21:120],nrow=10))
  print(i)
}

### manually sorting the fits that are undecided ###
####################################################

# just change i and the appropriate row-indices (inds) in Af and af
sorted$manual

i=94

Af=matrix(as.numeric(M[i,(2*dims+1):(2*dims+dims^2)]),nrow=dims)
af=as.numeric(M[i,1:dims])
Af

# check the appropriate elements of the columns to put on the diagonal and put the indices in inds
inds=c(7,3,1,4,6,2,5,8)

Af[1:dims,]=Af[inds,]
af=af[inds]
negdiag=which(diag(Af)<0)
Af[negdiag,]=-Af[negdiag,]
af[negdiag]=1-af[negdiag]
M[i,1:dims]=af
M[i,(2*dims+1):(2*dims+dims^2)]=Af
sorted$manual
# optional save for later use
save(M,file=paste0("~/D",dims,"S",sampsize,"sorted.Rdata"))






### Skewness measures for the different models ###
##################################################

# uncomment the desired model to calculate the different skewness measures below

### models introduction plot
############################

### M1 ###
d=2
basefunc=c("normal","logistic")
alpha=c(0.25,0.65)
df=c(NA,NA)
mu=c(20,20)
A=matrix(c(12,-5,4,8),nrow=d)
# 
# ### M2 ###
# d=2
# basefunc=c("normal","logistic")
# alpha=c(0.25,0.65)
# df=c(NA,NA)
# mu=c(20,20)
# A=matrix(c(7,0,-6,3),nrow=d)
# 
# ### M3 ###
# d=2
# basefunc=c("normal","t")
# alpha=c(0.25,0.65)
# df=c(NA,5)
# mu=c(20,20)
# A=matrix(c(12,0,0,8),nrow=d)
# 
# ### simulation models ###
# #########################
# 
# ### setting 1 ###
# d=2
# A=matrix(c(4,-3,1,4),nrow=d,ncol=d) 
# alpha=c(0.35,0.7)
# basefunc=c("normal","t")
# df=c(NA,6)
# 
# ### setting 2 ###
# d=6
# A=matrix(
#   c(10,0 ,5 ,0 ,1 ,0 ,
#     0 ,10,1 ,0 ,-4,2 ,
#     -5,-1,10,0 ,6 ,0 ,
#     0 ,0 ,0 ,10,0 ,-2,
#     -1,4 ,-6,0 ,10,0 ,
#     0 ,-2,0 ,2 ,0 ,10),nrow=d,ncol=d,byrow=T) 
# mu=1:d
# alpha=seq(0.2,0.4,length.out = d)
# basefunc=rep("laplace",d)
# df=NULL
# 
# ### fitted data example models ###
# 
# ### pokemon ###
# load("~/fittedmodels_pokemon.Rdata")
# pokemonfits=modelfit
# d=6
# A=matrix(pokemonfits$A[75,],nrow=6)
# alpha=pokemonfits$alpha[75,]
# mu=pokemonfits$mu[75,]
# basefunc=pokemonfits$basefunc[75,]
# df=pokemonfits$df[75,]
# 
# 
# ### AIS logistic-normal ###
# d=2
# alpha=c(0.2246, 0.3005)
# mu=c(20.1003, 54.5137)
# A=matrix(c(0.4262, 0.6305, 0.5402, 5.1894),nrow=2,ncol=2)
# df=NULL
# basefunc=c("logistic","normal")
# 
# 
# ### AIS normal-normal ###
# d=2
# alpha=c(0.2178, 0.3020)
# mu=c(20.0532, 54.6272)
# A=matrix(c(0.7490, 0.6493, 0.9029, 5.2251),nrow=2,ncol=2)
# df=NULL
# basefunc=c("normal","normal")



### preliminary calculations ###
################################

# calculation of E[Z]
meanZ=rep(NA,d)
for(j in 1:d){
  meanZ[j]=(1-2*alpha[j])/(alpha[j]*(1-alpha[j]))*switch(basefunc[j]
                                                         ,"normal"=sqrt(2/pi)
                                                         ,"laplace"=1
                                                         ,"logistic"=2*log(2)
                                                         ,"t"=sqrt(df[j]/pi)*(gamma((df[j]-1)/2)/gamma(df[j]/2))
  )
}

# inverse of A
Ainv=solve(A)

# calculation of var(Z) and var(X)
VarZ=rep(NA,d)
for(i in 1:d){
  VarZ[i]=switch(basefunc[i],"normal"=((1-2*alpha[i])^2*(1-2/pi)+alpha[i]*(1-alpha[i])*1)/(alpha[i]*(1-alpha[i]))^2
                 ,"laplace"=((1-2*alpha[i])^2+2*alpha[i]*(1-alpha[i]))/(alpha[i]*(1-alpha[i]))^2
                 ,"logistic"=((1-2*alpha[i])^2*(pi^2/3-4*log(2)^2)+alpha[i]*(1-alpha[i])*pi^2/3)/(alpha[i]*(1-alpha[i]))^2
                 ,"t"=((1-2*alpha[i])^2*(df[i]/(df[i]-2)-df[i]/pi*(gamma((df[i]-1)/2)/gamma(df[i]/2))^2)+alpha[i]*(1-alpha[i])*df[i]/(df[i]-2))/(alpha[i]*(1-alpha[i]))^2
  )
}
VarZ=diag(VarZ)

VarX=t(A)%*%VarZ%*%A

# calculation of Sigma inverse
SigmaInv=solve(VarX)

# calculation of Sigma^(-1/2)*A^T
SIA=expm::sqrtm(SigmaInv)%*%t(A)

# calculation of third order central moments of Z
ThirdCentralMomentZ=rep(NA,d)
skewZ=rep(NA,d)
for(i in 1:d){
  skewZ[i]=switch(basefunc[i]
                  ,"normal"= QBAsyDist::skewAND(alpha = alpha[i])
                  ,"laplace"= QBAsyDist::skewALaD(alpha = alpha[i])
                  ,"logistic"= QBAsyDist::skewALoD(alpha = alpha[i])
                  ,"t"= QBAsyDist::skewATD(alpha = alpha[i],nu = df[i])
  )
}
ThirdCentralMomentZ=skewZ*diag(VarZ)^(3/2)



### mardia's index ###
######################

MARDIA=sum(skewZ^2)


### Mori-Rohatgi-Szekely measure ###
####################################

MRS=matrix(0,nrow=d,ncol=1)
for(m in 1:d){
  for(i in 1:d){
    MRS[m]=MRS[m]+(SIA[i,]^2*SIA[m,])%*%ThirdCentralMomentZ
  }
}


### Kollo measure ###
#####################

KOLLO=matrix(0,nrow=d,ncol=1)
for(k in 1:d){
  KOLLO=KOLLO+sum(SIA[,k]%*%t(SIA[,k]))*ThirdCentralMomentZ[k]*SIA[,k]
}




### empirical versions of the measures ###
##########################################

### AIS data
ais=ais
X=as.matrix(ais[,c(6,9)])
n=nrow(X)
d=ncol(X)


ColMeanX=colMeans(X)
Sinv=solve((n-1)/n*var(X))

Y=sweep(x=X,MARGIN = 2,STATS = ColMeanX)%*%expm::sqrtm(solve(var(X)))
C=sweep(x=X,MARGIN = 2,STATS = ColMeanX)

MARDIA_emp=1/n^2*sum((C%*%Sinv%*%t(C))^3)

MRS_emp=1/n*colSums(sweep(x = Y,MARGIN = 1,STATS = rowSums(Y^2),FUN = "*"))

temp=matrix(NA,nrow=d,ncol=n)
for(i in 1:n){
  temp[,i]=sum(Y[i,]%*%t(Y[i,]))*Y[i,]
}
KOLLO_emp=1/n*rowSums(temp)

### Pokemon data
Pokemon <- read.csv("~/Pokemon.csv", header=FALSE, comment.char="#")
X=as.matrix(Pokemon[,6:11])
n=nrow(X)
d=ncol(X)

ColMeanX=colMeans(X)
Sinv=solve((n-1)/n*var(X))

Y=sweep(x=X,MARGIN = 2,STATS = ColMeanX)%*%expm::sqrtm(solve(var(X)))
C=sweep(x=X,MARGIN = 2,STATS = ColMeanX)

MARDIA_emp=1/n^2*sum((C%*%Sinv%*%t(C))^3)

MRS_emp=1/n*colSums(sweep(x = Y,MARGIN = 1,STATS = rowSums(Y^2),FUN = "*"))

temp=matrix(NA,nrow=d,ncol=n)
for(i in 1:n){
  temp[,i]=sum(Y[i,]%*%t(Y[i,]))*Y[i,]
}
KOLLO_emp=1/n*rowSums(temp)






### skew-normal on AIS data
ais=ais

X=as.matrix(ais[,c(6,9)])
n=1000
d=ncol(X)

fit=msn.mle(x=rep(1,length(X[,1])),y=X)

MARDIA_emp=rep(NA,1000)
MRS_emp=matrix(NA,nrow=1000,ncol=2)
KOLLO_emp=matrix(NA,nrow=1000,ncol=2)
for(i in 1:1000){
  S=rmsn(n = n,dp = fit$dp)
  
  ColMeanS=colMeans(S)
  Sinv=solve((n-1)/n*var(S))
  
  Y=sweep(x=S,MARGIN = 2,STATS = ColMeanS)%*%expm::sqrtm(solve(var(S)))
  C=sweep(x=S,MARGIN = 2,STATS = ColMeanS)
  
  MARDIA_emp[i]=1/n^2*sum((C%*%Sinv%*%t(C))^3)
  
  MRS_emp[i,]=1/n*colSums(sweep(x = Y,MARGIN = 1,STATS = rowSums(Y^2),FUN = "*"))
  
  temp=matrix(NA,nrow=d,ncol=n)
  for(j in 1:n){
    temp[,j]=sum(Y[j,]%*%t(Y[j,]))*Y[j,]
  }
  KOLLO_emp[i,]=1/n*rowSums(temp)
  print(i)
}
mean(MARDIA_emp)
colMeans(MRS_emp)
colMeans(KOLLO_emp)




### skew-t on Pokemon data
Pokemon <- read.csv("~/Pokemon.csv", header=FALSE, comment.char="#")
X=as.matrix(Pokemon[,6:11])
n=1000
d=ncol(X)
r=1000

fit=mst.mple(x=rep(1,length(X[,1])),y=X)

set.seed(12)

# main loop 
cores=detectCores()
cl <- makeCluster(4)
registerDoParallel(cl)

dir.create(paste0("~/MSM"))
result=foreach(i=1:r,.packages=c('sn','alabama'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 S=rmst(n = n,dp = fit$dp)
                 
                 ColMeanS=colMeans(S)
                 Sinv=solve((n-1)/n*var(S))
                 
                 Y=sweep(x=S,MARGIN = 2,STATS = ColMeanS)%*%expm::sqrtm(solve(var(S)))
                 C=sweep(x=S,MARGIN = 2,STATS = ColMeanS)
                 
                 # Mardia
                 MARDIA_emp=1/n^2*sum((C%*%Sinv%*%t(C))^3)
                 
                 # MRS
                 MRS_emp=1/n*colSums(sweep(x = Y,MARGIN = 1,STATS = rowSums(Y^2),FUN = "*"))
                 
                 # Kollo
                 temp=matrix(NA,nrow=d,ncol=n)
                 for(j in 1:n){
                   temp[,j]=sum(Y[j,]%*%t(Y[j,]))*Y[j,]
                 }
                 KOLLO_emp=1/n*rowSums(temp)
                 
                 # writing the output to the specified directory
                 output=list('MAR'=MARDIA_emp,"MRS"=MRS_emp,"KOL"=KOLLO_emp)
                 filename=paste0("~/MSM/run",i,".Rdata")
                 save(output,file=filename)
               }
stopCluster(cl)
x11()

MARDIA_emp=rep(NA,r)
MRS_emp=matrix(NA,nrow=r,ncol=6)
KOLLO_emp=matrix(NA,nrow=r,ncol=6)
for(i in 1:r){
  load(paste0("~/MSM/run",i,".Rdata"))
  MARDIA_emp[i]=output$MAR
  MRS_emp[i,]=output$MRS
  KOLLO_emp[i,]=output$KOL
}

mean(MARDIA_emp)
colMeans(MRS_emp)
colMeans(KOLLO_emp)



