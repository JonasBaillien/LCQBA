### Skewness measures for the different models mentioned in the paper ###
#########################################################################

# uncomment the desired model to calculate the different skewness measures below


### models introduction plot ###
################################

# ### M1 ###
# d=2
# basefunc=c("normal","logistic")
# alpha=c(0.25,0.65)
# df=c(NA,NA)
# mu=c(20,20)
# A=matrix(c(12,-5,4,8),nrow=d)

# ### M2 ###
# d=2
# basefunc=c("normal","logistic")
# alpha=c(0.25,0.65)
# df=c(NA,NA)
# mu=c(20,20)
# A=matrix(c(7,0,-6,3),nrow=d)

# ### M3 ###
# d=2
# basefunc=c("normal","t")
# alpha=c(0.25,0.65)
# df=c(NA,5)
# mu=c(20,20)
# A=matrix(c(12,0,0,8),nrow=d)
 
 
# ### simulation models ###
# #########################
 
# ### setting 1 ###
# d=2
# A=matrix(c(4,-3,1,4),nrow=d,ncol=d) 
# alpha=c(0.35,0.7)
# basefunc=c("normal","t")
# df=c(NA,6)
 
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
 
 
 
# ### fitted data example models ###
# ################################## 

# ### pokemon ###
# load("~/fittedmodels_pokemon.Rdata") # these are the fits obtained from "FitAllOptions", parameters can be set manually
# pokemonfits=modelfit
# d=6
# A=matrix(pokemonfits$A[75,],nrow=6)
# alpha=pokemonfits$alpha[75,]
# mu=pokemonfits$mu[75,]
# basefunc=pokemonfits$basefunc[75,]
# df=pokemonfits$df[75,]
 

# ### AIS logistic-normal ###
# d=2
# alpha=c(0.2246, 0.3005)
# mu=c(20.1003, 54.5137)
# A=matrix(c(0.4262, 0.6305, 0.5402, 5.1894),nrow=2,ncol=2)
# df=NULL
# basefunc=c("logistic","normal")
 
 
# ### AIS normal-normal ###
# d=2
# alpha=c(0.2178, 0.3020)
# mu=c(20.0532, 54.6272)
# A=matrix(c(0.7490, 0.6493, 0.9029, 5.2251),nrow=2,ncol=2)
# df=NULL
# basefunc=c("normal","normal")




#####################################
### calculating skewness measures ###
#####################################



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







### skewness from the skew-elliptical models ###
################################################


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
