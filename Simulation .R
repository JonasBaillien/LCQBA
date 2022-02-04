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
d=2                        # dimensions
A=matrix(c(4,-3,1,4),nrow=d,ncol=d) # transpose of A is used to generate samples according to the definition
location=c(0,0)               # mu
alpha=c(0.35,0.7)             # alpha
basefunc=c("normal","t")      # type of base function
tpars=c(6)                    # degrees of freedom


# # Simulation setting 2
# reps=400
# d=6
# sampsize=400 #100, 200, 400, 800
# A=matrix(
#   c(10,0 ,5 ,0 ,1 ,0 ,
#     0 ,10,1 ,0 ,-4,2 ,
#     -5,-1,10,0 ,6 ,0 ,
#     0 ,0 ,0 ,10,0 ,-2,
#     -1,4 ,-6,0 ,10,0 ,
#     0 ,-2,0 ,2 ,0 ,10),nrow=d,ncol=d,byrow=T) # of 20 op hoofddiagonaal
# location=1:d
# alpha=seq(0.2,0.4,length.out = d)
# basefunc=rep("laplace",d)
# tpars=NULL





## model for testing the impact of the sample size:

# # model dimensions 6
# reps=400
# d=6
# sampsize=400 #200, 400, 800
# A=matrix(
#   c(10,0 ,5 ,0 ,1 ,0 ,
#     0 ,10,1 ,0 ,-4,2 ,
#     -5,-1,10,0 ,6 ,0 ,
#     0 ,0 ,0 ,10,0 ,-2,
#     -1,4 ,-6,0 ,10,0 ,
#     0 ,-2,0 ,2 ,0 ,10),nrow=d,ncol=d,byrow=T) # of 20 op hoofddiagonaal
# location=1:d
# alpha=seq(0.2,0.4,length.out = d)
# basefunc=rep("laplace",d)
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
# d=2
# A=matrix(
#   c(20,0,
#     0 ,20),nrow=d,ncol=d,byrow=T)
# location=1:d
# alpha=seq(0.2,0.4,length.out = d)
# basefunc=rep("laplace",d)
# tpars=NULL

# ## model dimensions 4
# reps=100
# sampsize=10000
# d=4
# A=matrix(
#   c(20,0 ,5 ,0 ,
#     0 ,20,1 ,0 ,
#     -5,-1,20,3 ,
#     0 ,0 ,-3,20),nrow=d,ncol=d,byrow=T)
# location=1:d
# alpha=seq(0.2,0.4,length.out = d)
# basefunc=rep("laplace",d)
# tpars=NULL



# ## model dimensions 6
# reps=100
# sampsize=10000
# d=6
# A=matrix(
#   c(20,0 ,5 ,0 ,1 ,0 ,
#     0 ,20,1 ,0 ,-4,2 ,
#     -5,-1,20,3 ,5 ,0 ,
#     0 ,0 ,-3,20,0 ,-2,
#     -1 ,4,-5,0 ,20,0 ,
#     0 ,-2,0 ,2 ,0 ,20),nrow=d,ncol=d,byrow=T)
# location=1:d
# alpha=seq(0.2,0.4,length.out = d)
# basefunc=rep("laplace",d)
# tpars=NULL


# ## model dimensions 8
# reps=100
# sampsize=10000
# d=8
# A=matrix(
#   c(20,0 ,5 ,0 ,1 ,0 ,2 ,0 ,
#     0 ,20,1 ,0 ,-4,2 ,0 ,1 ,
#     -5,-1,20,3 ,5 ,0 ,4 ,1 ,
#     0 ,0 ,-3,20,0 ,-2,-3,2 ,
#     -1 ,4,-5,0 ,20,0 ,0 ,0 ,
#     0 ,-2,0 ,2 ,0 ,20,1 ,-3,
#     -2,0 ,-4,3 ,0 ,-1,20,0 ,
#     0 ,-1,-1,-2,0 ,3 ,0 ,20),nrow=d,ncol=d,byrow=T)
# location=1:d
# alpha=seq(0.2,0.4,length.out = d)
# basefunc=rep("laplace",d)
# tpars=NULL


# ## model dimensions 10
# reps=100
# sampsize=10000
# d=10
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
#     0 ,-3,4 ,0 ,-2,0 ,-5,3 ,0 ,20),nrow=d,ncol=d,byrow=T)
# location=1:d
# alpha=seq(0.2,0.4,length.out = d)
# basefunc=rep("laplace",d)
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
dir.create(paste0("~/D",d,"S",sampsize))

result=foreach(i=1:reps,.packages=c('combinat','mrfDepth','nloptr',"optimx"),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 # required files
                 source(file="Code part 1 Functions.R")

                 # generating a sample
                 sample=Xsample(A=A,location=location,basefunc=basefunc,alpha=alpha,sampsize=sampsize,d=d,tpars=tpars,seed=seed-i+1)
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
                 filename=paste0("~/D",d,"S",sampsize,"/run",i,".csv")
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
    temp_data = t(read.csv(paste0("~/D",d,"S",sampsize,"/run",i,".csv"), row.names=1))
    
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
  colnames(dataset)=c(paste0(a,1:d),paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")),"log likelihood","time taken (s)")
} else {
  colnames(dataset)=c(paste0(a,1:d),paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")),paste0("df",indt),"log likelihood","time taken (s)")
}
save(dataset,file=paste0("~/D",d,"S",sampsize,".Rdata"))



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

Af=matrix(as.numeric(M[i,(2*d+1):(2*d+d^2)]),nrow=d)
af=as.numeric(M[i,1:d])
Af

# check the appropriate elements of the columns to put on the diagonal and put the indices in inds
inds=c(7,3,1,4,6,2,5,8)

Af[1:d,]=Af[inds,]
af=af[inds]
negdiag=which(diag(Af)<0)
Af[negdiag,]=-Af[negdiag,]
af[negdiag]=1-af[negdiag]
M[i,1:d]=af
M[i,(2*d+1):(2*d+d^2)]=Af
sorted$manual
# optional save for later use
save(M,file=paste0("~/D",d,"S",sampsize,"sorted.Rdata"))
