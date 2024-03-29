###########################################################
### script to fit all possible models on a real dataset ###
###########################################################

# please select the correct data and change the directories where output is written, or choose your own data.
# the fitter can also be changed accordingly in the code below.
# take care in changing the directories of output and sourced files.


### load in the  data ###
#########################

# load ais data
library(DAAG)
dat=as.matrix(ais[,c(6,9)])



### matrix of all possible combinations of base functions ###
#############################################################
d=length(dat[1,])
bf=c("normal","t","laplace","logistic")
expandlist=list()
for(i in 1:d){
  expandlist[[i]]=bf
}
basefuncs=as.matrix(expand.grid(expandlist))
basefuncs=apply(basefuncs,1,sort)
basefuncs=t(basefuncs)
basefuncs=unique(basefuncs,margin=1)

nmods=length(basefuncs[,1])


### create folder to store output in (put correct name here) ###
################################################################
setwd(dir = choose.dir())
dir.create(paste0("~/name"))

# also change the file locations further below to write in this directory



### fitting the models to the data ###
######################################

# for the fitter
seed=126
maxiter=50000
tol=10^-5
maxstarts=50
reqfits=25

# main loop 
cores=detectCores()
cl <- makeCluster(4)
registerDoParallel(cl)

result=foreach(i=1:nmods,.packages=c('nloptr','optimx'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 # Required functions and correct path
                 source("names")
                 
                 # basefunctions of the model fitted to the data
                 basefunc=basefuncs[i,]
                 
                 # fitting the model to the data
                 fit=LCQBAFD_fit2(data=dat,basefunc=basefunc,seed=seed,maxiter=maxiter,tol=tol,maxstarts=maxstarts,reqfits=reqfits)

                 output=c(as.numeric(fit$`fitted parameters`),fit$`log likelihood fit`)
                 
                 # storing the output, put the right directory here (change the last file location)
                 filename=paste0("~/name/M",i,".csv")
                 write.csv(output,file=filename)
               }
stopCluster(cl)
x11()



### combining output in list containing matrices for alpha, mu, A, ###
### df, log-likelihood and base functions                          ###
######################################################################

alpha_fit=matrix(NA,nrow=nmods,ncol=d)
mu_fit=matrix(NA,nrow=nmods,ncol=d)
A_fit=matrix(NA,nrow=nmods,ncol=d^2)
df_fit=matrix(NA,nrow=nmods,ncol=d)
loglikelihood_fit=rep(NA,nmods)
basefunctions_fit=basefuncs

for (i in 1:nmods){
  # put the right directory here (change the last file location)
  temp_data = t(read.csv(paste0("~/name/M",i,".csv"), row.names=1))
  
  alpha_fit[i,] = as.numeric(temp_data[1:d])
  
  mu_fit[i,] = as.numeric(temp_data[(d+1):(2*d)])
  
  A_fit[i,] = as.numeric(temp_data[(2*d+1):(d^2+2*d)])
  
  if(is.element("t",basefuncs[i,])){
    indt=which(basefuncs[i,]=="t")
    nt=length(indt)
    df_fit[i,indt] = as.numeric(temp_data[(d^2+2*d+1):(d^2+2*d+nt)])
  } else {
    nt=0
  }
  
  loglikelihood_fit[i] = as.numeric(temp_data[d^2+2*d+nt+1])
}

# combine output
modelfit=list("alpha"=alpha_fit,"mu"=mu_fit,"A"=A_fit,"df"=df_fit,"logl"=loglikelihood_fit,"basefunc"=basefunctions_fit)

# save output in correct location (change directory name if needed)
save(modelfit,file="~/name/fittedmodels.Rdata")
