### libraries ###
#################
library(DAAG)
library(sn)
library(mrfDepth)
library(fields)
library(doParallel)
library(combinat)
library(nloptr)
library(optimx)
library(utils)
library(mvnTest)
library(mixtools)
library(LaplacesDemon)
library(QBAsyDist)
library(expm)
library(alabama)



### density functions of used 1D QBA-distributions ###
######################################################

# Two-piece normal density with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecenormaldensity=function(x,alpha){
  f=(x<=0)*alpha*(1-alpha)*sqrt(2/pi)*exp( -1/2*(1-alpha)^2*(x)^2)+
    (x>0)*alpha*(1-alpha)*sqrt(2/pi)*exp( -1/2*alpha^2*(x)^2)
  return(f)
}

# Two-piece laplace density with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecelaplacedensity=function(x,alpha){
  f=(x<=0)*alpha*(1-alpha)*exp( (1-alpha)*x)+
    (x>0)*alpha*(1-alpha)*exp( -alpha*x)
  return(f)
}

# Two-piece student-t density with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecestudentdensity=function(x,alpha,nu){
  f=(x<=0)*2*alpha*(1-alpha)/(sqrt(nu)*beta(1/2,nu/2))*(1+(1-alpha)^2/nu*(x)^2)^(-(nu+1)/2) +
    (x>0)*2*alpha*(1-alpha)/(sqrt(nu)*beta(1/2,nu/2))*(1+alpha^2/nu*(x)^2)^(-(nu+1)/2)
  return(f)
}

# Two-piece logistic density with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecelogisticdensity=function(x,alpha){
  f=(x<=0)*2*alpha*(1-alpha)*exp((1-alpha)*x)/(1+exp((1-alpha)*x))^2+
    (x>0)*2*alpha*(1-alpha)*exp(-alpha*x)/(1+exp(-alpha*x))^2
  return(f)
}


### Cumulative distribution functions ###
#########################################

# Two-piece normal distribution with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecenormaldistribution=function(x,alpha){
  f=(x<=0)*2*alpha*pnorm(q=(1-alpha)*x,mean=0,sd=1)+
    (x>0)*2*alpha-1+2*(1-alpha)*pnorm(q=alpha*x,mean=0,sd=1)
  return(f)
}

# Two-piece laplace distribution with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecelaplacedistribution=function(x,alpha){
  f=(x<=0)*2*alpha*plaplace(q=(1-alpha)*x,location=0,scale=1)+
    (x>0)*2*alpha-1+2*(1-alpha)*plaplace(q=alpha*x,location=0,scale=1)
  return(f)
}

# Two-piece student-t distribution with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
# df=nu
twopiecestudentdistribution=function(x,alpha,nu){
  f=(x<=0)*2*alpha*pt(q=(1-alpha)*x,df=nu)+
    (x>0)*2*alpha-1+2*(1-alpha)*pt(q=alpha*x,df=nu)
  return(f)
}

# Two-piece logistic distribution with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecelogisticdistribution=function(x,alpha){
  f=(x<=0)*2*alpha*plogis(q=(1-alpha)*x,location=0,scale=1)+
    (x>0)*2*alpha-1+2*(1-alpha)*plogis(q=alpha*x,location=0,scale=1)
  return(f)
}



### Quantile functions ###
##########################


# Two-piece normal beta quantile (numeric vector) with skewing parameter alpha (location = 0, scale = 1)
twopiecenormalquantile=function(beta,alpha){
  f=rep(NA,length(beta))
  for(i in 1:length(beta)){
    if(beta[i]<=alpha){
      f[i]=1/(1-alpha)*qnorm(p=beta[i]/(2*alpha),mean=0,sd=1)
    } else {
      f[i]=1/alpha*qnorm(p=(1+beta[i]-2*alpha)/(2*(1-alpha)),mean=0,sd=1)
    }
  }
  return(f)
}

# Two-piece laplace beta quantile (numeric vector) with skewing parameter alpha (location = 0, scale = 1)
twopiecelaplacequantile=function(beta,alpha){
  f=rep(NA,length(beta))
  for(i in 1:length(beta)){
    if(beta[i]<=alpha){
      f[i]=1/(1-alpha)*qlaplace(p=(beta[i]/(2*alpha)),location=0,scale=1)
    } else {
      f[i]=1/alpha*qlaplace(p=(1+beta[i]-2*alpha)/(2*(1-alpha)),location=0,scale=1)
    }
  }
  return(f)
}

# Two-piece student-t beta quantile (numeric vector) with skewing parameter alpha (location = 0, scale = 1)
# df=nu
twopiecestudentquantile=function(beta,alpha,nu){
  f=rep(NA,length(beta))
  for(i in 1:length(beta)){
    if(beta[i]<=alpha){
      f[i]=1/(1-alpha)*qt(p=beta[i]/(2*alpha),df=nu)
    } else {
      f[i]=1/alpha*qt(p=(1+beta[i]-2*alpha)/(2*(1-alpha)),df=nu)
    }
  }
  return(f)
}

# Two-piece logistic beta quantile (numeric vector) with skewing parameter alpha (location = 0, scale = 1)
twopiecelogisticquantile=function(beta,alpha){
  f=rep(NA,length(beta))
  for(i in 1:length(beta)){
    if(beta[i]<=alpha){
      f[i]=1/(1-alpha)*qlogis(p=beta[i]/(2*alpha),location=0,scale=1)
    } else {
      f[i]=1/alpha*qlogis(p=(1+beta[i]-2*alpha)/(2*(1-alpha)),location=0,scale=1)
    }
  }
  return(f)
}


### f'()/f() ###
################

# Two-piece normal density with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecenormalderiv=function(x,alpha){
  f=-(x<=0)*(1-alpha)^2*x+(x>0)*alpha^2*x
  return(f)
}

# Two-piece laplace density with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecelaplacederiv=function(x,alpha){
  f=(x<=0)*(1-alpha)/2-(x>0)*alpha/2
  return(f)
}

# Two-piece student-t density with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecestudentderiv=function(x,alpha,nu){
  f=(x<=0)*((nu+1)/nu*(1-alpha)*x)/(1+(1-alpha)^2*x^2/nu)+
    (x>0)*(-(nu+1)/nu*alpha*x)/(1+alpha^2*x^2/nu)
  return(f)
}

# Two-piece logistic density with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecelogisticderiv=function(x,alpha){
  f=(x<=0)*(1-alpha)*(1-exp((1-alpha)*x))/(1+exp((1-alpha)*x))
  -(x>0)*alpha*(1-exp(-alpha*x))/(1+exp(-alpha*x))
  return(f)
}





### Function for calculating the density of a linear combination ###
### of QBA-distributions                                         ###
####################################################################

dMQBA=function(x,alpha,mu,A,basefunc,tpars,LOG=FALSE){
  # x is a nxd matrix containing the points in which to evaluate the density
  # alpha is a d-vector containing the skewing parameters between 0 an 1
  # mu is a d-vector containing the location parameter mu_a
  # A is a dxd matrix containing the mixing matrix
  # basefunc is a character vector with the reference densities, can be "normal", "t", "logistic" or "laplace"
  # tpars is a d-vector containing the degrees of freedom for any basefunc=="t" on the correct indices
  # LOG is a logical whether the log-density should be returned or no
  
  
  
  # putting data in correct format in case it is a vector
  x=as.matrix(x)
  if(ncol(x)==1){
    x=t(x)
  }
  
  # calculation of density
  d=length(x[1,])
  B=solve(matrix(A,nrow=d))
  Y=sweep(x,2,mu)%*%B
  
  dens=abs(det(B))
  
  for(i in 1:d){
    dens=dens*switch(basefunc[i],
                     "normal"=QBAsyDist::dAND(y = Y[,i],mu = 0,phi = 1,alpha = alpha[i]),
                     "laplace"=QBAsyDist::dALaD(y = Y[,i],mu = 0,phi = 1,alpha = alpha[i]),
                     "logistic"=QBAsyDist::dALoD(y = Y[,i],mu = 0,phi = 1,alpha = alpha[i]),
                     "t"=QBAsyDist::dATD(y = Y[,i],mu = 0,phi = 1,alpha = alpha[i],nu = tpars[i])
    )
  }
  
  if(LOG==T){
    return(log(dens))
  }
  
  return(dens)
}






### function for fitting the linear combination of QBA-distributions model ###
##############################################################################

MVQBAFD_fit=function(data,basefunc,seed=NULL,maxiter=10^6,tol=10^-6,numstarts=20){
  # function for parameter estimation in a MultiVariate Quantile Based Asymmetric Family of Distributions
  # model
  
  # data is a nxd-matrix which contains n-observations of d-variate points which form the data
  
  # basefunc is a vector  of lengthe d which containing the names of univariate functions which 
  # are linearly combined to form the data. choises are "normal", "laplace", "logistic" or "t"
  
  # seed is a numeric seed to ensure consistents results on the same data
  
  # maxiter is a numeric value for the underlying optimizer which sets the maximum number of iterations
  # for the optimizer
  
  # tol is a numeric value for the underlying optimizer which sets the relative accuracy of the 
  # optimization step
  
  # numstarts is a numeric value which determines the number of different starting points for the optimizer
  
  X=na.omit(data)
  
  # dimensions of X
  n=length(X[,1])
  d=length(X[1,])
  
  # check if correct number of base functions is supplied
  if(length(basefunc)!=d){
    stop("incorrect number of base functions supplied")
  }
  
  # check if basefuncs are ok
  possibilities=c("t","normal","logistic","laplace")
  test=!is.element(basefunc,possibilities)
  if(sum(test)>0){
    stop("incorrect basis function detected")
  }
  
  
  
  # generate bounds on the parmaters (alpha between 0 and 1; mu in the range of the data)
  
  # lower bounds for parameters in optimization
  loweralpha=rep(0,d)
  lowermu=apply(X,2,min)
  lowerA=matrix(-max(sd(t(X))),nrow=d,ncol=d)
  lowerbounds=as.vector(cbind(loweralpha,lowermu,lowerA))
  
  
  # upper bounds for parameters in optimization
  upperalpha=rep(1,d)
  uppermu=apply(X,2,max)
  upperA=matrix(max(sd(t(X))),nrow=d,ncol=d)
  upperbounds=as.vector(cbind(upperalpha,uppermu,upperA))
  
  
  
  
  # generate starting values
  
  if(!is.null(seed)){
    set.seed(seed)
  } else {
    seed=sample(1:10^6,1)
    set.seed(seed)
  }
  
  
  # for alpha
  startalpha=matrix(runif(d*numstarts),nrow=numstarts,ncol=d,byrow=T)
  # for mu
  startmu=matrix(runif(d*numstarts,apply(X,2,min),apply(X,2,max)),nrow=numstarts,ncol=d,byrow=T)
  # for A 
  startA=matrix(runif(numstarts*d^2,min=-max(sd(t(X))),max=max(sd(t(X)))),nrow=numstarts,ncol=d^2,byrow=T)
  
  
  # combine them in matrix
  x0=cbind(startalpha,startmu,startA)
  
  
  # degrees of freedom if student t distribution is passed down
  
  if(is.element("t",basefunc)){
    
    # which basis functions are t-distributed
    indt=which(basefunc=="t")
    ldf=length(indt)
    
    # generate starting values for degrees of freedom as well as bounds
    lowerdf=rep(1,ldf)
    upperdf=rep(10000,ldf)
    startdf=matrix(runif(numstarts*ldf,min=2,max=100),ncol=ldf,nrow=numstarts,byrow=T)
    
    # pasting these to other starting values
    x0=cbind(x0,startdf)
    lowerbounds=c(lowerbounds,lowerdf)
    upperbounds=c(upperbounds,upperdf)
    
  } else {
    # in case no t-distribution is present
    indt=NULL
  }
  
  
  # holding matrix for best found parameters (one row is one set)
  paramsfit=data.frame(matrix(NA,ncol=length(x0[1,]),nrow=numstarts))
  a=intToUtf8(945) # character for alpha
  m=intToUtf8(956) # character for mu
  # tau=intToUtf8(964)
  if(is.null(indt)){
    colnames(paramsfit)=c(paste0(a,1:d),paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")))
  } else {
    colnames(paramsfit)=c(paste0(a,1:d),paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")),paste0("df",indt))
  }
  
  
  # vector with log-likelihood in found optimum 
  loglfit=rep(NA,numstarts)
  
  
  
  # main loop for parameter estimation using the bobyqa function from the nloptr package
  for(i in 1:numstarts){
    
    try({
      # set seed
      set.seed(seed+i)
      # starting values for parameters
      parstart=x0[i,]
      
      # optimization
      output=nloptr(x0=parstart,eval_f=opt_fun,lb=lowerbounds,ub=upperbounds,opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=maxiter,xtol_rel=tol,ftol_abs=tol,print_level=0),X=X,basefunc=basefunc,indt=indt)
      paramsfit[i,]=output$solution
      loglfit[i]=output$objective
    },silent=T)
  }
  
  # only return starting value which provides best fit (lowest minus log-likelihood)
  indmin=which.min(loglfit)
  paramsfit=paramsfit[indmin,]
  loglfit=-loglfit[indmin]
  
  return(list("fitted parameters"=paramsfit,"log likelihood fit"=loglfit))
}



# minus log-likelihood function which is to be minimized
opt_fun=function(parms,X,basefunc,indt=NULL){
  # parms is a vector containing the parameter values. In order these are: alpha, mu, A (by column) and 
  # possibly degrees of freedom for student-t distributions. 
  
  # X is a nxd matrix containing the data
  
  # basefunc is a vector  of lengthe d which containing the names of univariate functions which 
  # are linearly combined to form the data. choises are "normal", "laplace", "logistic" or "t". 
  # indices of t-distributions within this vector are supplied via indt. This is used to generate
  # a vector of length d which contains the degrees of freedom with the corresponding t-distribution in
  # in basefunc
  
  # dimensions
  n=length(X[,1])
  d=length(X[1,])
  
  # parameter values
  alpha=parms[1:d]
  mu=parms[(d+1):(2*d)]
  A=matrix(parms[(2*d+1):(d^2+2*d)],nrow=d)
  B=solve(A)
  
  # creating vector with degrees of freedom
  if(!is.null(indt)){
    tpars=rep(NA,d)
    tpars[indt]=parms[(d^2+2*d+1):(d^2+2*d+length(indt))]
  }
  
  # X*A^{-1}
  V=X%*%B
  # mu*A^{-1}
  VV=mu%*%B
  
  # calculation of minus log-likelihood, max ensures that no inf's of NA's are produced
  value=-n*log(abs(det(B)))
  for(i in 1:d){
    if(basefunc[i]=="laplace"){
      value=value-sum(log(twopiecelaplacedensity(V[,i]-VV[i],alpha[i])))
    } else if(basefunc[i]=="normal"){
      value=value-sum(log(twopiecenormaldensity(V[,i]-VV[i],alpha[i])))
    } else if(basefunc[i]=="logistic"){
      value=value-sum(log(twopiecelogisticdensity(V[,i]-VV[i],alpha[i])))
    } else {
      value=value-sum(log(twopiecestudentdensity(V[,i]-VV[i],alpha[i],nu=tpars[i])))
    }
  }
  
  return(value)
}


### alternative function when a set number of good fits needs to be obtained
MVQBAFD_fit2=function(data,basefunc,seed=NULL,maxiter=10^8,tol=10^-3,maxstarts=100,reqfits=20){
  # function for parameter estimation in a MultiVariate Quantile Based Asymmetric Family of Distributions
  # model
  
  # data is a nxd-matrix which contains n-observations of d-variate points which form the data
  
  # basefunc is a vector  of lengthe d which containing the names of univariate functions which 
  # are linearly combined to form the data. choises are "normal", "laplace", "logistic" or "t"
  
  # seed is a numeric seed to ensure consistents results on the same data
  
  # maxiter is a numeric value for the underlying optimizer which sets the maximum number of iterations
  # for the optimizer
  
  # tol is a numeric value for the underlying optimizer which sets the relative accuracy of the 
  # optimization step
  
  # maxstarts is a numeric value which determines the maximum number of different starting points for the optimizer
  
  # reqfits is a numeric value which determines the number of valid fits to be obtained, the final fit will be
  # the best one of these
  
  
  X=na.omit(data)
  
  # dimensions of X
  n=length(X[,1])
  d=length(X[1,])
  
  # check if correct number of base functions is supplied
  if(length(basefunc)!=d){
    stop("incorrect number of base functions supplied")
  }
  
  # check if basefuncs are ok
  possibilities=c("t","normal","logistic","laplace")
  test=!is.element(basefunc,possibilities)
  if(sum(test)>0){
    stop("incorrect basis function detected")
  }
  
  
  
  # generate bounds on the parmaters (alpha between 0 and 1; mu in the range of the data)
  
  # lower bounds for parameters in optimization
  loweralpha=rep(0,d)
  lowermu=apply(X,2,min)
  lowerA=matrix(-max(apply(X,2,sd)),nrow=d,ncol=d)
  
  lowerbounds=as.vector(cbind(loweralpha,lowermu,lowerA))
  
  
  # upper bounds for parameters in optimization
  upperalpha=rep(1,d)
  uppermu=apply(X,2,max)
  upperA=matrix(max(apply(X,2,sd)),nrow=d,ncol=d)
  
  upperbounds=as.vector(cbind(upperalpha,uppermu,upperA))
  
  
  
  
  # generate starting values
  
  if(!is.null(seed)){
    set.seed(seed)
  } else {
    seed=sample(1:10^6,1)
    set.seed(seed)
  }
  
  # degrees of freedom if student t distribution is passed down
  
  if(is.element("t",basefunc)){
    
    # which basis functions are t-distributed
    indt=which(basefunc=="t")
    ldf=length(indt)
    
    # generate starting values for degrees of freedom as well as bounds
    lowerdf=rep(2,ldf)
    upperdf=rep(1000,ldf)
    
    # pasting these to other starting values
    lowerbounds=c(lowerbounds,lowerdf)
    upperbounds=c(upperbounds,upperdf)
    
  } else {
    # in case no t-distribution is present
    indt=NULL
    upperdf=NULL
    ldf=0
  }
  
  
  # reference number of iterations required when no convergence is achieved
  output=nloptr(x0=c(rep(0.05,d),lowermu,upperA-runif(d^2,min=upperA/20,max=upperA/10),upperdf/2),eval_f=opt_fun,lb=lowerbounds,ub=upperbounds,opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=maxiter,xtol_rel=tol,ftol_abs=tol,print_level=0),X=X,basefunc=basefunc,indt=indt)
  refiters=output$iterations
  
  j=1
  i=0
  
  paramsfit=matrix(NA,nrow=reqfits,ncol=2*d+d^2+ldf)
  loglfit=rep(NA,reqfits)
  a=intToUtf8(945) # character for alpha
  m=intToUtf8(956) # character for mu
  # tau=intToUtf8(964)
  if(is.null(indt)){
    colnames(paramsfit)=c(paste0(a,1:d),paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")))
  } else {
    colnames(paramsfit)=c(paste0(a,1:d),paste0(m,1:d),paste0("A",apply(expand.grid(1:d, 1:d), 1, paste, collapse="")),paste0("df",indt))
  }
  
  
  # main loop for parameter estimation using the bobyqa function from the nloptr package
  while(i<maxstarts & j<=reqfits){
    
    try({
      # set seed
      set.seed(seed+i)
      
      # starting values for parameters
      
      # for alpha
      startalpha=runif(d)
      # for mu
      startmu=runif(d,apply(X,2,min),apply(X,2,max))
      # for A (ensure positive diagonal elements)
      startA=runif(d^2,min=-max(apply(X,2,sd)),max=max(apply(X,2,sd)))
      startA[seq(1,d^2,by=d+1)]=abs(startA[seq(1,d^2,by=d+1)])
      # for df
      startdf=runif(ldf,min=2,max=100)
      
      parstart=c(startalpha,startmu,startA,startdf)
      
      
      # optimization
      output=nloptr(x0=parstart,eval_f=opt_fun,lb=lowerbounds,ub=upperbounds,opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=maxiter,xtol_rel=tol,ftol_abs=tol,print_level=0),X=X,basefunc=basefunc,indt=indt)
      iters=output$iterations
      if(iters>refiters*1.3){
        paramsfit[j,]=output$solution
        loglfit[j]=output$objective
        j=j+1
      } 
    },silent=T)
    i=i+1
  }
  
  indmin=which.min(loglfit)
  paramsfit=paramsfit[indmin,]
  loglfit=-loglfit[indmin]
  
  return(list("fitted parameters"=paramsfit,"log likelihood fit"=loglfit))
}






### Functions to calculate the fisher information matrix ###
############################################################

### helper functions
# function to calculate kappa_{k,1}
expectation=function(alpha,basefunc,tpars=NULL){
  d=length(alpha)
  m=rep(NA,d)
  for(i in 1:d){
    m[i]=switch(basefunc[i],normal=sqrt(2/pi),
                laplace=1,
                logistic=2*log(2),
                t=(sqrt(tpars[i])*gamma((tpars[i]-1)/2))/(sqrt(pi)*gamma(tpars[i]/2))   )
  }
  out=(1-2*alpha)/(alpha*(1-alpha))*m
  return(out)
}

# function to calculate kappa_{k,2}
expectationsquare=function(alpha,basefunc,tpars=NULL){
  d=length(alpha)
  m=rep(NA,d)
  for(i in 1:d){
    m[i]=switch(basefunc[i],normal=1,
                laplace=2,
                logistic=pi^2/3,
                t=tpars[i]/(tpars[i]-2)   )
  }
  out=((1-alpha)^3+alpha^3)/(alpha^2*(1-alpha)^2)*m
  return(out)
}

# function to calculate gamma_{k,1}
gamma1=function(d,basefunc,tpars=NULL){
  g=rep(NA,d)
  for(i in 1:d){
    g[i]=switch(basefunc[i],normal=1/2,
                laplace=1/2,
                logistic=1/6,
                t=(tpars[i]+1)/(2*(tpars[i]+3))   )
  }
  return(g)
}

# function to calculate gamma_{k,2}
gamma2=function(d,basefunc,tpars=NULL){
  g=rep(NA,d)
  for(i in 1:d){
    g[i]=switch(basefunc[i],normal=sqrt(2/pi),
                laplace=1/2,
                logistic=1/6+log(2)/3,
                t=4*gamma(tpars[i]/2+3/2)/(sqrt(pi*tpars[i])*(tpars[i]+3)*gamma(tpars[i]/2))   )
  }
  return(g)
}

# function to calculate gamma_{k,3}
gamma3=function(d,basefunc,tpars=NULL){
  g=rep(NA,d)
  for(i in 1:d){
    g[i]=switch(basefunc[i],normal=3/2,
                laplace=1,
                logistic=2/3+pi^2/18,
                t=3*(tpars[i]+1)/(2*(tpars[i]+3))   )
  }
  return(g)
}

# function to calculate B_{n,j}^{k,l}
Bn=function(A,n,j,k,l){
  # according to the formula for B_n
  # A is the mixing matrix
  # j,k,l are the indices to be used in the function
  # k and l are the row resp. column index of the element of A to which the derivative is taken
  # j must be different from k and is the index of the basefunction
  
  # j must be different from k
  if(j==k){
    stop("j must be different from k")
  }
  
  # range of the sum
  d=length(A[1,])
  
  # calculation
  if(d==2){
    for(i in (1:d)[-l]){
      out=(-1)^(i+j+k+l+(j>k)+(i>l))*A[n,i]/det(A)
    }
  } else {
    out=0
    for(i in (1:d)[-l]){
      out=out+(-1)^(i+j+k+l+(j>k)+(i>l))*A[n,i]*det(matrix(A[-c(j,k),-c(i,l)],ncol=d-2,nrow=d-2))/det(A)
    }
  }
  
  # output
  return(out)
}



### main function for theoretical asymptotic Fisher information matrix 
FIM=function(alpha,mu,A,basefunc,tpars=NULL){
  # alpha is a numeric vector of length d containing the skewing parameters
  
  # mu is a numeric vector of length d containing the location parameters
  
  # A is vector of length d^2 or a dxd-matrix containing the mixing parameters
  
  # basefunc is a vector of lenght d containing the names of the basefunctions (options are "normal",
  # "laplace","logistic" or "t")
  
  # tpars is a vector containing the degrees of freedom for possible t-distributions
  
  
  ### checks and dimensions ###
  #############################
  
  # which basefunctions are t-distributed
  indt=which(basefunc=="t")
  
  # dimensions
  d=length(alpha)
  end=2*d+d^2+length(indt)
  
  # formatting A to correct type
  A=matrix(A,nrow=d,ncol=d)
  
  # checks for correct input
  if(length(basefunc)!=d){
    stop("incorrect number of basefunctions supplied")
  }
  
  possibilities=c("t","normal","logistic","laplace")
  test=!is.element(basefunc,possibilities)
  if(sum(test)>0){
    stop("incorrect basis function detected")
  }
  
  if(length(mu)!=d){
    stop("incorrect dimensions for mu")
  }
  
  if(length(A[,1])!=d | length(A[1,])!=d){
    stop("incorrect dimensions for A")
  }
  
  if(length(tpars)!=length(indt)){
    stop("incorrect number of degrees of freedom supplied")
  }
  
  if(det(A)==0){
    stop("A is singular")
  }
  
  if(sum(diag(A)<=0)!=0){
    stop("Diagonal elements of A must be positive")
  }
  
  # formatting tpars to correct type
  df=rep(NA,d)
  df[indt]=tpars
  
  
  ### used quantities in the FI ###
  #################################
  
  # kappas
  K1=expectation(alpha=alpha,basefunc=basefunc,tpars=df)
  K2=expectationsquare(alpha=alpha,basefunc=basefunc,tpars=df)
  
  # gammas
  g1=gamma1(d=d,basefunc=basefunc,tpars=df)
  g2=gamma2(d=d,basefunc=basefunc,tpars=df)
  g3=gamma3(d=d,basefunc=basefunc,tpars=df)
  
  # inverse of A
  Ainv=solve(A)
  
  # B_{n,j}^{k,l}'s
  B=array(NA,dim = rep(d,4))
  for(n in 1:d){
    for(j in 1:d){
      for(k in (1:d)[-j]){
        for(l in 1:d){
          B[n,j,k,l]=Bn(A = A,n = n,j = j,k = k,l = l)
        }
      }
    }
  }
  
  
  
  ### constructing the FI matrix ###
  ##################################
  
  # note: index letters are the same as used in the paper to avoid confusion
  # elements of A are ordered column by column, so A_11, A_21, A31, ....
  
  
  # empty FI matrix with correct dimensions
  I=matrix(0,nrow=end,ncol=end)
  
  
  # Variances-covariances for alpha-alpha
  diag(I[1:d,1:d])=(2*(alpha^3+(1-alpha)^3)*g3-(1-2*alpha)^2)/(alpha^2*(1-alpha)^2)
  
  # Variances-covariances for alpha-mu
  I[1:d,(d+1):(2*d)]=-2*t(Ainv)*g2
  I[(d+1):(2*d),1:d]=t(I[1:d,(d+1):(2*d)])
  
  # Variances-covariances for alpha-A
  temp=matrix(0,nrow=d,ncol=d^2)
  for(k in 1:d){
    # if k=l
    for(m in 1:d){
      temp[k,k+(m-1)*d]=(-1)^(k+m+1)*(1-2*alpha[k])*(2*g3[k]-1)*det(matrix(A[-k,-m],nrow=d-1,ncol=d-1))/(alpha[k]*(1-alpha[k])*det(A))
    }
    
    # if k!=l
    for(l in (1:d)[-k]){
      for(m in 1:d){
        tempvalue=0
        for(n in (1:d)[-k]){
          tempvalue=tempvalue+B[n,k,l,m]*K1[n]
        }
        tempvalue=tempvalue*2*g2[k]
        temp[k,(m-1)*d+l]=tempvalue
      }
    }
  }
  
  I[1:d,(2*d+1):(2*d+d^2)]=temp
  I[(2*d+1):(2*d+d^2),1:d]=t(temp)
  
  
  # Variances-covariances for mu-mu
  temp=matrix(0,nrow=d,ncol=d)
  for(k in 1:d){
    for(l in 1:d){
      tempvalue=0
      for(j in 1:d){
        tempvalue=tempvalue+2*alpha[j]*(1-alpha[j])*Ainv[k,j]*Ainv[l,j]*g1[j]
      }
      temp[k,l]=tempvalue
    }
  }
  I[(d+1):(2*d),(d+1):(2*d)]=temp
  
  
  # Variances-covariances for mu-A
  temp=matrix(0,nrow=d,ncol=d^2)
  for(k in 1:d){
    for(l in 1:d){
      for(m in 1:d){
        tempvalue1=0
        for(j in (1:d)[-l]){
          tempvalue2=0
          for(n in (1:d)[-j]){
            tempvalue2=tempvalue2+B[n,j,l,m]*K1[n]
          }
          tempvalue2=alpha[j]*(1-alpha[j])*Ainv[k,j]*g1[j]*tempvalue2
          tempvalue1=tempvalue1+tempvalue2
        }
        temp[k,(m-1)*d+l]=-2*tempvalue1
      }
    }
  }
  
  I[(d+1):(2*d),(2*d+1):(2*d+d^2)]=temp
  I[(2*d+1):(2*d+d^2),(d+1):(2*d)]=t(temp)
  
  # Variances-covariances for A-A
  temp=matrix(0,nrow=d^2,ncol=d^2)
  for(k in 1:d){
    for(l in 1:d){
      for(r in 1:d){
        for(s in 1:d){
          tempvalue1=0
          for(j in (1:d)[-c(k,r)]){
            tempvalue2=0
            for(m in (1:d)[-j]){
              for(n in (1:d)[-c(m,j)]){
                tempvalue2=tempvalue2+B[m,j,k,l]*K1[m]*B[n,j,r,s]*K1[n]
              }
            }
            for(g in (1:d)[-j]){
              tempvalue2=tempvalue2+B[g,j,k,l]*B[g,j,r,s]*K2[g]
            }
            tempvalue2=2*alpha[j]*(1-alpha[j])*g1[j]*tempvalue2
            tempvalue1=tempvalue1+tempvalue2
          }
          tempvalue3=0
          for(q in (1:d)[-k]){
            for(j in (1:d)[-c(q,r)]){
              tempvalue3=tempvalue3+B[q,j,r,s]*B[j,q,k,l]
            }
          }
          temp[(l-1)*d+k,(s-1)*d+r]=tempvalue1+tempvalue3
        }
      }
    }
  }
  for(k in 1:d){
    for(l in 1:d){
      for(s in 1:d){
        temp[(l-1)*d+k,(s-1)*d+k]=temp[(l-1)*d+k,(s-1)*d+k]+(-1)^(l+s)*det(matrix(A[-k,-l],nrow=d-1,ncol=d-1))*det(matrix(A[-k,-s],nrow=d-1,ncol=d-1))*(2*g3[k]-1)/(det(A))^2
      }
      
    }
  }
  I[(2*d+1):(2*d+d^2),(2*d+1):(2*d+d^2)]=temp
  
  
  
  # for degrees of freedom if present
  if(is.element("t",basefunc)){
    # Variances-covariances for alpha-df
    temp=matrix(0,nrow=d,ncol=length(indt))
    for(i in 1:length(indt)){
      for(j in 1:d){
        if(indt[i]==j){
          temp[j,i]=2*(1-2*alpha[j])/(alpha[j]*(1-alpha[j])*(df[j]+1)*(df[j]+3))
        }
      }
    }
    I[1:d,(2*d+d^2+1):end]=temp
    I[(2*d+d^2+1):end,1:d]=t(temp)
    
    # Variances-covariances for mu-df are 0
    
    # Variances-covariances for A-df
    temp=matrix(0,nrow=d^2,ncol=length(indt))
    for(i in 1:length(indt)){
      for(j in 1:d){
        for(l in 1:d){
          if(indt[i]==j){
            temp[(l-1)*d+j,i]=2*(-1)^(j+l+1)*det(matrix(A[-j,-l],nrow=d-1,ncol=d-1))/((df[j]+1)*(df[j]+3)*det(A))
            
          }
        }
      }
    }
    
    I[(2*d+1):(2*d+d^2),(2*d+d^2+1):end]=temp
    I[(2*d+d^2+1):end,(2*d+1):(2*d+d^2)]=t(temp)
    
    # Variances-covariances for df-df
    temp=matrix(0,nrow=length(indt),ncol=length(indt))
    diag(temp)=na.omit(1/4*(trigamma(df/2)-trigamma((df+1)/2)) - (df+5)/(2*df*(df+1)*(df+3)))
    I[(2*d+d^2+1):end,(2*d+d^2+1):end]=temp
  }
  
  return(I)
}




### Function for generating a sample from the linear combination of QBA-distributions ###
#########################################################################################

# generate sample from the untransformed, independent random variables
Zsample=function(sampsize,dims,basefunc,alpha,seed=NULL,tpars=NULL){
  
  # generating uniform samples
  if(!is.null(seed)){
    set.seed(seed)
  }
  uniformsample=matrix(runif(sampsize*dims),nrow=sampsize,ncol=dims)
  
  
  # obtaining sample from independent random variables
  Z=matrix(NA,nrow=sampsize,ncol=dims)
  
  for(i in 1:dims){
    if(basefunc[i]=="laplace"){
      for(j in 1:sampsize){
        Z[j,i]=twopiecelaplacequantile(uniformsample[j,i],alpha[i])
      }
    } else if(basefunc[i]=="normal"){
      for(j in 1:sampsize){
        Z[j,i]=twopiecenormalquantile(uniformsample[j,i],alpha[i])
      }
    } else if(basefunc[i]=="logistic"){
      for(j in 1:sampsize){
        Z[j,i]=twopiecelogisticquantile(uniformsample[j,i],alpha[i])
      }
    } else {
      for(j in 1:sampsize){
        Z[j,i]=twopiecestudentquantile(uniformsample[j,i],alpha[i],tpars[i])
      }
    } 
  }
  return(Z)
}

# generate a sample from the target distribution
Xsample=function(A,location,basefunc,alpha,sampsize,dims,seed=NULL,tpars=NULL){
  # input:
  # A: dxd mixing matrix for affine transformation
  # location: vector of lengthe d with location of the mode 
  # basefunc: a character vector of lengthe d containing the underlying basis distrubution.
  # (options are: "laplace", "normal", "logistic" and "t")
  # alpha: vector of lengthe d with values for skewing parameter (in [0,1])
  # sampsize: the sample size (number of samples needed)
  # dims: the dimensionality of the data (=d>1)
  # seed: a seed for random number generation
  # tpars: a vector containing the degrees of freedom for possible student-t distributions
  
  # check if dims is strict larger than 1
  if(dims<=1){
    stop("dimensions must be 2 or larger")
  }
  
  
  # check if matrix A has correct dimensions and is invertible
  if(sum(dim(A)==c(dims,dims))!=2){
    stop("A has wrong dimensions")
  }
  if(det(A)==0){
    stop("A is not invertible")
  }
  
  # check if location has correct length
  if(length(location)!=dims){
    stop("location vector has incorrect length")
  }
  
  # check if correct number of basis function are given
  if(length(basefunc)!=dims){
    stop("dimensions and number of basis function are not equal")
  }
  
  # check if basefuncs are ok
  possibilities=c("t","normal","logistic","laplace")
  test=!is.element(basefunc,possibilities)
  if(sum(test)>0){
    stop("incorrect basis function detected")
  }
  
  # check if correct amount of different df are given
  tpars=na.omit(tpars)
  if(is.element("t",basefunc)){
    if(is.null(tpars) | length(tpars)!=sum(basefunc=="t")){
      stop("no/incorrect number of degrees of freedom for student-t distribution provided")
    } else { # create helpvector with degrees of freedom
      index=which(basefunc=="t")
      df=as.matrix(rep(NA,dims),nrow=1)
      df[index]=tpars
    }
  }
  
  # check if alpha has correct length and values
  testalpha=c(alpha>=0 & alpha<=1)
  if(sum(testalpha)!=dims){
    stop("invalid or incorrect number of alpha parameters provided")
  }
  
  Z=Zsample(basefunc = basefunc, alpha = alpha, sampsize = sampsize,dims = dims, seed = seed, tpars = df)
  
  X=matrix(NA,nrow=sampsize,ncol=dims)
  for(i in 1:sampsize){
    X[i,]=t(A)%*%Z[i,]+location
  }
  return(list(Z,X))
}
















### Function for plotting the heatmap of a correlation or covariance matrix ###
###############################################################################

varcov_MQBA=function(alpha,A,basefunc,tpars){
  # function for the plotting of the heatmap of the correlation matrix of a
  # multivariate quantile based model with mixing matrix A (numeric dxd matrix),
  # skewing parameters alpha (numeric vector of length d), base functions
  # basefunc (character vector of length d) and tpars a vector with the degrees
  # of freedom for possible student-t base functions
  
  
  # dimension of the problem
  d=length(alpha)
  
  if(length(A[1,])!=d | length(A[,1])!=d | length(basefunc)!=d){
    stop("incorrect dimensions")
  }
  
  # check for correct basefunctions
  possibilities=c("t","normal","logistic","laplace")
  test=!is.element(basefunc,possibilities)
  if(sum(test)>0){
    stop("incorrect basis function detected")
  }
  
  indt=which(basefunc=="t")
  if(is.null(indt)){
    df=NULL
  } else {
    df=rep(NA,d)
    df[indt]=stats::na.omit(tpars)
  }
  
  # mu1 for the contruction of diag(var(Z_i))
  mu1=rep(NA,d)
  for(i in 1:d){
    mu1[i]=switch(basefunc[i],normal=sqrt(2/pi),
                  laplace=1,
                  logistic=2*log(2),
                  t=(sqrt(df[i])*gamma((df[i]-1)/2))/(sqrt(pi)*gamma(df[i]/2))
    )
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
  covmat=t(A)%*%varZ%*%A
  
  # standard deviations used to obtain the correlation
  sdmat=sqrt(diag(covmat))%*%t(sqrt(diag(covmat)))
  
  # correlation matrix of data according to fitted model
  cormat=covmat/sdmat
  
  
  # row and colnames
  rownames(cormat)=paste0("X",1:d)
  colnames(cormat)=paste0("X",1:d)
  
  
  return(list("correlation"=cormat,"covariance"=covmat))
}