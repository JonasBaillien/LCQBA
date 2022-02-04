### function for fitting the linear combination of QBA-distributions model ###
##############################################################################

LCQBAFD_fit=function(data,basefunc,seed=NULL,maxiter=10^6,tol=10^-6,numstarts=20){
  # function for parameter estimation in a Linear Combinations of Quantile Based Asymmetric Family of Distributions
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
LCQBAFD_fit2=function(data,basefunc,seed=NULL,maxiter=10^8,tol=10^-3,maxstarts=100,reqfits=20){
  # function for parameter estimation in a Linear Combination of Quantile Based Asymmetric Family of Distributions
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

