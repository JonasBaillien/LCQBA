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
