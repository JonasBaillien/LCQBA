### Function for plotting the heatmap of a correlation or covariance matrix ###
###############################################################################

varcov_LCQBA=function(alpha,A,basefunc,tpars){
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

# ###############
# ### example ###
# ###############
# # model parameters
# alpha=c(0.35,0.7)
# mu=c(0,0)
# A=matrix(c(4,-3,1,4),nrow=2,ncol=2)
# basefunc=c("normal","laplace")
# tpars=c(NA,NA)


# # variance-covariance structure of the model
# varcov=varcov_LCQBA(alpha = alpha,A = A,basefunc = basefunc,tpars = tpars)
# varcov$covariance
# varcov$correlation
