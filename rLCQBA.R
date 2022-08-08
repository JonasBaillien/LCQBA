### Function for generating a sample from the linear combination of QBA-distributions ###
#########################################################################################

# generate sample from the untransformed, independent random variables
Zsample=function(sampsize,d,basefunc,alpha,seed=NULL,tpars=NULL){
  # input:
  # basefunc: a character vector of length d containing the underlying reference distribution.
  # (options are: "laplace", "normal", "logistic" and "t")
  # alpha: vector of lengthe d with values for skewing parameter (in [0,1])
  # sampsize: the sample size (number of samples needed)
  # d: the dimensionality of the data (=d>1)
  # seed: a seed for random number generation
  # tpars: a vector containing the degrees of freedom for possible student-t distributions
  
  # generating uniform samples
  if(!is.null(seed)){
    set.seed(seed)
  }
  uniformsample=matrix(runif(sampsize*d),nrow=sampsize,ncol=d)
  
  
  # obtaining sample from independent random variables
  Z=matrix(NA,nrow=sampsize,ncol=d)
  
  for(i in 1:d){
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
rLCQBA=function(A,mu,basefunc,alpha,sampsize,d,seed=NULL,tpars=NULL){
  # input:
  # A: dxd mixing matrix for affine transformation
  # mu: vector of length d with location of the mode 
  # basefunc: a character vector of length d containing the underlying reference distribution.
  # (options are: "laplace", "normal", "logistic" and "t")
  # alpha: vector of lengthe d with values for skewing parameter (in [0,1])
  # sampsize: the sample size (number of samples needed)
  # d: the dimensionality of the data (=d>1)
  # seed: a seed for random number generation
  # tpars: a vector containing the degrees of freedom for possible student-t distributions
  
  # check if d is strict larger than 1
  if(d<=1){
    stop("dimensions must be 2 or larger")
  }
  
  
  # check if matrix A has correct dimensions and is invertible
  if(sum(dim(A)==c(d,d))!=2){
    stop("A has wrong dimensions")
  }
  if(det(A)==0){
    stop("A is not invertible")
  }
  
  # check if location has correct length
  if(length(mu)!=d){
    stop("mu has incorrect length")
  }
  
  # check if correct number of basis function are given
  if(length(basefunc)!=d){
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
      df=as.matrix(rep(NA,d),nrow=1)
      df[index]=tpars
    }
  }
  
  # check if alpha has correct length and values
  testalpha=c(alpha>=0 & alpha<=1)
  if(sum(testalpha)!=d){
    stop("invalid or incorrect number of alpha parameters provided")
  }
  
  Z=Zsample(basefunc = basefunc, alpha = alpha, sampsize = sampsize,d = d, seed = seed, tpars = df)
  
  X=matrix(NA,nrow=sampsize,ncol=d)
  for(i in 1:sampsize){
    X[i,]=t(A)%*%Z[i,]+mu
  }
  return(list(Z,X))
}


# ###############
# ### Example ###
# ###############
# # model parameters
# alpha=c(0.35,0.7)
# mu=c(0,0)
# A=matrix(c(4,-3,1,4),nrow=2,ncol=2)
# basefunc=c("normal","laplace")
# tpars=c(NA,NA)


# # generate a sample
# samp=rLCQBA(alpha = alpha,mu = mu,A = A,basefunc = basefunc,
                 tpars = tpars,sampsize = 500,seed = 248)
# Z=samp$Z
# X=samp$X

# # visualize unmixed and final data
# plot(X,col="red",xlab="",ylab="")
# points(Z,col="black")
# legend("bottomleft",legend=c('Z','X'),col=c(1,2),pch=c(1,1))
