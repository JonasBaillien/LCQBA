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
  f=(x<=0)*(2*alpha*pnorm(q=(1-alpha)*x,mean=0,sd=1))+
    (x>0)*(2*alpha-1+2*(1-alpha)*pnorm(q=alpha*x,mean=0,sd=1))
  return(f)
}

# Two-piece laplace distribution with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecelaplacedistribution=function(x,alpha){
  f=(x<=0)*(2*alpha*plaplace(q=(1-alpha)*x,location=0,scale=1))+
    (x>0)*(2*alpha-1+2*(1-alpha)*plaplace(q=alpha*x,location=0,scale=1))
  return(f)
}

# Two-piece student-t distribution with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
# df=nu
twopiecestudentdistribution=function(x,alpha,nu){
  f=(x<=0)*(2*alpha*pt(q=(1-alpha)*x,df=nu))+
    (x>0)*(2*alpha-1+2*(1-alpha)*pt(q=alpha*x,df=nu))
  return(f)
}

# Two-piece logistic distribution with skewing parameter alpha in numeric vector x (location = 0, scale = 1)
twopiecelogisticdistribution=function(x,alpha){
  f=(x<=0)*(2*alpha*plogis(q=(1-alpha)*x,location=0,scale=1))+
    (x>0)*(2*alpha-1+2*(1-alpha)*plogis(q=alpha*x,location=0,scale=1))
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


###############
### Example ###
###############
# x=seq(-5,5,by=0.01)
# y=twopiecelogisticdensity(x,0.4)
# plot(x,y,type="l")

# x=seq(-5,5,by=0.01)
# y=twopiecelaplacedistribuiton(x,0.6)
# plot(x,y,type="l")
