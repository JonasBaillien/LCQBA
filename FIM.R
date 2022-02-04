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
