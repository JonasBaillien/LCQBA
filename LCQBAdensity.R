### Function for calculating the density of a linear combination ###
### of QBA-distributions                                         ###
####################################################################

dLCQBA=function(x,alpha,mu,A,basefunc,tpars,LOG=FALSE){
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

