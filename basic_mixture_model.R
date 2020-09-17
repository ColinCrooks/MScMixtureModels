
niter<-100

X<-as.matrix(cbind(final_data$ic_beta,final_data$gw_beta))
psi<-as.matrix(cbind(final_data$ic_se^2,final_data$gw_se^2))


# SNP level f(x | m )
fx<-function(X_,psi_,m,position)
{
  px<-1
  for(j in c(1,2))
  {
    px=px*exp((-(X_-m)^2)/(2*psi_))/sqrt(2*pi*psi_)
  }
  return(px)
}

# SNP level (f(x | m0))
fx0<-function(X_,psi_,position)
{
  px0<-1
  for(j in c(1,2))
  {
    px0=px0*exp(-(X_^2)/(2*psi_))/sqrt(2*pi*psi_)
  }
  return(px0)
}

# SNP level expectation of phi
ppost<-function(Phi,X_,psi_,m,position)
{
  
  numneg<-fx(X_,psi_,m[1],position)*Phi[1]
  numpos<-fx(X_,psi_,m[2],position)*Phi[2]
  den<-numneg + numpos + (fx0(X_,psi_,position)*(1-Phi[1] - Phi[2]))
  return(c(numneg,numpos)/den)
}

niter <- 100

Phipost<-matrix(rep(1/3, (niter+1)*2), nrow = niter+1 )
mpost<-matrix(c(rep(-0.03, (niter+1)),rep(0.03, (niter+1))), nrow = niter+1 )

for (iter in 1:niter)
{
  # Initialise parameters
  pmpost<-matrix(rep(1/3, nrow(X)*2), nrow = nrow(X) )
  m <-  matrix(c(rep(-0.03, nrow(X)),rep(0.03, nrow(X))), nrow = nrow(X) )
  for(position in seq(1,nrow(X)))
  {
    #Expectation at SNP level
    pmpost[position,]<-ppost(Phi = Phipost[iter,],X_ = X[position,1],psi_ = psi[position,1],m = mpost[iter,],position) + ppost(Phi = Phipost[iter,],X_ = X[position,2],psi_ = psi[position,2],m = mpost[iter,],position)
    m[position,]<- ppost(Phi = Phipost[iter,],X_ = X[position,1],psi_ = psi[position,1],m = mpost[iter,],position)*X[position,1] + ppost(Phi = Phipost[iter,],X_= X[position,2],psi_ = psi[position,2],m = mpost[iter,],position)*X[position,2]
  }
  # Maximise mixture distributions at overall level
  Phipost[iter+1,]<-(colSums(pmpost))/(nrow(X)*2)
  mpost[iter+1,] <- colSums(m)/colSums(pmpost)
}  