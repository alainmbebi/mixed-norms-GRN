# '@ Alain Mbebi
set.seed(143)
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# '@ libraries
library(data.table)
library(MASS)
library(corpcor)
library(Matrix)           
library(lattice)
library(glasso)           # This will help to estimate sparse Omega_Mixt using Graphical Lasso
library(glassoFast)       # Faster alternative to glasso
library(RColorBrewer)
library(matrixcalc)
library(igraph)
library(precrec)
library(reshape2)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the joint estimation for B and Omega in the mixed l1-l21 model
Mixedl2l21<-function(X, Y, lambdaO, lambdaB, eps=1e-5, zeta=1e-6, tol.out=1e-6){
  #lambdaO regularization parameter on the precition matrix
  #lambdaB regularization parameter on the regression coefficient matrix
  #X is the matrix of predictors (TFs)
  #Y is the matrix of responses (TGs)
  #zeta is the treshold value when computing the l21norm to avoid division by zero
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # t=0 initialization
  diag_C_0=c(rep(1, ncol(Y)))      
  invdiag_C_0=1/diag_C_0 
  P_0=diag(ncol(Y))
  Omega_0=(1/2*lambdaO)*(mpower(P_0^2 + 8*lambdaO*diag(ncol(Y)), 1/2))-(1/2*lambdaO)*P_0

  #----------------------
  SVD_X=svd(X)
  U1=SVD_X$u
  Gammma=diag(SVD_X$d) + diag(zeta, ncol(U1)) #nxn
  V1=SVD_X$v
  
  #----------------------
  Btilde0= matrix(0, nrow=ncol(V1), ncol=ncol(Y))
  K0=diag_C_0*P_0
  S0=t(V1)%*%(t(X))%*%Y
  for(j in 1:nrow(S0)){
    print(j)
    Btilde0[j,]=S0[j,]%*%solve((Gammma[j,j]^2)*diag(ncol(Y))+ (2*nrow(Y)*lambdaB)*K0)
  }
  B_0=V1%*%Btilde0
  
  #----------------------
  #compute the l21norm for t(B_0) that will help to update invdiag_C_1
  l21B_0=c(rep(0, ncol(Y)))
  for(i in 1:ncol(Y)){
    l21B_0[i]=2*norm((t(B_0)[i,]),type = "2")# + zeta
  }
  
  #----------------------
  P_1=round(cov.shrink(Y-X%*%B_0),9)
  diag_C_1=l21B_0
  invdiag_C_1=1/l21B_0
  Omega_1=(1/2*lambdaO)*(mpower(P_1^2 + 8*lambdaO*diag(ncol(Y)), 1/2))-(1/2*lambdaO)*P_1
  
  #----------------------
  #----------------------
  Btilde1= matrix(0, nrow=ncol(V1), ncol=ncol(Y))
  K1=diag_C_1*P_1
  S1=t(V1)%*%(t(X))%*%Y
  for(j in 1:nrow(S1)){
    print(j)
    Btilde1[j,]=S1[j,]%*%solve((Gammma[j,j]^2)*diag(ncol(Y))+ (2*nrow(Y)*lambdaB)*K1)
  }
  B_1=V1%*%Btilde1
  
  #----------------------
  l21B_1=c(rep(0, ncol(Y)))
  for(i in 1:ncol(Y)){
    l21B_1[i]=2*norm((t(B_1)[i,]),type = "2")# + zeta
  }
  
  #----------------------
  #t=2 update should follow in this way...
  
  tcontmixtl2l21=0
  if((sum(l21B_1)<=sum(l21B_0))==TRUE && tcontmixtl2l21<=itermax){
    l21B_0=l21B_1
    diag_C_0=diag_C_1
    invdiag_C_0=invdiag_C_1
    P_0=P_1
    Omega_0=Omega_1 
    
    #----------------------
    K0=K1
    S0=S1
    Btilde0= Btilde1
    B_0=B_1

    #----------------------
    P_1=round(cov.shrink(Y-X%*%B_0),9)
    diag_C_1=l21B_0
    invdiag_C_1=1/l21B_0
    Omega_1=(1/2*lambdaO)*(mpower(P_1^2 + 8*lambdaO*diag(ncol(Y)), 1/2))-(1/2*lambdaO)*P_1
    #----------------------
    Btilde1= matrix(0, nrow=ncol(V1), ncol=ncol(Y))
    K1=diag_C_1*P_1
    S1=t(V1)%*%(t(X))%*%Y
    for(j in 1:nrow(S1)){
      print(j)
      Btilde1[j,]=S1[j,]%*%solve((Gammma[j,j]^2)*diag(ncol(Y))+ (2*nrow(Y)*lambdaB)*K1)
    }
    B_1=V1%*%Btilde1
    
    
    #----------------------
    
    tcontmixtl2l21 <- sum(tcontmixtl2l21, 1)
    print(tcontmixtl2l21)
    
  }
  
  return(list(Bhat_l2l21=B_0, Omegahat_l2l21=Omega_0, Sigmahat_l2l21=P_0, iteropt=tcontmixtl2l21))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


