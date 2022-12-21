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
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# '@functions

#function for the joint estimation for B and Omega in the mixed l1-l21 model
Mixedl1l21_Gauss<-function(X, Y, lambdaO, lambdaB, zeta=1e-8){
  #lambdaO regularization parameter on the precition matrix
  #lambdaB regularization parameter on the regression coefficient matrix
  #X is the matrix of predictors (TFs)
  #Y is the matrix of responses (TGs)
  #zeta is the treshold value when computing the l21norm to avoid division by zero
  #----------------------

  Sigma_0=diag(ncol(Y))
  Omega_0=glassoFast(Sigma_0, rho=lambdaO,thr=1.0e-4, maxIt=1e4)$wi
  P0=t(chol(Sigma_0) + diag(zeta, ncol(Y)))               # to get the lower triangular matrix P (notice the transpose otherwise upper triangular)
  
  #----------------------
  SVD_P0=svd(P0)
  U2_0=SVD_P0$u
  Psi_0=diag(SVD_P0$d) + diag(zeta, ncol(Y)) #sxs
  V2_0=SVD_P0$v
  
  #----------------------
  SVD_X=svd(X)
  U1=SVD_X$u
  Gammma=diag(SVD_X$d) + diag(zeta, ncol(U1)) #nxn
  V1=SVD_X$v
  
  #----------------------
  Btilde0= matrix(0, nrow=ncol(V1), ncol=ncol(Y))
  S0=t(V1)%*%(t(X))%*%Y%*%U2_0
  for(i in 1:nrow(S0)){
    for(j in 1:ncol(S0)){
      Btilde0[i,j]=S0[i,j]/((Gammma[i,i])^2 + nrow(Y)*lambdaB*(Psi_0[j,j])^2)
    }
  }
  B_0=V1%*%Btilde0%*%t(U2_0)
  
  #----------------------
  #compute the l21norm for t(B_0) that will help to update invdiag_C_1
  l21B_0=c(rep(0, ncol(Y)))
  for(i in 1:ncol(Y)){
    l21B_0[i]=2*norm((t(B_0)[i,]),type = "2")# + zeta
  }
  
  #----------------------
  #t=1 update
  Sigma_1=cov.shrink(Y-X%*%B_0)  # estimate the cov as the shrinkage covariance estimator using cov.shrink(X)   # shrinkage estimate
  Omega_1=glassoFast(Sigma_1, rho=lambdaO,thr=1.0e-4, maxIt=1e4)$wi    #get the precision matrix at t=1 from glasso
  P1=t(chol(Sigma_1) + diag(zeta, ncol(Y)))              
  
  #----------------------
  SVD_P1=svd(P1)
  U2_1=SVD_P1$u
  Psi_1=diag(SVD_P1$d) + diag(zeta, ncol(Y)) #sxs
  V2_1=SVD_P1$v

  #----------------------
  Btilde1= matrix(0, nrow=ncol(V1), ncol=ncol(Y))
  S1=t(V1)%*%(t(X))%*%Y%*%U2_1
  for(i in 1:nrow(S1)){
    for(j in 1:ncol(S1)){
      Btilde1[i,j]=S1[i,j]/((Gammma[i,i])^2 + nrow(Y)*lambdaB*(Psi_1[j,j])^2)
    }
  }
  B_1=V1%*%Btilde1%*%t(U2_1)
  
  #----------------------
  #compute the l21norm for t(B_0) that will help to update invdiag_C_1
  l21B_1=c(rep(0, ncol(Y)))
  for(i in 1:ncol(Y)){
    l21B_1[i]=2*norm((t(B_1)[i,]),type = "2")
  }

  #----------------------
  
  tcontmixtl1l21_Gauss=0

  if((sum(l21B_1)<=sum(l21B_0))==TRUE && tcontmixtl1l21_Gauss<=itermax){
    l21B_0=l21B_1
    Sigma_0=Sigma_1
    Omega_0=Omega_1 
    P0=P1
    
    #----------------------
    SVD_P0=SVD_P1
    U2_0=U2_1
    Psi_0=Psi_1
    V2_0=V2_1
    
    #----------------------
    S0=S1
    Btilde0= Btilde1
    B_0=B_1

    #----------------------
    #Now we replace all values above
    Sigma_1=cov.shrink(Y-X%*%B_0)  
    Omega_1=glassoFast(Sigma_1, rho=lambdaO,thr=1.0e-4, maxIt=1e4)$wi    
    P1=t(chol(Sigma_1) + diag(zeta, ncol(Y)))              
    
    #----------------------
    SVD_P1=svd(P1)
    U2_1=SVD_P1$u
    Psi_1=diag(SVD_P1$d) + diag(zeta, ncol(Y)) 
    V2_1=SVD_P1$v
    
    #----------------------
    Btilde1= matrix(0, nrow=ncol(V1), ncol=ncol(Y))
    S1=t(V1)%*%(t(X))%*%Y%*%U2_1
    for(i in 1:nrow(S1)){
      for(j in 1:ncol(S1)){
        Btilde1[i,j]=S1[i,j]/((Gammma[i,i])^2 + nrow(Y)*lambdaB*(Psi_1[j,j])^2)
      }
    }
    B_1=V1%*%Btilde1%*%t(U2_1)
    
    #----------------------
    
    tcontmixtl1l21_Gauss <- sum(tcontmixtl1l21_Gauss, 1)
    print(tcontmixtl1l21_Gauss)
    
  }
  
  return(list(Bhat_l1l21_Gauss=B_0, Omegahat_l1l21_Gauss=Omega_0, Sigmahat_l1l21_Gauss=Sigma_0, iteropt=tcontmixtl1l21_Gauss))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

