# '@ Alain Mbebi
set.seed(143)
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# '@ libraries
library(data.table)
library(MASS)
library(corpcor)
library(Matrix)           # use Diagonal() when computing kronecker instead of diag()
library(lattice)
library(mvtnorm)
library(glasso)           # This will help to estimate sparse Omega_Mixt using Graphical Lasso
library(glassoFast)
library(RColorBrewer)
library(matrixcalc)
library(glmnet)
library(igraph)
library(precrec)
library(reshape2)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# '@functions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#function for the joint estimation for B and Omega in the mixed l1-l21 model
Mixedl2l21_Gauss<-function(X, Y, lambdaO, lambdaB, eps=1e-5, zeta=1e-6, tol.out=1e-6){
  #eps=1e-5, the algorithm will terminate if the minimum diagonal entry of the current iteration for
  #residual sample covariance is less than eps.
  #This may need adjustment depending on the scales of the variables.
  #zeta=1e-8 is the treshold value when computing the l21norm to avoid division by zero
  #----------------------
  #***********************************In fact Sigma_ here is P in the mixed l2l21_Gauss algorithm
  Sigma_0=diag(ncol(Y))
  Omega_0=(1/2*lambdaO)*(mpower(Sigma_0^2 + 8*lambdaO*diag(ncol(Y)), 1/2))-(1/2*lambdaO)*Sigma_0
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
      Btilde0[i,j]=S0[i,j]/((Gammma[i,i])^2 + 2*nrow(Y)*lambdaB*(Psi_0[j,j])^2)
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
  Sigma_1=round(cov.shrink(Y-X%*%B_0),9)
  Omega_1=(1/2*lambdaO)*(mpower(Sigma_1^2 + 8*lambdaO*diag(ncol(Y)), 1/2))-(1/2*lambdaO)*Sigma_1
  P1=t(chol(Sigma_1) + diag(zeta, ncol(Y)))               # to get the lower triangular matrix P (notice the transpose otherwise upper triangular)
  
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
      Btilde1[i,j]=S1[i,j]/((Gammma[i,i])^2 + 2*nrow(Y)*lambdaB*(Psi_1[j,j])^2)
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
  tcontmixtl2l21_Gauss=0

  if((sum(l21B_1)<=sum(l21B_0))==TRUE && tcontmixtl2l21_Gauss<=itermax){
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
    Sigma_1=round(cov.shrink(Y-X%*%B_0),9)
    Omega_1=(1/2*lambdaO)*(mpower(Sigma_1^2 + 8*lambdaO*diag(ncol(Y)), 1/2))-(1/2*lambdaO)*Sigma_1
    P1=t(chol(Sigma_1) + diag(zeta, ncol(Y)))               # to get the lower triangular matrix P (notice the transpose otherwise upper triangular)
    
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
        Btilde1[i,j]=S1[i,j]/((Gammma[i,i])^2 + 2*nrow(Y)*lambdaB*(Psi_1[j,j])^2)
      }
    }
    B_1=V1%*%Btilde1%*%t(U2_1)
    
    #----------------------
    
    tcontmixtl2l21_Gauss <- sum(tcontmixtl2l21_Gauss, 1)
    print(tcontmixtl2l21_Gauss)
    
  }
  
  return(list(Bhat_l2l21_Gauss=B_0, Omegahat_l2l21_Gauss=Omega_0, Sigmahat_l2l21_Gauss=Sigma_0, iteropt=tcontmixtl2l21_Gauss))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Mixedl2l21_Gauss cross validation
CV_Mixedl2l21_Gauss<-function( X, Y, lambdaO.vect, lambdaB.vect, kfold=3){
  set.seed(123)
  lambdaO.vect = sort(lambdaO.vect)
  lambdaB.vect = sort(lambdaB.vect)
  n = nrow(Y)
  CV_errors = matrix(0, nrow=length(lambdaO.vect), ncol=length(lambdaB.vect))
  ind=sample(n)
  for (k in 1:kfold){
    
    leave.out = ind[(1 + floor((k - 1) * n/kfold)):floor(k * n/kfold)]
    
    # training set
    Y.train = Y[-leave.out,]
    X.train = X[-leave.out,]
    Y.valid = Y[leave.out, ]
    X.valid = X[leave.out, ]

    # loop over all tuning parameters
         for(i in 1:length(lambdaO.vect)){
         for(j in 1:length(lambdaB.vect)){
          Estimes.train = Mixedl2l21_Gauss(X.train, Y.train, lambdaO.vect[i], lambdaB.vect[j]) 
          B_hat.train = Estimes.train[[1]]
          Omega_hat_train =  Estimes.train[[2]]
          Sigma_hat_train =  Estimes.train[[3]]
          CV_errors[i,j] = mean(mse.matrix(X.valid%*%B_hat.train,Y.valid))#eventually use the RV coef
        }
       }
      }
    # determine optimal tuning parameters
  AVG = apply(CV_errors, 1, mean)
  error = min(AVG)
  opt = which.min(CV_errors) %% (dim(CV_errors)[1])
  opt = (opt != 0)*opt + (opt == 0)*(dim(CV_errors)[1])
  opt.i = opt
  opt.j = which.min(CV_errors[opt,])
  best.lamO = lambdaO.vect[opt.i]
  best.lamB = lambdaB.vect[opt.j]

  #output on all data
  Estimes.all = Mixedl2l21_Gauss(X, Y, best.lamO, best.lamB) 
  B_hat_final = Estimes.all[[1]]
  Omega_hat_final =  Estimes.all[[2]]
  Sigma_hat_final =  Estimes.all[[3]]
  # return best values in lambdaO.vect, lambdaB.vect and other meaningful values
  return(list(best.lamO.final=best.lamO, best.lamB.final=best.lamB,
              B_hat_final=B_hat_final,Omega_hat_final=Omega_hat_final, Sigma_hat_final=Sigma_hat_final, min.error = error, avg.error = AVG, cv.err=CV_errors)) 
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

