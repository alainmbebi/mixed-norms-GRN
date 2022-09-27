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
library(glassoFast)
library(matrixcalc)
library(igraph)
library(precrec)
library(reshape2)
library(minet)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# '@functions
#function for the joint estimation for B and Omega in the mixed l2l21 model
Mixedl2l21<-function(X, Y, lambdaO, lambdaB, eps=1e-5, zeta=1e-6, tol.out=1e-6){  
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
  l21B_0=c(rep(0, ncol(Y)))
  for(i in 1:ncol(Y)){
    l21B_0[i]=2*norm((t(B_0)[i,]),type = "2")# + zeta
  }
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
  l21B_1=c(rep(0, ncol(Y)))
  for(i in 1:ncol(Y)){
    l21B_1[i]=2*norm((t(B_1)[i,]),type = "2")# + zeta
  }
  
  #----------------------  
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
#Mixedl2l21 cross validation
CV_Mixedl2l21<-function( X, Y, lambdaO.vect, lambdaB.vect, kfold=3){
  set.seed(123)
  lambdaO.vect = sort(lambdaO.vect)
  lambdaB.vect = sort(lambdaB.vect)
  n = nrow(Y)
  CV_errors = matrix(0, nrow=length(lambdaO.vect), ncol=length(lambdaB.vect))
  ind=sample(n)
  for (k in 1:kfold){
    
    leave.out = ind[(1 + floor((k - 1) * n/kfold)):floor(k * n/kfold)]
    
    Y.train = Y[-leave.out,]
    X.train = X[-leave.out,]
    Y.valid = Y[leave.out, ]
    X.valid = X[leave.out, ]
         for(i in 1:length(lambdaO.vect)){
         for(j in 1:length(lambdaB.vect)){
          Estimes.train = Mixedl2l21(X.train, Y.train, lambdaO.vect[i], lambdaB.vect[j]) 
          B_hat.train = Estimes.train[[1]]
          Omega_hat_train =  Estimes.train[[2]]
          Sigma_hat_train =  Estimes.train[[3]]
          CV_errors[i,j] = mean(mse.matrix(X.valid%*%B_hat.train,Y.valid))#eventually use the RV coef
        }
       }
      }
  AVG = apply(CV_errors, 1, mean)
  error = min(AVG)
  opt = which.min(CV_errors) %% (dim(CV_errors)[1])
  opt = (opt != 0)*opt + (opt == 0)*(dim(CV_errors)[1])
  opt.i = opt
  opt.j = which.min(CV_errors[opt,])
  best.lamO = lambdaO.vect[opt.i]
  best.lamB = lambdaB.vect[opt.j]

  Estimes.all = Mixedl2l21(X, Y, best.lamO, best.lamB) 
  B_hat_final = Estimes.all[[1]]
  Omega_hat_final =  Estimes.all[[2]]
  Sigma_hat_final =  Estimes.all[[3]]
  return(list(best.lamO.final=best.lamO, best.lamB.final=best.lamB,
              B_hat_final=B_hat_final,Omega_hat_final=Omega_hat_final, Sigma_hat_final=Sigma_hat_final, min.error = error, avg.error = AVG, cv.err=CV_errors)) 
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
