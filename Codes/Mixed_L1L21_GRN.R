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
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#function for the joint estimation for B and Omega in the mixed l1-l21 model
Mixedl1l21<-function(X, Y, lambdaO, lambdaB, eps=1e-5, zeta=1e-6, tol.out=1e-6){
  #eps=1e-5, the algorithm will terminate if the minimum diagonal entry of the current iteration for
  #residual sample covariance is less than eps.
  #This may need adjustment depending on the scales of the variables.
  #zeta=1e-8 is the treshold value when computing the l21norm to avoid division by zero
  #----------------------
  
  # t=0 initialization
  diag_C_0=c(rep(1, ncol(Y)))       # Initialize C as the identity matrix of dim the number nrow of W_true  
  invdiag_C_0=1/diag_C_0 
  Sigma_0=diag(ncol(Y))
  #Sigma_0=cov(Y) + diag(zeta, ncol(Y))            # Initialize Omega as the sample covariance matrix of dim the number ncol of Y
  #Sigma_0=cov.shrink(Y)  # estimate the cov as the shrinkage covariance estimator using cov.shrink(X)   # shrinkage estimate 
  Omega_0=glassoFast(Sigma_0, rho=lambdaO,thr=1.0e-4, maxIt=1e4)$wi
  #Omega_0= mpower(Sigma_0, -1)
  
  #----------------------
  SVD_X=svd(X)
  U1=SVD_X$u
  Gammma=diag(SVD_X$d) + diag(zeta, ncol(U1)) #nxn
  V1=SVD_X$v
  
  #----------------------
  Btilde0= matrix(0, nrow=ncol(V1), ncol=ncol(Y))
  K0=diag_C_0*Sigma_0
  S0=t(V1)%*%(t(X))%*%Y
  #library(Rlinsolve)
  for(j in 1:nrow(S0)){
    print(j)
    Btilde0[j,]=S0[j,]%*%solve((Gammma[j,j]^2)*diag(ncol(Y)) + 2*nrow(Y)*lambdaB*K0)
  }
  B_0=V1%*%Btilde0
  
  #----------------------
  #compute the l21norm for t(B_0) that will help to update invdiag_C_1
  l21B_0=c(rep(0, ncol(Y)))
  for(i in 1:ncol(Y)){
    l21B_0[i]=2*norm((t(B_0)[i,]),type = "2")# + zeta
  }
  
  #----------------------
  
  #t=1 update
  Sigma_1=cov.shrink(Y-X%*%B_0)
  diag_C_1=l21B_0
  invdiag_C_1=1/l21B_0
  Omega_1=glassoFast(Sigma_1, rho=lambdaO,thr=1.0e-4, maxIt=1e4)$wi    #get the precision matrix at t=1 from glasso
  
  #----------------------
  Btilde1= matrix(0, nrow=ncol(V1), ncol=ncol(Y))
  K1=diag_C_1*Sigma_1
  S1=t(V1)%*%(t(X))%*%Y
  #library(Rlinsolve)
  for(j in 1:nrow(S1)){
    print(j)
    Btilde1[j,]=S1[j,]%*%solve((Gammma[j,j]^2)*diag(ncol(Y))+ 2*nrow(Y)*lambdaB*K1)
  }
  B_1=V1%*%Btilde1
  
  #----------------------
  #compute the l21norm for t(B_0) that will help to update invdiag_C_1
  l21B_1=c(rep(0, ncol(Y)))
  for(i in 1:ncol(Y)){
    l21B_1[i]=2*norm((t(B_1)[i,]),type = "2")# + zeta
  }
  
  #----------------------
  #t=2 update should follow in this way...
  
  tcontmixtl1l21=0
  if((sum(l21B_1)<=sum(l21B_0))==TRUE && tcontmixtl1l21<=itermax){
    l21B_0=l21B_1
    diag_C_0=diag_C_1
    invdiag_C_0=invdiag_C_1
    Sigma_0=Sigma_1
    Omega_0=Omega_1 
    
    #----------------------
    K0=K1
    S0=S1
    Btilde0= Btilde1
    B_0=B_1

    #----------------------
    Sigma_1=cov.shrink(Y-X%*%B_0)
    diag_C_1=l21B_0
    invdiag_C_1=1/l21B_0
    Omega_1=glassoFast(Sigma_1, rho=lambdaO,thr=1.0e-4, maxIt=1e4)$wi    #get the precision matrix at t=1 from glasso
    #----------------------
    Btilde1= matrix(0, nrow=ncol(V1), ncol=ncol(Y))
    K1=diag_C_1*Sigma_1
    S1=t(V1)%*%(t(X))%*%Y
    for(j in 1:nrow(S1)){
      print(j)
      Btilde1[j,]=S1[j,]%*%solve((Gammma[j,j]^2)*diag(ncol(Y))+ 2*nrow(Y)*lambdaB*K1)
    }
    B_1=V1%*%Btilde1    
    
    #----------------------
    
    tcontmixtl1l21 <- sum(tcontmixtl1l21, 1)
    print(tcontmixtl1l21)
    
  }
  
  return(list(Bhat_l1l21=B_0, Omegahat_l1l21=Omega_0, Sigmahat_l1l21=Sigma_0, iteropt=tcontmixtl1l21))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Mixedl1l21 cross validation
CV_Mixedl1l21<-function( X, Y, lambdaO.vect, lambdaB.vect, kfold=3){
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
          Estimes.train = Mixedl1l21(X.train, Y.train, lambdaO.vect[i], lambdaB.vect[j]) 
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
  Estimes.all = Mixedl1l21(X, Y, best.lamO, best.lamB) 
  B_hat_final = Estimes.all[[1]]
  Omega_hat_final =  Estimes.all[[2]]
  Sigma_hat_final =  Estimes.all[[3]]
  # return best values in lambdaO.vect, lambdaB.vect and other meaningful values
  return(list(best.lamO.final=best.lamO, best.lamB.final=best.lamB,
              B_hat_final=B_hat_final,Omega_hat_final=Omega_hat_final, Sigma_hat_final=Sigma_hat_final, min.error = error, avg.error = AVG, cv.err=CV_errors)) 
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
