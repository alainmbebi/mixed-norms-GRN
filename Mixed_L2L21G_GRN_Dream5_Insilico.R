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
library(CVTuningCov)
library(glmnet)
library(hydroGOF)
library(igraph)
#library(precrec)
library(reshape2)
#library(minet)
#library(expm)             # to find the unique psd sqrt of matrices with the function sqrtm().
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# '@functions
#function to center data matrices using 'colMeans()' note that we center the col because they represent variables
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function to compute the L2-norm of a given vector x
L2norm <- function(x) {
  sqrt(crossprod(x))
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function to compute the ridge estimate
ridgeEst<-function(lambdaB, X, Y){
  Bhat_ridge= round(t(X)%*%solve(X%*%t(X) + lambdaB*diag(nrow(Y)))%*%Y ,9)
  return(Bhat_ridge)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
  #Sigma_0=cov(Y) + diag(zeta, ncol(Y))            # Initialize Omega as the sample covariance matrix of dim the number ncol of Y 
  #Sigma_0=round(cov.shrink(Y),9)
  #Omega_0=(1/2*lambdaO)*(sqrtm(Sigma_0^2 + 8*lambdaO*diag(ncol(Y))))-(1/2*lambdaO)*Sigma_0
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
 #Sigma_1=cov(Y-X%*%B_0) + diag(zeta, ncol(Y))
  Sigma_1=round(cov.shrink(Y-X%*%B_0),9)
  #Omega_1=(1/2*lambdaO)*(sqrtm(Sigma_1^2 + 8*lambdaO*diag(ncol(Y))))-(1/2*lambdaO)*Sigma_1
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
  #t=2 update should follow in this way...
  
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
    #Now we replace all values above
    #Sigma_1=cov(Y-X%*%B_0) + diag(zeta, ncol(Y))
    Sigma_1=round(cov.shrink(Y-X%*%B_0),9)
    #Omega_1=(1/2*lambdaO)*(sqrtm(Sigma_1^2 + 8*lambdaO*diag(ncol(Y))))-(1/2*lambdaO)*Sigma_1
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
          # compute the joint penalized regression matrix estimate with the iner training set
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
##remove nonTF-nonTF interactions
rm_NTF_int <- function(net, tfs)
{
  tmp <- net
  net$from <- net$from%in%tfs
  net$to <- net$to%in%tfs
  return(tmp[which(net$from | net$to),])
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### read expression data
#expr <- as.matrix(read.table(file = "net1_expression_data.tsv",header = T,sep = "\t",stringsAsFactors = F))
expr<-t(as.matrix(read.csv("net1_expression_data_log.csv", head=FALSE, row.names = 1, sep=",")) )
head(expr[,1:6])
dim(expr)
expr=scale(expr)
##### read transcription factors
TFs <- read.table(file = "net1_transcription_factors.tsv",header = F,sep = "\t",stringsAsFactors = F)
TFs <- TFs$V1
head(TFs)
length(TFs)
##### read Gold standard network
GS <- read.table(file = "DREAM5_NetworkInference_GoldStandard_Network1_Insilico.tsv",header = F,sep = "\t",stringsAsFactors = F)
GS <- GS[,c(2,1,3)]
colnames(GS) <- c( "from","to",paste("weight","GS",sep = "_"))
head(GS)
dim(GS)
ordernet <- function(net,TFs)
{
  print(dim(net))
  
  net$from <- toupper(trimws(net$from))
  net$to <- toupper(trimws(net$to))
  
  wrong_orders <- which(!net$from%in%TFs)
  tmp <- net[wrong_orders,1]
  net[wrong_orders,1] <- net[wrong_orders,2]
  net[wrong_orders,2] <- tmp
  
  ############## (5) Merge duplicates and remove duplicates and self interactions
  # sort based on from and to columns
  net <- net[order(net$from,net$to),]
  
  # remove self interactions
  net <- net[which(net$from!=net$to),]
  
  print(dim(net))
  
  net <- net[!duplicated(net[,1:2]),]
  
  print(dim(net))
  
  return(net)
}
GS <- ordernet(GS,TFs)
head(GS)
dim(GS)

GS_genes <- unique(c(GS$from,GS$to))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### select expression profiles for all genes
#TGsexpr <- expr
#TFsexpr <- expr[,colnames(expr)%in%TFs]
### select expression profiles for only GS_genes
TGsexpr <- expr[,which(colnames(expr)%in%GS_genes)]
TFsexpr <- TGsexpr[,colnames(TGsexpr)%in%TFs]
TFs=colnames(TFsexpr)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

itermax=5
#Yl21 = scale(responses, scale = TRUE) #dim n x s
#Xl21 = scale(predictors, scale = TRUE) #dim n x p
Yl21=TGsexpr #expr
Xl21=TFsexpr
### remove duplicate colnames in a matrix Yl21 and Xl21 with
Yl21 <- Yl21[, !duplicated(colnames(Yl21))]
Xl21 <- Xl21[, !duplicated(colnames(Xl21))]

ttt=Mixedl2l21_Gauss(Xl21, Yl21, lambdaO=2^(-6), lambdaB=.5) 
#in CV set: lambdaO=(2^(-8),2^(-7)) and lambdaB=seq(0.2, 0.5,0.1)
#when all genes are used 
#best tuning for old model, i.e. C is identity
#lambdaO=2^(-4), lambdaB=2^(-1)  AUPR=0.014  and AUROC=0.634 for all expr genes and bad  precision matrix with scale(expr)
#lambdaO=2^(-5), lambdaB=2^(1)   AUPR=0.012  and AUROC=0.681 for all expr genes and bad  precision matrix with scale(expr)
#lambdaO=2^(-5), lambdaB=2^(2)   AUPR=0.010  and AUROC=0.683 for all expr genes and bad  precision matrix with scale(expr)
#with shrinkage covariance estimator 
#lambdaO=2^(-4), lambdaB=2^(2)   AUPR=0.087  and AUROC=0.639 for all expr genes and good  precision matrix and shrinked initial cov mat
#lambdaO=2^(-6), lambdaB=2^(4)   AUPR=0.061  and AUROC=0.672 for all expr genes and good  precision matrix and shrinked initial cov mat

#Bhat_l2l21_Gauss_D5_Insilico<-t(as.matrix(read.csv("Bhat_l2l21_Gauss_D5_Insilico.csv", head=TRUE, row.names = 1, sep=","))) 
#Omegahat_l2l21_Gauss_D5_Insilico<-as.matrix(read.csv("Omegahat_l2l21_Gauss_D5_Insilico.csv", head=TRUE, row.names = 1, sep=",")) 
#Sigmahat_l2l21_Gauss_D5_Insilico<-as.matrix(read.csv("Sigmahat_l2l21_Gauss_D5_Insilico.csv", head=TRUE, row.names = 1, sep=",")) 
#Bhat_l2l21_Gauss_D5_Insilico[abs(Bhat_l2l21_Gauss_D5_Insilico)<0.00009]=0
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bhat_l2l21_Gauss_D5_Insilico=ttt$Bhat_l2l21_Gauss
rownames(Bhat_l2l21_Gauss_D5_Insilico) <- colnames(Xl21)
colnames(Bhat_l2l21_Gauss_D5_Insilico) <- colnames(Yl21)

Omegahat_l2l21_Gauss_D5_Insilico=ttt$Omegahat_l2l21_Gauss
rownames(Omegahat_l2l21_Gauss_D5_Insilico) <- colnames(Yl21)
colnames(Omegahat_l2l21_Gauss_D5_Insilico) <- colnames(Yl21)

Sigmahat_l2l21_Gauss_D5_Insilico=ttt$Sigmahat_l2l21_Gauss
rownames(Sigmahat_l2l21_Gauss_D5_Insilico) <- colnames(Yl21)
colnames(Sigmahat_l2l21_Gauss_D5_Insilico) <- colnames(Yl21)

write.csv(t(Bhat_l2l21_Gauss_D5_Insilico), "Bhat_l2l21_Gauss_D5_Insilico.csv")
write.csv(Omegahat_l2l21_Gauss_D5_Insilico, "Omegahat_l2l21_Gauss_D5_Insilico.csv")
write.csv(Sigmahat_l2l21_Gauss_D5_Insilico, "Sigmahat_l2l21_Gauss_D5_Insilico.csv")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Rescale each column to range between 0 and 1
#Bhat_l2l21_Gauss_D5_Insilico=t(apply(abs(Bhat_l2l21_Gauss_D5_Insilico), 1, function(x)(x-min(x))/(max(x)-min(x))))
#Bhat_l2l21_Gauss_D5_Insilico=t(apply(Bhat_l2l21_Gauss_D5_Insilico, 1, function(x)(abs(x))/abs(max(x))))
mixedl2l21_Gauss_D5_Insilico <- reshape2::melt(Bhat_l2l21_Gauss_D5_Insilico)
mixedl2l21_Gauss_D5_Insilico <- mixedl2l21_Gauss_D5_Insilico[,c(2,1,3)]
colnames(mixedl2l21_Gauss_D5_Insilico) <- c("from","to","weight_mixedl2l21_Gauss_D5_Insilico")
mixedl2l21_Gauss_D5_Insilico=rm_NTF_int(mixedl2l21_Gauss_D5_Insilico,TFs)
mixedl2l21_Gauss_D5_Insilico <- ordernet(mixedl2l21_Gauss_D5_Insilico,TFs)

mixedl2l21_Gauss_D5_Insilico <- mixedl2l21_Gauss_D5_Insilico[which(abs(mixedl2l21_Gauss_D5_Insilico$weight_mixedl2l21_Gauss_D5_Insilico)>0.00001),]
mixedl2l21_Gauss_D5_Insilico$score <- abs(mixedl2l21_Gauss_D5_Insilico$weight_mixedl2l21_Gauss_D5_Insilico)/max(abs(mixedl2l21_Gauss_D5_Insilico$weight_mixedl2l21_Gauss_D5_Insilico))
mixedl2l21_Gauss_D5_Insilico <- mixedl2l21_Gauss_D5_Insilico[order(mixedl2l21_Gauss_D5_Insilico$score,decreasing = T),]
mixedl2l21_Gauss_D5_Insilico=head(mixedl2l21_Gauss_D5_Insilico,100000)
colnames(mixedl2l21_Gauss_D5_Insilico)[3:4] <-  c("weight_mixedl2l21_Gauss_D5_Insilico","score_mixedl2l21_Gauss_D5_Insilico")

nl2l21_Gauss_D5_Insilico <- merge(mixedl2l21_Gauss_D5_Insilico,GS,by = c("from","to"),all.y = T)
nl2l21_Gauss_D5_Insilico$weight_mixedl2l21_Gauss_D5_Insilico[which(is.na(nl2l21_Gauss_D5_Insilico$weight_mixedl2l21_Gauss_D5_Insilico)==T)] <- 0
#GS=head(GS,100000)
library(precrec)
sscurves_nl2l21_Gauss_D5_Insilico <- evalmod(scores = nl2l21_Gauss_D5_Insilico$weight_mixedl2l21_Gauss_D5_Insilico, labels = GS$weight_GS)
plot(sscurves_nl2l21_Gauss_D5_Insilico)
auc(sscurves_nl2l21_Gauss_D5_Insilico)

nl2l21_Gauss_D5_Insilico_net <- as.matrix(get.adjacency(graph.data.frame(mixedl2l21_Gauss_D5_Insilico,directed=T),attr="weight_mixedl2l21_Gauss_D5_Insilico"))
nl2l21_Gauss_D5_Insilico_net <- cbind(nl2l21_Gauss_D5_Insilico_net,matrix(0,nrow(nl2l21_Gauss_D5_Insilico_net),length(setdiff(GS_genes,colnames(nl2l21_Gauss_D5_Insilico_net))),dimnames=list(rownames(nl2l21_Gauss_D5_Insilico_net), setdiff(GS_genes,colnames(nl2l21_Gauss_D5_Insilico_net)))))
nl2l21_Gauss_D5_Insilico_net <- rbind(nl2l21_Gauss_D5_Insilico_net,matrix(0,length(setdiff(GS_genes,rownames(nl2l21_Gauss_D5_Insilico_net))),ncol(nl2l21_Gauss_D5_Insilico_net),dimnames=list(setdiff(GS_genes,rownames(nl2l21_Gauss_D5_Insilico_net)),colnames(nl2l21_Gauss_D5_Insilico_net))))

GS_net <- as.matrix(get.adjacency(graph.data.frame(GS,directed=T),attr="weight_GS"))
#GS_net<-GS_net[rownames(nl2l21_Gauss_D5_Insilico_net),colnames(nl2l21_Gauss_D5_Insilico_net)] #because the generated network is bigger than the GS, we must change the subsetting
nl2l21_Gauss_D5_Insilico_net<-nl2l21_Gauss_D5_Insilico_net[rownames(GS_net),colnames(GS_net)]
#### convert the predicted network into a binary matrix 
#newnl2l21_Gauss_D5_Insilico_net=as.matrix((nl2l21_Gauss_D5_Insilico_net>0)+0)
library(minet)
table_l2l21_Gauss_D5_Insilico <- validate(inet=nl2l21_Gauss_D5_Insilico_net,tnet=GS_net)
auc.roc(table_l2l21_Gauss_D5_Insilico)
auc.pr(table_l2l21_Gauss_D5_Insilico)

dev <- show.roc(table_l2l21_Gauss_D5_Insilico,col="orange", type="l",lty=2 ,lwd=4)
dev <- show.pr(table_l2l21_Gauss_D5_Insilico,col="orange", type="l",lty=2 ,lwd=4)
#count the number of non zero entries
nnzero(Omegahat_l2l21_Gauss_D5_Insilico, na.counted = NA)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Needed for genes selection with cond indep
final_Omegahat_l2l21_Gauss_D5_Insilico=Omegahat_l2l21_Gauss_D5_Insilico[c(rownames(Bhat_l2l21_Gauss_D5_Insilico)),]
Sparse_Bhat_l2l21_Gauss_D5_Insilico=Bhat_l2l21_Gauss_D5_Insilico
for (i in 1:nrow(Sparse_Bhat_l2l21_Gauss_D5_Insilico)) {
  for (j in 1:ncol(Sparse_Bhat_l2l21_Gauss_D5_Insilico)) {
    if(!(Sparse_Bhat_l2l21_Gauss_D5_Insilico[i,j] == 0) && (final_Omegahat_l2l21_Gauss_D5_Insilico[i,j])==0){
      Sparse_Bhat_l2l21_Gauss_D5_Insilico[i,j]<-0
    }
  }
}

write.csv(t(Sparse_Bhat_l2l21_Gauss_D5_Insilico), "Sparse_Bhat_l2l21_Gauss_D5_Insilico.csv")
write.csv(t(final_Omegahat_l2l21_Gauss_D5_Insilico), "final_Omegahat_l2l21_Gauss_D5_Insilico.csv")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##plot with package pROC
library(pROC)
pROC_obj <- roc(GS$weight_GS, nl2l21_Gauss_D5_Insilico$weight_mixedl2l21_Gauss_D5_Insilico,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")

## Warning in plot.ci.se(sens.ci, type = "shape", col = "lightblue"): Low
## definition shape.

plot(sens.ci, type="bars")
### remove duplicate colnames in a matrix Yl21 with
Yl21 <- Yl21[, !duplicated(colnames(Yl21))]

gl2l21_Gauss_D5_Insilico_net <- as.matrix(get.adjacency(graph.data.frame(genie_all,directed=T),attr="weight_mixedl2l21_Gauss_D5_Insilico"))
gl2l21_Gauss_D5_Insilico_net <- cbind(gl2l21_Gauss_D5_Insilico_net,matrix(0,nrow(gl2l21_Gauss_D5_Insilico_net),length(setdiff(GS_genes,colnames(gl2l21_Gauss_D5_Insilico_net))),dimnames=list(rownames(gl2l21_Gauss_D5_Insilico_net), setdiff(GS_genes,colnames(gl2l21_Gauss_D5_Insilico_net)))))
gl2l21_Gauss_D5_Insilico_net <- rbind(gl2l21_Gauss_D5_Insilico_net,matrix(0,length(setdiff(GS_genes,rownames(gl2l21_Gauss_D5_Insilico_net))),ncol(gl2l21_Gauss_D5_Insilico_net),dimnames=list(setdiff(GS_genes,rownames(gl2l21_Gauss_D5_Insilico_net)),colnames(gl2l21_Gauss_D5_Insilico_net))))


#test code to set self interadtion to zero

Bhat_l2l21_Gauss_D5_Insilico_minself=Bhat_l2l21_Gauss_D5_Insilico
for (i in 1:nrow(Bhat_l2l21_Gauss_D5_Insilico_minself)) {
  for (j in 1:ncol(Bhat_l2l21_Gauss_D5_Insilico_minself)) {
    Bhat_l2l21_Gauss_D5_Insilico_minself[i,i]<-0
  }
}

##B_0=solve(X%*%Omega_0%*%t(X) + nrow(Y)*lambdaB*diag(nrow(X)))%*%X%*%Omega_0%*%t(Y)
#the new l21 norm willl be

