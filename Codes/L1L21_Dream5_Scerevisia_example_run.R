#This is an example run for L1L21 solution with S. cerevisiae data (network 4 from DREAM5 challenge)
source("Mixed_L1L21_GRN.R")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### read expression data
expr <- as.matrix(read.table(file = "net4_expression_data.tsv",header = T,sep = "\t",stringsAsFactors = F))
expr=scale(expr)
##### read transcription factors
TFs <- read.table(file = "net4_transcription_factors.tsv",header = F,sep = "\t",stringsAsFactors = F)
TFs <- TFs$V1
##### read Gold standard network
GS <- read.table(file = "DREAM5_NetworkInference_GoldStandard_Network4 - S. cerevisiae.tsv",header = F,sep = "\t",stringsAsFactors = F)
colnames(GS) <- c( "from","to",paste("weight","GS",sep = "_"))
ordernet <- function(net,TFs)
{
  print(dim(net))
  
  net$from <- toupper(trimws(net$from))
  net$to <- toupper(trimws(net$to))
  
  wrong_orders <- which(!net$from%in%TFs)
  tmp <- net[wrong_orders,1]
  net[wrong_orders,1] <- net[wrong_orders,2]
  net[wrong_orders,2] <- tmp
  
  net <- net[order(net$from,net$to),]
  net <- net[which(net$from!=net$to),]
  net <- net[!duplicated(net[,1:2]),]
  return(net)
}
GS <- ordernet(GS,TFs)
GS_genes <- unique(c(GS$from,GS$to))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGsexpr <- expr[,which(colnames(expr)%in%GS_genes)]
TFsexpr <- TGsexpr[,colnames(TGsexpr)%in%TFs]
TFs=colnames(TFsexpr)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

itermax=50 
Y=TGsexpr #the response matrix 
X=TFsexpr #the predictor matrix 
#lambdaO=2^(-6) and lambdaB=0.8 are the optimal hyperparameters learned from 10-fold CV. This will also reproduce results (for Scerevisiae) in the manuscript.
ttt=Mixedl1l21(X, Y, lambdaO=2^(-6), lambdaB=0.8) 

#Get the estimated regression coefficient matrix
Bhat_l1l21_D5_Scerevisiae=ttt$Bhat_l1l21
rownames(Bhat_l1l21_D5_Scerevisiae) <- colnames(X)
colnames(Bhat_l1l21_D5_Scerevisiae) <- colnames(Y)

#Get the estimated precision matrix
Omegahat_l1l21_D5_Scerevisiae=ttt$Omegahat_l1l21
rownames(Omegahat_l1l21_D5_Scerevisiae) <- colnames(Y)
colnames(Omegahat_l1l21_D5_Scerevisiae) <- colnames(Y)

#Get the estimated covariance matrix
Sigmahat_l1l21_D5_Scerevisiae=ttt$Sigmahat_l1l21
rownames(Sigmahat_l1l21_D5_Scerevisiae) <- colnames(Y)
colnames(Sigmahat_l1l21_D5_Scerevisiae) <- colnames(Y)

