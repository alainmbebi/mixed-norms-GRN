#This is an example run for L1L21 solution with S. cerevisiae data (network 4 from DREAM5 challenge)
source("Mixed_L1L21_GRN.R")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### read expression data
expr <- as.matrix(read.table(file = "net4_expression_data.tsv",header = T,sep = "\t",stringsAsFactors = F))
head(expr[,1:6])
dim(expr)
expr=scale(expr)
##### read transcription factors
TFs <- read.table(file = "net4_transcription_factors.tsv",header = F,sep = "\t",stringsAsFactors = F)
TFs <- TFs$V1
head(TFs)
length(TFs)
##### read Gold standard network
GS <- read.table(file = "DREAM5_NetworkInference_GoldStandard_Network4 - S. cerevisiae.tsv",header = F,sep = "\t",stringsAsFactors = F)
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
  
  net <- net[order(net$from,net$to),]
  net <- net[which(net$from!=net$to),]
  net <- net[!duplicated(net[,1:2]),]
  return(net)
}
GS <- ordernet(GS,TFs)
head(GS)
dim(GS)

GS_genes <- unique(c(GS$from,GS$to))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### If the goal is to learn new interactions, then select expression profiles for all genes
#TGsexpr <- expr
#TFsexpr <- expr[,colnames(expr)%in%TFs]

### Otherwise (e.g. comparative analysis) select expression profiles for only GS genes only
TGsexpr <- expr[,which(colnames(expr)%in%GS_genes)]
TFsexpr <- TGsexpr[,colnames(TGsexpr)%in%TFs]
TFs=colnames(TFsexpr)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

itermax=50 
Y=TGsexpr #the response matrix 
X=TFsexpr #the predictor matrix 
#lambdaO=2^(-6) and lambdaB=0.8 are the optimal hyperparameters learned from 10-fold CV. This will also reproduce results in the manuscript.
ttt=Mixedl1l21(X, Y, lambdaO=2^(-6), lambdaB=0.8) 
#when all genes are used 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bhat_l1l21_D5_Scerevisiae=ttt$Bhat_l1l21
rownames(Bhat_l1l21_D5_Scerevisiae) <- colnames(X)
colnames(Bhat_l1l21_D5_Scerevisiae) <- colnames(Y)

Omegahat_l1l21_D5_Scerevisiae=ttt$Omegahat_l1l21
rownames(Omegahat_l1l21_D5_Scerevisiae) <- colnames(Y)
colnames(Omegahat_l1l21_D5_Scerevisiae) <- colnames(Y)

Sigmahat_l1l21_D5_Scerevisiae=ttt$Sigmahat_l1l21
rownames(Sigmahat_l1l21_D5_Scerevisiae) <- colnames(Y)
colnames(Sigmahat_l1l21_D5_Scerevisiae) <- colnames(Y)

