
source("Mixed_L1L21_GRN.R")
source("Mixed_L1L21G_GRN.R")
source("Mixed_L2L21_GRN.R")
source("Mixed_L2L21G_GRN.R")

##### read the expression profile matrix
expr <- as.matrix(read.table(file = "net1_expression_data.tsv",header = T,sep = "\t",stringsAsFactors = F))
expr=scale(expr)
##### read transcription factors
TFs <- read.table(file = "net1_transcription_factors.tsv",header = F,sep = "\t",stringsAsFactors = F)
TFs <- TFs$V1)
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
#uncomment if the goal is to learn new interactions
#TGsexpr <- expr
#TFsexpr <- expr[,colnames(expr)%in%TFs]

### for inference only: select expression profiles for only GS_genes
TGsexpr <- expr[,which(colnames(expr)%in%GS_genes)]
TFsexpr <- TGsexpr[,colnames(TGsexpr)%in%TFs]
TFs=colnames(TFsexpr)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

itermax=50
Yl21=TGsexpr 
Xl21=TFsexpr
### remove duplicate colnames (if any) in a matrix Yl21 and Xl21 with
Yl21 <- Yl21[, !duplicated(colnames(Yl21))]
Xl21 <- Xl21[, !duplicated(colnames(Xl21))]
##lambdaO=2^(-6) and lambdaB=.6 reproduce the results in the manuscript and they are learned from CV
ttt=Mixedl1l21(Xl21, Yl21, lambdaO=2^(-6), lambdaB=.6) 
#Get the regression
Bhat_l1l21_D5_Insilico=ttt$Bhat_l1l21
rownames(Bhat_l1l21_D5_Insilico) <- colnames(Xl21)
colnames(Bhat_l1l21_D5_Insilico) <- colnames(Yl21)
#Get the precision matrix
Omegahat_l1l21_D5_Insilico=ttt$Omegahat_l1l21
rownames(Omegahat_l1l21_D5_Insilico) <- colnames(Yl21)
colnames(Omegahat_l1l21_D5_Insilico) <- colnames(Yl21)
#Get the covariance matrix
Sigmahat_l1l21_D5_Insilico=ttt$Sigmahat_l1l21
rownames(Sigmahat_l1l21_D5_Insilico) <- colnames(Yl21)
colnames(Sigmahat_l1l21_D5_Insilico) <- colnames(Yl21)
