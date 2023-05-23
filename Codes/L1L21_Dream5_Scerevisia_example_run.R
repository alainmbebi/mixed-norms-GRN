#This is an example run for L1L21 solution with S. cerevisiae data (network 4 from DREAM5 challenge)
#the function "ordernet" is adapted from the source code in Fused LASSO (https://github.com/omranian/inference-of-GRN-using-Fused-LASSO#start-of-content)
source("Mixed_L1L21_GRN.R")

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
expr <- as.matrix(read.table(file = "net4_expression_data.tsv",header = T,sep = "\t",stringsAsFactors = F))
head(expr[,1:6])
dim(expr)
expr=scale(expr)
##### read transcription factors
TFs <- read.table(file = "net4_transcription_factors.tsv",header = F,sep = "\t",stringsAsFactors = F)
TFs <- TFs$V1
head(TFs)
length(TFs)
##### read Gold standard network with all interactions
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
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Rescale each row of the regression coefficient matrix to range between 0 and 1
Bhat_l1l21_D5_Scerevisiae=t(apply(abs(Bhat_l1l21_D5_Scerevisiae), 1, function(x)(x-min(x))/(max(x)-min(x))))
mixedl1l21_D5_Scerevisiae <- reshape2::melt(Bhat_l1l21_D5_Scerevisiae)
mixedl1l21_D5_Scerevisiae <- mixedl1l21_D5_Scerevisiae[,c(2,1,3)]
colnames(mixedl1l21_D5_Scerevisiae) <- c("from","to","weight_mixedl1l21_D5_Scerevisiae")
mixedl1l21_D5_Scerevisiae=rm_NTF_int(mixedl1l21_D5_Scerevisiae,TFs)
mixedl1l21_D5_Scerevisiae <- ordernet(mixedl1l21_D5_Scerevisiae,TFs)
mixedl1l21_D5_Scerevisiae$score <- abs(mixedl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae)/max(abs(mixedl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae))
mixedl1l21_D5_Scerevisiae <- mixedl1l21_D5_Scerevisiae[order(mixedl1l21_D5_Scerevisiae$score,decreasing = T),]
mixedl1l21_D5_Scerevisiae=head(mixedl1l21_D5_Scerevisiae,100000)

colnames(mixedl1l21_D5_Scerevisiae)[3:4] <-  c("weight_mixedl1l21_D5_Scerevisiae","score_mixedl1l21_D5_Scerevisiae")
write.table(mixedl1l21_D5_Scerevisiae, file = "mixedl1l21_D5_Scerevisiae.txt",row.names = FALSE, col.names = TRUE)

nl1l21_D5_Scerevisiae <- merge(mixedl1l21_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
nl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae[which(is.na(nl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae)==T)] <- 0
nl1l21_D5_Scerevisiae$score_mixedl1l21_D5_Scerevisiae[which(is.na(nl1l21_D5_Scerevisiae$score_mixedl1l21_D5_Scerevisiae)==T)] <- 0
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# get the metrics (AUROC, AUPR) with precrec package
library(precrec)
sscurves_nl1l21_D5_Scerevisiae <- evalmod(scores = nl1l21_D5_Scerevisiae$score_mixedl1l21_D5_Scerevisiae, labels = GS$weight_GS)
plot(sscurves_nl1l21_D5_Scerevisiae)
auc(sscurves_nl1l21_D5_Scerevisiae)

############################################
####Compute the early PR
# make sure there is no self-loops in gold standard network
GS <- GS[GS$from != GS$to,]
# Find the number of edges in gold standard network
k <- 3940 #nrow(GS) #with edges only
GS <- GS[order(GS$weight_GS, decreasing = T),]
GS <- GS[1:k,]
#####EPR
############################################l1l21
# sort predicted network by decreasing interaction strength
mixedl1l21_D5_Scerevisiae <- mixedl1l21_D5_Scerevisiae[order(mixedl1l21_D5_Scerevisiae$score_mixedl1l21_D5_Scerevisiae ,decreasing = T),]
# keep only top-k predicted edges
mixedl1l21_D5_Scerevisiae <- mixedl1l21_D5_Scerevisiae[1:k,]
# count number of true positives in top-k predicted edges
num_true_positives_mixedl1l21_D5_Scerevisiae <- sum(mixedl1l21_D5_Scerevisiae$from %in% GS$from & mixedl1l21_D5_Scerevisiae$to %in% GS$to)
# compute early precision ratio
early_precision_ratio_mixedl1l21_D5_Scerevisiae<- (num_true_positives_mixedl1l21_D5_Scerevisiae / k)*100

############################################
######## compute the normalized Discounted Cumulative Gain (nDCG)
# load required libraries
library(dplyr)
#Load the gold-standard network (with only positive interactions) and inferred network data frames
positive_interactions <- read.table(file = "DREAM5_NetworkInference_Edges_Network4.tsv",header = F,sep = "\t",stringsAsFactors = F)
colnames(positive_interactions) <- c( "from","to",paste("weight","GS",sep = "_"))

# filter the inferred network to only include links present in the gold-standard network
mixedl1l21_D5_Scerevisiae <- mixedl1l21_D5_Scerevisiae[,c(1,2,4)]   #select the from, to and score columns
mixedl1l21_D5_Scerevisiae <- mixedl1l21_D5_Scerevisiae %>% filter(mixedl1l21_D5_Scerevisiae$from %in% positive_interactions$from & to %in% positive_interactions$to)

# compute the nDCG
relevant_genes <- unique(positive_interactions$from) # identify the relevant genes in the gold-standard network
ndcg_mixedl1l21_D5_Scerevisiae <- data.frame(gene = character(), ndcg_mixedl1l21_D5_Scerevisiae = double(), stringsAsFactors = FALSE) # initialize the data frame to store the NDCG scores
for (gene in relevant_genes) {
  # subset the gold-standard and inferred networks to the relevant gene
  goldstandard_sub <- positive_interactions %>% filter(positive_interactions$from == gene)
  inferred_sub <- mixedl1l21_D5_Scerevisiae %>% filter(mixedl1l21_D5_Scerevisiae$from == gene)
  
  # join the two data frames on the "to" column to get the relevance scores
  joined <- left_join(inferred_sub, goldstandard_sub, by = "to")
  joined$weight_GS[is.na(joined$weight_GS)] <- 0 # set missing values to 0
  joined$relevance <- joined$weight_GS # set the relevance scores to the presence values
  
  # calculate the DCG and IDCG scores
  joined <- joined %>% arrange(desc(score_mixedl1l21_D5_Scerevisiae))
  joined$discount <- log2(1 + seq_along(joined$score_mixedl1l21_D5_Scerevisiae))
  dcg_mixedl1l21_D5_Scerevisiae <- sum(joined$relevance / joined$discount)
  idcg_mixedl1l21_D5_Scerevisiae <- sum(sort(joined$relevance, decreasing = TRUE) / log2(1 + seq_along(joined$relevance)))
  
  # calculate the NDCG score
  if (idcg_mixedl1l21_D5_Scerevisiae == 0) {
    ndcg_mixedl1l21_D5_Scerevisiae_score <- 0
  } else {
    ndcg_mixedl1l21_D5_Scerevisiae_score <- dcg_mixedl1l21_D5_Scerevisiae / idcg_mixedl1l21_D5_Scerevisiae
  }
  
  # add the NDCG score to the data frame
  ndcg_mixedl1l21_D5_Scerevisiae <- ndcg_mixedl1l21_D5_Scerevisiae %>% add_row(gene = gene, ndcg_mixedl1l21_D5_Scerevisiae = ndcg_mixedl1l21_D5_Scerevisiae_score)
}

# print the NDCG scores
print(ndcg_mixedl1l21_D5_Scerevisiae)
print(ndcg_mixedl1l21_D5_Scerevisiae_score)

############################################
####fonction to compute tau index
#This function is adapted from the code provided in https://doi.org/10.1093/bib/bbw008 (A benchmark of gene expression tissue-specificity metrics)
fTau_matrix <- function(x) {
  res <- numeric(nrow(x))
  for (i in 1:nrow(x)) {
    if(all(!is.na(x[i,]))) {
      if(min(x[i,], na.rm=TRUE) >= 0) {
        if(max(x[i,]) != 0) {
          x_norm <- (1 - x[i,] / max(x[i,]))
          res[i] <- sum(x_norm, na.rm=TRUE) / (length(x_norm) - 1)
        } else {
          res[i] <- 0
        }
      } else {
        res[i] <- NA
        #print("Expression values have to be positive!")
      }
    } else {
      res[i] <- NA
      #print("No data for this gene avalable.")
    }
  }
  names(res) <- rownames(x)
  return(res)
}

################################To obtain the figure in the manuscript

#L1L21
Bhat_l1l21_D5_Scerevisiae<-t(as.matrix(read.csv("Bhat_l1l21_D5_Scerevisiae.csv", head=TRUE, row.names = 1, sep=","))) 
Bhat_l1l21_D5_Scerevisiae=t(apply(abs(Bhat_l1l21_D5_Scerevisiae), 1, function(x)(x-min(x))/(max(x)-min(x))))
mixedl1l21_D5_Scerevisiae <- reshape2::melt(Bhat_l1l21_D5_Scerevisiae)
mixedl1l21_D5_Scerevisiae <- mixedl1l21_D5_Scerevisiae[,c(2,1,3)]
colnames(mixedl1l21_D5_Scerevisiae) <- c("from","to","weight_mixedl1l21_D5_Scerevisiae")
mixedl1l21_D5_Scerevisiae=rm_NTF_int(mixedl1l21_D5_Scerevisiae,TFs)
mixedl1l21_D5_Scerevisiae <- ordernet(mixedl1l21_D5_Scerevisiae,TFs)
mixedl1l21_D5_Scerevisiae$score <- abs(mixedl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae)/max(abs(mixedl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae))
mixedl1l21_D5_Scerevisiae <- mixedl1l21_D5_Scerevisiae[order(mixedl1l21_D5_Scerevisiae$score,decreasing = T),]
mixedl1l21_D5_Scerevisiae=head(mixedl1l21_D5_Scerevisiae,100000)

colnames(mixedl1l21_D5_Scerevisiae)[3:4] <-  c("weight_mixedl1l21_D5_Scerevisiae","score_mixedl1l21_D5_Scerevisiae")

nl1l21_D5_Scerevisiae <- merge(mixedl1l21_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
nl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae[which(is.na(nl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae)==T)] <- 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#L1L21G
Bhat_l1l21_Gauss_D5_Scerevisiae<-t(as.matrix(read.csv("Bhat_l1l21_Gauss_D5_Scerevisiae.csv", head=TRUE, row.names = 1, sep=","))) 
Bhat_l1l21_Gauss_D5_Scerevisiae=t(apply(abs(Bhat_l1l21_Gauss_D5_Scerevisiae), 1, function(x)(x-min(x))/(max(x)-min(x))))
mixedl1l21_Gauss_D5_Scerevisiae <- reshape2::melt(Bhat_l1l21_Gauss_D5_Scerevisiae)
mixedl1l21_Gauss_D5_Scerevisiae <- mixedl1l21_Gauss_D5_Scerevisiae[,c(2,1,3)]
colnames(mixedl1l21_Gauss_D5_Scerevisiae) <- c("from","to","weight_mixedl1l21_Gauss_D5_Scerevisiae")
mixedl1l21_Gauss_D5_Scerevisiae=rm_NTF_int(mixedl1l21_Gauss_D5_Scerevisiae,TFs)
mixedl1l21_Gauss_D5_Scerevisiae <- ordernet(mixedl1l21_Gauss_D5_Scerevisiae,TFs)
mixedl1l21_Gauss_D5_Scerevisiae$score <- abs(mixedl1l21_Gauss_D5_Scerevisiae$weight_mixedl1l21_Gauss_D5_Scerevisiae)/max(abs(mixedl1l21_Gauss_D5_Scerevisiae$weight_mixedl1l21_Gauss_D5_Scerevisiae))
mixedl1l21_Gauss_D5_Scerevisiae <- mixedl1l21_Gauss_D5_Scerevisiae[order(mixedl1l21_Gauss_D5_Scerevisiae$score,decreasing = T),]
mixedl1l21_Gauss_D5_Scerevisiae=head(mixedl1l21_Gauss_D5_Scerevisiae,100000)
colnames(mixedl1l21_Gauss_D5_Scerevisiae)[3:4] <-  c("weight_mixedl1l21_Gauss_D5_Scerevisiae","score_mixedl1l21_Gauss_D5_Scerevisiae")
nl1l21_Gauss_D5_Scerevisiae <- merge(mixedl1l21_Gauss_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
nl1l21_Gauss_D5_Scerevisiae$weight_mixedl1l21_Gauss_D5_Scerevisiae[which(is.na(nl1l21_Gauss_D5_Scerevisiae$weight_mixedl1l21_Gauss_D5_Scerevisiae)==T)] <- 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#L2L21
Bhat_l2l21_D5_Scerevisiae<-t(as.matrix(read.csv("Bhat_l2l21_D5_Scerevisiae.csv", head=TRUE, row.names = 1, sep=","))) 
Bhat_l2l21_D5_Scerevisiae=t(apply(abs(Bhat_l2l21_D5_Scerevisiae), 1, function(x)(x-min(x))/(max(x)-min(x))))
mixedl2l21_D5_Scerevisiae <- reshape2::melt(Bhat_l2l21_D5_Scerevisiae)
mixedl2l21_D5_Scerevisiae <- mixedl2l21_D5_Scerevisiae[,c(2,1,3)]
colnames(mixedl2l21_D5_Scerevisiae) <- c("from","to","weight_mixedl2l21_D5_Scerevisiae")
mixedl2l21_D5_Scerevisiae=rm_NTF_int(mixedl2l21_D5_Scerevisiae,TFs)
mixedl2l21_D5_Scerevisiae <- ordernet(mixedl2l21_D5_Scerevisiae,TFs)
mixedl2l21_D5_Scerevisiae$score <- abs(mixedl2l21_D5_Scerevisiae$weight_mixedl2l21_D5_Scerevisiae)/max(abs(mixedl2l21_D5_Scerevisiae$weight_mixedl2l21_D5_Scerevisiae))
mixedl2l21_D5_Scerevisiae <- mixedl2l21_D5_Scerevisiae[order(mixedl2l21_D5_Scerevisiae$score,decreasing = T),]
mixedl2l21_D5_Scerevisiae=head(mixedl2l21_D5_Scerevisiae,100000)
colnames(mixedl2l21_D5_Scerevisiae)[3:4] <-  c("weight_mixedl2l21_D5_Scerevisiae","score_mixedl2l21_D5_Scerevisiae")
nl2l21_D5_Scerevisiae <- merge(mixedl2l21_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
nl2l21_D5_Scerevisiae$weight_mixedl2l21_D5_Scerevisiae[which(is.na(nl2l21_D5_Scerevisiae$weight_mixedl2l21_D5_Scerevisiae)==T)] <- 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#L2L21G
Bhat_l2l21_Gauss_D5_Scerevisiae<-t(as.matrix(read.csv("Bhat_l2l21_Gauss_D5_Scerevisiae.csv", head=TRUE, row.names = 1, sep=","))) 
Bhat_l2l21_Gauss_D5_Scerevisiae=t(apply(abs(Bhat_l2l21_Gauss_D5_Scerevisiae), 1, function(x)(x-min(x))/(max(x)-min(x))))
mixedl2l21_Gauss_D5_Scerevisiae <- reshape2::melt(Bhat_l2l21_Gauss_D5_Scerevisiae)
mixedl2l21_Gauss_D5_Scerevisiae <- mixedl2l21_Gauss_D5_Scerevisiae[,c(2,1,3)]
colnames(mixedl2l21_Gauss_D5_Scerevisiae) <- c("from","to","weight_mixedl2l21_Gauss_D5_Scerevisiae")
mixedl2l21_Gauss_D5_Scerevisiae=rm_NTF_int(mixedl2l21_Gauss_D5_Scerevisiae,TFs)
mixedl2l21_Gauss_D5_Scerevisiae <- ordernet(mixedl2l21_Gauss_D5_Scerevisiae,TFs)
mixedl2l21_Gauss_D5_Scerevisiae$score <- abs(mixedl2l21_Gauss_D5_Scerevisiae$weight_mixedl2l21_Gauss_D5_Scerevisiae)/max(abs(mixedl2l21_Gauss_D5_Scerevisiae$weight_mixedl2l21_Gauss_D5_Scerevisiae))
mixedl2l21_Gauss_D5_Scerevisiae <- mixedl2l21_Gauss_D5_Scerevisiae[order(mixedl2l21_Gauss_D5_Scerevisiae$score,decreasing = T),]
mixedl2l21_Gauss_D5_Scerevisiae=head(mixedl2l21_Gauss_D5_Scerevisiae,100000)
colnames(mixedl2l21_Gauss_D5_Scerevisiae)[3:4] <-  c("weight_mixedl2l21_Gauss_D5_Scerevisiae","score_mixedl2l21_Gauss_D5_Scerevisiae")
nl2l21_Gauss_D5_Scerevisiae <- merge(mixedl2l21_Gauss_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
nl2l21_Gauss_D5_Scerevisiae$weight_mixedl2l21_Gauss_D5_Scerevisiae[which(is.na(nl2l21_Gauss_D5_Scerevisiae$weight_mixedl2l21_Gauss_D5_Scerevisiae)==T)] <- 0
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#GENIE3
Genie3_D5_Scerevisiae <- read.table(file = "Genie3_D5_Scerevisiae.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(Genie3_D5_Scerevisiae) <- c("from","to","weight_Genie3_D5_Scerevisiae")
Genie3_D5_Scerevisiae$score <- abs(Genie3_D5_Scerevisiae$weight_Genie3_D5_Scerevisiae)/max(abs(Genie3_D5_Scerevisiae$weight_Genie3_D5_Scerevisiae))
Genie3_D5_Scerevisiae <- Genie3_D5_Scerevisiae[order(Genie3_D5_Scerevisiae$score,decreasing = T),]
colnames(Genie3_D5_Scerevisiae)[3:4] <-  c("weight_Genie3_D5_Scerevisiae","score_Genie3_D5_Scerevisiae")

Genie3_D5_Scerevisiae <- merge(Genie3_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
Genie3_D5_Scerevisiae$weight_Genie3_D5_Scerevisiae[which(is.na(Genie3_D5_Scerevisiae$weight_Genie3_D5_Scerevisiae)==T)] <- 0
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ANOVA
Anova_D5_Scerevisiae <- read.table(file = "Anova_D5_Scerevisiae.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(Anova_D5_Scerevisiae) <- c("from","to","weight_Anova_D5_Scerevisiae")
Anova_D5_Scerevisiae$score <- abs(Anova_D5_Scerevisiae$weight_Anova_D5_Scerevisiae)/max(abs(Anova_D5_Scerevisiae$weight_Anova_D5_Scerevisiae))
Anova_D5_Scerevisiae <- Anova_D5_Scerevisiae[order(Anova_D5_Scerevisiae$score,decreasing = T),]
colnames(Anova_D5_Scerevisiae)[3:4] <-  c("weight_Anova_D5_Scerevisiae","score_Anova_D5_Scerevisiae")
Anova_D5_Scerevisiae <- merge(Anova_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
Anova_D5_Scerevisiae$weight_Anova_D5_Scerevisiae[which(is.na(Anova_D5_Scerevisiae$weight_Anova_D5_Scerevisiae)==T)] <- 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TIGRESS
Tigress_D5_Scerevisiae <- read.table(file = "Tigress_D5_Scerevisiae.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(Tigress_D5_Scerevisiae) <- c("from","to","weight_Tigress_D5_Scerevisiae")
Tigress_D5_Scerevisiae$score <- abs(Tigress_D5_Scerevisiae$weight_Tigress_D5_Scerevisiae)/max(abs(Tigress_D5_Scerevisiae$weight_Tigress_D5_Scerevisiae))
Tigress_D5_Scerevisiae <- Tigress_D5_Scerevisiae[order(Tigress_D5_Scerevisiae$score,decreasing = T),]
colnames(Tigress_D5_Scerevisiae)[3:4] <-  c("weight_Tigress_D5_Scerevisiae","score_Tigress_D5_Scerevisiae")
Tigress_D5_Scerevisiae <- merge(Tigress_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
Tigress_D5_Scerevisiae$weight_Tigress_D5_Scerevisiae[which(is.na(Tigress_D5_Scerevisiae$weight_Tigress_D5_Scerevisiae)==T)] <- 0
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#PORTIA
Portia_D5_Scerevisiae <- read.table(file = "Portia_D5_Scerevisiae.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(Portia_D5_Scerevisiae) <- c("from","to","weight_Portia_D5_Scerevisiae")
Portia_D5_Scerevisiae$score <- abs(Portia_D5_Scerevisiae$weight_Portia_D5_Scerevisiae)/max(abs(Portia_D5_Scerevisiae$weight_Portia_D5_Scerevisiae))
Portia_D5_Scerevisiae <- Portia_D5_Scerevisiae[order(Portia_D5_Scerevisiae$score,decreasing = T),]
colnames(Portia_D5_Scerevisiae)[3:4] <-  c("weight_Portia_D5_Scerevisiae","score_Portia_D5_Scerevisiae")
Portia_D5_Scerevisiae <- merge(Portia_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
Portia_D5_Scerevisiae$weight_Portia_D5_Scerevisiae[which(is.na(Portia_D5_Scerevisiae$weight_Portia_D5_Scerevisiae)==T)] <- 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
etePortia_D5_Scerevisiae <- read.table(file = "etePortia_D5_Scerevisiae.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(etePortia_D5_Scerevisiae) <- c("from","to","weight_etePortia_D5_Scerevisiae")
etePortia_D5_Scerevisiae$score <- abs(etePortia_D5_Scerevisiae$weight_etePortia_D5_Scerevisiae)/max(abs(etePortia_D5_Scerevisiae$weight_etePortia_D5_Scerevisiae))
etePortia_D5_Scerevisiae <- etePortia_D5_Scerevisiae[order(etePortia_D5_Scerevisiae$score,decreasing = T),]
colnames(etePortia_D5_Scerevisiae)[3:4] <-  c("weight_etePortia_D5_Scerevisiae","score_etePortia_D5_Scerevisiae")
etePortia_D5_Scerevisiae <- merge(etePortia_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
etePortia_D5_Scerevisiae$weight_etePortia_D5_Scerevisiae[which(is.na(etePortia_D5_Scerevisiae$weight_etePortia_D5_Scerevisiae)==T)] <- 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Ennet_D5_Scerevisiae <- read.table(file = "Ennet_D5_Scerevisiae.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(Ennet_D5_Scerevisiae) <- c("from","to","weight_Ennet_D5_Scerevisiae")
Ennet_D5_Scerevisiae$score <- abs(Ennet_D5_Scerevisiae$weight_Ennet_D5_Scerevisiae)/max(abs(Ennet_D5_Scerevisiae$weight_Ennet_D5_Scerevisiae))
Ennet_D5_Scerevisiae <- Ennet_D5_Scerevisiae[order(Ennet_D5_Scerevisiae$score,decreasing = T),]
colnames(Ennet_D5_Scerevisiae)[3:4] <-  c("weight_Ennet_D5_Scerevisiae","score_Ennet_D5_Scerevisiae")
Ennet_D5_Scerevisiae <- merge(Ennet_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
Ennet_D5_Scerevisiae$weight_Ennet_D5_Scerevisiae[which(is.na(Ennet_D5_Scerevisiae$weight_Ennet_D5_Scerevisiae)==T)] <- 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Plsnet_D5_Scerevisiae <- read.table(file = "Plsnet_D5_Scerevisiae.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(Plsnet_D5_Scerevisiae) <- c("from","to","weight_Plsnet_D5_Scerevisiae")
Plsnet_D5_Scerevisiae$score <- abs(Plsnet_D5_Scerevisiae$weight_Plsnet_D5_Scerevisiae)/max(abs(Plsnet_D5_Scerevisiae$weight_Plsnet_D5_Scerevisiae))
Plsnet_D5_Scerevisiae <- Plsnet_D5_Scerevisiae[order(Plsnet_D5_Scerevisiae$score,decreasing = T),]
colnames(Plsnet_D5_Scerevisiae)[3:4] <-  c("weight_Plsnet_D5_Scerevisiae","score_Plsnet_D5_Scerevisiae")
Plsnet_D5_Scerevisiae <- merge(Plsnet_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
Plsnet_D5_Scerevisiae$weight_Plsnet_D5_Scerevisiae[which(is.na(Plsnet_D5_Scerevisiae$weight_Plsnet_D5_Scerevisiae)==T)] <- 0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
d3grn_D5_Scerevisiae <- read.table(file = "d3grn_D5_Scerevisiae.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(d3grn_D5_Scerevisiae) <- c("from","to","weight_d3grn_D5_Scerevisiae")
d3grn_D5_Scerevisiae$score <- abs(d3grn_D5_Scerevisiae$weight_d3grn_D5_Scerevisiae)/max(abs(d3grn_D5_Scerevisiae$weight_d3grn_D5_Scerevisiae))
d3grn_D5_Scerevisiae <- d3grn_D5_Scerevisiae[order(d3grn_D5_Scerevisiae$score,decreasing = T),]
colnames(d3grn_D5_Scerevisiae)[3:4] <-  c("weight_d3grn_D5_Scerevisiae","score_d3grn_D5_Scerevisiae")
d3grn_D5_Scerevisiae <- merge(d3grn_D5_Scerevisiae,GS,by = c("from","to"),all.y = T)
d3grn_D5_Scerevisiae$weight_d3grn_D5_Scerevisiae[which(is.na(d3grn_D5_Scerevisiae$weight_d3grn_D5_Scerevisiae)==T)] <- 0
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##plot with package pROC
library(pROC)
library("RColorBrewer")
library(viridis)
library(ggplot2)
#par(bg="transparent")
pROC_obj_l1l21  <-roc(GS$weight_GS, nl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae)
pROC_obj_l1l21G <-roc(GS$weight_GS, nl1l21_Gauss_D5_Scerevisiae$weight_mixedl1l21_Gauss_D5_Scerevisiae)
pROC_obj_l2l21  <-roc(GS$weight_GS, nl2l21_D5_Scerevisiae$weight_mixedl2l21_D5_Scerevisiae)
pROC_obj_l2l21G <-roc(GS$weight_GS, nl2l21_Gauss_D5_Scerevisiae$weight_mixedl2l21_Gauss_D5_Scerevisiae)
pROC_obj_genie3 <-roc(GS$weight_GS, Genie3_D5_Scerevisiae$weight_Genie3_D5_Scerevisiae)
pROC_obj_anova  <-roc(GS$weight_GS, Anova_D5_Scerevisiae$weight_Anova_D5_Scerevisiae)
pROC_obj_tigress<-roc(GS$weight_GS, Tigress_D5_Scerevisiae$weight_Tigress_D5_Scerevisiae)
pROC_obj_Portia<-roc(GS$weight_GS, Portia_D5_Scerevisiae$weight_Portia_D5_Scerevisiae)
pROC_obj_etePortia<-roc(GS$weight_GS, etePortia_D5_Scerevisiae$weight_etePortia_D5_Scerevisiae)
pROC_obj_Ennet<-roc(GS$weight_GS, Ennet_D5_Scerevisiae$weight_Ennet_D5_Scerevisiae)
pROC_obj_Plsnet<-roc(GS$weight_GS, Plsnet_D5_Scerevisiae$weight_Plsnet_D5_Scerevisiae)
pROC_obj_d3grn<-roc(GS$weight_GS, d3grn_D5_Scerevisiae$weight_d3grn_D5_Scerevisiae)
#-----------------------------------------------------------------------------------
rocobj <- list(Data_1 = pROC_obj_l1l21, Data_2 = pROC_obj_l1l21G, Data_3 = pROC_obj_l2l21, Data_4 = pROC_obj_l2l21G,
               Data_5=pROC_obj_genie3, Data_6=pROC_obj_anova, Data_7=pROC_obj_tigress, Data_8=pROC_obj_Portia,
               Data_9=pROC_obj_etePortia, Data_10=pROC_obj_Ennet, Data_11=pROC_obj_Plsnet, Data_12=pROC_obj_d3grn)

##############################################drawing plot
#set colors
Alain12 <- c(
  "maroon", 
  "blue1",
  "steelblue4",
  "darkturquoise", 
  "yellow3",
  "brown",
  "dodgerblue2", 
  "#E31A1C",
  "green4",
  "gold1",
  "skyblue2", 
  "darkorange4"
)
tiff("ROC_Scerevisiae_Dream5.tiff", units="in", width=5, height=5, res=300)
ggroc(rocobj,aes=c("linetype", "color"), linetype = 2, size = .2, alpha=1, show.legend = FALSE)+
  labs(x = "False positive rate", y = "") + 
  geom_abline(intercept = 1, slope = 1, color='grey',size = 0.5,linetype = 1) + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(.8, 0.45),
        legend.box = "vertical",
        legend.key = element_rect(colour = NA, fill = NA),
        legend.title = element_text(colour="black",size=16),
        axis.line = element_line(colour = "black", 
                                 size = 0.3, linetype = "solid"),
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="white"),
        legend.text=element_text(size=16),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(colour="black",size=14),
        axis.text.y = element_text(colour="black",size=14),
        axis.title.x = element_text(colour="black",size=20),
        axis.title.y = element_text(colour="black",size=20)) +
  scale_y_continuous(expand = c(0.01, 0.001),breaks= seq(from = 0 , to = 1 ,by = 0.2))+
  scale_x_reverse(expand = c(0.01, 0.001),breaks= seq(from = 0 , to = 1 ,by = 0.2),limits = c(1,-0.01),
                  labels = c("1.0","0.8","0.6","0.4","0.2","0.0")) +
  scale_color_manual(name = "", values=Alain12,
                     labels = c("L1L21", "L1L21G", "L2L21", "L2L21G", "GENIE3", "ANOVerence", "TIGRESS", "PORTA", "etePORTA", "ENNET", "PLSNET", "D3GRN"))
dev.off()
#-----------------------------------------------------------------------------------------------------------
#for PR curve use 
pPR_obj_l1l21  <-roc(GS$weight_GS, nl1l21_D5_Scerevisiae$weight_mixedl1l21_D5_Scerevisiae)
pPR_obj_l1l21G <-roc(GS$weight_GS, nl1l21_Gauss_D5_Scerevisiae$weight_mixedl1l21_Gauss_D5_Scerevisiae)
pPR_obj_l2l21  <-roc(GS$weight_GS, nl2l21_D5_Scerevisiae$weight_mixedl2l21_D5_Scerevisiae)
pPR_obj_l2l21G <-roc(GS$weight_GS, nl2l21_Gauss_D5_Scerevisiae$weight_mixedl2l21_Gauss_D5_Scerevisiae)
pPR_obj_genie3 <-roc(GS$weight_GS, Genie3_D5_Scerevisiae$weight_Genie3_D5_Scerevisiae)
pPR_obj_anova  <-roc(GS$weight_GS, Anova_D5_Scerevisiae$weight_Anova_D5_Scerevisiae)
pPR_obj_tigress<-roc(GS$weight_GS, Tigress_D5_Scerevisiae$weight_Tigress_D5_Scerevisiae)
pPR_obj_portia<-roc(GS$weight_GS, Portia_D5_Scerevisiae$weight_Portia_D5_Scerevisiae)
pPR_obj_eteportia<-roc(GS$weight_GS, etePortia_D5_Scerevisiae$weight_etePortia_D5_Scerevisiae)
pPR_obj_ennet<-roc(GS$weight_GS, Ennet_D5_Scerevisiae$weight_Ennet_D5_Scerevisiae)
pPR_obj_plsnet<-roc(GS$weight_GS, Plsnet_D5_Scerevisiae$weight_Plsnet_D5_Scerevisiae)
pPR_obj_d3grn<-roc(GS$weight_GS, d3grn_D5_Scerevisiae$weight_d3grn_D5_Scerevisiae)

coor_pPR_obj_l1l21=coords(pPR_obj_l1l21, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_l1l21G=coords(pPR_obj_l1l21G, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_l2l21=coords(pPR_obj_l2l21, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_l2l21G=coords(pPR_obj_l2l21G, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_genie3=coords(pPR_obj_genie3, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_anova=coords(pPR_obj_anova, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_tigress=coords(pPR_obj_tigress, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_portia=coords(pPR_obj_portia, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_eteportia=coords(pPR_obj_eteportia, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_ennet=coords(pPR_obj_ennet, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_plsnet=coords(pPR_obj_plsnet, "all", ret = c("threshold", "precision", "recall"))
coor_pPR_obj_d3grn=coords(pPR_obj_d3grn, "all", ret = c("threshold", "precision", "recall"))
#coord can also give  ret = c("threshold", "sensitivity", "specificity", "precision", "recall")
###ploting
tiff("PR_Scerevisiae_Dream5.tiff", units="in", width=5, height=5, res=300)
plot(coor_pPR_obj_l1l21$precision ~ coor_pPR_obj_l1l21$recall, t(coords(pPR_obj_l1l21, "all", 
                                                                        ret = c("recall", "precision"))), lty = 2, type="l", xlab="", ylab="",
     main = expression(italic("\n \n S. cerevisiae")), col="maroon", bty="l", cex.axis = 1.2, cex.main = 2)
title(xlab = "Recall", ylab = "", line=2.5, cex.lab = 1.5, mgp=c(1,1,0))

lines(coor_pPR_obj_l1l21G$precision ~ coor_pPR_obj_l1l21G$recall, lty = 2, col="blue1")
lines(coor_pPR_obj_l2l21$precision ~ coor_pPR_obj_l2l21$recall, lty = 2, col="steelblue4")
lines(coor_pPR_obj_l2l21G$precision ~ coor_pPR_obj_l2l21G$recall, lty = 2, col="darkturquoise")
lines(coor_pPR_obj_genie3$precision ~ coor_pPR_obj_genie3$recall, lty = 2, col="yellow3")
lines(coor_pPR_obj_anova$precision ~ coor_pPR_obj_anova$recall, lty = 2, col="brown")
lines(coor_pPR_obj_tigress$precision ~ coor_pPR_obj_tigress$recall, lty = 2, col="dodgerblue2")
lines(coor_pPR_obj_portia$precision ~ coor_pPR_obj_portia$recall, lty = 2, col="#E31A1C")
lines(coor_pPR_obj_eteportia$precision ~ coor_pPR_obj_eteportia$recall, lty = 2, col="green4")
lines(coor_pPR_obj_ennet$precision ~ coor_pPR_obj_ennet$recall, lty = 2, col="gold1")
lines(coor_pPR_obj_plsnet$precision ~ coor_pPR_obj_plsnet$recall, lty = 2, col="skyblue2")
lines(coor_pPR_obj_d3grn$precision ~ coor_pPR_obj_d3grn$recall, lty = 2, col="darkorange4")

dev.off()


