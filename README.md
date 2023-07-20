# mixed-norms-GRN
This is the repository for the manuscript "Gene regulatory network inference using mixed-norms regularized multivariate model with covariance selection" by Alain J. Mbebi & Zoran Nikoloski.

# Organisation

1. The folder Codes contains the following R scripts with the K-folds cross-validation option to learn the hyperparameters:
  * Mixed_L1L21_GRN.R which computes L1L21-solution 
  * Mixed_L1L21G_GRN.R which computes L1L21G-solution
  * Mixed_L2L21_GRN.R which computes L2L21-solution
  * Mixed_L2L21G_GRN.R which computes L2L21G-solution
  * L1L21_Dream5_Scerevisiae_example_run.R is an example run using the L1L21-solution with S. cerevisiae data (Network 4 in DREAM5 challenge) 
  All files needed to successfully run "L1L21_Dream5_Scerevisiae_example_run" are locaded in the folder Codes.

2. The folder Figures contains all figures in the manuscript.

3. The folder Inferred-networks contains all network objects for each dataset and each inference methods in the comparative analysis.

4. Due to the space limitation on this repository, all data along with predicted networks for the contending approaches have been made publicly available and can be accessed from https://zenodo.org/record/7965949. 

# Dependencies and required packages
The following packages are required for the contending approaches in the comparative analysis: "devtools", "foreach", "plyr", "glmnet" and "randomForest".

# GENIE3
The GENIE3 package can be installed from: http://bioconductor.org/packages/release/bioc/html/GENIE3.html

# TIGRESS
The TIGRESS repository can be obtained from: https://github.com/jpvert/tigress

# ENNET
The ENNET repository can be obtained from: https://github.com/slawekj/ennet

# PLSNET
The Matlab source code of PLSNET can be obtained from: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1398-6#Sec17

# PORTIA 
The PORTIA repository can be obtained from: https://github.com/AntoinePassemiers/PORTIA

# D3GRN
The Matlab source code of D3GRN can be obtained from: https://github.com/chenxofhit/D3GRN

# Fused-LASSO
The fused-LASSO repository can be obtained from: https://github.com/omranian/inference-of-GRN-using-Fused-LASSO

# ANOVerence
Because of some technical issues (e.g code's accessibility:  http://www2.bio.ifi.lmu.de/˜kueffner/anova.tar.gz), we were not able to reproduce ANOVerence results and used the inferred network from DREAM5 challenge instead.

4. Although the codes here were tested on Fedora 29 (Workstation Edition) using R (version 4.2.2), they can run under any Linux or Windows OS distributions, as long as all the required packages are compatible with the desired R version.
5. All data underlying this publication are publicly available and their corresponding references provided in the article. They can also be downloaded from: https://doi.org/10.5281/zenodo.7965949
