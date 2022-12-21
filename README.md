# mixed-norms-GRN
This is the repository for the manuscript "Gene regulatory network inference using mixed-norms regularized multivariate model with covariance selection" by Alain J. Mbebi & Zoran Nikoloski.

# Organisation

1. The folder Codes contains the following R scripts with the K-folds cross-validation option to learn the hyperparameters:
  * Mixed_L1L21_GRN.R which computes L1L21-solution 
  * Mixed_L1L21G_GRN.R which computes L1L21G-solution
  * Mixed_L2L21_GRN.R which computes L2L21-solution
  * Mixed_L2L21G_GRN.R which computes L2L21G-solution
  * L1L21_Dream5_Scerevisiae_example_run.R is an example run using the L1L21-solution with S. cerevisiae data (Network 4 in DREAM5 challenge) 

2. The folder Figures contains all figures in the manuscript.

3. The folder Inferred-networks contains all network objects for each dataset and each inference methods in the comparative analysis.

# Dependencies and required packages
The following packages are required for the contending approaches in the comparative analysis: "devtools", "foreach", "plyr", "doRNG", "glmnet" and "randomForest".

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
# ANOVerence
Because of some technical issue (e.g code's accessibility), we were not able to reproduce ANOVerence results and used the inferred network from DREAM5 challenge instead.

4. Although the codes here were tested on Fedora 29 (Workstation Edition) using R (version 3.6.1), they can run under any Linux or Windows OS distributions, as long as all the required packages are compatible with the desired R version.

