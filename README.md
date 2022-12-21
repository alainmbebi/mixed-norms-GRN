This is the source code for the manuscript "Gene regulatory network inference using mixed-norms regularized multivariate model with covariance selection", by Alain J. Mbebi & Zoran Nikoloski.
# mixed-norms-GRN

R functions to implement all algorithms described in Gene regulatory network inference using mixed-norms regularized multivariate model with covariance selection.

1. The folder Codes contains the following R scripts with the K-folds cross-validation option to learn the hyperparameters:
  * Mixed_L1L21_GRN.R which computes L1L21-solution 
  * Mixed_L1L21G_GRN.R which computes L1L21G-solution
  * Mixed_L2L21_GRN.R which computes L2L21-solution
  * Mixed_L2L21G_GRN.R which computes L2L21G-solution

2. L1L21_Dream5_Scerevisiae_example_run.R is an example run using the L1L21-solution with S. cerevisiae data (Network 4 in DREAM5 challenge) 


3. Although the codes here were tested on Fedora 29 (Workstation Edition) using R (version 3.6.1), they can run under any Linux or Windows OS distributions, as long as all the required packages are compatible with the desired R version.

