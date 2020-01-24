# CmmD: Constant Multiplex network Module Detector
CmmD is an tool for R that runs the MolTi community detection algorithm (Didier et al. 2018) and computes the "CmmD" analysis for the network community structures obtained. This "CmmD" analysis just looks for all the nodes of the multiplex network that are consistently found in the same community while changing MolTi's reslution parameter, a value that can be set in order to tune the number of communities found when running the algortihm.

# How to Install CmmD

# 1. Donwload and install molti-console in your system, available at: https://github.com/gilles-didier/MolTi-DREAM

# 2. Install package "devtools" in R:
      install.packages('devtools')
 
# 3.Install package "AnnotationDbi" using Bioconductor:
      install.packages('BiocManager') #Only if Bioconductor is not installed in your system
      BiocManager::install('AnnotationDbi')

# 4. Install CmmD package using devtools:
      devtools::install_github("ikernunezca/CmmD)

# Author and Mantainer: Iker Nuñez-Carpintero <iker.nunez@bsc.es>
