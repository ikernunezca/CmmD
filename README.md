# CmmD: Continual Multiplex network Module Detector
CmmD is an R tool that runs MolTi community detection algorithm (Didier et al. 2018) and computes a multilayer community trajectory analysis for the network community structures obtained. a multilayer community trajectory corresponds for all the nodes of the multiplex network that are consistently found in the same community while changing MolTi's reslution parameter, a value that can be set in order to tune the number of communities found when running the algorithm.

When using CmmD, please cite: Núñez-Carpintero, I. et al. (2021) ‘The multilayer community structure of medulloblastoma’, iScience, 24(4). doi: 10.1016/j.isci.2021.102365.

## A) How to Install CmmD

### 1. Download and install molti-console in your system, available at:
https://github.com/gilles-didier/MolTi-DREAM

* Make sure to add molti-console as an executable program from any directory by adding it to your ~/.bashrc file. (This file may change depending on your OS ditribution)

### 2. Install package "devtools" in R:
      install.packages('devtools')
 
### 3. Install package "AnnotationDbi" using Bioconductor:
      install.packages('BiocManager') #Only if Bioconductor is not installed in your system
      BiocManager::install('AnnotationDbi')

### 4. Install CmmD package using devtools:
      devtools::install_github('ikernunezca/CmmD')
      
## B) Using CmmD package functions:
Example 1: Let's imagine that we want to run MolTi and analyze the community structures given by MolTi at a determined range of resolution parameter: 0 to 30, in intervals of 0.5 (0, 0.5, 1, 1.5, 2, 2.5, etc). In order to do so, we call the function CmmD():
      
      #First, we generate a vector with the paths where our network files are saved: 
      nets <- c("Input_networks/net1.csv","Input_networks/net2.csv","Input_networks/net3.csv")
      
      #In order to run MolTi and analyze the outputs given by the program, we call the function CmmD():
      #If we want all the nodes of the network to be returned in the final output:
      CmmD(nodelist=NULL, input_layers= nets, resolution_start= 0, resolution_end= 30, interval= 0.5, destfile_community_analysis= "Output/",distmethod= "hamming",threads = 4)
      
      #If we want only a fixed set of nodes of the network to be returned in the final output:
      nodelist <- c("your_node1","your_node2,"your_node3")
      CmmD(nodelist=nets, input_layers= nets, resolution_start= 0, resolution_end= 30, interval= 0.5, destfile_community_analysis= "Output/",distmethod= "hamming",threads = 4)
      
Example 2: Now, suppose that we want to perform the very same analysis, but we have already precomputed community structures given by MolTi (or that have the same output file structure). In order to do so, we call the function CmmD_from_community_structures():
      
      # 1. Create a vector with the paths where the community structures are located:
      structures <- c("Molti_Output/0.5.csv","Molti_Output/1.csv")
      
      # 2. In order to analyze the range of resolution parameter 0 to 30, at intervals of 0.5 (e.g: c(0.5,1,1.5,2,2.5,3...,30)):
      CmmD_from_community_structures(nodelist=NULL, community_structures=estructuras,resolution_start=0,resolution_end=4, interval=0.5,distmethod= "hamming",threads = 4)
      
      #Again, if we want only a fixed set of nodes of the network to be returned in the final output:
      nodes_to_analyze <- c("A","B","C","D")
      CmmD_from_community_structures(nodelist=nodes_to_analyze, community_structures=structures,resolution_start=0,resolution_end=4, interval=0.5,distmethod= "hamming",threads = 4)

You will find more information in how to use the functions and the output returned from them in the documentation of the functions.

#### References: 
Núñez-Carpintero, I., Petrizzelli, M., Zinovyev, A., Cirillo, D. & Valencia, A. The multilayer community structure of medulloblastoma. iScience 24, (2021). DOI: 10.1016/j.isci.2021.102365

Didier, G., Valdeolivas, A., Baudot, A., 2018. Identifying communities from multiplex biological networks by randomized optimization of modularity. F1000Res 7. https://doi.org/10.12688/f1000research.15486.2

#### Contact: Iker Nuñez-Carpintero, Barcelona Supercomputing Center <iker.nunez@bsc.es>
