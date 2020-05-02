# Project Description

This repository is for replicating parts of the results in the paper "Baron, Maayan, Adrian Veres, Samuel L. Wolock, Aubrey L. Faust, Renaud Gaujoux, Amedeo Vetere, Jennifer Hyoje Ryu, et al. 2016. “A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-Cell Population Structure.” Cell Systems 3 (4): 346–60.e4". This project is divided into three parts: data curator, programmer and analyst. It contains the bash code and R code that used to pre-processing the data and analyze the down-stream data.

# Contributors

Analyst part: Jinghan Huang (zx1an95) jh50@bu.edu

# Repository Contents

The analyst R code was used to identify marker genes for each cluster, label clusters as a cell type based on marker genes, visualize the clustered cells using UMAP, visualize the top marker genes per cluster and find novel marker genes. The functions were all from Seurat package. 
