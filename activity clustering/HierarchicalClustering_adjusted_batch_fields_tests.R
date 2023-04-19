#Author: Jan Hirtz 2020 - 2023

# this code makes use of the Dynamic tree cut algorithm as pubslihed in Langfelder et al. (2008): Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. Bioinformatics 24:719-720.10.1093/bioinformatics/btm563

#install.packages("cluster")
#install.packages("dynamicTreeCut")
#install.packages("igraph")
#install.packages("ggplot2")

library(cluster)
library(dynamicTreeCut)
library(igraph)
library(ggplot2)

username <- Sys.getenv("USERNAME")
exchange_directory <- paste0("fill_in_path")

##hierarchiccal clustering

  
for (n in 0:3) {

f <- paste0("fill_in_path/R_cell_distance_matrix_field_",n,".txt")
if (file.exists(f)) {
  z <- read.table(paste0(exchange_directory,"R_cell_distance_matrix_field_",n,".txt"),  row.names=1, fill = TRUE, stringsAsFactors=FALSE)
  z <- as.matrix(z)
  
  distCo <- as.dist(z, diag = TRUE)
  hclCo <- hclust(distCo, method = "average")
  
  plot(hclCo, labels = FALSE)
  
  clust_hclCo_dyn <- dynamicTreeCut::cutreeDynamic(
    hclCo, distM = z, method = "hybrid",
    minClusterSize = 2, deepSplit = 0.5, pamStage = FALSE
  )
  table(clust_hclCo_dyn)
  
  write.table(clust_hclCo_dyn, file=paste0(exchange_directory, "export_R_cell_clusters_field_",n), sep=",") }

f <- paste0("fill_in_path/R_sound_distance_matrix_field_",n,".txt")
if (file.exists(f)) {
  zs <- read.table(paste0(exchange_directory, "R_sound_distance_matrix_field_",n,".txt"),  row.names=1, fill = TRUE, stringsAsFactors=FALSE)
  zs <- as.matrix(zs)
  
  distCos <- as.dist(zs, diag = TRUE)
  hclCos <- hclust(distCos, method = "average")
  
  plot(hclCos, labels = FALSE)
  
  clust_hclCo_dyns <- dynamicTreeCut::cutreeDynamic(
    hclCos, distM = zs, method = "hybrid",
    minClusterSize = 2, deepSplit = 2.5, pamStage = FALSE
  )
  table(clust_hclCo_dyns)
  
  write.table(clust_hclCo_dyns, file=paste0(exchange_directory, "export_R_sound_clusters_field_",n), sep=",")} }
