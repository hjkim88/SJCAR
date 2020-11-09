###
#   File name : Read10xOutputsToSeurat.R
#   Author    : Hyunjin Kim
#   Date      : Oct 19, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : We would like to read 10x outputs and combine them to one Seurat object.
#               Other preprocessing will be performed and it will be saved as a RDATA file.
#
#   Instruction
#               1. Source("Read10xOutputsToSeurat.R")
#               2. Run the function "read10x_and_make_seuratobj" - specify the input directory and the output path
#               3. The results will be generated under the output path
#
#   Example
#               > source("The_directory_of_Read10xOutputsToSeurat.R/Read10xOutputsToSeurat.R")
#               > read10x_and_make_seuratobj(ten_x_dir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/JCC/JCC212_SJCAR19/JCC212_SJCAR19_C1AggregOct2020ns/filtered_feature_bc_matrix/",
#                                            outputPath="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/SJCAR19_Seurat_Obj_Oct2020.RDS")
###

read10x_and_make_seuratobj <- function(ten_x_dir="C:/Users/hkim8/SJ/SJCAR19/JCC212_SJCAR19_C1AggregOct2020ns/filtered_feature_bc_matrix/",
                                       outputPath="./data/SJCAR19_Seurat_Obj_Oct2020.RDS") {
  
  ### load libraries
  if(!require(dplyr, quietly = TRUE)) {
    install.packages("dplyr")
    require(dplyr, quietly = TRUE)
  }
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(patchwork, quietly = TRUE)) {
    install.packages("patchwork")
    require(patchwork, quietly = TRUE)
  }
  
  ### load the SJCAR19 dataset
  sjcar19.data <- Read10X(data.dir = ten_x_dir)
  
  ### create a Seurat object
  ### min.cells : include genes detected in at least this many cells
  ### min.features : include cells where at least this many features are detected
  SJCAR19_Oct2020_Seurat_Obj <- CreateSeuratObject(counts = sjcar19.data,
                                                   project = "SJCAR19_Oct2020",
                                                   min.cells = 3,
                                                   min.features = 200)
  
  ### MT percentage
  SJCAR19_Oct2020_Seurat_Obj[["percent.mt"]] <- PercentageFeatureSet(SJCAR19_Oct2020_Seurat_Obj, pattern = "^MT-")
  
  ### normalization
  SJCAR19_Oct2020_Seurat_Obj <- NormalizeData(SJCAR19_Oct2020_Seurat_Obj,
                                              normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  SJCAR19_Oct2020_Seurat_Obj <- FindVariableFeatures(SJCAR19_Oct2020_Seurat_Obj,
                                                     selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  SJCAR19_Oct2020_Seurat_Obj <- ScaleData(SJCAR19_Oct2020_Seurat_Obj)
  
  ### PCA
  SJCAR19_Oct2020_Seurat_Obj <- RunPCA(SJCAR19_Oct2020_Seurat_Obj,
                                       features = VariableFeatures(object = SJCAR19_Oct2020_Seurat_Obj))
  
  ### UMAP
  SJCAR19_Oct2020_Seurat_Obj <- RunUMAP(SJCAR19_Oct2020_Seurat_Obj, dims = 1:15)
  
  ### clustering
  SJCAR19_Oct2020_Seurat_Obj <- FindNeighbors(SJCAR19_Oct2020_Seurat_Obj, dims = 1:10)
  SJCAR19_Oct2020_Seurat_Obj <- FindClusters(SJCAR19_Oct2020_Seurat_Obj, resolution = 0.5)
  
  ### save the Seurat object as RDS file
  saveRDS(SJCAR19_Oct2020_Seurat_Obj, file = outputPath)
  
}
