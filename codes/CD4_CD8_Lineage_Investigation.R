###
#   File name : CD4_CD8_Lineage_Investigation.R
#   Author    : Hyunjin Kim
#   Date      : Nov 8, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Annotate CD4/CD8 to the cells and identify CD4/CD8 - associated lineages.
#
#   Instruction
#               1. Source("CD4_CD8_Lineage_Investigation.R")
#               2. Run the function "cd4_cd8_investigation" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_CD4_CD8_Lineage_Investigation.R/CD4_CD8_Lineage_Investigation.R")
#               > cd4_cd8_investigation(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/SJCAR19_Oct2020_Seurat_Obj3.RDS",
#                                       clonotype_lineage_info_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/SJCAR19_Clonotype_Lineages.RDS",
#                                       outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/")
###

cd4_cd8_investigation <- function(Seurat_RObj_path="./data/SJCAR19_Oct2020_Seurat_Obj3.RDS",
                                  clonotype_lineage_info_path="./data/SJCAR19_Clonotype_Lineages.RDS",
                                  outputDir="./results/New/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(SingleR, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("SingleR")
    require(SingleR, quietly = TRUE)
  }
  if(!require(SingleCellExperiment, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("SingleCellExperiment")
    require(SingleCellExperiment, quietly = TRUE)
  }
  
  ### create outputDir
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### load Seurat object
  Seurat_Obj <- readRDS(Seurat_RObj_path)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
  
  ### active assay = "RNA"
  Seurat_Obj@active.assay <- "RNA"
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  ### set time points
  total_time_points <- c("PreTrans", "PreTransB", "Wk-1", "Wk-1Run1", "Wk-1Run2", "Wk0", "GMP", "GMP-redo",
                         "Wk1", "Wk1b", "Wk2", "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo")
  
  ### download the hematopoietic cell population data
  NHD <- NovershternHematopoieticData()    
  
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$library)
  subset_Seurat_Obj <- subset(Seurat_Obj, idents = "SJCAR19-05_GMP_GMP")
  
  ### get gene expressions from the Seurat object
  target_mat <- as.SingleCellExperiment(subset_Seurat_Obj, assay = "RNA")
  
  ###
  NHD.main <- SingleR(test = target_mat, ref = list(NHD), labels = list(NHD$label.main))
  table(NHD.main$labels)
  
  
  
  NHD.fine <- SingleR(test = target_mat, ref = list(NHD), labels = list(NHD$label.fine ))
  table(C1_GMP_sub_NHD.fine$labels)
  
  
  
  
}
