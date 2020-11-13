###
#   File name : IFITM3_PI3K_Pathway_Investigation.R
#   Author    : Hyunjin Kim
#   Date      : Nov 12, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Investigate associations of IFITM3 expression and PI3K signaling pathway.
#               1. (Combined) Correlation plots of IFITM3 - PI3K signaling genes based on expression
#               2. Heatmap of the IFITM3 & PI3K signaling genes with various annotation labels
#
#   Instruction
#               1. Source("IFITM3_PI3K_Pathway_Investigation.R")
#               2. Run the function "ifitm3_pi3k" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_IFITM3_PI3K_Pathway_Investigation.R/IFITM3_PI3K_Pathway_Investigation.R")
#               > ifitm3_pi3k(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                             outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results2/")
###

ifitm3_pi3k <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
                        outputDir="./results/New2/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(msigdbr, quietly = TRUE)) {
    install.packages("msigdbr")
    library(msigdbr, quietly = TRUE)
  }
  
  ### new output directory
  outputDir2 <- paste0(outputDir, "IFITM3_PI3K/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
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
  
  ### db preparation
  # MSIGDB - human
  m_df <- msigdbr(species = "Homo sapiens") 
  m_list <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  ### PI3K associated pathways
  pi3k_pathways <- m_list[names(m_list)[grep("PI3K", names(m_list))][c(9:13, 16:17, 21:23)]]
  
  ### 
  for(pathway in names(pi3k_pathways)) {
    
  }
  
  
}
