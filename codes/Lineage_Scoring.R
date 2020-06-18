###
#   File name : Lineage_Scoring.R
#   Author    : Hyunjin Kim
#   Date      : Jun 17, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Give scores (or weights) to clones that indicate the probability of the clone will be
#               lasting after GMP infusion.
#
#   * If GMP prob == Wk2 prob is the null hypothesis,
#     When GMP prob >> Wk2 prob, the lineage represents a phenotype not good at lasting.
#     When GMP prob << Wk2 prob, the lineage represents a phenotype that is good at lasting.
#     We need to incorporate the pre-transduction population dynamics
#
#   Instruction
#               1. Source("Lineage_Scoring.R")
#               2. Run the function "give_lineage_score" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Lineage_Scoring.R/Lineage_Scoring.R")
#               > give_lineage_score(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                                    outputDir="./results/PROTO/")
###

give_lineage_score <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                               outputDir="./results/PROTO/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(DESeq2, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("DESeq2")
    require(DESeq2, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(metaseqR, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("metaseqR")
    require(metaseqR, quietly = TRUE)
  }
  if(!require(ggalluvial, quietly = TRUE)) {
    install.packages("ggalluvial")
    require(ggalluvial, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  if(!require(viridis, quietly = TRUE)) {
    install.packages("viridis")
    require(viridis, quietly = TRUE)
  }
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  
  
}
