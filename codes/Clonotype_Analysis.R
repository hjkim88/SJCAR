###
#   File name : Clonotype_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : May 5, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Now we have combined GEX and TCR info (GE, clonotypes, and lineages) of the SJCAR19 data
#               I would like to see if there is any clonotype that appeared frequently in the CAR+ cells
#               across PATIENTS, and also want to see gene expression patterns of the CAR+ cells that
#               have a lineage, and their pathways (biological functions).
#
#   Instruction
#               1. Source("Clonotype_Analysis.R")
#               2. Run the function "clonotype_analysis" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Clonotype_Analysis.R/Clonotype_Analysis.R")
#               > clonotype_analysis(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                                    outputDir="./results/PROTO/DEEP/")
###

clonotype_analysis <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                               outputDir="./results/PROTO/DEEP/") {
  
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
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### global clonotypes in the Seurat objects
  global_clonotypes <- colnames(Seurat_Obj@meta.data)[grep("global_clonotype", colnames(Seurat_Obj@meta.data), fixed = TRUE)]
  
  
  
}
