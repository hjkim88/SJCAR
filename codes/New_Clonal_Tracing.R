###
#   File name : New_Clonal_Tracing.R
#   Author    : Hyunjin Kim
#   Date      : Apr 5, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Define various clonotypes and evaluate the results
#
#   Instruction
#               1. Source("New_Clonal_Tracing.R")
#               2. Run the function "new_clonal_tracing" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_New_Clonal_Tracing.R/New_Clonal_Tracing.R")
#               > new_clonal_tracing(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                    outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/New_Tracing/")
###

new_clonal_tracing <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
                               outputDir="./results/New3/New_Tracing/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  
  
  
  
  
  
  
  
    
}
