###
#   File name : Epitope_Spreading.R
#   Author    : Hyunjin Kim
#   Date      : Dec 9, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Analyses about epitope spreading in SJCAR19 data.
#
#   Instruction
#               1. Source("Epitope_Spreading.R")
#               2. Run the function "epitope_spreading_investigation" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Epitope_Spreading.R/Epitope_Spreading.R")
#               > epitope_spreading_investigation(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                                 clonotype_lineage_info_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Clonotype_Lineages.RDS",
#                                                 outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results2/")
###

epitope_spreading_investigation <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
                                            clonotype_lineage_info_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Clonotype_Lineages.RDS",
                                            outputDir="./results/New2/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  
  
  
  
  
  
}
