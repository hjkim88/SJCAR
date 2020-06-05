###
#   File name : CD4_CD8_Ratio.R
#   Author    : Hyunjin Kim
#   Date      : Jun 4, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : With CD4, CD8A, and CD8B expression, we determine whether a given cell is a CD4 cell
#               or a CD8 cell. Then we summarize how many CD4 & CD8 cells are in each library.
#
#   Instruction
#               1. Source("CD4_CD8_Ratio.R")
#               2. Run the function "get_cd4_cd8_ratio" - specify the input Seurat object and the output directory
#               3. The ratio summary table will be generated under the output directory
#
#   Example
#               > source("The_directory_of_CD4_CD8_Ratio.R/CD4_CD8_Ratio.R")
#               > get_cd4_cd8_ratio(Seurat_Obj=Seurat_Obj,
#                                   outputDir="./results/")
###

get_cd4_cd8_ratio <- function(Seurat_Obj=Seurat_Obj,
                              outputDir="./results/") {
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  
  ### get gene expressions of CD4/CD8A/CD8B
  cd4_exps <- Seurat_Obj@assays$RNA@counts["CD4",]
  cd8a_exps <- Seurat_Obj@assays$RNA@counts["CD8A",]
  cd8b_exps <- Seurat_Obj@assays$RNA@counts["CD8B",]
  
  ###
  
  
  
    
}
