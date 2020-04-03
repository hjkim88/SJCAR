###
#   File name : UMAP_with_GEX_TCR.R
#   Author    : Hyunjin Kim
#   Date      : Apr 3, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Draw a UMAP plot that shows how TCR clonality changes over time
#               and besides that, how GEX clusters change over time.
#
#   Instruction
#               1. Source("UMAP_with_GEX_TCR.R")
#               2. Run the function "umap_gex_tcr" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_UMAP_with_GEX_TCR.R/UMAP_with_GEX_TCR.R")
#               > umap_gex_tcr(Seurat_RObj_path="./data/JCC212_Px5_TCR_combined.Robj",
#                              outputDir="./results/")
###

umap_gex_tcr <- function(Seurat_RObj_path="./data/JCC212_Px5_TCR_combined.Robj",
                         outputDir="./results/") {
  
  ### load library
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  
  ### load the object
  load(Seurat_RObj_path)
  
  ###
  
  
  
}
