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
#               2. Run the function "shared_tcrs" - specify the input directory and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Read10xOutputsToSeurat.R/Read10xOutputsToSeurat.R")
#               > shared_tcrs(ten_x_dir="C:/Users/hkim8/SJ/SJCAR19/JCC212_SJCAR19_C1AggregOct2020ns/filtered_feature_bc_matrix/",
#                             outputDir="./data/")
###

read10x_and_make_seuratobj <- function(ten_x_dir="C:/Users/hkim8/SJ/SJCAR19/JCC212_SJCAR19_C1AggregOct2020ns/filtered_feature_bc_matrix/",
                                       outputDir="./data/") {
  ### SJCAR19_Seurat_Obj_Oct2020.RDATA
  
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
  
  
  
  
  
  
  
  
}
