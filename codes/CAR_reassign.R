###
#   File name : CAR_reassign.R
#   Author    : Hyunjin Kim
#   Date      : Apr 11, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Currently, the SJCAR19 Seurat object's meta data has the CAR column
#               with 3 entries: GMP, CARpos, and CARneg. But I would like to make them
#               only have either CARpos or CARneg. Reassign the GMP rows based on the
#               CARcount column. Any rows that have CAR UMI more than 0 will be assigned
#               as CARpos. Otherwise, CARneg.
#
#   Instruction
#               1. Source("CAR_reassign.R")
#               2. Run the function "car_column_reassign" - specify the input file path and the output file path
#               3. The modified Seurat object will be generated in the output file path
#
#   Example
#               > source("The_directory_of_CAR_reassign.R/CAR_reassign.R")
#               > car_column_reassign(Seurat_RObj_path="./data/JCC212_Px5_TCR_clonotyped.Robj",
#                                     outRobjPath="./data/JCC212_Px5_TCR_clonotyped2.Robj")
###

car_column_reassign <- function(Seurat_RObj_path="./data/JCC212_Px5_TCR_clonotyped.Robj",
                                outRobjPath="./data/JCC212_Px5_TCR_clonotyped2.Robj") {
  
  ### load library
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### reassign
  Seurat_Obj@meta.data$CAR <- "CARneg"
  Seurat_Obj@meta.data$CAR[which(Seurat_Obj@meta.data$CARcount > 0)] <- "CARpos"
  
  ### change the object name to the original one
  assign(obj_name, Seurat_Obj)
  rm(Seurat_Obj)
  gc()
  
  ### save the combined Seurat object
  save(list = c(obj_name), file = outRobjPath)
  
}
