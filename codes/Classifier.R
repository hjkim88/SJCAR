###
#   File name : Classifier.R
#   Author    : Hyunjin Kim
#   Date      : Jun 15, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Building some interesting classifier with the SJCAR19 data to classify such as:
#               1. GMP CAR+ cells that will last or not
#               2. Responder vs Non-responder
#               3. Cytokine levels
#
#   Instruction
#               1. Source("Classifier.R")
#               2. Run the function "classifier_run" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Classifier.R/Classifier.R")
#               > classifier_run(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                                outputDir="./results/PROTO/DEEP/")
###

classifier_run <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                           outputDir="./results/PROTO/DEEP/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(caret, quietly = TRUE)) {
    install.packages("caret")
    require(caret, quietly = TRUE)
  }
  if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
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
  
  ### set new result directory
  outputDir2 <- paste0(outputDir, "Classifier/")
  dir.create(path = outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  
  
    
}
