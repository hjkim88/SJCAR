###
#   File name : TCR_Clonotyping.R
#   Author    : Hyunjin Kim
#   Date      : Apr 3, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : After combining the GEX and TCR data, the "clonotype_id" represents
#               clonotypes in each time point. Therefore, I would like to do
#               a global clnotypying with all the time points.
#   Instruction
#               1. Source("TCR_Clonotyping.R")
#               2. Run the function "clonotyping" - specify the input file path and the output file path
#               3. The result Robj file will be generated in the output file path
#
#   Example
#               > source("The_directory_of_TCR_Clonotyping.R/TCR_Clonotyping.R")
#               > clonotyping(Seurat_RObj_path="./data/JCC212_Px5_TCR_combined.Robj",
#                             outRobjPath="./data/JCC212_Px5_TCR_clonotyped.Robj")
###

clonotyping <- function(Seurat_RObj_path="./data/JCC212_Px5_TCR_combined.Robj",
                        outRobjPath="./data/JCC212_Px5_TCR_clonotyped.Robj") {
  
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
  
  ### give the strict version of global clonotypes
  ### since the variables in each of the "cdr3_aa" are ordered,
  ### if the "cdr3_aa" are exactly the same, they are the same clonotype
  ### just ignore rows that are unique across all the rows - we are not interested in those
  Seurat_Obj@meta.data$global_clonotype_strict <- NA
  dup_idx <- intersect(which(duplicated(Seurat_Obj@meta.data$cdr3_aa)), which(!is.na(Seurat_Obj@meta.data$cdr3_aa)))
  ### the unique "cdr3_aa" sequences that exist more than 1 
  unique_dup_seqs <- unique(Seurat_Obj@meta.data$cdr3_aa[dup_idx])
  ### give clonotypes
  for(i in 1:length(unique_dup_seqs)) {
    idx <- which(Seurat_Obj@meta.data$cdr3_aa == unique_dup_seqs[i])
    Seurat_Obj@meta.data$global_clonotype_strict[idx] <- paste0("clonotype", i)
  }
  
  ### change the object name to the original one
  assign(obj_name, Seurat_Obj)
  rm(Seurat_Obj)
  gc()
  
  ### save the combined Seurat object
  save(list = c(obj_name), file = outRobjPath)
  
}
