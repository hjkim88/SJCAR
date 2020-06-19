###
#   File name : Lineage_Scoring.R
#   Author    : Hyunjin Kim
#   Date      : Jun 17, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Give p-values to the GMP CAR+ clones that represent how significant a given CAR+ clone size is
#               at each time point in terms of GMP lasting capability.
#               * A matrix of rows: All GMP CAR+ Clones
#                           & cols: All the GMP+ time points + Total p-value
#               * One-tailed Fisher's exact test + Fisher's method
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
  
  ### get patient ids (dir names) from the result directory
  f <- list.dirs(outputDir, full.names = FALSE, recursive = FALSE)
  f <- f[grep("SJCAR19", f)]
  
  ### for each patient, perform the analysis
  for(px in f) {
    ### print progress
    writeLines(paste(px))
    
    ### get unique GMP CAR+ clones
    unique_gmp_carpos_clones <- unique(Seurat_Obj@meta.data$global_clonotype_ab_strict0[intersect(intersect(which(Seurat_Obj@meta.data$Px == px),
                                                                                                            which(Seurat_Obj@meta.data$Time == "GMP")),
                                                                                                  which(Seurat_Obj@meta.data$CAR == "CARpos"))])
    
    ### remove NA, "NA", "" from the clone list
    unique_gmp_carpos_clones <- unique_gmp_carpos_clones[which(!is.na(unique_gmp_carpos_clones))]
    unique_gmp_carpos_clones <- unique_gmp_carpos_clones[which(unique_gmp_carpos_clones != "NA")]
    unique_gmp_carpos_clones <- unique_gmp_carpos_clones[which(unique_gmp_carpos_clones != "")]
    
    ### get CAR+ time points
    time_points <- 
    
    if(length(unique_gmp_carpos_clones) > 0) {
      ### make an empty matrix
      result_mat <- matrix(0, nrow = length(unique_gmp_carpos_clones), ncol = )
    }
  }
  
  
  
  
  
  
  
}
