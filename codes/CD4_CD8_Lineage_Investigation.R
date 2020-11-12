###
#   File name : CD4_CD8_Lineage_Investigation.R
#   Author    : Hyunjin Kim
#   Date      : Nov 8, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Annotate CD4/CD8 to the cells and identify CD4/CD8 - associated lineages.
#
#   Instruction
#               1. Source("CD4_CD8_Lineage_Investigation.R")
#               2. Run the function "cd4_cd8_investigation" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_CD4_CD8_Lineage_Investigation.R/CD4_CD8_Lineage_Investigation.R")
#               > cd4_cd8_investigation(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                       clonotype_lineage_info_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Clonotype_Lineages.RDS",
#                                       outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results2/")
###

cd4_cd8_investigation <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
                                  clonotype_lineage_info_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Clonotype_Lineages.RDS",
                                  outputDir="./results/New2/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(SingleR, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("SingleR")
    require(SingleR, quietly = TRUE)
  }
  if(!require(SingleCellExperiment, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("SingleCellExperiment")
    require(SingleCellExperiment, quietly = TRUE)
  }
  
  ### create outputDir
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### load Seurat object
  Seurat_Obj <- readRDS(Seurat_RObj_path)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
  
  ### active assay = "RNA"
  Seurat_Obj@active.assay <- "RNA"
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  ### set time points
  total_time_points <- c("PreTrans", "PreTransB", "Wk-1", "Wk-1Run1", "Wk-1Run2", "Wk0", "GMP", "GMP-redo",
                         "Wk1", "Wk1b", "Wk2", "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo")
  
  ### CD4/CD8 by expression
  exp_quotient <- Seurat_Obj@assays$RNA@counts["CD4",] / Seurat_Obj@assays$RNA@counts["CD8A",]
  cd4_idx <- which(exp_quotient > 1)
  cd8_idx <- which(exp_quotient < 1)
  both_idx <- which(exp_quotient == 1)
  Seurat_Obj@meta.data$CD4_CD8_by_Exp <- NA
  Seurat_Obj@meta.data$CD4_CD8_by_Exp[cd4_idx] <- "CD4"
  Seurat_Obj@meta.data$CD4_CD8_by_Exp[cd8_idx] <- "CD8"
  Seurat_Obj@meta.data$CD4_CD8_by_Exp[both_idx] <- "BOTH"
  rm(exp_quotient)
  gc()
  
  ### download the hematopoietic cell population data
  NHD <- NovershternHematopoieticData()
  
  ### set idents with the libary
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$library)
  
  ### for each library add annotations of CD4/CD8 using singleR
  Seurat_Obj@meta.data$CD4_CD8_by_NHD_main <- NA
  Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine <- NA
  for(lib in unique(Seurat_Obj@meta.data$library)) {
    ### print progress
    writeLines(paste(lib))
    
    ### get gene expressions from the subset of the Seurat object
    target_mat <- as.SingleCellExperiment(subset(Seurat_Obj, idents = lib), assay = "RNA")
    
    ### get CD4/CD8 annotation using the main labels
    NHD.main <- SingleR(test = target_mat, ref = list(NHD), labels = list(NHD$label.main))$labels
    
    ### get CD4/CD8 annotation using the fine labels
    NHD.fine <- SingleR(test = target_mat, ref = list(NHD), labels = list(NHD$label.fine))$labels
    
    ### add annotions to the meta data
    Seurat_Obj@meta.data[colnames(target_mat),"CD4_CD8_by_NHD_main"] <- NHD.main
    Seurat_Obj@meta.data[colnames(target_mat),"CD4_CD8_by_NHD_fine"] <- NHD.fine
    
    gc()
  }
  
  ### CD4/CD8 Consensus
  Seurat_Obj@meta.data$CD4_CD8_by_Consensus <- NA
  for(i in 1:nrow(Seurat_Obj@meta.data)) {
    cd4_cnt <- 0
    cd8_cnt <- 0
    if(Seurat_Obj@meta.data$CD4_CD8_by_Exp == "CD4") {
      cd4_cnt <- cd4_cnt + 1
    } else if (Seurat_Obj@meta.data$CD4_CD8_by_Exp == "CD8") {
      cd8_cnt <- cd8_cnt + 1
    }
    if(grepl("CD4", Seurat_Obj@meta.data$CD4_CD8_by_NHD_main[i], fixed = TRUE)) {
      cd4_cnt <- cd4_cnt + 1
    } else if(grepl("CD8", Seurat_Obj@meta.data$CD4_CD8_by_NHD_main[i], fixed = TRUE)) {
      cd8_cnt <- cd8_cnt + 1
    }
    if(grepl("CD4", Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine[i], fixed = TRUE)) {
      cd4_cnt <- cd4_cnt + 1
    } else if(grepl("CD8", Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine[i], fixed = TRUE)) {
      cd8_cnt <- cd8_cnt + 1
    }
    
    ### decide with the consensus
    if(cd4_cnt > cd8_cnt) {
      Seurat_Obj@meta.data$CD4_CD8_by_Consensus[i] <- "CD4"
    } else if(cd4_cnt < cd8_cnt) {
      Seurat_Obj@meta.data$CD4_CD8_by_Consensus[i] <- "CD8"
    } else if((cd4_cnt != 0) && (cd4_cnt == cd8_cnt)) {
      Seurat_Obj@meta.data$CD4_CD8_by_Consensus[i] <- "BOTH"
    }
    
    ### print progress
    if(i %% 10000 == 0) {
      writeLines(paste(i, "/", nrow(Seurat_Obj@meta.data)))
    }
  }
  
  ### load clonotype lineages info
  SJCAR19_Clonotype_Frequency <- readRDS(clonotype_lineage_info_path)
  
  
}
