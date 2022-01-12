###
#   File name : CD_Revision.R
#   Author    : Hyunjin Kim
#   Date      : Jan 11, 2022
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. a) DE genes between cluster 3 & 8
#                  b) DE genes between GMP that ends up with cluster 3 vs GMP that ends up with cluster 8
#               2. a) DE genes between CAR > 0 VS CAR = 0
#                  b) DE genes between CAR > 0; HIGH vs CAR > 0; low
#               
#   Instruction
#               1. Source("CD_Revision.R")
#               2. Run the function "manuscript_revision" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_CD_Revision.R/CD_Revision.R")
#               > manuscript_revision(Seurat_RObj_path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds",
#                                     outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Manuscript/Revision/")
###

manuscript_revision <- function(Seurat_RObj_path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds",
                                outputDir="./results/New3/Manuscript/Revision/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### create outputDir
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### load Jeremy's object
  JCC_Seurat_Obj <- readRDS(file = "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds")
  
  ### check whether the orders are the same
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
  print(identical(names(Idents(object = JCC_Seurat_Obj)), rownames(JCC_Seurat_Obj@meta.data)))
  
  ###
  ### 1. a) DE genes between cluster 3 & 8
  #      b) DE genes between GMP that ends up with cluster 3 vs GMP that ends up with cluster 8
  ###
  
  ### a) DE genes between PI cluster 3 & 8
  
  ### check UMAP
  DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
          group.by = "AllSeuratClusters", label = TRUE,
          pt.size = 0.5, raster = TRUE)
  
  ### set new column for PI cluster 3 and PI cluster 8
  JCC_Seurat_Obj$PI_Cluster3_8 <- "Others"
  JCC_Seurat_Obj$PI_Cluster3_8[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                         which(JCC_Seurat_Obj$AllSeuratClusters == "3"))] <- "PI_Cluster3"
  JCC_Seurat_Obj$PI_Cluster3_8[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                         which(JCC_Seurat_Obj$AllSeuratClusters == "8"))] <- "PI_Cluster8"
  
  ### set idents
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$PI_Cluster3_8)
  
  ### cluster3 vs cluster8
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "PI_Cluster3",
                           ident.2 = "PI_Cluster8",
                           min.pct = 0.2,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/1_DE_PI_Clstr3_vs_8.xlsx"),
              sheetName = "1_DE_PI_Clstr3_vs_8", row.names = FALSE)
  
  
  ### b) DE genes between GMP that ends up with cluster 3 vs GMP that ends up with cluster 8
  
  ### PI subsisters in the cluster3 & 8
  cluster3_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                   which(JCC_Seurat_Obj$AllSeuratClusters == "3")),
                                                                         which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
  cluster8_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                   which(JCC_Seurat_Obj$AllSeuratClusters == "8")),
                                                                         which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
  cluster3_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster3_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
  cluster8_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster8_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
  
  ### remove NA clones
  cluster3_pi_subsister_clones <- cluster3_pi_subsister_clones[which(!is.na(cluster3_pi_subsister_clones))]
  cluster8_pi_subsister_clones <- cluster8_pi_subsister_clones[which(!is.na(cluster8_pi_subsister_clones))]
  
  ### get the GMP subsister indicies that end up in PI cluster 3 & 8
  GMP_Subsisters_PI_Cluster3_idx <- intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                              which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% cluster3_pi_subsister_clones))
  GMP_Subsisters_PI_Cluster8_idx <- intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                              which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% cluster8_pi_subsister_clones))
  
  ### set new column for GMP cells that end up in PI cluster 3 and in PI cluster 8
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8 <- "Others"
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8[GMP_Subsisters_PI_Cluster3_idx] <- "GMP_End_In_Cluster3"
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8[GMP_Subsisters_PI_Cluster8_idx] <- "GMP_End_In_Cluster8"
  
  ### set idents
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8)
  
  ### cluster3 vs cluster8
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "GMP_End_In_Cluster3",
                           ident.2 = "GMP_End_In_Cluster8",
                           min.pct = 0.2,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/1_DE_GMP_End_In_Clstr3_vs_8.xlsx"),
              sheetName = "1_DE_GMP_End_In_Clstr3_vs_8", row.names = FALSE)
  
  ### there can be GMP cells that will end up in both cluster3 and in cluster 8
  ### remove those and rerun the DE analysis
  
  ### set new column for GMP cells that end up in PI cluster 3 and in PI cluster 8
  ### but no duplicates this time
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8 <- "Others"
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8[setdiff(GMP_Subsisters_PI_Cluster3_idx,
                                                  GMP_Subsisters_PI_Cluster8_idx)] <- "GMP_End_In_Cluster3"
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8[setdiff(GMP_Subsisters_PI_Cluster8_idx,
                                                  GMP_Subsisters_PI_Cluster3_idx)] <- "GMP_End_In_Cluster8"
  
  ### set idents
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8)
  
  ### cluster3 vs cluster8
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "GMP_End_In_Cluster3",
                           ident.2 = "GMP_End_In_Cluster8",
                           min.pct = 0.2,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/1_DE_GMP_End_In_Clstr3_vs_8(2).xlsx"),
              sheetName = "1_DE_GMP_End_In_Clstr3_vs_8", row.names = FALSE)
  
  
  ###
  ### 2. a) DE genes between CAR > 0 VS CAR = 0
  #      b) DE genes between CAR > 0; HIGH vs CAR > 0; low
  ###
  ###    This should be done in 2 versions - 1. GMP only, 2. PI only
  ###    To see how CAR expression correlates with differentiation, exhaustion and/or apoptosis
  ###    Because the reviewer suspect that high CAR can drive those.
  ###
  
  ### The current JCC_Seurat_Obj is only with CAR+ cells
  ### need to load the full dataset
  ### make sure this is the right data - PB/BM all included with the same patient filtered
  ### why "Analyses_Figures_And_Tables.R" firstly use this file but later used JCC object?
  ### CAR? PB/BM? CD4/CD8? Maybe we wanna compare
  Seurat_Obj <- readRDS(file = "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj_Total2.rds")
  
  
  ### check UMAP based on CAR
  DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
          group.by = "CAR", label = TRUE,
          pt.size = 0.5, raster = TRUE)
  
  ### set new column for PI cluster 3 and PI cluster 8
  JCC_Seurat_Obj$PI_Cluster3_8 <- "Others"
  JCC_Seurat_Obj$PI_Cluster3_8[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                         which(JCC_Seurat_Obj$AllSeuratClusters == "3"))] <- "PI_Cluster3"
  JCC_Seurat_Obj$PI_Cluster3_8[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                         which(JCC_Seurat_Obj$AllSeuratClusters == "8"))] <- "PI_Cluster8"
  
  ### set idents
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$PI_Cluster3_8)
  
  ### cluster3 vs cluster8
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "PI_Cluster3",
                           ident.2 = "PI_Cluster8",
                           min.pct = 0.2,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  
  
  
  
  
  
}