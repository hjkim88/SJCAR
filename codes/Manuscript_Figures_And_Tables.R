###
#   File name : Manuscript_Figures_And_Tables.R
#   Author    : Hyunjin Kim
#   Date      : Apr 15, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. Alluvial plots of lineages - CAR+ a) all patients into one b) separated by patient
#               2. Number of lineages & clonotypes table CD4+, CD8+, CAR+, CAR-, per patient
#                  (alpha only, beta only, one from each, strict)
#               3. CAR+ % for each library - proportional bar graph & table - CAR+, CAR-, CD4+, CD8+
#               4. Clone size between CAR+ lineages vs non-lineage CAR+ after infusion
#               5. a) PCA/UMAP plot of CAR+s over time (coloring based on time)
#                  b) Clustering + UMAP
#                  c) Cluster 0,1,4,8 vs others
#                  d) After infusion CAR+ subsister vs non-subsisters to find CAR is differentially expressed
#                  e) Where are the subsisters that close to non-subsisters located in the after infusion CAR+ UMAP?
#                  f) Where are the outliers in the after infusion CAR+ UMAP located in the GMP PCA?
#                  g) Where are the CD4+ CAR+ subsister cells lie in the UMAP?
#                  h) Pseudotime analysis on PCA
#                  i) PCA/UMAP plot of lineages with size=1 vs lineages with size > 1 (coloring differently)
#                  j) Find all markers based on the clustering
#                  k) UMAP plot of GMP subsisters
#               6. Visualize some interesting genes on UMAP
#               7. CAR+ (>0, >1, >2, etc.) numbers for each patient between "From Sorting" and "From scRNA-Seq"
#               8. a) If sampling from GMP and sampling from an after infusion time point, how many matches do we see?
#                  Using GMP - clone size as a background to estimate a selection factor
#                  And show distribution shift
#                  b) Use the % of CD4/CD8 subsisters vs Total CD4/CD8 in GMP to estimate how many CD8 subsisters were infused
#                  c) Use the % of Cluster01 vs Total CD4/CD8 in GMP to estimate how many subsister-like cells were infused
#                  d) Correlation between [subsister infusion amount (per kg)] and [PeakCAR]
#               9. The number of cells and clones in each cluster
#               10. Gini index of the CAR+ lineages over time
#               11. DE result - Violin & Dot plot
#               12. BM cells in UMAP & their info + compare repertoires and GEx profiles- are BM and PB CAR+ cells similar?
#               13. Clustering on GMP UMAP and find subsister clusters + calculate dose (related to 5(k))
#               14. What are the DEGs between cells from the late time (after six weeks) points vs. early time points?
#               15. Feature plots and DE list with specific genes from Tan paper
#               16. Trajectory analysis of differentiation of CD8 CAR+ (start with all time points and perhaps move to only post-infusion depending on how it looks)
#               17. Comaprison to a reference data set of well characterized t cell differentiation gene sets
#               18. CAR+ cells with high CAR expression vs CAR+ cells with low CAR expression (DE, pathway, & subsister difference)
#               19. Look at all GMP (CAR+ & CAR-) and all CAR+ GMP cells whether they are two separated clusters
#               20. Tay's request to look at some genes of Px11
#               21. Multiple regression to estimate peakcar using both dose level & tumor burden
#               22. Classifier - reperform with the new data and also with the subsister cluster info (not only based on the subsisters all the cluster cells)
#               23. PCA/UMAP/MNN comparison with the original GMP CARpos cells
#               24. Mapping GMP clusters to PI clusters
#               25. 06/14/2021 - Re-analyze everything with the PB-Px filtered data with different CD4/CD8 definition
#               26. New Task - 3D plane mapping between GMP & PI clusters
#               27. 06/21/2021 - Pseudotime (Slingshot & Monocle2) analysis on Jeremy's object
#               28. Comparison of DE genes between "GMP CAR+ S vs NS" & "After infusion CAR+ S vs NS"
#
#   Instruction
#               1. Source("Manuscript_Figures_And_Tables.R")
#               2. Run the function "manuscript_prep" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Manuscript_Figures_And_Tables.R/Manuscript_Figures_And_Tables.R")
#               > manuscript_prep(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                 gmp_carpos_cd8_de_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/GMP_CARpos_CD8_Persisters_vs_NonPersisters.xlsx",
#                                 px_result_dir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/",
#                                 outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Manuscript/")
###

manuscript_prep <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj_Total2.RDS",
                            gmp_carpos_cd8_de_path="./results/New3/Manuscript/New/GMP_CARpos_CD8_Persisters_vs_NonPersisters.xlsx",
                            px_result_dir="./results/New3/PB_ONLY/",
                            outputDir="./results/New3/Manuscript/New/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
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
  if(!require(grid, quietly = TRUE)) {
    install.packages("grid")
    require(grid, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(slingshot, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("slingshot")
    require(slingshot, quietly = TRUE)
  }
  if(!require(tradeSeq, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("tradeSeq")
    require(tradeSeq, quietly = TRUE)
  }
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    require(scales, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(DescTools, quietly = TRUE)) {
    install.packages("DescTools")
    require(DescTools, quietly = TRUE)
  }
  if(!require(ggrepel, quietly = TRUE)) {
    install.packages("ggrepel")
    require(ggrepel, quietly = TRUE)
  }
  if(!require(monocle, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("monocle")
    require(monocle, quietly = TRUE)
  }
  if(!require(gage, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("gage")
    require(gage, quietly = TRUE)
  }
  if(!require(msigdbr, quietly = TRUE)) {
    install.packages("msigdbr")
    library(msigdbr, quietly = TRUE)
  }
  if(!require(ggthemes, quietly = TRUE)) {
    install.packages("ggthemes")
    require(ggthemes, quietly = TRUE)
  }
  if(!require(DelayedMatrixStats, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("DelayedMatrixStats")
    require(DelayedMatrixStats, quietly = TRUE)
  }
  if(!require(scran, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("scran")
    require(scran, quietly = TRUE)
  }
  if(!require(remotes, quietly = TRUE)) {
    install.packages("remotes")
    require(remotes, quietly = TRUE)
  }
  if(!require(SeuratWrappers, quietly = TRUE)) {
    remotes::install_github('satijalab/seurat-wrappers')
    require(SeuratWrappers, quietly = TRUE)
  }
  if(!require(e1071, quietly = TRUE)) {
    install.packages("e1071")
    require(e1071, quietly = TRUE)
  }
  if(!require(hydroGOF, quietly = TRUE)) {
    install.packages("hydroGOF")
    require(hydroGOF, quietly = TRUE)
  }
  if(!require(caret, quietly = TRUE)) {
    install.packages("caret")
    require(caret, quietly = TRUE)
  }
  if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
  }
  if(!require(scDblFinder, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("scDblFinder")
    require(scDblFinder, quietly = TRUE)
  }
  if(!require(plot3D, quietly = TRUE)) {
    install.packages("plot3D")
    require(plot3D, quietly = TRUE)
  }
  if(!require(plot3Drgl, quietly = TRUE)) {
    install.packages("plot3Drgl")
    require(plot3Drgl, quietly = TRUE)
  }
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    require(gplots, quietly = TRUE)
  }
  if(!require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    require(RColorBrewer, quietly = TRUE)
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
  
  ### combine some separated time points into one
  Seurat_Obj@meta.data$time2 <- Seurat_Obj@meta.data$time
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "GMP-redo")] <- "GMP"
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "PreTransB")] <- "PreTrans"
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "Wk1b")] <- "Wk1"
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "Wk-1Run1")] <- "Wk-1"
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "Wk-1Run2")] <- "Wk-1"
  
  ### set time points
  total_time_points <- c("PreTrans", "Wk-1", "Wk0", "GMP", "Wk1", "Wk2", 
                         "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo")
  gmp_after_time_points <- c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", 
                             "Wk6", "Wk8", "3mo", "6mo", "9mo")
  
  ### save the BM cells and only retain PB cells
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$tissue)
  BM_Seurat_obj <- subset(Seurat_Obj, idents = c("GMP", "BM"))
  Seurat_Obj <- subset(Seurat_Obj, idents = c("GMP", "PB"))
  
  ### clonotype types
  clonotype_names <- c("Alpha", "Beta", "One_From_Each", "Strict")
  clonotype_types <- c("clonotype_id_by_patient_alpha", "clonotype_id_by_patient_beta",
                       "clonotype_id_by_patient_one_alpha_beta", "clonotype_id_by_patient")
  names(clonotype_names) <- clonotype_types
  names(clonotype_types) <- clonotype_names
  
  ### for each clonotype type to get lineages
  total_lineages <- vector("list", length = length(clonotype_names))
  names(total_lineages) <- clonotype_names
  total_gmp_lineages <- vector("list", length = length(clonotype_names))
  names(total_gmp_lineages) <- clonotype_names
  for(type in clonotype_names) {
    total_lineages[[type]] <- NULL
    total_gmp_lineages[[type]] <- NULL
    ### run for each patient
    for(patient in unique(Seurat_Obj@meta.data$px)) {
      
      ### print progress
      writeLines(paste(patient))
      
      ### load the file
      target_file <- read.xlsx2(file = paste0(px_result_dir, patient, "/car_clonotype_frequency_over_time_", patient, ".xlsx"),
                                sheetName = paste0("CARpos_Clonotype_Frequency_", substr(type, 1, 4)), stringsAsFactors = FALSE, check.names = FALSE,
                                row.names = 1)
      
      ### numerize the table
      for(i in 1:ncol(target_file)) {
        target_file[,i] <- as.numeric(target_file[,i])
      }
      
      ### combine some redundant time points to one
      if(length(which(colnames(target_file) == "GMP-redo")) > 0) {
        target_file[,"GMP"] <- target_file[,"GMP"] + target_file[,"GMP-redo"]
        target_file <- target_file[,-which(colnames(target_file) == "GMP-redo")]
      }
      if(length(which(colnames(target_file) == "PreTransB")) > 0) {
        if(length(which(colnames(target_file) == "PreTrans")) > 0) {
          target_file[,"PreTrans"] <- target_file[,"PreTrans"] + target_file[,"PreTransB"]
          target_file <- target_file[,-which(colnames(target_file) == "PreTransB")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "PreTransB")] <- "PreTrans"
        }
      }
      if(length(which(colnames(target_file) == "Wk1b")) > 0) {
        if(length(which(colnames(target_file) == "Wk1")) > 0) {
          target_file[,"Wk1"] <- target_file[,"Wk1"] + target_file[,"Wk1b"]
          target_file <- target_file[,-which(colnames(target_file) == "Wk1b")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "Wk1b")] <- "Wk1"
        }
      }
      if(length(which(colnames(target_file) == "Wk-1Run1")) > 0) {
        if(length(which(colnames(target_file) == "Wk-1")) > 0) {
          target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run1"]
          target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run1")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "Wk-1Run1")] <- "Wk-1"
        }
      }
      if(length(which(colnames(target_file) == "Wk-1Run2")) > 0) {
        if(length(which(colnames(target_file) == "Wk-1")) > 0) {
          target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run2"]
          target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run2")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "Wk-1Run2")] <- "Wk-1"
        }
      }
      
      ### remove all zero time points
      time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
      
      ### get time points except the Total
      time_points <- setdiff(time_points, c("Total"))
      
      ### if there are at least two time points
      if(length(time_points) > 1) {
        ###  get lineages
        lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
          return(length(which(x > 0)) > 1)  
        })),time_points]
        lineages <- rownames(lineage_table)
        
        ### get gmp lineages
        gmp_lineage_table <- lineage_table[which(lineage_table$GMP > 0),]
        gmp_lineages <- rownames(gmp_lineage_table)
        
        ### combine lineages
        if(length(lineages) > 0) {
          ### combine
          if(is.null(total_lineages[[type]])) {
            total_lineages[[type]] <- lineages
          } else {
            total_lineages[[type]] <- c(total_lineages[[type]], lineages)
          }
        }
        
        ### combine gmp lineages
        if(length(gmp_lineages) > 0) {
          ### combine
          if(is.null(total_gmp_lineages[[type]])) {
            total_gmp_lineages[[type]] <- gmp_lineages
          } else {
            total_gmp_lineages[[type]] <- c(total_gmp_lineages[[type]], gmp_lineages)
          }
        }
      }
      
      gc()
    }
  }
  
  ### some pre-calculated indices
  carpos_idx <- which(Seurat_Obj@meta.data$CAR == "CARpos")
  gmp_carpos_idx <- intersect(carpos_idx,
                              which(Seurat_Obj@meta.data$time2 == "GMP"))
  cd4_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4"),
                              carpos_idx)
  cd8_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"),
                              carpos_idx)
  gmp_cd8_carpos_idx <- intersect(cd8_carpos_idx,
                                  which(Seurat_Obj@meta.data$time2 == "GMP"))
  no_tcr_idx <- which(is.na(Seurat_Obj@meta.data$cdr3_aa))
  lineages_idx <- vector("list", length = length(total_lineages))
  names(lineages_idx) <- names(total_lineages)
  for(type in names(lineages_idx)) {
    lineages_idx[[type]] <- which(Seurat_Obj@meta.data[,clonotype_types[type]] %in% total_lineages[[type]])
  }
  gmp_lineages_idx <- vector("list", length = length(total_gmp_lineages))
  names(gmp_lineages_idx) <- names(total_gmp_lineages)
  for(type in names(gmp_lineages_idx)) {
    gmp_lineages_idx[[type]] <- which(Seurat_Obj@meta.data[,clonotype_types[type]] %in% total_gmp_lineages[[type]])
  }
  
  #
  ### define lineages
  #
  
  ### all the CAR+ lineages appeared in at least any two time points
  ### including all the time points
  ### "YES": All CAR+ Persister Cells
  ### "NO": CAR+ cells that are not in the lineages (non-persisters)
  ### NA: CAR- cells and CAR+ cells that don't have TCR info
  Seurat_Obj@meta.data$ALL_CARpos_Persister <- NA
  Seurat_Obj@meta.data$ALL_CARpos_Persister[intersect(lineages_idx[["One_From_Each"]],
                                                      carpos_idx)] <- "YES"
  Seurat_Obj@meta.data$ALL_CARpos_Persister[setdiff(setdiff(carpos_idx,
                                                            lineages_idx[["One_From_Each"]]),
                                                    no_tcr_idx)] <- "NO"
  
  ### CAR+ lineages that have GMP time point (and at least one more post-infusion time point)
  ### Cells in the GMP time point ONLY
  ### "YES": GMP CAR+ Persister cells
  ### "NO": GMP CAR+ cells that are non-persisters
  ### NA: Others including cells that don't have TCR info
  Seurat_Obj@meta.data$GMP_CARpos_Persister <- NA
  Seurat_Obj@meta.data$GMP_CARpos_Persister[intersect(gmp_lineages_idx[["One_From_Each"]],
                                                      gmp_carpos_idx)] <- "YES"
  Seurat_Obj@meta.data$GMP_CARpos_Persister[setdiff(setdiff(gmp_carpos_idx,
                                                            gmp_lineages_idx[["One_From_Each"]]),
                                                    no_tcr_idx)] <- "NO"
  
  ### CD8 cells in CAR+ lineages that have GMP time point (and at least one more post-infusion time point)
  ### CD8 Cells in the GMP time point ONLY
  ### "YES": GMP CAR+ CD8 Persister cells
  ### "NO": GMP CAR+ CD8 cells that are non-persisters
  ### NA: Others including cells that don't have TCR info
  Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister <- NA
  Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister[intersect(gmp_lineages_idx[["One_From_Each"]],
                                                          gmp_cd8_carpos_idx)] <- "YES"
  Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister[setdiff(setdiff(gmp_cd8_carpos_idx,
                                                                gmp_lineages_idx[["One_From_Each"]]),
                                                        no_tcr_idx)] <- "NO"
  
  ### lineages
  all_gmp_persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")])
  
  ### CAR+ lineages that have GMP time point (and at least one more post-infusion time point)
  ### including all the time points
  ### "YES": CAR+ Persister cells in GMP lineages (all time)
  ### "NO": CAR+ cells that are not in the GMP lineages
  ### NA: Others including cells that don't have TCR info
  Seurat_Obj@meta.data$ALL_GMP_CARpos_Persister <- NA
  Seurat_Obj@meta.data$ALL_GMP_CARpos_Persister[intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% all_gmp_persister_clones),
                                                          carpos_idx)] <- "YES"
  Seurat_Obj@meta.data$ALL_GMP_CARpos_Persister[setdiff(setdiff(carpos_idx,
                                                                which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% all_gmp_persister_clones)),
                                                        no_tcr_idx)] <- "NO"
  
  ### filter out some insufficient patients
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$px)
  Seurat_Obj <- subset(Seurat_Obj, idents = c("SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                                              "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                              "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### remove some trash meta.data columns
  Seurat_Obj@meta.data$RNA_snn_res.0.5 <- NULL
  Seurat_Obj@meta.data$seurat_clusters <- NULL
  
  ### theme that draws dotted lines for each y-axis ticks
  ### this function is from "immunarch" package
  theme_cleveland2 <- function(rotate = TRUE) {
    if (rotate) {
      theme(
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(
          colour = "grey70",
          linetype = "dashed"
        )
      )
    }
    else {
      theme(
        panel.grid.major.x = element_line(
          colour = "grey70",
          linetype = "dashed"
        ), panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )
    }
  }
  
  ### a function for color brewer
  cell_pal <- function(cell_vars, pal_fun) {
    if (is.numeric(cell_vars)) {
      pal <- pal_fun(100)
      return(pal[cut(cell_vars, breaks = 100)])
    } else {
      categories <- sort(unique(cell_vars))
      pal <- setNames(pal_fun(length(categories)), categories)
      return(pal[cell_vars])
    }
  }
  
  #' @title Plot Slingshot output
  #' @name plot-SlingshotDataSet
  #' @aliases plot-SlingshotDataSet plot,SlingshotDataSet,ANY-method
  #'
  #' @description Tools for visualizing lineages inferred by \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{"curves"}, or \code{"both"} (by partial matching),
  #'   see Details for more.
  #' @param linInd integer, an index indicating which lineages should be plotted
  #'   (default is to plot all lineages). If \code{col} is a vector, it will be
  #'   subsetted by \code{linInd}.
  #' @param show.constraints logical, whether or not the user-specified initial
  #'   and terminal clusters should be specially denoted by green and red dots,
  #'   respectively.
  #' @param add logical, indicates whether the output should be added to an
  #'   existing plot.
  #' @param dims numeric, which dimensions to plot (default is \code{1:2}).
  #' @param asp numeric, the y/x aspect ratio, see \code{\link{plot.window}}.
  #' @param cex numeric, amount by which points should be magnified, see
  #'   \code{\link{par}}.
  #' @param lwd numeric, the line width, see \code{\link{par}}.
  #' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
  #' @param ... additional parameters to be passed to \code{\link{lines}}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' plot(sds, type = 'b')
  #'
  #' # add to existing plot
  #' plot(rd, col = 'grey50')
  #' lines(sds, lwd = 3)
  #'
  #' @import graphics
  #' @export
  setMethod(
    f = "plot",
    signature = signature(x = "SlingshotDataSet"),
    definition = function(x, type = NULL,
                          linInd = NULL,
                          show.constraints = FALSE,
                          constraints.col = NULL,
                          add = FALSE,
                          dims = seq_len(2),
                          asp = 1,
                          cex = 2,
                          lwd = 2,
                          col = 1,
                          ...) {
      col <- rep(col, length(slingLineages(x)))
      curves <- FALSE
      lineages <- FALSE
      if(is.null(type)){
        if(length(slingCurves(x)) > 0){
          type <- 'curves'
        }else if(length(slingLineages(x)) > 0){
          type <- 'lineages'
        }else{
          stop('No lineages or curves detected.')
        }
      }else{
        type <- c('curves','lineages','both')[pmatch(type,
                                                     c('curves','lineages','both'))]
        if(is.na(type)){
          stop('Unrecognized type argument.')
        }
      }
      
      if(type %in% c('lineages','both')){
        lineages <- TRUE
      }
      if(type %in% c('curves','both')){
        curves <- TRUE
      }
      
      if(lineages & (length(slingLineages(x))==0)){
        stop('No lineages detected.')
      }
      if(curves & (length(slingCurves(x))==0)){
        stop('No curves detected.')
      }
      
      if(is.null(linInd)){
        linInd <- seq_along(slingLineages(x))
      }else{
        linInd <- as.integer(linInd)
        if(!all(linInd %in% seq_along(slingLineages(x)))){
          if(any(linInd %in% seq_along(slingLineages(x)))){
            linInd.removed <-
              linInd[! linInd %in% seq_along(slingLineages(x))]
            linInd <-
              linInd[linInd %in% seq_along(slingLineages(x))]
            message('Unrecognized lineage indices (linInd): ',
                    paste(linInd.removed, collapse = ", "))
          }else{
            stop('None of the provided lineage indices',
                 ' (linInd) were found.')
          }
        }
      }
      
      if(lineages){
        X <- reducedDim(x)
        clusterLabels <- slingClusterLabels(x)
        connectivity <- slingMST(x)
        clusters <- rownames(connectivity)
        nclus <- nrow(connectivity)
        centers <- t(vapply(clusters,function(clID){
          w <- clusterLabels[,clID]
          return(apply(X, 2, weighted.mean, w = w))
        }, rep(0,ncol(X))))
        rownames(centers) <- clusters
        X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
        clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                       drop = FALSE]
        linC <- slingParams(x)
        clus2include <- unique(unlist(slingLineages(x)[linInd]))
      }
      
      if(!add){
        xs <- NULL
        ys <- NULL
        if(lineages){
          xs <- c(xs, centers[,dims[1]])
          ys <- c(ys, centers[,dims[2]])
        }
        if(curves){
          npoints <- nrow(slingCurves(x)[[1]]$s)
          xs <- c(xs, as.numeric(vapply(slingCurves(x),
                                        function(c){ c$s[,dims[1]] }, rep(0,npoints))))
          ys <- c(ys, as.numeric(vapply(slingCurves(x),
                                        function(c){ c$s[,dims[2]] }, rep(0,npoints))))
        }
        plot(x = NULL, y = NULL, asp = asp,
             xlim = range(xs), ylim = range(ys),
             xlab = colnames(reducedDim(x))[dims[1]],
             ylab = colnames(reducedDim(x))[dims[2]])
      }
      
      if(lineages){
        for(i in seq_len(nclus-1)){
          for(j in seq(i+1,nclus)){
            if(connectivity[i,j]==1 &
               all(clusters[c(i,j)] %in% clus2include)){
              lines(centers[c(i,j), dims],
                    lwd = lwd, col = col[1], ...)
            }
          }
        }
        points(centers[clusters %in% clus2include, dims],
               cex = cex+1, pch = 16, col = col[1])
        if(show.constraints && !is.null(constraints.col)){
          for(const in names(constraints.col)) {
            points(centers[clusters %in% const, dims,
                           drop=FALSE], cex = cex-0.5,
                   col = constraints.col[const], pch = 16)
            # text(x = centers[clusters %in% const, dims[1]]+0,
            #      y = centers[clusters %in% const, dims[2]]+2,
            #      labels = const,
            #      font = 2,
            #      cex = cex-0.5,
            #      col = "black")
          }
        }
      }
      if(curves){
        for(ii in seq_along(slingCurves(x))[linInd]){
          c <- slingCurves(x)[[ii]]
          lines(c$s[c$ord, dims], lwd = lwd, col = col[ii], ...)
        }
      }
      invisible(NULL)
    }
  )
  
  #' @rdname plot-SlingshotDataSet
  #' @importFrom graphics lines
  #' @export
  lines.SlingshotDataSet <- function(x,
                                     type = NULL,
                                     dims = seq_len(2),
                                     ...) {
    plot(x, type = type, add = TRUE, dims = dims, ...)
    invisible(NULL)
  }
  
  #' @title Pairs plot of Slingshot output
  #' @name pairs-SlingshotDataSet
  #'
  #' @description A tool for quickly visualizing lineages inferred by
  #'   \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
  #'   Details for more.
  #' @param show.constraints logical, whether or not the user-specified initial
  #'   and terminal clusters should be specially denoted by green and red dots,
  #'   respectively.
  #' @param col character, color vector for points.
  #' @param pch integer or character specifying the plotting symbol, see
  #'   \code{\link{par}}.
  #' @param cex numeric, amount by which points should be magnified, see
  #'   \code{\link{par}}.
  #' @param lwd numeric, the line width, see \code{\link{par}}.
  #' @param ... additional parameters for \code{plot} or \code{axis}, see
  #'   \code{\link{pairs}}.
  #' @param labels character, the names of the variables, see \code{\link{pairs}}.
  #' @param horInd see \code{\link{pairs}}.
  #' @param verInd see \code{\link{pairs}}.
  #' @param lower.panel see \code{\link{pairs}}.
  #' @param upper.panel see \code{\link{pairs}}.
  #' @param diag.panel see \code{\link{pairs}}.
  #' @param text.panel see \code{\link{pairs}}.
  #' @param label.pos see \code{\link{pairs}}.
  #' @param line.main see \code{\link{pairs}}.
  #' @param cex.labels see \code{\link{pairs}}.
  #' @param font.labels see \code{\link{pairs}}.
  #' @param row1attop see \code{\link{pairs}}.
  #' @param gap see \code{\link{pairs}}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' pairs(sds, type = 'curves')
  #'
  #' @export
  pairs.SlingshotDataSet <-
    function (x, type = NULL, show.constraints = FALSE, col = NULL,
              constraints.col = NULL,
              pch = 16, cex=1, lwd=2, ...,
              labels, horInd = seq_len(nc), verInd = seq_len(nc),
              lower.panel = FALSE, upper.panel = TRUE,
              diag.panel = NULL, text.panel = textPanel,
              label.pos = 0.5 + has.diag/3, line.main = 3,
              cex.labels = NULL, font.labels = 1,
              row1attop = TRUE, gap = 1,
              xlim=NULL, ylim=NULL) {
      #####
      lp.sling <- lower.panel
      up.sling <- upper.panel
      panel <- points
      if(!up.sling){
        upper.panel <- NULL
      }else{
        upper.panel <- panel
      }
      if(!lower.panel){
        lower.panel <- NULL
      }else{
        lower.panel <- panel
      }
      log = ""
      sds <- x
      x <- reducedDim(sds)
      curves <- FALSE
      lineages <- FALSE
      if(is.null(type)){
        if(length(slingCurves(sds)) > 0){
          type <- 'curves'
        }else if(length(slingLineages(sds)) > 0){
          type <- 'lineages'
        }else{
          stop('No lineages or curves detected.')
        }
      }else{
        type <- c('curves','lineages','both')[pmatch(type,
                                                     c('curves','lineages',
                                                       'both'))]
        if(is.na(type)){
          stop('Unrecognized type argument.')
        }
      }
      if(type %in% c('lineages','both')){
        lineages <- TRUE
      }
      if(type %in% c('curves','both')){
        curves <- TRUE
      }
      if(lineages & (length(slingLineages(sds))==0)){
        stop('No lineages detected.')
      }
      if(curves & (length(slingCurves(sds))==0)){
        stop('No curves detected.')
      }
      if(lineages){
        forest <- slingAdjacency(sds)
        clusters <- rownames(forest)
        nclus <- nrow(forest)
        centers <- t(vapply(clusters,function(clID){
          w <- slingClusterLabels(sds)[,clID]
          return(apply(x, 2, weighted.mean, w = w))
        }, rep(0,ncol(reducedDim(sds)))))
        rownames(centers) <- clusters
        linC <- slingParams(sds)
      }
      range.max <- max(apply(x,2,function(xi){
        r <- range(xi, na.rm = TRUE)
        return(abs(r[2] - r[1]))
      }))
      plot.ranges <- apply(x,2,function(xi){
        mid <- (max(xi,na.rm = TRUE) + min(xi,na.rm = TRUE))/2
        return(c(mid - range.max/2, mid + range.max/2))
      })
      if(is.null(col)){
        if(requireNamespace("RColorBrewer", quietly = TRUE)) {
          cc <- c(RColorBrewer::brewer.pal(9, "Set1")[-c(1,3,6)],
                  RColorBrewer::brewer.pal(7, "Set2")[-2],
                  RColorBrewer::brewer.pal(6, "Dark2")[-5],
                  RColorBrewer::brewer.pal(8, "Set3")[-c(1,2)])
        } else {
          cc <- seq_len(100)
        }
        col <- cc[apply(slingClusterLabels(sds),1,which.max)]
      }
      #####
      if(doText <- missing(text.panel) || is.function(text.panel))
        textPanel <-
        function(x = 0.5, y = 0.5, txt, cex, font)
          text(x, y, txt, cex = cex, font = font)
      
      localAxis <- function(side, x, y, xpd, bg, col=NULL, lwd=NULL, main,
                            oma, ...) {
        ## Explicitly ignore any color argument passed in as
        ## it was most likely meant for the data points and
        ## not for the axis.
        xpd <- NA
        if(side %% 2L == 1L && xl[j]) xpd <- FALSE
        if(side %% 2L == 0L && yl[i]) xpd <- FALSE
        if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
      }
      
      localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
      localLowerPanel <- function(..., main, oma, font.main, cex.main)
        lower.panel(...)
      localUpperPanel <- function(..., main, oma, font.main, cex.main)
        upper.panel(...)
      localDiagPanel <- function(..., main, oma, font.main, cex.main)
        diag.panel(...)
      
      dots <- list(...); nmdots <- names(dots)
      if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for(i in seq_along(names(x))) {
          if(is.factor(x[[i]]) || is.logical(x[[i]]))
            x[[i]] <- as.numeric(x[[i]])
          if(!is.numeric(unclass(x[[i]])))
            stop("non-numeric argument to 'pairs'")
        }
      } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")
      panel <- match.fun(panel)
      if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel <- match.fun(lower.panel)
      if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel <- match.fun(upper.panel)
      if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))
        diag.panel <- match.fun( diag.panel)
      
      if(row1attop) {
        tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
        tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
      }
      
      nc <- ncol(x)
      if (nc < 2L) stop("only one column in the argument to 'pairs'")
      if(!all(horInd >= 1L & horInd <= nc))
        stop("invalid argument 'horInd'")
      if(!all(verInd >= 1L & verInd <= nc))
        stop("invalid argument 'verInd'")
      if(doText) {
        if (missing(labels)) {
          labels <- colnames(x)
          if (is.null(labels)) labels <- paste("var", 1L:nc)
        }
        else if(is.null(labels)) doText <- FALSE
      }
      oma <- if("oma" %in% nmdots) dots$oma
      main <- if("main" %in% nmdots) dots$main
      if (is.null(oma))
        oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)
      opar <- par(mfrow = c(length(horInd), length(verInd)),
                  mar = rep.int(gap/2, 4), oma = oma)
      on.exit(par(opar))
      dev.hold(); on.exit(dev.flush(), add = TRUE)
      
      xl <- yl <- logical(nc)
      if (is.numeric(log)) xl[log] <- yl[log] <- TRUE
      else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}
      for (i in if(row1attop) verInd else rev(verInd))
        for (j in horInd) {
          l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
          
          if(is.null(xlim) & !is.null(ylim))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = plot.ranges[,j], ylim=ylim)
          else if(!is.null(xlim) & is.null(ylim))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = xlim, ylim = plot.ranges[,i])
          else if(!is.null(xlim) & !is.null(ylim))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = xlim, ylim=ylim)
          else
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = plot.ranges[,j], ylim = plot.ranges[,i])
          
          if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
            box()
            if(i == 1  && (!(j %% 2L) || !has.upper || !has.lower ))
              localAxis(1L + 2L*row1attop, x[, j], x[, i], ...)
            if(i == nc && (  j %% 2L  || !has.upper || !has.lower ))
              localAxis(3L - 2L*row1attop, x[, j], x[, i], ...)
            if(j == 1  && (!(i %% 2L) || !has.upper || !has.lower ))
              localAxis(2L, x[, j], x[, i], ...)
            if(j == nc && (  i %% 2L  || !has.upper || !has.lower ))
              localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if(i == j) {
              if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
              if (doText) {
                par(usr = c(0, 1, 0, 1))
                if(is.null(cex.labels)) {
                  l.wid <- strwidth(labels, "user")
                  cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                }
                xlp <- if(xl[i]) 10^0.5 else 0.5
                ylp <- if(yl[j]) 10^label.pos else label.pos
                text.panel(xlp, ylp, labels[i],
                           cex = cex.labels, font = font.labels)
              }
            } else if(i < j){
              if(up.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        seg.col <- 1
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col, ...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex=2*cex)
                  if(show.constraints && is.null(constraints.col)){
                    if(any(linC$start.given)){
                      st.ind <- clusters %in%
                        linC$start.clus[linC$start.given]
                      points(centers[st.ind,j],
                             centers[st.ind,i], cex = cex,
                             col = 'green3',
                             pch = pch)
                    }
                    if(any(linC$end.given)){
                      en.ind <- clusters %in%
                        linC$end.clus[linC$end.given]
                      points(centers[en.ind,j],
                             centers[en.ind,i], cex = cex,
                             col = 'red2', pch = pch)
                    }
                  } else if(show.constraints && !is.null(constraints.col)){
                    for(const in names(constraints.col)) {
                      points(centers[clusters %in% const, j, drop=FALSE],
                             centers[clusters %in% const, i, drop=FALSE],
                             cex = cex, pch = 16,
                             col = constraints.col[const])
                    }
                  }
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)], lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            else{
              if(lp.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        if(clusters[ii] %in%
                           linC$start.clus |
                           clusters[jj] %in%
                           linC$start.clus){
                          seg.col <- 'green3'
                        }else if(clusters[ii] %in%
                                 linC$end.clus[
                                   linC$end.given] |
                                 clusters[jj] %in%
                                 linC$end.clus[
                                   linC$end.given]){
                          seg.col <- 'red2'
                        }else{
                          seg.col <- 1
                        }
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col,...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex = 2*cex)
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)],lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            if (any(par("mfg") != mfg))
              stop("the 'panel' function made a new plot")
          } else par(new = FALSE)
          
        }
      if (!is.null(main)) {
        font.main <- if("font.main" %in% nmdots){
          dots$font.main
        }else par("font.main")
        cex.main <- if("cex.main" %in%
                       nmdots) dots$cex.main else par("cex.main")
        mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main,
              font = font.main)
      }
      invisible(NULL)
    }
  
  #' @name plot3d-SlingshotDataSet
  #' @title Plot Slingshot output in 3D
  #'
  #' @description Tools for visualizing lineages inferred by \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
  #'   Details for more.
  #' @param linInd integer, an index indicating which lineages should be plotted
  #'   (default is to plot all lineages). If \code{col} is a vector, it will be
  #'   subsetted by \code{linInd}.
  #' @param add logical, indicates whether the output should be added to an
  #'   existing plot.
  #' @param dims numeric, which dimensions to plot (default is \code{1:3}).
  #' @param aspect either a logical indicating whether to adjust the aspect ratio
  #'   or a new ratio, see \code{\link[rgl:plot3d]{plot3d}}.
  #' @param size numeric, size of points for MST (default is \code{10}), see
  #'   \code{\link[rgl:plot3d]{plot3d}}.
  #' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
  #' @param ... additional parameters to be passed to \code{lines3d}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' \dontrun{
  #' library(rgl)
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' rd <- cbind(rd, rnorm(nrow(rd)))
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' plot3d(sds, type = 'b')
  #'
  #' # add to existing plot
  #' plot3d(rd, col = 'grey50', aspect = 'iso')
  #' plot3d(sds, lwd = 3, add = TRUE)
  #' }
  # #' @importFrom rgl plot3d
  #' @export
  plot3d.SlingshotDataSet <- function(x,
                                      type = NULL,
                                      linInd = NULL,
                                      add = FALSE,
                                      dims = seq_len(3),
                                      aspect = 'iso',
                                      size = 10,
                                      col = 1,
                                      col2 = NULL,
                                      ...){
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("Package 'rgl' is required for 3D plotting.",
           call. = FALSE)
    }
    col <- rep(col, length(slingLineages(x)))
    n <- nrow(reducedDim(x))
    curves <- FALSE
    lineages <- FALSE
    if(is.null(type)){
      if(length(slingCurves(x)) > 0){
        type <- 'curves'
      }else if(length(slingLineages(x)) > 0){
        type <- 'lineages'
      }else{
        stop('No lineages or curves detected.')
      }
    }else{
      type <- c('curves','lineages','both')[pmatch(type,c('curves','lineages',
                                                          'both'))]
      if(is.na(type)){
        stop('Unrecognized type argument.')
      }
    }
    
    if(type %in% c('lineages','both')){
      lineages <- TRUE
    }
    if(type %in% c('curves','both')){
      curves <- TRUE
    }
    
    if(lineages & (length(slingLineages(x))==0)){
      stop('No lineages detected.')
    }
    if(curves & (length(slingCurves(x))==0)){
      stop('No curves detected.')
    }
    
    if(is.null(linInd)){
      linInd <- seq_along(slingLineages(x))
    }else{
      linInd <- as.integer(linInd)
      if(!all(linInd %in% seq_along(slingLineages(x)))){
        if(any(linInd %in% seq_along(slingLineages(x)))){
          linInd.removed <-
            linInd[! linInd %in% seq_along(slingLineages(x))]
          linInd <-
            linInd[linInd %in% seq_along(slingLineages(x))]
          message('Unrecognized lineage indices (linInd): ',
                  paste(linInd.removed, collapse = ", "))
        }else{
          stop('None of the provided lineage indices',
               ' (linInd) were found.')
        }
      }
    }
    
    if(lineages){
      X <- reducedDim(x)
      clusterLabels <- slingClusterLabels(x)
      connectivity <- slingAdjacency(x)
      clusters <- rownames(connectivity)
      nclus <- nrow(connectivity)
      centers <- t(vapply(clusters,function(clID){
        w <- clusterLabels[,clID]
        return(apply(X, 2, weighted.mean, w = w))
      }, rep(0,ncol(X))))
      rownames(centers) <- clusters
      X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
      clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                     drop = FALSE]
      clus2include <- unique(unlist(slingLineages(x)[linInd]))
    }
    
    if(!add){
      xs <- NULL
      ys <- NULL
      zs <- NULL
      if(lineages){
        xs <- c(xs, centers[,dims[1]])
        ys <- c(ys, centers[,dims[2]])
        zs <- c(zs, centers[,dims[3]])
      }
      if(curves){
        npoints <- nrow(slingCurves(x)[[1]]$s)
        xs <- c(xs, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[1]] }, rep(0,npoints))))
        ys <- c(ys, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[2]] }, rep(0,npoints))))
        zs <- c(zs, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[3]] }, rep(0,npoints))))
      }
      rgl::plot3d(x = NULL, y = NULL, z = NULL, aspect = aspect,
                  xlim = range(xs), ylim = range(ys), zlim = range(zs),
                  xlab = colnames(reducedDim(x))[dims[1]],
                  ylab = colnames(reducedDim(x))[dims[2]],
                  zlab = colnames(reducedDim(x))[dims[3]])
    }
    
    if(lineages){
      for(i in seq_len(nclus-1)){
        for(j in seq(i+1,nclus)){
          if(connectivity[i,j]==1 &
             all(clusters[c(i,j)] %in% clus2include)){
            rgl::lines3d(x = centers[c(i,j),dims[1]],
                         y = centers[c(i,j),dims[2]],
                         z = centers[c(i,j),dims[3]],
                         col = col[1], ...)
          }
        }
      }
      rgl::points3d(centers[clusters %in% clus2include, dims],
                    size = size/2, col = col2[clusters[clusters %in% clus2include]])
      rgl::points3d(centers[clusters %in% clus2include, dims],
                    size = size, col = col[1])
    }
    if(curves){
      for(ii in seq_along(slingCurves(x))[linInd]){
        c <- slingCurves(x)[[ii]]
        rgl::lines3d(c$s[c$ord,dims], col = col[ii], ...)
      }
    }
    invisible(NULL)
  }
  
  ### the plot3d.SlingshotDataSet of the slingshot package is incomplete and too simple,
  ### so, i'm implementing a 3d plot function myself
  slingshot_3d_lineages <- function(slingshot_obj, color, title,
                                    print=FALSE, outputDir=NULL,
                                    width=1200, height=800) {
    
    ### load libraries
    if(!require(Seurat, quietly = TRUE)) {
      install.packages("Seurat")
      require(Seurat, quietly = TRUE)
    }
    if(!require(slingshot, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("slingshot")
      require(slingshot, quietly = TRUE)
    }
    if(!require(rgl, quietly = TRUE)) {
      install.packages("rgl")
      require(rgl, quietly = TRUE)
    }
    if(!require(rmarkdown, quietly = TRUE)) {
      install.packages("rmarkdown")
      require(rmarkdown, quietly = TRUE)
    }
    
    #
    ### 3D Slingshot
    #
    
    ### draw 3D PCA
    par3d(windowRect = c(50, 50, width+50, height+50))
    plot3d.SlingshotDataSet(slingshot_obj, dims = 1:3, col = "black", col2 = color, type = "lineages", add = TRUE)
    plot3d(slingshot_obj@reducedDim, col = apply(slingshot_obj@clusterLabels, 1, function(x) color[names(x)[which(x == 1)]]),
           size = 5, alpha = 0.5, aspect = FALSE, add = TRUE)
    axes3d(edges=c("x+-", "y+-", "z++"), lwd = 2,
           labels=TRUE, tick = FALSE, nticks = 3, box = TRUE, expand = 1.05)
    mtext3d(text = expression(bold("PC1")), edge="x+-", line = -2, at = min(slingshot_obj@reducedDim[,1]), pos = NA)
    mtext3d(text = expression(bold("PC2")), edge="y+-", line = -2, at = min(slingshot_obj@reducedDim[,2]), pos = NA)
    mtext3d(text = expression(bold("PC3")), edge="z++", line = -2, at = max(slingshot_obj@reducedDim[,3]), pos = NA)
    decorate3d(xlim = NULL, ylim = NULL, zlim = NULL, 
               xlab = "", ylab = "", zlab = "", 
               box = FALSE, axes = FALSE, main = title, sub = NULL,
               top = TRUE, aspect = FALSE, expand = 1.05)
    legend3d("topright", legend = names(color), title = "Clusters",
             col = color, pch = 19, cex=2)
    if(print) {
      writeWebGL(dir=outputDir, filename = paste0(outputDir, title, ".html"),
                 width=width, height = height)
    }
    
  }
  
  # ******************************************************************************************
  # Pathway Analysis with clusterProfiler package
  # Input: geneList     = a vector of gene Entrez IDs for pathway analysis [numeric or character]
  #        org          = organism that will be used in the analysis ["human" or "mouse"]
  #                       should be either "human" or "mouse"
  #        database     = pathway analysis database (KEGG or GO) ["KEGG" or "GO"]
  #        title        = title of the pathway figure [character]
  #        pv_threshold = pathway analysis p-value threshold (not DE analysis threshold) [numeric]
  #        displayNum   = the number of pathways that will be displayed [numeric]
  #                       (If there are many significant pathways show the few top pathways)
  #        imgPrint     = print a plot of pathway analysis [TRUE/FALSE]
  #        dir          = file directory path of the output pathway figure [character]
  #
  # Output: Pathway analysis results in figure - using KEGG and GO pathways
  #         The x-axis represents the number of DE genes in the pathway
  #         The y-axis represents pathway names
  #         The color of a bar indicates adjusted p-value from the pathway analysis
  #         For Pathview Result, all colored genes are found DE genes in the pathway,
  #         and the color indicates log2(fold change) of the DE gene from DE analysis
  # ******************************************************************************************
  pathwayAnalysis_CP <- function(geneList,
                                 org,
                                 database,
                                 title="Pathway_Results",
                                 pv_threshold=0.05,
                                 displayNum=Inf,
                                 imgPrint=TRUE,
                                 dir="./") {
    
    ### load library
    if(!require(clusterProfiler, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("clusterProfiler")
      require(clusterProfiler, quietly = TRUE)
    }
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    
    ### collect gene list (Entrez IDs)
    geneList <- geneList[which(!is.na(geneList))]
    
    if(!is.null(geneList)) {
      ### make an empty list
      p <- list()
      
      if(database == "KEGG") {
        ### KEGG Pathway
        kegg_enrich <- enrichKEGG(gene = geneList, organism = org, pvalueCutoff = pv_threshold)
        
        if(is.null(kegg_enrich)) {
          writeLines("KEGG Result does not exist")
          return(NULL)
        } else {
          kegg_enrich@result <- kegg_enrich@result[which(kegg_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(kegg_enrich@result) <= displayNum)) {
              result <- kegg_enrich@result
              description <- kegg_enrich@result$Description
            } else {
              result <- kegg_enrich@result[1:displayNum,]
              description <- kegg_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(kegg_enrich) > 0) {
              p[[1]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 30) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title)) +
                theme(axis.text = element_text(size = 30))
              
              ggsave(file = paste0(dir, "KEGG_", title, "_CB.png"), plot = p[[1]], width = 35, height = 10, dpi = 350)
            } else {
              writeLines("KEGG Result does not exist")
            }
          }
          
          return(kegg_enrich@result)
        }
      } else if(database == "GO") {
        ### GO Pathway
        if(org == "human") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else if(org == "mouse") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else {
          go_enrich <- NULL
          writeLines(paste("Unknown org variable:", org))
        }
        
        if(is.null(go_enrich)) {
          writeLines("GO Result does not exist")
          return(NULL)
        } else {
          go_enrich@result <- go_enrich@result[which(go_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(go_enrich@result) <= displayNum)) {
              result <- go_enrich@result
              description <- go_enrich@result$Description
            } else {
              result <- go_enrich@result[1:displayNum,]
              description <- go_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(go_enrich) > 0) {
              p[[2]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 30) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("GO ", title)) +
                theme(axis.text = element_text(size = 30))
              
              ggsave(file = paste0(dir, "GO_", title, "_CB.png"), plot = p[[2]], width = 35, height = 10, dpi = 350)
            } else {
              writeLines("GO Result does not exist")
            }
          }
          
          return(go_enrich@result)
        }
      } else {
        stop("database prameter should be \"GO\" or \"KEGG\"")
      }
    } else {
      writeLines("geneList = NULL")
    }
  }
  
  #
  ### GSEA with the important genes of the PC1
  #
  
  #'****************************************************************************************
  #' Gene Set Enrichment Analysis function
  #' 
  #' It receives gene list (character vector) and signature profiles (named numeric vector)
  #' as inputs, performs GSEA and returns a table of GSEA result table and draws
  #' a GSEA plot. It is basically a statistical significance test to check how the
  #' given given genes are biased on the signature profiles.
  #' 
  #' Whether there are multiple gene sets or multiple signatures,
  #' multiple testing (FDR computation) is performed.
  #' But if the input gene set and the input signature are both lists with multiple
  #' items (The length of the two are both more than 1) then we return an error message.
  #' 
  #' The plot file names will be determined by names(gene_list) or names(signature)
  #' If length(gene_list) > 1, then names(gene_list) will be used and
  #' if length(signature) > 1, then names(signature) will be used as file names.
  #' If there is no list names, then file names will be "GSEA_Plot_i.png".
  #' Here, i indicates that the plot is from the i-th row of the GSEA result table.
  #' 
  #' * Some plot drawing codes were from Rtoolbox/R/ReplotGSEA.R written by Thomas Kuilman. 
  #'****************************************************************************************
  #' @title	run_gsea
  #' 
  #' @param gene_list   A list of character vectors containing gene names to be tested
  #' @param signature   A list of named numeric vectors of signature values for GSEA. The gene_list
  #'                    should be included in the names(signature)
  #' @param printPlot   If TRUE, it also generates GSEA plot of the results
  #'                    (Default = FALSE)
  #' @param fdr_cutoff  When printing GSEA plots, print them with the FDR < fdr_cutoff only
  #'                    (Default = 0.05)
  #' @param heatmap_color_type  when 'relative', the heatmap of the GSEA colors the bottom half of the
  #'                            absolute range of the signature as blue and the upper half as red
  #'                            when 'absolute', the heatmap of GSEA colors the negative signature as blue
  #'                            and the positives as red
  #' @param printPath   When printing GSEA plots, print them in the designated path
  #'                    (Default = "./")
  #' @param width       The width of the plot file
  #'                    (Default = 2000)
  #' @param height      The height of the plot file
  #'                    (Default = 1200)
  #' @param res         The resolution of the plot file
  #'                    (Default = 130)
  #' 
  #' @return 	          It tests bias of the "gene_list" on the "signature" range and
  #'                    returns a table including p-values and FDRs (adjusted p-values)
  #'                    If fdr_cutoff == TRUE, it also generates a GSEA plot with the result
  #' 
  run_gsea <- function(gene_list,
                       signature,
                       printPlot = FALSE,
                       fdr_cutoff = 0.05,
                       heatmap_color_type = c("relative", "absolute"),
                       width = 2000,
                       height = 1200,
                       res = 350,
                       printPath = "./") {
    
    ### load required libraries
    if(!require("fgsea", quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("fgsea")
      require("fgsea", quietly = TRUE)
    }
    if(!require("limma", quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("limma")
      require("limma", quietly = TRUE)
    } 
    if(!require("checkmate", quietly = TRUE)) {
      install.packages("checkmate")
      require("checkmate", quietly = TRUE)
    }
    
    ### argument checking
    assertList(gene_list)
    assertList(signature)
    assertLogical(printPlot)
    assertNumeric(fdr_cutoff)
    assertIntegerish(width)
    assertIntegerish(height)
    assertIntegerish(res)
    assertString(printPath)
    if(length(gene_list) > 1 && length(signature) > 1) {
      stop("ERROR: \"gene_list\" and \"signature\" cannot be both \"list\"")
    }
    
    ### set random seed
    set.seed(1234)
    
    ### run GSEA
    ### if there are more than one signatures
    if(length(signature) > 1) {
      ### combine GSEA results of every signature inputs
      for(i in 1:length(signature)) {
        temp <- data.frame(fgseaMultilevel(pathways = gene_list, stats = signature[[i]]))
        if(i == 1) {
          gsea_result <- temp
        } else {
          gsea_result <- rbind(gsea_result, temp)
        }
      }
      
      ### compute FDRs
      corrected_gsea_result <- gsea_result[order(gsea_result$pval),]
      corrected_gsea_result$padj <- p.adjust(corrected_gsea_result$pval, method = "BH")
      gsea_result <- corrected_gsea_result[rownames(gsea_result),]
    } else {
      ### if there are more than one gene sets
      gsea_result <- data.frame(fgseaMultilevel(pathways = gene_list, stats = signature[[1]]))
    }
    
    ### print GSEA plot
    sIdx <- which(gsea_result$padj < fdr_cutoff)
    if(printPlot && length(sIdx) > 0) {
      for(i in sIdx) {
        ### get required values ready
        if(length(signature) > 1) {
          geneset <- gene_list[[1]]
          stats <- signature[[i]]
          stats <- stats[order(-stats)]
          fileName <- names(signature)[i]
        } else {
          geneset <- gene_list[[i]]
          stats <- signature[[1]]
          stats <- stats[order(-stats)]
          fileName <- names(gene_list)[i]
        }
        if(is.null(fileName)) {
          fileName <- paste0("GSEA_Plot_", i)
        }
        stats <- stats[!is.na(stats)]
        gsea.hit.indices <- which(names(stats) %in% geneset)
        es.temp <- calcGseaStat(stats, gsea.hit.indices, returnAllExtremes = TRUE)
        if(es.temp$res >= 0) {
          gsea.es.profile <- es.temp$tops
        } else {
          gsea.es.profile <- es.temp$bottoms
        }
        enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
        metric.range <- c(min(stats), max(stats))
        gsea.p.value <- round(gsea_result$pval[i] ,5)
        gsea.fdr <- round(gsea_result$padj[i] ,5)
        gsea.enrichment.score <- round(gsea_result$ES[i], 5)
        gsea.normalized.enrichment.score <- round(gsea_result$NES[i], 5)
        
        ### print GSEA result plot
        png(paste0(printPath, fileName, ".png"), width = width, height = height, res = res)
        
        ### set layout
        layout.show(layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2)))
        
        ### draw the GSEA plot
        par(mar = c(0, 5, 2, 2))
        plot(c(1, gsea.hit.indices, length(stats)),
             c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
             xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
             ylim = enrichment.score.range,
             main = list(fileName, font = 1, cex = 1),
             panel.first = {
               abline(h = seq(round(enrichment.score.range[1], digits = 1),
                              enrichment.score.range[2], 0.1),
                      col = "gray95", lty = 2)
               abline(h = 0, col = "gray50", lty = 2)
             }
        )
        
        ### add informative text to the GSEA plot
        plot.coordinates <- par("usr")
        if(es.temp$res < 0) {
          text(length(stats) * 0.01, plot.coordinates[3] * 0.98,
               paste("P-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
                     gsea.enrichment.score, "\nNormalized ES:",
                     gsea.normalized.enrichment.score), adj = c(0, 0))
        } else {
          text(length(stats) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
               paste("P-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
                     gsea.enrichment.score, "\nNormalized ES:",
                     gsea.normalized.enrichment.score), adj = c(1, 1))
        }
        
        ### draw hit indices
        par(mar = c(0, 5, 0, 2))
        plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
             ylab = "", xlim = c(1, length(stats)))
        abline(v = gsea.hit.indices, lwd = 0.75)
        
        ### create color palette for the heatmap
        par(mar = c(0, 5, 0, 2))
        if(heatmap_color_type[1] == "relative") {
          rank.colors <- stats - metric.range[1]
          rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
          rank.colors <- ceiling(rank.colors * 511 + 1)
          rank.colors <- colorRampPalette(c("blue", "white", "red"))(512)[rank.colors]
        } else {
          rank.colors1 <- stats[which(stats >= 0)]
          rank.colors1 <- rank.colors1 - min(rank.colors1)
          rank.colors1 <- rank.colors1 / (max(rank.colors1) - min(rank.colors1))
          rank.colors1 <- ceiling(rank.colors1 * 255 + 1)
          rank.colors1 <- colorRampPalette(c("white", "red"))(256)[rank.colors1]
          rank.colors2 <- stats[which(stats < 0)]
          rank.colors2 <- rank.colors2 - min(rank.colors2)
          rank.colors2 <- rank.colors2 / (max(rank.colors2) - min(rank.colors2))
          rank.colors2 <- ceiling(rank.colors2 * 255 + 1)
          rank.colors2 <- colorRampPalette(c("blue", "white"))(256)[rank.colors2]
          rank.colors <- c(rank.colors1, rank.colors2)
        }
        
        ### draw the heatmap
        rank.colors <- rle(rank.colors)
        barplot(matrix(rank.colors$lengths), col = rank.colors$values,
                border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(stats)))
        box()
        text(length(stats) / 2, 0.7,
             labels = "Signature")
        text(length(stats) * 0.01, 0.7, "Largest", adj = c(0, NA))
        text(length(stats) * 0.99, 0.7, "Smallest", adj = c(1, NA))
        
        ### draw signature values
        par(mar = c(5, 5, 0, 2))
        rank.metric <- rle(round(stats, digits = 2))
        plot(stats, type = "n", xaxs = "i",
             xlab = "Rank in ordered gene list", xlim = c(0, length(stats)),
             ylim = metric.range, yaxs = "i",
             ylab = "Signature values",
             panel.first = abline(h = seq(metric.range[1] / 2,
                                          metric.range[2] - metric.range[1] / 4,
                                          metric.range[2] / 2), col = "gray95", lty = 2))
        
        barplot(rank.metric$values, col = "lightgrey", lwd = 0.1,
                xlim = c(0, length(stats)), ylim = c(-1, 1),
                width = rank.metric$lengths, border = NA,
                space = 0, add = TRUE, xaxt = "n")
        box()
        
        ### print out the file
        dev.off()
      }
    }
    
    return(gsea_result)
    
  }
  
  
  #
  ### 1. Alluvial plots of lineages - CAR+ a) all patients into one b) separated by patient
  #
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/1/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### run for each patient
  total_plot_df <- NULL
  p <- vector("list", length= length(unique(Seurat_Obj@meta.data$px)))
  names(p) <- unique(Seurat_Obj@meta.data$px)
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    
    ### print progress
    writeLines(paste(patient))
    
    ### load the file
    target_file <- read.xlsx2(file = paste0(px_result_dir, patient, "/car_clonotype_frequency_over_time_", patient, ".xlsx"),
                              sheetName = paste0("CARpos_Clonotype_Frequency_One_"), stringsAsFactors = FALSE, check.names = FALSE,
                              row.names = 1)
    
    ### numerize the table
    for(i in 1:ncol(target_file)) {
      target_file[,i] <- as.numeric(target_file[,i])
    }
    
    ### combine some redundant time points to one
    if(length(which(colnames(target_file) == "GMP-redo")) > 0) {
      target_file[,"GMP"] <- target_file[,"GMP"] + target_file[,"GMP-redo"]
      target_file <- target_file[,-which(colnames(target_file) == "GMP-redo")]
    }
    if(length(which(colnames(target_file) == "PreTransB")) > 0) {
      if(length(which(colnames(target_file) == "PreTrans")) > 0) {
        target_file[,"PreTrans"] <- target_file[,"PreTrans"] + target_file[,"PreTransB"]
        target_file <- target_file[,-which(colnames(target_file) == "PreTransB")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "PreTransB")] <- "PreTrans"
      }
    }
    if(length(which(colnames(target_file) == "Wk1b")) > 0) {
      if(length(which(colnames(target_file) == "Wk1")) > 0) {
        target_file[,"Wk1"] <- target_file[,"Wk1"] + target_file[,"Wk1b"]
        target_file <- target_file[,-which(colnames(target_file) == "Wk1b")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "Wk1b")] <- "Wk1"
      }
    }
    if(length(which(colnames(target_file) == "Wk-1Run1")) > 0) {
      if(length(which(colnames(target_file) == "Wk-1")) > 0) {
        target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run1"]
        target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run1")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "Wk-1Run1")] <- "Wk-1"
      }
    }
    if(length(which(colnames(target_file) == "Wk-1Run2")) > 0) {
      if(length(which(colnames(target_file) == "Wk-1")) > 0) {
        target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run2"]
        target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run2")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "Wk-1Run2")] <- "Wk-1"
      }
    }
    
    ### remove all zero time points
    time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
    
    ### get time points except the Total
    time_points <- setdiff(time_points, c("Total"))
    
    ### remove before GMP time points
    time_points <- intersect(time_points, gmp_after_time_points)
    
    ### draw when there are at least two time points
    if(length(time_points) > 1) {
      ###  get lineages
      lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
        return(length(which(x > 0)) > 1)  
      })),time_points]
      
      ### get an input data frame for the alluvial plot
      total_rows <- length(which(lineage_table[,time_points] > 0))
      if(total_rows > 0) {
        plot_df <- data.frame(Time=rep("", total_rows),
                              Clone_Size=rep(0, total_rows),
                              Clone=rep("", total_rows),
                              CDR3=rep("", total_rows))
        cnt <- 1
        for(i in 1:nrow(lineage_table)) {
          for(tp in time_points) {
            if(lineage_table[i,tp] > 0) {
              plot_df[cnt,] <- c(tp,
                                 lineage_table[i,tp],
                                 rownames(lineage_table)[i],
                                 "CDR3")
              cnt <- cnt + 1
            }
          }
        }
        plot_df$Time <- factor(plot_df$Time, levels = intersect(time_points, unique(plot_df$Time)))
        
        ### numerize the clone_size column
        plot_df$Clone_Size <- as.numeric(plot_df$Clone_Size)
        
        ### draw an alluvial plot
        p[[patient]] <- ggplot(plot_df,
               aes(x = Time, stratum = Clone, alluvium = Clone,
                   y = Clone_Size,
                   fill = Clone, label = Clone)) +
          ggtitle(paste(patient)) +
          geom_flow() +
          geom_stratum(alpha = 1) +
          # geom_text(stat = "stratum", size = 2) +
          rotate_x_text(90) +
          theme_pubr(legend = "none") +
          # theme_cleveland2() +
          scale_fill_viridis(discrete = T) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
          theme_classic(base_size = 36) +
          theme(axis.text.x = element_text(size = 30),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 30),
                legend.position = "none")
        ggsave(file = paste0(outputDir2, "Car+_Clonal_Tracing_", patient, ".png"),
               plot = p[[patient]],
               width = 20, height = 10, dpi = 350)
        
        ### combine the plot data into one
        if(is.null(total_plot_df)) {
          total_plot_df <- plot_df
        } else {
          total_plot_df <- rbind(total_plot_df, plot_df)
        }
      }
    }
    
    gc()
  }
  
  ### draw an alluvial plot with all-patient-combined plot data
  total_plot_df$Time <- factor(total_plot_df$Time, levels = intersect(total_time_points, unique(total_plot_df$Time)))
  ggplot(total_plot_df,
         aes(x = Time, stratum = Clone, alluvium = Clone,
             y = Clone_Size,
             fill = Clone, label = Clone)) +
    ggtitle(paste("All patients - Px00-Px15")) +
    geom_flow() +
    geom_stratum(alpha = 1) +
    # geom_text(stat = "stratum", size = 2) +
    rotate_x_text(90) +
    theme_pubr(legend = "none") +
    # theme_cleveland2() +
    scale_fill_viridis(discrete = T) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30),
          legend.position = "none")
  ggsave(file = paste0(outputDir2, "All_Car+_Clonal_Tracing.png"), width = 20, height = 10, dpi = 350)
  
  ### combine plots of the selected 9 patients into one
  g <- arrangeGrob(grobs = p[c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                               "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                               "SJCAR19-09", "SJCAR19-10", "SJCAR19-11")],
                   nrow = 3,
                   ncol = 3,
                   top = "")
  ggsave(file = paste0(outputDir2, "Car+_Clonal_Tracing_9pxs_in_one.png"), g, width = 40, height = 15, dpi = 350)
  
  
  ### 2. Number of lineages & clonotypes table CD4+, CD8+, CAR+, CAR-, per patient
  ### (alpha only, beta only, one from each, strict)
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/2/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### set table colnames & rownames
  table_colnames <- c("Alpha", "Beta", "One_From_Each", "Strict")
  clonotype_types <- c("clonotype_id_by_patient_alpha", "clonotype_id_by_patient_beta",
                       "clonotype_id_by_patient_one_alpha_beta", "clonotype_id_by_patient")
  names(table_colnames) <- clonotype_types
  names(clonotype_types) <- table_colnames
  table_rownames <- c("CD4+/CAR+ Cell #", "CD4+/CAR- Cell #",
                      "CD8+/CAR+ Cell #", "CD8+/CAR- Cell #",
                      "CD4+/CAR+ Clonotype #", "CD4+/CAR- Clonotype #",
                      "CD8+/CAR+ Clonotype #", "CD8+/CAR- Clonotype #",
                      "CD4+/CAR+ Lineage #", "CD4+/CAR- Lineage #",
                      "CD8+/CAR+ Lineage #", "CD8+/CAR- Lineage #",
                      "CD4+/CAR+ GMP Lineage #", "CD8+/CAR+ GMP Lineage #")
  px_table_rownames <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) paste(x, table_rownames)))
  
  ### make an empty table
  result_table <- matrix(0, length(px_table_rownames), length(table_colnames))
  colnames(result_table) <- table_colnames
  rownames(result_table) <- px_table_rownames
  
  ### for each clonotype type to get lineages
  total_lineages <- vector("list", length = length(table_colnames))
  names(total_lineages) <- table_colnames
  total_gmp_lineages <- vector("list", length = length(table_colnames))
  names(total_gmp_lineages) <- table_colnames
  for(type in table_colnames) {
    total_lineages[[type]] <- NULL
    total_gmp_lineages[[type]] <- NULL
    ### run for each patient
    for(patient in unique(Seurat_Obj@meta.data$px)) {
      
      ### print progress
      writeLines(paste(patient))
      
      ### load the file
      target_file <- read.xlsx2(file = paste0(px_result_dir, patient, "/car_clonotype_frequency_over_time_", patient, ".xlsx"),
                                sheetName = paste0("CARpos_Clonotype_Frequency_", substr(type, 1, 4)), stringsAsFactors = FALSE, check.names = FALSE,
                                row.names = 1)
      
      ### numerize the table
      for(i in 1:ncol(target_file)) {
        target_file[,i] <- as.numeric(target_file[,i])
      }
      
      ### combine some redundant time points to one
      if(length(which(colnames(target_file) == "GMP-redo")) > 0) {
        target_file[,"GMP"] <- target_file[,"GMP"] + target_file[,"GMP-redo"]
        target_file <- target_file[,-which(colnames(target_file) == "GMP-redo")]
      }
      if(length(which(colnames(target_file) == "PreTransB")) > 0) {
        if(length(which(colnames(target_file) == "PreTrans")) > 0) {
          target_file[,"PreTrans"] <- target_file[,"PreTrans"] + target_file[,"PreTransB"]
          target_file <- target_file[,-which(colnames(target_file) == "PreTransB")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "PreTransB")] <- "PreTrans"
        }
      }
      if(length(which(colnames(target_file) == "Wk1b")) > 0) {
        if(length(which(colnames(target_file) == "Wk1")) > 0) {
          target_file[,"Wk1"] <- target_file[,"Wk1"] + target_file[,"Wk1b"]
          target_file <- target_file[,-which(colnames(target_file) == "Wk1b")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "Wk1b")] <- "Wk1"
        }
      }
      if(length(which(colnames(target_file) == "Wk-1Run1")) > 0) {
        if(length(which(colnames(target_file) == "Wk-1")) > 0) {
          target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run1"]
          target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run1")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "Wk-1Run1")] <- "Wk-1"
        }
      }
      if(length(which(colnames(target_file) == "Wk-1Run2")) > 0) {
        if(length(which(colnames(target_file) == "Wk-1")) > 0) {
          target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run2"]
          target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run2")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "Wk-1Run2")] <- "Wk-1"
        }
      }
      
      ### remove all zero time points
      time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
      
      ### get time points except the Total
      time_points <- setdiff(time_points, c("Total"))
      
      ### draw when there are at least two time points
      if(length(time_points) > 1) {
        ###  get lineages
        lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
          return(length(which(x > 0)) > 1)  
        })),time_points]
        lineages <- rownames(lineage_table)
        
        ### get gmp lineages
        gmp_lineage_table <- lineage_table[which(lineage_table$GMP > 0),]
        gmp_lineages <- rownames(gmp_lineage_table)
        
        ### combine lineages
        if(length(lineages) > 0) {
          ### combine
          if(is.null(total_lineages[[type]])) {
            total_lineages[[type]] <- lineages
          } else {
            total_lineages[[type]] <- c(total_lineages[[type]], lineages)
          }
        }
        
        ### combine gmp lineages
        if(length(gmp_lineages) > 0) {
          ### combine
          if(is.null(total_gmp_lineages[[type]])) {
            total_gmp_lineages[[type]] <- gmp_lineages
          } else {
            total_gmp_lineages[[type]] <- c(total_gmp_lineages[[type]], gmp_lineages)
          }
        }
      }
      
      gc()
    }
  }
  
  ### some pre-calculated indices
  cd4_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4"),
                              which(Seurat_Obj@meta.data$CAR == "CARpos"))
  cd4_carneg_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4"),
                              which(Seurat_Obj@meta.data$CAR == "CARneg"))
  cd8_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"),
                              which(Seurat_Obj@meta.data$CAR == "CARpos"))
  cd8_carneg_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"),
                              which(Seurat_Obj@meta.data$CAR == "CARneg"))
  lineages_idx <- vector("list", length = length(total_lineages))
  names(lineages_idx) <- names(total_lineages)
  for(type in names(lineages_idx)) {
    lineages_idx[[type]] <- which(Seurat_Obj@meta.data[,clonotype_types[type]] %in% total_lineages[[type]])
  }
  gmp_lineages_idx <- vector("list", length = length(total_gmp_lineages))
  names(gmp_lineages_idx) <- names(total_gmp_lineages)
  for(type in names(gmp_lineages_idx)) {
    gmp_lineages_idx[[type]] <- which(Seurat_Obj@meta.data[,clonotype_types[type]] %in% total_gmp_lineages[[type]])
  }
  
  ### fill out the table
  for(i in 1:length(table_colnames)) {
    for(j in 1:length(unique(Seurat_Obj@meta.data$px))) {
      
      ### fill out for each column and for each patient
      px <- unique(Seurat_Obj@meta.data$px)[j]
      result_table[paste(px, "CD4+/CAR+ Cell #"),
                   table_colnames[i]] <- length(intersect(which(Seurat_Obj@meta.data$px == px),
                                                          cd4_carpos_idx))
      result_table[paste(px, "CD4+/CAR- Cell #"),
                   table_colnames[i]] <- length(intersect(which(Seurat_Obj@meta.data$px == px),
                                                          cd4_carneg_idx))
      result_table[paste(px, "CD8+/CAR+ Cell #"),
                   table_colnames[i]] <- length(intersect(which(Seurat_Obj@meta.data$px == px),
                                                          cd8_carpos_idx))
      result_table[paste(px, "CD8+/CAR- Cell #"),
                   table_colnames[i]] <- length(intersect(which(Seurat_Obj@meta.data$px == px),
                                                          cd8_carneg_idx))
      result_table[paste(px, "CD4+/CAR+ Clonotype #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                      cd4_carpos_idx),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD4+/CAR- Clonotype #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                      cd4_carneg_idx),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR+ Clonotype #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                      cd8_carpos_idx),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR- Clonotype #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                      cd8_carneg_idx),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD4+/CAR+ Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd4_carpos_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD4+/CAR- Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd4_carneg_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR+ Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd8_carpos_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR- Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd8_carneg_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD4+/CAR+ GMP Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(gmp_lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd4_carpos_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR+ GMP Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(gmp_lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd8_carpos_idx)),
                                                                            clonotype_types[i]]))
      
    }
  }
  
  ### save the result table
  write.xlsx2(result_table,
              sheetName = "Lineage_Statistics_Table",
              file = paste0(outputDir2, "Lineage_Statistics_Table_Per_Px.xlsx"))
  
  ### make a px-combined table
  specific_pxs <- c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                    "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                    "SJCAR19-09", "SJCAR19-10", "SJCAR19-11")
  specifix_pxs_idx <- as.vector(sapply(specific_pxs, function(x) grep(x, rownames(result_table), fixed = TRUE)))
  total_result_table <- matrix(0, length(table_rownames), length(table_colnames))
  rownames(total_result_table) <- table_rownames
  colnames(total_result_table) <- table_colnames
  for(r in table_rownames) {
    for(c in table_colnames) {
      total_result_table[r,c] <- sum(result_table[intersect(specifix_pxs_idx,
                                                            grep(r, rownames(result_table), fixed = TRUE)),
                                                  c])
    }
  }
  
  ### save the total result table
  write.xlsx2(total_result_table,
              sheetName = "Total_Lineage_Statistics_Table",
              file = paste0(outputDir2, "Lineage_Statistics_Table_Total.xlsx"))
  
  
  #
  ### 3. CAR+ % for each library - proportional bar graph & table - CAR+, CAR-, CD4+, CD8+
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/3/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### some pre-calculated indices
  cd4_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4"),
                              which(Seurat_Obj@meta.data$CAR == "CARpos"))
  cd4_carneg_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4"),
                              which(Seurat_Obj@meta.data$CAR == "CARneg"))
  cd8_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"),
                              which(Seurat_Obj@meta.data$CAR == "CARpos"))
  cd8_carneg_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"),
                              which(Seurat_Obj@meta.data$CAR == "CARneg"))
  
  ### examine the proportion of cells in each patient - barplot 
  plot_df <- data.frame(LIB=as.vector(sapply(unique(Seurat_Obj@meta.data$library), function(x) rep(x,4))),
                        TYPE=rep(c("CD4+ CAR+", "CD4+ CAR-", "CD8+ CAR+", "CD8+ CAR-"),
                                length(unique(Seurat_Obj@meta.data$library))),
                        NUM=0,
                        PCNT=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:nrow(plot_df)) {
    type <- strsplit(plot_df$TYPE[i], split = " ", fixed = TRUE)[[1]]
    if(type[1] == "CD4+" && type[2] == "CAR+") {
      target_idx <- cd4_carpos_idx
    } else if(type[1] == "CD4+" && type[2] == "CAR-") {
      target_idx <- cd4_carneg_idx
    } else if(type[1] == "CD8+" && type[2] == "CAR+") {
      target_idx <- cd8_carpos_idx
    } else {
      target_idx <- cd8_carneg_idx
    }
    plot_df$NUM[i] <- length(intersect(which(Seurat_Obj@meta.data$library == plot_df$LIB[i]),
                                       target_idx))
    plot_df$PCNT[i] <- round(plot_df$NUM[i]*100/length(which(Seurat_Obj@meta.data$library == plot_df$LIB[i])), digits = 0)
  }
  ### pcnt < 5 -> ""
  plot_df$PCNT[which(as.numeric(plot_df$PCNT) < 5)] <- ""
  plot_df$PCNT <- as.character(plot_df$PCNT)
  ggplot(data=plot_df, aes_string(x="LIB", y="NUM", fill="TYPE", label="PCNT")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Cell %") +
    geom_text(size = 3, position = position_stack(vjust = 1)) +
    coord_flip() +
    scale_y_continuous(expand = c(0,300)) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 10),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "Cell_Distribution_In_Each_Library.png"), width = 20, height = 10, dpi = 400)
  
  
  #
  ### 4. Clone size between CAR+ lineages vs non-lineage CAR+ after infusion
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/4/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### AFTER-INFUSION - GMP SUBSISTER CAR+ VS THE REST CAR+
  
  ### get AI indicies
  car_pi_indicies <- intersect(which(Seurat_Obj@meta.data$time2 %in% c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo")),
                               which(Seurat_Obj@meta.data$CAR == "CARpos"))
  
  ### target lineages
  subsister_clones_pi <- unique(intersect(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")],
                                          Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[car_pi_indicies]))
  rest_carpos_clones_pi <- unique(setdiff(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[car_pi_indicies],
                                          subsister_clones_pi))
  rest_carpos_clones_pi <- rest_carpos_clones_pi[which(!is.na(rest_carpos_clones_pi))]
  
  ### get persister vs rest carpos clone sizes in PI time points
  subsister_clone_sizes_pi <- rep(1, length(subsister_clones_pi))
  names(subsister_clone_sizes_pi) <- subsister_clones_pi
  for(clone in subsister_clones_pi) {
    subsister_clone_sizes_pi[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[car_pi_indicies] == clone))
  }
  
  rest_carpos_clone_sizes_pi <- rep(0, length(rest_carpos_clones_pi))
  names(rest_carpos_clone_sizes_pi) <- rest_carpos_clones_pi
  for(clone in rest_carpos_clones_pi) {
    rest_carpos_clone_sizes_pi[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[car_pi_indicies] == clone))
  }
  
  ### prepare a data frame for the plot
  plot_df <- data.frame(Clone_Size=c(subsister_clone_sizes_pi,
                                     rest_carpos_clone_sizes_pi),
                        Group=c(rep("Post_Infusion_CARpos_Subsisters", length(subsister_clone_sizes_pi)),
                                rep("Post_Infusion_CARpos_Rest", length(rest_carpos_clone_sizes_pi))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df$Group <- factor(plot_df$Group, levels = c("Post_Infusion_CARpos_Subsisters", "Post_Infusion_CARpos_Rest"))
  
  ### draw a violin plot
  ggplot(plot_df, aes_string(x="Group", y="Clone_Size", fill = "Group")) +
    geom_violin(trim=FALSE)+
    geom_text(data = data.frame(Group=c("Post_Infusion_CARpos_Subsisters", "Post_Infusion_CARpos_Rest"),
                                Median=c(median(subsister_clone_sizes_pi),
                                         median(rest_carpos_clone_sizes_pi)),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Median"),
              size = 10, hjust = -2, vjust = -2) +
    geom_text(data = data.frame(Group=c("Post_Infusion_CARpos_Subsisters", "Post_Infusion_CARpos_Rest"),
                                Median=c(median(subsister_clone_sizes_pi),
                                         median(rest_carpos_clone_sizes_pi)),
                                Length=c(paste("n =", length(subsister_clone_sizes_pi)),
                                         paste("n =", length(rest_carpos_clone_sizes_pi))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Length"),
              size = 5, hjust = 0.5, vjust = 3) +
    labs(title="", x="", y = "Clone Size", fill = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    stat_summary(fun="median") +
    # stat_compare_means(size = 8) +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 25),
          axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 30))
  ggsave(paste0(outputDir2, "/", "PI_CARpos_Violin_Clone_Size_S_vs_R.png"), width = 22, height = 12, dpi = 500)
  
  ### draw a violin plot
  ggplot(plot_df, aes_string(x="Group", y="Clone_Size", fill = "Group")) +
    geom_violin(trim=FALSE)+
    ylim(c(-0.5, 2)) +
    geom_text(data = data.frame(Group=c("Post_Infusion_CARpos_Subsisters", "Post_Infusion_CARpos_Rest"),
                                Median=c(median(subsister_clone_sizes_pi),
                                         median(rest_carpos_clone_sizes_pi)),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Median"),
              size = 10, hjust = -2, vjust = -2) +
    geom_text(data = data.frame(Group=c("Post_Infusion_CARpos_Subsisters", "Post_Infusion_CARpos_Rest"),
                                Median=c(median(subsister_clone_sizes_pi),
                                         median(rest_carpos_clone_sizes_pi)),
                                Length=c(paste("n =", length(subsister_clone_sizes_pi)),
                                         paste("n =", length(rest_carpos_clone_sizes_pi))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Length"),
              size = 5, hjust = 0.5, vjust = 35) +
    labs(title="", x="", y = "Clone Size", fill = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    stat_summary(fun="median") +
    # stat_compare_means(size = 8) +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 25),
          axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 30))
  ggsave(paste0(outputDir2, "/", "PI_CARpos_Violin_Clone_Size_S_vs_R_CUT.png"), width = 22, height = 12, dpi = 500)
  
  ### density plot
  ggplot(plot_df, aes_string(x="Clone_Size", col="Group")) +
    geom_density(size = 2) +
    xlim(c(0,20)) +
    geom_text(data = data.frame(Group=c("Post_Infusion_CARpos_Subsisters", "Post_Infusion_CARpos_Rest"),
                                x = c(18, 18),
                                y = c(0.3, 0.2),
                                Length=c(paste("n =", length(subsister_clone_sizes_pi)),
                                         paste("n =", length(rest_carpos_clone_sizes_pi))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "x", y = "y", label = "Length"),
              size = 5, hjust = 0.5, vjust = 0.5) +
    labs(title="", x="Clone Size", y = "Density", col = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_classic(base_size = 36) +
    theme(legend.text = element_text(size = 35))
  ggsave(paste0(outputDir2, "/", "PI_CARpos_Density_Clone_Size_S_vs_R.png"), width = 18, height = 12, dpi = 500)
  
  ### target lineages
  subsister_clones_gmp <- unique(target_Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  rest_carpos_clones_gmp <- unique(target_Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "NO")])
  rest_carpos_clones_gmp <- rest_carpos_clones_gmp[which(!is.na(rest_carpos_clones_gmp))]
  
  ### gmp indicies
  gmp_carpos_cd8_indicies <- intersect(intersect(which(Seurat_Obj@meta.data$time2 == "GMP"),
                                                 which(Seurat_Obj@meta.data$CAR == "CARpos")),
                                       which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### get persister vs rest carpos clone sizes in GMP
  subsister_clone_sizes_gmp <- rep(1, length(subsister_clones_gmp))
  names(subsister_clone_sizes_gmp) <- subsister_clones_gmp
  for(clone in subsister_clones_gmp) {
    subsister_clone_sizes_gmp[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[gmp_carpos_cd8_indicies] == clone))
  }
  
  rest_carpos_clone_sizes_gmp <- rep(0, length(rest_carpos_clones_gmp))
  names(rest_carpos_clone_sizes_gmp) <- rest_carpos_clones_gmp
  for(clone in rest_carpos_clones_gmp) {
    rest_carpos_clone_sizes_gmp[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[gmp_carpos_cd8_indicies] == clone))
  }
  
  ### prepare a data frame for the plot
  plot_df <- data.frame(Clone_Size=c(subsister_clone_sizes_gmp,
                                     rest_carpos_clone_sizes_gmp),
                        Group=c(rep("GMP_CARpos_CD8_Subsisters", length(subsister_clone_sizes_gmp)),
                                rep("GMP_CARpos_CD8_Rest", length(rest_carpos_clone_sizes_gmp))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df$Group <- factor(plot_df$Group, levels = c("GMP_CARpos_CD8_Subsisters", "GMP_CARpos_CD8_Rest"))
  
  ### draw a violin plot
  ggplot(plot_df, aes_string(x="Group", y="Clone_Size", fill = "Group")) +
    geom_violin(trim=FALSE)+
    geom_text(data = data.frame(Group=c("GMP_CARpos_CD8_Subsisters", "GMP_CARpos_CD8_Rest"),
                                Median=c(median(subsister_clone_sizes_gmp),
                                         median(rest_carpos_clone_sizes_gmp)),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Median"),
              size = 10, hjust = -2, vjust = -2) +
    geom_text(data = data.frame(Group=c("GMP_CARpos_CD8_Subsisters", "GMP_CARpos_CD8_Rest"),
                                Median=c(median(subsister_clone_sizes_gmp),
                                         median(rest_carpos_clone_sizes_gmp)),
                                Length=c(paste("n =", length(subsister_clone_sizes_gmp)),
                                         paste("n =", length(rest_carpos_clone_sizes_gmp))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Length"),
              size = 5, hjust = 0.5, vjust = 3) +
    labs(title="", x="", y = "Clone Size", fill = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    stat_summary(fun="median") +
    # stat_compare_means(size = 8) +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 25),
          axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 30))
  ggsave(paste0(outputDir2, "/", "GMP_CARpos_CD8_Violin_Clone_Size_S_vs_R.png"), width = 22, height = 12, dpi = 500)
  
  ### draw a violin plot
  ggplot(plot_df, aes_string(x="Group", y="Clone_Size", fill = "Group")) +
    geom_violin(trim=FALSE)+
    ylim(c(0.5, 2)) +
    geom_text(data = data.frame(Group=c("GMP_CARpos_CD8_Subsisters", "GMP_CARpos_CD8_Rest"),
                                Median=c(median(subsister_clone_sizes_gmp),
                                         median(rest_carpos_clone_sizes_gmp)),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Median"),
              size = 10, hjust = -2, vjust = -2) +
    geom_text(data = data.frame(Group=c("GMP_CARpos_CD8_Subsisters", "GMP_CARpos_CD8_Rest"),
                                Median=c(median(subsister_clone_sizes_gmp),
                                         median(rest_carpos_clone_sizes_gmp)),
                                Length=c(paste("n =", length(subsister_clone_sizes_gmp)),
                                         paste("n =", length(rest_carpos_clone_sizes_gmp))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Length"),
              size = 5, hjust = 0.5, vjust = 15) +
    labs(title="", x="", y = "Clone Size", fill = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    stat_summary(fun="median") +
    # stat_compare_means(size = 8) +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 25),
          axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 30))
  ggsave(paste0(outputDir2, "/", "GMP_CARpos_CD8_Violin_Clone_Size_S_vs_R_CUT.png"), width = 22, height = 12, dpi = 500)
  
  ### density plot
  ggplot(plot_df, aes_string(x="Clone_Size", col="Group")) +
    geom_density(size = 2) +
    xlim(c(0,20)) +
    geom_text(data = data.frame(Group=c("GMP_CARpos_CD8_Subsisters", "GMP_CARpos_CD8_Rest"),
                                x = c(18, 18),
                                y = c(0.4, 0.2),
                                Length=c(paste("n =", length(subsister_clone_sizes_gmp)),
                                         paste("n =", length(rest_carpos_clone_sizes_gmp))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "x", y = "y", label = "Length"),
              size = 5, hjust = 0.5, vjust = 0.5) +
    labs(title="", x="Clone Size", y = "Density", col = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_classic(base_size = 36) +
    theme(legend.text = element_text(size = 35))
  ggsave(paste0(outputDir2, "/", "GMP_CARpos_CD8_Density_Clone_Size_S_vs_R.png"), width = 18, height = 12, dpi = 500)
  
  ### draw a cumulative plot
  max_clone_size <- Reduce(max, c(subsister_clone_sizes_gmp, rest_carpos_clone_sizes_gmp,
                                  subsister_clone_sizes_pi, rest_carpos_clone_sizes_pi))
  png(filename = paste0(outputDir2, "/", "Cumulative_Clone_Sizes.png"), width = 2000, height = 1500, res = 350)
  plot(ecdf(subsister_clone_sizes_gmp), col="#F2C03E", lwd = 1,
       main = "Cumulative Clone Sizes",
       sub = paste0("KS test p-value: ",
                    "GMP=", formatC(ks.test(subsister_clone_sizes_gmp, rest_carpos_clone_sizes_gmp)$p.value, format = "e", digits = 2), ", ",
                    "PI=", formatC(ks.test(subsister_clone_sizes_pi, rest_carpos_clone_sizes_pi)$p.value, format = "e", digits = 2)),
       xlab = "Clone Size",
       ylab = "Proportion",
       xaxs="i", xlim=c(0,max_clone_size))
  lines(ecdf(rest_carpos_clone_sizes_gmp), col="#045282", lwd = 1)
  lines(ecdf(subsister_clone_sizes_pi), col="#D21414", lwd = 1)
  lines(ecdf(rest_carpos_clone_sizes_pi), col="#039076", lwd = 1)
  legend("bottomright", 
         legend=c("GMP CAR+ CD8 Subsisters (n=198)", "GMP CAR+ CD8 Rest (n=47994)",
                  "PI CAR Subsister (n=198)","PI CAR Rest (n=45448)"),
         col=c("#F2C03E", "#045282","#D21414","#039076"),
         pch=15)
  dev.off()
  
  ### with downsampled ones
  set.seed(1234)
  rest_carpos_clone_sizes_gmp_ds <- sample(rest_carpos_clone_sizes_gmp, length(subsister_clone_sizes_gmp))
  rest_carpos_clone_sizes_pi_ds <- sample(rest_carpos_clone_sizes_pi, length(subsister_clone_sizes_pi))
  
  ### draw a cumulative plot
  max_clone_size <- Reduce(max, c(subsister_clone_sizes_gmp, rest_carpos_clone_sizes_gmp_ds,
                                  subsister_clone_sizes_pi, rest_carpos_clone_sizes_pi_ds))
  png(filename = paste0(outputDir2, "/", "Cumulative_Clone_Sizes_Downsampled.png"), width = 2000, height = 1500, res = 350)
  plot(ecdf(subsister_clone_sizes_gmp), col="#F2C03E", lwd = 1,
       main = "Cumulative Clone Sizes",
       sub = paste0("KS test p-value: ",
                    "GMP=", formatC(ks.test(subsister_clone_sizes_gmp, rest_carpos_clone_sizes_gmp_ds)$p.value, format = "e", digits = 2), ", ",
                    "PI=", formatC(ks.test(subsister_clone_sizes_pi, rest_carpos_clone_sizes_pi_ds)$p.value, format = "e", digits = 2)),
       xlab = "Clone Size",
       ylab = "Proportion",
       xaxs="i", xlim=c(0,max_clone_size))
  lines(ecdf(rest_carpos_clone_sizes_gmp_ds), col="#045282", lwd = 1)
  lines(ecdf(subsister_clone_sizes_pi), col="#D21414", lwd = 1)
  lines(ecdf(rest_carpos_clone_sizes_pi_ds), col="#039076", lwd = 1)
  legend("bottomright", 
         legend=c("GMP CAR+ CD8 Subsisters (n=198)", "GMP CAR+ CD8 Rest (n=198)",
                  "PI CAR Subsister (n=198)","PI CAR Rest (n=198)"),
         col=c("#F2C03E", "#045282","#D21414","#039076"),
         pch=15)
  dev.off()
  
  
  #
  ### 5. a) PCA/UMAP plot of CAR+s over time (coloring based on time)
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/5/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get CARpos-only seurat object
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$CAR)
  sub_seurat_obj <- subset(Seurat_Obj, idents = c("CARpos"))
  
  ### after gmp time points only
  after_gmp_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6",
                             "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$time2)
  sub_seurat_obj <- subset(sub_seurat_obj, idents = intersect(after_gmp_time_points,
                                                              unique(sub_seurat_obj@meta.data$time2)))
  
  #
  ### run pca on the sub_seurat_obj
  #
  ### normalization
  sub_seurat_obj <- NormalizeData(sub_seurat_obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  sub_seurat_obj <- FindVariableFeatures(sub_seurat_obj,
                                         selection.method = "vst", nfeatures = 2000)
  ### scaling
  sub_seurat_obj <- ScaleData(sub_seurat_obj,
                              vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  ### PCA
  sub_seurat_obj <- RunPCA(sub_seurat_obj,
                           features = VariableFeatures(object = sub_seurat_obj),
                           npcs = 15)
  ### UMAP
  sub_seurat_obj <- RunUMAP(sub_seurat_obj, dims = 1:15)
  
  ### get seurat object for some specific patients
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$px)
  sub_seurat_obj2 <- subset(sub_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                       "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### factorize the time2 column
  sub_seurat_obj2@meta.data$time2 <- factor(sub_seurat_obj2@meta.data$time2,
                                            levels = intersect(total_time_points, sub_seurat_obj2@meta.data$time2))
  
  ### PCA with time by each patient
  p <- DimPlot(object = sub_seurat_obj2, reduction = "pca",
               group.by = "time2", split.by = "px",
               pt.size = 3, ncol = 3) +
    ggtitle("") +
    labs(color="Time") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30)) +
    scale_colour_brewer(palette = "OrRd")
  # p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "PCA_CARpos_Time.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### PCA with time with all patients
  p <- DimPlot(object = sub_seurat_obj2, reduction = "pca",
               group.by = "time2",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Time") +
    scale_colour_brewer(palette = "OrRd") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "PCA_CARpos_Time_ALL.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  ### change UMAP coordinates consistency with the previous result
  ### doesn't affect the result itself
  sub_seurat_obj2@reductions$umap@cell.embeddings <- -sub_seurat_obj2@reductions$umap@cell.embeddings
  sub_seurat_obj2@reductions$umap@feature.loadings <- -sub_seurat_obj2@reductions$umap@feature.loadings
  
  ### UMAP with time by each patient
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "time2", split.by = "px",
               pt.size = 3, ncol = 3) +
    ggtitle("") +
    labs(color="Time") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30)) +
    scale_colour_brewer(palette = "OrRd")
  # p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "UMAP_CARpos_Time.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### UMAP with time with all patients
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "time2",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Time") +
    scale_colour_brewer(palette = "OrRd") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "UMAP_CARpos_Time_ALL.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  #
  ### 5. b) Clustering + UMAP
  #
  
  ### perform clustering
  sub_seurat_obj2 <- FindNeighbors(sub_seurat_obj2, dims = 1:10)
  sub_seurat_obj2 <- FindClusters(sub_seurat_obj2, resolution = 0.5)
  
  ### save the clustering result to meta.data
  sub_seurat_obj2@meta.data$clusters <- Idents(sub_seurat_obj2)
  
  ### UMAP with clusters by each patient
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "clusters", split.by = "px",
               pt.size = 3, ncol = 3) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "UMAP_CARpos_Clusters.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### UMAP with clusters with all patients
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "clusters",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "UMAP_CARpos_Clusters_ALL.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  ### cluster0 - UMAP
  sub_seurat_obj2@meta.data$pi_cluster0 <- as.character(sub_seurat_obj2@meta.data$clusters)
  sub_seurat_obj2@meta.data$pi_cluster0[which(sub_seurat_obj2@meta.data$pi_cluster0 != "0")] <- "Others"
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "pi_cluster0",
               pt.size = 5, cols = c("0" = "red", "Others" = "lightgray"),
               order = c("0", "Others")) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "UMAP_CARpos_Cluster0_ALL.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  #
  ### 5. c) Cluster0 vs others
  #
  
  ### set column for cluster0 and the others
  sub_seurat_obj2@meta.data$pi_cluster0 <- as.character(sub_seurat_obj2@meta.data$clusters)
  sub_seurat_obj2@meta.data$pi_cluster0[which(sub_seurat_obj2@meta.data$pi_cluster0 != "0")] <- "Others"
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$pi_cluster0)
  
  ### DE analysis
  de_result <- FindMarkers(sub_seurat_obj2,
                           ident.1 = "0",
                           ident.2 = "Others",
                           min.pct = 0.5,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_Cluster0_vs_Others.xlsx"),
              sheetName = "CARpos_Cluster0_DE_Result", row.names = FALSE)
  
  #
  ### 5. d) After infusion CAR+ subsister vs non-subsisters to find CAR is differentially expressed
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/5/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get CARpos-only seurat object
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$CAR)
  sub_seurat_obj <- subset(Seurat_Obj, idents = c("CARpos"))
  
  ### after gmp time points only
  after_gmp_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6",
                             "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$time2)
  sub_seurat_obj <- subset(sub_seurat_obj, idents = intersect(after_gmp_time_points,
                                                              unique(sub_seurat_obj@meta.data$time2)))
  
  ### normalization
  sub_seurat_obj <- NormalizeData(sub_seurat_obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  sub_seurat_obj <- FindVariableFeatures(sub_seurat_obj,
                                         selection.method = "vst", nfeatures = 2000)
  ### scaling
  sub_seurat_obj <- ScaleData(sub_seurat_obj,
                              vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  ### PCA
  sub_seurat_obj <- RunPCA(sub_seurat_obj,
                           features = VariableFeatures(object = sub_seurat_obj),
                           npcs = 15)
  ### UMAP
  sub_seurat_obj <- RunUMAP(sub_seurat_obj, dims = 1:15)
  
  ### get seurat object for some specific patients
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$px)
  sub_seurat_obj2 <- subset(sub_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                       "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### factorize the time2 column
  sub_seurat_obj2@meta.data$time2 <- factor(sub_seurat_obj2@meta.data$time2,
                                            levels = intersect(total_time_points, sub_seurat_obj2@meta.data$time2))
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  
  ### define persisters
  sub_seurat_obj2@meta.data$CD8_Persisters <- "Non-subsisters"
  sub_seurat_obj2@meta.data$CD8_Persisters[which(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones)] <- "Subsisters"
  
  ### set ident with the persistency info
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$CD8_Persisters)
  
  ### DE analysis
  de_result <- FindMarkers(sub_seurat_obj2,
                           ident.1 = "Subsisters",
                           ident.2 = "Non-subsisters",
                           min.pct = 0.5,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_Post_Infusion_Subsister_vs_Non-Subsisters.xlsx"),
              sheetName = "CARpos_PI_S_vs_NS", row.names = FALSE)
  
  ### now after infusion CAR+ subsister vs cluster (1,2,3,5,6,7,9,10)
  
  ### perform clustering
  sub_seurat_obj2 <- FindNeighbors(sub_seurat_obj2, dims = 1:10)
  sub_seurat_obj2 <- FindClusters(sub_seurat_obj2, resolution = 0.5)
  
  ### save the clustering result to meta.data
  sub_seurat_obj2@meta.data$clusters <- Idents(sub_seurat_obj2)
  
  ### set column for cluster0 and the others
  sub_seurat_obj2@meta.data$pi_cluster0 <- as.character(sub_seurat_obj2@meta.data$clusters)
  sub_seurat_obj2@meta.data$pi_cluster0[which(sub_seurat_obj2@meta.data$pi_cluster0 != "0")] <- "Others"
  
  ### combine subsister column and the cluster0 column
  ### if overlapped, overwirte it with Cluster12345678910
  sub_seurat_obj2@meta.data$new_group <- NA
  sub_seurat_obj2@meta.data$new_group[which(sub_seurat_obj2@meta.data$CD8_Persisters == "Subsisters")] <- "Subsisters"
  sub_seurat_obj2@meta.data$new_group[which(sub_seurat_obj2@meta.data$pi_cluster0 == "Others")] <- "Cluster12345678910"
  
  ### set ident with the persistency info
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$new_group)
  
  ### DE analysis
  de_result <- FindMarkers(sub_seurat_obj2,
                           ident.1 = "Subsisters",
                           ident.2 = "Cluster12345678910",
                           min.pct = 0.5,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_Post_Infusion_Subsister_vs_C12345678910.xlsx"),
              sheetName = "CARpos_PI_S_vs_C12345678910", row.names = FALSE)
  
  
  #
  ### 5. e) Where are the subsisters that close to non-subsister located in the after infusion CAR+ UMAP?
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/5/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### only get the GMP persisters and non-persisters
  target_Seurat_Obj2 <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[union(which(Seurat_Obj$GMP_CARpos_CD8_Persister == "YES"),
                                                                                        which(Seurat_Obj$GMP_CARpos_CD8_Persister == "NO"))])
  
  ### only the patients we are interested
  target_Seurat_Obj2 <- SetIdent(object = target_Seurat_Obj2,
                                 cells = rownames(target_Seurat_Obj2@meta.data),
                                 value = target_Seurat_Obj2@meta.data$px)
  target_Seurat_Obj2 <- subset(target_Seurat_Obj2, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                              "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                              "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### get the same number of samples (persister: non-persister - same ratio) and draw a UMAP plot again
  set.seed(1234)
  target_Seurat_Obj2@meta.data$New_Persistency <- NA
  multiplier_k <- 1
  for(px in unique(target_Seurat_Obj2@meta.data$px)) {
    ### get specific indicies
    px_gmp_last <- intersect(which(target_Seurat_Obj2@meta.data$px == px),
                             which(target_Seurat_Obj2@meta.data$GMP_CARpos_CD8_Persister == "YES"))
    px_gmp_not_last <- intersect(which(target_Seurat_Obj2@meta.data$px == px),
                                 which(target_Seurat_Obj2@meta.data$GMP_CARpos_CD8_Persister == "NO"))
    
    ### sampling
    if((length(px_gmp_last) > 0) && (length(px_gmp_not_last) > length(px_gmp_last))) {
      px_gmp_not_last <- sample(px_gmp_not_last, size = length(px_gmp_last)*multiplier_k)
    }
    
    ### annotate new persistency info
    target_Seurat_Obj2@meta.data$New_Persistency[px_gmp_last] <- "YES"
    target_Seurat_Obj2@meta.data$New_Persistency[px_gmp_not_last] <- "NO"
  }
  
  ### set idents with the new info
  target_Seurat_Obj2 <- SetIdent(object = target_Seurat_Obj2,
                                 cells = rownames(target_Seurat_Obj2@meta.data),
                                 value = target_Seurat_Obj2@meta.data$New_Persistency)
  
  ### only using the specific cells
  target_Seurat_Obj2 <- subset(target_Seurat_Obj2, idents = c("YES", "NO"))
  
  
  ### get DE genes
  target_Seurat_Obj2 <- SetIdent(object = target_Seurat_Obj2,
                                 cells = rownames(target_Seurat_Obj2@meta.data),
                                 value = target_Seurat_Obj2@meta.data$New_Persistency)
  de_result <- FindMarkers(target_Seurat_Obj2,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.5,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### retain the DE genes only
  de_genes <- rownames(de_result)[which(de_result$p_val_adj < 0.01)]
  target_Seurat_Obj2@assays$RNA@counts <- target_Seurat_Obj2@assays$RNA@counts[de_genes,]
  target_Seurat_Obj2@assays$RNA@data <- target_Seurat_Obj2@assays$RNA@data[de_genes,]
  target_Seurat_Obj2@assays$RNA@var.features <- de_genes
  
  ### normalization
  target_Seurat_Obj2 <- NormalizeData(target_Seurat_Obj2,
                                      normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj2 <- FindVariableFeatures(target_Seurat_Obj2,
                                             selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj2 <- ScaleData(target_Seurat_Obj2,
                                  vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj2 <- RunPCA(target_Seurat_Obj2,
                               features = VariableFeatures(object = target_Seurat_Obj2),
                               npcs = 15)
  target_Seurat_Obj2 <- RunUMAP(target_Seurat_Obj2, dims = 1:15)
  
  ### PCA with all those info
  # p <- list()
  # p[[1]] <- DimPlot(object = target_Seurat_Obj2, reduction = "pca",
  #                   group.by = "px", shape.by = "New_Persistency",
  #                   pt.size = 3, order = c("YES", "NO")) +
  #   ggtitle("PCA of SJCAR19 Data") +
  #   theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
  #   labs(color="Patient",
  #        shape="Is Persistent")
  # p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  # 
  # ### pca for each px
  # target_Seurat_Obj2 <- SetIdent(object = target_Seurat_Obj2,
  #                                cells = rownames(target_Seurat_Obj2@meta.data),
  #                                value = target_Seurat_Obj2@meta.data$px)
  # g <- ggplot_build(p[[1]])
  # color_code <- data.frame(colours = unique(g$data[[1]]["colour"]), 
  #                          label = levels(g$plot$data[, "px"]),
  #                          stringsAsFactors = FALSE, check.names = FALSE)
  # rownames(color_code) <- color_code$label
  # unique_px <- c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
  #                "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
  #                "SJCAR19-09", "SJCAR19-10", "SJCAR19-11")
  
  # for(i in 1:length(unique_px)) {
  #   target_Seurat_Obj2 <- subset(target_Seurat_Obj, idents = c(unique_px[i]))
  #   p[[i+1]] <- DimPlot(object = target_Seurat_Obj2, reduction = "pca",
  #                       group.by = "New_Persistency",
  #                       pt.size = 3, cols = c("YES" = "red", "NO" = "lightgray"),
  #                       order = c("YES", "NO")) +
  #     ggtitle(unique_px[i]) +
  #     lims(x = g$layout$panel_scales_x[[1]]$range$range,
  #          y = g$layout$panel_scales_y[[1]]$range$range) +
  #     theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
  #     labs(color="Is Persistent")
  #   # p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  # }
  
  # ### arrange the plots and save
  # fName <- paste0("PCA_Plot_Persistence_per_Px_Sampled_DE_Genes_Only")
  # rowNum <- 3
  # colNum <- 3
  # g <- arrangeGrob(grobs = p[2:length(p)],
  #                  nrow = rowNum,
  #                  ncol = colNum,
  #                  top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  # ggsave(file = paste0(outputDir2, fName, ".png"), g, width = 25, height = 15, dpi = 300)
  
  ### draw PCA with GMP time point only
  p <- DimPlot(object = target_Seurat_Obj2, reduction = "pca",
               group.by = "New_Persistency", split.by = "px",
               cols = c("NO" = "lightgray", "YES" = "red"),
               order = c("YES", "NO"),
               pt.size = 5, ncol = 3) +
    ggtitle("") +
    labs(color="") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.8
  ggsave(paste0(outputDir2, "PCA_Plot_Persistence_per_Px_Sampled_DE_Genes_Only.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### get subsister cells with PC1 < 0
  pca_map <- Embeddings(target_Seurat_Obj2, reduction = "pca")[rownames(target_Seurat_Obj2@meta.data), 1:10]
  target_idx <- intersect(which(target_Seurat_Obj2@meta.data$New_Persistency == "YES"),
                          which(pca_map[,"PC_1"] < 0))
  target_cells <- rownames(target_Seurat_Obj2@meta.data)[target_idx]
  target_clonotypes <- unique(target_Seurat_Obj2@meta.data$clonotype_id_by_patient_one_alpha_beta[target_idx])
  
  ### factorize the column
  target_Seurat_Obj2@meta.data$New_Persistency <- factor(target_Seurat_Obj2@meta.data$New_Persistency, levels = c("YES", "NO", "MID"))
  
  ### re-define the subsisters
  target_Seurat_Obj2@meta.data$New_Persistency[target_idx] <- "MID"
  
  # ### redraw the pca plots
  # g <- ggplot_build(p[[1]])
  # for(i in 1:length(unique_px)) {
  #   temp_seurat_obj <- subset(target_Seurat_Obj2, idents = c(unique_px[i]))
  #   p[[i+1]] <- DimPlot(object = temp_seurat_obj, reduction = "pca",
  #                       group.by = "New_Persistency",
  #                       pt.size = 3, cols = c("YES" = "red", "NO" = "lightgray", "MID" = "orange"),
  #                       order = c("YES", "MID", "NO")) +
  #     ggtitle(unique_px[i]) +
  #     lims(x = g$layout$panel_scales_x[[1]]$range$range,
  #          y = g$layout$panel_scales_y[[1]]$range$range) +
  #     theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
  #     labs(color="Is Persistent")
  #   p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  # }
  # 
  # ### arrange the plots and save
  # fName <- paste0("PCA_Plot_Persistence_per_Px_Sampled_DE_Genes_Only")
  # rowNum <- 3
  # colNum <- 3
  # g <- arrangeGrob(grobs = p[2:length(p)],
  #                  nrow = rowNum,
  #                  ncol = colNum,
  #                  top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  # ggsave(file = paste0(outputDir2, fName, "2.png"), g, width = 25, height = 15, dpi = 300)
  
  ### draw PCA with GMP time point only
  p <- DimPlot(object = target_Seurat_Obj2, reduction = "pca",
               group.by = "New_Persistency", split.by = "px",
               cols = c("YES" = "red", "NO" = "lightgray", "MID" = "orange"),
               order = c("NO", "MID", "YES"),
               pt.size = 5, ncol = 3) +
    ggtitle("") +
    labs(color="") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.8
  ggsave(paste0(outputDir2, "PCA_Plot_Persistence_per_Px_Sampled_DE_Genes_Only2.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  
  ### add MID annotation to after infusion CAR+ object
  sub_seurat_obj2@meta.data$CD8_Persisters[which(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% target_clonotypes)] <- "MID"
  
  ### draw UMAP only with after infusion CAR+ cells coloring with subsisters vs non-subsisters
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "CD8_Persisters", split.by = "px",
               cols = c("Non-subsisters" = "lightgray", "Subsisters" = "red", "MID" = "orange"),
               order = c("Subsisters", "MID", "Non-subsisters"),
               pt.size = 5, ncol = 3) +
    ggtitle("") +
    labs(color="Is_Persistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.8
  ggsave(paste0(outputDir2, "UMAP_CARpos_Subsister_MID.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  #
  ### 5. f) Where are the outliers (those not overlapping with Cluster 0 & 1) in the after infusion CAR+ UMAP located in the GMP PCA?
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/5/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  
  ### define persisters
  sub_seurat_obj2@meta.data$CD8_Persisters <- "Non-subsisters"
  sub_seurat_obj2@meta.data$CD8_Persisters[which(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones)] <- "Subsisters"
  
  ### get the outlier clones
  outlier_idx <- intersect(which(sub_seurat_obj2@meta.data$CD8_Persisters == "Subsisters"),
                           which(!sub_seurat_obj2@meta.data$clusters %in% c("0", "1", "4")))
  outlier_clones <- unique(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta[outlier_idx])
  
  ### make outlier column
  sub_seurat_obj2@meta.data$CD8_Persisters2 <- sub_seurat_obj2@meta.data$CD8_Persisters
  sub_seurat_obj2@meta.data$CD8_Persisters2[outlier_idx] <- "Outliers"
  
  ### draw UMAP only with after infusion CAR+ cells coloring with subsisters vs non-subsisters
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "CD8_Persisters2", split.by = "px",
               cols = c("Non-subsisters" = "lightgray", "Subsisters" = "red", "Outliers" = "black"),
               order = c("Outliers", "Subsisters", "Non-subsisters"),
               pt.size = 5, ncol = 3) +
  ggtitle("") +
  labs(color="") +
  theme_classic(base_size = 64) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
        axis.text.x = element_text(size = 48),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 48),
        legend.title = element_text(size = 36),
        legend.text = element_text(size = 36)) +
  guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.8
  ggsave(paste0(outputDir2, "UMAP_CARpos_Subsister_Outliers.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### look at the PCA now
  ### set new column for the outliers
  target_Seurat_Obj2@meta.data$New_Persistency2 <- as.character(target_Seurat_Obj2@meta.data$New_Persistency)
  target_Seurat_Obj2@meta.data$New_Persistency2[which(target_Seurat_Obj2@meta.data$New_Persistency2 == "MID")] <- "YES"
  target_Seurat_Obj2@meta.data$New_Persistency2[which(target_Seurat_Obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% outlier_clones)] <- "Outliers"
  
  ### draw PCA with GMP time point only
  p <- DimPlot(object = target_Seurat_Obj2, reduction = "pca",
               group.by = "New_Persistency2", split.by = "px",
               cols = c("NO" = "lightgray", "YES" = "red", "Outliers" = "black"),
               order = c("NO", "Outliers", "YES"),
               pt.size = 5, ncol = 3) +
    ggtitle("") +
    labs(color="") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.8
  ggsave(paste0(outputDir2, "PCA_GMP_CARpos_Subsister_Outliers.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  #
  ### 5. i) PCA/UMAP plot of lineages with size=1 vs lineages with size > 1 (coloring differently)
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/5/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get CARpos-only seurat object
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$CAR)
  sub_seurat_obj <- subset(Seurat_Obj, idents = c("CARpos"))
  
  ### after gmp time points only
  after_gmp_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6",
                             "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$time2)
  sub_seurat_obj <- subset(sub_seurat_obj, idents = intersect(after_gmp_time_points,
                                                              unique(sub_seurat_obj@meta.data$time2)))
  
  #
  ### run pca on the sub_seurat_obj
  #
  ### normalization
  sub_seurat_obj <- NormalizeData(sub_seurat_obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  sub_seurat_obj <- FindVariableFeatures(sub_seurat_obj,
                                         selection.method = "vst", nfeatures = 2000)
  ### scaling
  sub_seurat_obj <- ScaleData(sub_seurat_obj,
                              vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  ### PCA
  sub_seurat_obj <- RunPCA(sub_seurat_obj,
                           features = VariableFeatures(object = sub_seurat_obj),
                           npcs = 15)
  ### UMAP
  sub_seurat_obj <- RunUMAP(sub_seurat_obj, dims = 1:15)
  
  ### get seurat object for some specific patients
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$px)
  sub_seurat_obj2 <- subset(sub_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                       "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### factorize the time2 column
  sub_seurat_obj2@meta.data$time2 <- factor(sub_seurat_obj2@meta.data$time2,
                                            levels = intersect(total_time_points, sub_seurat_obj2@meta.data$time2))
  
  
  ### change UMAP coordinates consistency with the previous result
  ### doesn't affect the result itself
  sub_seurat_obj2@reductions$umap@cell.embeddings <- -sub_seurat_obj2@reductions$umap@cell.embeddings
  sub_seurat_obj2@reductions$umap@feature.loadings <- -sub_seurat_obj2@reductions$umap@feature.loadings
  
  ### perform clustering
  sub_seurat_obj2 <- FindNeighbors(sub_seurat_obj2, dims = 1:10)
  sub_seurat_obj2 <- FindClusters(sub_seurat_obj2, resolution = 0.5)
  
  ### save the clustering result to meta.data
  sub_seurat_obj2@meta.data$clusters <- Idents(sub_seurat_obj2)
  
  ### Cluster0 & others column
  sub_seurat_obj2@meta.data$pi_cluster0 <- as.character(sub_seurat_obj2@meta.data$clusters)
  sub_seurat_obj2@meta.data$pi_cluster0[which(sub_seurat_obj2@meta.data$pi_cluster0 != "0")] <- "Others"
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  
  ### define persisters
  sub_seurat_obj2@meta.data$CD8_Persisters <- "Non-Subsisters"
  sub_seurat_obj2@meta.data$CD8_Persisters[which(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones)] <- "Subsisters"
  
  ### define super persisters (lineage size > 1)
  super_persister_clones <- sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta[intersect(which(sub_seurat_obj2@meta.data$ALL_GMP_CARpos_Persister == "YES"),
                                                                                                       which(sub_seurat_obj2@meta.data$CD4_CD8_by_Consensus == "CD8"))]
  super_persister_clones <- unique(super_persister_clones[which(duplicated(super_persister_clones))])
  sub_seurat_obj2@meta.data$Super_Persisters <- sub_seurat_obj2@meta.data$CD8_Persisters
  sub_seurat_obj2@meta.data$Super_Persisters[which(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% super_persister_clones)] <- "Super-Subsisters"
  
  ### UMAP with subsister info split by each patient
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "Super_Persisters", split.by = "px",
               pt.size = 5, ncol = 3,
               cols = c("Non-Subsisters" = "lightgray", "Subsisters" = "red", "Super-Subsisters" = "magenta"),
               order = c("Super-Subsisters", "Subsisters", "Non-Subsisters")) +
    ggtitle("") +
    labs(color="Is_Subsistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_CARpos_Super_Subsister.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  #
  ### 5. j) Find all markers based on the clustering
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/5/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### perform clustering
  sub_seurat_obj2 <- FindNeighbors(sub_seurat_obj2, dims = 1:10)
  sub_seurat_obj2 <- FindClusters(sub_seurat_obj2, resolution = 0.5)
  
  ### save the clustering result to meta.data
  sub_seurat_obj2@meta.data$clusters <- Idents(sub_seurat_obj2)
  
  ### set cluster info as idents
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$clusters)
  
  ### DE analysis
  de_result <- FindAllMarkers(sub_seurat_obj2,
                              min.pct = 0.5,
                              logfc.threshold = 0.2,
                              test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_Clusters_AllMarkers.xlsx"),
              sheetName = "CARpos_Clusters_AllMarkers_DE_Result", row.names = FALSE)
  
  #
  ### 5. k) UMAP plot of GMP subsisters
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/5/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### only get the GMP persisters and non-persisters
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[union(which(Seurat_Obj$GMP_CARpos_CD8_Persister == "YES"),
                                                                                       which(Seurat_Obj$GMP_CARpos_CD8_Persister == "NO"))])
  
  ### only the patients we are interested
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  target_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                            "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                            "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### normalization
  target_Seurat_Obj <- NormalizeData(target_Seurat_Obj,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj <- FindVariableFeatures(target_Seurat_Obj,
                                            selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj <- ScaleData(target_Seurat_Obj,
                                 vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj <- RunPCA(target_Seurat_Obj,
                              features = VariableFeatures(object = target_Seurat_Obj),
                              npcs = 15)
  target_Seurat_Obj <- RunUMAP(target_Seurat_Obj, dims = 1:15)
  
  ### draw a UMAP with GMP CAR+ CD8 cells only by patient
  p <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
               group.by = "GMP_CARpos_CD8_Persister", split.by = "px",
               cols = c("YES" = "red", "NO" = "lightgray"),
               order = c("YES", "NO"),
               pt.size = 5, ncol = 3) +
    ggtitle("") +
    labs(color="Is_Subsistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.8
  ggsave(paste0(outputDir2, "UMAP_GMP_CARpos_CD8_Subsister_vs_Non-Subsister.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### draw a UMAP with GMP CAR+ CD8 cells only ALL patients
  p <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
               group.by = "GMP_CARpos_CD8_Persister",
               cols = c("YES" = "red", "NO" = "lightgray"),
               order = c("YES", "NO"),
               pt.size = 5) +
    ggtitle("") +
    labs(color="Is_Subsistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.8
  ggsave(paste0(outputDir2, "UMAP_GMP_CARpos_CD8_Subsister_vs_Non-Subsister_ALL.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### perform clustering
  target_Seurat_Obj <- FindNeighbors(target_Seurat_Obj, dims = 1:10)
  target_Seurat_Obj <- FindClusters(target_Seurat_Obj, resolution = 0.01)
  
  ### save the clustering result to meta.data
  target_Seurat_Obj@meta.data$gmp_clusters <- Idents(target_Seurat_Obj)
  
  ### UMAP
  p <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
               group.by = "gmp_clusters", pt.size = 3,
               order = c("1", "0")) +
    ggtitle("") +
    labs(color="") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.8
  ggsave(paste0(outputDir2, "UMAP_GMP_CARpos_CD8_Two_Clusters.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### DE analysis
  de_result_01 <- FindMarkers(target_Seurat_Obj,
                              ident.1 = "0",
                              ident.2 = "1",
                              min.pct = 0.2,
                              logfc.threshold = 0.2,
                              test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result_01),
                         de_result_01,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CARpos_CD8_Two_Clusters_DE_Result.xlsx"),
              sheetName = "GMP_CARpos_CD8_Two_Clusters_DE_Result", row.names = FALSE)
  
  
  #
  ### 6. Visualize some interesting genes on UMAP
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/6/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get CARpos-only seurat object
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$CAR)
  sub_seurat_obj <- subset(Seurat_Obj, idents = c("CARpos"))
  
  ### after gmp time points only
  after_gmp_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6",
                             "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$time2)
  sub_seurat_obj <- subset(sub_seurat_obj, idents = intersect(after_gmp_time_points,
                                                              unique(sub_seurat_obj@meta.data$time2)))
  
  #
  ### run pca on the sub_seurat_obj
  #
  ### normalization
  sub_seurat_obj <- NormalizeData(sub_seurat_obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  sub_seurat_obj <- FindVariableFeatures(sub_seurat_obj,
                                         selection.method = "vst", nfeatures = 2000)
  ### scaling
  sub_seurat_obj <- ScaleData(sub_seurat_obj,
                              vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  ### PCA
  sub_seurat_obj <- RunPCA(sub_seurat_obj,
                           features = VariableFeatures(object = sub_seurat_obj),
                           npcs = 15)
  ### UMAP
  sub_seurat_obj <- RunUMAP(sub_seurat_obj, dims = 1:15)
  
  ### get seurat object for some specific patients
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$px)
  sub_seurat_obj2 <- subset(sub_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                       "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### factorize the time2 column
  sub_seurat_obj2@meta.data$time2 <- factor(sub_seurat_obj2@meta.data$time2,
                                            levels = intersect(total_time_points, sub_seurat_obj2@meta.data$time2))
  
  ### UMAP with px info
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "px",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Patient") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "UMAP_CARpos_Px_ALL.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  
  ### define persisters
  sub_seurat_obj2@meta.data$CD8_Persisters <- "Non-subsisters"
  sub_seurat_obj2@meta.data$CD8_Persisters[which(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones)] <- "Subsisters"
  
  ### perform clustering
  sub_seurat_obj2 <- FindNeighbors(sub_seurat_obj2, dims = 1:10)
  sub_seurat_obj2 <- FindClusters(sub_seurat_obj2, resolution = 0.5)
  
  ### save the clustering result to meta.data
  sub_seurat_obj2@meta.data$clusters <- Idents(sub_seurat_obj2)
  
  ### UMAP with subsister info
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "CD8_Persisters",
               pt.size = 5,
               cols = c("Non-subsisters" = "lightgray", "Subsisters" = "red"),
               order = c("Subsisters", "Non-subsisters")) +
    ggtitle("") +
    labs(color="Is_Persistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_CARpos_Subsister_ALL.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### UMAP with subsister info split by each patient
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "CD8_Persisters", split.by = "px",
               pt.size = 5, ncol = 3,
               cols = c("Non-subsisters" = "lightgray", "Subsisters" = "red"),
               order = c("Subsisters", "Non-subsisters")) +
    ggtitle("") +
    labs(color="Is_Persistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_CARpos_Subsister.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### load DE genes
  gmp_carpos_cd8_de_genes <- read.xlsx2(gmp_carpos_cd8_de_path, sheetIndex = 1,
                                        stringsAsFactors = FALSE, check.names = FALSE)
  rownames(gmp_carpos_cd8_de_genes) <- gmp_carpos_cd8_de_genes$Gene
  
  ### See gene expressions on UMAP with the top 1:9 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[1:9], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(1).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 10:18 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[10:18], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(2).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 19:27 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[19:27], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(3).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 28:36 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[28:36], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(4).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 37:45 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[37:45], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(5).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### set column for cluster 0,1,4,8 and the others
  sub_seurat_obj2@meta.data$pi_cluster0 <- as.character(sub_seurat_obj2@meta.data$clusters)
  sub_seurat_obj2@meta.data$pi_cluster0[which(sub_seurat_obj2@meta.data$pi_cluster0 != "0")] <- "Others"
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$pi_cluster0)
  
  ### DE analysis
  de_result <- FindMarkers(sub_seurat_obj2,
                           ident.1 = "0",
                           ident.2 = "Others",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### See gene expressions on UMAP with the top 1:9 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[1:9], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0_GEXP_9_DE_Genes(1).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 10:18 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[10:18], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0_GEXP_9_DE_Genes(2).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 19:27 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[19:27], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0_GEXP_9_DE_Genes(3).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 28:36 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[28:36], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0_GEXP_9_DE_Genes(4).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 37:45 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[37:45], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0_GEXP_9_DE_Genes(5).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### get entrez ids for the genes
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[which(de_result$p_val_adj < 0.01)],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("Pathways_with_DE_Genes_PI_Cluster0"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir2)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathways_with_DE_Genes_PI_Cluster0"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir2)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir2, "GO_Pathways_with_DE_genes_PI_Cluster0.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, "KEGG_Pathways_with_DE_genes_PI_Cluster0.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  # CD3D - t cells
  # CD4 - cd4 t cells
  # CD8A - cd8 t cells
  # KLRF1 - NK cells
  # CD14 - monocytes/macrophages
  # FCGR3A - NK cells
  # FOXP3 - t regs
  # STAT4 - th1 cells
  # STAT6 - th2 cells
  # STAT3 - th17 cells
  # CD44 - memory
  # SELL(CD62L) - central memory cells
  # IL7R - central memory cells
  # CCR7 - central memory
  # CD27 - effector cells
  # KLRG1 - SLEC (short-lived effector cells)
  # CD69 - tissue resident memory
  # CD127 - MPEC (memory precursor effector cells)
  # CD45RO - effector memory t cells
  
  ### 1
  # CD4 - cd4 t cells
  # CD8A - cd8 t cells
  # KLRF1 - NK cells
  # CD14 - monocytes/macrophages
  # FCGR3A - NK cells
  # FOXP3 - t regs
  # STAT4 - th1 cells
  # STAT6 - th2 cells
  # STAT3 - th17 cells
  gene_set <- c("CD4", "CD8A", "KLRF1",
                "CD14", "FCGR3A", "FOXP3",
                "STAT4", "STAT6", "STAT3")
  names(gene_set) <- c("CD4", "CD8", "NK",
                       "Mono/Mcrphg", "NK", "Treg",
                       "Th1", "Th2", "Th17")
  p <- FeaturePlot(sub_seurat_obj2, features = gene_set,
              cols = c("lightgray", "red"))
  for(i in 1:length(gene_set)) {
    p[[i]]$labels$title <- paste(gene_set[i], "-", names(gene_set)[i])
  }
  ggsave(paste0(outputDir2, "CARpos_9_Interesting_Genes(1).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### 2
  # CD4 - cd4 t cells
  # CD8A - cd8 t cells
  # CD44 - memory
  # SELL(CD62L) - central memory cells
  # IL7R - central memory cells
  # CCR7 - central memory
  # CD27 - effector cells
  # KLRG1 - SLEC (short-lived effector cells)
  # CD69 - tissue resident memory
  gene_set <- c("CD4", "CD8A", "CD44",
                "SELL", "IL7R", "CCR7",
                "CD27", "KLRG1", "CD69")
  names(gene_set) <- c("CD4", "CD8", "Memory",
                       "(CD62L) Central Memory", "Centeral Memory", "Central Memory",
                       "Effector", "SLEC", "Tissue Resident Memory")
  p <- FeaturePlot(sub_seurat_obj2, features = gene_set,
                   cols = c("lightgray", "red"), max.cutoff = 3)
  for(i in 1:length(gene_set)) {
    p[[i]]$labels$title <- paste(gene_set[i], "-", names(gene_set)[i])
  }
  ggsave(paste0(outputDir2, "CARpos_9_Interesting_Genes(2).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  #
  ### 7. CAR+ (>0, >1, >2, etc.) numbers for each patient between "From Sorting" and "From scRNA-Seq"
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/7/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### this sould be TRUE for further analysis
  print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
  
  ### change the CAR exp threshold
  n <- 9
  car_num_mat <- matrix(0, length(unique(Seurat_Obj@meta.data$px))+1, n)
  rownames(car_num_mat) <- c(unique(Seurat_Obj@meta.data$px), "Total")
  colnames(car_num_mat) <- paste0("CAR_TRANSCRIPT>", 0:(n-1))
  gmp_car_num_mat <- car_num_mat
  for(i in 1:n) {
    ### get the indicies for CAR+ cells based on the threshold
    car_idx <- which(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",] > (i-1))
    gmp_car_idx <- intersect(which(Seurat_Obj@meta.data$time2 == "GMP"), car_idx)
    
    ### get CAR+ numbers for each patient
    for(px in unique(Seurat_Obj@meta.data$px)) {
      ### indicies for the given patient
      px_idx <- which(Seurat_Obj@meta.data$px == px)
      
      ### get CAR+ numbers
      car_num_mat[px, i] <- length(intersect(car_idx, px_idx))
      gmp_car_num_mat[px, i] <- length(intersect(gmp_car_idx, px_idx))
    }
    
    ### get total numbers
    car_num_mat["Total", i] <- sum(car_num_mat[1:(nrow(car_num_mat)-1),i])
    gmp_car_num_mat["Total", i] <- sum(gmp_car_num_mat[1:(nrow(gmp_car_num_mat)-1),i])
  }
  
  ### write out the result
  write.xlsx2(data.frame(Px=rownames(car_num_mat),
                         car_num_mat,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CAR_Numbers_per_px.xlsx"),
              sheetName = "CAR_Numbers_per_px", row.names = FALSE)
  
  ### write out the result
  write.xlsx2(data.frame(Px=rownames(gmp_car_num_mat),
                         gmp_car_num_mat,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CAR_Numbers_per_px.xlsx"),
              sheetName = "GMP_CAR_Numbers_per_px", row.names = FALSE)
  
  ### CAR+ percentages
  car_pcnt_mat <- car_num_mat
  gmp_car_pcnt_mat <- gmp_car_num_mat
  for(px in unique(Seurat_Obj@meta.data$px)) {
    car_pcnt_mat[px,] <- round(car_pcnt_mat[px,] * 100 / length(which(Seurat_Obj@meta.data$px == px)), 2)
    gmp_car_pcnt_mat[px,] <- round(gmp_car_pcnt_mat[px,] * 100 / length(intersect(which(Seurat_Obj@meta.data$time2 == "GMP"),
                                                                                  which(Seurat_Obj@meta.data$px == px))), 2)
  }
  car_pcnt_mat["Total",] <- round(car_pcnt_mat["Total",] * 100 / nrow(Seurat_Obj@meta.data), 2)
  gmp_car_pcnt_mat["Total",] <- round(gmp_car_pcnt_mat["Total",] * 100 / length(which(Seurat_Obj@meta.data$time2 == "GMP")), 2)
  
  ### write out the result
  write.xlsx2(data.frame(Px=rownames(car_pcnt_mat),
                         car_pcnt_mat,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CAR_Numbers_per_px.xlsx"),
              sheetName = "CAR_Percentages_per_px", row.names = FALSE, append = TRUE)
  
  ### write out the result
  write.xlsx2(data.frame(Px=rownames(gmp_car_pcnt_mat),
                         gmp_car_pcnt_mat,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CAR_Numbers_per_px.xlsx"),
              sheetName = "GMP_CAR_Percentages_per_px", row.names = FALSE, append = TRUE)
  
  ### draw density plot for all the time points
  plot_df <- data.frame(CAR_EXP=Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",rownames(Seurat_Obj@meta.data)],
                        Px=Seurat_Obj@meta.data$px,
                        stringsAsFactors = FALSE, check.names = FALSE)
  p <- ggplot(plot_df, aes_string(x="CAR_EXP", col="Px")) +
    geom_density(size = 2) +
    labs(color="") +
    ylab("") +
    scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
    theme_classic(base_size = 30)
  ggsave(paste0(outputDir2, "CAR_Numbers_Density.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### draw density plot for GMP only
  plot_df <- data.frame(GMP_CAR_EXP=Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",rownames(Seurat_Obj@meta.data)[which(Seurat_Obj@meta.data$time2 == "GMP")]],
                        Px=Seurat_Obj@meta.data$px[which(Seurat_Obj@meta.data$time2 == "GMP")],
                        stringsAsFactors = FALSE, check.names = FALSE)
  p <- ggplot(plot_df, aes_string(x="GMP_CAR_EXP", col="Px")) +
    geom_density(size = 2) +
    labs(color="") +
    ylab("") +
    scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
    theme_classic(base_size = 30)
  ggsave(paste0(outputDir2, "GMP_CAR_Numbers_Density.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  
  #
  ### 8. a) If sampling from GMP and sampling from an after infusion time point, how many matches do we see?
  #      Using GMP - clone size as a background to estimate a selection factor
  #      And show distribution shift
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/8/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  
  ### after infusion time points
  after_gmp_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6",
                             "Wk8", "3mo", "6mo", "9mo")
  
  ### make empty tables
  selection_table <- matrix(0, 1, 11)
  colnames(selection_table) <- c("Clone", "Patient", "Time", "GMP CD8 CAR+ #", "GMP Clone Size",
                                 "Given Time CD8 CAR+ #", "Given Time Expected Clone Size", "Given Time Clone size",
                                 "Selection %", "Odds Ratio", "P-value")
  selection_table2 <- matrix(0, 1, 10)
  colnames(selection_table2) <- c("Clone", "Patient", "GMP CD8 CAR+ #", "GMP Clone Size",
                                  "Post-Infusion CD8 CAR+ #", "Post-Infusion Expected Clone Size", "Post-Infusion Clone size",
                                  "Selection %", "Odds Ratio", "P-value")
  
  ### make data frame
  selection_table <- data.frame(selection_table,
                                stringsAsFactors = FALSE, check.names = FALSE)
  selection_table2 <- data.frame(selection_table2,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  
  ### for each lineage what's the selection factor?
  for(lin in persister_clones) {
    
    ### get px info
    px <- strsplit(lin, split = "_", fixed = TRUE)[[1]][1]
    
    ### get time points that the given lineage appeared
    tps <- unique(Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta == lin)])
    tps <- intersect(after_gmp_time_points, tps)
    
    ### get GMP CD8 CAR+ # of the patient
    gmp_cd8_carpos_numbers <- length(intersect(intersect(which(Seurat_Obj@meta.data$time2 == "GMP"),
                                                         which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8")),
                                               intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                         which(Seurat_Obj@meta.data$px == px))))
    
    ### get the clone size in GMP CD8 CAR+ of the patient
    gmp_clone_size <- length(intersect(intersect(which(Seurat_Obj@meta.data$time2 == "GMP"),
                                                 which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta == lin)),
                                       intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                 which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))))
    
    ### get info of all the post-infusion time points
    pi_cd8_carpos_numbers <- length(intersect(intersect(which(Seurat_Obj@meta.data$time2 %in% tps),
                                                        which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8")),
                                              intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                        which(Seurat_Obj@meta.data$px == px))))
    pi_clone_size <- length(intersect(intersect(which(Seurat_Obj@meta.data$time2 %in% tps),
                                                which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta == lin)),
                                      intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))))
    pi_expected_clone_size <- pi_cd8_carpos_numbers * gmp_clone_size / gmp_cd8_carpos_numbers
    pi_selection_pcnt <- round(pi_clone_size * 100 / pi_expected_clone_size, 2)
    
    ### calculate p-value
    ### Fisher's exact test
    ###
    ###           TP Clone   No TP (GMP) Clone
    ###          ----------------------------
    ###    Clone |   X              Y
    ### No Clone |   Z              W
    X <- pi_clone_size
    Y <- gmp_clone_size
    Z <- pi_cd8_carpos_numbers - pi_clone_size
    W <- gmp_cd8_carpos_numbers - gmp_clone_size
    
    pi_odds_ratio <- fisher.test(matrix(c(X, Z, Y, W), 2, 2), alternative = "greater")$estimate
    pi_p_value <- fisher.test(matrix(c(X, Z, Y, W), 2, 2), alternative = "greater")$p.value
    
    ### add new row to the table
    selection_table2 <- rbind(selection_table2, c(lin, px, gmp_cd8_carpos_numbers, gmp_clone_size,
                                                  pi_cd8_carpos_numbers, pi_expected_clone_size,
                                                  pi_clone_size, pi_selection_pcnt, pi_odds_ratio, pi_p_value))
    
    
    ### for each time point
    for(tp in tps) {
      ### number for the given time point
      given_tp_cd8_carpos_numbers <- length(intersect(intersect(which(Seurat_Obj@meta.data$time2 == tp),
                                                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8")),
                                                      intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                                which(Seurat_Obj@meta.data$px == px))))
      given_tp_clone_size <- length(intersect(intersect(which(Seurat_Obj@meta.data$time2 == tp),
                                                        which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta == lin)),
                                              intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                        which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))))
      
      ### additional info (selection %, odds ratio, p-value)
      expected_clone_size <- given_tp_cd8_carpos_numbers * gmp_clone_size / gmp_cd8_carpos_numbers
      selection_pcnt <- round(given_tp_clone_size * 100 / expected_clone_size, 2)
      
      ### calculate p-value
      ### Fisher's exact test
      ###
      ###           TP Clone   No TP (GMP) Clone
      ###          ----------------------------
      ###    Clone |   X              Y
      ### No Clone |   Z              W
      X <- given_tp_clone_size
      Y <- gmp_clone_size
      Z <- given_tp_cd8_carpos_numbers - given_tp_clone_size
      W <- gmp_cd8_carpos_numbers - gmp_clone_size
      
      odds_ratio <- fisher.test(matrix(c(X, Z, Y, W), 2, 2), alternative = "greater")$estimate
      p_value <- fisher.test(matrix(c(X, Z, Y, W), 2, 2), alternative = "greater")$p.value
      
      ### add new row to the table
      selection_table <- rbind(selection_table, c(lin, px, tp, gmp_cd8_carpos_numbers, gmp_clone_size,
                                                  given_tp_cd8_carpos_numbers, expected_clone_size,
                                                  given_tp_clone_size, selection_pcnt, odds_ratio, p_value))
    }
    
  }
  
  ### remove the fake (first) row
  selection_table <- selection_table[-1,]
  selection_table2 <- selection_table2[-1,]
  
  ### numerize some columns
  selection_table[,4:11] <- sapply(selection_table[,4:11], as.numeric)
  selection_table2[,3:10] <- sapply(selection_table2[,3:10], as.numeric)
  
  ### multiple testing correction
  selection_table$FDR <- p.adjust(selection_table$`P-value`, method = "BH")
  selection_table2$FDR <- p.adjust(selection_table2$`P-value`, method = "BH")
  
  ### order them based on FDR
  selection_table <- selection_table[order(selection_table$FDR),]
  selection_table2 <- selection_table2[order(selection_table2$FDR),]
  
  ### write out the result
  write.xlsx2(selection_table,
              file = paste0(outputDir2, "/Subsister_Selection_Fator_Foreach_Time.xlsx"),
              sheetName = "Subsister_Selection_Time", row.names = FALSE)
  write.xlsx2(selection_table2,
              file = paste0(outputDir2, "/Subsister_Selection_Fator_All_PostInfusion.xlsx"),
              sheetName = "Subsister_Selection_All_Post_Infusion", row.names = FALSE)
  
  
  #
  ### clone based not cell based
  #
  
  ### target patients
  target_pxs <- unique(Seurat_Obj@meta.data$px[which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones)])
  
  ### make an empty table
  selection_table3 <- matrix(0, 1, 10)
  colnames(selection_table3) <- c("Patient", "Time", "GMP CAR+ Clone #", "GMP Subsister Clone #",
                                  "Given Time CAR+ Clone #", "Given Time Expected Subsister Clone #",
                                  "Given Time Observed Subsister Clone #",
                                  "Selection %", "Odds Ratio", "P-value")
  
  ### for each patient
  for(patient in target_pxs) {
    
    ### get time points of the patients
    tps <- unique(Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$px == patient)])
    tps <- intersect(after_gmp_time_points, tps)
    
    ### get info of GMP
    gmp_carpos_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[intersect(intersect(which(Seurat_Obj@meta.data$px == patient),
                                                                                                                which(Seurat_Obj@meta.data$time2 == "GMP")),
                                                                                                      which(Seurat_Obj@meta.data$CAR == "CARpos"))])
    ### remove NA in the clones
    gmp_carpos_clones <- gmp_carpos_clones[which(!is.na(gmp_carpos_clones))]
    
    gmp_carpos_clone_num <- length(gmp_carpos_clones)
    gmp_subsister_clone_num <- length(intersect(gmp_carpos_clones,
                                                persister_clones))
    
    ### for each time point
    for(tp in tps) {
      ### number for the given time point
      given_time_carpos_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[intersect(intersect(which(Seurat_Obj@meta.data$px == patient),
                                                                                                                         which(Seurat_Obj@meta.data$time2 == tp)),
                                                                                                               which(Seurat_Obj@meta.data$CAR == "CARpos"))])
      given_time_carpos_clone_num <- length(given_time_carpos_clones)
      given_time_subsister_clone_num <- length(intersect(given_time_carpos_clones,
                                                         persister_clones))
      
      ### additional info (selection %, odds ratio, p-value)
      expected_clone_num <- given_time_carpos_clone_num * gmp_subsister_clone_num / gmp_carpos_clone_num
      selection_pcnt <- round(given_time_subsister_clone_num * 100 / expected_clone_num, 2)
      
      ### calculate p-value
      ### Fisher's exact test
      ###
      ###           TP Clone   No TP (GMP) Clone
      ###          ----------------------------
      ###    Clone |   X              Y
      ### No Clone |   Z              W
      X <- given_time_subsister_clone_num
      Y <- gmp_subsister_clone_num
      Z <- given_time_carpos_clone_num - given_time_subsister_clone_num
      W <- gmp_carpos_clone_num - gmp_subsister_clone_num
      
      odds_ratio <- fisher.test(matrix(c(X, Z, Y, W), 2, 2), alternative = "greater")$estimate
      p_value <- fisher.test(matrix(c(X, Z, Y, W), 2, 2), alternative = "greater")$p.value
      
      ### add new row to the table
      selection_table3 <- rbind(selection_table3, c(patient, tp, gmp_carpos_clone_num, gmp_subsister_clone_num,
                                                    given_time_carpos_clone_num, expected_clone_num,
                                                    given_time_subsister_clone_num, selection_pcnt, odds_ratio, p_value))
    }
    
  }
  
  ### remove the fake (first) row
  selection_table3 <- selection_table3[-1,]
  
  ### make data frame
  selection_table3 <- data.frame(selection_table3,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  
  ### numerize some columns
  selection_table3[,3:10] <- sapply(selection_table3[,3:10], as.numeric)
  
  ### multiple testing correction
  selection_table3$FDR <- p.adjust(selection_table3$`P-value`, method = "BH")
  
  ### order them based on FDR
  selection_table3 <- selection_table3[order(selection_table3$FDR),]
  
  ### write out the result
  write.xlsx2(selection_table3,
              file = paste0(outputDir2, "/Subsister_Selection_Fator_Clone.xlsx"),
              sheetName = "Subsister_Selection_Fator_Clone", row.names = FALSE)
  
  
  
  #
  ### 8. b) Use the % of CD4/CD8 subsisters vs Total CD4/CD8 in GMP to estimate how many CD8 subsisters were infused
  #         Do not use GMP-redo here
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/8/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  # ### the percentage of CD4/CD8 in the CAR+ GMP product
  # gmp_product_cd4_pcnt <- rep(0, length(unique(Seurat_Obj@meta.data$px)))
  # names(gmp_product_cd4_pcnt) <- unique(Seurat_Obj@meta.data$px)
  # gmp_product_cd8_pcnt <- rep(0, length(unique(Seurat_Obj@meta.data$px)))
  # names(gmp_product_cd8_pcnt) <- unique(Seurat_Obj@meta.data$px)
  
  # gmp_product_cd4_pcnt <- c(63.7, 79.1, 49.4, 57.5, 61.7, 39.2, 89.5, 62.0, 53.2, NA, NA, NA, NA, NA, NA, NA)
  # gmp_product_cd8_pcnt <- c(14,5, 10.2, 48.1, 22.4, 36.4, 59.8, 8.1, 7.5, 29.1, NA, NA, NA, NA, NA, NA, NA)
  # 
  # ### cd3 percentage of the GMP product
  # gmp_cd3_pcnt <- c(0.962, 0.999, 0.995, 0.996, 0.990, 0.999, 0.993, 0.982, 0.983,
  #                   0.999, 0.998, 0.993, 0.996, 0.999, 0.997, 0.996, 0.995, 0.996, 99.810)
  
  ### the info we've got
  Real_GMP_CARpos_pcnt <- c(41.90, 37.00, 41.90, 59.10, 33.30, 23.20, 26.70, 14.70,
                            58.30, 33.60, 42.10, 32.00, 27.40, 32.60, NA, NA)
  Real_GMP_CAR_CD4_pcnt <- c(60.60, 93.60, 60.60, 75.50, 77.50, 61.40, 95.40, 85.00,
                             66.90, 52.10, 63.50, 58.60, 70.70, 71.20, NA, NA)
  Real_GMP_CAR_CD8_pcnt <- c(37.20, 5.50, 37.20, 23.60, 21.00, 36.80, 2.90, 2.40,
                             30.80, 46.70, 35.00, 39.90, 28.00, 27.70, NA, NA)
  names(Real_GMP_CARpos_pcnt) <- unique(Seurat_Obj@meta.data$px)
  names(Real_GMP_CAR_CD4_pcnt) <- unique(Seurat_Obj@meta.data$px)
  names(Real_GMP_CAR_CD8_pcnt) <- unique(Seurat_Obj@meta.data$px)
  
  ### calculate some
  Real_GMP_CAR_CD4_pcnt2 <- Real_GMP_CARpos_pcnt * Real_GMP_CAR_CD4_pcnt / 10000
  Real_GMP_CAR_CD8_pcnt2 <- Real_GMP_CARpos_pcnt * Real_GMP_CAR_CD8_pcnt / 10000
  
  ### the cell count of the GMP product
  gmp_cell_cnt <- c(1887000000, 3920000000, 2620000000, 4000000000, 3870000000,
                    4030000000, 3620000000, 2820000000, 2700000000, 4060000000,
                    4180000000, 2930000000, 5150000000, 3200000000, 3010000000, 4590000000)
  names(gmp_cell_cnt) <- unique(Seurat_Obj@meta.data$px)
  
  ### dose info
  dose_levels <- c(1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000,
                   3000000, 3000000, 3000000, 3000000, 3000000, 3000000, 3000000, 3000000, 3000000)
  
  ### result_table
  result_table <- matrix(0, length(unique(Seurat_Obj@meta.data$px)), 16)
  rownames(result_table) <- unique(Seurat_Obj@meta.data$px)
  colnames(result_table) <- c("GMP CAR+ Cell %", "CD4 % of CAR+", "CD8 % of CAR+", "GMP Total Cell #",
                              "10x GMP CAR+ CD4 Cell #", "10x GMP CAR+ CD8 Cell #",
                              "10x GMP CAR+ CD4 Subsisters #", "10x GMP CAR+ CD8 Subsisters #",
                              "GMP CD4 Subsister %", "GMP CD8 Subsister %",
                              "Estimated GMP CD4 Subsisters #", "Estimated GMP CD8 Subsisters #",
                              "Dose (per Kg)",
                              "Adjusted Dose (per Kg)",
                              "Estimated Injected GMP CD4 Subsisters # (per kg)",
                              "Estimated Injected GMP CD8 Subsisters # (per kg)")
  
  ### data frame
  result_table <- data.frame(result_table, stringsAsFactors = FALSE, check.names = FALSE)
  result_table$`GMP CAR+ Cell %` <- Real_GMP_CARpos_pcnt
  result_table$`CD4 % of CAR+` <- Real_GMP_CAR_CD4_pcnt
  result_table$`CD8 % of CAR+` <- Real_GMP_CAR_CD8_pcnt
  result_table$`GMP Total Cell #` <- gmp_cell_cnt
  result_table$`Dose (per Kg)` <- dose_levels
  result_table$`Adjusted Dose (per Kg)` <- round(dose_levels * 100 / Real_GMP_CARpos_pcnt, 0)
  result_table$`Estimated Injected GMP CD4 Subsisters # (per kg)` <- NA
  result_table$`Estimated Injected GMP CD8 Subsisters # (per kg)` <- NA
  
  ### for each patient
  for(px in unique(Seurat_Obj@meta.data$px)) {
    ### some pre-computed indicies (exclude GMP-redo)
    gmp_carpos_cd4_idx <- intersect(intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                              which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4")),
                                    intersect(which(Seurat_Obj@meta.data$time == "GMP"),
                                              which(Seurat_Obj@meta.data$px == px)))
    gmp_carpos_cd8_idx <- intersect(intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                              which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8")),
                                    intersect(which(Seurat_Obj@meta.data$time == "GMP"),
                                              which(Seurat_Obj@meta.data$px == px)))
    
    ### CD4/CD8 GMP CARpos subsister cells
    ### the numbers can different because we are only considering GMP excluding GMP-redo
    sc_gmp_carpos_cd4_subsister_cell_num <- length(intersect(which(Seurat_Obj@meta.data$ALL_GMP_CARpos_Persister == "YES"),
                                                             gmp_carpos_cd4_idx))
    sc_gmp_carpos_cd8_subsister_cell_num <- length(intersect(which(Seurat_Obj@meta.data$ALL_GMP_CARpos_Persister == "YES"),
                                                             gmp_carpos_cd8_idx))
    
    ### single cell GMP CAR+ CD4/CD8 total vs CD4/CD8 subsister ratio
    sc_cd4_subsister_ratio <- sc_gmp_carpos_cd4_subsister_cell_num / length(gmp_carpos_cd4_idx)
    sc_cd8_subsister_ratio <- sc_gmp_carpos_cd8_subsister_cell_num / length(gmp_carpos_cd8_idx)
    
    ### fill the table
    result_table[px, "10x GMP CAR+ CD4 Cell #"] <- length(gmp_carpos_cd4_idx)
    result_table[px, "10x GMP CAR+ CD8 Cell #"] <- length(gmp_carpos_cd8_idx)
    result_table[px, "10x GMP CAR+ CD4 Subsisters #"] <- sc_gmp_carpos_cd4_subsister_cell_num
    result_table[px, "10x GMP CAR+ CD8 Subsisters #"] <- sc_gmp_carpos_cd8_subsister_cell_num
    result_table[px, "Estimated GMP CD4 Subsisters #"] <- round(gmp_cell_cnt[px] * Real_GMP_CAR_CD4_pcnt2[px] * sc_cd4_subsister_ratio, 2)
    result_table[px, "Estimated GMP CD8 Subsisters #"] <- round(gmp_cell_cnt[px] * Real_GMP_CAR_CD8_pcnt2[px] * sc_cd8_subsister_ratio, 2)
    result_table[px, "GMP CD4 Subsister %"] <- round(Real_GMP_CAR_CD4_pcnt2[px] * sc_cd4_subsister_ratio * 100, 5)
    result_table[px, "GMP CD8 Subsister %"] <- round(Real_GMP_CAR_CD8_pcnt2[px] * sc_cd8_subsister_ratio * 100, 5)
  }
  
  ### calculate some more info
  result_table[,"Estimated Injected GMP CD4 Subsisters # (per kg)"] <- round(result_table[,"Adjusted Dose (per Kg)"] * result_table[,"GMP CD4 Subsister %"] / 100, 0)
  result_table[,"Estimated Injected GMP CD8 Subsisters # (per kg)"] <- round(result_table[,"Adjusted Dose (per Kg)"] * result_table[,"GMP CD8 Subsister %"] / 100, 0)
  
  ### write out the result_table
  write.xlsx2(data.frame(Patient=rownames(result_table),
                         result_table,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/Estimated_CD4_CD8_Subsister_#_In_GMP.xlsx"),
              sheetName = "Estimated_CD4_CD8_Subsister_#_GMP", row.names = FALSE)
  
  ### bar graph
  plot_df <- data.frame(Cell_Num=c(result_table$`Estimated Injected GMP CD4 Subsisters # (per kg)`,
                                   result_table$`Estimated Injected GMP CD8 Subsisters # (per kg)`),
                        Patient=c(rownames(result_table), rownames(result_table)),
                        Cell_Type=c(rep("CD4", nrow(result_table)), rep("CD8", nrow(result_table))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df <- plot_df[which(!is.na(plot_df$Cell_Num)),]
  plot_df$Cell_Num2 <- as.character(plot_df$Cell_Num)
  plot_df$Cell_Num2[which(plot_df$Cell_Num2 == "0")] <- ""
  plot_df <- plot_df[-which(plot_df$Patient %in% c("SJCAR19-00", "SJCAR19-01", "SJCAR19-03", "SJCAR19-12")),]
  temp <- plot_df[order(-plot_df$Cell_Num),]
  temp <- temp[which(temp$Cell_Type == "CD8"),]
  plot_df$Patient <- factor(plot_df$Patient, levels = unique(temp$Patient))
  plot_df$Cell_Num[which(plot_df$Cell_Num > 10000)] <- 10000
  p <- ggplot(data=plot_df, aes_string(x="Patient", y="Cell_Num", fill="Cell_Type", label="Cell_Num2")) +
    geom_bar(position = position_dodge(), stat = "identity") +
    ggtitle("Infused GMP Subsisters # (per kg)") +
    geom_text(size = 7, position = position_dodge(1)) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5, size = 20))
  ggsave(file = paste0(outputDir2, "Estimated_CD4_CD8_Subsister_#_In_GMP.png"), plot = p, width = 20, height = 10, dpi = 400)
  
  temp_df <- plot_df[which(plot_df$Cell_Type == "CD4"),]
  temp_df <- temp_df[order(-temp_df$Cell_Num),]
  temp_df <- temp_df[which(temp_df$Cell_Num > 0),]
  temp_df$Patient <- factor(temp_df$Patient, levels = unique(temp_df$Patient))
  p <- ggplot(data=temp_df, aes_string(x="Patient", y="Cell_Num", label="Cell_Num2")) +
    geom_bar(position = position_dodge(), stat = "identity", fill="#F8766D") +
    ggtitle("Infused GMP CD4 Subsisters # (per kg)") +
    geom_text(size = 6, position = position_dodge(1)) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5, size = 20))
  ggsave(file = paste0(outputDir2, "Estimated_CD4_Subsister_#_In_GMP.png"), plot = p, width = 25, height = 10, dpi = 400)
  
  temp_df <- plot_df[which(plot_df$Cell_Type == "CD8"),]
  temp_df <- temp_df[order(-temp_df$Cell_Num),]
  temp_df <- temp_df[which(temp_df$Cell_Num > 0),]
  temp_df$Patient <- factor(temp_df$Patient, levels = unique(temp_df$Patient))
  p <- ggplot(data=temp_df, aes_string(x="Patient", y="Cell_Num", label="Cell_Num2")) +
    geom_bar(position = position_dodge(), stat = "identity", fill="#00BFC4") +
    ggtitle("Infused GMP CD8 Subsisters # (per kg)") +
    geom_text(size = 6, position = position_dodge(1)) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5, size = 20))
  ggsave(file = paste0(outputDir2, "Estimated_CD8_Subsister_#_In_GMP.png"), plot = p, width = 25, height = 10, dpi = 400)
  
  
  #
  ### 8. c) Use the % of Cluster01 vs Total CD4/CD8 in GMP to estimate how many subsister-like cells were infused
  #         Do not use GMP-redo here
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/8/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### cluster01 clones
  cluster01_clones <- unique(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta[which(sub_seurat_obj2@meta.data$clusters %in% c(0, 1))])
  
  ### remove NA from cluster01 clones
  cluster01_clones <- cluster01_clones[which(!is.na(cluster01_clones))]
  
  ### add cluster01 column to the original seurat object
  Seurat_Obj@meta.data$Cluster01 <- "NO"
  Seurat_Obj@meta.data$Cluster01[which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% cluster01_clones)] <- "YES"
  
  ### result_table
  result_table2 <- matrix(0, length(unique(Seurat_Obj@meta.data$px)), 16)
  rownames(result_table2) <- unique(Seurat_Obj@meta.data$px)
  colnames(result_table2) <- c("GMP CAR+ Cell %", "CD4 % of CAR+", "CD8 % of CAR+", "GMP Total Cell #",
                              "10x GMP CAR+ CD4 Cell #", "10x GMP CAR+ CD8 Cell #",
                              "10x GMP CAR+ CD4 Subsisters #", "10x GMP CAR+ CD8 Subsisters #",
                              "GMP CD4 Subsister %", "GMP CD8 Subsister %",
                              "Estimated GMP CD4 Subsisters #", "Estimated GMP CD8 Subsisters #",
                              "Dose (per Kg)",
                              "Adjusted Dose (per Kg)",
                              "Estimated Injected GMP CD4 Subsisters # (per kg)",
                              "Estimated Injected GMP CD8 Subsisters # (per kg)")
  
  ### data frame
  result_table2 <- data.frame(result_table2, stringsAsFactors = FALSE, check.names = FALSE)
  result_table2$`GMP CAR+ Cell %` <- Real_GMP_CARpos_pcnt
  result_table2$`CD4 % of CAR+` <- Real_GMP_CAR_CD4_pcnt
  result_table2$`CD8 % of CAR+` <- Real_GMP_CAR_CD8_pcnt
  result_table2$`GMP Total Cell #` <- gmp_cell_cnt
  result_table2$`Dose (per Kg)` <- dose_levels
  result_table2$`Adjusted Dose (per Kg)` <- round(dose_levels * 100 / Real_GMP_CARpos_pcnt, 0)
  result_table2$`Estimated Injected GMP CD4 Subsisters # (per kg)` <- NA
  result_table2$`Estimated Injected GMP CD8 Subsisters # (per kg)` <- NA
  
  ### for each patient
  for(px in unique(Seurat_Obj@meta.data$px)) {
    ### some pre-computed indicies (exclude GMP-redo)
    gmp_carpos_cd4_idx <- intersect(intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                              which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4")),
                                    intersect(which(Seurat_Obj@meta.data$time == "GMP"),
                                              which(Seurat_Obj@meta.data$px == px)))
    gmp_carpos_cd8_idx <- intersect(intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                              which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8")),
                                    intersect(which(Seurat_Obj@meta.data$time == "GMP"),
                                              which(Seurat_Obj@meta.data$px == px)))
    
    ### CD4/CD8 GMP CARpos subsister cells
    ### the numbers can different because we are only considering GMP excluding GMP-redo
    sc_gmp_carpos_cd4_subsister_cell_num <- length(intersect(which(Seurat_Obj@meta.data$Cluster01 == "YES"),
                                                             gmp_carpos_cd4_idx))
    sc_gmp_carpos_cd8_subsister_cell_num <- length(intersect(which(Seurat_Obj@meta.data$Cluster01 == "YES"),
                                                             gmp_carpos_cd8_idx))
    
    ### single cell GMP CAR+ CD4/CD8 total vs CD4/CD8 subsister ratio
    sc_cd4_subsister_ratio <- sc_gmp_carpos_cd4_subsister_cell_num / length(gmp_carpos_cd4_idx)
    sc_cd8_subsister_ratio <- sc_gmp_carpos_cd8_subsister_cell_num / length(gmp_carpos_cd8_idx)
    
    ### fill the table
    result_table2[px, "10x GMP CAR+ CD4 Cell #"] <- length(gmp_carpos_cd4_idx)
    result_table2[px, "10x GMP CAR+ CD8 Cell #"] <- length(gmp_carpos_cd8_idx)
    result_table2[px, "10x GMP CAR+ CD4 Subsisters #"] <- sc_gmp_carpos_cd4_subsister_cell_num
    result_table2[px, "10x GMP CAR+ CD8 Subsisters #"] <- sc_gmp_carpos_cd8_subsister_cell_num
    result_table2[px, "Estimated GMP CD4 Subsisters #"] <- round(gmp_cell_cnt[px] * Real_GMP_CAR_CD4_pcnt2[px] * sc_cd4_subsister_ratio, 2)
    result_table2[px, "Estimated GMP CD8 Subsisters #"] <- round(gmp_cell_cnt[px] * Real_GMP_CAR_CD8_pcnt2[px] * sc_cd8_subsister_ratio, 2)
    result_table2[px, "GMP CD4 Subsister %"] <- round(Real_GMP_CAR_CD4_pcnt2[px] * sc_cd4_subsister_ratio * 100, 5)
    result_table2[px, "GMP CD8 Subsister %"] <- round(Real_GMP_CAR_CD8_pcnt2[px] * sc_cd8_subsister_ratio * 100, 5)
  }
  
  ### calculate some more info
  result_table2[,"Estimated Injected GMP CD4 Subsisters # (per kg)"] <- round(result_table2[,"Adjusted Dose (per Kg)"] * result_table2[,"GMP CD4 Subsister %"] / 100, 0)
  result_table2[,"Estimated Injected GMP CD8 Subsisters # (per kg)"] <- round(result_table2[,"Adjusted Dose (per Kg)"] * result_table2[,"GMP CD8 Subsister %"] / 100, 0)
  
  ### write out the result_table
  write.xlsx2(data.frame(Patient=rownames(result_table2),
                         result_table2,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/Estimated_CD4_CD8_Subsister(Cluster01)_#_In_GMP.xlsx"),
              sheetName = "Estimated_CD4_CD8_Subsister_#_GMP", row.names = FALSE)
  
  ### bar graph
  plot_df <- data.frame(Cell_Num=c(result_table2$`Estimated Injected GMP CD4 Subsisters # (per kg)`,
                                   result_table2$`Estimated Injected GMP CD8 Subsisters # (per kg)`),
                        Patient=c(rownames(result_table2), rownames(result_table2)),
                        Cell_Type=c(rep("CD4", nrow(result_table2)), rep("CD8", nrow(result_table2))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df <- plot_df[which(!is.na(plot_df$Cell_Num)),]
  plot_df$Cell_Num2 <- as.character(plot_df$Cell_Num)
  plot_df$Cell_Num2[which(plot_df$Cell_Num2 == "0")] <- ""
  temp <- plot_df[order(-plot_df$Cell_Num),]
  temp <- temp[which(temp$Cell_Type == "CD8"),]
  plot_df$Patient <- factor(plot_df$Patient, levels = unique(temp$Patient))
  plot_df$Cell_Num[which(plot_df$Cell_Num > 200000)] <- 200000
  p <- ggplot(data=plot_df, aes_string(x="Patient", y="Cell_Num", fill="Cell_Type", label="Cell_Num2")) +
    geom_bar(position = position_dodge(), stat = "identity") +
    ggtitle("Infused GMP Subsisters (Cluster01) # (per kg)") +
    geom_text(size = 6, position = position_dodge(1)) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5, size = 20))
  ggsave(file = paste0(outputDir2, "Estimated_CD4_CD8_Subsister(Cluster01)_#_In_GMP.png"), plot = p, width = 25, height = 10, dpi = 400)
  
  temp_df <- plot_df[which(plot_df$Cell_Type == "CD4"),]
  temp_df <- temp_df[order(-temp_df$Cell_Num),]
  temp_df <- temp_df[which(temp_df$Cell_Num > 0),]
  temp_df$Patient <- factor(temp_df$Patient, levels = unique(temp_df$Patient))
  p <- ggplot(data=temp_df, aes_string(x="Patient", y="Cell_Num", label="Cell_Num2")) +
    geom_bar(position = position_dodge(), stat = "identity", fill="#F8766D") +
    ggtitle("Infused GMP CD4 Subsisters (Cluster01) # (per kg)") +
    geom_text(size = 6, position = position_dodge(1)) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5, size = 20))
  ggsave(file = paste0(outputDir2, "Estimated_CD4_Subsister(Cluster01)_#_In_GMP.png"), plot = p, width = 25, height = 10, dpi = 400)
  
  temp_df <- plot_df[which(plot_df$Cell_Type == "CD8"),]
  temp_df <- temp_df[order(-temp_df$Cell_Num),]
  temp_df <- temp_df[which(temp_df$Cell_Num > 0),]
  temp_df$Patient <- factor(temp_df$Patient, levels = unique(temp_df$Patient))
  p <- ggplot(data=temp_df, aes_string(x="Patient", y="Cell_Num", label="Cell_Num2")) +
    geom_bar(position = position_dodge(), stat = "identity", fill="#00BFC4") +
    ggtitle("Infused GMP CD8 Subsisters (Cluster01) # (per kg)") +
    geom_text(size = 6, position = position_dodge(1)) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5, size = 20))
  ggsave(file = paste0(outputDir2, "Estimated_CD8_Subsister(Cluster01)_#_In_GMP.png"), plot = p, width = 25, height = 10, dpi = 400)
  
  
  #
  ### 8. d) Correlation between [subsister infusion amount (per kg)] and [PeakCAR]
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/8/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### PeakCAR info
  peakcar <- c(4215, 6178, 28867, 33667, 16709)
  names(peakcar) <- c("SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05", "SJCAR19-06")
  
  ### make the plot data frame
  plot_df <- data.frame(X=c(result_table[names(peakcar),"Estimated Injected GMP CD4 Subsisters # (per kg)"],
                            result_table[names(peakcar), "Estimated Injected GMP CD8 Subsisters # (per kg)"],
                            result_table2[names(peakcar),"Estimated Injected GMP CD4 Subsisters # (per kg)"],
                            result_table2[names(peakcar), "Estimated Injected GMP CD8 Subsisters # (per kg)"]),
                        Y=rep(peakcar, 4),
                        Group=c(rep("CD4 Subsisters", length(peakcar)),
                                rep("CD8 Subsisters", length(peakcar)),
                                rep("CD4 Cluster01", length(peakcar)),
                                rep("CD8 Cluster01", length(peakcar))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### spearman correlation
  s_cor <- c(round(cor(result_table[names(peakcar),"Estimated Injected GMP CD4 Subsisters # (per kg)"], peakcar, method = "spearman"), 2),
             round(cor(result_table[names(peakcar), "Estimated Injected GMP CD8 Subsisters # (per kg)"], peakcar, method = "spearman"), 2),
             round(cor(result_table2[names(peakcar),"Estimated Injected GMP CD4 Subsisters # (per kg)"], peakcar, method = "spearman"), 2),
             round(cor(result_table2[names(peakcar), "Estimated Injected GMP CD8 Subsisters # (per kg)"], peakcar, method = "spearman"), 2))
  names(s_cor) <- c("CD4 Subsisters", "CD8 Subsisters", "CD4 Cluster01", "CD8 Cluster01")
  
  ### draw the correlation plot
  p <- ggplot(data = plot_df, aes(x=X, y=Y)) +
    geom_point(aes_string(color="Group"), size = 8) +
    labs(title = paste0("Spearman Correlation\n", paste0(names(s_cor), ": ", s_cor, collapse = ", "))) +
    xlab("Estimated Injected GMP # (per kg)") +
    ylab("PeakCAR") +
    geom_smooth(method = lm, color="black", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24))
  ggsave(file = paste0(outputDir2, "Correlation_GMP_Dose_vs_PeakCAR.png"), plot = p, width = 18, height = 10, dpi = 400)
  
  
  #
  ### 9. The number of cells and clones in each cluster
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/9/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  
  ### make an empty table
  out_table <- data.frame(matrix(0, length(unique(sub_seurat_obj2@meta.data$clusters)), 4),
                          stringsAsFactors = FALSE, check.names = FALSE)
  rownames(out_table) <- unique(sub_seurat_obj2@meta.data$clusters)
  colnames(out_table) <- c("Cell #", "Clone #", "Duplicated Cell #", "Subsister Cell #")
  
  ### ordering
  out_table <- out_table[order(as.numeric(rownames(out_table))),]
  
  ### for each cluster
  for(cluster in rownames(out_table)) {
    ### some indicies
    cluster_idx <- which(sub_seurat_obj2@meta.data$clusters == cluster)
    
    ### fill out the table
    out_table[cluster, "Cell #"] <- length(cluster_idx)
    out_table[cluster, "Clone #"] <- length(unique(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta[cluster_idx]))
    out_table[cluster, "Subsister Cell #"] <- length(intersect(which(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones),
                                                                cluster_idx))
    out_table[cluster, "Duplicated Cell #"] <- length(which(duplicated(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta[cluster_idx])))
  }
  
  ### save the result table
  write.xlsx2(data.frame(Cluster=rownames(out_table),
                         out_table,
                         stringsAsFactors = FALSE, check.names = FALSE),
              sheetName = "Numbers_In_PI_Clusters",
              file = paste0(outputDir2, "Numbers_In_PI_Clusters.xlsx"),
              row.names = FALSE)
  
  ### draw a pie chart to show the percentage of subsisters in each cluster
  plot_df <- data.frame(Cluster=rownames(out_table),
                        Numbers=out_table$`Subsister Cell #`,
                        Pcnt=round(out_table$`Subsister Cell #`*100/sum(out_table$`Subsister Cell #`),1),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df <- plot_df[order(as.numeric(plot_df$Cluster)),]
  plot_df$Cluster <- factor(plot_df$Cluster, levels = unique(plot_df$Cluster))
  p <- ggplot(data = plot_df,
         aes(x = "", y = Numbers, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta="y") +
    labs(x = NULL, y = NULL, title = "Subsister # in Post-Infusion Clusters") +
    scale_fill_discrete(name = "Cluster", labels = paste0(plot_df$Cluster, ": ",
                                                          plot_df$Numbers, " (",
                                                          plot_df$Pcnt, "%)")) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", size = 48),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  ggsave(paste0(outputDir2, "Subsister_Numbers_In_PI_Clusters.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  #
  ### 10. Gini index of the CAR+ lineages over time
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/10/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get CARpos-only seurat object
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$CAR)
  sub_seurat_obj3 <- subset(Seurat_Obj, idents = c("CARpos"))
  
  ### after gmp time points only
  gmp_after_time_points <- c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", 
                             "Wk6", "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj3 <- SetIdent(object = sub_seurat_obj3,
                              cells = rownames(sub_seurat_obj3@meta.data),
                              value = sub_seurat_obj3@meta.data$time2)
  sub_seurat_obj3 <- subset(sub_seurat_obj3, idents = intersect(gmp_after_time_points,
                                                                unique(sub_seurat_obj3@meta.data$time2)))
  
  ### get seurat object for some specific patients
  sub_seurat_obj3 <- SetIdent(object = sub_seurat_obj3,
                             cells = rownames(sub_seurat_obj3@meta.data),
                             value = sub_seurat_obj3@meta.data$px)
  sub_seurat_obj3 <- subset(sub_seurat_obj3, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                       "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### cells only with the TCR info
  sub_seurat_obj3 <- SetIdent(object = sub_seurat_obj3,
                              cells = rownames(sub_seurat_obj3@meta.data),
                              value = sub_seurat_obj3@meta.data$ALL_GMP_CARpos_Persister)
  sub_seurat_obj3 <- subset(sub_seurat_obj3, idents = c("YES", "NO"))
  
  ### get time points
  tps <- unique(sub_seurat_obj3@meta.data$time2)
  
  ### for each time point
  for(tp in tps) {
    
    ### time point index
    tp_idx <- which(sub_seurat_obj3@meta.data$time2 == tp)
    
    ### unique clones for the time point
    clones_tp <- sub_seurat_obj3@meta.data$clonotype_id_by_patient_one_alpha_beta[tp_idx]
    clones_tp <- clones_tp[which(!is.na(clones_tp))]
    unique_clones_tp <- unique(clones_tp)
    
    ### make an empty table
    out_table <- data.frame(Clone_Name=unique_clones_tp,
                            Clone_Size=1,
                            TP=tp,
                            stringsAsFactors = FALSE, check.names = FALSE)
    rownames(out_table) <- unique_clones_tp
    
    ### compute the clone size
    dups <- clones_tp[which(duplicated(clones_tp))]
    unique_dups <- unique(dups)
    for(dup in unique_dups) {
      out_table[dup,"Clone_Size"] <- length(intersect(tp_idx,
                                                      which(sub_seurat_obj3@meta.data$clonotype_id_by_patient_one_alpha_beta == dup)))
    }
    
    ### combine the table
    if(tp == tps[1]) {
      plot_df <- out_table
    } else {
      plot_df <- rbind(plot_df, out_table)
    }
  }
  
  ### factorize the time column
  plot_df$TP <- factor(plot_df$TP,
                       levels = tps)
  
  ### gini index
  gini_idx <- sapply(tps, function(x) {
    return(Gini(plot_df$Clone_Size[plot_df$TP == x]))
  })
  
  ### draw a density plot
  p <- ggplot(plot_df, aes_string(x="Clone_Size", col="TP")) +
    xlim(c(0,3)) +
    geom_density(size = 2) +
    labs(title = "Density - All CAR+", color="") +
    ylab("") +
    theme_classic(base_size = 30) +
    scale_color_discrete(name = "Gini Idx", labels = paste0(names(gini_idx), ": ", round(gini_idx, 2))) +
    theme(legend.title = element_text(size = 50),
          legend.text = element_text(size = 40))
  ggsave(paste0(outputDir2, "CAR+_Evenness_Gini.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### per patient
  plot_df$Px <- sapply(plot_df$Clone_Name, function(x) {
    return(strsplit(x, split = "_", fixed = TRUE)[[1]][1])
  })
  p <- vector("list", length = length(unique(plot_df$Px)))
  names(p) <- unique(plot_df$Px)
  for(px in unique(plot_df$Px)) {
    temp_plot_df <- plot_df[which(plot_df$Px == px),]
    gini_idx1 <- sapply(tps, function(x) {
      return(Gini(temp_plot_df$Clone_Size[temp_plot_df$TP == x]))
    })
    gini_idx1 <- gini_idx1[which(!is.nan(gini_idx1))]
    temp_plot_df <- temp_plot_df[which(temp_plot_df$TP %in% names(gini_idx1)),]
    p[[px]] <- ggplot(temp_plot_df, aes_string(x="Clone_Size", col="TP")) +
      xlim(c(0,3)) +
      geom_density(size = 2) +
      labs(title = px, color="") +
      ylab("") +
      theme_classic(base_size = 30) +
      scale_color_discrete(name = "Gini Idx", labels = paste0(names(gini_idx1), ": ", round(gini_idx1, 2))) +
      theme(legend.title = element_text(size = 50),
            legend.text = element_text(size = 40))
    
  }
  
  ### combine the plots 
  g <- arrangeGrob(grobs = p,
                   ncol = 3,
                   top = "")
  ggsave(file = paste0(outputDir2, "CAR+_Evenness_Gini_per_px.png"), g, width = 40, height = 20, dpi = 350)
  
  
  #
  ### only withe CD8 CAR+ lineages
  #
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  
  ### get CD8 CAR+ lineages only
  plot_df2 <- plot_df[which(plot_df$Clone_Name %in% persister_clones),]
  
  ### gini index
  gini_idx2 <- sapply(tps, function(x) {
    return(Gini(plot_df2$Clone_Size[plot_df2$TP == x]))
  })
  
  ### draw a density plot
  p <- ggplot(plot_df2, aes_string(x="Clone_Size", col="TP")) +
    xlim(c(0,5)) +
    geom_density(size = 2) +
    labs(title = "Density - CAR+ Lineages", color="") +
    ylab("") +
    theme_classic(base_size = 30) +
    scale_color_discrete(name = "Gini Idx", labels = paste0(names(gini_idx2), ": ", round(gini_idx2, 2))) +
    theme(legend.title = element_text(size = 50),
          legend.text = element_text(size = 40))
  ggsave(paste0(outputDir2, "CAR+_Lineages_Evenness_Gini.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### per patient
  plot_df2$Px <- sapply(plot_df2$Clone_Name, function(x) {
    return(strsplit(x, split = "_", fixed = TRUE)[[1]][1])
  })
  p <- vector("list", length = length(unique(plot_df2$Px)))
  names(p) <- unique(plot_df2$Px)
  for(px in unique(plot_df2$Px)) {
    temp_plot_df <- plot_df2[which(plot_df2$Px == px),]
    gini_idx3 <- sapply(tps, function(x) {
      return(Gini(temp_plot_df$Clone_Size[temp_plot_df$TP == x]))
    })
    gini_idx3 <- gini_idx3[which(!is.nan(gini_idx3))]
    temp_plot_df <- temp_plot_df[which(temp_plot_df$TP %in% names(gini_idx3)),]
    p[[px]] <- ggplot(temp_plot_df, aes_string(x="Clone_Size", col="TP")) +
      xlim(c(0,5)) +
      geom_density(size = 2) +
      labs(title = px, color="") +
      ylab("") +
      theme_classic(base_size = 30) +
      scale_color_discrete(name = "Gini Idx", labels = paste0(names(gini_idx3), ": ", round(gini_idx3, 2))) +
      theme(legend.title = element_text(size = 50),
            legend.text = element_text(size = 40))
    
  }
  
  ### combine the plots 
  g <- arrangeGrob(grobs = p,
                   ncol = 3,
                   top = "")
  ggsave(file = paste0(outputDir2, "CAR+_Lineages_Evenness_Gini_per_px.png"), g, width = 40, height = 20, dpi = 350)
  
  
  #
  ### 11. DE result - Violin & Dot plot
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/11/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### only get the GMP persisters and non-persisters
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[union(which(Seurat_Obj$GMP_CARpos_CD8_Persister == "YES"),
                                                                                       which(Seurat_Obj$GMP_CARpos_CD8_Persister == "NO"))])
  
  ### only the patients we are interested
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  target_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                            "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                            "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### get the same number of samples (persister: non-persister - same ratio)
  set.seed(1234)
  target_Seurat_Obj@meta.data$New_Persistency <- NA
  multiplier_k <- 1
  for(px in unique(target_Seurat_Obj@meta.data$px)) {
    ### get specific indicies
    px_gmp_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                             which(target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES"))
    px_gmp_not_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                                 which(target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "NO"))
    
    ### sampling
    if((length(px_gmp_last) > 0) && (length(px_gmp_not_last) > length(px_gmp_last))) {
      px_gmp_not_last <- sample(px_gmp_not_last, size = length(px_gmp_last)*multiplier_k)
    }
    
    ### annotate new persistency info
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_last] <- "YES"
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_not_last] <- "NO"
  }
  
  ### set idents with the new info
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$New_Persistency)
  
  ### only using the specific cells
  target_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("YES", "NO"))
  
  ### change the labels
  target_Seurat_Obj@meta.data$New_Persistency[which(target_Seurat_Obj@meta.data$New_Persistency == "YES")] <- "Subsisters"
  target_Seurat_Obj@meta.data$New_Persistency[which(target_Seurat_Obj@meta.data$New_Persistency == "NO")] <- "Non-Subsisters"
  
  ### violin plot
  Idents(target_Seurat_Obj) <- target_Seurat_Obj@meta.data$New_Persistency
  p <- VlnPlot(target_Seurat_Obj, features = c("MAL", "IFITM3", "CD27", "SELL", "TIGIT",
                                               "HLA-DQA1", "CCL4", "GZMH", "IFNG"),
               pt.size = 0)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
      stat_compare_means(size = 8) +
      stat_summary(fun=mean, geom="point", size=5, color="red") +
      theme_classic(base_size = 40) +
      theme(legend.position = "none",
            axis.title.x = element_blank())
  }
  
  ### save the violin plot
  ggsave(file = paste0(outputDir2, "GMP_CD8_CARpos_Subsisters_vs_Non-Subsisters_Downsampled.png"), plot = p, width = 25, height = 25, dpi = 350)
  
  ### violin plot
  Idents(target_Seurat_Obj) <- target_Seurat_Obj@meta.data$New_Persistency
  p <- VlnPlot(target_Seurat_Obj, features = c("CAPG", "MAL", "CD52", "IFITM2", "IFITM3",
                                               "SPOCK2", "SELL", "CD27", "CD7"),
               pt.size = 0)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
      stat_compare_means(size = 8) +
      stat_summary(fun=mean, geom="point", size=5, color="red") +
      theme_classic(base_size = 40) +
      theme(legend.position = "none",
            axis.title.x = element_blank())
  }
  
  ### save the violin plot
  ggsave(file = paste0(outputDir2, "GMP_CD8_CARpos_Subsisters_vs_Non-Subsisters_Downsampled(2).png"), plot = p, width = 25, height = 25, dpi = 350)
  
  ### violin plot
  Idents(target_Seurat_Obj) <- target_Seurat_Obj@meta.data$New_Persistency
  p <- VlnPlot(target_Seurat_Obj, features = c("HLA-DRB5", "HLA-DPA1", "LAG3", "HLA-DQB1", "NPM1",
                                               "CCL5", "FABP5", "CD74", "CD70"),
               pt.size = 0)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
      stat_compare_means(size = 8) +
      stat_summary(fun=mean, geom="point", size=5, color="red") +
      theme_classic(base_size = 40) +
      theme(legend.position = "none",
            axis.title.x = element_blank())
  }
  
  ### save the violin plot
  ggsave(file = paste0(outputDir2, "GMP_CD8_CARpos_Subsisters_vs_Non-Subsisters_Downsampled(3).png"), plot = p, width = 25, height = 25, dpi = 350)
  
  # ### change the labels
  # target_Seurat_Obj@meta.data$Cluster01[which(target_Seurat_Obj@meta.data$Cluster01 == "YES")] <- "Cluster01"
  # target_Seurat_Obj@meta.data$Cluster01[which(target_Seurat_Obj@meta.data$Cluster01 == "NO")] <- "Other_Clusters"
  # 
  # ### violin plot
  # Idents(target_Seurat_Obj) <- target_Seurat_Obj@meta.data$Cluster01
  # p <- VlnPlot(target_Seurat_Obj, features = c("MAL", "IFITM3", "CD27", "SELL", "TIGIT",
  #                                              "HLA-DQA1", "CCL4", "GZMH", "IFNG"),
  #              pt.size = 0)
  # for(i in 1:9) {
  #   p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
  #     stat_compare_means(size = 8) +
  #     stat_summary(fun=mean, geom="point", size=5, color="red") +
  #     theme_classic(base_size = 40) +
  #     theme(legend.position = "none",
  #           axis.title.x = element_blank())
  # }
  # 
  # ### save the violin plot
  # ggsave(file = paste0(outputDir2, "GMP_CD8_CARpos_Cluster01_vs_Others_Downsampled.png"), plot = p, width = 25, height = 25, dpi = 350)
  # 
  
  ### get CARpos-only seurat object
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$CAR)
  sub_seurat_obj <- subset(Seurat_Obj, idents = c("CARpos"))
  
  ### after gmp time points only
  after_gmp_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6",
                             "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$time2)
  sub_seurat_obj <- subset(sub_seurat_obj, idents = intersect(after_gmp_time_points,
                                                              unique(sub_seurat_obj@meta.data$time2)))
  
  ### normalization
  sub_seurat_obj <- NormalizeData(sub_seurat_obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  sub_seurat_obj <- FindVariableFeatures(sub_seurat_obj,
                                         selection.method = "vst", nfeatures = 2000)
  ### scaling
  sub_seurat_obj <- ScaleData(sub_seurat_obj,
                              vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### get seurat object for some specific patients
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$px)
  sub_seurat_obj2 <- subset(sub_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                       "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  
  ### define persisters
  sub_seurat_obj2@meta.data$CD8_Persisters <- "Non-subsisters"
  sub_seurat_obj2@meta.data$CD8_Persisters[which(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones)] <- "Subsisters"
  sub_seurat_obj2@meta.data$CD8_Persisters <- factor(sub_seurat_obj2@meta.data$CD8_Persisters, levels = c("Subsisters", "Non-subsisters"))
  
  ### set ident with the persistency info
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$CD8_Persisters)
  
  ### violin plot
  p <- VlnPlot(sub_seurat_obj2, features = c("MAL", "IFITM3", "CD27", "SELL", "TIGIT",
                                               "HLA-DQA1", "CCL4", "GZMH", "IFNG"),
               pt.size = 0)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
      stat_compare_means(size = 8) +
      stat_summary(fun=mean, geom="point", size=5, color="red") +
      theme_classic(base_size = 40) +
      theme(legend.position = "none",
            axis.title.x = element_blank())
  }
  
  ### save the violin plot
  ggsave(file = paste0(outputDir2, "Post-Infusion_CD8_CARpos_Subsisters_vs_Non-Subsisters.png"), plot = p, width = 25, height = 25, dpi = 350)
  
  # ### change the labels
  # sub_seurat_obj2@meta.data$Cluster01[which(sub_seurat_obj2@meta.data$Cluster01 == "YES")] <- "Cluster01"
  # sub_seurat_obj2@meta.data$Cluster01[which(sub_seurat_obj2@meta.data$Cluster01 == "NO")] <- "Other_Clusters"
  # 
  # ### set ident with the persistency info
  # sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
  #                             cells = rownames(sub_seurat_obj2@meta.data),
  #                             value = sub_seurat_obj2@meta.data$Cluster01)
  # 
  # ### violin plot
  # p <- VlnPlot(sub_seurat_obj2, features = c("MAL", "IFITM3", "CD27", "SELL", "TIGIT",
  #                                            "HLA-DQA1", "CCL4", "GZMH", "IFNG"),
  #              pt.size = 0)
  # for(i in 1:9) {
  #   p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
  #     stat_compare_means(size = 8) +
  #     stat_summary(fun=mean, geom="point", size=5, color="red") +
  #     theme_classic(base_size = 40) +
  #     theme(legend.position = "none",
  #           axis.title.x = element_blank())
  # }
  # 
  # ### save the violin plot
  # ggsave(file = paste0(outputDir2, "Post-Infusion_CD8_CARpos_Cluster01_vs_Others.png"), plot = p, width = 25, height = 25, dpi = 350)
  
  
  #
  ### 12. BM cells in UMAP & their info + compare repertoires and GEx profiles- are BM and PB CAR+ cells similar?
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/12/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### pre indicies
  bm_gmp_idx <- which(BM_Seurat_obj@meta.data$time2 == "GMP")
  bm_ai_idx <- which(BM_Seurat_obj@meta.data$time2 %in% c("Wk4", "Wk8", "3mo"))
  bm_carpos_idx <- which(BM_Seurat_obj@meta.data$CAR == "CARpos")
  bm_ai_carpos_idx <- intersect(bm_ai_idx, bm_carpos_idx)
  bm_gmp_carpos_idx <- intersect(bm_gmp_idx, bm_carpos_idx)
  bm_gmp_carpos_cd8_idx <- intersect(bm_gmp_carpos_idx,
                                     which(BM_Seurat_obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  bm_gmp_carpos_clones <- unique(BM_Seurat_obj$clonotype_id_by_patient_one_alpha_beta[bm_gmp_carpos_idx])
  bm_gmp_carpos_clones <- bm_gmp_carpos_clones[which(!is.na(bm_gmp_carpos_clones))]
  bm_ai_carpos_clones <- unique(BM_Seurat_obj$clonotype_id_by_patient_one_alpha_beta[bm_ai_carpos_idx])
  bm_ai_carpos_clones <- bm_ai_carpos_clones[which(!is.na(bm_ai_carpos_clones))]
  bm_gmp_ai_carpos_clones <- intersect(bm_gmp_carpos_clones, bm_ai_carpos_clones)
  
  ### define subsisters in BM
  
  ### remove ALL_CARpos_Persister
  BM_Seurat_obj@meta.data$ALL_CARpos_Persister <- NULL
  
  ### BM GMP CARpos lineages at least one in GMP and others in post-infusion
  ### including all time points
  ### "YES": CAR+ Persister cells in GMP lineages (all time)
  ### "NO": CAR+ cells that are not in the GMP lineages
  ### NA: Others
  BM_Seurat_obj@meta.data$ALL_GMP_CARpos_Persister <- NA
  BM_Seurat_obj@meta.data$ALL_GMP_CARpos_Persister[which(BM_Seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% bm_gmp_ai_carpos_clones)] <- "YES"
  BM_Seurat_obj@meta.data$ALL_GMP_CARpos_Persister[setdiff(bm_carpos_idx,
                                                           which(BM_Seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% bm_gmp_ai_carpos_clones))] <- "NO"
  
  ### BM GMP CARpos lineages at least one in GMP and others in post-infusion
  ### including GMP CELLS ONLY
  ### "YES": GMP CAR+ Persister cells
  ### "NO": GMP CAR+ cells that are non-persisters
  ### NA: Others
  BM_Seurat_obj@meta.data$GMP_CARpos_Persister <- NA
  BM_Seurat_obj@meta.data$GMP_CARpos_Persister[intersect(which(BM_Seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% bm_gmp_ai_carpos_clones),
                                                         bm_gmp_carpos_idx)] <- "YES"
  BM_Seurat_obj@meta.data$GMP_CARpos_Persister[setdiff(bm_gmp_carpos_idx,
                                                       which(BM_Seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% bm_gmp_ai_carpos_clones))] <- "NO"
  
  ### BM CD8 cells in GMP CARpos lineages at least one in GMP and others in post-infusion
  ### including GMP CELLS ONLY
  ### "YES": GMP CAR+ CD8 Persister cells
  ### "NO": GMP CAR+ CD8 cells that are non-persisters
  ### NA: Others
  BM_Seurat_obj@meta.data$GMP_CARpos_CD8_Persister <- NA
  BM_Seurat_obj@meta.data$GMP_CARpos_CD8_Persister[intersect(which(BM_Seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% bm_gmp_ai_carpos_clones),
                                                             bm_gmp_carpos_cd8_idx)] <- "YES"
  BM_Seurat_obj@meta.data$GMP_CARpos_CD8_Persister[setdiff(bm_gmp_carpos_cd8_idx,
                                                           which(BM_Seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% bm_gmp_ai_carpos_clones))] <- "NO"
  
  
  ### get CARpos-only seurat object
  BM_Seurat_obj <- SetIdent(object = BM_Seurat_obj,
                         cells = rownames(BM_Seurat_obj@meta.data),
                         value = BM_Seurat_obj@meta.data$CAR)
  BM_sub_seurat_obj <- subset(BM_Seurat_obj, idents = c("CARpos"))
  
  ### after gmp time points only
  total_time_points <- c("PreTrans", "Wk-1", "Wk0", "GMP", "Wk1", "Wk2", 
                         "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo")
  after_gmp_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6",
                             "Wk8", "3mo", "6mo", "9mo")
  BM_sub_seurat_obj <- SetIdent(object = BM_sub_seurat_obj,
                                cells = rownames(BM_sub_seurat_obj@meta.data),
                                value = BM_sub_seurat_obj@meta.data$time2)
  BM_sub_seurat_obj <- subset(BM_sub_seurat_obj, idents = intersect(after_gmp_time_points,
                                                                    unique(BM_sub_seurat_obj@meta.data$time2)))
  
  #
  ### run pca on the BM_sub_seurat_obj
  #
  ### normalization
  BM_sub_seurat_obj <- NormalizeData(BM_sub_seurat_obj,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  BM_sub_seurat_obj <- FindVariableFeatures(BM_sub_seurat_obj,
                                            selection.method = "vst", nfeatures = 2000)
  ### scaling
  BM_sub_seurat_obj <- ScaleData(BM_sub_seurat_obj,
                                 vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  ### PCA
  BM_sub_seurat_obj <- RunPCA(BM_sub_seurat_obj,
                              features = VariableFeatures(object = BM_sub_seurat_obj),
                              npcs = 15)
  ### UMAP
  BM_sub_seurat_obj <- RunUMAP(BM_sub_seurat_obj, dims = 1:15)
  
  ### get seurat object for some specific patients
  BM_sub_seurat_obj <- SetIdent(object = BM_sub_seurat_obj,
                                cells = rownames(BM_sub_seurat_obj@meta.data),
                                value = BM_sub_seurat_obj@meta.data$px)
  BM_sub_seurat_obj2 <- subset(BM_sub_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                             "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                             "SJCAR19-10", "SJCAR19-11"))
  
  ### factorize the time2 column
  BM_sub_seurat_obj2@meta.data$time2 <- factor(BM_sub_seurat_obj2@meta.data$time2,
                                               levels = intersect(total_time_points, BM_sub_seurat_obj2@meta.data$time2))
  
  ### perform clustering
  BM_sub_seurat_obj2 <- FindNeighbors(BM_sub_seurat_obj2, dims = 1:10)
  BM_sub_seurat_obj2 <- FindClusters(BM_sub_seurat_obj2, resolution = 0.5)
  
  ### save the clustering result to meta.data
  BM_sub_seurat_obj2@meta.data$clusters <- Idents(BM_sub_seurat_obj2)
  
  ### UMAP with subsisters in post-infusion
  p <- DimPlot(object = BM_sub_seurat_obj2, reduction = "umap",
               group.by = "ALL_GMP_CARpos_Persister",
               pt.size = 5,
               cols = c("NO" = "lightgray", "YES" = "red"),
               order = c("YES", "NO")) +
    ggtitle("") +
    labs(color="Is_Subsistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_BM_PI_CARpos_Subsister_ALL.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### UMAP with subsister info split by each patient
  p <- DimPlot(object = BM_sub_seurat_obj2, reduction = "umap",
               group.by = "ALL_GMP_CARpos_Persister", split.by = "px",
               pt.size = 5, ncol = 3,
               cols = c("NO" = "lightgray", "YES" = "red"),
               order = c("YES", "NO")) +
    ggtitle("") +
    labs(color="Is_Subsistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_BM_PI_CARpos_Subsister.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### DE analysis - BM Post Infusion subsister vs non-subsister
  BM_sub_seurat_obj2 <- SetIdent(object = BM_sub_seurat_obj2,
                                 cells = rownames(BM_sub_seurat_obj2@meta.data),
                                 value = BM_sub_seurat_obj2@meta.data$ALL_GMP_CARpos_Persister)
  de_result <- FindMarkers(BM_sub_seurat_obj2,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/BM_PI_CARpos_Subsister_vs_Non-subsister.xlsx"),
              sheetName = "BM_PI_CARpos_Subsister_vs_Non-subsister", row.names = FALSE)
  
  ### get CARpos-only seurat object
  BM_Seurat_obj <- SetIdent(object = BM_Seurat_obj,
                            cells = rownames(BM_Seurat_obj@meta.data),
                            value = BM_Seurat_obj@meta.data$CAR)
  BM_sub_seurat_obj3 <- subset(BM_Seurat_obj, idents = c("CARpos"))
  
  ### get GMP seurat
  BM_sub_seurat_obj3 <- SetIdent(object = BM_sub_seurat_obj3,
                                 cells = rownames(BM_sub_seurat_obj3@meta.data),
                                 value = BM_sub_seurat_obj3@meta.data$time2)
  BM_sub_seurat_obj3 <- subset(BM_sub_seurat_obj3, idents = c("GMP"))
  
  ### normalization
  BM_sub_seurat_obj3 <- NormalizeData(BM_sub_seurat_obj3,
                                      normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  BM_sub_seurat_obj3 <- FindVariableFeatures(BM_sub_seurat_obj3,
                                             selection.method = "vst", nfeatures = 2000)
  ### scaling
  BM_sub_seurat_obj3 <- ScaleData(BM_sub_seurat_obj3,
                                  vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  ### PCA
  BM_sub_seurat_obj3 <- RunPCA(BM_sub_seurat_obj3,
                               features = VariableFeatures(object = BM_sub_seurat_obj3),
                               npcs = 15)
  ### UMAP
  BM_sub_seurat_obj3 <- RunUMAP(BM_sub_seurat_obj3, dims = 1:15)
  
  ### get seurat object for some specific patients
  BM_sub_seurat_obj3 <- SetIdent(object = BM_sub_seurat_obj3,
                                 cells = rownames(BM_sub_seurat_obj3@meta.data),
                                 value = BM_sub_seurat_obj3@meta.data$px)
  BM_sub_seurat_obj3 <- subset(BM_sub_seurat_obj3, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                             "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                             "SJCAR19-10", "SJCAR19-11"))
  
  ### DE analysis - BM GMP subsister vs non-subsister
  BM_sub_seurat_obj3 <- SetIdent(object = BM_sub_seurat_obj3,
                                 cells = rownames(BM_sub_seurat_obj3@meta.data),
                                 value = BM_sub_seurat_obj3@meta.data$GMP_CARpos_CD8_Persister)
  de_result <- FindMarkers(BM_sub_seurat_obj3,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/BM_GMP_CARpos_Subsister_vs_Non-subsister.xlsx"),
              sheetName = "BM_GMP_CARpos_Subsister_vs_Non-subsister", row.names = FALSE)
  
  ### PB GMP CARpos CD8 Subsisters vs BM GMP CARpos CD8 Subsisters
  pb_gmp_carpos_cd8_persister_idx <- which(rownames(Seurat_Obj@meta.data) %in% rownames(target_Seurat_Obj@meta.data)[which(target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  bm_gmp_carpos_cd8_persister_idx <- which(rownames(Seurat_Obj@meta.data) %in% rownames(BM_sub_seurat_obj3@meta.data)[which(BM_sub_seurat_obj3@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  ### 23 clones intersected
  intersect_gmp_carpos_cd8_persister_idx <- intersect(pb_gmp_carpos_cd8_persister_idx,
                                                      bm_gmp_carpos_cd8_persister_idx)
  intersect_gmp_carpos_cd8_persister_clones <- Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[intersect_gmp_carpos_cd8_persister_idx]
  Seurat_Obj@meta.data$PB_BM_GMP_CD8_Persister <- NA
  Seurat_Obj@meta.data$PB_BM_GMP_CD8_Persister[setdiff(pb_gmp_carpos_cd8_persister_idx,
                                                       bm_gmp_carpos_cd8_persister_idx)] <- "PB_GMP_CARpos_CD8_Subsisters"
  Seurat_Obj@meta.data$PB_BM_GMP_CD8_Persister[setdiff(bm_gmp_carpos_cd8_persister_idx,
                                                       pb_gmp_carpos_cd8_persister_idx)] <- "BM_GMP_CARpos_CD8_Subsisters"
  
  ### DE analysis - PB GMP CARpos CD8 Subsisters vs BM GMP CARpos CD8 Subsisters
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$PB_BM_GMP_CD8_Persister)
  de_result <- FindMarkers(Seurat_Obj,
                           ident.1 = "PB_GMP_CARpos_CD8_Subsisters",
                           ident.2 = "BM_GMP_CARpos_CD8_Subsisters",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CARpos_CD8_Subsister_PB_vs_BM.xlsx"),
              sheetName = "GMP_CARpos_CD8_Subsister_PB_vs_BM", row.names = FALSE)
  
  ### get GMP cells only
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$time2)
  pb_bm_gmp_seurat_obj <- subset(Seurat_Obj, idents = c("GMP"))
  
  ### get seurat object for some specific patients
  pb_bm_gmp_seurat_obj <- SetIdent(object = pb_bm_gmp_seurat_obj,
                                   cells = rownames(pb_bm_gmp_seurat_obj@meta.data),
                                   value = pb_bm_gmp_seurat_obj@meta.data$px)
  pb_bm_gmp_seurat_obj <- subset(pb_bm_gmp_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                                  "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                                  "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### get seurat object for specific cells of interests
  pb_bm_gmp_seurat_obj <- subset(pb_bm_gmp_seurat_obj, cells = rownames(pb_bm_gmp_seurat_obj@meta.data)[intersect(which(pb_bm_gmp_seurat_obj@meta.data$CAR == "CARpos"),
                                                                                                                  which(pb_bm_gmp_seurat_obj@meta.data$CD4_CD8_by_Consensus == "CD8"))])
  
  ### normalization
  pb_bm_gmp_seurat_obj <- NormalizeData(pb_bm_gmp_seurat_obj,
                                        normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  pb_bm_gmp_seurat_obj <- FindVariableFeatures(pb_bm_gmp_seurat_obj,
                                               selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  pb_bm_gmp_seurat_obj <- ScaleData(pb_bm_gmp_seurat_obj,
                                    vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  pb_bm_gmp_seurat_obj <- RunPCA(pb_bm_gmp_seurat_obj,
                                 features = VariableFeatures(object = pb_bm_gmp_seurat_obj),
                                 npcs = 15)
  pb_bm_gmp_seurat_obj <- RunUMAP(pb_bm_gmp_seurat_obj, dims = 1:15)
  
  ### NA to "Others"
  pb_bm_gmp_seurat_obj@meta.data$PB_BM_GMP_CD8_Persister[which(is.na(pb_bm_gmp_seurat_obj@meta.data$PB_BM_GMP_CD8_Persister))] <- "Others"
  
  ### change UMAP coordinates consistency with the previous result
  ### doesn't affect the result itself
  pb_bm_gmp_seurat_obj@reductions$umap@cell.embeddings <- -pb_bm_gmp_seurat_obj@reductions$umap@cell.embeddings
  pb_bm_gmp_seurat_obj@reductions$umap@feature.loadings <- -pb_bm_gmp_seurat_obj@reductions$umap@feature.loadings
  
  ### GMP UMAP difference between PB & BM subsisters
  p <- DimPlot(object = pb_bm_gmp_seurat_obj, reduction = "umap",
               group.by = "PB_BM_GMP_CD8_Persister",
               pt.size = 5,
               cols = c("PB_GMP_CARpos_CD8_Subsisters" = "red", "BM_GMP_CARpos_CD8_Subsisters" = "blue", "Others" = "gray"),
               order = c("BM_GMP_CARpos_CD8_Subsisters", "PB_GMP_CARpos_CD8_Subsisters", "Others")) +
    ggtitle("") +
    labs(color="Is_Subsistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_PB_BM_GMP_CARpos_Subsister_ALL.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### UMAP with subsister info split by each patient
  p <- DimPlot(object = pb_bm_gmp_seurat_obj, reduction = "umap",
               group.by = "PB_BM_GMP_CD8_Persister", split.by = "px",
               pt.size = 5, ncol = 3,
               cols = c("PB_GMP_CARpos_CD8_Subsisters" = "red", "BM_GMP_CARpos_CD8_Subsisters" = "blue", "Others" = "gray"),
               order = c("BM_GMP_CARpos_CD8_Subsisters", "PB_GMP_CARpos_CD8_Subsisters", "Others")) +
    ggtitle("") +
    labs(color="Is_Subsistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_PB_BM_GMP_CARpos_Subsister.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  
  #
  ### 13. Clustering on GMP UMAP and find subsister clusters + calculate dose (related to 5(k))
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/13/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### only get the GMP persisters and non-persisters
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[union(which(Seurat_Obj$GMP_CARpos_CD8_Persister == "YES"),
                                                                                       which(Seurat_Obj$GMP_CARpos_CD8_Persister == "NO"))])
  
  ### only the patients we are interested
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  target_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                            "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                            "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### normalization
  target_Seurat_Obj <- NormalizeData(target_Seurat_Obj,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj <- FindVariableFeatures(target_Seurat_Obj,
                                            selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj <- ScaleData(target_Seurat_Obj,
                                 vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj <- RunPCA(target_Seurat_Obj,
                              features = VariableFeatures(object = target_Seurat_Obj),
                              npcs = 15)
  target_Seurat_Obj <- RunUMAP(target_Seurat_Obj, dims = 1:15)
  
  ### perform clustering
  target_Seurat_Obj <- FindNeighbors(target_Seurat_Obj, dims = 1:10)
  target_Seurat_Obj <- FindClusters(target_Seurat_Obj, resolution = 1.5)
  
  ### save the clustering result to meta.data
  target_Seurat_Obj@meta.data$gmp_seurat_clusters <- Idents(target_Seurat_Obj)
  
  ### specific cluster
  DimPlot(object = target_Seurat_Obj, reduction = "umap",
          group.by = "gmp_seurat_clusters",
          pt.size = 1)
  target_Seurat_Obj@meta.data$gmp_seurat_clusters2 <- as.character(target_Seurat_Obj@meta.data$gmp_seurat_clusters)
  target_Seurat_Obj@meta.data$gmp_seurat_clusters2[which(target_Seurat_Obj@meta.data$gmp_seurat_clusters2 != "22")] <- "Others"
  DimPlot(object = target_Seurat_Obj, reduction = "umap",
          group.by = "gmp_seurat_clusters2",
          pt.size = 1,
          cols = c("Others" = "lightgray", "22" = "red"),
          order = c("22", "Others"))
  
  ### GMP UMAP difference between PB & BM subsisters
  p <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
               group.by = "gmp_seurat_clusters",
               pt.size = 3) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_PB_GMP_CARpos_CD8_Clusters_ALL.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### UMAP with subsister info split by each patient
  p <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
               group.by = "gmp_seurat_clusters", split.by = "px",
               pt.size = 3, ncol = 3) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_PB_GMP_CARpos_CD8_Clusters.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### UMAP cluster22 and the others
  p <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
          group.by = "gmp_seurat_clusters2",
          pt.size = 3,
          cols = c("Others" = "lightgray", "22" = "red"),
          order = c("22", "Others")) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "UMAP_PB_GMP_CARpos_CD8_Cluster22_and_Others.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### DE analysis
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$gmp_seurat_clusters2)
  de_result <- FindMarkers(target_Seurat_Obj,
                           ident.1 = "22",
                           ident.2 = "Others",
                           min.pct = 0.5,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/PB_GMP_CARpos_CD8_Cluster22_vs_Others.xlsx"),
              sheetName = "PB_GMP_CARpos_CD8_Cluster22_DE_Result", row.names = FALSE)
  
  ### cluster22 clones
  cluster22_clones <- unique(target_Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(target_Seurat_Obj@meta.data$gmp_seurat_clusters2 == "22")])
  
  ### remove NA from cluster01 clones
  cluster22_clones <- cluster22_clones[which(!is.na(cluster22_clones))]
  
  ### add cluster22 column to the original seurat object
  Seurat_Obj@meta.data$Cluster22 <- "NO"
  Seurat_Obj@meta.data$Cluster22[which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% cluster22_clones)] <- "YES"
  
  ### the info we've got
  Real_GMP_CARpos_pcnt <- c(41.90, 37.00, 41.90, 59.10, 33.30, 23.20, 26.70, 14.70,
                            58.30, 33.60, 42.10, 32.00, 27.40, 32.60, NA, NA)
  Real_GMP_CAR_CD4_pcnt <- c(60.60, 93.60, 60.60, 75.50, 77.50, 61.40, 95.40, 85.00,
                             66.90, 52.10, 63.50, 58.60, 70.70, 71.20, NA, NA)
  Real_GMP_CAR_CD8_pcnt <- c(37.20, 5.50, 37.20, 23.60, 21.00, 36.80, 2.90, 2.40,
                             30.80, 46.70, 35.00, 39.90, 28.00, 27.70, NA, NA)
  names(Real_GMP_CARpos_pcnt) <- unique(Seurat_Obj@meta.data$px)
  names(Real_GMP_CAR_CD4_pcnt) <- unique(Seurat_Obj@meta.data$px)
  names(Real_GMP_CAR_CD8_pcnt) <- unique(Seurat_Obj@meta.data$px)
  
  ### calculate some
  Real_GMP_CAR_CD4_pcnt2 <- Real_GMP_CARpos_pcnt * Real_GMP_CAR_CD4_pcnt / 10000
  Real_GMP_CAR_CD8_pcnt2 <- Real_GMP_CARpos_pcnt * Real_GMP_CAR_CD8_pcnt / 10000
  
  ### the cell count of the GMP product
  gmp_cell_cnt <- c(1887000000, 3920000000, 2620000000, 4000000000, 3870000000,
                    4030000000, 3620000000, 2820000000, 2700000000, 4060000000,
                    4180000000, 2930000000, 5150000000, 3200000000, 3010000000, 4590000000)
  names(gmp_cell_cnt) <- unique(Seurat_Obj@meta.data$px)
  
  ### dose info
  dose_levels <- c(1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000,
                   3000000, 3000000, 3000000, 3000000, 3000000, 3000000, 3000000, 3000000, 3000000)
  
  ### result_table
  result_table2 <- matrix(0, length(unique(Seurat_Obj@meta.data$px)), 16)
  rownames(result_table2) <- unique(Seurat_Obj@meta.data$px)
  colnames(result_table2) <- c("GMP CAR+ Cell %", "CD4 % of CAR+", "CD8 % of CAR+", "GMP Total Cell #",
                               "10x GMP CAR+ CD4 Cell #", "10x GMP CAR+ CD8 Cell #",
                               "10x GMP CAR+ CD4 Subsisters #", "10x GMP CAR+ CD8 Subsisters #",
                               "GMP CD4 Subsister %", "GMP CD8 Subsister %",
                               "Estimated GMP CD4 Subsisters #", "Estimated GMP CD8 Subsisters #",
                               "Dose (per Kg)",
                               "Adjusted Dose (per Kg)",
                               "Estimated Injected GMP CD4 Subsisters # (per kg)",
                               "Estimated Injected GMP CD8 Subsisters # (per kg)")
  
  ### data frame
  result_table2 <- data.frame(result_table2, stringsAsFactors = FALSE, check.names = FALSE)
  result_table2$`GMP CAR+ Cell %` <- Real_GMP_CARpos_pcnt
  result_table2$`CD4 % of CAR+` <- Real_GMP_CAR_CD4_pcnt
  result_table2$`CD8 % of CAR+` <- Real_GMP_CAR_CD8_pcnt
  result_table2$`GMP Total Cell #` <- gmp_cell_cnt
  result_table2$`Dose (per Kg)` <- dose_levels
  result_table2$`Adjusted Dose (per Kg)` <- round(dose_levels * 100 / Real_GMP_CARpos_pcnt, 0)
  result_table2$`Estimated Injected GMP CD4 Subsisters # (per kg)` <- NA
  result_table2$`Estimated Injected GMP CD8 Subsisters # (per kg)` <- NA
  
  ### for each patient
  for(px in unique(Seurat_Obj@meta.data$px)) {
    ### some pre-computed indicies (exclude GMP-redo)
    gmp_carpos_cd4_idx <- intersect(intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                              which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4")),
                                    intersect(which(Seurat_Obj@meta.data$time == "GMP"),
                                              which(Seurat_Obj@meta.data$px == px)))
    gmp_carpos_cd8_idx <- intersect(intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                              which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8")),
                                    intersect(which(Seurat_Obj@meta.data$time == "GMP"),
                                              which(Seurat_Obj@meta.data$px == px)))
    
    ### CD4/CD8 GMP CARpos subsister cells
    ### the numbers can different because we are only considering GMP excluding GMP-redo
    sc_gmp_carpos_cd4_subsister_cell_num <- length(intersect(which(Seurat_Obj@meta.data$Cluster22 == "YES"),
                                                             gmp_carpos_cd4_idx))
    sc_gmp_carpos_cd8_subsister_cell_num <- length(intersect(which(Seurat_Obj@meta.data$Cluster22 == "YES"),
                                                             gmp_carpos_cd8_idx))
    
    ### single cell GMP CAR+ CD4/CD8 total vs CD4/CD8 subsister ratio
    sc_cd4_subsister_ratio <- sc_gmp_carpos_cd4_subsister_cell_num / length(gmp_carpos_cd4_idx)
    sc_cd8_subsister_ratio <- sc_gmp_carpos_cd8_subsister_cell_num / length(gmp_carpos_cd8_idx)
    
    ### fill the table
    result_table2[px, "10x GMP CAR+ CD4 Cell #"] <- length(gmp_carpos_cd4_idx)
    result_table2[px, "10x GMP CAR+ CD8 Cell #"] <- length(gmp_carpos_cd8_idx)
    result_table2[px, "10x GMP CAR+ CD4 Subsisters #"] <- sc_gmp_carpos_cd4_subsister_cell_num
    result_table2[px, "10x GMP CAR+ CD8 Subsisters #"] <- sc_gmp_carpos_cd8_subsister_cell_num
    result_table2[px, "Estimated GMP CD4 Subsisters #"] <- round(gmp_cell_cnt[px] * Real_GMP_CAR_CD4_pcnt2[px] * sc_cd4_subsister_ratio, 2)
    result_table2[px, "Estimated GMP CD8 Subsisters #"] <- round(gmp_cell_cnt[px] * Real_GMP_CAR_CD8_pcnt2[px] * sc_cd8_subsister_ratio, 2)
    result_table2[px, "GMP CD4 Subsister %"] <- round(Real_GMP_CAR_CD4_pcnt2[px] * sc_cd4_subsister_ratio * 100, 5)
    result_table2[px, "GMP CD8 Subsister %"] <- round(Real_GMP_CAR_CD8_pcnt2[px] * sc_cd8_subsister_ratio * 100, 5)
  }
  
  ### calculate some more info
  result_table2[,"Estimated Injected GMP CD4 Subsisters # (per kg)"] <- round(result_table2[,"Adjusted Dose (per Kg)"] * result_table2[,"GMP CD4 Subsister %"] / 100, 0)
  result_table2[,"Estimated Injected GMP CD8 Subsisters # (per kg)"] <- round(result_table2[,"Adjusted Dose (per Kg)"] * result_table2[,"GMP CD8 Subsister %"] / 100, 0)
  
  ### write out the result_table
  write.xlsx2(data.frame(Patient=rownames(result_table2),
                         result_table2,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/Estimated_CD4_CD8_Subsister(Cluster22)_#_In_GMP.xlsx"),
              sheetName = "Estimated_CD4_CD8_Subsister_#_GMP", row.names = FALSE)
  
  ### bar graph
  plot_df <- data.frame(Cell_Num=c(result_table2$`Estimated Injected GMP CD4 Subsisters # (per kg)`,
                                   result_table2$`Estimated Injected GMP CD8 Subsisters # (per kg)`),
                        Patient=c(rownames(result_table2), rownames(result_table2)),
                        Cell_Type=c(rep("CD4", nrow(result_table2)), rep("CD8", nrow(result_table2))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df <- plot_df[which(!is.na(plot_df$Cell_Num)),]
  plot_df$Cell_Num2 <- as.character(plot_df$Cell_Num)
  plot_df$Cell_Num2[which(plot_df$Cell_Num2 == "0")] <- ""
  temp <- plot_df[order(-plot_df$Cell_Num),]
  temp <- temp[which(temp$Cell_Type == "CD8"),]
  plot_df$Patient <- factor(plot_df$Patient, levels = unique(temp$Patient))
  plot_df$Cell_Num[which(plot_df$Cell_Num > 10000)] <- 10000
  p <- ggplot(data=plot_df, aes_string(x="Patient", y="Cell_Num", fill="Cell_Type", label="Cell_Num2")) +
    geom_bar(position = position_dodge(), stat = "identity") +
    ggtitle("Infused GMP Subsisters (Cluster22) # (per kg)") +
    geom_text(size = 6, position = position_dodge(1)) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5, size = 20))
  ggsave(file = paste0(outputDir2, "Estimated_CD4_CD8_Subsister(Cluster22)_#_In_GMP.png"), plot = p, width = 25, height = 10, dpi = 400)
  
  temp_df <- plot_df[which(plot_df$Cell_Type == "CD4"),]
  temp_df <- temp_df[order(-temp_df$Cell_Num),]
  temp_df <- temp_df[which(temp_df$Cell_Num > 0),]
  temp_df$Patient <- factor(temp_df$Patient, levels = unique(temp_df$Patient))
  p <- ggplot(data=temp_df, aes_string(x="Patient", y="Cell_Num", label="Cell_Num2")) +
    geom_bar(position = position_dodge(), stat = "identity", fill="#F8766D") +
    ggtitle("Infused GMP CD4 Subsisters (Cluster22) # (per kg)") +
    geom_text(size = 6, position = position_dodge(1)) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5, size = 20))
  ggsave(file = paste0(outputDir2, "Estimated_CD4_Subsister(Cluster22)_#_In_GMP.png"), plot = p, width = 25, height = 10, dpi = 400)
  
  temp_df <- plot_df[which(plot_df$Cell_Type == "CD8"),]
  temp_df <- temp_df[order(-temp_df$Cell_Num),]
  temp_df <- temp_df[which(temp_df$Cell_Num > 0),]
  temp_df$Patient <- factor(temp_df$Patient, levels = unique(temp_df$Patient))
  p <- ggplot(data=temp_df, aes_string(x="Patient", y="Cell_Num", label="Cell_Num2")) +
    geom_bar(position = position_dodge(), stat = "identity", fill="#00BFC4") +
    ggtitle("Infused GMP CD8 Subsisters (Cluster22) # (per kg)") +
    geom_text(size = 6, position = position_dodge(1)) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = -45, hjust = 0.5, vjust = 0.5, size = 20))
  ggsave(file = paste0(outputDir2, "Estimated_CD8_Subsister(Cluster22)_#_In_GMP.png"), plot = p, width = 25, height = 10, dpi = 400)
  
  ### correlation between adjusted dose & subsister #
  plot_df <- data.frame(Adjusted_Dose=as.numeric(result_table2$`Adjusted Dose (per Kg)`),
                        Injected_Subsisters=as.numeric(result_table2$`Estimated Injected GMP CD8 Subsisters # (per kg)`),
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw the correlation plot
  s_cor <- round(cor(plot_df$Adjusted_Dose, plot_df$Injected_Subsisters, method = "spearman", use = "complete.obs"), 2)
  p <- ggplot(data = plot_df, aes(x=Adjusted_Dose, y=Injected_Subsisters)) +
    geom_point(col = "blue", size = 8) +
    labs(title = paste0("Spearman Correlation:", s_cor)) +
    xlab("Injected Dose (per kg)") +
    ylab("Injected Subsisters (per kg)") +
    geom_smooth(method = lm, color="black", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24))
  ggsave(file = paste0(outputDir2, "Correlation_GMP_Dose_vs_Subsister#.png"), plot = p, width = 18, height = 10, dpi = 400)
  
  
  #
  ### 14. What are the DEGs between cells from the late time (after six weeks) points vs. early time points?
  #   - CARpos CD8 after infusion [After 6 weeks] vs [Early time points; not GMP]
  #   - CARpos CD8 subsisters [After 6 weeks] vs [Early time points; not GMP]
  #   DE analysis and pathway analysis
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/14/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  #
  ### get seurat object for CAR+ CD8 post-infusion time points
  #
  ### get CARpos-only seurat object
  carpos_cd8_cells <- rownames(Seurat_Obj@meta.data)[intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                               which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))]
  sub_seurat_obj5 <- subset(Seurat_Obj, cells = carpos_cd8_cells)
  
  ### after gmp time points only
  pi_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", 
                      "Wk6", "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj5 <- SetIdent(object = sub_seurat_obj5,
                              cells = rownames(sub_seurat_obj5@meta.data),
                              value = sub_seurat_obj5@meta.data$time2)
  sub_seurat_obj5 <- subset(sub_seurat_obj5, idents = intersect(pi_time_points,
                                                                unique(sub_seurat_obj5@meta.data$time2)))
  
  ### get seurat object for some specific patients
  sub_seurat_obj5 <- SetIdent(object = sub_seurat_obj5,
                              cells = rownames(sub_seurat_obj5@meta.data),
                              value = sub_seurat_obj5@meta.data$px)
  sub_seurat_obj5 <- subset(sub_seurat_obj5, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                        "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                        "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### normalization
  sub_seurat_obj5 <- NormalizeData(sub_seurat_obj5,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  sub_seurat_obj5 <- FindVariableFeatures(sub_seurat_obj5,
                                          selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  sub_seurat_obj5 <- ScaleData(sub_seurat_obj5,
                               vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  sub_seurat_obj5 <- RunPCA(sub_seurat_obj5,
                            features = VariableFeatures(object = sub_seurat_obj5),
                            npcs = 15)
  sub_seurat_obj5 <- RunUMAP(sub_seurat_obj5, dims = 1:15)
  
  ### set early vs late
  sub_seurat_obj5@meta.data$early_late <- "Early"
  sub_seurat_obj5@meta.data$early_late[which(sub_seurat_obj5@meta.data$time2 %in% c("Wk8", "3mo", "6mo", "9mo"))] <- "Late"
  
  ### DE analysis
  sub_seurat_obj5 <- SetIdent(object = sub_seurat_obj5,
                              cells = rownames(sub_seurat_obj5@meta.data),
                              value = sub_seurat_obj5@meta.data$early_late)
  de_result <- FindMarkers(sub_seurat_obj5,
                           ident.1 = "Early",
                           ident.2 = "Late",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/DE_PI_CARpos_CD8_Early_vs_Late.xlsx"),
              sheetName = "PI_CARpos_CD8_Early_vs_Late_DE_Result", row.names = FALSE)
  
  ### get entrez ids for the genes
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[which(de_result$p_val_adj < 0.01)],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("Pathways_with_DE_Genes_PI_CARpos_CD8_Early_vs_Late"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir2)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathways_with_DE_Genes_PI_CARpos_CD8_Early_vs_Late"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir2)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir2, "GO_Pathways_with_DE_genes_PI_CARpos_CD8_Early_vs_Late.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, "KEGG_Pathways_with_DE_genes_PI_CARpos_CD8_Early_vs_Late.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  
  ### CARpos CD8 subsisters [After 6 weeks] vs [Early time points; not GMP]
  sub_seurat_obj5@meta.data$early_late2 <- NA
  sub_seurat_obj5@meta.data$early_late2[intersect(which(sub_seurat_obj5@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones),
                                                  which(sub_seurat_obj5@meta.data$early_late == "Early"))] <- "Early_Subsisters"
  sub_seurat_obj5@meta.data$early_late2[intersect(which(sub_seurat_obj5@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones),
                                                  which(sub_seurat_obj5@meta.data$early_late == "Late"))] <- "Late_Subsisters"
  
  ### DE analysis
  sub_seurat_obj5 <- SetIdent(object = sub_seurat_obj5,
                              cells = rownames(sub_seurat_obj5@meta.data),
                              value = sub_seurat_obj5@meta.data$early_late2)
  de_result <- FindMarkers(sub_seurat_obj5,
                           ident.1 = "Early_Subsisters",
                           ident.2 = "Late_Subsisters",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/DE_PI_CARpos_CD8_Subsisters_Early_vs_Late.xlsx"),
              sheetName = "PI_CARpos_CD8_Subsisters_Early_vs_Late_DE_Result", row.names = FALSE)
  
  ### get entrez ids for the genes
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[which(de_result$p_val_adj < 0.01)],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("Pathways_with_DE_Genes_PI_CARpos_CD8_Subsisters_Early_vs_Late"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir2)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathways_with_DE_Genes_PI_CARpos_CD8_Subsisters_Early_vs_Late"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir2)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir2, "GO_Pathways_with_DE_genes_PI_CARpos_CD8_Subsisters_Early_vs_Late.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, "KEGG_Pathways_with_DE_genes_PI_CARpos_CD8_Subsisters_Early_vs_Late.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  #
  ### 15. Feature plots and DE list with specific genes from Tan paper
  #   DE result with the specific genes
  #   & Feature plot of those genes with GMP & PI UMAP
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/15/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### define gene signatures from Tan paper
  tan_poor_persistency_genes <- c("IRF7", "RSAD2", "MX1", "ISG15", "OASL", "IFIT3")
  tan_good_persistency_genes <- c("TCF7", "LEF1", "TBX21", "LDHA", "GLUT1", "GAPDH", "MKI67", "PRDM1", "ZEB2")
  
  ### check whether they are in the data
  tan_poor_persistency_genes <- intersect(tan_poor_persistency_genes,
                                          rownames(Seurat_Obj@assays$RNA@counts))
  tan_good_persistency_genes <- intersect(tan_good_persistency_genes,
                                          rownames(Seurat_Obj@assays$RNA@counts))
  
  #
  ### get DE result of the genes
  #
  ### PB GMP subsister vs non-subsister
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister)
  gmp_de_result <- FindMarkers(target_Seurat_Obj,
                               ident.1 = "YES",
                               ident.2 = "NO",
                               min.pct = 0.1,
                               logfc.threshold = 0.2,
                               test.use = "wilcox")
  ### PB PI subsister vs non-subsister
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$ALL_GMP_CARpos_Persister)
  pi_de_result <- FindMarkers(sub_seurat_obj2,
                              ident.1 = "YES",
                              ident.2 = "NO",
                              min.pct = 0.1,
                              logfc.threshold = 0.2,
                              test.use = "wilcox")
  
  ### only look at the specific genes and save them
  ### GMP
  temp <- gmp_de_result[c(tan_poor_persistency_genes, tan_good_persistency_genes),]
  temp <- temp[which(!is.na(temp$p_val)),]
  write.xlsx2(data.frame(Gene=rownames(temp),
                         temp,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/PB_GMP_CARpos_CD8_S_vs_NS_Tan_DE_Result.xlsx"),
              sheetName = "PB_GMP_CARpos_CD8_DE_Result", row.names = FALSE)
  ### PI
  temp <- pi_de_result[c(tan_poor_persistency_genes, tan_good_persistency_genes),]
  temp <- temp[which(!is.na(temp$p_val)),]
  write.xlsx2(data.frame(Gene=rownames(temp),
                         temp,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/PB_PI_CARpos_CD8_S_vs_NS_Tan_DE_Result.xlsx"),
              sheetName = "PB_PI_CARpos_CD8_DE_Result", row.names = FALSE)
  
  ### Feature plot of the genes
  ### GMP
  p <- FeaturePlot(target_Seurat_Obj, features = tan_poor_persistency_genes, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "PB_GMP_CARpos_CD8_Tan_Poor_Persistency.png"), plot = p, width = 15, height = 10, dpi = 350)
  p <- FeaturePlot(target_Seurat_Obj, features = tan_good_persistency_genes, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "PB_GMP_CARpos_CD8_Tan_Good_Persistency.png"), plot = p, width = 15, height = 10, dpi = 350)
  ### PI
  p <- FeaturePlot(sub_seurat_obj2, features = tan_poor_persistency_genes, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "PB_PI_CARpos_CD8_Tan_Poor_Persistency.png"), plot = p, width = 15, height = 10, dpi = 350)
  p <- FeaturePlot(sub_seurat_obj2, features = tan_good_persistency_genes, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "PB_PI_CARpos_CD8_Tan_Good_Persistency.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  #
  ### 16. Trajectory analysis of differentiation of CD8 CAR+ (start with all time points and perhaps move to only post-infusion depending on how it looks)
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/16/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get CARpos-only seurat object
  carpos_cd8_cells <- rownames(Seurat_Obj@meta.data)[intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                               which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))]
  sub_seurat_obj4 <- subset(Seurat_Obj, cells = carpos_cd8_cells)
  
  ### after gmp time points only
  gmp_after_time_points <- c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", 
                             "Wk6", "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$time2)
  sub_seurat_obj4 <- subset(sub_seurat_obj4, idents = intersect(gmp_after_time_points,
                                                                unique(sub_seurat_obj4@meta.data$time2)))
  
  ### get seurat object for some specific patients
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$px)
  sub_seurat_obj4 <- subset(sub_seurat_obj4, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                        "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                        "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### run mnn
  sub_seurat_obj4.list <- SplitObject(sub_seurat_obj4, split.by = "px")
  sub_seurat_obj4 <- RunFastMNN(object.list = sub_seurat_obj4.list)
  rm(sub_seurat_obj4.list)
  gc()
  
  ### normalization
  sub_seurat_obj4 <- NormalizeData(sub_seurat_obj4,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  sub_seurat_obj4 <- FindVariableFeatures(sub_seurat_obj4,
                                          selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  sub_seurat_obj4 <- ScaleData(sub_seurat_obj4,
                               vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  sub_seurat_obj4 <- RunPCA(sub_seurat_obj4,
                            features = VariableFeatures(object = sub_seurat_obj4),
                            npcs = 15)
  sub_seurat_obj4 <- RunUMAP(sub_seurat_obj4, dims = 1:15)
  
  ### PCA map
  pca_map <- Embeddings(sub_seurat_obj4, reduction = "pca")[rownames(sub_seurat_obj4@meta.data),1:2]
  
  ### get slingshot object
  slingshot_obj <- slingshot(pca_map,
                             clusterLabels = sub_seurat_obj4@meta.data$time2,
                             start.clus = "GMP",
                             end.clus = "6mo")
  slingshot_obj <- as.SlingshotDataSet(slingshot_obj)
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(sub_seurat_obj4@meta.data$time2), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, "CARpos_CD8_Trajectory_Inference_Time_PCA.png"), width = 5000, height = 3000, res = 350)
  par(mar=c(7, 7, 7, 1), mgp=c(4,1,0))
  plot(reducedDim(slingshot_obj),
       main=paste("CAR+ CD8 Trajectory Inference"),
       col = cell_colors_clust[as.character(sub_seurat_obj4@meta.data$time2)],
       pch = 19, cex = 1, cex.lab = 3, cex.main = 3, cex.axis = 2)
  # title(xlab="PC1", mgp=c(1,1,0), cex.lab=3)
  # title(ylab="PC2", mgp=c(1,1,0), cex.lab=3)
  lines(slingshot_obj, lwd = 4, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, cex = 1.5)
  dev.off()
  
  ### MNN map
  mnn_map <- Embeddings(sub_seurat_obj4, reduction = "mnn")[rownames(sub_seurat_obj4@meta.data),1:2]
  
  ### get slingshot object
  slingshot_obj_mnn <- slingshot(mnn_map,
                                 clusterLabels = sub_seurat_obj4@meta.data$time2,
                                 start.clus = "GMP",
                                 end.clus = "6mo")
  slingshot_obj_mnn <- as.SlingshotDataSet(slingshot_obj_mnn)
  
  ### Trajectory inference
  png(paste0(outputDir2, "CARpos_CD8_Trajectory_Inference_Time_MNN.png"), width = 5000, height = 3000, res = 350)
  par(mar=c(7, 7, 7, 1), mgp=c(4,1,0))
  plot(reducedDim(slingshot_obj_mnn),
       main=paste("CAR+ CD8 Trajectory Inference"),
       col = cell_colors_clust[as.character(sub_seurat_obj4@meta.data$time2)],
       pch = 19, cex = 1, cex.lab = 3, cex.main = 3, cex.axis = 2)
  # title(xlab="PC1", mgp=c(1,1,0), cex.lab=3)
  # title(ylab="PC2", mgp=c(1,1,0), cex.lab=3)
  lines(slingshot_obj_mnn, lwd = 4, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, cex = 1.5)
  dev.off()
  
  
  ### cell # for each time point
  cell_num <- sapply(unique(sub_seurat_obj4@meta.data$time2), function(x) {
    return(length(which(sub_seurat_obj4@meta.data$time2 == x)))
  })
  
  ### remove 6mo data since there are only 3 cells
  cell_num <- cell_num[-8]
  
  ### the data is too big (can't run for monocle), so we down-sample the cells in each time point
  set.seed(1234)
  fixed_min_cell_num <- min(cell_num)
  sub_seurat_obj4@meta.data$downsampled <- "NO"
  for(tp in names(cell_num)) {
    sub_seurat_obj4@meta.data$downsampled[sample(which(sub_seurat_obj4@meta.data$time2 == tp), fixed_min_cell_num)] <- "YES"
  }
  
  ### PCA map
  pca_map2 <- Embeddings(sub_seurat_obj4, reduction = "pca")[rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$downsampled == "YES")],1:2]
  
  ### get slingshot object
  slingshot_obj2 <- slingshot(pca_map2,
                              clusterLabels = sub_seurat_obj4@meta.data$time2[which(sub_seurat_obj4@meta.data$downsampled == "YES")],
                              start.clus = "GMP",
                              end.clus = "3mo")
  slingshot_obj2 <- as.SlingshotDataSet(slingshot_obj2)
  
  ### get colors for the clustering result
  cell_colors_clust2 <- cell_pal(names(cell_num), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, "CARpos_CD8_Trajectory_Inference_Time_PCA_Downsampled.png"), width = 5000, height = 3000, res = 350)
  par(mar=c(7, 7, 7, 1), mgp=c(4,1,0))
  plot(reducedDim(slingshot_obj2),
       main=paste("CAR+ CD8 Trajectory Inference"),
       col = cell_colors_clust2[as.character(sub_seurat_obj4@meta.data$time2[which(sub_seurat_obj4@meta.data$downsampled == "YES")])],
       pch = 19, cex = 1, cex.lab = 3, cex.main = 3, cex.axis = 2)
  # title(xlab="PC1", mgp=c(1,1,0), cex.lab=3)
  # title(ylab="PC2", mgp=c(1,1,0), cex.lab=3)
  lines(slingshot_obj2, lwd = 4, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust2)
  legend("bottomleft", legend = names(cell_colors_clust2), col = cell_colors_clust2,
         pch = 19, cex = 1.5)
  dev.off()
  
  ### MNN map
  mnn_map2 <- Embeddings(sub_seurat_obj4, reduction = "mnn")[rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$downsampled == "YES")],1:2]
  
  ### get slingshot object
  slingshot_obj2_mnn <- slingshot(mnn_map2,
                                  clusterLabels = sub_seurat_obj4@meta.data$time2[which(sub_seurat_obj4@meta.data$downsampled == "YES")],
                                  start.clus = "GMP",
                                  end.clus = "3mo")
  slingshot_obj2_mnn <- as.SlingshotDataSet(slingshot_obj2_mnn)
  
  ### Trajectory inference
  png(paste0(outputDir2, "CARpos_CD8_Trajectory_Inference_Time_MNN_Downsampled.png"), width = 5000, height = 3000, res = 350)
  par(mar=c(7, 7, 7, 1), mgp=c(4,1,0))
  plot(reducedDim(slingshot_obj2_mnn),
       main=paste("CAR+ CD8 Trajectory Inference"),
       col = cell_colors_clust2[as.character(sub_seurat_obj4@meta.data$time2[which(sub_seurat_obj4@meta.data$downsampled == "YES")])],
       pch = 19, cex = 1, cex.lab = 3, cex.main = 3, cex.axis = 2)
  # title(xlab="PC1", mgp=c(1,1,0), cex.lab=3)
  # title(ylab="PC2", mgp=c(1,1,0), cex.lab=3)
  lines(slingshot_obj2_mnn, lwd = 4, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust2)
  legend("bottomleft", legend = names(cell_colors_clust2), col = cell_colors_clust2,
         pch = 19, cex = 1.5)
  dev.off()
  
  ### run mnn
  sub_seurat_obj4.list <- SplitObject(sub_seurat_obj4, split.by = "library")
  sub_seurat_obj4 <- RunFastMNN(object.list = sub_seurat_obj4.list)
  rm(sub_seurat_obj4.list)
  gc()
  
  ### normalization
  sub_seurat_obj4 <- NormalizeData(sub_seurat_obj4,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  sub_seurat_obj4 <- FindVariableFeatures(sub_seurat_obj4,
                                          selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  sub_seurat_obj4 <- ScaleData(sub_seurat_obj4,
                               vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  sub_seurat_obj4 <- RunPCA(sub_seurat_obj4,
                            features = VariableFeatures(object = sub_seurat_obj4),
                            npcs = 15)
  sub_seurat_obj4 <- RunUMAP(sub_seurat_obj4, dims = 1:15)
  
  ### MNN map2
  mnn_map <- Embeddings(sub_seurat_obj4, reduction = "mnn")[rownames(sub_seurat_obj4@meta.data),1:2]
  
  ### get slingshot object
  slingshot_obj_mnn <- slingshot(mnn_map,
                                 clusterLabels = sub_seurat_obj4@meta.data$time2,
                                 start.clus = "GMP",
                                 end.clus = "6mo")
  slingshot_obj_mnn <- as.SlingshotDataSet(slingshot_obj_mnn)
  
  ### Trajectory inference
  png(paste0(outputDir2, "CARpos_CD8_Trajectory_Inference_Time_MNN2.png"), width = 5000, height = 3000, res = 350)
  par(mar=c(7, 7, 7, 1), mgp=c(4,1,0))
  plot(reducedDim(slingshot_obj_mnn),
       main=paste("CAR+ CD8 Trajectory Inference"),
       col = cell_colors_clust[as.character(sub_seurat_obj4@meta.data$time2)],
       pch = 19, cex = 1, cex.lab = 3, cex.main = 3, cex.axis = 2)
  # title(xlab="PC1", mgp=c(1,1,0), cex.lab=3)
  # title(ylab="PC2", mgp=c(1,1,0), cex.lab=3)
  lines(slingshot_obj_mnn, lwd = 4, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, cex = 1.5)
  dev.off()
  
  ### MNN map
  mnn_map2 <- Embeddings(sub_seurat_obj4, reduction = "mnn")[rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$downsampled == "YES")],1:2]
  
  ### get slingshot object
  slingshot_obj2_mnn <- slingshot(mnn_map2,
                                  clusterLabels = sub_seurat_obj4@meta.data$time2[which(sub_seurat_obj4@meta.data$downsampled == "YES")],
                                  start.clus = "GMP",
                                  end.clus = "3mo")
  slingshot_obj2_mnn <- as.SlingshotDataSet(slingshot_obj2_mnn)
  
  ### Trajectory inference
  png(paste0(outputDir2, "CARpos_CD8_Trajectory_Inference_Time_MNN2_Downsampled.png"), width = 5000, height = 3000, res = 350)
  par(mar=c(7, 7, 7, 1), mgp=c(4,1,0))
  plot(reducedDim(slingshot_obj2_mnn),
       main=paste("CAR+ CD8 Trajectory Inference"),
       col = cell_colors_clust2[as.character(sub_seurat_obj4@meta.data$time2[which(sub_seurat_obj4@meta.data$downsampled == "YES")])],
       pch = 19, cex = 1, cex.lab = 3, cex.main = 3, cex.axis = 2)
  # title(xlab="PC1", mgp=c(1,1,0), cex.lab=3)
  # title(ylab="PC2", mgp=c(1,1,0), cex.lab=3)
  lines(slingshot_obj2_mnn, lwd = 4, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust2)
  legend("bottomleft", legend = names(cell_colors_clust2), col = cell_colors_clust2,
         pch = 19, cex = 1.5)
  dev.off()
  
  ### see whether the subsisters are still in the same cluster or not
  sub_seurat_obj4 <- RunUMAP(sub_seurat_obj4, reduction = "mnn", dims = 1:15)
  sub_seurat_obj4 <- FindNeighbors(sub_seurat_obj4, reduction = "mnn", dims = 1:15)
  sub_seurat_obj4 <- FindClusters(sub_seurat_obj4)
  
  p <- DimPlot(object = sub_seurat_obj4, reduction = "mnn",
               group.by = "ALL_GMP_CARpos_Persister", split.by = "time2",
               pt.size = 2, ncol = 3,
               cols = c("NO" = "lightgray", "YES" = "red"),
               order = c("YES", "NO")) +
    ggtitle("") +
    labs(color="Is_Subsistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "MNN2_CARpos_CD8_Subsisters.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  p <- DimPlot(object = sub_seurat_obj4, reduction = "umap",
               group.by = "ALL_GMP_CARpos_Persister", split.by = "time2",
               pt.size = 2, ncol = 3,
               cols = c("NO" = "lightgray", "YES" = "red"),
               order = c("YES", "NO")) +
    ggtitle("") +
    labs(color="Is_Subsistent") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir2, "MNN2_UMAP_CARpos_CD8_Subsisters.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  
  ### Construct a monocle cds
  monocle_metadata <- sub_seurat_obj4@meta.data[rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$downsampled == "YES")],]
  monocle_metadata$time2 <- factor(monocle_metadata$time2, levels = unique(monocle_metadata$time2))
  monocle_cds <- newCellDataSet(as(as.matrix(sub_seurat_obj4@assays$RNA@data[,rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$downsampled == "YES")]]), 'sparseMatrix'),
                                phenoData = new('AnnotatedDataFrame', data = monocle_metadata),
                                featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(sub_seurat_obj4@assays$RNA@data),
                                                                                          row.names = row.names(sub_seurat_obj4@assays$RNA@data),
                                                                                          stringsAsFactors = FALSE, check.names = FALSE)),
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  
  ### run monocle
  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)
  monocle_cds <- reduceDimension(monocle_cds, reduction_method = "DDRTree")
  monocle_cds <- orderCells(monocle_cds)
  
  ### determine the beginning state
  plot_cell_trajectory(monocle_cds, color_by = "time2") + geom_point(alpha=0.1)
  plot_cell_trajectory(monocle_cds, color_by = "State")
  plot_complex_cell_trajectory(monocle_cds, color_by = "time2")
  plot_complex_cell_trajectory(monocle_cds, color_by = "State")
  monocle_cds <- orderCells(monocle_cds, root_state = "5")
  
  ### draw monocle plots
  p <- plot_cell_trajectory(monocle_cds, color_by = "time2", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_Time_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_complex_cell_trajectory(monocle_cds, color_by = "time2", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_Time_Complex_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_complex_cell_trajectory(monocle_cds, color_by = "State", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_State_Complex_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  ### gene expression changes
  interesting_genes <- c("SELL", "GZMH", "GZMK", "CCL4", "KLRD1", "IFITM3")
  cds_subset <- monocle_cds[interesting_genes,]
  p <- plot_genes_in_pseudotime(cds_subset, color_by = "time2",
                           cell_size = 2, ncol = 2) +
    labs(color = "Time") +
    theme_classic(base_size = 20) +
    theme(legend.title = element_text(size = 30),
          legend.text = element_text(size = 20))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_GEX.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  ### what are the gene expression differences in GMP State 1 vs State 5?
  sub_seurat_obj4@meta.data$State <- NA
  sub_seurat_obj4@meta.data[rownames(monocle_cds@phenoData@data),"State"] <- monocle_cds@phenoData@data$State
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$State)
  de_result_st1_vs_5 <- FindMarkers(sub_seurat_obj4,
                                    ident.1 = "1",
                                    ident.2 = "5",
                                    min.pct = 0.2,
                                    logfc.threshold = 0.2,
                                    test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result_st1_vs_5),
                         de_result_st1_vs_5,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_CD8_St1_vs_St5.xlsx"),
              sheetName = "CARpos_CD8_St1_vs_St5", row.names = FALSE)
  
  ### get entrez ids for the genes
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result_st1_vs_5)[which(de_result_st1_vs_5$p_val_adj < 0.01)],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("Pathways_with_DE_Genes_St1_vs_St5"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir2)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathways_with_DE_Genes_St1_vs_St5"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir2)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir2, "GO_Pathways_with_DE_genes_St1_vs_St5.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, "KEGG_Pathways_with_DE_genes_St1_vs_St5.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  ### the data is too big (can't run for monocle), so we down-sample the cells in each time point
  ### BUT THIS TIME, INCLUDE ALL THE GMP SUBSISTERS (ADDITION TO THE ORIGINAL DOWNSAMPLED ONES)
  sub_seurat_obj4@meta.data$downsampled2 <- sub_seurat_obj4@meta.data$downsampled
  sub_seurat_obj4@meta.data$downsampled2[which(sub_seurat_obj4@meta.data$GMP_CARpos_CD8_Persister == "YES")] <- "YES"
  
  ### Construct a monocle cds
  monocle_metadata <- sub_seurat_obj4@meta.data[rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$downsampled2 == "YES")],]
  monocle_metadata$time2 <- factor(monocle_metadata$time2, levels = unique(monocle_metadata$time2))
  monocle_cds <- newCellDataSet(as(as.matrix(sub_seurat_obj4@assays$RNA@data[,rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$downsampled2 == "YES")]]), 'sparseMatrix'),
                                phenoData = new('AnnotatedDataFrame', data = monocle_metadata),
                                featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(sub_seurat_obj4@assays$RNA@data),
                                                                                          row.names = row.names(sub_seurat_obj4@assays$RNA@data),
                                                                                          stringsAsFactors = FALSE, check.names = FALSE)),
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  
  ### run monocle
  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)
  monocle_cds <- reduceDimension(monocle_cds, reduction_method = "DDRTree")
  monocle_cds <- orderCells(monocle_cds)
  
  ### determine the beginning state
  plot_cell_trajectory(monocle_cds, color_by = "time2") + geom_point(alpha=0.1)
  plot_cell_trajectory(monocle_cds, color_by = "State")
  plot_complex_cell_trajectory(monocle_cds, color_by = "time2")
  plot_complex_cell_trajectory(monocle_cds, color_by = "State")
  monocle_cds <- orderCells(monocle_cds, root_state = "1")
  
  ### draw monocle plots
  p <- plot_cell_trajectory(monocle_cds, color_by = "time2", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_Time_Monocle2_Subsisters.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_complex_cell_trajectory(monocle_cds, color_by = "time2", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_Time_Complex_Monocle2_Subsisters.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_cell_trajectory(monocle_cds, color_by = "GMP_CARpos_CD8_Persister", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="GMP Subsisters") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_Persistency_Complex_Monocle2_Subsisters.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  ### State1:  241, State 6: 111
  print(length(intersect(which(monocle_cds$State == 1),
                         which(monocle_cds$GMP_CARpos_CD8_Persister == "YES"))))
  print(length(intersect(which(monocle_cds$State == 6),
                         which(monocle_cds$GMP_CARpos_CD8_Persister == "YES"))))
  print(unique(monocle_cds$px[intersect(which(monocle_cds$State == 1),
                                        which(monocle_cds$GMP_CARpos_CD8_Persister == "YES"))]))
  print(unique(monocle_cds$px[intersect(which(monocle_cds$State == 6),
                                        which(monocle_cds$GMP_CARpos_CD8_Persister == "YES"))]))
  
  ### the data is too big (can't run for monocle), so we down-sample the cells in each time point
  ### BUT THIS TIME, INCLUDE ALL THE SUBSISTERS (ADDITION TO THE ORIGINAL DOWNSAMPLED ONES)
  sub_seurat_obj4@meta.data$downsampled3 <- sub_seurat_obj4@meta.data$downsampled2
  sub_seurat_obj4@meta.data$downsampled3[which(sub_seurat_obj4@meta.data$ALL_GMP_CARpos_Persister == "YES")] <- "YES"
  
  ### Construct a monocle cds
  monocle_metadata <- sub_seurat_obj4@meta.data[rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$downsampled3 == "YES")],]
  monocle_metadata$time2 <- factor(monocle_metadata$time2, levels = unique(monocle_metadata$time2))
  monocle_cds <- newCellDataSet(as(as.matrix(sub_seurat_obj4@assays$RNA@data[,rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$downsampled3 == "YES")]]), 'sparseMatrix'),
                                phenoData = new('AnnotatedDataFrame', data = monocle_metadata),
                                featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(sub_seurat_obj4@assays$RNA@data),
                                                                                          row.names = row.names(sub_seurat_obj4@assays$RNA@data),
                                                                                          stringsAsFactors = FALSE, check.names = FALSE)),
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  
  ### run monocle
  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)
  monocle_cds <- reduceDimension(monocle_cds, reduction_method = "DDRTree")
  monocle_cds <- orderCells(monocle_cds)
  
  ### make state consistent with the one before adding subsisters
  monocle_cds@phenoData@data$State <- as.character(monocle_cds@phenoData@data$State)
  monocle_cds@phenoData@data$State[which(monocle_cds@phenoData@data$State == "5")] <- "99"
  monocle_cds@phenoData@data$State[which(monocle_cds@phenoData@data$State == "1")] <- "5"
  monocle_cds@phenoData@data$State[which(monocle_cds@phenoData@data$State == "7")] <- "1"
  monocle_cds@phenoData@data$State[which(monocle_cds@phenoData@data$State == "99")] <- "7"
  monocle_cds@phenoData@data$State <- factor(monocle_cds@phenoData@data$State, levels = unique(monocle_cds@phenoData@data$State[order(as.numeric(monocle_cds@phenoData@data$State))]))
  
  ### determine the beginning state
  plot_cell_trajectory(monocle_cds, color_by = "time2") + geom_point(alpha=0.1)
  plot_cell_trajectory(monocle_cds, color_by = "State")
  plot_complex_cell_trajectory(monocle_cds, color_by = "time2")
  plot_complex_cell_trajectory(monocle_cds, color_by = "State")
  monocle_cds <- orderCells(monocle_cds, root_state = "1")
  
  ### draw monocle plots
  p <- plot_cell_trajectory(monocle_cds, color_by = "time2", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_Time_Monocle2_Subsisters2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_complex_cell_trajectory(monocle_cds, color_by = "time2", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_Time_Complex_Monocle2_Subsisters2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_cell_trajectory(monocle_cds, color_by = "ALL_GMP_CARpos_Persister", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="ALL Subsisters") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CAR_CD8_Trajectory_Inference_Persistency_Complex_Monocle2_Subsisters2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  ### get state info and find all markers for the states
  sub_seurat_obj4@meta.data$State <- NA
  sub_seurat_obj4@meta.data[rownames(monocle_cds@phenoData@data),"State"] <- monocle_cds@phenoData@data$State
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$State)
  de_result_st_all <- FindAllMarkers(sub_seurat_obj4,
                                     min.pct = 0.2,
                                     logfc.threshold = 0.2,
                                     test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result_st_all),
                         de_result_st_all,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_CD8_State_AllMarkers.xlsx"),
              sheetName = "CARpos_CD8_State_AllMarkers", row.names = FALSE)
  
  ### what are the gene expression differences in GMP State 1 vs State 5?
  sub_seurat_obj4@meta.data$State <- NA
  sub_seurat_obj4@meta.data[rownames(monocle_cds@phenoData@data),"State"] <- monocle_cds@phenoData@data$State
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$State)
  de_result_st1_vs_5 <- FindMarkers(sub_seurat_obj4,
                                    ident.1 = "1",
                                    ident.2 = "5",
                                    min.pct = 0.2,
                                    logfc.threshold = 0.2,
                                    test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result_st1_vs_5),
                         de_result_st1_vs_5,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_CD8_St1_vs_St5_Including_Subsisters.xlsx"),
              sheetName = "CARpos_CD8_St1_vs_St5_Including_Subsisters", row.names = FALSE)
  
  ### State1 subsisters vs State5 subsisters
  sub_seurat_obj4@meta.data$state_1_5_subsisters <- NA
  sub_seurat_obj4@meta.data$state_1_5_subsisters[intersect(which(sub_seurat_obj4@meta.data$State == "1"),
                                                           which(sub_seurat_obj4@meta.data$GMP_CARpos_CD8_Persister == "YES"))] <- "State1_Subsisters"
  sub_seurat_obj4@meta.data$state_1_5_subsisters[intersect(which(sub_seurat_obj4@meta.data$State == "5"),
                                                           which(sub_seurat_obj4@meta.data$GMP_CARpos_CD8_Persister == "YES"))] <- "State5_Subsisters"
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$state_1_5_subsisters)
  de_result_st1_vs_5 <- FindMarkers(sub_seurat_obj4,
                                    ident.1 = "State1_Subsisters",
                                    ident.2 = "State5_Subsisters",
                                    min.pct = 0.2,
                                    logfc.threshold = 0.2,
                                    test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result_st1_vs_5),
                         de_result_st1_vs_5,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CARpos_CD8_Subsisters_St1_vs_St5.xlsx"),
              sheetName = "GMP_CARpos_CD8_Subsisters_St1_vs_St5", row.names = FALSE)
  
  #
  ### Proportion of cells in each state - proportional bar graph (time, subsisters, GMP-subsisters, clone_size)
  #
  
  ### time
  plot_df <- data.frame(State=as.vector(sapply(levels(monocle_cds$State), function(x) rep(x, length(levels(monocle_cds$time2))))),
                        Time=rep(levels(monocle_cds$time2), length(levels(monocle_cds$State))),
                        Freq=0,
                        Pcnt=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### fill out the plot_df
  for(i in 1:length(levels(monocle_cds$State))) {
    for(j in 1:length(levels(monocle_cds$time2))) {
      plot_df[(i-1)*length(levels(monocle_cds$time2))+j,"Freq"] <- length(intersect(which(monocle_cds$State == levels(monocle_cds$State)[i]),
                                       which(monocle_cds$time2 == levels(monocle_cds$time2)[j])))
    }
  }
  
  ### remove 0 rows
  plot_df <- plot_df[which(plot_df$Freq != 0),]
  
  ### percentage calculation
  state_sum <- rep(0, length(levels(monocle_cds$State)))
  names(state_sum) <- levels(monocle_cds$State)
  for(i in 1:length(levels(monocle_cds$State))) {
    state_sum[i] <- sum(plot_df[which(plot_df$State == levels(monocle_cds$State)[i]),"Freq"])
    plot_df$Pcnt[which(plot_df$State == levels(monocle_cds$State)[i])] <- round(plot_df$Freq[which(plot_df$State == levels(monocle_cds$State)[i])] * 100 / state_sum[i], 1)
  }
  
  ### factorize the columns
  plot_df$State <- factor(plot_df$State, levels = levels(monocle_cds$State))
  plot_df$Time <- factor(plot_df$Time, levels = levels(monocle_cds$time2))
  
  ### draw a proportional bar plot
  ### pcnt < 10 -> ""
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) < 10)] <- ""
  plot_df$Pcnt <- as.character(plot_df$Pcnt)
  ### state_sum < 20 -> ""
  blank_state <- names(state_sum)[which(state_sum < 20)]
  plot_df$Pcnt[which(plot_df$State %in% blank_state)] <- ""
  p <- ggplot(data=plot_df, aes_string(x="State", y="Freq", fill="Time", label="Pcnt")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Proportion of Cells - Time") +
    geom_text(size = 5, position = position_stack(vjust = 1)) +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "State_Cell_Proportion_Time.png"), plot = p,
         width = 20, height = 10, dpi = 350)
  
  ### subsisters
  plot_df <- data.frame(State=as.vector(sapply(levels(monocle_cds$State), function(x) rep(x, length(unique(monocle_cds$ALL_GMP_CARpos_Persister))))),
                        Subsister=rep(unique(monocle_cds$ALL_GMP_CARpos_Persister), length(levels(monocle_cds$State))),
                        Freq=0,
                        Pcnt=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### fill out the plot_df
  for(i in 1:length(levels(monocle_cds$State))) {
    for(j in 1:length(unique(monocle_cds$ALL_GMP_CARpos_Persister))) {
      plot_df[(i-1)*length(unique(monocle_cds$ALL_GMP_CARpos_Persister))+j,"Freq"] <- length(intersect(which(monocle_cds$State == levels(monocle_cds$State)[i]),
                                                                                    which(monocle_cds$ALL_GMP_CARpos_Persister == unique(monocle_cds$ALL_GMP_CARpos_Persister)[j])))
    }
  }
  
  ### remove 0 rows
  plot_df <- plot_df[which(plot_df$Freq != 0),]
  
  ### percentage calculation
  state_sum <- rep(0, length(levels(monocle_cds$State)))
  names(state_sum) <- levels(monocle_cds$State)
  for(i in 1:length(levels(monocle_cds$State))) {
    state_sum[i] <- sum(plot_df[which(plot_df$State == levels(monocle_cds$State)[i]),"Freq"])
    plot_df$Pcnt[which(plot_df$State == levels(monocle_cds$State)[i])] <- round(plot_df$Freq[which(plot_df$State == levels(monocle_cds$State)[i])] * 100 / state_sum[i], 1)
  }
  
  ### factorize the columns
  plot_df$State <- factor(plot_df$State, levels = levels(monocle_cds$State))
  plot_df$Subsister <- factor(plot_df$Subsister, levels = unique(plot_df$Subsister))
  
  ### draw a proportional bar plot
  ### pcnt < 10 -> ""
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) < 10)] <- ""
  plot_df$Pcnt <- as.character(plot_df$Pcnt)
  ### state_sum < 20 -> ""
  blank_state <- names(state_sum)[which(state_sum < 20)]
  plot_df$Pcnt[which(plot_df$State %in% blank_state)] <- ""
  p <- ggplot(data=plot_df, aes_string(x="State", y="Freq", fill="Subsister", label="Pcnt")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Proportion of Cells - Subsister Group") +
    geom_text(size = 5, position = position_stack(vjust = 1)) +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "State_Cell_Proportion_Subsister_Group.png"), plot = p,
         width = 20, height = 10, dpi = 350)
  
  ### GMP subsisters
  plot_df <- data.frame(State=as.vector(sapply(levels(monocle_cds$State), function(x) rep(x, length(unique(monocle_cds$GMP_CARpos_CD8_Persister))))),
                        Subsister=rep(unique(monocle_cds$GMP_CARpos_CD8_Persister), length(levels(monocle_cds$State))),
                        Freq=0,
                        Pcnt=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### fill out the plot_df
  for(i in 1:length(levels(monocle_cds$State))) {
    for(j in 1:length(unique(monocle_cds$GMP_CARpos_CD8_Persister))) {
      plot_df[(i-1)*length(unique(monocle_cds$GMP_CARpos_CD8_Persister))+j,"Freq"] <- length(intersect(which(monocle_cds$State == levels(monocle_cds$State)[i]),
                                                                                                       which(monocle_cds$GMP_CARpos_CD8_Persister == unique(monocle_cds$GMP_CARpos_CD8_Persister)[j])))
    }
  }
  
  ### remove 0 rows
  plot_df <- plot_df[which(plot_df$Freq != 0),]
  
  ### percentage calculation
  state_sum <- rep(0, length(levels(monocle_cds$State)))
  names(state_sum) <- levels(monocle_cds$State)
  for(i in 1:length(levels(monocle_cds$State))) {
    state_sum[i] <- sum(plot_df[which(plot_df$State == levels(monocle_cds$State)[i]),"Freq"])
    plot_df$Pcnt[which(plot_df$State == levels(monocle_cds$State)[i])] <- round(plot_df$Freq[which(plot_df$State == levels(monocle_cds$State)[i])] * 100 / state_sum[i], 1)
  }
  
  ### factorize the columns
  plot_df$State <- factor(plot_df$State, levels = levels(monocle_cds$State))
  plot_df$Subsister <- factor(plot_df$Subsister, levels = unique(plot_df$Subsister))
  
  ### draw a proportional bar plot
  ### pcnt < 10 -> ""
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) < 10)] <- ""
  plot_df$Pcnt <- as.character(plot_df$Pcnt)
  ### state_sum < 20 -> ""
  blank_state <- names(state_sum)[which(state_sum < 20)]
  plot_df$Pcnt[which(plot_df$State %in% blank_state)] <- ""
  p <- ggplot(data=plot_df, aes_string(x="State", y="Freq", fill="Subsister", label="Pcnt")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Proportion of Cells - GMP Subsisters") +
    geom_text(size = 15, position = position_stack(vjust = 1)) +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "State_Cell_Proportion_GMP_Subsisters.png"), plot = p,
         width = 20, height = 10, dpi = 350)
  

  ### save CAR+ clone size info
  SJCAR19_Clonotype_Frequency <- vector("list", length = length(unique(Seurat_Obj@meta.data$px)))
  names(SJCAR19_Clonotype_Frequency) <- unique(Seurat_Obj@meta.data$px)
  type <- "One_From_Each"
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    
    ### print progress
    writeLines(paste(patient))
    
    ### load the file
    target_file <- read.xlsx2(file = paste0(px_result_dir, patient, "/car_clonotype_frequency_over_time_", patient, ".xlsx"),
                              sheetName = paste0("CARpos_Clonotype_Frequency_", substr(type, 1, 4)), stringsAsFactors = FALSE, check.names = FALSE,
                              row.names = 1)
    
    ### numerize the table
    for(i in 1:ncol(target_file)) {
      target_file[,i] <- as.numeric(target_file[,i])
    }
    
    ### combine some redundant time points to one
    if(length(which(colnames(target_file) == "GMP-redo")) > 0) {
      target_file[,"GMP"] <- target_file[,"GMP"] + target_file[,"GMP-redo"]
      target_file <- target_file[,-which(colnames(target_file) == "GMP-redo")]
    }
    if(length(which(colnames(target_file) == "PreTransB")) > 0) {
      if(length(which(colnames(target_file) == "PreTrans")) > 0) {
        target_file[,"PreTrans"] <- target_file[,"PreTrans"] + target_file[,"PreTransB"]
        target_file <- target_file[,-which(colnames(target_file) == "PreTransB")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "PreTransB")] <- "PreTrans"
      }
    }
    if(length(which(colnames(target_file) == "Wk1b")) > 0) {
      if(length(which(colnames(target_file) == "Wk1")) > 0) {
        target_file[,"Wk1"] <- target_file[,"Wk1"] + target_file[,"Wk1b"]
        target_file <- target_file[,-which(colnames(target_file) == "Wk1b")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "Wk1b")] <- "Wk1"
      }
    }
    if(length(which(colnames(target_file) == "Wk-1Run1")) > 0) {
      if(length(which(colnames(target_file) == "Wk-1")) > 0) {
        target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run1"]
        target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run1")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "Wk-1Run1")] <- "Wk-1"
      }
    }
    if(length(which(colnames(target_file) == "Wk-1Run2")) > 0) {
      if(length(which(colnames(target_file) == "Wk-1")) > 0) {
        target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run2"]
        target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run2")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "Wk-1Run2")] <- "Wk-1"
      }
    }
    
    SJCAR19_Clonotype_Frequency[[patient]] <- target_file
    
  }
  
  ### CLONE SIZE
  Seurat_Obj@meta.data$car_clone_size <- 0
  Seurat_Obj@meta.data$car_clone_size[intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                which(!is.na(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta)))] <- 1
  for(patient in names(SJCAR19_Clonotype_Frequency)) {
    ### per patient per clone give clone size info
    for(cln in rownames(SJCAR19_Clonotype_Frequency[[patient]])) {
      Seurat_Obj@meta.data$car_clone_size[which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta == cln)] <- as.numeric(SJCAR19_Clonotype_Frequency[[patient]][cln,"Total"])
    }
  }
  
  ### plot_df - clone size
  monocle_cds$car_clone_size <- Seurat_Obj@meta.data[rownames(monocle_cds@phenoData@data), "car_clone_size"]
  plot_df <- data.frame(State=as.vector(sapply(levels(monocle_cds$State), function(x) rep(x, length(unique(monocle_cds$car_clone_size))))),
                        Clone_Size=rep(unique(monocle_cds$car_clone_size), length(levels(monocle_cds$State))),
                        Freq=0,
                        Pcnt=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### fill out the plot_df
  for(i in 1:length(levels(monocle_cds$State))) {
    for(j in 1:length(unique(monocle_cds$car_clone_size))) {
      plot_df[(i-1)*length(unique(monocle_cds$car_clone_size))+j,"Freq"] <- length(intersect(which(monocle_cds$State == levels(monocle_cds$State)[i]),
                                                                                                       which(monocle_cds$car_clone_size == unique(monocle_cds$car_clone_size)[j])))
    }
  }
  
  ### remove 0 rows
  plot_df <- plot_df[which(plot_df$Freq != 0),]
  
  ### percentage calculation
  state_sum <- rep(0, length(levels(monocle_cds$State)))
  names(state_sum) <- levels(monocle_cds$State)
  for(i in 1:length(levels(monocle_cds$State))) {
    state_sum[i] <- sum(plot_df[which(plot_df$State == levels(monocle_cds$State)[i]),"Freq"])
    plot_df$Pcnt[which(plot_df$State == levels(monocle_cds$State)[i])] <- round(plot_df$Freq[which(plot_df$State == levels(monocle_cds$State)[i])] * 100 / state_sum[i], 1)
  }
  
  ### factorize the columns
  plot_df$State <- factor(plot_df$State, levels = levels(monocle_cds$State))
  plot_df$Clone_Size <- factor(plot_df$Clone_Size, levels = unique(plot_df$Clone_Size)[order(unique(plot_df$Clone_Size))])
  
  ### draw a proportional bar plot
  ### pcnt < 10 -> ""
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) < 10)] <- ""
  plot_df$Pcnt <- as.character(plot_df$Pcnt)
  ### state_sum < 20 -> ""
  blank_state <- names(state_sum)[which(state_sum < 20)]
  plot_df$Pcnt[which(plot_df$State %in% blank_state)] <- ""
  p <- ggplot(data=plot_df, aes_string(x="State", y="Freq", fill="Clone_Size", label="Pcnt")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Proportion of Cells - Clone_Size") +
    geom_text(size = 5, position = position_stack(vjust = 1)) +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "State_Cell_Proportion_Clone_Size.png"), plot = p,
         width = 20, height = 10, dpi = 350)
  
  
  ### add state info to the target_seurat_obj - GMP CAR+ CD8
  target_Seurat_Obj@meta.data$State <- NA
  target_Seurat_Obj@meta.data[rownames(monocle_cds@phenoData@data),"State"] <- monocle_cds@phenoData@data$State
  
  
  #
  ### 17. Comaprison to a reference data set of well characterized t cell differentiation gene sets
  #       T cell differentiation: GO:0030217 (Gene Ontology)
  #       T cell exhaustion:  GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN (msigdb)
  #
  #       CD8 effector vs. memory T cells (downregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN.html
  #       CD8 effector vs. memory T cells (upregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP.html
  #       CD8 naive vs. effector T cells (downregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_DN.html
  #       CD8 naive vs. effector T cells (upregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP.html
  #       CD8 central vs. effector memory T cells (downregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GSE23321_CENTRAL_VS_EFFECTOR_MEMORY_CD8_TCELL_DN.html
  #       CD8 central vs. effector memory T cells (upregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GSE23321_CENTRAL_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.html
  #       CD8 stem cell vs. central memory T cells (downregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GSE23321_CD8_STEM_CELL_MEMORY_VS_CENTRAL_MEMORY_CD8_TCELL_DN.html
  #       CD8 stem cell vs. central memory T cells (upregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GSE23321_CD8_STEM_CELL_MEMORY_VS_CENTRAL_MEMORY_CD8_TCELL_UP.html
  #       CD8 stem cell vs. effector memory T cells (downregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_DN.html
  #       CD8 stem cell vs. effector memory T cells (upregulated genes): https://www.gsea-msigdb.org/gsea/msigdb/cards/GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP.html
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/17/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get t cell differentiation genes
  go.hs <- go.gsets(species="human")
  t_cell_diff_genes <- go.hs[["go.sets"]][["GO:0030217 T cell differentiation"]]
  t_cell_diff_genes <- mapIds(org.Hs.eg.db,
                              t_cell_diff_genes,
                              "SYMBOL", "ENTREZID")
  t_cell_diff_genes <- t_cell_diff_genes[!is.na(t_cell_diff_genes)]
  
  ### get t cell exhaustion genes
  m_df <- msigdbr(species = "Homo sapiens")
  m_list <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  t_cell_exhaustion_genes <- m_list[["GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN"]]
  
  ### get t_cell_diff_list based on Tay's suggestion
  t_cell_diff_list <- list(m_list[["GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN"]],
                           m_list[["GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP"]],
                           m_list[["GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_DN"]],
                           m_list[["GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP"]],
                           m_list[["GSE23321_CENTRAL_VS_EFFECTOR_MEMORY_CD8_TCELL_DN"]],
                           m_list[["GSE23321_CENTRAL_VS_EFFECTOR_MEMORY_CD8_TCELL_UP"]],
                           m_list[["GSE23321_CD8_STEM_CELL_MEMORY_VS_CENTRAL_MEMORY_CD8_TCELL_DN"]],
                           m_list[["GSE23321_CD8_STEM_CELL_MEMORY_VS_CENTRAL_MEMORY_CD8_TCELL_UP"]],
                           m_list[["GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_DN"]],
                           m_list[["GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_UP"]])
  names(t_cell_diff_list) <- c("CD8 effector vs. memory T cells (downregulated genes)",
                               "CD8 effector vs. memory T cells (upregulated genes)",
                               "CD8 naive vs. effector T cells (downregulated genes)",
                               "CD8 naive vs. effector T cells (upregulated genes)",
                               "CD8 central vs. effector memory T cells (downregulated genes)",
                               "CD8 central vs. effector memory T cells (upregulated genes)",
                               "CD8 stem cell vs. central memory T cells (downregulated genes)",
                               "CD8 stem cell vs. central memory T cells (upregulated genes)",
                               "CD8 stem cell vs. effector memory T cells (downregulated genes)",
                               "CD8 stem cell vs. effector memory T cells (upregulated genes)")
  
  #
  ### GMP CARpos CD8 subsister vs non-subsister
  #
  ### DE analysis
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister)
  de_result <- FindMarkers(target_Seurat_Obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### subset
  de_result_diff <- de_result[intersect(rownames(de_result), t_cell_diff_genes),]
  de_result_exhaustion <- de_result[intersect(rownames(de_result), t_cell_exhaustion_genes),]
  
  ### write out the DE results
  write.xlsx2(data.frame(Gene=rownames(de_result_diff),
                         de_result_diff,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CARpos_CD8_Subsisters_vs_Non-subsisters_diff.xlsx"),
              sheetName = "GMP_CARpos_CD8_Subsisters_vs_Non-subsisters_diff", row.names = FALSE)
  
  write.xlsx2(data.frame(Gene=rownames(de_result_exhaustion),
                         de_result_exhaustion,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CARpos_CD8_Subsisters_vs_Non-subsisters_exhaustion.xlsx"),
              sheetName = "GMP_CARpos_CD8_Subsisters_vs_Non-subsisters_exhaustion", row.names = FALSE)
  
  ### dot plot
  p <- DotPlot(target_Seurat_Obj,
               features = rownames(de_result_diff),
               cols = c("blue", "red"),
               group.by = "GMP_CARpos_CD8_Persister") +
    scale_size(range = c(2, 15)) +
    coord_flip() +
    xlab("T Cell Differentiation Genes") +
    ylab("Is_Subsistent") +
    theme_calc(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir2, "GMP_CARpos_CD8_Subsisters_vs_Non-subsisters_diff.png"),
         plot = p, width = 15, height = 15, dpi = 350)
  
  p <- DotPlot(target_Seurat_Obj,
               features = rownames(de_result_exhaustion),
               cols = c("blue", "red"),
               group.by = "GMP_CARpos_CD8_Persister") +
    scale_size(range = c(2, 15)) +
    coord_flip() +
    xlab("T Cell Exhaustion Genes") +
    ylab("Is_Subsistent") +
    theme_calc(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir2, "GMP_CARpos_CD8_Subsisters_vs_Non-subsisters_exhaustion.png"),
         plot = p, width = 15, height = 15, dpi = 350)
  
  ### signature preparation
  # signat <- sign(de_result$avg_log2FC) * -log10(de_result$p_val_adj)
  # names(signat) <- rownames(de_result)
  subsister_idx <- which(target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")
  non_subsister_idx <- which(target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "NO")
  ### add 0.01 to avoid Inf, -Inf, & NaN
  s_mat <- data.frame(target_Seurat_Obj@assays$RNA@counts[,subsister_idx]) + 0.01
  ns_mat <- data.frame(target_Seurat_Obj@assays$RNA@counts[,non_subsister_idx]) + 0.01
  s_mat <- apply(s_mat, 1, mean)
  ns_mat <- apply(ns_mat, 1, mean)
  signat <- log2(s_mat / ns_mat)
  
  ### run GSEA
  GSEA_result <- run_gsea(gene_list = t_cell_diff_list, signature = list(signat),
                          fdr_cutoff = 0.05,
                          printPlot = TRUE, printPath = outputDir2)
  
  ### write out the result file
  write.xlsx2(GSEA_result, file = paste0(outputDir2, "/GSEA_GMP_CARpos_CD8_Subsisters_vs_Non-subsisters.xlsx"),
              sheetName = "GSEA_GMP_CARpos_CD8_Subsisters_vs_Non-subsisters", row.names = FALSE)
  
  
  #
  ### CARpos post-infusion subsister vs non-subsister
  #
  ### DE analysis
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$CD8_Persisters)
  de_result <- FindMarkers(sub_seurat_obj2,
                           ident.1 = "Subsisters",
                           ident.2 = "Non-subsisters",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### subset
  de_result_diff <- de_result[intersect(rownames(de_result), t_cell_diff_genes),]
  de_result_exhaustion <- de_result[intersect(rownames(de_result), t_cell_exhaustion_genes),]
  
  ### write out the DE results
  write.xlsx2(data.frame(Gene=rownames(de_result_diff),
                         de_result_diff,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/PI_CARpos_Subsisters_vs_Non-subsisters_diff.xlsx"),
              sheetName = "PI_CARpos_Subsisters_vs_Non-subsisters_diff", row.names = FALSE)
  
  write.xlsx2(data.frame(Gene=rownames(de_result_exhaustion),
                         de_result_exhaustion,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/PI_CARpos_Subsisters_vs_Non-subsisters_exhaustion.xlsx"),
              sheetName = "PI_CARpos_Subsisters_vs_Non-subsisters_exhaustion", row.names = FALSE)
  
  ### dot plot
  p <- DotPlot(sub_seurat_obj2,
          features = rownames(de_result_diff),
          cols = c("blue", "red"),
          group.by = "CD8_Persisters") +
    scale_size(range = c(2, 15)) +
    coord_flip() +
    xlab("T Cell Differentiation Genes") +
    ylab("Is_Subsistent") +
    theme_calc(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir2, "PI_CARpos_Subsisters_vs_Non-subsisters_diff.png"),
         plot = p, width = 15, height = 15, dpi = 350)
  
  p <- DotPlot(sub_seurat_obj2,
               features = rownames(de_result_exhaustion),
               cols = c("blue", "red"),
               group.by = "CD8_Persisters") +
    scale_size(range = c(2, 15)) +
    coord_flip() +
    xlab("T Cell Exhaustion Genes") +
    ylab("Is_Subsistent") +
    theme_calc(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir2, "PI_CARpos_Subsisters_vs_Non-subsisters_exhaustion.png"),
         plot = p, width = 15, height = 15, dpi = 350)
  
  ### signature preparation
  # signat <- sign(de_result$avg_log2FC) * -log10(de_result$p_val_adj)
  # names(signat) <- rownames(de_result)
  subsister_idx <- which(sub_seurat_obj2@meta.data$CD8_Persisters == "Subsisters")
  non_subsister_idx <- which(sub_seurat_obj2@meta.data$CD8_Persisters == "Non-subsisters")
  ### add 0.01 to avoid Inf, -Inf, & NaN
  s_mat <- data.frame(sub_seurat_obj2@assays$RNA@counts[,subsister_idx]) + 0.01
  ns_mat <- data.frame(sub_seurat_obj2@assays$RNA@counts[,non_subsister_idx]) + 0.01
  s_mat <- apply(s_mat, 1, mean)
  ns_mat <- apply(ns_mat, 1, mean)
  signat <- log2(s_mat / ns_mat)
  
  ### run GSEA
  GSEA_result <- run_gsea(gene_list = t_cell_diff_list, signature = list(signat),
                          fdr_cutoff = 0.05,
                          printPlot = TRUE, printPath = outputDir2)
  
  ### write out the result file
  write.xlsx2(GSEA_result, file = paste0(outputDir2, "/PI_CARpos_Subsisters_vs_Non-subsisters.xlsx"),
              sheetName = "PI_CARpos_Subsisters_vs_Non-subsisters", row.names = FALSE)
  
  
  #
  ### 18. GMP CAR+ CD8 cells with high CAR expression vs low CAR expression (DE, pathway, & subsister difference)
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/18/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### GMP
  ### CAR EXP < 2: 4208, CAR EXP > 50: 3677
  car_cnt <- target_Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",]
  target_Seurat_Obj@meta.data$CAR2 <- "MID"
  target_Seurat_Obj@meta.data$CAR2[which(car_cnt[rownames(target_Seurat_Obj@meta.data)] < 2)] <- "LOW"
  target_Seurat_Obj@meta.data$CAR2[which(car_cnt[rownames(target_Seurat_Obj@meta.data)] > 50)] <- "HIGH"
  
  ### UMAP with time with all patients
  p <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
               group.by = "CAR2",
               pt.size = 2,
               cols = c("HIGH" = "red", "LOW" = "blue", "MID" = "lightgray"),
               order = c("HIGH", "LOW", "MID")) +
    ggtitle("") +
    labs(color="CAR EXP") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "UMAP_GMP_CAR+_CD8_CAR_EXP.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  ### DE analysis
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$CAR2)
  de_result <- FindMarkers(target_Seurat_Obj,
                           ident.1 = "HIGH",
                           ident.2 = "LOW",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE results
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CARpos_CD8_CAR_EXP_HIGH_vs_Low.xlsx"),
              sheetName = "GMP_CARpos_CD8_CAR_EXP_HIGH_vs_Low", row.names = FALSE)
  
  ### get entrez ids for the genes
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[which(de_result$p_val_adj < 0.01)],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("Pathways_with_DE_Genes_GMP_CARpos_CD8_CAR_EXP_HIGH_vs_Low"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir2)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathways_with_DE_Genes_GMP_CARpos_CD8_CAR_EXP_HIGH_vs_Low"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir2)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir2, "GO_Pathways_with_DE_genes_GMP_CARpos_CD8_CAR_EXP_HIGH_vs_Low.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, "KEGG_Pathways_with_DE_genes_GMP_CARpos_CD8_CAR_EXP_HIGH_vs_Low.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  ### Ridge plot
  features <- rownames(de_result)[1:25]
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$CAR2)
  
  ### CAR2
  p <- RidgePlot(target_Seurat_Obj, features = features, ncol = 5,
                 cols = c("blue", "orange", "red"))
  for(i in 1:25) {
    p[[i]] <- p[[i]] + labs(y = "CAR EXP")
  }
  ggsave(paste0(outputDir2, "GMP_CARpos_CD8_CAR_EXP_DEGs.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  features <- rownames(de_result)[1:9]
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$gmp_seurat_clusters)
  cluster_col <- rainbow(length(levels(target_Seurat_Obj@meta.data$gmp_seurat_clusters)))
  names(cluster_col) <- levels(target_Seurat_Obj@meta.data$gmp_seurat_clusters)
  cluster_col["22"] <- "black"
  ### Clusters
  p <- RidgePlot(target_Seurat_Obj, features = features, cols = cluster_col, ncol = 3)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "GMP_CARpos_CD8_CAR_EXP_DEGs_Clusters.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### Time2
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$time2)
  p <- RidgePlot(target_Seurat_Obj, features = features, ncol = 3)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + labs(y = "Time")
  }
  ggsave(paste0(outputDir2, "GMP_CARpos_CD8_CAR_EXP_DEGs_Time.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### State
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$State)
  p <- RidgePlot(target_Seurat_Obj, features = features, ncol = 3)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + labs(y = "State")
  }
  ggsave(paste0(outputDir2, "GMP_CARpos_CD8_CAR_EXP_DEGs_State.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  
  ### Post infusion
  ### CAR EXP < 2: 8515, CAR EXP > 30: 1443
  car_cnt <- sub_seurat_obj2@assays$RNA@counts["JCC-SJCAR19short",]
  sub_seurat_obj2@meta.data$CAR2 <- "MID"
  sub_seurat_obj2@meta.data$CAR2[which(car_cnt[rownames(sub_seurat_obj2@meta.data)] < 2)] <- "LOW"
  sub_seurat_obj2@meta.data$CAR2[which(car_cnt[rownames(sub_seurat_obj2@meta.data)] > 30)] <- "HIGH"
  
  ### UMAP with time with all patients
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "CAR2",
               pt.size = 2,
               cols = c("HIGH" = "red", "LOW" = "blue", "MID" = "lightgray"),
               order = c("HIGH", "LOW", "MID")) +
    ggtitle("") +
    labs(color="CAR EXP") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "UMAP_PI_CAR+_CAR_EXP.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  ### DE analysis
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$CAR2)
  de_result <- FindMarkers(sub_seurat_obj2,
                           ident.1 = "HIGH",
                           ident.2 = "LOW",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE results
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/PI_CARpos_CAR_EXP_HIGH_vs_Low.xlsx"),
              sheetName = "PI_CARpos_CAR_EXP_HIGH_vs_Low", row.names = FALSE)
  
  ### get entrez ids for the genes
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[which(de_result$p_val_adj < 0.01)],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("Pathways_with_DE_Genes_PI_CARpos_CAR_EXP_HIGH_vs_Low"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir2)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathways_with_DE_Genes_PI_CARpos_CAR_EXP_HIGH_vs_Low"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir2)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir2, "GO_Pathways_with_DE_genes_PI_CARpos_CAR_EXP_HIGH_vs_Low.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, "KEGG_Pathways_with_DE_genes_PI_CARpos_CAR_EXP_HIGH_vs_Low.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  ### Ridge plot
  features <- rownames(de_result)[1:25]
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$CAR2)
  p <- RidgePlot(sub_seurat_obj2, features = features, ncol = 5,
                 cols = c("blue", "orange", "red"))
  for(i in 1:25) {
    p[[i]] <- p[[i]] + labs(y = "CAR EXP")
  }
  ggsave(paste0(outputDir2, "PI_CARpos_CAR_EXP_DEGs.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  features <- rownames(de_result)[1:9]
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$clusters)
  cluster_col <- rainbow(length(levels(sub_seurat_obj2@meta.data$clusters)))
  names(cluster_col) <- levels(sub_seurat_obj2@meta.data$clusters)
  cluster_col["0"] <- "black"
  p <- RidgePlot(sub_seurat_obj2, features = features, cols = cluster_col, ncol = 3)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "PI_CARpos_CAR_EXP_DEGs_Clusters.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### Time2
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$time2)
  p <- RidgePlot(sub_seurat_obj2, features = features, ncol = 3)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + labs(y = "Time")
  }
  ggsave(paste0(outputDir2, "PI_CARpos_CAR_EXP_DEGs_Time.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### State
  common_cell_names <- intersect(rownames(sub_seurat_obj2@meta.data), rownames(sub_seurat_obj4@meta.data))
  sub_seurat_obj2@meta.data$State <- NA
  sub_seurat_obj2@meta.data[common_cell_names,"State"] <- sub_seurat_obj4@meta.data[common_cell_names,"State"]
  sub_seurat_obj2@meta.data$State <- factor(sub_seurat_obj2@meta.data$State, levels = unique(sub_seurat_obj2@meta.data$State[order(sub_seurat_obj2@meta.data$State)]))
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$State)
  p <- RidgePlot(sub_seurat_obj2, features = features, ncol = 3)
  for(i in 1:9) {
    p[[i]] <- p[[i]] + labs(y = "State")
  }
  ggsave(paste0(outputDir2, "PI_CARpos_CAR_EXP_DEGs_State.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  
  #
  ### 19. Look at all GMP (CAR+ & CAR-) and all CAR+ GMP cells whether they are two separated clusters
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/19/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### subset with GMP cells only (including both CAR+ & CAR-)
  target_Seurat_Obj3 <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[which(Seurat_Obj@meta.data$time2 == "GMP")])
  
  ### only the patients we are interested
  target_Seurat_Obj3 <- SetIdent(object = target_Seurat_Obj3,
                                 cells = rownames(target_Seurat_Obj3@meta.data),
                                 value = target_Seurat_Obj3@meta.data$px)
  target_Seurat_Obj3 <- subset(target_Seurat_Obj3, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                              "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                              "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### normalization
  target_Seurat_Obj3 <- NormalizeData(target_Seurat_Obj3,
                                      normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj3 <- FindVariableFeatures(target_Seurat_Obj3,
                                             selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj3 <- ScaleData(target_Seurat_Obj3,
                                  vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj3 <- RunPCA(target_Seurat_Obj3,
                               features = VariableFeatures(object = target_Seurat_Obj3),
                               npcs = 15)
  target_Seurat_Obj3 <- RunUMAP(target_Seurat_Obj3, dims = 1:15)
  
  ### factorize the CAR
  target_Seurat_Obj3@meta.data$CAR <- factor(target_Seurat_Obj3@meta.data$CAR, levels = c("CARpos", "CARneg"))
  
  ### draw a UMAP plot
  p <- DimPlot(object = target_Seurat_Obj3, reduction = "umap",
          group.by = "CAR",
          cols = c("CARpos" = "blue", "CARneg" = "red"),
          order = c("CARpos", "CARneg"),
          pt.size = 1) +
    ggtitle("") +
    labs(color="CAR") +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "UMAP_GMP_CAR1.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- DimPlot(object = target_Seurat_Obj3, reduction = "umap",
               group.by = "CAR",
               cols = c("CARpos" = "blue", "CARneg" = "red"),
               order = c("CARneg", "CARpos"),
               pt.size = 1) +
    ggtitle("") +
    labs(color="CAR") +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "UMAP_GMP_CAR2.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### GMP CAR+
  target_Seurat_Obj3 <- subset(target_Seurat_Obj3, cells = rownames(target_Seurat_Obj3@meta.data)[which(target_Seurat_Obj3@meta.data$CAR == "CARpos")])
  
  ### normalization
  target_Seurat_Obj3 <- NormalizeData(target_Seurat_Obj3,
                                      normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj3 <- FindVariableFeatures(target_Seurat_Obj3,
                                             selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj3 <- ScaleData(target_Seurat_Obj3,
                                  vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj3 <- RunPCA(target_Seurat_Obj3,
                               features = VariableFeatures(object = target_Seurat_Obj3),
                               npcs = 15)
  target_Seurat_Obj3 <- RunUMAP(target_Seurat_Obj3, dims = 1:15)
  
  ### draw a UMAP plot
  p <- DimPlot(object = target_Seurat_Obj3, reduction = "umap",
               group.by = "time",
               cols = c("GMP" = "blue", "GMP-redo" = "skyblue"),
               order = c("GMP", "GMP-redo"),
               pt.size = 1) +
    ggtitle("") +
    labs(color="Time") +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "UMAP_GMP_CARpos_redo1.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- DimPlot(object = target_Seurat_Obj3, reduction = "umap",
               group.by = "time",
               cols = c("GMP" = "blue", "GMP-redo" = "skyblue"),
               order = c("GMP-redo", "GMP"),
               pt.size = 1) +
    ggtitle("") +
    labs(color="Time") +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "UMAP_GMP_CARpos_redo2.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### GMP cells only (not GMP-redo)
  target_Seurat_Obj3 <- subset(target_Seurat_Obj3, cells = rownames(target_Seurat_Obj3@meta.data)[which(target_Seurat_Obj3@meta.data$time == "GMP")])
  
  ### normalization
  target_Seurat_Obj3 <- NormalizeData(target_Seurat_Obj3,
                                      normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj3 <- FindVariableFeatures(target_Seurat_Obj3,
                                             selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj3 <- ScaleData(target_Seurat_Obj3,
                                  vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj3 <- RunPCA(target_Seurat_Obj3,
                               features = VariableFeatures(object = target_Seurat_Obj3),
                               npcs = 15)
  target_Seurat_Obj3 <- RunUMAP(target_Seurat_Obj3, dims = 1:15)
  
  ### draw a UMAP plot
  p <- DimPlot(object = target_Seurat_Obj3, reduction = "umap",
               group.by = "px", pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Patient") +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "UMAP_GMP_Only_CARpos_Px.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### draw a PCA plot
  p <- DimPlot(object = target_Seurat_Obj3, reduction = "pca",
               group.by = "px", pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Patient") +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "PCA_GMP_Only_CARpos_Px.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  #
  ### 20. Tay's request to look at some genes of Px11
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/20/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### genes of interest
  genes_of_interest <- c("TIGIT", "KLRD1", "CD86", "IL2RA", "CD70",
                         "LAG3", "CD7", "SELL", "CD27", "IL7R")
  genes_of_interest <- intersect(genes_of_interest,
                                 rownames(target_Seurat_Obj@assays$RNA@counts))
  
  ### get GMP CARpos CD8 Px11 only data
  px11_gmp_seurat_obj <- subset(target_Seurat_Obj, cells = rownames(target_Seurat_Obj@meta.data)[which(target_Seurat_Obj@meta.data$px == "SJCAR19-11")])
  
  ### de analysis
  px11_gmp_seurat_obj <- SetIdent(object = px11_gmp_seurat_obj,
                                  cells = rownames(px11_gmp_seurat_obj@meta.data),
                                  value = px11_gmp_seurat_obj@meta.data$GMP_CARpos_CD8_Persister)
  px_11_de_result <- FindMarkers(px11_gmp_seurat_obj,
                                 ident.1 = "YES",
                                 ident.2 = "NO",
                                 min.pct = 0.2,
                                 logfc.threshold = 0.2,
                                 test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(px_11_de_result),
                         px_11_de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CARpos_CD8_Px11_Subsisters_vs_Non-Subsisters.xlsx"),
              sheetName = "GMP_CARpos_CD8_Px11_Subsisters_vs_Non-Subsisters", row.names = FALSE)
  
  ### DE result - dot plot
  px11_gmp_seurat_obj@meta.data$GMP_CARpos_CD8_Persister <- factor(px11_gmp_seurat_obj@meta.data$GMP_CARpos_CD8_Persister,
                                                                   levels = c("YES", "NO"))
  p <- DotPlot(px11_gmp_seurat_obj,
               features = genes_of_interest,
               cols = c("blue", "red"),
               group.by = "GMP_CARpos_CD8_Persister") +
    scale_size(range = c(2, 15)) +
    coord_flip() +
    xlab("Interesting Genes") +
    ylab("Is Subsistent") +
    theme_calc(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir2, "Dotplot_GMP_CARpos_CD8_Px11_Subsisters_vs_Non-subsisters.png"),
         plot = p, width = 15, height = 15, dpi = 350)
  
  ### DE result - ridge plot
  p <- RidgePlot(px11_gmp_seurat_obj, features = genes_of_interest, ncol = 5,
                 cols = c("#D21414", "#039076"))
  for(i in 1:length(genes_of_interest)) {
    p[[i]] <- p[[i]] + labs(y = "Is Subsistent")
  }
  ggsave(paste0(outputDir2, "Ridgeplot_GMP_CARpos_CD8_Px11_Subsisters_vs_Non-subsisters.png"), plot = p, width = 15, height = 7, dpi = 350)
  
  
  #
  ### 21. Multiple regression to estimate peakcar & wk1car using both dose level & tumor burden
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/21/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### set values
  # peakcar <- c(4215, 6178, 28867, 33667, 16709)
  peakcar <- c(199054, 42149, 61777, 288670, 224445, 167092, 4806, 142422,
               301705, 356424, 64212, 5835)
  wk1car_ug <- c(3114, 13912, 61777, 3745, 224445, 1011, 657, 46046,
                 268449, 356424, 8071, 5835)
  wk1car_ml <- c(9341, 296787, 947250, 24967, 1122227, 50700, 7561, 314646,
                 268449, 1960334, 236747, 66126)
  names(peakcar) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                      "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                      "SJCAR19-11", "SJCAR19-12")
  names(wk1car_ug) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                        "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                        "SJCAR19-11", "SJCAR19-12")
  names(wk1car_ml) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                        "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                        "SJCAR19-11", "SJCAR19-12")
  dose_level_table <- read.xlsx2(file = paste0(outputDir, "/8/Estimated_CD4_CD8_Subsister_#_In_GMP.xlsx"),
                                 sheetIndex = 1, row.names = 1,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  cd8_subsister_dose_level <- as.numeric(dose_level_table[names(peakcar),"Estimated Injected GMP CD8 Subsisters # (per kg)"])
  dose_level_table <- read.xlsx2(file = paste0(outputDir, "/13/Estimated_CD4_CD8_Subsister(Cluster22)_#_In_GMP.xlsx"),
                                 sheetIndex = 1, row.names = 1,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  cd8_cluster22_dose_level <- as.numeric(dose_level_table[names(peakcar),"Estimated Injected GMP CD8 Subsisters # (per kg)"])
  tumor_burden <- c(98, 10, 1, 0, 12, 1, 80, 72, 78, 84, 0, 2)
  names(tumor_burden) <- c("SJCAR19-01", "SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                           "SJCAR19-06", "SJCAR19-07", "SJCAR19-08", "SJCAR19-09", "SJCAR19-10",
                           "SJCAR19-11", "SJCAR19-12")
  
  ### make a data frame
  plot_df <- data.frame(CD8_Subsister_Dose_Level=cd8_subsister_dose_level,
                        CD8_Cluster22_Dose_Level=cd8_cluster22_dose_level,
                        PeakCAR=peakcar,
                        Wk1CAR_ug=wk1car_ug,
                        Wk1CAR_ml=wk1car_ml,
                        tumor_burden,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### remove NA rows
  plot_df <- plot_df[which(!is.na(plot_df$CD8_Subsister_Dose_Level)),]
  
  ### correlation plot
  png(filename = paste0(outputDir2, "Correlation_Between_Factors.png"), width = 2500, height = 1500, res = 350)
  pairs(data=plot_df,
        ~PeakCAR + tumor_burden + CD8_Subsister_Dose_Level + CD8_Cluster22_Dose_Level,
        pch = 19)
  dev.off()
  
  ### multiple regression - subsisters
  fit <- lm(PeakCAR ~ tumor_burden + CD8_Subsister_Dose_Level, data=plot_df)
  smr <- summary(fit)
  f <- smr$fstatistic
  pv <- pf(f[1],f[2],f[3],lower.tail=F)
  
  ### make the plot data frame
  new_plot_df <- data.frame(plot_df,
                            Residuals=fit$residuals,
                            Fitted_Values=fit$fitted.values,
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw the correlation plot
  p <- list()
  p[[1]] <- ggplot(data = new_plot_df, aes_string(x="PeakCAR", y="Fitted_Values")) +
    geom_point(col = "black", size = 8) +
    geom_abline(intercept = 0, slope = 1, col = "red", size = 2) +
    labs(title = paste0("Spearman Correlation = ", round(cor(new_plot_df$PeakCAR, new_plot_df$Fitted_Values, method = "spearman"), 2))) +
    xlab("PeakCAR") +
    ylab("Predicted Values") +
    geom_smooth(method = lm, color="blue", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24, color = "blue"))
  p[[2]] <- ggplot(data = new_plot_df, aes_string(x="Fitted_Values", y="Residuals")) +
    geom_point(col = "black", size = 8) +
    geom_line(size = 3) +
    labs(title = paste0("R2 = ", round(smr$r.squared, 2),
                        ", Adjusted R2 = ", round(smr$adj.r.squared, 2),
                        ", P-value = ", round(pv, 2))) +
    xlab("Predicted Values") +
    ylab("Residuals") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24))
  g <- arrangeGrob(grobs = p,
                   nrow = 1,
                   ncol = 2)
  ggsave(file = paste0(outputDir2, "PeakCAR_Multiple_Regression_Subsisters.png"), g, width = 25, height = 10, dpi = 400)
  
  ### multiple regression - cluster22
  fit <- lm(PeakCAR ~ tumor_burden + CD8_Cluster22_Dose_Level, data=plot_df)
  smr <- summary(fit)
  f <- smr$fstatistic
  pv <- pf(f[1],f[2],f[3],lower.tail=F)
  
  ### make the plot data frame
  new_plot_df <- data.frame(plot_df,
                            Residuals=fit$residuals,
                            Fitted_Values=fit$fitted.values,
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw the correlation plot
  p <- list()
  p[[1]] <- ggplot(data = new_plot_df, aes_string(x="PeakCAR", y="Fitted_Values")) +
    geom_point(col = "black", size = 8) +
    geom_abline(intercept = 0, slope = 1, col = "red", size = 2) +
    labs(title = paste0("Spearman Correlation = ", round(cor(new_plot_df$PeakCAR, new_plot_df$Fitted_Values, method = "spearman"), 2))) +
    xlab("PeakCAR") +
    ylab("Predicted Values") +
    geom_smooth(method = lm, color="blue", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24, color = "blue"))
  p[[2]] <- ggplot(data = new_plot_df, aes_string(x="Fitted_Values", y="Residuals")) +
    geom_point(col = "black", size = 8) +
    geom_line(size = 3) +
    labs(title = paste0("R2 = ", round(smr$r.squared, 2),
                        ", Adjusted R2 = ", round(smr$adj.r.squared, 2),
                        ", P-value = ", round(pv, 2))) +
    xlab("Predicted Values") +
    ylab("Residuals") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24))
  g <- arrangeGrob(grobs = p,
                   nrow = 1,
                   ncol = 2)
  ggsave(file = paste0(outputDir2, "PeakCAR_Multiple_Regression_Cluster22.png"), g, width = 25, height = 10, dpi = 400)
  
  
  ### r-squared function
  rsq <- function (x, y) round(cor(x, y)^2, 4)
  adj_rsq <- function(x, y, p, n) round(1-((1-cor(x, y)^2) * (n-1) / (n-p-1)), 4)
  
  ### f-statistic function
  fstat <- function(x, y, df1, df2) round((cor(x, y)^2 / (1-cor(x, y)^2)) * (df2 / df1), 4)
  fstat_pv <- function(x, y, df1, df2) round(pf((cor(x, y)^2 / (1-cor(x, y)^2)) * (df2 / df1), df1, df2, lower.tail = FALSE), 4)
  
  
  ### SVM regression - subsisters
  svm_model = svm(PeakCAR ~ tumor_burden + CD8_Subsister_Dose_Level, data=plot_df)
  predictYsvm <- predict(svm_model, data=plot_df)
  
  ### make the plot data frame
  new_plot_df <- data.frame(plot_df,
                            Residuals=predictYsvm-plot_df$PeakCAR,
                            Fitted_Values=predictYsvm,
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw the correlation plot
  p <- list()
  p[[1]] <- ggplot(data = new_plot_df, aes_string(x="PeakCAR", y="Fitted_Values")) +
    geom_point(col = "black", size = 8) +
    geom_abline(intercept = 0, slope = 1, col = "red", size = 2) +
    labs(title = paste0("Spearman Correlation = ", round(cor(new_plot_df$PeakCAR, new_plot_df$Fitted_Values, method = "spearman"), 2))) +
    xlab("PeakCAR") +
    ylab("Predicted Values") +
    geom_smooth(method = lm, color="blue", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24, color = "blue"))
  p[[2]] <- ggplot(data = new_plot_df, aes_string(x="Fitted_Values", y="Residuals")) +
    geom_point(col = "black", size = 8) +
    geom_line(size = 3) +
    labs(title = paste0("R2 = ", rsq(new_plot_df$PeakCAR, new_plot_df$Fitted_Values),
                        ", Adjusted R2 = ", adj_rsq(new_plot_df$PeakCAR, new_plot_df$Fitted_Values, 2, length(plot_df$PeakCAR)),
                        ", P-value = ", fstat_pv(new_plot_df$PeakCAR, new_plot_df$Fitted_Values, 2, 9))) +
    xlab("Predicted Values") +
    ylab("Residuals") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24))
  g <- arrangeGrob(grobs = p,
                   nrow = 1,
                   ncol = 2)
  ggsave(file = paste0(outputDir2, "PeakCAR_SVM_Regression_Subsisters.png"), g, width = 25, height = 10, dpi = 400)
  
  ### SVM regression - cluster22
  svm_model = svm(PeakCAR ~ tumor_burden + CD8_Cluster22_Dose_Level, data=plot_df)
  predictYsvm <- predict(svm_model, data=plot_df)
  
  ### make the plot data frame
  new_plot_df <- data.frame(plot_df,
                            Residuals=predictYsvm-plot_df$PeakCAR,
                            Fitted_Values=predictYsvm,
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw the correlation plot
  p <- list()
  p[[1]] <- ggplot(data = new_plot_df, aes_string(x="PeakCAR", y="Fitted_Values")) +
    geom_point(col = "black", size = 8) +
    geom_abline(intercept = 0, slope = 1, col = "red", size = 2) +
    labs(title = paste0("Spearman Correlation = ", round(cor(new_plot_df$PeakCAR, new_plot_df$Fitted_Values, method = "spearman"), 2))) +
    xlab("PeakCAR") +
    ylab("Predicted Values") +
    geom_smooth(method = lm, color="blue", se=TRUE) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24, color = "blue"))
  p[[2]] <- ggplot(data = new_plot_df, aes_string(x="Fitted_Values", y="Residuals")) +
    geom_point(col = "black", size = 8) +
    geom_line(size = 3) +
    labs(title = paste0("R2 = ", rsq(new_plot_df$PeakCAR, new_plot_df$Fitted_Values),
                        ", Adjusted R2 = ", adj_rsq(new_plot_df$PeakCAR, new_plot_df$Fitted_Values, 2, length(plot_df$PeakCAR)),
                        ", P-value = ", fstat_pv(new_plot_df$PeakCAR, new_plot_df$Fitted_Values, 2, 9))) +
    xlab("Predicted Values") +
    ylab("Residuals") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0, vjust = 0.5, size = 24))
  g <- arrangeGrob(grobs = p,
                   nrow = 1,
                   ncol = 2)
  ggsave(file = paste0(outputDir2, "PeakCAR_SVM_Regression_Cluster22.png"), g, width = 25, height = 10, dpi = 400)
  
  
  #
  ### 22. Classifier - reperform with the new data and also with the subsister cluster info (not only based on the subsisters all the cluster cells)
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/22/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  #'******************************************************************************
  #' A function to transform RNA-Seq data with VST in DESeq2 package
  #' readCount: RNA-Seq rawcounts in a matrix or in a data frame form
  #'            Rows are genes and columns are samples
  #' filter_thresh: The function filters out genes that have at least one sample
  #'                with counts larger than the 'filter_thresh' value
  #'                e.g., if the 'filter_thresh' = 1, then it removes genes
  #'                that have counts <= 1 across all the samples
  #'                if 0, then there will be no filtering
  #'******************************************************************************
  normalizeRNASEQwithVST <- function(readCount, filter_thresh=1) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filter_thresh > 0) {
      ### Remove rubbish rows - this will decrease the number of rows
      keep = apply(counts(deSeqData), 1, function(r){
        return(sum(r > filter_thresh) > 0)
      })
      deSeqData <- deSeqData[keep,]
    }
    
    ### VST
    vsd <- varianceStabilizingTransformation(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  ### only get the GMP persisters and non-persisters
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[union(which(Seurat_Obj$GMP_CARpos_CD8_Persister == "YES"),
                                                                                       which(Seurat_Obj$GMP_CARpos_CD8_Persister == "NO"))])
  
  ### only the patients we are interested
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  target_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                            "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                            "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### normalization
  target_Seurat_Obj <- NormalizeData(target_Seurat_Obj,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj <- FindVariableFeatures(target_Seurat_Obj,
                                            selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj <- ScaleData(target_Seurat_Obj,
                                 vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj <- RunPCA(target_Seurat_Obj,
                              features = VariableFeatures(object = target_Seurat_Obj),
                              npcs = 15)
  target_Seurat_Obj <- RunUMAP(target_Seurat_Obj, dims = 1:15)
  
  ### perform clustering
  target_Seurat_Obj <- FindNeighbors(target_Seurat_Obj, dims = 1:10)
  target_Seurat_Obj <- FindClusters(target_Seurat_Obj, resolution = 1.5)
  
  ### save the clustering result to meta.data
  target_Seurat_Obj@meta.data$gmp_seurat_clusters <- Idents(target_Seurat_Obj)
  
  ### specific cluster - cluster 22
  DimPlot(object = target_Seurat_Obj, reduction = "umap",
          group.by = "gmp_seurat_clusters",
          pt.size = 1)
  target_Seurat_Obj@meta.data$gmp_seurat_clusters2 <- as.character(target_Seurat_Obj@meta.data$gmp_seurat_clusters)
  target_Seurat_Obj@meta.data$gmp_seurat_clusters2[which(target_Seurat_Obj@meta.data$gmp_seurat_clusters2 != "22")] <- "Others"
  DimPlot(object = target_Seurat_Obj, reduction = "umap",
          group.by = "gmp_seurat_clusters2",
          pt.size = 1,
          cols = c("Others" = "lightgray", "22" = "red"),
          order = c("22", "Others"))
  
  ### because of the imbalance of the two cluster sizes, we randomly choose
  ### the same number of samples in each class and iteratively build the classifier
  iteration <- 10
  set.seed(2990)
  featureSelectionNum <- 100
  sampleNum <- 100
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "parRF", "glmboost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "Linear_Model", "K-NN")
  train_control <- trainControl(method="LOOCV", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
  
  ### a small number that will be added to the counts for normalization
  ### log transform should not have 0 values
  log_trans_add <- 1
  
  ### de result - subsisters
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister)
  gmp_de_result <- FindMarkers(target_Seurat_Obj,
                               ident.1 = "YES",
                               ident.2 = "NO",
                               min.pct = 0.2,
                               logfc.threshold = 0.2,
                               test.use = "wilcox")
  
  ### build the classifier 10 times - Subsisters
  for(i in 1:iteration) {
    
    ### get random samples for each condition
    cond1_samps <- rownames(target_Seurat_Obj@meta.data)[sample(which(target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES"), sampleNum)]
    cond2_samps <- rownames(target_Seurat_Obj@meta.data)[sample(which(target_Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "NO"), sampleNum)]
    
    ### normalize the read counts
    ### before the normalization, only keep the samples that will be used in the classifier
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(target_Seurat_Obj@assays$RNA@counts[rownames(gmp_de_result)[1:featureSelectionNum],
                                                                                             c(cond1_samps,
                                                                                               cond2_samps)] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE),
                                         filter_thresh = 0)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(c(rep("GMP_Last", sampleNum),
                                 rep("GMP_Not_Last", sampleNum)),
                               levels = c("GMP_Last", "GMP_Not_Last"))
    
    ### build classifier and test
    ### LOOCV
    p <- list()
    acc <- NULL
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=input_data, trControl=train_control, method=methodTypes[j])
      roc <- roc(model$pred$obs, model$pred$GMP_Last)
      acc <- c(acc, round(mean(model$results$Accuracy), 3))
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using DE Genes\n",
                                           "Accuracy =", acc[j]),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      gc()
    }
    
    ### draw ROC curves
    png(paste0(outputDir2, "Classifier_GMP_CARpos_CD8_Subsisters_vs_Non-Subsisters_", featureSelectionNum, "_(", i, ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", acc[j]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
    
    gc()
  }
  
  ### de result - cluster22
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$gmp_seurat_clusters2)
  gmp_de_result_c22 <- FindMarkers(target_Seurat_Obj,
                                   ident.1 = "22",
                                   ident.2 = "Others",
                                   min.pct = 0.2,
                                   logfc.threshold = 0.2,
                                   test.use = "wilcox")
  
  ### build the classifier 10 times - cluster22
  set.seed(2990)
  for(i in 1:iteration) {
    
    ### get random samples for each condition
    cond1_samps <- rownames(target_Seurat_Obj@meta.data)[sample(which(target_Seurat_Obj@meta.data$gmp_seurat_clusters2 == "22"), sampleNum)]
    cond2_samps <- rownames(target_Seurat_Obj@meta.data)[sample(which(target_Seurat_Obj@meta.data$gmp_seurat_clusters2 == "Others"), sampleNum)]
    
    ### normalize the read counts
    ### before the normalization, only keep the samples that will be used in the classifier
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(target_Seurat_Obj@assays$RNA@counts[rownames(gmp_de_result_c22)[1:featureSelectionNum],
                                                                                                    c(cond1_samps,
                                                                                                      cond2_samps)] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE),
                                         filter_thresh = 0)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(c(rep("GMP_Last", sampleNum),
                                 rep("GMP_Not_Last", sampleNum)),
                               levels = c("GMP_Last", "GMP_Not_Last"))
    
    ### build classifier and test
    ### LOOCV
    p <- list()
    acc <- NULL
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=input_data, trControl=train_control, method=methodTypes[j])
      roc <- roc(model$pred$obs, model$pred$GMP_Last)
      acc <- c(acc, round(mean(model$results$Accuracy), 3))
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using DE Genes\n",
                                           "Accuracy =", acc[j]),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      gc()
    }
    
    ### draw ROC curves
    png(paste0(outputDir2, "Classifier_GMP_CARpos_CD8_C22_vs_Others_", featureSelectionNum, "_(", i, ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", acc[j]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
    
    gc()
  }
  
  #
  ### 23. PCA/UMAP/MNN comparison with the original GMP CARpos cells
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/23/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### only get the GMP CARpos cells (no GMP-redo)
  target_Seurat_Obj3 <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[union(which(Seurat_Obj$time == "GMP"),
                                                                                        which(Seurat_Obj$CAR == "CARpos"))])
  
  ### only the patients we are interested
  target_Seurat_Obj3 <- SetIdent(object = target_Seurat_Obj3,
                                 cells = rownames(target_Seurat_Obj3@meta.data),
                                 value = target_Seurat_Obj3@meta.data$px)
  target_Seurat_Obj3 <- subset(target_Seurat_Obj3, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                              "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                              "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### run mnn
  target_Seurat_Obj3.list <- SplitObject(target_Seurat_Obj3, split.by = "px")
  target_Seurat_Obj3 <- RunFastMNN(object.list = target_Seurat_Obj3.list)
  rm(target_Seurat_Obj3.list)
  gc()
  
  ### normalization
  target_Seurat_Obj3 <- NormalizeData(target_Seurat_Obj3,
                                      normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj3 <- FindVariableFeatures(target_Seurat_Obj3,
                                             selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj3 <- ScaleData(target_Seurat_Obj3,
                                  vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj3 <- RunPCA(target_Seurat_Obj3,
                               features = VariableFeatures(object = target_Seurat_Obj3),
                               npcs = 15)
  target_Seurat_Obj3 <- RunUMAP(target_Seurat_Obj3, dims = 1:15)
  
  ### draw a PCA plot
  p <- DimPlot(object = target_Seurat_Obj3, reduction = "pca",
               group.by = "px", pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Patient") +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "PCA_Original_GMP_Only_CARpos_Px.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### draw a UMAP plot
  p <- DimPlot(object = target_Seurat_Obj3, reduction = "umap",
               group.by = "px", pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Patient") +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "UMAP_Original_GMP_Only_CARpos_Px.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### draw a MNN plot
  p <- DimPlot(object = target_Seurat_Obj3, reduction = "mnn",
               group.by = "px", pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Patient") +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "MNN_Original_GMP_Only_CARpos_Px.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  #
  ### 24. Mapping GMP clusters to PI clusters
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/24/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get CARpos-only seurat object
  carpos_cd8_cells <- rownames(Seurat_Obj@meta.data)[intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                               which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))]
  sub_seurat_obj4 <- subset(Seurat_Obj, cells = carpos_cd8_cells)
  
  ### after gmp time points only
  after_gmp_time_points <- c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", 
                             "Wk6", "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$time2)
  sub_seurat_obj4 <- subset(sub_seurat_obj4, idents = intersect(after_gmp_time_points,
                                                                unique(sub_seurat_obj4@meta.data$time2)))
  
  ### get seurat object for some specific patients
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$px)
  sub_seurat_obj4 <- subset(sub_seurat_obj4, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                        "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                        "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### run mnn
  sub_seurat_obj4.list <- SplitObject(sub_seurat_obj4, split.by = "library")
  sub_seurat_obj4 <- RunFastMNN(object.list = sub_seurat_obj4.list)
  rm(sub_seurat_obj4.list)
  gc()
  
  ### normalization
  sub_seurat_obj4 <- NormalizeData(sub_seurat_obj4,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  sub_seurat_obj4 <- FindVariableFeatures(sub_seurat_obj4,
                                          selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  sub_seurat_obj4 <- ScaleData(sub_seurat_obj4,
                               vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap & clustering
  sub_seurat_obj4 <- RunPCA(sub_seurat_obj4,
                            features = VariableFeatures(object = sub_seurat_obj4),
                            npcs = 15)
  sub_seurat_obj4 <- RunUMAP(sub_seurat_obj4, reduction = "mnn", dims = 1:15)
  sub_seurat_obj4 <- FindNeighbors(sub_seurat_obj4, reduction = "mnn", dims = 1:15)
  sub_seurat_obj4 <- FindClusters(sub_seurat_obj4, resolution = 0.5)
  
  ### save the clustering result to meta.data
  sub_seurat_obj4@meta.data$clusters <- Idents(sub_seurat_obj4)
  
  ### UMAP with clusters by each time
  p <- DimPlot(object = sub_seurat_obj4, reduction = "umap",
               group.by = "clusters", split.by = "time2",
               pt.size = 3, ncol = 3) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "MNN_UMAP_CARpos_CD8_Clusters.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### add GMP - PI distinguishable column
  sub_seurat_obj4@meta.data$GMP_PI <- "PI"
  sub_seurat_obj4@meta.data$GMP_PI[which(sub_seurat_obj4@meta.data$time2 == "GMP")] <- "GMP"
  
  ### UMAP with mapping for each cluster (GMP - PI)
  p <- vector("list", length(levels(sub_seurat_obj4@meta.data$clusters)))
  names(p) <- levels(sub_seurat_obj4@meta.data$clusters)
  for(clstr in levels(sub_seurat_obj4@meta.data$clusters)) {
    ### get seurat object for the given cluster
    temp_seurat_obj <- subset(sub_seurat_obj4, cells = rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$clusters == clstr)])
    
    ### draw a UMAP within the given cluster - coloring with GMP & PI
    p[[clstr]] <- DimPlot(object = temp_seurat_obj, reduction = "umap",
                          group.by = "GMP_PI", cols = c("GMP" = "blue", "PI" = "green"),
                          pt.size = 3) +
      ggtitle(paste("Cluster", clstr)) +
      labs(color="") +
      theme_classic(base_size = 36) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    ### transpency
    p[[clstr]][[1]]$layers[[1]]$aes_params$alpha <- 0.3
  }
  
  g <- arrangeGrob(grobs = p,
                   nrow = 4,
                   ncol = 3)
  ggsave(file = paste0(outputDir2, "MNN_UMAP_GMP_PI_within_Clusters.png"), g, width = 25, height = 15, dpi = 350)
  
  ### how many subsisters have lineageS between GMP and PI in each cluster?
  ### but this can't say whether two cells are in one lineage
  sub_seurat_obj4@meta.data$GMP_PI2 <- sub_seurat_obj4@meta.data$GMP_PI
  sub_seurat_obj4@meta.data$GMP_PI2[intersect(which(sub_seurat_obj4@meta.data$GMP_PI == "GMP"),
                                              which(sub_seurat_obj4@meta.data$ALL_GMP_CARpos_Persister == "YES"))] <- "GMP_Subsister"
  sub_seurat_obj4@meta.data$GMP_PI2[intersect(which(sub_seurat_obj4@meta.data$GMP_PI == "PI"),
                                              which(sub_seurat_obj4@meta.data$ALL_GMP_CARpos_Persister == "YES"))] <- "PI_Subsister"
  
  for(clstr in levels(sub_seurat_obj4@meta.data$clusters)) {
    ### get seurat object for the given cluster
    temp_seurat_obj <- subset(sub_seurat_obj4, cells = rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$clusters == clstr)])
    
    ### get subsister lineages in the given cluster
    temp_subsister_clones <- unique(temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(temp_seurat_obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
    temp_subsister_clones <- temp_subsister_clones[which(!is.na(temp_subsister_clones))]
    temp_subsister_clones <- intersect(temp_subsister_clones, temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(temp_seurat_obj@meta.data$GMP_PI == "PI")])
    
    ### define subsisters in the same lineage in the given cluster
    temp_seurat_obj@meta.data$GMP_PI3 <- temp_seurat_obj@meta.data$GMP_PI
    temp_seurat_obj@meta.data$GMP_PI3[intersect(which(temp_seurat_obj@meta.data$GMP_PI == "GMP"),
                                                which(temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% temp_subsister_clones))] <- "GMP_Subsister"
    temp_seurat_obj@meta.data$GMP_PI3[intersect(which(temp_seurat_obj@meta.data$GMP_PI == "PI"),
                                                which(temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% temp_subsister_clones))] <- "PI_Subsister"
    
    
    ### draw a UMAP within the given cluster - coloring with GMP & PI
    p[[clstr]] <- DimPlot(object = temp_seurat_obj, reduction = "umap",
                          group.by = "GMP_PI3", cols = c("GMP" = "blue", "PI" = "green", "GMP_Subsister" = "red", "PI_Subsister" = "orange"),
                          order = c("PI_Subsister", "GMP_Subsister", "PI", "GMP"),
                          pt.size = 3) +
      ggtitle(paste("Cluster", clstr, "-",
                    length(temp_subsister_clones), "Out of",
                    length(intersect(unique(sub_seurat_obj4@meta.data$clonotype_id_by_patient_one_alpha_beta[which(sub_seurat_obj4@meta.data$GMP_PI2 == "PI_Subsister")]),
                                     unique(temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(temp_seurat_obj@meta.data$GMP_CARpos_CD8_Persister == "YES")]))))) +
      labs(color="") +
      theme_classic(base_size = 36) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    ### transpency
    p[[clstr]][[1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  g <- arrangeGrob(grobs = p,
                   nrow = 4,
                   ncol = 3)
  ggsave(file = paste0(outputDir2, "MNN_UMAP_GMP_PI_within_Clusters_Subsisters.png"), g, width = 30, height = 15, dpi = 350)
  
  #
  ### look at GMP subsisters in each cluster and find in which clusters they are in the PI
  ### proportional bar graph
  #
  
  ### make a data frame for the plot
  plot_df <- data.frame(GMP_Subsister_Cluster=as.character(sapply(levels(sub_seurat_obj4@meta.data$clusters),
                                                                  function(x) rep(x, length(levels(sub_seurat_obj4@meta.data$clusters))))),
                        PI_Subsister_Cluster=rep(levels(sub_seurat_obj4@meta.data$clusters), length(levels(sub_seurat_obj4@meta.data$clusters))),
                        Lineage_Num = 0,
                        Lineage_Pct = 0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### fill out the plot_df
  i <- 1
  for(clstr1 in levels(sub_seurat_obj4@meta.data$clusters)) {
    ### clstr1 indicies
    clstr1_index <- which(sub_seurat_obj4@meta.data$clusters == clstr1)
    clstr1_gmp_subsister_index <- intersect(clstr1_index,
                                            which(sub_seurat_obj4@meta.data$GMP_CARpos_CD8_Persister == "YES"))
    
    ### get clstr1 gmp subsister clones
    clstr1_gmp_subsister_clones <- unique(sub_seurat_obj4@meta.data$clonotype_id_by_patient_one_alpha_beta[clstr1_gmp_subsister_index])
    clstr1_gmp_subsister_clones <- clstr1_gmp_subsister_clones[which(!is.na(clstr1_gmp_subsister_clones))]
    
    ### for clstr2 find GMP - PI lineage in the given clusters
    for(clstr2 in levels(sub_seurat_obj4@meta.data$clusters)) {
      ### clstr2 indicies
      clstr2_index <- which(sub_seurat_obj4@meta.data$clusters == clstr2)
      clstr2_pi_index <- intersect(clstr2_index,
                                   which(sub_seurat_obj4@meta.data$GMP_PI == "PI"))
      clstr1_clstr2_pi_subsister_index <- intersect(clstr2_pi_index,
                                                    which(sub_seurat_obj4@meta.data$clonotype_id_by_patient_one_alpha_beta %in% clstr1_gmp_subsister_clones))
      
      ### get the GMP - PI lineage numbers
      plot_df$Lineage_Num[i] <- length(unique(sub_seurat_obj4@meta.data$clonotype_id_by_patient_one_alpha_beta[clstr1_clstr2_pi_subsister_index]))
      
      ### iteration count ++
      i <- i+1
    }
    
    ### calculate the percentage
    lineage_sum <- sum(plot_df$Lineage_Num[which(plot_df$GMP_Subsister_Cluster == clstr1)])
    plot_df$Lineage_Pct[which(plot_df$GMP_Subsister_Cluster == clstr1)] <- round(plot_df$Lineage_Num[which(plot_df$GMP_Subsister_Cluster == clstr1)] * 100 / lineage_sum, 1)
  }
  
  ### draw the proportional bar graph
  
  ### remove 0 rows
  plot_df <- plot_df[which(plot_df$Lineage_Num != 0),]
  
  ### factorize the columns
  plot_df$GMP_Subsister_Cluster <- factor(plot_df$GMP_Subsister_Cluster, levels = levels(sub_seurat_obj4@meta.data$clusters))
  plot_df$PI_Subsister_Cluster <- factor(plot_df$PI_Subsister_Cluster, levels = levels(sub_seurat_obj4@meta.data$clusters))
  
  ### draw a proportional bar plot
  ### pcnt < 30 -> ""
  plot_df$Lineage_Pct[which(as.numeric(plot_df$Lineage_Pct) < 30)] <- ""
  plot_df$Lineage_Pct <- as.character(plot_df$Lineage_Pct)
  p <- ggplot(data=plot_df, aes_string(x="GMP_Subsister_Cluster", y="Lineage_Num", fill="PI_Subsister_Cluster", label="Lineage_Pct")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Subsister Lineages Across Clusters") +
    geom_text(size = 5, position = position_stack(vjust = 1)) +
    xlab("GMP_Subsister_Cluster") +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "GMP_CARpos_CD8_Lineages_Across_Clusters.png"), plot = p,
         width = 25, height = 12, dpi = 350)
  
  
  ### try with different clustering
  sub_seurat_obj4 <- FindClusters(sub_seurat_obj4, resolution = 0.1)
  sub_seurat_obj4@meta.data$clusters2 <- Idents(sub_seurat_obj4)
  
  ### UMAP with clusters by each patient
  p <- DimPlot(object = sub_seurat_obj4, reduction = "umap",
               group.by = "clusters2", split.by = "time2",
               pt.size = 3, ncol = 3) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "MNN_UMAP_CARpos_CD8_Clusters2.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### UMAP with mapping for each cluster (GMP - PI)
  p <- vector("list", length(levels(sub_seurat_obj4@meta.data$clusters2)))
  names(p) <- levels(sub_seurat_obj4@meta.data$clusters2)
  for(clstr in levels(sub_seurat_obj4@meta.data$clusters2)) {
    ### get seurat object for the given cluster
    temp_seurat_obj <- subset(sub_seurat_obj4, cells = rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$clusters2 == clstr)])
    
    ### draw a UMAP within the given cluster - coloring with GMP & PI
    p[[clstr]] <- DimPlot(object = temp_seurat_obj, reduction = "umap",
                          group.by = "GMP_PI", cols = c("GMP" = "blue", "PI" = "green"),
                          pt.size = 3) +
      ggtitle(paste("Cluster2", clstr)) +
      labs(color="") +
      theme_classic(base_size = 36) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    ### transpency
    p[[clstr]][[1]]$layers[[1]]$aes_params$alpha <- 0.3
  }
  
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   ncol = 2)
  ggsave(file = paste0(outputDir2, "MNN_UMAP_GMP_PI_within_Clusters2.png"), g, width = 25, height = 15, dpi = 350)
  
  for(clstr in levels(sub_seurat_obj4@meta.data$clusters2)) {
    ### get seurat object for the given cluster
    temp_seurat_obj <- subset(sub_seurat_obj4, cells = rownames(sub_seurat_obj4@meta.data)[which(sub_seurat_obj4@meta.data$clusters2 == clstr)])
    
    ### get subsister lineages in the given cluster
    temp_subsister_clones <- unique(temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(temp_seurat_obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
    temp_subsister_clones <- temp_subsister_clones[which(!is.na(temp_subsister_clones))]
    temp_subsister_clones <- intersect(temp_subsister_clones, temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(temp_seurat_obj@meta.data$GMP_PI == "PI")])
    
    ### define subsisters in the same lineage in the given cluster
    temp_seurat_obj@meta.data$GMP_PI3 <- temp_seurat_obj@meta.data$GMP_PI
    temp_seurat_obj@meta.data$GMP_PI3[intersect(which(temp_seurat_obj@meta.data$GMP_PI == "GMP"),
                                                which(temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% temp_subsister_clones))] <- "GMP_Subsister"
    temp_seurat_obj@meta.data$GMP_PI3[intersect(which(temp_seurat_obj@meta.data$GMP_PI == "PI"),
                                                which(temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% temp_subsister_clones))] <- "PI_Subsister"
    
    
    ### draw a UMAP within the given cluster - coloring with GMP & PI
    p[[clstr]] <- DimPlot(object = temp_seurat_obj, reduction = "umap",
                          group.by = "GMP_PI3", cols = c("GMP" = "blue", "PI" = "green", "GMP_Subsister" = "red", "PI_Subsister" = "orange"),
                          order = c("PI_Subsister", "GMP_Subsister", "PI", "GMP"),
                          pt.size = 3) +
      ggtitle(paste("Cluster", clstr, "-",
                    length(temp_subsister_clones), "Out of",
                    length(intersect(unique(sub_seurat_obj4@meta.data$clonotype_id_by_patient_one_alpha_beta[which(sub_seurat_obj4@meta.data$GMP_PI2 == "PI_Subsister")]),
                                     unique(temp_seurat_obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(temp_seurat_obj@meta.data$GMP_CARpos_CD8_Persister == "YES")]))))) +
      labs(color="") +
      theme_classic(base_size = 36) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    ### transpency
    p[[clstr]][[1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   ncol = 2)
  ggsave(file = paste0(outputDir2, "MNN_UMAP_GMP_PI_within_Clusters2_Subsisters.png"), g, width = 30, height = 15, dpi = 350)
  
  #
  ### look at GMP subsisters in each cluster and find in which clusters they are in the PI
  ### proportional bar graph
  #
  
  ### make a data frame for the plot
  plot_df <- data.frame(GMP_Subsister_Cluster=as.character(sapply(levels(sub_seurat_obj4@meta.data$clusters2),
                                                                  function(x) rep(x, length(levels(sub_seurat_obj4@meta.data$clusters2))))),
                        PI_Subsister_Cluster=rep(levels(sub_seurat_obj4@meta.data$clusters2), length(levels(sub_seurat_obj4@meta.data$clusters2))),
                        Lineage_Num = 0,
                        Lineage_Pct = 0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### fill out the plot_df
  i <- 1
  for(clstr1 in levels(sub_seurat_obj4@meta.data$clusters2)) {
    ### clstr1 indicies
    clstr1_index <- which(sub_seurat_obj4@meta.data$clusters2 == clstr1)
    clstr1_gmp_subsister_index <- intersect(clstr1_index,
                                            which(sub_seurat_obj4@meta.data$GMP_CARpos_CD8_Persister == "YES"))
    
    ### get clstr1 gmp subsister clones
    clstr1_gmp_subsister_clones <- unique(sub_seurat_obj4@meta.data$clonotype_id_by_patient_one_alpha_beta[clstr1_gmp_subsister_index])
    clstr1_gmp_subsister_clones <- clstr1_gmp_subsister_clones[which(!is.na(clstr1_gmp_subsister_clones))]
    
    ### for clstr2 find GMP - PI lineage in the given clusters
    for(clstr2 in levels(sub_seurat_obj4@meta.data$clusters2)) {
      ### clstr2 indicies
      clstr2_index <- which(sub_seurat_obj4@meta.data$clusters2 == clstr2)
      clstr2_pi_index <- intersect(clstr2_index,
                                   which(sub_seurat_obj4@meta.data$GMP_PI == "PI"))
      clstr1_clstr2_pi_subsister_index <- intersect(clstr2_pi_index,
                                                    which(sub_seurat_obj4@meta.data$clonotype_id_by_patient_one_alpha_beta %in% clstr1_gmp_subsister_clones))
      
      ### get the GMP - PI lineage numbers
      plot_df$Lineage_Num[i] <- length(unique(sub_seurat_obj4@meta.data$clonotype_id_by_patient_one_alpha_beta[clstr1_clstr2_pi_subsister_index]))
      
      ### iteration count ++
      i <- i+1
    }
    
    ### calculate the percentage
    lineage_sum <- sum(plot_df$Lineage_Num[which(plot_df$GMP_Subsister_Cluster == clstr1)])
    plot_df$Lineage_Pct[which(plot_df$GMP_Subsister_Cluster == clstr1)] <- round(plot_df$Lineage_Num[which(plot_df$GMP_Subsister_Cluster == clstr1)] * 100 / lineage_sum, 1)
  }
  
  ### draw the proportional bar graph
  
  ### remove 0 rows
  plot_df <- plot_df[which(plot_df$Lineage_Num != 0),]
  
  ### factorize the columns
  plot_df$GMP_Subsister_Cluster <- factor(plot_df$GMP_Subsister_Cluster, levels = levels(sub_seurat_obj4@meta.data$clusters2))
  plot_df$PI_Subsister_Cluster <- factor(plot_df$PI_Subsister_Cluster, levels = levels(sub_seurat_obj4@meta.data$clusters2))
  
  ### draw a proportional bar plot
  ### pcnt < 20 -> ""
  plot_df$Lineage_Pct[which(as.numeric(plot_df$Lineage_Pct) < 20)] <- ""
  plot_df$Lineage_Pct <- as.character(plot_df$Lineage_Pct)
  p <- ggplot(data=plot_df, aes_string(x="GMP_Subsister_Cluster", y="Lineage_Num", fill="PI_Subsister_Cluster", label="Lineage_Pct")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Subsister Lineages Across Clusters") +
    geom_text(size = 5, position = position_stack(vjust = 1)) +
    xlab("GMP_Subsister_Cluster") +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "GMP_CARpos_CD8_Lineages_Across_Clusters2.png"), plot = p,
         width = 25, height = 12, dpi = 350)
  
  #
  ### separate the GMP & PI
  #
  
  ### separate the seurat object
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$GMP_PI)
  gmp_carpos_cd8_seurat_obj <- subset(sub_seurat_obj4, idents = c("GMP"))
  pi_carpos_cd8_seurat_obj <- subset(sub_seurat_obj4, idents = c("PI"))
  
  ### preprocess each seurat object
  
  ### run mnn
  gmp_carpos_cd8_seurat_obj.list <- SplitObject(gmp_carpos_cd8_seurat_obj, split.by = "library")
  gmp_carpos_cd8_seurat_obj <- RunFastMNN(object.list = gmp_carpos_cd8_seurat_obj.list)
  rm(gmp_carpos_cd8_seurat_obj.list)
  gc()
  pi_carpos_cd8_seurat_obj.list <- SplitObject(pi_carpos_cd8_seurat_obj, split.by = "library")
  pi_carpos_cd8_seurat_obj <- RunFastMNN(object.list = pi_carpos_cd8_seurat_obj.list)
  rm(pi_carpos_cd8_seurat_obj.list)
  gc()
  
  ### normalization
  gmp_carpos_cd8_seurat_obj <- NormalizeData(gmp_carpos_cd8_seurat_obj,
                                             normalization.method = "LogNormalize", scale.factor = 10000)
  pi_carpos_cd8_seurat_obj <- NormalizeData(pi_carpos_cd8_seurat_obj,
                                            normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  gmp_carpos_cd8_seurat_obj <- FindVariableFeatures(gmp_carpos_cd8_seurat_obj,
                                                    selection.method = "vst", nfeatures = 2000)
  pi_carpos_cd8_seurat_obj <- FindVariableFeatures(pi_carpos_cd8_seurat_obj,
                                                   selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  gmp_carpos_cd8_seurat_obj <- ScaleData(gmp_carpos_cd8_seurat_obj,
                                         vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  pi_carpos_cd8_seurat_obj <- ScaleData(pi_carpos_cd8_seurat_obj,
                                        vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap & clustering
  gmp_carpos_cd8_seurat_obj <- RunPCA(gmp_carpos_cd8_seurat_obj,
                                      features = VariableFeatures(object = gmp_carpos_cd8_seurat_obj),
                                      npcs = 15)
  gmp_carpos_cd8_seurat_obj <- RunUMAP(gmp_carpos_cd8_seurat_obj, reduction = "mnn", dims = 1:15)
  gmp_carpos_cd8_seurat_obj <- FindNeighbors(gmp_carpos_cd8_seurat_obj, reduction = "mnn", dims = 1:15)
  pi_carpos_cd8_seurat_obj <- RunPCA(pi_carpos_cd8_seurat_obj,
                                     features = VariableFeatures(object = pi_carpos_cd8_seurat_obj),
                                     npcs = 15)
  pi_carpos_cd8_seurat_obj <- RunUMAP(pi_carpos_cd8_seurat_obj, reduction = "mnn", dims = 1:15)
  pi_carpos_cd8_seurat_obj <- FindNeighbors(pi_carpos_cd8_seurat_obj, reduction = "mnn", dims = 1:15)
  
  gmp_carpos_cd8_seurat_obj <- FindClusters(gmp_carpos_cd8_seurat_obj, resolution = 0.3)
  pi_carpos_cd8_seurat_obj <- FindClusters(pi_carpos_cd8_seurat_obj, resolution = 0.3)
  
  ### save the clustering result to meta.data
  gmp_carpos_cd8_seurat_obj@meta.data$clusters <- Idents(gmp_carpos_cd8_seurat_obj)
  pi_carpos_cd8_seurat_obj@meta.data$clusters <- Idents(pi_carpos_cd8_seurat_obj)
  
  ### MNN UMAP with clusters
  p <- DimPlot(object = gmp_carpos_cd8_seurat_obj, reduction = "umap",
               group.by = "clusters",
               pt.size = 3) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "MNN_UMAP_GMP_CARpos_CD8_Clusters.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- DimPlot(object = pi_carpos_cd8_seurat_obj, reduction = "umap",
               group.by = "clusters",
               pt.size = 3) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "MNN_UMAP_PI_CARpos_CD8_Clusters.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### find all markers
  
  ### set cluster info as idents
  gmp_carpos_cd8_seurat_obj <- SetIdent(object = gmp_carpos_cd8_seurat_obj,
                                        cells = rownames(gmp_carpos_cd8_seurat_obj@meta.data),
                                        value = gmp_carpos_cd8_seurat_obj@meta.data$clusters)
  pi_carpos_cd8_seurat_obj <- SetIdent(object = pi_carpos_cd8_seurat_obj,
                                       cells = rownames(pi_carpos_cd8_seurat_obj@meta.data),
                                       value = pi_carpos_cd8_seurat_obj@meta.data$clusters)
  
  ### DE analysis
  de_result <- FindAllMarkers(gmp_carpos_cd8_seurat_obj,
                              min.pct = 0.2,
                              logfc.threshold = 0.2,
                              test.use = "wilcox")
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CARpos_CD8_Clusters_AllMarkers.xlsx"),
              sheetName = "GMP_CARpos_CD8_Clusters_AllMarkers_DE_Result", row.names = FALSE)
  de_result <- FindAllMarkers(pi_carpos_cd8_seurat_obj,
                              min.pct = 0.2,
                              logfc.threshold = 0.2,
                              test.use = "wilcox")
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/PI_CARpos_CD8_Clusters_AllMarkers.xlsx"),
              sheetName = "PI_CARpos_CD8_Clusters_AllMarkers_DE_Result", row.names = FALSE)
  
  ### see feature plot of the markers genes determined by Tay
  effector_markers <- c("LAG3", "PRF1", "GZMB", "GZMK", "GZMH", "IFNG", "TBX21", "STAT4", "KLRD1", "GNLY", "KLRG1")
  treg_markers <- c("FOXP3", "IL2RA", "TIGIT", "CTLA4")
  central_memory_markers <- c("SELL", "IL7R", "EOMES", "CCR7")
  effector_memory_markers <- c("IL7R", effector_markers)
  tissue_resident_memory_markers <- c("CD69", "ITGAE", "CXCR6", "ITGA1")
  stem_cell_memory_markers <- c("SELL", "CCR7", "FAS", "CXCR3", "IL2RB")
  exhaustion_markers <- c("TOX", "PDCD1", "EOMES", "LAG3", "LEF1", "TCF7", "CX3CR1")
  
  ### filter with the existing genes
  effector_markers <- intersect(effector_markers, rownames(sub_seurat_obj4@assays$RNA@counts))
  treg_markers <- intersect(treg_markers, rownames(sub_seurat_obj4@assays$RNA@counts))
  central_memory_markers <- intersect(central_memory_markers, rownames(sub_seurat_obj4@assays$RNA@counts))
  effector_memory_markers <- intersect(effector_memory_markers, rownames(sub_seurat_obj4@assays$RNA@counts))
  tissue_resident_memory_markers <- intersect(tissue_resident_memory_markers, rownames(sub_seurat_obj4@assays$RNA@counts))
  stem_cell_memory_markers <- intersect(stem_cell_memory_markers, rownames(sub_seurat_obj4@assays$RNA@counts))
  exhaustion_markers <- intersect(exhaustion_markers, rownames(sub_seurat_obj4@assays$RNA@counts))
  
  ### See gene expressions on UMAP with the effector markers
  p <- FeaturePlot(gmp_carpos_cd8_seurat_obj, features = effector_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_GMP_CARpos_CD8_Effector_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  p <- FeaturePlot(pi_carpos_cd8_seurat_obj, features = effector_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_PI_CARpos_CD8_Effector_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the tregs markers
  p <- FeaturePlot(gmp_carpos_cd8_seurat_obj, features = treg_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_GMP_CARpos_CD8_Treg_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  p <- FeaturePlot(pi_carpos_cd8_seurat_obj, features = treg_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_PI_CARpos_CD8_Treg_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the central memory markers
  p <- FeaturePlot(gmp_carpos_cd8_seurat_obj, features = central_memory_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_GMP_CARpos_CD8_Central_Memory_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  p <- FeaturePlot(pi_carpos_cd8_seurat_obj, features = central_memory_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_PI_CARpos_CD8_Central_Memory_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the effector memory markers
  p <- FeaturePlot(gmp_carpos_cd8_seurat_obj, features = effector_memory_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_GMP_CARpos_CD8_Effector_Memory_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  p <- FeaturePlot(pi_carpos_cd8_seurat_obj, features = effector_memory_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_PI_CARpos_CD8_Effector_Memory_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the tissue resident memory markers
  p <- FeaturePlot(gmp_carpos_cd8_seurat_obj, features = tissue_resident_memory_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_GMP_CARpos_CD8_Tissue_Resident_Memory_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  p <- FeaturePlot(pi_carpos_cd8_seurat_obj, features = tissue_resident_memory_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_PI_CARpos_CD8_Tissue_Resident_Memory_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the stem cell memory markers
  p <- FeaturePlot(gmp_carpos_cd8_seurat_obj, features = stem_cell_memory_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_GMP_CARpos_CD8_Stem_Cell_Memory_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  p <- FeaturePlot(pi_carpos_cd8_seurat_obj, features = stem_cell_memory_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_PI_CARpos_CD8_Stem_Cell_Memory_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the exhaustion markers
  p <- FeaturePlot(gmp_carpos_cd8_seurat_obj, features = exhaustion_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_GMP_CARpos_CD8_Exhaustion_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  p <- FeaturePlot(pi_carpos_cd8_seurat_obj, features = exhaustion_markers, cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "FeaturePlot_PI_CARpos_CD8_Exhaustion_Markers.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### Ridge plot
  gmp_carpos_cd8_seurat_obj <- SetIdent(object = gmp_carpos_cd8_seurat_obj,
                                        cells = rownames(gmp_carpos_cd8_seurat_obj@meta.data),
                                        value = gmp_carpos_cd8_seurat_obj@meta.data$clusters)
  pi_carpos_cd8_seurat_obj <- SetIdent(object = pi_carpos_cd8_seurat_obj,
                                       cells = rownames(pi_carpos_cd8_seurat_obj@meta.data),
                                       value = pi_carpos_cd8_seurat_obj@meta.data$clusters)
  
  ### ridge plot with effector markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = effector_markers)
  for(i in 1:length(effector_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Effector_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = effector_markers)
  for(i in 1:length(effector_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Effector_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### ridge plot with tregs markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = treg_markers)
  for(i in 1:length(treg_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Treg_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = treg_markers)
  for(i in 1:length(treg_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Treg_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### ridge plot with central memory markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = central_memory_markers)
  for(i in 1:length(central_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Central_Memory_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = central_memory_markers)
  for(i in 1:length(central_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Central_Memory_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### ridge plot with effector memory markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = effector_memory_markers)
  for(i in 1:length(effector_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Effector_Memory_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = effector_memory_markers)
  for(i in 1:length(effector_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Effector_Memory_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### ridge plot with tissue resident memory markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = tissue_resident_memory_markers)
  for(i in 1:length(tissue_resident_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Tissue_Resident_Memory_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = tissue_resident_memory_markers)
  for(i in 1:length(tissue_resident_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Tissue_Resident_Memory_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### ridge plot with stem cell memory markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = stem_cell_memory_markers)
  for(i in 1:length(stem_cell_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Stem_Cell_Memory_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = stem_cell_memory_markers)
  for(i in 1:length(stem_cell_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Stem_Cell_Memory_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### ridge plot with exhaustion markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = exhaustion_markers)
  for(i in 1:length(exhaustion_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Exhaustion_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = exhaustion_markers)
  for(i in 1:length(exhaustion_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Exhaustion_Markers.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### look at all the GMP & PI subsister cells and look at them closely
  ### test some genes (Ridgeplot) between subsisters and non-subsisters
  gmp_carpos_cd8_seurat_obj <- SetIdent(object = gmp_carpos_cd8_seurat_obj,
                                        cells = rownames(gmp_carpos_cd8_seurat_obj@meta.data),
                                        value = gmp_carpos_cd8_seurat_obj@meta.data$GMP_CARpos_CD8_Persister)
  pi_carpos_cd8_seurat_obj <- SetIdent(object = pi_carpos_cd8_seurat_obj,
                                       cells = rownames(pi_carpos_cd8_seurat_obj@meta.data),
                                       value = pi_carpos_cd8_seurat_obj@meta.data$ALL_GMP_CARpos_Persister)
  
  ### ridge plot with effector memory markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = effector_memory_markers)
  for(i in 1:length(effector_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Is_Subsistent")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Effector_Memory_Markers_Subsister.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = effector_memory_markers)
  for(i in 1:length(effector_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Is_Subsistent")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Effector_Memory_Markers_Subsister.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### ridge plot with central memory markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = central_memory_markers, ncol = 2)
  for(i in 1:length(central_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Is_Subsistent")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Central_Memory_Markers_Subsister.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = central_memory_markers, ncol = 2)
  for(i in 1:length(central_memory_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Is_Subsistent")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Central_Memory_Markers_Subsister.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  ### ridge plot with exhaustion markers
  p <- RidgePlot(gmp_carpos_cd8_seurat_obj, features = exhaustion_markers)
  for(i in 1:length(exhaustion_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Is_Subsistent")
  }
  ggsave(paste0(outputDir2, "RidgePlot_GMP_CARpos_CD8_Exhaustion_Markers_Subsister.png"), plot = p, width = 20, height = 15, dpi = 350)
  p <- RidgePlot(pi_carpos_cd8_seurat_obj, features = exhaustion_markers)
  for(i in 1:length(exhaustion_markers)) {
    p[[i]] <- p[[i]] + labs(y = "Is_Subsistent")
  }
  ggsave(paste0(outputDir2, "RidgePlot_PI_CARpos_CD8_Exhaustion_Markers_Subsister.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  
  ### AFTER CHATTING WITH JEREMY,
  ### INCLUDE ALL THE PATIENTS except pxs with insufficient samples & INCREASE THE CLUSTER #
  
  ### get CARpos-only seurat object
  carpos_cd8_cells <- rownames(Seurat_Obj@meta.data)[intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                               which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))]
  sub_seurat_obj4 <- subset(Seurat_Obj, cells = carpos_cd8_cells)
  
  ### pi time points only
  pi_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", 
                      "Wk6", "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$time2)
  sub_seurat_obj4 <- subset(sub_seurat_obj4, idents = intersect(pi_time_points,
                                                                unique(sub_seurat_obj4@meta.data$time2)))
  
  ### get seurat object for some specific patients
  sub_seurat_obj4 <- SetIdent(object = sub_seurat_obj4,
                              cells = rownames(sub_seurat_obj4@meta.data),
                              value = sub_seurat_obj4@meta.data$px)
  sub_seurat_obj4 <- subset(sub_seurat_obj4, idents = c("SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                                                        "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                        "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### normalization
  sub_seurat_obj4 <- NormalizeData(sub_seurat_obj4,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  sub_seurat_obj4 <- FindVariableFeatures(sub_seurat_obj4,
                                          selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  sub_seurat_obj4 <- ScaleData(sub_seurat_obj4,
                               vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run mnn
  sub_seurat_obj4.list <- SplitObject(sub_seurat_obj4, split.by = "library")
  sub_seurat_obj4 <- RunFastMNN(object.list = sub_seurat_obj4.list)
  rm(sub_seurat_obj4.list)
  gc()
  
  ### normalization
  sub_seurat_obj4 <- NormalizeData(sub_seurat_obj4,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  sub_seurat_obj4 <- FindVariableFeatures(sub_seurat_obj4,
                                          selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  sub_seurat_obj4 <- ScaleData(sub_seurat_obj4,
                               vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap & clustering
  sub_seurat_obj4 <- RunPCA(sub_seurat_obj4,
                            features = VariableFeatures(object = sub_seurat_obj4),
                            npcs = 15)
  sub_seurat_obj4 <- RunUMAP(sub_seurat_obj4, reduction = "mnn", dims = 1:15)
  sub_seurat_obj4 <- FindNeighbors(sub_seurat_obj4, reduction = "mnn", dims = 1:15)
  sub_seurat_obj4 <- FindClusters(sub_seurat_obj4)
  
  ### save the clustering result to meta.data
  sub_seurat_obj4@meta.data$clusters <- Idents(sub_seurat_obj4)
  
  ### UMAP with clusters
  p <- DimPlot(object = sub_seurat_obj4, reduction = "umap",
               group.by = "clusters",
               pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "PI_MNN_UMAP_CARpos_CD8_13_Clusters.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### UMAP with clusters
  p <- DimPlot(object = sub_seurat_obj4, reduction = "umap",
               group.by = "px",
               pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "PI_MNN_UMAP_CARpos_CD8_Px.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  #
  ### Jeremy's workflow
  #
  
  ### get CARpos-only seurat object
  carpos_cd8_cells <- rownames(Seurat_Obj@meta.data)[intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                               which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))]
  gmp_carpos_cd8_seurat_obj <- subset(Seurat_Obj, cells = carpos_cd8_cells)
  pi_carpos_cd8_seurat_obj <- subset(Seurat_Obj, cells = carpos_cd8_cells)
  
  ### pi time points only
  gmp_time_points <- c("GMP")
  pi_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", 
                      "Wk6", "Wk8", "3mo", "6mo", "9mo")
  gmp_carpos_cd8_seurat_obj <- SetIdent(object = gmp_carpos_cd8_seurat_obj,
                                        cells = rownames(gmp_carpos_cd8_seurat_obj@meta.data),
                                        value = gmp_carpos_cd8_seurat_obj@meta.data$time2)
  gmp_carpos_cd8_seurat_obj <- subset(gmp_carpos_cd8_seurat_obj, idents = intersect(gmp_time_points,
                                                                                    unique(gmp_carpos_cd8_seurat_obj@meta.data$time2)))
  pi_carpos_cd8_seurat_obj <- SetIdent(object = pi_carpos_cd8_seurat_obj,
                                       cells = rownames(pi_carpos_cd8_seurat_obj@meta.data),
                                       value = pi_carpos_cd8_seurat_obj@meta.data$time2)
  pi_carpos_cd8_seurat_obj <- subset(pi_carpos_cd8_seurat_obj, idents = intersect(pi_time_points,
                                                                                  unique(pi_carpos_cd8_seurat_obj@meta.data$time2)))
  
  ### get seurat object for some specific patients
  gmp_carpos_cd8_seurat_obj <- SetIdent(object = gmp_carpos_cd8_seurat_obj,
                                        cells = rownames(gmp_carpos_cd8_seurat_obj@meta.data),
                                        value = gmp_carpos_cd8_seurat_obj@meta.data$px)
  gmp_carpos_cd8_seurat_obj <- subset(gmp_carpos_cd8_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                                                                            "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                                            "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  pi_carpos_cd8_seurat_obj <- SetIdent(object = pi_carpos_cd8_seurat_obj,
                                       cells = rownames(pi_carpos_cd8_seurat_obj@meta.data),
                                       value = pi_carpos_cd8_seurat_obj@meta.data$px)
  pi_carpos_cd8_seurat_obj <- subset(pi_carpos_cd8_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-03", "SJCAR19-04", "SJCAR19-05",
                                                                          "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                                          "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### normalization
  gmp_carpos_cd8_seurat_obj <- NormalizeData(gmp_carpos_cd8_seurat_obj,
                                             normalization.method = "LogNormalize", scale.factor = 10000)
  pi_carpos_cd8_seurat_obj <- NormalizeData(pi_carpos_cd8_seurat_obj,
                                            normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### Lets pull out Alpha-Beta and Gamma-Delta genes and IG genes
  ### GMP
  TRAVgenes <- grep(pattern = "TRAV", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRBVgenes <- grep(pattern = "TRBV", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRAJgenes <- grep(pattern = "TRAJ", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRBJgenes <- grep(pattern = "TRBJ", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRACgenes <- grep(pattern = "TRAC", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRBCgenes <- grep(pattern = "TRBC", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  
  AlphaBetaGenes <- c(TRAVgenes, TRBVgenes, TRAJgenes, TRBJgenes, TRACgenes, TRBCgenes)
  
  percent.AlphaBeta <- Matrix::colSums(GetAssayData(object = gmp_carpos_cd8_seurat_obj, slot = 'counts')[AlphaBetaGenes,]) / Matrix::colSums(GetAssayData(object = gmp_carpos_cd8_seurat_obj, slot = 'counts'))
  gmp_carpos_cd8_seurat_obj[['percent.AlphaBeta']] <- percent.AlphaBeta
  
  
  TRGVgenes <- grep(pattern = "TRGV", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRDVgenes <- grep(pattern = "TRDV", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRGJgenes <- grep(pattern = "TRGJ", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRDJgenes <- grep(pattern = "TRDJ", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRGCgenes <- grep(pattern = "TRGC", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRGCgenes2 <- grep(pattern = "TCRG", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  TRDCgenes <- grep(pattern = "TRDC", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  
  GammaDeltagenes <- c(TRGVgenes, TRDVgenes, TRGJgenes, TRDVgenes, 
                       TRGCgenes, TRDCgenes, TRGCgenes2)
  
  percent.GammaDelta <- Matrix::colSums(GetAssayData(object = gmp_carpos_cd8_seurat_obj, slot = 'counts')[GammaDeltagenes,]) / Matrix::colSums(GetAssayData(object = gmp_carpos_cd8_seurat_obj, slot = 'counts'))
  gmp_carpos_cd8_seurat_obj[['percent.GammaDelta']] <- percent.GammaDelta
  
  
  IGKVgenes <- grep(pattern = "IGKV", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  IGKJgenes <- grep(pattern = "IGKJ", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  IGKDgenes <- grep(pattern = "IGKD", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  IGKCgenes <- grep(pattern = "IGKC", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  IGJchaingenes <- grep(pattern = "JCHAIN", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  
  IGHVgenes <- grep(pattern = "IGHV", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  IGHJgenes <- grep(pattern = "IGHJ", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  IGHDgenes <- grep(pattern = "IGHD", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  IGHCgenes <- grep(pattern = "IGHC", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  
  IGLVgenes <- grep(pattern = "IGLV", x = rownames(gmp_carpos_cd8_seurat_obj), value = TRUE)
  
  
  IGgenes <- c(IGKVgenes, IGKJgenes, IGKDgenes, IGKCgenes, IGHVgenes, IGHJgenes, IGHDgenes, IGHCgenes, IGJchaingenes,
               IGLVgenes)
  
  percent.IG <- Matrix::colSums(GetAssayData(object = gmp_carpos_cd8_seurat_obj, slot = 'counts')[IGgenes,]) / Matrix::colSums(GetAssayData(object = gmp_carpos_cd8_seurat_obj, slot = 'counts'))
  gmp_carpos_cd8_seurat_obj[['percent.IG']] <- percent.IG
  
  
  gmp_markers.remove <- c(AlphaBetaGenes, GammaDeltagenes, IGgenes)
  
  ### PI
  TRAVgenes <- grep(pattern = "TRAV", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRBVgenes <- grep(pattern = "TRBV", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRAJgenes <- grep(pattern = "TRAJ", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRBJgenes <- grep(pattern = "TRBJ", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRACgenes <- grep(pattern = "TRAC", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRBCgenes <- grep(pattern = "TRBC", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  
  AlphaBetaGenes <- c(TRAVgenes, TRBVgenes, TRAJgenes, TRBJgenes, TRACgenes, TRBCgenes)
  
  percent.AlphaBeta <- Matrix::colSums(GetAssayData(object = pi_carpos_cd8_seurat_obj, slot = 'counts')[AlphaBetaGenes,]) / Matrix::colSums(GetAssayData(object = pi_carpos_cd8_seurat_obj, slot = 'counts'))
  pi_carpos_cd8_seurat_obj[['percent.AlphaBeta']] <- percent.AlphaBeta
  
  
  TRGVgenes <- grep(pattern = "TRGV", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRDVgenes <- grep(pattern = "TRDV", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRGJgenes <- grep(pattern = "TRGJ", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRDJgenes <- grep(pattern = "TRDJ", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRGCgenes <- grep(pattern = "TRGC", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRGCgenes2 <- grep(pattern = "TCRG", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  TRDCgenes <- grep(pattern = "TRDC", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  
  GammaDeltagenes <- c(TRGVgenes, TRDVgenes, TRGJgenes, TRDVgenes, 
                       TRGCgenes, TRDCgenes, TRGCgenes2)
  
  percent.GammaDelta <- Matrix::colSums(GetAssayData(object = pi_carpos_cd8_seurat_obj, slot = 'counts')[GammaDeltagenes,]) / Matrix::colSums(GetAssayData(object = pi_carpos_cd8_seurat_obj, slot = 'counts'))
  pi_carpos_cd8_seurat_obj[['percent.GammaDelta']] <- percent.GammaDelta
  
  
  IGKVgenes <- grep(pattern = "IGKV", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  IGKJgenes <- grep(pattern = "IGKJ", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  IGKDgenes <- grep(pattern = "IGKD", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  IGKCgenes <- grep(pattern = "IGKC", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  IGJchaingenes <- grep(pattern = "JCHAIN", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  
  IGHVgenes <- grep(pattern = "IGHV", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  IGHJgenes <- grep(pattern = "IGHJ", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  IGHDgenes <- grep(pattern = "IGHD", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  IGHCgenes <- grep(pattern = "IGHC", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  
  IGLVgenes <- grep(pattern = "IGLV", x = rownames(pi_carpos_cd8_seurat_obj), value = TRUE)
  
  
  IGgenes <- c(IGKVgenes, IGKJgenes, IGKDgenes, IGKCgenes, IGHVgenes, IGHJgenes, IGHDgenes, IGHCgenes, IGJchaingenes,
               IGLVgenes)
  
  percent.IG <- Matrix::colSums(GetAssayData(object = pi_carpos_cd8_seurat_obj, slot = 'counts')[IGgenes,]) / Matrix::colSums(GetAssayData(object = pi_carpos_cd8_seurat_obj, slot = 'counts'))
  pi_carpos_cd8_seurat_obj[['percent.IG']] <- percent.IG
  
  
  pi_markers.remove <- c(AlphaBetaGenes, GammaDeltagenes, IGgenes)
  
  ### Detect variable features
  ### Will use 'vst' method, which Seurat3 claims is superior to the mean.var.plot method
  gmp_carpos_cd8_seurat_obj <- FindVariableFeatures(object = gmp_carpos_cd8_seurat_obj, selection.method = 'vst')
  pi_carpos_cd8_seurat_obj <- FindVariableFeatures(object = pi_carpos_cd8_seurat_obj, selection.method = 'vst')
  
  ##Let's get rid of the TCR genes from Variable Genes so that clonotype info isn't affecting anything at this stage
  #This is important because of the way that the alignment works. If a read maps to more than one gene equally well,
  #then it is discounted entirely. Some TCR gene segments are more distinct than others, so this could cause bias.
  NumFeaturesToRemove <- length(VariableFeatures(gmp_carpos_cd8_seurat_obj)[VariableFeatures(gmp_carpos_cd8_seurat_obj) %in% gmp_markers.remove])
  gmp_carpos_cd8_seurat_obj <- FindVariableFeatures(object = gmp_carpos_cd8_seurat_obj, selection.method = 'vst', nfeatures = (2000 +NumFeaturesToRemove) )
  gmp_carpos_cd8_seurat_obj@assays$RNA@var.features <- gmp_carpos_cd8_seurat_obj@assays$RNA@var.features[!(gmp_carpos_cd8_seurat_obj@assays$RNA@var.features %in% gmp_markers.remove)]
  
  NumFeaturesToRemove <- length(VariableFeatures(pi_carpos_cd8_seurat_obj)[VariableFeatures(pi_carpos_cd8_seurat_obj) %in% pi_markers.remove])
  pi_carpos_cd8_seurat_obj <- FindVariableFeatures(object = pi_carpos_cd8_seurat_obj, selection.method = 'vst', nfeatures = (2000 +NumFeaturesToRemove) )
  pi_carpos_cd8_seurat_obj@assays$RNA@var.features <- pi_carpos_cd8_seurat_obj@assays$RNA@var.features[!(pi_carpos_cd8_seurat_obj@assays$RNA@var.features %in% pi_markers.remove)]
  
  ### exclude the weird cluster (13) in PI
  
  
  
  ### run mnn
  gmp_carpos_cd8_seurat_obj.list <- SplitObject(gmp_carpos_cd8_seurat_obj, split.by = "library")
  gmp_carpos_cd8_seurat_obj <- RunFastMNN(object.list = gmp_carpos_cd8_seurat_obj.list)
  rm(gmp_carpos_cd8_seurat_obj.list)
  gc()
  
  pi_carpos_cd8_seurat_obj.list <- SplitObject(pi_carpos_cd8_seurat_obj, split.by = "library")
  pi_carpos_cd8_seurat_obj <- RunFastMNN(object.list = pi_carpos_cd8_seurat_obj.list)
  rm(pi_carpos_cd8_seurat_obj.list)
  gc()
  
  ### run umap and find clusters
  gmp_carpos_cd8_seurat_obj <- RunUMAP(gmp_carpos_cd8_seurat_obj, reduction = "mnn", dims = 1:50)
  gmp_carpos_cd8_seurat_obj <- FindNeighbors(gmp_carpos_cd8_seurat_obj, reduction = "mnn", dims = 1:50)
  gmp_carpos_cd8_seurat_obj <- FindClusters(gmp_carpos_cd8_seurat_obj)
  
  pi_carpos_cd8_seurat_obj <- RunUMAP(pi_carpos_cd8_seurat_obj, reduction = "mnn", dims = 1:50)
  pi_carpos_cd8_seurat_obj <- FindNeighbors(pi_carpos_cd8_seurat_obj, reduction = "mnn", dims = 1:50)
  pi_carpos_cd8_seurat_obj <- FindClusters(pi_carpos_cd8_seurat_obj)
  
  ### save the clustering result to meta.data
  gmp_carpos_cd8_seurat_obj@meta.data$clusters <- Idents(gmp_carpos_cd8_seurat_obj)
  pi_carpos_cd8_seurat_obj@meta.data$clusters <- Idents(pi_carpos_cd8_seurat_obj)
  
  ### UMAP with clusters
  p <- DimPlot(object = gmp_carpos_cd8_seurat_obj, reduction = "umap",
               group.by = "clusters",
               pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "GMP_MNN_UMAP_CARpos_CD8_Clusters.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- DimPlot(object = pi_carpos_cd8_seurat_obj, reduction = "umap",
               group.by = "clusters",
               pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "PI_MNN_UMAP_CARpos_CD8_Clusters.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### find doublets & multiplets
  NumFeaturesToRemove <- length(VariableFeatures(gmp_carpos_cd8_seurat_obj)[VariableFeatures(gmp_carpos_cd8_seurat_obj) %in% gmp_markers.remove])
  gmp_carpos_cd8_seurat_obj <- FindVariableFeatures(object = gmp_carpos_cd8_seurat_obj, selection.method = 'vst', nfeatures = (2000 +NumFeaturesToRemove) )
  gmp_carpos_cd8_seurat_obj@assays$RNA@var.features <- gmp_carpos_cd8_seurat_obj@assays$RNA@var.features[!(gmp_carpos_cd8_seurat_obj@assays$RNA@var.features %in% gmp_markers.remove)]
  
  NumFeaturesToRemove <- length(VariableFeatures(pi_carpos_cd8_seurat_obj)[VariableFeatures(pi_carpos_cd8_seurat_obj) %in% pi_markers.remove])
  pi_carpos_cd8_seurat_obj <- FindVariableFeatures(object = pi_carpos_cd8_seurat_obj, selection.method = 'vst', nfeatures = (2000 +NumFeaturesToRemove) )
  pi_carpos_cd8_seurat_obj@assays$RNA@var.features <- pi_carpos_cd8_seurat_obj@assays$RNA@var.features[!(pi_carpos_cd8_seurat_obj@assays$RNA@var.features %in% pi_markers.remove)]
  
  set.seed(1234)
  gmp_carpos_cd8_seurat_obj@meta.data$doublet_score <- computeDoubletDensity(as.SingleCellExperiment(gmp_carpos_cd8_seurat_obj), subset.row=gmp_carpos_cd8_seurat_obj@assays$RNA@var.features)
  pi_carpos_cd8_seurat_obj@meta.data$doublet_score <- computeDoubletDensity(as.SingleCellExperiment(pi_carpos_cd8_seurat_obj), subset.row=pi_carpos_cd8_seurat_obj@assays$RNA@var.features)
  
  p <- FeaturePlot(gmp_carpos_cd8_seurat_obj, features = "doublet_score", cols = c("lightgray", "red")) +
    ggtitle("") +
    labs(color="Doublet Score") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "GMP_MNN_UMAP_CARpos_CD8_Doublet_Score.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- FeaturePlot(pi_carpos_cd8_seurat_obj, features = "doublet_score", cols = c("lightgray", "red")) +
    ggtitle("") +
    labs(color="Doublet Score") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "PI_MNN_UMAP_CARpos_CD8_Doublet_Score.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### cluster markers
  gmp_carpos_cd8_seurat_obj <- SetIdent(object = gmp_carpos_cd8_seurat_obj,
                                        cells = rownames(gmp_carpos_cd8_seurat_obj@meta.data),
                                        value = gmp_carpos_cd8_seurat_obj@meta.data$clusters)
  de_result <- FindAllMarkers(gmp_carpos_cd8_seurat_obj,
                              min.pct = 0.2,
                              logfc.threshold = 0.2,
                              test.use = "wilcox")
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/GMP_CARpos_CD8_Clusters_AllMarkers2.xlsx"),
              sheetName = "GMP_CARpos_CD8_Clusters_AllMarkers_DE_Result", row.names = FALSE)
  
  pi_carpos_cd8_seurat_obj <- SetIdent(object = pi_carpos_cd8_seurat_obj,
                                       cells = rownames(pi_carpos_cd8_seurat_obj@meta.data),
                                       value = pi_carpos_cd8_seurat_obj@meta.data$clusters)
  de_result <- FindAllMarkers(pi_carpos_cd8_seurat_obj,
                              min.pct = 0.2,
                              logfc.threshold = 0.2,
                              test.use = "wilcox")
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/PI_CARpos_CD8_Clusters_AllMarkers2.xlsx"),
              sheetName = "PI_CARpos_CD8_Clusters_AllMarkers_DE_Result", row.names = FALSE)
  
  
  #
  ### 25. 06/14/2021 - Re-analyze everything with the PB-Px filtered data with different CD4/CD8 definition
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/25/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### load the seurat object that I gave to Jeremy
  final_seurat_obj <- readRDS(file = "Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Jun2021_PB_Px_Filtered.RDS")
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  final_seurat_obj@meta.data <- final_seurat_obj@meta.data[colnames(final_seurat_obj@assays$RNA@counts),]
  print(identical(rownames(final_seurat_obj@meta.data), colnames(final_seurat_obj@assays$RNA@counts)))
  
  ### active assay = "RNA"
  final_seurat_obj@active.assay <- "RNA"
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = final_seurat_obj)), rownames(final_seurat_obj@meta.data)))
  
  
  #
  ### 26. New Task - 3D plane mapping between GMP & PI clusters
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/26/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### load the seurat object that I gave to Jeremy
  final_seurat_obj <- readRDS(file = "Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Jun2021_PB_Px_Filtered.RDS")
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  final_seurat_obj@meta.data <- final_seurat_obj@meta.data[colnames(final_seurat_obj@assays$RNA@counts),]
  print(identical(rownames(final_seurat_obj@meta.data), colnames(final_seurat_obj@assays$RNA@counts)))
  
  ### active assay = "RNA"
  final_seurat_obj@active.assay <- "RNA"
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = final_seurat_obj)), rownames(final_seurat_obj@meta.data)))
  
  ### only use the time points we are interested
  gmp_after_time_points <- c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", 
                             "Wk6", "Wk8", "3mo", "6mo", "9mo")
  final_seurat_obj <- SetIdent(object = final_seurat_obj,
                               cells = rownames(final_seurat_obj@meta.data),
                               value = final_seurat_obj@meta.data$time2)
  final_seurat_obj <- subset(final_seurat_obj, idents = intersect(gmp_after_time_points,
                                                                  unique(final_seurat_obj@meta.data$time2)))
  
  ### add GMP - PI distinguishable column
  final_seurat_obj@meta.data$GMP_PI <- "PI"
  final_seurat_obj@meta.data$GMP_PI[which(final_seurat_obj@meta.data$time2 == "GMP")] <- "GMP"
  
  ### normalization
  final_seurat_obj <- NormalizeData(final_seurat_obj,
                                    normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  final_seurat_obj <- FindVariableFeatures(final_seurat_obj,
                                           selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  final_seurat_obj <- ScaleData(final_seurat_obj,
                                vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run mnn
  final_seurat_obj.list <- SplitObject(final_seurat_obj, split.by = "library")
  final_seurat_obj <- RunFastMNN(object.list = final_seurat_obj.list)
  rm(final_seurat_obj.list)
  gc()
  
  ### normalization
  final_seurat_obj <- NormalizeData(final_seurat_obj,
                                    normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  final_seurat_obj <- FindVariableFeatures(final_seurat_obj,
                                           selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  final_seurat_obj <- ScaleData(final_seurat_obj,
                                vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run umap and find clusters
  final_seurat_obj <- RunUMAP(final_seurat_obj, reduction = "mnn", dims = 1:30)
  final_seurat_obj <- FindNeighbors(final_seurat_obj, reduction = "mnn", dims = 1:30)
  final_seurat_obj <- FindClusters(final_seurat_obj)
  
  ### save the clustering result to meta.data
  final_seurat_obj@meta.data$clusters <- Idents(final_seurat_obj)
  
  # ### temporary save
  # saveRDS(object = final_seurat_obj, file = "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/final_seurat_obj_mnn_performed.rds")
  
  ### only pick CARpos cells
  final_seurat_obj <- subset(final_seurat_obj, cells = rownames(final_seurat_obj@meta.data)[which(final_seurat_obj@meta.data$CAR == "CARpos")])
  
  ### color generation for the clusters
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  color_scale <- gg_color_hue(length(levels(final_seurat_obj@meta.data$clusters)))
  names(color_scale) <- levels(final_seurat_obj@meta.data$clusters)
  show_col(color_scale)
  
  ### Combined UMAP plot
  p <- DimPlot(object = final_seurat_obj, reduction = "umap",
               group.by = "clusters",
               cols = color_scale[as.character(final_seurat_obj@meta.data$clusters)],
               pt.size = 1, label = TRUE) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "MNN_UMAP_CARpos_Combined_Clusters.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### split gmp & pi umap
  gmp_umap <- Embeddings(final_seurat_obj, reduction = "umap")[rownames(final_seurat_obj@meta.data)[which(final_seurat_obj@meta.data$GMP_PI == "GMP")], 1:2]
  pi_umap <- Embeddings(final_seurat_obj, reduction = "umap")[rownames(final_seurat_obj@meta.data)[which(final_seurat_obj@meta.data$GMP_PI == "PI")], 1:2]
  
  ### draw in 3D
  shift_x <- 10
  shift_z <- 5
  plot_df <- data.frame(x=c(gmp_umap[,1], pi_umap[,1]+shift_x),
                        y=c(gmp_umap[,2], pi_umap[,2]),
                        z=c(rep(0, nrow(gmp_umap)), rep(shift_z, nrow(pi_umap))),
                        cluster=as.character(c(final_seurat_obj@meta.data[rownames(gmp_umap),"clusters"],
                                               final_seurat_obj@meta.data[rownames(pi_umap),"clusters"])),
                        clone=as.character(c(final_seurat_obj@meta.data[rownames(gmp_umap),"clonotype_id_by_patient_one_alpha_beta"],
                                             final_seurat_obj@meta.data[rownames(pi_umap),"clonotype_id_by_patient_one_alpha_beta"])),
                        color=c(rep("red", nrow(gmp_umap)), rep("blue", nrow(pi_umap))),
                        gmp_pi=c(rep("GMP", nrow(gmp_umap)), rep("PI", nrow(pi_umap))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df$color <- color_scale[as.character(plot_df$cluster)]
  rownames(plot_df) <- c(rownames(gmp_umap), rownames(pi_umap))
  
  ### find connections between GMP & PI
  gmp_pi_connection_clones <- intersect(final_seurat_obj@meta.data[rownames(gmp_umap),"clonotype_id_by_patient_one_alpha_beta"],
                                        final_seurat_obj@meta.data[rownames(pi_umap),"clonotype_id_by_patient_one_alpha_beta"])
  gmp_pi_connection_clones <- gmp_pi_connection_clones[which(!is.na(gmp_pi_connection_clones))]
  
  ### define the subsister column
  ### one gmp -> one pi
  gmp_subsisters_idx <- which(final_seurat_obj@meta.data[rownames(gmp_umap),"clonotype_id_by_patient_one_alpha_beta"] %in% gmp_pi_connection_clones)
  temp <- final_seurat_obj@meta.data[rownames(gmp_umap),"clonotype_id_by_patient_one_alpha_beta"][gmp_subsisters_idx]
  gmp_subsisters_name <- rownames(gmp_umap)[gmp_subsisters_idx[which(!duplicated(temp))]]
  names(gmp_subsisters_name) <- final_seurat_obj@meta.data[gmp_subsisters_name,"clonotype_id_by_patient_one_alpha_beta"]
  gmp_subsisters_name <- gmp_subsisters_name[order(names(gmp_subsisters_name))]
  
  pi_subsisters_idx <- which(final_seurat_obj@meta.data[rownames(pi_umap),"clonotype_id_by_patient_one_alpha_beta"] %in% gmp_pi_connection_clones)
  temp <- final_seurat_obj@meta.data[rownames(pi_umap),"clonotype_id_by_patient_one_alpha_beta"][pi_subsisters_idx]
  pi_subsisters_name <- rownames(pi_umap)[pi_subsisters_idx[which(!duplicated(temp))]]
  names(pi_subsisters_name) <- final_seurat_obj@meta.data[pi_subsisters_name,"clonotype_id_by_patient_one_alpha_beta"]
  pi_subsisters_name <- pi_subsisters_name[order(names(pi_subsisters_name))]
  
  plot_df2 <- data.frame(x0 = plot_df[gmp_subsisters_name, "x"],
                         y0 = plot_df[gmp_subsisters_name, "y"],
                         z0 = plot_df[gmp_subsisters_name, "z"]-0.1,
                         x1 = plot_df[pi_subsisters_name, "x"],
                         y1 = plot_df[pi_subsisters_name, "y"],
                         z1 = plot_df[pi_subsisters_name, "z"]+0.1,
                         stringsAsFactors = FALSE, check.names = FALSE)
  set.seed(1234)
  plot_df2 <- plot_df2[sample(nrow(plot_df2), 5),]
  plot_df2$color <- c(rep("black", 3), rep("black", 2))
  
  ### 2D
  png(filename = paste0(outputDir2, "/", "MNN_UMAP_CARpos_GMP_PI_Mapping_2D.png"), width = 2500, height = 1500, res = 350)
  plot(c(gmp_umap[,1], pi_umap[,1]+50), plot_df$y, col = plot_df$color, pch = 19,
       main = "GMP-PI MNN-UMAP in 2D", xlab = "UMAP1", ylab = "UMAP2")
  dev.off()
  
  ### 3D
  png(filename = paste0(outputDir2, "/", "MNN_UMAP_CARpos_GMP_PI_Mapping_3D.png"), width = 2000, height = 1500, res = 350)
  scatter3D(plot_df$x, plot_df$y, plot_df$z,
            colvar = NULL,
            col = plot_df$color,
            pch = ".",  theta = 0, phi = 0,
            main = "GMP-PI MNN-UMAP IN 3D", xlab = "UMAP1",
            ylab ="UMAP2", zlab = "",
            colkey = FALSE, bty = "n")
  arrows3D(plot_df2$x0, plot_df2$y0, plot_df2$z0,
           plot_df2$x1, plot_df2$y1, plot_df2$z1,
           colvar = NULL, col = plot_df2$color,
           lwd = 1, d = 1,
           main = "", ticktype = "simple",
           add = TRUE, bty = "n")
  dev.off()
  
  ### make it interactive
  plotrgl(lighting = TRUE, smooth = TRUE)
  htmlwidgets::saveWidget(rglwidget(width = 1200, height = 1200), 
                          file = paste0(outputDir2, "MNN_UMAP_CARpos_GMP_PI_Mapping_3D_INTERACTIVE.html"),
                          selfcontained = TRUE)
  
  
  #
  ### 27. 06/21/2021 - Pseudotime (Slingshot & Monocle2) analysis on Jeremy's object
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/27/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### load Jeremy's object
  JCC_Seurat_Obj <- readRDS(file = "./data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds")
  
  ### check whether the orders are the same
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
  print(identical(names(Idents(object = JCC_Seurat_Obj)), rownames(JCC_Seurat_Obj@meta.data)))
  
  ### draw a UMAP to confirm that the UMAP I see is the same as his
  DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
          group.by = "AllSeuratClusters", label = TRUE,
          pt.size = 0.5)
  
  #
  ### perform monocle2 first (since this is faster than the slingshot)
  #
  
  ### cell # for each time point
  cell_num <- sapply(unique(JCC_Seurat_Obj@meta.data$time2), function(x) {
    return(length(which(JCC_Seurat_Obj@meta.data$time2 == x)))
  })
  
  ### remove wk6 & 6mo data since there are only 1 & 7 cells
  cell_num <- cell_num[-c(6,9)]
  
  ### the data is too big (can't run for monocle), so we down-sample the cells in each time point
  set.seed(1234)
  fixed_min_cell_num <- min(cell_num)
  JCC_Seurat_Obj@meta.data$downsampled <- "NO"
  for(tp in names(cell_num)) {
    JCC_Seurat_Obj@meta.data$downsampled[sample(which(JCC_Seurat_Obj@meta.data$time2 == tp), fixed_min_cell_num)] <- "YES"
  }
  
  ### Construct a monocle cds
  monocle_metadata <- JCC_Seurat_Obj@meta.data[rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj@meta.data$downsampled == "YES")],]
  monocle_metadata$time2 <- factor(monocle_metadata$time2, levels = unique(monocle_metadata$time2))
  monocle_cds <- newCellDataSet(as(as.matrix(JCC_Seurat_Obj@assays$RNA@data[,rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj@meta.data$downsampled == "YES")]]), 'sparseMatrix'),
                                phenoData = new('AnnotatedDataFrame', data = monocle_metadata),
                                featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(JCC_Seurat_Obj@assays$RNA@data),
                                                                                          row.names = row.names(JCC_Seurat_Obj@assays$RNA@data),
                                                                                          stringsAsFactors = FALSE, check.names = FALSE)),
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  
  ### run monocle
  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)
  monocle_cds <- reduceDimension(monocle_cds, reduction_method = "DDRTree")
  monocle_cds <- orderCells(monocle_cds)
  
  ### determine the beginning state
  plot_cell_trajectory(monocle_cds, color_by = "time2")
  plot_cell_trajectory(monocle_cds, color_by = "State")
  plot_complex_cell_trajectory(monocle_cds, color_by = "time2")
  plot_complex_cell_trajectory(monocle_cds, color_by = "State")
  monocle_cds <- orderCells(monocle_cds, root_state = "6")
  
  ### draw monocle plots
  p <- plot_cell_trajectory(monocle_cds, color_by = "time2", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_Time_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 16))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_Pseudotime_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_cell_trajectory(monocle_cds, color_by = "State", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_State_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_cell_trajectory(monocle_cds, color_by = "AllSeuratClusters", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_Cluster_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_cell_trajectory(monocle_cds, color_by = "AllSeuratClusters", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30)) +
    facet_wrap(~AllSeuratClusters, nrow = 5, ncol = 6)
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_Cluster2_Monocle2.png"),
         plot = p,
         width = 15, height = 20, dpi = 350)
  
  p <- plot_complex_cell_trajectory(monocle_cds, color_by = "time2", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_Time_Complex_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_complex_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 10))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_Pseudotime_Complex_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_complex_cell_trajectory(monocle_cds, color_by = "State", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_State_Complex_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_complex_cell_trajectory(monocle_cds, color_by = "AllSeuratClusters", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 30))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_Cluster_Complex_Monocle2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  ### get DE genes based on pseudotime
  degs <- differentialGeneTest(monocle_cds,
                               fullModelFormulaStr = "~sm.ns(Pseudotime)",
                               cores = 4)
  degs <- degs[order(degs$qval),]
  
  ### write out the DE genes
  write.xlsx2(data.frame(Gene=rownames(degs),
                         degs,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_Trajectory_Inference_Pseudotime_DEGs.xlsx"),
              sheetName = "CARpos_Trajectory_Inference_Pseudotime_DEGs", row.names = FALSE)
  
  ### gene expression plots
  interesting_genes <- degs$gene_short_name[1:9]
  cds_subset <- monocle_cds[interesting_genes,]
  p <- plot_genes_in_pseudotime(cds_subset, color_by = "time2",
                                cell_size = 5, ncol = 3) +
    labs(color = "Time") +
    theme_classic(base_size = 20) +
    theme(legend.title = element_text(size = 30),
          legend.text = element_text(size = 20))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_Pseudotime_DEGs.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  p <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime",
                                cell_size = 5, ncol = 3) +
    labs(color = "Time") +
    theme_classic(base_size = 20) +
    theme(legend.title = element_text(size = 30),
          legend.text = element_text(size = 20))
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_Pseudotime_DEGs2.png"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  #
  ### Now it's time for Slingshot
  #
  
  ### MNN map
  mnn_map <- Embeddings(JCC_Seurat_Obj, reduction = "mnn")[rownames(JCC_Seurat_Obj@meta.data),1:2]
  
  ### remove wk6 & 6mo data since there are only 1 & 7 cells
  usable_idx <- which(JCC_Seurat_Obj@meta.data$time2 %in% c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", "Wk8", "3mo"))
  mnn_map <- mnn_map[rownames(JCC_Seurat_Obj@meta.data)[usable_idx],]
  mnn_map_time <- JCC_Seurat_Obj@meta.data$time2[usable_idx]
  mnn_map_exp <- JCC_Seurat_Obj@assays$RNA@counts[,rownames(JCC_Seurat_Obj@meta.data)[usable_idx]]
  
  ### get slingshot object
  slingshot_obj_mnn <- slingshot(mnn_map,
                                 clusterLabels = mnn_map_time,
                                 start.clus = "GMP",
                                 end.clus = "3mo")
  slingshot_obj_mnn <- as.SlingshotDataSet(slingshot_obj_mnn)
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(mnn_map_time), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, "CARpos_Trajectory_Inference_Time_MNN.png"), width = 5000, height = 3000, res = 350)
  par(mar=c(7, 7, 7, 1), mgp=c(4,1,0))
  plot(reducedDim(slingshot_obj_mnn),
       main=paste("CAR+ Trajectory Inference"),
       col = cell_colors_clust[as.character(mnn_map_time)],
       pch = 19, cex = 1, cex.lab = 3, cex.main = 3, cex.axis = 2)
  lines(slingshot_obj_mnn, lwd = 4, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, cex = 1.5)
  dev.off()
  
  ### the matrix is too large
  ### so we are down-sampling them
  mnn_map <- mnn_map[rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj@meta.data$downsampled == "YES")],]
  mnn_map_time <- JCC_Seurat_Obj@meta.data$time2[which(JCC_Seurat_Obj@meta.data$downsampled == "YES")]
  mnn_map_exp <- JCC_Seurat_Obj@assays$RNA@counts[,rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj@meta.data$downsampled == "YES")]]
  
  ### get slingshot object
  slingshot_obj_mnn <- slingshot(mnn_map,
                                 clusterLabels = mnn_map_time,
                                 start.clus = "GMP",
                                 end.clus = "3mo")
  slingshot_obj_mnn <- as.SlingshotDataSet(slingshot_obj_mnn)
  
  ### Trajectory inference
  png(paste0(outputDir2, "CARpos_Trajectory_Inference_Time_MNN_DOWNSAMPLED.png"), width = 5000, height = 3000, res = 350)
  par(mar=c(7, 7, 7, 1), mgp=c(4,1,0))
  plot(reducedDim(slingshot_obj_mnn),
       main=paste("CAR+ Trajectory Inference"),
       col = cell_colors_clust[as.character(mnn_map_time)],
       pch = 19, cex = 1, cex.lab = 3, cex.main = 3, cex.axis = 2)
  lines(slingshot_obj_mnn, lwd = 4, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, cex = 1.5)
  dev.off()
  
  ### find genes that change their expression over the course of development
  ###  If "consecutive", then consecutive points along each lineage will be used as contrasts
  sce <- fitGAM(counts = mnn_map_exp, sds = slingshot_obj_mnn)
  ATres <- associationTest(sce, contrastType = "consecutive")
  
  ### give FDR
  ATres <- ATres[order(ATres$pvalue),]
  ATres$FDR <- p.adjust(p = ATres$pvalue, method = "BH")
  
  ### save the result
  write.xlsx2(data.frame(Gene=rownames(ATres),
                         ATres,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_Trajectory_Inference_Pseudotime_DEGs_Slingshot.xlsx"),
              sheetName = "CARpos_Trajectory_Inference_Pseudotime_DEGs_Slingshot", row.names = FALSE)
  
  ### get those genes from the Slingshot
  topgenes <- rownames(ATres[order(ATres$FDR), ])[1:100]
  pst.ord <- order(sce$slingshot$pseudotime.Lineage1, na.last = NA)
  heatdata <- assays(sce)$counts[topgenes, pst.ord]
  heatclus <- JCC_Seurat_Obj@meta.data[colnames(heatdata),"time2"]
  
  ### draw the heatmap
  png(paste0(outputDir2, "CARpos_Trajectory_Inference_Pseudotime_DEGs_Slingshot_Heatmap.png"),
      width = 3000, height = 3000, res = 350)
  par(oma=c(0,3,0,3))
  heatmap.2(log1p(heatdata), col = colorpanel(24, low = "blue", high = "red"),
            scale = "none", dendrogram = "row", trace = "none",
            cexRow = 0.5, key.title = "", main = "Top 100 Genes Associated With The Pseudotime",
            Colv = FALSE, labCol = FALSE,  key.xlab = "log(Count+1)", key.ylab = "Frequency",
            ColSideColors = cell_colors_clust[heatclus])
  legend("left", inset = -0.1,
         box.lty = 0, cex = 0.8,
         title = "Time", xpd = TRUE,
         legend=names(cell_colors_clust),
         col=cell_colors_clust,
         pch=15)
  dev.off()
  
  ### compare DEGs of Monocle2 and Slingshot
  print(length(intersect(degs$gene_short_name[which(degs$qval < 0.01)],
                         rownames(ATres)[which(ATres$FDR < 0.01)])))
  print(length(which(degs$qval < 0.01)))
  print(length(which(ATres$FDR < 0.01)))
  
  
  ### draw a pie chart to show the percentage of time points in each state
  plot_df <- data.frame(State=as.character(sapply(levels(monocle_cds@phenoData@data$State), function(x) rep(x, length(levels(monocle_cds@phenoData@data$time2))))),
                        Time_Point=rep(levels(monocle_cds@phenoData@data$time2), length(levels(monocle_cds@phenoData@data$State))),
                        Numbers=0,
                        Pcnt=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### calculate numbers
  cnt <- 1
  for(state in levels(monocle_cds@phenoData@data$State)) {
    for(tp in levels(monocle_cds@phenoData@data$time2)) {
      plot_df$Numbers[cnt] <- length(intersect(which(monocle_cds@phenoData@data$State == state),
                                               which(monocle_cds@phenoData@data$time2 == tp)))
      cnt <- cnt + 1
    }
  }
  
  ### calculate percentages
  state_sum <- rep(0, length(levels(monocle_cds$State)))
  names(state_sum) <- levels(monocle_cds$State)
  for(i in 1:length(levels(monocle_cds$State))) {
    state_sum[i] <- sum(plot_df[which(plot_df$State == levels(monocle_cds$State)[i]),"Numbers"])
    plot_df$Pcnt[which(plot_df$State == levels(monocle_cds$State)[i])] <- round(plot_df$Numbers[which(plot_df$State == levels(monocle_cds$State)[i])] * 100 / state_sum[i], 1)
  }
  
  ### remove NaN rows
  plot_df <- plot_df[which(!is.nan(plot_df$Pcnt)),]
  
  ### factorize the time point & state
  plot_df$Time_Point <- factor(plot_df$Time_Point, levels = levels(monocle_cds$time2))
  plot_df$State <- factor(plot_df$State, levels = unique(plot_df$State))
  
  ### draw the pie charts
  plot_df <- plot_df[order(as.numeric(plot_df$State)),]
  p <- vector("list", length = length(levels(plot_df$State)))
  names(p) <- levels(plot_df$State)
  for(state in levels(plot_df$State)) {
    p[[state]] <- ggplot(data = plot_df[which(plot_df$State == state),],
                aes(x = "", y = Numbers, fill = Time_Point)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta="y") +
      labs(x = NULL, y = NULL, title = paste("State", state)) +
      scale_fill_discrete(name = "Time", labels = paste0(as.character(plot_df$Time_Point[which(plot_df$State == state)]), ": ",
                                                          plot_df$Numbers[which(plot_df$State == state)], " (",
                                                          plot_df$Pcnt[which(plot_df$State == state)], "%)")) +
      theme_classic(base_size = 48) +
      theme(plot.title = element_text(hjust = 0.5, color = "black", size = 48),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank())
  }
  
  ### combine plots into one
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   ncol = 2,
                   top = "")
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_DOWNSAMPLED_Monocle_Pie.png"), g, width = 25, height = 20, dpi = 350)
  
  ### draw a proportional bar plot
  ### pcnt < 10 -> ""
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) < 10)] <- ""
  plot_df$Pcnt <- as.character(plot_df$Pcnt)
  ### state_sum < 100 -> ""
  blank_state <- names(state_sum)[which(state_sum < 100)]
  plot_df$Pcnt[which(plot_df$State %in% blank_state)] <- ""
  p <- ggplot(data=plot_df, aes_string(x="State", y="Numbers", fill="Time_Point", label="Pcnt")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Proportion of Cells") +
    xlab("State") + ylab("Cell #") +
    geom_text(size = 5, position = position_stack(vjust = 1)) +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_DOWNSAMPLED_Monocle_Bar.png"), plot = p,
         width = 20, height = 10, dpi = 350)
  
  
  ### Now in the inverted way
  ### what states are in each time point
  
  ### draw a pie chart to show the percentage of time points in each state
  plot_df <- data.frame(State=as.character(sapply(levels(monocle_cds@phenoData@data$State), function(x) rep(x, length(levels(monocle_cds@phenoData@data$time2))))),
                        Time_Point=rep(levels(monocle_cds@phenoData@data$time2), length(levels(monocle_cds@phenoData@data$State))),
                        Numbers=0,
                        Pcnt=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### calculate numbers
  cnt <- 1
  for(state in levels(monocle_cds@phenoData@data$State)) {
    for(tp in levels(monocle_cds@phenoData@data$time2)) {
      plot_df$Numbers[cnt] <- length(intersect(which(monocle_cds@phenoData@data$State == state),
                                               which(monocle_cds@phenoData@data$time2 == tp)))
      cnt <- cnt + 1
    }
  }
  
  ### calculate percentages
  tp_sum <- rep(0, length(levels(monocle_cds$time2)))
  names(tp_sum) <- levels(monocle_cds$time2)
  for(i in 1:length(levels(monocle_cds$time2))) {
    tp_sum[i] <- sum(plot_df[which(plot_df$Time_Point == levels(monocle_cds$time2)[i]),"Numbers"])
    plot_df$Pcnt[which(plot_df$Time_Point == levels(monocle_cds$time2)[i])] <- round(plot_df$Numbers[which(plot_df$Time_Point == levels(monocle_cds$time2)[i])] * 100 / tp_sum[i], 1)
  }
  
  ### remove NaN rows
  plot_df <- plot_df[which(!is.nan(plot_df$Pcnt)),]
  
  ### factorize the time point & state
  plot_df$Time_Point <- factor(plot_df$Time_Point, levels = levels(monocle_cds$time2))
  plot_df$State <- factor(plot_df$State, levels = unique(plot_df$State))
  
  ### draw the pie charts
  plot_df <- plot_df[order(plot_df$Time_Point),]
  p <- vector("list", length = length(levels(plot_df$Time_Point)))
  names(p) <- levels(plot_df$Time_Point)
  for(tp in levels(plot_df$Time_Point)) {
    p[[tp]] <- ggplot(data = plot_df[which(plot_df$Time_Point == tp),],
                         aes(x = "", y = Numbers, fill = State)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta="y") +
      labs(x = NULL, y = NULL, title = paste("", tp)) +
      scale_fill_discrete(name = "State", labels = paste0(as.character(plot_df$State[which(plot_df$Time_Point == tp)]), ": ",
                                                          plot_df$Numbers[which(plot_df$Time_Point == tp)], " (",
                                                          plot_df$Pcnt[which(plot_df$Time_Point == tp)], "%)")) +
      theme_classic(base_size = 48) +
      theme(plot.title = element_text(hjust = 0.5, color = "black", size = 48),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank())
  }
  
  ### combine plots into one
  g <- arrangeGrob(grobs = p,
                   nrow = 4,
                   ncol = 2,
                   top = "")
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_DOWNSAMPLED_Monocle_Pie2.png"), g, width = 25, height = 22, dpi = 350)
  
  ### draw a proportional bar plot
  ### pcnt < 10 -> ""
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) < 10)] <- ""
  plot_df$Pcnt <- as.character(plot_df$Pcnt)
  p <- ggplot(data=plot_df, aes_string(x="Time_Point", y="Numbers", fill="State", label="Pcnt")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Proportion of Cells") +
    xlab("Time") + ylab("Cell #") +
    geom_text(size = 5, position = position_stack(vjust = 1)) +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "CARpos_Trajectory_Inference_DOWNSAMPLED_Monocle_Bar2.png"), plot = p,
         width = 20, height = 10, dpi = 350)
  
  
  #
  ### 28. Comparison of DE genes between "GMP CAR+ S vs NS" & "After infusion CAR+ S vs NS"
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/28/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### only get the GMP persisters and non-persisters
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[union(which(Seurat_Obj$GMP_CARpos_CD8_Persister == "YES"),
                                                                                       which(Seurat_Obj$GMP_CARpos_CD8_Persister == "NO"))])
  
  ### only the patients we are interested
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  target_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                            "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                            "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### normalization
  target_Seurat_Obj <- NormalizeData(target_Seurat_Obj,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj <- FindVariableFeatures(target_Seurat_Obj,
                                            selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj <- ScaleData(target_Seurat_Obj,
                                 vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj <- RunPCA(target_Seurat_Obj,
                              features = VariableFeatures(object = target_Seurat_Obj),
                              npcs = 15)
  target_Seurat_Obj <- RunUMAP(target_Seurat_Obj, dims = 1:15)
  
  ### get CARpos-only seurat object
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$CAR)
  sub_seurat_obj <- subset(Seurat_Obj, idents = c("CARpos"))
  
  ### after gmp time points only
  after_gmp_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6",
                             "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$time2)
  sub_seurat_obj <- subset(sub_seurat_obj, idents = intersect(after_gmp_time_points,
                                                              unique(sub_seurat_obj@meta.data$time2)))
  
  ### normalization
  sub_seurat_obj <- NormalizeData(sub_seurat_obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  sub_seurat_obj <- FindVariableFeatures(sub_seurat_obj,
                                         selection.method = "vst", nfeatures = 2000)
  ### scaling
  sub_seurat_obj <- ScaleData(sub_seurat_obj,
                              vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  ### PCA
  sub_seurat_obj <- RunPCA(sub_seurat_obj,
                           features = VariableFeatures(object = sub_seurat_obj),
                           npcs = 15)
  ### UMAP
  sub_seurat_obj <- RunUMAP(sub_seurat_obj, dims = 1:15)
  
  ### get seurat object for some specific patients
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$px)
  sub_seurat_obj2 <- subset(sub_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                       "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### factorize the time2 column
  sub_seurat_obj2@meta.data$time2 <- factor(sub_seurat_obj2@meta.data$time2,
                                            levels = intersect(total_time_points, sub_seurat_obj2@meta.data$time2))
  
  
  
}
