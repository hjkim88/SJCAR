###
#   File name : GMP_Persistency_Marker_Discovery.R
#   Author    : Hyunjin Kim
#   Date      : Dec 29, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. What markers predict whether a lineage (CAR+/CAR-) will persist after infusion?
#               - DE analysis -> Consider patient effect in the design matrix
#               - GO/Pathway analysis
#               - Comparison of various factors between persisters & non-persisters (violin plot?)
#               - Using the DE genes to make a classifier
#               - To ensure an individual does not drive the found pattern?
#                 a) if the classifier (LOOCV) worked well, it is a proof
#                 b) One patient wouldn't be very powerful, because I will consider the patient info
#                    as a factor in the design matrix (with a toy example to show)
#               - Revered PCA/UMAP (gene-patient conversion) to show whether those genes are separated by persistency
#
#   Instruction
#               1. Source("GMP_Persistency_Marker_Discovery.R")
#               2. Run the function "persistency_study" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_GMP_Persistency_Marker_Discovery.R/GMP_Persistency_Marker_Discovery.R")
#               > persistency_study(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                   clonotype_lineage_info_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Clonotype_Lineages.RDS",
#                                   outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Persistency/")
###

persistency_study <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
                              clonotype_lineage_info_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Clonotype_Lineages.RDS",
                              outputDir="./results/New2/Persistency/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  options(java.parameters = "-Xmx10240m")
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
  if(!require(DESeq2, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("DESeq2")
    require(DESeq2, quietly = TRUE)
  }
  if(!require(zinbwave, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("zinbwave")
    require(zinbwave, quietly = TRUE)
  }
  if(!require(BiocParallel, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("BiocParallel")
    require(BiocParallel, quietly = TRUE)
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
  
  ### load clonotype lineages info
  SJCAR19_Clonotype_Frequency <- readRDS(clonotype_lineage_info_path)
  
  ### set idents with the libary
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$GMP_CARpos_Persister)
  
  ### get GMP CAR+ cells as a SingleCellExperiment object
  GMP_CARpos_sce <- as.SingleCellExperiment(subset(Seurat_Obj, idents = c("YES", "NO")), assay = "RNA")
  
  ### for memory efficiency, remove the large Seurat object
  rm(Seurat_Obj)
  gc()
  
  ### QC
  ### 26027 genes & 109614 cells
  row_sum <- rowSums(counts(GMP_CARpos_sce) > 0)
  col_sum <- colSums(counts(GMP_CARpos_sce) > 0)
  plot(density(row_sum))
  plot(density(col_sum))
  
  ### set parameter
  ### gene.min.pct: in a gene, at least n % of each condition (cells) should be expressed
  ### cell.min.num: in a cell, at least n genes should be expressed (determine after looking at the density plot above)
  gene.min.pct <- 0.2
  cell.min.num <- 2000
  cell.max.num <- 5000
  
  ### low count filter for genes
  yes_num <- length(which(GMP_CARpos_sce$ident == "YES"))
  no_num <- length(which(GMP_CARpos_sce$ident == "NO"))
  yes_sum <- rowSums(counts(GMP_CARpos_sce)[,which(GMP_CARpos_sce$ident == "YES")] > 0)
  no_sum <- rowSums(counts(GMP_CARpos_sce)[,which(GMP_CARpos_sce$ident == "NO")] > 0)
  yes_exp_pct <- yes_sum / yes_num
  no_exp_pct <- no_sum / no_num
  keep <- intersect(which(yes_exp_pct >= gene.min.pct),
                    which(no_exp_pct >= gene.min.pct))
  GMP_CARpos_sce <- GMP_CARpos_sce[keep,]
  
  ### low count filter for cells
  keep <- intersect(which(colSums(counts(GMP_CARpos_sce) > 0) >= cell.min.num),
                    which(colSums(counts(GMP_CARpos_sce) > 0) <= cell.max.num))
  GMP_CARpos_sce <- GMP_CARpos_sce[,keep]
  
  
  ### there are too many cells (especially in the second condition), so we will randomly choose some
  ### to lower the computational complexity
  set.seed(1234)
  keep <- sample(which(GMP_CARpos_sce$ident == "NO"), length(which(GMP_CARpos_sce$ident == "YES")))
  keep <- union(keep, which(GMP_CARpos_sce$ident == "YES"))
  GMP_CARpos_sce <- GMP_CARpos_sce[,keep]
  
  ### "counts" should be in the first place in the assayNames(GMP_CARpos_sce)
  nms <- c("counts", setdiff(assayNames(GMP_CARpos_sce), "counts"))
  assays(GMP_CARpos_sce) <- assays(GMP_CARpos_sce)[nms]
  
  ### epsilon setting as recommended by the ZINB-WaVE integration paper
  system.time({
    zinb <- zinbwave(GMP_CARpos_sce, K=0, observationalWeights=TRUE,
                     BPPARAM=SerialParam(), epsilon=1e12)
  })
  
  
  dds <- DESeqDataSet(zinb, design=~condition)
  dds <- estimateSizeFactors(dds, type="poscounts")
  library(scran)
  scr <- computeSumFactors(dds)
  dat <- data.frame(true=dds$ExpLibSize,
                    pos=sizeFactors(dds),
                    sum=sizeFactors(scr))
  dat$true <- dat$true / exp(mean(log(dat$true)))
  panel.scatter <- function(x,y,...) {
    points(x,y,...)
    abline(0,1,col="red",lwd=2)
    legend("topleft", legend=round(cor(x,y),3))
  }
  pairs(dat, panel=panel.scatter)
  
  
  
  
}
