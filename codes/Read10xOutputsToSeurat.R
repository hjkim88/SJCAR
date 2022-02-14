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
#               2. Run the function "read10x_and_make_seuratobj" - specify the input directory and the output path
#               3. The results will be generated under the output path
#
#   Example
#               > source("The_directory_of_Read10xOutputsToSeurat.R/Read10xOutputsToSeurat.R")
#               > read10x_and_make_seuratobj(ten_x_dir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/JCC/JCC212_SJCAR19/JCC212_SJCAR19_C1AggregOct2020ns/filtered_feature_bc_matrix/",
#                                            outputPath="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/SJCAR19_Oct2020_Seurat_Obj5.RDS")
###

read10x_and_make_seuratobj <- function(ten_x_dir="C:/Users/hkim8/SJ/SJCAR19/JCC212_SJCAR19_C1AggregOct2020ns/filtered_feature_bc_matrix/",
                                       outputPath="./data/SJCAR19_Oct2020_Seurat_Obj5.RDS") {
  
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
  if(!require(SeuratDisk, quietly = TRUE)) {
    remotes::install_github("mojaveazure/seurat-disk")
    require(SeuratDisk, quietly = TRUE)
  }
  if(!require(reticulate, quietly = TRUE)) {
    install.packages("reticulate")
    require(reticulate, quietly = TRUE)
  }
  if(!require(zellkonverter, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("zellkonverter")
    require(zellkonverter, quietly = TRUE)
  }
  if(!require(SummarizedExperiment, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("SummarizedExperiment")
    require(SummarizedExperiment, quietly = TRUE)
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
  
  ### if the matrix is too big,
  # Convert(source = "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/filtered_feature_bc_matrix/mat.h5ad", dest = "h5seurat")
  # SJCAR19_Oct2020_Seurat_Obj <- LoadH5Seurat("Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/filtered_feature_bc_matrix/mat.h5seurat")
  # ad <- import("anndata", convert = FALSE)
  # SJCAR19_Oct2020_Seurat_Obj <- ad$read_h5ad("Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/filtered_feature_bc_matrix/mat.h5ad")
  # SJCAR19_Oct2020_Seurat_Obj <- Convert(SJCAR19_Oct2020_Seurat_Obj, to = "seurat")
  # a <- SJCAR19_Oct2020_Seurat_Obj$X
  
  
  # ### HEALTHY DONOR ONLY - TEMPORARY
  # sce <- readH5AD(file = "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/filtered_feature_bc_matrix/mat_healthy_donors.h5ad")
  # SJCAR19_Oct2020_healthy_Seurat_Obj <- CreateSeuratObject(counts = assay(sce, "X"),
  #                                                          project = "SJCAR19",
  #                                                          min.cells = 3,
  #                                                          min.features = 200)
  # 
  # ### dissect the combined barcodes
  # SJCAR19_Oct2020_healthy_Seurat_Obj$orig.barcode <- sapply(rownames(SJCAR19_Oct2020_healthy_Seurat_Obj@meta.data), function(x) {
  #   return(strsplit(x, split = "-", fixed = TRUE)[[1]][1])
  # })
  # SJCAR19_Oct2020_healthy_Seurat_Obj$lib.ident <- sapply(rownames(SJCAR19_Oct2020_healthy_Seurat_Obj@meta.data), function(x) {
  #   return(strsplit(x, split = "-", fixed = TRUE)[[1]][2])
  # })
  # 
  # ### jacard index
  # cal_jaccard <- function(x, y) {
  #   return(length(intersect(x, y)) / length(union(x, y)))
  # }
  # 
  # ### annotate the cells which are which
  # ### the barcodes don't match exactly because when creating seurat, there was a filtering process
  # bar_path <- "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/"
  # bar_f <- list.files(path = bar_path, pattern = "*barcodes.tsv.gz$")
  # SJCAR19_Oct2020_healthy_Seurat_Obj$library <- NA
  # for(sample_name in bar_f) {
  #   
  #   temp_barcodes <- as.character(read.csv(gzfile(paste0(bar_path, sample_name)),
  #                                          header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)[,1])
  #   temp_barcodes <- sapply(temp_barcodes, function(x) {
  #     return(strsplit(x, split = "-", fixed = TRUE)[[1]][1])
  #   }, USE.NAMES = FALSE)
  #   
  #   jaccard_list <- rep(0, length(unique(SJCAR19_Oct2020_healthy_Seurat_Obj$lib.ident)))
  #   names(jaccard_list) <- unique(SJCAR19_Oct2020_healthy_Seurat_Obj$lib.ident)
  #   for(lib in unique(SJCAR19_Oct2020_healthy_Seurat_Obj$lib.ident)) {
  #     lib_barcodes <- SJCAR19_Oct2020_healthy_Seurat_Obj$orig.barcode[which(SJCAR19_Oct2020_healthy_Seurat_Obj$lib.ident == lib)]
  #     
  #     jaccard_list[lib] <- cal_jaccard(temp_barcodes, lib_barcodes)
  #   }
  #   
  #   if(length(which(jaccard_list > 0.8)) > 1) {
  #     writeLines(paste("ERROR: There are multiple samples matched with this barcodes:", sample_name))
  #   }
  #   
  #   SJCAR19_Oct2020_healthy_Seurat_Obj$library[which(SJCAR19_Oct2020_healthy_Seurat_Obj$lib.ident == names(jaccard_list)[which(jaccard_list > 0.8)])] <- paste(strsplit(sample_name, split = "_", fixed = TRUE)[[1]][3:4], collapse = "_")
  #   
  # }
  # 
  # ### set Other columns
  # SJCAR19_Oct2020_healthy_Seurat_Obj$px <- ""
  # SJCAR19_Oct2020_healthy_Seurat_Obj$px[which(SJCAR19_Oct2020_healthy_Seurat_Obj$library %in% c("Donor33_PreTransB", "GMPdonor33_short"))] <- "Donor33"
  # SJCAR19_Oct2020_healthy_Seurat_Obj$px[which(SJCAR19_Oct2020_healthy_Seurat_Obj$library %in% c("Donor30_PreTransB", "GMPdonor30_short"))] <- "Donor30"
  # SJCAR19_Oct2020_healthy_Seurat_Obj$px[which(SJCAR19_Oct2020_healthy_Seurat_Obj$library %in% c("Donor32_PreTransB", "GMPdonor32_short"))] <- "Donor32"
  # 
  # SJCAR19_Oct2020_healthy_Seurat_Obj$time <- ""
  # SJCAR19_Oct2020_healthy_Seurat_Obj$time[which(SJCAR19_Oct2020_healthy_Seurat_Obj$library %in% c("Donor33_PreTransB", "Donor30_PreTransB", "Donor32_PreTransB"))] <- "PreTransB"
  # SJCAR19_Oct2020_healthy_Seurat_Obj$time[which(SJCAR19_Oct2020_healthy_Seurat_Obj$library %in% c("GMPdonor30_short", "GMPdonor33_short", "GMPdonor32_short"))] <- "GMP"
  # 
  # ### CD4/CD8
  # FeaturePlot(object = SJCAR19_Oct2020_healthy_Seurat_Obj, features = "CD4")
  # FeaturePlot(object = SJCAR19_Oct2020_healthy_Seurat_Obj, features = "CD8A")
  # FeaturePlot(object = SJCAR19_Oct2020_healthy_Seurat_Obj, features = "CD8B")
  # SJCAR19_Oct2020_healthy_Seurat_Obj$CD4_CD8_by_Exp <- NA
  # SJCAR19_Oct2020_healthy_Seurat_Obj$CD4_CD8_by_Exp[which(SJCAR19_Oct2020_healthy_Seurat_Obj@assays$RNA@counts["CD4",] > SJCAR19_Oct2020_healthy_Seurat_Obj@assays$RNA@counts["CD8A",])] <- "CD4"
  # SJCAR19_Oct2020_healthy_Seurat_Obj$CD4_CD8_by_Exp[which(SJCAR19_Oct2020_healthy_Seurat_Obj@assays$RNA@counts["CD8A",] > SJCAR19_Oct2020_healthy_Seurat_Obj@assays$RNA@counts["CD4",])] <- "CD8"
  # DimPlot(object = SJCAR19_Oct2020_healthy_Seurat_Obj, reduction = "umap", raster = TRUE,
  #         group.by = "CD4_CD8_by_Exp")
  # 
  # ### CAR
  # SJCAR19_Oct2020_healthy_Seurat_Obj@meta.data$CAR <- sapply(SJCAR19_Oct2020_healthy_Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",rownames(SJCAR19_Oct2020_healthy_Seurat_Obj@meta.data)],
  #                                    function(x) {
  #                                      if(x > 0) {
  #                                        return("CARpos")
  #                                      } else {
  #                                        return("CARneg")
  #                                      }
  #                                    })
  # 
  # ### to follow the other procedure in the below,
  # SJCAR19_Oct2020_Seurat_Obj <- SJCAR19_Oct2020_healthy_Seurat_Obj
  
  
  ### MT percentage
  SJCAR19_Oct2020_Seurat_Obj[["percent.mt"]] <- PercentageFeatureSet(SJCAR19_Oct2020_Seurat_Obj, pattern = "^MT-")
  
  ### Visualize QC metrics as a violin plot
  VlnPlot(SJCAR19_Oct2020_Seurat_Obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot(density(SJCAR19_Oct2020_Seurat_Obj@meta.data$nFeature_RNA))
  plot(density(SJCAR19_Oct2020_Seurat_Obj@meta.data$percent.mt))
  
  ### Filter cells that have unique feature counts over 5000 or less than 300
  ### Filter cells that have > 10% mitochondrial counts
  SJCAR19_Oct2020_Seurat_Obj <- subset(SJCAR19_Oct2020_Seurat_Obj, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)
  
  ### normalization
  SJCAR19_Oct2020_Seurat_Obj <- NormalizeData(SJCAR19_Oct2020_Seurat_Obj,
                                              normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### Cell cycle score (will be used later for regression out)
  SJCAR19_Oct2020_Seurat_Obj <- CellCycleScoring(object = SJCAR19_Oct2020_Seurat_Obj,
                                                 g2m.features = cc.genes$g2m.genes,
                                                 s.features = cc.genes$s.genes)
  
  ### find variable genes
  SJCAR19_Oct2020_Seurat_Obj <- FindVariableFeatures(SJCAR19_Oct2020_Seurat_Obj,
                                                     selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  SJCAR19_Oct2020_Seurat_Obj <- ScaleData(SJCAR19_Oct2020_Seurat_Obj,
                                          vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### PCA
  SJCAR19_Oct2020_Seurat_Obj <- RunPCA(SJCAR19_Oct2020_Seurat_Obj,
                                       features = VariableFeatures(object = SJCAR19_Oct2020_Seurat_Obj))
  
  ### UMAP
  SJCAR19_Oct2020_Seurat_Obj <- RunUMAP(SJCAR19_Oct2020_Seurat_Obj, dims = 1:15)
  
  ### Elbow plot
  ElbowPlot(SJCAR19_Oct2020_Seurat_Obj, ndims = 15, reduction = "pca")
  
  ### clustering
  SJCAR19_Oct2020_Seurat_Obj <- FindNeighbors(SJCAR19_Oct2020_Seurat_Obj, dims = 1:15)
  SJCAR19_Oct2020_Seurat_Obj <- FindClusters(SJCAR19_Oct2020_Seurat_Obj, resolution = 0.5)
  
  ### draw a umap
  DimPlot(object = SJCAR19_Oct2020_Seurat_Obj, reduction = "umap", raster = TRUE,
          group.by = "seurat_clusters")
  
  ### save the Seurat object as RDS file
  saveRDS(SJCAR19_Oct2020_Seurat_Obj, file = outputPath)
  
}
