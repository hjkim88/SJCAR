###
#   File name : Post_Scenic_Process.R
#   Author    : Hyunjin Kim
#   Date      : Feb 27, 2022
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : I ran PyScenic on the SJCAR19 data, now it's time to load those results into R,
#               and visualize the results
#
#   Instruction
#               1. Source("Post_Scenic_Process.R")
#               2. Run the function "scenic_process" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Post_Scenic_Process.R/Post_Scenic_Process.R")
#               > scenic_process(Seurat_RObj_path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds",
#                                pyscenic_result_dir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/",
#                                outputDir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/")
###

scenic_process <- function(Seurat_RObj_path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds",
                           pyscenic_result_dir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/",
                           outputDir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(reticulate, quietly = TRUE)) {
    install.packages("reticulate")
    require(reticulate, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(shadowtext, quietly = TRUE)) {
    install.packages("shadowtext")
    require(shadowtext, quietly = TRUE)
  }
  
  ### load Jeremy's object
  JCC_Seurat_Obj <- readRDS(file = Seurat_RObj_path)
  
  ### check whether the orders are the same
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
  print(identical(names(Idents(object = JCC_Seurat_Obj)), rownames(JCC_Seurat_Obj@meta.data)))
  
  ### normalization
  JCC_Seurat_Obj <- NormalizeData(JCC_Seurat_Obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### Cell cycle score (will be used later for regression out)
  JCC_Seurat_Obj <- CellCycleScoring(object = JCC_Seurat_Obj,
                                     g2m.features = cc.genes$g2m.genes,
                                     s.features = cc.genes$s.genes)
  
  ### find variable genes
  JCC_Seurat_Obj <- FindVariableFeatures(JCC_Seurat_Obj,
                                         selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  JCC_Seurat_Obj <- ScaleData(JCC_Seurat_Obj,
                              vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
  
  ### PI subsisters in the cluster3&8
  cluster3_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                   which(JCC_Seurat_Obj$AllSeuratClusters == "3")),
                                                                         which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
  cluster8_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                   which(JCC_Seurat_Obj$AllSeuratClusters == "8")),
                                                                         which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
  cluster38_pi_subsisters <- c(cluster3_pi_subsisters, cluster8_pi_subsisters)
  
  ### PI subsister clones in the cluster3&8
  cluster3_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster3_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
  cluster8_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster8_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
  cluster38_pi_subsister_clones <- unique(c(cluster3_pi_subsister_clones, cluster8_pi_subsister_clones))
  
  ### remove NA clones
  cluster3_pi_subsister_clones <- cluster3_pi_subsister_clones[which(!is.na(cluster3_pi_subsister_clones))]
  cluster8_pi_subsister_clones <- cluster8_pi_subsister_clones[which(!is.na(cluster8_pi_subsister_clones))]
  cluster38_pi_subsister_clones <- cluster38_pi_subsister_clones[which(!is.na(cluster38_pi_subsister_clones))]
  
  ### get the GMP subsister indicies that end up in PI cluster 3 & 8
  GMP_Subsisters_PI_Cluster38_idx <- intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                               which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% cluster38_pi_subsister_clones))
  
  ### add a column for the info
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 <- "Others"
  JCC_Seurat_Obj@meta.data[cluster38_pi_subsisters, "GMP_Subsisters_End_Up_In_Cluster38"] <- "PI_Subsisters_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38[GMP_Subsisters_PI_Cluster38_idx] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38[intersect(which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "Others"),
                                                              which(JCC_Seurat_Obj$GMP_CARpos_Persister == "YES"))] <- "Other_GMP_Subsisters"
  
  ### GMP subsisters end up in cluster3 and 8 vs other CD8 GMP subsisters
  JCC_Seurat_Obj@meta.data$GMP_Subsisters_End_Up_In_Cluster38_CD8 <- "Others"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_CD8[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_CD8[intersect(which(JCC_Seurat_Obj$AllSeuratClusters %in% c("1", "3", "5", "6", "7", "8", "12", "13", "16", "17", "19", "20")),
                                                                  which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "Other_GMP_Subsisters"))] <- "Other_CD8_GMP_Subsisters"
  
  ### GMP subsisters end up in cluster3 and 8 vs all other GMPs
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 <- "Others"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2[setdiff(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                              which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8"))] <- "Other_GMPs"
  
  ### GMP subsisters end up in cluster3 and 8 vs all other CD8 GMPs
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 <- "Others"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[intersect(which(JCC_Seurat_Obj$AllSeuratClusters %in% c("1", "3", "5", "6", "7", "8", "12", "13", "16", "17", "19", "20")),
                                                                    which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 == "Other_GMPs"))] <- "Other_CD8_GMPs"
  
  ### NEW FUNCTIONAL ANNOTATION BY TAY - 09/27/2021
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters <- "Others"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "1", "5", "7", "9", "10", "11", "12", "15", "19", "21"))] <- "Proliferating"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("2", "4", "6", "17"))] <- "Transitioning"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("3", "8", "14"))] <- "Functional Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("13", "20"))] <- "Dysfunctional"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("16"))] <- "Early Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("18"))] <- "Metabolically Active"
  
  
  ### a function to retrieve word between two specific characters
  get_words <- function(line, specific_chr = "'") {
    temp <- strsplit(line, split = specific_chr, fixed = TRUE)[[1]]
    temp_len <- length(temp)
    total_word_num <- floor(temp_len / 2)
    even_idx <- (1:total_word_num) * 2
    
    return(temp[even_idx])
  }
  
  ### we found that *_regulons.p files are actually an organized form of *_motifs.csv
  ### *_motifs.csv has all the TFs with their target genes,
  ### and if we organize those to unique sets, then it's *_regulons.p
  
  ### get pyscenic file list
  f <- list.files(path = pyscenic_result_dir, pattern = "*.p$")
  
  ### load regulon files and load them into R
  pd <- import("pandas")
  regulon_lists <- vector("list", length = length(f))
  names(regulon_lists) <- f
  for(file in f) {
    ### load the *.p pickle data
    pickle_data <- pd$read_pickle(paste0(pyscenic_result_dir, "/", file))
    names(pickle_data) <- sapply(pickle_data, function(x) x$transcription_factor)
    
    ### make an empty regulon list
    regulon_lists[[file]] <- vector("list", length(pickle_data))
    names(regulon_lists[[file]]) <- names(pickle_data)
    
    ### save target genes to the list
    for(tf in names(pickle_data)) {
      regulon_lists[[file]][[tf]] <- unlist(pickle_data[[tf]]$genes)
    }
    
    ### garbage collection
    gc()
  }
  
  ### load AUCELL results
  f2 <- list.files(path = pyscenic_result_dir, pattern = "*_aucell.csv$")
  aucell_results <- vector("list", length = length(f2))
  names(aucell_results) <- f2
  for(file in f2) {
    aucell_results[[file]] <- read.csv(file = paste0(pyscenic_result_dir, "/", file),
                                       row.names = 1,
                                       stringsAsFactors = FALSE, check.names = FALSE)
    
    ### cell names have additional prefix, so remove those to keep consistency with the seurat object
    if(file == "combined_aucell.csv") {
      rownames(aucell_results[[file]]) <- sapply(rownames(aucell_results[[file]]), function(x) {
        return(paste(strsplit(x, split = "_", fixed = TRUE)[[1]][2:3], collapse = "_"))
      })
    }
  }
  
  ### subset seurat objects based on the cells
  precursor_seurat <- subset(JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data)[which(rownames(JCC_Seurat_Obj@meta.data) %in% rownames(aucell_results[["precursor_aucell.csv"]]))])
  non_precursor_seurat <- subset(JCC_Seurat_Obj,
                                 cells = rownames(JCC_Seurat_Obj@meta.data)[which(rownames(JCC_Seurat_Obj@meta.data) %in% rownames(aucell_results[["non_precursor_subset_aucell.csv"]]))])
  combined_seurat <- subset(JCC_Seurat_Obj,
                            cells = rownames(JCC_Seurat_Obj@meta.data)[which(rownames(JCC_Seurat_Obj@meta.data) %in% rownames(aucell_results[["combined_aucell.csv"]]))])
  cluster3_8_seurat <- subset(JCC_Seurat_Obj,
                              cells = rownames(JCC_Seurat_Obj@meta.data)[which(rownames(JCC_Seurat_Obj@meta.data) %in% rownames(aucell_results[["cluster3_8_aucell.csv"]]))])
  
  ### save the aucell results to each of the object
  precursor_seurat[["Scenic"]] <- CreateAssayObject(counts = t(aucell_results[["precursor_aucell.csv"]]))
  non_precursor_seurat[["Scenic"]] <- CreateAssayObject(counts = t(aucell_results[["non_precursor_subset_aucell.csv"]]))
  combined_seurat[["Scenic"]] <- CreateAssayObject(counts = t(aucell_results[["combined_aucell.csv"]]))
  cluster3_8_seurat[["Scenic"]] <- CreateAssayObject(counts = t(aucell_results[["cluster3_8_aucell.csv"]]))
  
  ### run UMAP based on the Scenic result
  combined_seurat <- FindVariableFeatures(combined_seurat,
                                          selection.method = "vst",
                                          nfeatures = nrow(combined_seurat@assays$Scenic@counts),
                                          assay = "Scenic")
  combined_seurat <- ScaleData(combined_seurat,
                               assay = "Scenic",
                               do.scale = FALSE,
                               do.center = FALSE)
  combined_seurat <- RunPCA(combined_seurat,
                            features = VariableFeatures(object = combined_seurat, assay = "Scenic"),
                            npcs = 15, assay = "Scenic",
                            reduction.name = "scenic_pca")
  ElbowPlot(combined_seurat, reduction = "scenic_pca", ndims = 15)
  combined_seurat <- RunUMAP(combined_seurat, dims = 1:15,
                             assay = "Scenic",
                             reduction = "scenic_pca",
                             reduction.name = "scenic_umap",
                             reduction.key = "ScenicUMAP_")
  DefaultAssay(combined_seurat) <- "Scenic"
  combined_seurat <- FindNeighbors(combined_seurat, assay = "Scenic",
                                   reduction = "scenic_umap", dims = 1:2)
  combined_seurat <- FindClusters(combined_seurat, resolution = 0.2)
  
  ### draw Scenic UMAPs
  p <- DimPlot(object = combined_seurat, reduction = "scenic_umap",
               group.by = "New_Functional_Annotation_Based_On_Clusters",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Functional Groups") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_text(size = 48),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir, "UMAP_Scenic_Combined_Functional_Groups.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  combined_seurat$class <- combined_seurat$GMP_Subsisters_End_Up_In_Cluster38_2_CD8
  combined_seurat$class[which(combined_seurat$class == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "Precursors"
  combined_seurat$class[which(combined_seurat$class == "Other_CD8_GMPs")] <- "Non-Precursors"
  combined_seurat$class <- factor(combined_seurat$class, levels = c("Precursors", "Non-Precursors"))
  
  p <- DimPlot(object = combined_seurat, reduction = "scenic_umap",
               group.by = "class",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Effector Precursor") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_text(size = 48),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 24),
          legend.position = "top") +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir, "UMAP_Scenic_Combined_Effector_Precursor.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  p <- DimPlot(object = combined_seurat, reduction = "scenic_umap",
               group.by = "Scenic_snn_res.0.2",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_text(size = 48),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 40)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir, "UMAP_Scenic_Combined_Clusters.png"), plot = p, width = 20, height = 12, dpi = 350)
  
  ### what are the TFs thave have different regulation activity among clusters? 
  combined_seurat <- SetIdent(object = combined_seurat,
                              cells = rownames(combined_seurat@meta.data),
                              value = combined_seurat$Scenic_snn_res.0.2)
  de_result <- FindAllMarkers(combined_seurat,
                              assay = "Scenic",
                              min.pct = 0.2,
                              logfc.threshold = 0.1,
                              test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/Scenic_Differentially_Activated_TFs_Among_Clusters.xlsx"),
              sheetName = "Differentially_Activated_TFs", row.names = FALSE)
  
  
  ### just look at the differentially activated TFs between precursors vs control
  combined_seurat <- SetIdent(object = combined_seurat,
                              cells = rownames(combined_seurat@meta.data),
                              value = combined_seurat$GMP_Subsisters_End_Up_In_Cluster38_2_CD8)
  de_result <- FindMarkers(combined_seurat,
                           assay = "Scenic",
                           ident.1 = "GMP_Subsisters_End_Up_In_Cluster_3_And_8",
                           ident.2 = "Other_CD8_GMPs",
                           min.pct = 0.2,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/Scenic_Differentially_Activated_TFs_between_Precursor_vs_Control.xlsx"),
              sheetName = "Differentially_Activated_TFs", row.names = FALSE)
  
  ### draw a violin plot between precursor vs others for the differentially activated regulons
  combined_seurat$class <- combined_seurat$GMP_Subsisters_End_Up_In_Cluster38_2_CD8
  combined_seurat$class[which(combined_seurat$class == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "Precursors"
  combined_seurat$class[which(combined_seurat$class == "Other_CD8_GMPs")] <- "Non-Precursors"
  combined_seurat <- SetIdent(object = combined_seurat,
                              cells = rownames(combined_seurat@meta.data),
                              value = combined_seurat$class)
  de_result <- de_result[order(de_result$avg_log2FC),]
  p <- VlnPlot(combined_seurat, features = rownames(de_result),
               pt.size = 0)
  for(i in 1:nrow(de_result)) {
    p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
      stat_summary(fun=mean, geom="point", size=3, color="black") +
      xlab("") +
      ylab("Activity Level") +
      theme_classic(base_size = 25) +
      theme(legend.key.size = unit(3, 'cm'),
            legend.position = "none",
            legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
            legend.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35, color = "black", face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 25, color = "black", face = "bold"),
            axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
            axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
  }
  ggsave(file = paste0(outputDir, "Violin_Scenic_Differentially_Activated_TFs_between_Precursor_vs_Control.png"),
         plot = p, width = 30, height = 20, dpi = 350)
  
  
  ###
  ### cluster 3+8's turn
  ###
  
  ### run UMAP based on the Scenic result
  cluster3_8_seurat <- FindVariableFeatures(cluster3_8_seurat,
                                            selection.method = "vst",
                                            nfeatures = nrow(cluster3_8_seurat@assays$Scenic@counts),
                                            assay = "Scenic")
  cluster3_8_seurat <- ScaleData(cluster3_8_seurat,
                                 assay = "Scenic",
                                 do.scale = FALSE,
                                 do.center = FALSE)
  cluster3_8_seurat <- RunPCA(cluster3_8_seurat,
                              features = VariableFeatures(object = cluster3_8_seurat, assay = "Scenic"),
                              npcs = 15, assay = "Scenic",
                              reduction.name = "scenic_pca")
  ElbowPlot(cluster3_8_seurat, reduction = "scenic_pca", ndims = 15)
  cluster3_8_seurat <- RunUMAP(cluster3_8_seurat, dims = 1:15,
                               assay = "Scenic",
                               reduction = "scenic_pca",
                               reduction.name = "scenic_umap",
                               reduction.key = "ScenicUMAP_")
  DefaultAssay(cluster3_8_seurat) <- "Scenic"
  cluster3_8_seurat <- FindNeighbors(cluster3_8_seurat, assay = "Scenic",
                                     reduction = "scenic_umap", dims = 1:2)
  cluster3_8_seurat <- FindClusters(cluster3_8_seurat, resolution = 0.2)
  
  
  ### draw Scenic UMAPs
  p <- DimPlot(object = cluster3_8_seurat, reduction = "scenic_umap",
               group.by = "AllSeuratClusters",
               pt.size = 3) +
    ggtitle("") +
    labs(color="Original Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_text(size = 48),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir, "UMAP_Scenic_Cluster3_8_AllSeuratClusters.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  p <- DimPlot(object = cluster3_8_seurat, reduction = "scenic_umap",
               group.by = "GMP",
               pt.size = 3) +
    ggtitle("") +
    labs(color="GMP/PI") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_text(size = 48),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir, "UMAP_Scenic_Cluster3_8_GMP.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  p <- DimPlot(object = cluster3_8_seurat, reduction = "scenic_umap",
               group.by = "New_Functional_Annotation_Based_On_Clusters",
               pt.size = 3) +
    ggtitle("") +
    labs(color="Functional Groups") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_text(size = 48),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 48),
          legend.text = element_text(size = 48)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir, "UMAP_Scenic_Cluster3_8_New_Functional_Annotation_Based_On_Clusters.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ###
  ### get the PI cells only
  ###
  cluster3_8_PI_seurat <- subset(cluster3_8_seurat,
                                 cells = rownames(cluster3_8_seurat@meta.data)[which(cluster3_8_seurat$GMP == "PI")])
  
  ### run UMAP based on the Scenic result
  cluster3_8_PI_seurat <- FindVariableFeatures(cluster3_8_PI_seurat,
                                               selection.method = "vst",
                                               nfeatures = nrow(cluster3_8_PI_seurat@assays$Scenic@counts),
                                               assay = "Scenic")
  cluster3_8_PI_seurat <- ScaleData(cluster3_8_PI_seurat,
                                    assay = "Scenic",
                                    do.scale = FALSE,
                                    do.center = FALSE)
  cluster3_8_PI_seurat <- RunPCA(cluster3_8_PI_seurat,
                                 features = VariableFeatures(object = cluster3_8_PI_seurat, assay = "Scenic"),
                                 npcs = 15, assay = "Scenic",
                                 reduction.name = "scenic_pca")
  ElbowPlot(cluster3_8_PI_seurat, reduction = "scenic_pca", ndims = 15)
  cluster3_8_PI_seurat <- RunUMAP(cluster3_8_PI_seurat, dims = 1:15,
                                  assay = "Scenic",
                                  reduction = "scenic_pca",
                                  reduction.name = "scenic_umap",
                                  reduction.key = "ScenicUMAP_")
  DefaultAssay(cluster3_8_PI_seurat) <- "Scenic"
  cluster3_8_PI_seurat <- FindNeighbors(cluster3_8_PI_seurat, assay = "Scenic",
                                        reduction = "scenic_umap", dims = 1:2)
  cluster3_8_PI_seurat <- FindClusters(cluster3_8_PI_seurat, resolution = 0.2)
  
  ### draw Scenic UMAPs
  p <- DimPlot(object = cluster3_8_PI_seurat, reduction = "scenic_umap",
               group.by = "AllSeuratClusters", label = TRUE,
               pt.size = 3) +
    ggtitle("") +
    labs(color="Original Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_text(size = 48),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 24)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p <- p + geom_shadowtext(data = p$layers[[2]]$data, aes(x = ScenicUMAP_1, y = ScenicUMAP_2, label=AllSeuratClusters),
                           size=6, color="cornsilk2", bg.color="black", bg.r=0.2)
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir, "UMAP_Scenic_Cluster3_8_PI_AllSeuratClusters.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  p <- DimPlot(object = cluster3_8_PI_seurat, reduction = "scenic_umap",
               group.by = "AllSeuratClusters", label = TRUE,
               pt.size = 3,
               cols = c("3" = "red", "8" = "blue")) +
    ggtitle("") +
    labs(color="Original Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_text(size = 48),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 24)) +
    guides(colour = guide_legend(override.aes = list(size=10)))
  p <- p + geom_shadowtext(data = p$layers[[2]]$data, aes(x = ScenicUMAP_1, y = ScenicUMAP_2, label=AllSeuratClusters),
                           size=6, color="cornsilk2", bg.color="black", bg.r=0.2)
  p[[1]]$layers[[1]]$aes_params$alpha <- 1
  ggsave(paste0(outputDir, "UMAP_Scenic_Cluster3_8_PI_AllSeuratClusters2.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### Differentially activated TFs between cluster 3 vs cluster 8 in PI
  cluster3_8_PI_seurat <- SetIdent(object = cluster3_8_PI_seurat,
                                   cells = rownames(cluster3_8_PI_seurat@meta.data),
                                   value = cluster3_8_PI_seurat$AllSeuratClusters)
  de_result <- FindMarkers(cluster3_8_PI_seurat,
                           assay = "Scenic",
                           ident.1 = "3",
                           ident.2 = "8",
                           min.pct = 0.1,
                           logfc.threshold = 0,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/Scenic_Differentially_Activated_TFs_between_Cluster3_vs_Cluster8.xlsx"),
              sheetName = "Differentially_Activated_TFs", row.names = FALSE)
  
  
}
