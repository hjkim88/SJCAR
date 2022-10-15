###
#   File name : Additional_Analyses_Apr.R
#   Author    : Hyunjin Kim
#   Date      : Apr 9, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. Is the clone size naturally bigger in persisters when compared to non-persisters?
#               2. How do gene expressions of persisters change over time?
#               3. Also do the analyses for persisters vs CAR- / CAR+ vs CAR-
#               4. Use alpha chain only to redo the analyses
#               5. TCRdist3 to find TCR clusters that are associated to persisters
#
#   Instruction
#               1. Source("Additional_Analyses_Apr.R")
#               2. Run the function "additional_apr" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Additional_Analyses_Apr.R/Additional_Analyses_Apr.R")
#               > additional_apr(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Additional_Apr/")
###

additional_apr <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj_Total2.RDS",
                           outputDir="./results/New3/Additional_Apr/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
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
  gmp_after_time_points <- c("GMP", "GMP-redo", "Wk1", "Wk1b", "Wk2", "Wk3",
                             "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo")
  
  ### 1. Is the clone size naturally bigger in persisters when compared to background?
  ###    And how about persister vs non-persisters in GMP?
  ### with violin plot
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  non_persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "NO")])
  
  ### get gmp indicies
  gmp_indicies <- union(which(Seurat_Obj@meta.data$time == "GMP"),
                        which(Seurat_Obj@meta.data$time == "GMP-redo"))
  
  ### get persister vs non-persister clone sizes in GMP
  persister_clone_sizes_gmp <- rep(0, length(persister_clones))
  names(persister_clone_sizes_gmp) <- persister_clones
  for(clone in persister_clones) {
    persister_clone_sizes_gmp[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[gmp_indicies] == clone))
  }
  
  non_persister_clone_sizes_gmp <- rep(0, length(non_persister_clones))
  names(non_persister_clone_sizes_gmp) <- non_persister_clones
  for(clone in non_persister_clones) {
    non_persister_clone_sizes_gmp[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[gmp_indicies] == clone))
  }
  
  ### prepare a data frame for the plot
  plot_df <- data.frame(Clone_Size=c(persister_clone_sizes_gmp,
                                     non_persister_clone_sizes_gmp),
                        Group=c(rep("Persisters", length(persister_clone_sizes_gmp)),
                                rep("Non-Persisters", length(non_persister_clone_sizes_gmp))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df$Group <- factor(plot_df$Group, levels = c("Persisters", "Non-Persisters"))
  
  ### draw a violin plot
  ggplot(plot_df, aes_string(x="Group", y="Clone_Size", fill = "Group")) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill="white") +
    ylim(c(-0.5, 2)) +
    geom_text(data = data.frame(Group=c("Persisters", "Non-Persisters"),
                                Median=c(median(persister_clone_sizes_gmp),
                                median(non_persister_clone_sizes_gmp)),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Median"),
              size = 10, hjust = -2, vjust = -2) +
    geom_text(data = data.frame(Group=c("Persisters", "Non-Persisters"),
                                Median=c(median(persister_clone_sizes_gmp),
                                         median(non_persister_clone_sizes_gmp)),
                                Length=c(paste("n =", length(persister_clone_sizes_gmp)),
                                         paste("n =", length(non_persister_clone_sizes_gmp))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Length"),
              size = 5, hjust = 0.5, vjust = 35) +
    labs(title="", x="", y = "Clone Size", fill = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 30))
  ggsave(paste0(outputDir, "/", "Violin_GMP_Clone_Size_P_vs_NP.png"), width = 20, height = 12, dpi = 500)
  
  ### density plot
  ggplot(plot_df, aes_string(x="Clone_Size", col="Group")) +
    geom_density(size = 2) +
    xlim(c(0,20)) +
    geom_text(data = data.frame(Group=c("Persisters", "Non-Persisters"),
                                x = c(18, 18),
                                y = c(0.4, 0.2),
                                Length=c(paste("n =", length(persister_clone_sizes_gmp)),
                                         paste("n =", length(non_persister_clone_sizes_gmp))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "x", y = "y", label = "Length"),
              size = 5, hjust = 0.5, vjust = 0.5) +
    labs(title="", x="Clone Size", y = "Density", col = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_classic(base_size = 36) +
    theme(legend.text = element_text(size = 35))
  ggsave(paste0(outputDir, "/", "Density_GMP_Clone_Size_P_vs_NP.png"), width = 15, height = 12, dpi = 500)
  
  
  #
  ### AFTER-INFUSION - PERSISTER VS BACKGROUND
  #
  
  ### get AI indicies
  ai_indicies <- which(Seurat_Obj@meta.data$time %in% c("Wk1", "Wk1b", "Wk2", "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo"))
  
  ### target lineages
  persister_clones_ai <- unique(intersect(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")],
                                          Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[ai_indicies]))
  background_clones_ai <- unique(intersect(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(is.na(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister))],
                                           Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[ai_indicies]))
  
  ### get persister vs background clone sizes in AI time points
  persister_clone_sizes_ai <- rep(0, length(persister_clones_ai))
  names(persister_clone_sizes_ai) <- persister_clones_ai
  for(clone in persister_clones_ai) {
    persister_clone_sizes_ai[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[ai_indicies] == clone))
  }
  
  background_clone_sizes_ai <- rep(0, length(background_clones_ai))
  names(background_clone_sizes_ai) <- background_clones_ai
  for(clone in background_clones_ai) {
    background_clone_sizes_ai[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[ai_indicies] == clone))
  }
  
  ### prepare a data frame for the plot
  plot_df <- data.frame(Clone_Size=c(persister_clone_sizes_ai,
                                     background_clone_sizes_ai),
                        Group=c(rep("Persisters_After_Infusion", length(persister_clone_sizes_ai)),
                                rep("Background_After_Infusion", length(background_clone_sizes_ai))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df$Group <- factor(plot_df$Group, levels = c("Persisters_After_Infusion", "Background_After_Infusion"))
  
  ### draw a violin plot
  ggplot(plot_df, aes_string(x="Group", y="Clone_Size", fill = "Group")) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill="white") +
    ylim(c(-0.5, 2)) +
    geom_text(data = data.frame(Group=c("Persisters_After_Infusion", "Background_After_Infusion"),
                                Median=c(median(persister_clone_sizes_ai),
                                         median(background_clone_sizes_ai)),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Median"),
              size = 10, hjust = -2, vjust = -2) +
    geom_text(data = data.frame(Group=c("Persisters_After_Infusion", "Background_After_Infusion"),
                                Median=c(median(persister_clone_sizes_ai),
                                         median(background_clone_sizes_ai)),
                                Length=c(paste("n =", length(persister_clone_sizes_ai)),
                                         paste("n =", length(background_clone_sizes_ai))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Length"),
              size = 5, hjust = 0.5, vjust = 35) +
    labs(title="", x="", y = "Clone Size", fill = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 30))
  ggsave(paste0(outputDir, "/", "Violin_After_Infusion_Clone_Size_P_vs_B.png"), width = 20, height = 12, dpi = 500)
  
  ### density plot
  ggplot(plot_df, aes_string(x="Clone_Size", col="Group")) +
    geom_density(size = 2) +
    xlim(c(0,20)) +
    geom_text(data = data.frame(Group=c("Persisters_After_Infusion", "Background_After_Infusion"),
                                x = c(18, 18),
                                y = c(0.3, 0.2),
                                Length=c(paste("n =", length(persister_clone_sizes_ai)),
                                         paste("n =", length(background_clone_sizes_ai))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "x", y = "y", label = "Length"),
              size = 5, hjust = 0.5, vjust = 0.5) +
    labs(title="", x="Clone Size", y = "Density", col = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_classic(base_size = 36) +
    theme(legend.text = element_text(size = 35))
  ggsave(paste0(outputDir, "/", "Density_After_Infusion_Clone_Size_P_vs_B.png"), width = 15, height = 12, dpi = 500)
  
  
  #
  ### how the DE gene expression from persisters and non-persisters are changing over time
  #
  
  ### get some DE genes
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister)
  de_result <- FindMarkers(Seurat_Obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### the number of DE genes to be looked into
  n <- 10
  de_result <- de_result[order(de_result$p_val_adj),]
  de_genes <- rownames(de_result)[1:n]
  
  ### persister cell indicies
  persister_clones_ai_idx <- which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones_ai)
  
  ### background cell indicies
  background_clone_sizes_ai_idx <- which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta %in% background_clones_ai)
  
  ### for each gene make a plot
  for(gene in de_genes) {
    
    ### prepare a data frame for the plot
    medians <- data.frame(Median=rep(0, length(gmp_after_time_points)*2),
                          Group=rep(c("Persisters", "Background"), length(gmp_after_time_points)),
                          Time=as.vector(sapply(gmp_after_time_points, function(x) rep(x, 2))),
                          stringsAsFactors = FALSE, check.names = FALSE)
    plot_df <- data.frame(Gex=c(Seurat_Obj@assays$RNA@data[gene,intersect(which(Seurat_Obj@meta.data$time == gmp_after_time_points[1]),
                                                                          persister_clones_ai_idx)],
                                Seurat_Obj@assays$RNA@data[gene,intersect(which(Seurat_Obj@meta.data$time == gmp_after_time_points[1]),
                                                                          background_clone_sizes_ai_idx)]),
                          Group=c(rep("Persisters", length(intersect(which(Seurat_Obj@meta.data$time == gmp_after_time_points[1]),
                                                                    persister_clones_ai_idx))),
                                 rep("Background", length(intersect(which(Seurat_Obj@meta.data$time == gmp_after_time_points[1]),
                                                                    background_clone_sizes_ai_idx)))),
                          Time=rep(gmp_after_time_points[1], length(intersect(which(Seurat_Obj@meta.data$time == gmp_after_time_points[1]),
                                                                          persister_clones_ai_idx)) + length(intersect(which(Seurat_Obj@meta.data$time == gmp_after_time_points[1]),
                                                                                                                       background_clone_sizes_ai_idx))),
                          stringsAsFactors = FALSE, check.names = FALSE)
    cnt <- 1
    medians$Median[cnt] <- median(Seurat_Obj@assays$RNA@data[gene,intersect(which(Seurat_Obj@meta.data$time == gmp_after_time_points[1]),
                                                                            persister_clones_ai_idx)])
    cnt <- cnt + 1
    medians$Median[cnt] <- median(Seurat_Obj@assays$RNA@data[gene,intersect(which(Seurat_Obj@meta.data$time == gmp_after_time_points[1]),
                                                                            background_clone_sizes_ai_idx)])
    cnt <- cnt + 1
    for(tp in gmp_after_time_points[-1]) {
      plot_df <- rbind(plot_df,
                       data.frame(Gex=c(Seurat_Obj@assays$RNA@data[gene,intersect(which(Seurat_Obj@meta.data$time == tp),
                                                                                  persister_clones_ai_idx)],
                                        Seurat_Obj@assays$RNA@data[gene,intersect(which(Seurat_Obj@meta.data$time == tp),
                                                                                  background_clone_sizes_ai_idx)]),
                                  Group=c(rep("Persisters", length(intersect(which(Seurat_Obj@meta.data$time == tp),
                                                                             persister_clones_ai_idx))),
                                          rep("Background", length(intersect(which(Seurat_Obj@meta.data$time == tp),
                                                                             background_clone_sizes_ai_idx)))),
                                  Time=rep(tp, length(intersect(which(Seurat_Obj@meta.data$time == tp),
                                                                persister_clones_ai_idx)) + length(intersect(which(Seurat_Obj@meta.data$time == tp),
                                                                                                             background_clone_sizes_ai_idx))),
                                  stringsAsFactors = FALSE, check.names = FALSE))
      medians$Median[cnt] <- median(Seurat_Obj@assays$RNA@data[gene,intersect(which(Seurat_Obj@meta.data$time == tp),
                                                                              persister_clones_ai_idx)])
      cnt <- cnt + 1
      medians$Median[cnt] <- median(Seurat_Obj@assays$RNA@data[gene,intersect(which(Seurat_Obj@meta.data$time == tp),
                                                                              background_clone_sizes_ai_idx)])
      cnt <- cnt + 1
    }
    plot_df$Group <- factor(plot_df$Group, levels = c("Persisters", "Background"))
    plot_df$Time <- factor(plot_df$Time, levels = gmp_after_time_points)
    
    ### draw a violin plot
    ggplot(plot_df, aes_string(x="Time", y="Gex", fill = "Group")) +
      geom_violin(trim=FALSE) +
      geom_text(data = data.frame(Group=c("Persisters_After_Infusion", "Background_After_Infusion"),
                                  Median=c(median(persister_clone_sizes_ai),
                                           median(background_clone_sizes_ai)),
                                  stringsAsFactors = FALSE, check.names = FALSE),
                aes_string(x = "Group", y = "Median", label = "Median"),
                size = 10, hjust = -2, vjust = -2) +
      labs(title="", x="", y = "Normalized Gene Expression", fill = "") +
      scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
      theme_classic(base_size = 36) +
      theme(axis.text.x = element_text(size = 30),
            axis.title.y = element_text(size = 30),
            legend.text = element_text(size = 30))
    ggsave(paste0(outputDir, "/", "Violin_Gex_P_vs_B_", gene, ".png"), width = 36, height = 12, dpi = 350, limitsize = FALSE)
    
  }
  
}
