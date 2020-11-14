###
#   File name : IFITM3_PI3K_Pathway_Investigation.R
#   Author    : Hyunjin Kim
#   Date      : Nov 12, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Investigate associations of IFITM3 expression and PI3K signaling pathway.
#               1. (Combined) Correlation plots of IFITM3 - PI3K signaling genes based on expression
#               2. Heatmap of the IFITM3 & PI3K signaling genes with various annotation labels
#
#   Instruction
#               1. Source("IFITM3_PI3K_Pathway_Investigation.R")
#               2. Run the function "ifitm3_pi3k" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_IFITM3_PI3K_Pathway_Investigation.R/IFITM3_PI3K_Pathway_Investigation.R")
#               > ifitm3_pi3k(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                             clonotype_lineage_info_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/NEW_SJCAR_SEURAT_OBJ//SJCAR19_Clonotype_Lineages.RDS",
#                             outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results2/")
###

ifitm3_pi3k <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
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
  if(!require(msigdbr, quietly = TRUE)) {
    install.packages("msigdbr")
    library(msigdbr, quietly = TRUE)
  }
  
  ### new output directory
  outputDir2 <- paste0(outputDir, "IFITM3_PI3K/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
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
  
  ### db preparation
  # MSIGDB - human
  m_df <- msigdbr(species = "Homo sapiens") 
  m_list <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  ### PI3K associated pathways
  pi3k_pathways <- m_list[names(m_list)[grep("PI3K", names(m_list))][c(9:13, 16:17, 21:23)]]
  
  
  ### load clonotype lineages info
  SJCAR19_Clonotype_Frequency <- readRDS(clonotype_lineage_info_path)
  
  ### GMP CAR+ persistent clones
  pClones <- NULL
  for(i in 1:length(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]])) {
    gmp_idx <- which(colnames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) == "GMP")
    gmp_redo_idx <- which(colnames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) == "GMP-redo")
    last_gmp_idx <- max(gmp_idx, gmp_redo_idx)
    
    ### if at least GMP or GMP-redo exist and there are at least one afterward-time point
    if((last_gmp_idx != -Inf) && (ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) - last_gmp_idx > 1)) {
      ### collect persistent clones that appeared in GMP and persist afterwards
      if(nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) > 0) {
        for(j in 1:nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])) {
          for(k in last_gmp_idx:(ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])-1)) {
            if((SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,"GMP"] > 0 ||
                SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,"GMP-redo"] > 0) &&
               SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,k] > 0) {
              pClones <- c(pClones, rownames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])[j])
              break;
            }
          }
        }
      }
    }
  }
  
  ### GMP CAR+ persistent cells
  pIdx <- intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient %in% pClones),
                    intersect(union(which(Seurat_Obj@meta.data$time == "GMP"),
                                    which(Seurat_Obj@meta.data$time == "GMP-redo")),
                              which(Seurat_Obj@meta.data$CAR == "CARpos")))
  
  ### GMP CAR+ non-persistent cells
  npIdx <- setdiff(intersect(union(which(Seurat_Obj@meta.data$time == "GMP"),
                                   which(Seurat_Obj@meta.data$time == "GMP-redo")),
                             which(Seurat_Obj@meta.data$CAR == "CARpos")),
                   which(Seurat_Obj@meta.data$clonotype_id_by_patient %in% pClones))
  
  ### because there are too many cells in non-persisters
  ### randomly select some from those and perform DE analysis
  set.seed(1234)
  GMP_CARpos_Persister <- Seurat_Obj@meta.data$GMP_CARpos_Persister
  GMP_CARpos_Persister[sample(npIdx, length(npIdx) - length(pIdx))] <- NA
  
  ### get target samples
  persister_samples <- rownames(Seurat_Obj@meta.data)[which(GMP_CARpos_Persister == "YES")]
  non_persister_samples <- rownames(Seurat_Obj@meta.data)[which(GMP_CARpos_Persister == "NO")]
  
  ### remove samples that have zero values -> remove raw counts = 0
  ### only use samples that have IFITM3 counts of larger than 0
  # persisters
  temp_mat <- Seurat_Obj@assays$RNA@data["IFITM3",persister_samples]
  persister_samples <- unlist(sapply(persister_samples, function(x) {
    if(temp_mat[x] > 0) {
      return(x)
    } else {
      return(NULL)
    }
  }), use.names = FALSE)
  # non-persisters
  temp_mat <- Seurat_Obj@assays$RNA@data["IFITM3",non_persister_samples]
  non_persister_samples <- unlist(sapply(non_persister_samples, function(x) {
    if(temp_mat[x] > 0) {
      return(x)
    } else {
      return(NULL)
    }
  }), use.names = FALSE)
  # total target samples
  target_samples <- c(persister_samples, non_persister_samples)
  
  ### color for the samples based on their classes
  colors <- c(rep("blue", length(persister_samples)),
              rep("red", length(non_persister_samples)))
  
  ### for each pathway, create a correlation plot and a heatmap
  for(pathway in names(pi3k_pathways)) {
    
    ### new output directory
    outputDir3 <- paste0(outputDir2, pathway, "/")
    dir.create(outputDir3, showWarnings = FALSE, recursive = TRUE)
    
    ### if there are more than 25 genes in the pathway, select 25 randomly
    set.seed(1234)
    if(length(pi3k_pathways[[pathway]]) > 25) {
      target_genes <- sample(pi3k_pathways[[pathway]], 25)
    } else {
      target_genes <- pi3k_pathways[[pathway]]
    }
    
    ### exp matrix
    exp_mat <- as.data.frame(Seurat_Obj@assays$RNA@data[c("IFITM3", target_genes),target_samples])
    
    ### draw a combined correlation plot
    png(paste0(outputDir3, pathway, "_Correlations_Between_IFITM3_and_Others.png"),
        width = 2200, height = 1200, res = 120)
    ### for each target genes, draw a correlation plot
    par(mfrow=c(4,5), oma = c(0,0,3,0))
    cnt <- 0
    for(gene in target_genes) {
      ### remove samples that have zero values -> remove raw counts = 0
      ### only use samples that have IFITM3 counts of larger than 0
      # persisters
      temp_mat <- exp_mat[gene,persister_samples]
      persister_samples2 <- unlist(sapply(persister_samples, function(x) {
        if(temp_mat[x] > 0) {
          return(x)
        } else {
          return(NULL)
        }
      }), use.names = FALSE)
      # non-persisters
      temp_mat <- exp_mat[gene,non_persister_samples]
      non_persister_samples2 <- unlist(sapply(non_persister_samples, function(x) {
        if(temp_mat[x] > 0) {
          return(x)
        } else {
          return(NULL)
        }
      }), use.names = FALSE)
      # total target samples
      target_samples2 <- c(persister_samples2, non_persister_samples2)
      
      ### color for the samples based on their classes
      colors <- c(rep("blue", length(persister_samples2)),
                  rep("red", length(non_persister_samples2)))
      
      ### draw the correlation plot
      if((length(persister_samples2) > 0) && (length(non_persister_samples2) > 0) && cnt < 20) {
        plot(as.numeric(exp_mat["IFITM3",target_samples2]),
             as.numeric(exp_mat[gene,target_samples2]),
             pch = 19,
             col = alpha(colors, 0.3),
             main = paste0("Persisters P.Cor = ", round(cor(as.numeric(exp_mat["IFITM3",persister_samples2]),
                                                            as.numeric(exp_mat[gene,persister_samples2]),
                                                            use = "pairwise.complete.obs"), 3),
                           ", p-value = ", signif(cor.test(as.numeric(exp_mat["IFITM3",persister_samples2]),
                                                           as.numeric(exp_mat[gene,persister_samples2]))$p.value, 3),
                           "\nNon-Persisters P.Cor = ", round(cor(as.numeric(exp_mat["IFITM3",non_persister_samples2]),
                                                                  as.numeric(exp_mat[gene,non_persister_samples2]),
                                                                  use = "pairwise.complete.obs"), 3),
                           ", p-value = ", signif(cor.test(as.numeric(exp_mat["IFITM3",non_persister_samples2]),
                                                           as.numeric(exp_mat[gene,non_persister_samples2]))$p.value, 3)),
             xlab = paste("Normalized Expression of", "IFITM3"),
             ylab = paste("Normalized Beta value of", gene))
        abline(lm(as.numeric(exp_mat[gene,persister_samples2])~as.numeric(exp_mat["IFITM3",persister_samples2])), col="blue", lwd=2)
        abline(lm(as.numeric(exp_mat[gene,non_persister_samples2])~as.numeric(exp_mat["IFITM3",non_persister_samples2])), col="red", lwd=2)
        legend("topright", legend = c("Persisters", "Non-Persisters"),
               col = c("blue", "red"), pch = 19,
               title = "Sample Groups", cex = 0.8)
        cnt <- cnt + 1
      }
    }
    mtext(paste("Expression Correlations - IFITM3 vs 20 Random Pathway Genes"), outer = TRUE, cex = 2)
    dev.off()
    gc()
    
  }
  
  
  ### heatmap
  
  
  
  
}
