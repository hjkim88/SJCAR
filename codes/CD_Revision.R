###
#   File name : CD_Revision.R
#   Author    : Hyunjin Kim
#   Date      : Jan 11, 2022
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. R3C2
#                  a) DE genes between cluster 3 & 8
#                  b) DE genes between GMP that ends up with cluster 3 vs GMP that ends up with cluster 8
#               2. R1C1
#                  a) DE genes between CAR > 0 VS CAR = 0
#                  b) DE genes between CAR > 0; HIGH vs CAR > 0; low
#                  c) Does the functional groups correlate with CAR expression?
#               3. R3C5
#                  PBMC and BMMC samples are treated identically. It would be interesting to understand how much cells
#                  from these different samples are driving some of the differences in the clusters.
#               4. R4C7 & R4C8
#                  expression of CASP8 was highest in state B relative to all other states (Fig 3C)". Is it significant?
#                  higher relative expression of TOX and LAG3. Is it significant? - State C
#               5. R4C11
#                  Why were TIGIT, SELL and CD27 selected? According to the supplementary data file, they are not the most
#                  differentially expressed genes between the groups. What's the prediction accuracy of using only these 3 genes
#                  as features in the SVM classifier analysis?
#               6. R3C4
#                  Findallmarkers among our functional groups to make sure those gene markers represent each group
#               7. R2C12
#                  Fig1D - change the x-axis to only the cluster numbers
#               8. R4C1
#                  New Figure - UMAP with cell cycle score
#               9. R4C6 - Figure 2D to a table
#              10. R4C9 - find out TCRs shared across patients
#              11. R4MC2 -  add frequency of CD4/CD8 to Fig1C
#              12. R4C3 - Let’s do a giant heatmap of the top 5 (defined by absolute log fold change) DEGs per cluster
#              13. R3C5 - UMAP with PB & BM (same time point)
#                         cluster proportions (or functional group proportions) for each BM sample and compare to PB samples from the same time points
#              14. UMAP - GMP effectors - effectors to the cluster3 first vs cluster8 first - coloring differently
#              15. GSEA on precursor TF regulons - with signature from cluster3 vs cluster8 (or GMP end up in cluster 3 vs in cluster 8)
#              16. Recluster GMP without proliferating clusters - then effector precursors might be clustered together
#              17. DE analysis for PI CD8s:
#                  - Compare PI effectors vs all other PIs
#                  - Compare PI effectors vs each CD8 dysfunctional/dying cluster
#                  - Compare each dysfunctional/dying cluster to each other
#              18. Ridge plot of Tay's genes comparing cluster 3, 8, 13, & 20
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
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(monocle, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("monocle")
    require(monocle, quietly = TRUE)
  }
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    require(scales, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(caret, quietly = TRUE)) {
    install.packages("caret")
    require(caret, quietly = TRUE)
  }
  if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
  }
  if(!require(shadowtext, quietly = TRUE)) {
    install.packages("shadowtext")
    require(shadowtext, quietly = TRUE)
  }
  if(!require(reticulate, quietly = TRUE)) {
    install.packages("reticulate")
    require(reticulate, quietly = TRUE)
  }
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    require(gplots, quietly = TRUE)
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
  
  ### create outputDir
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### load Jeremy's object
  JCC_Seurat_Obj <- readRDS(file = Seurat_RObj_path)
  
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
  ###    b) DE genes between CAR > 0; HIGH vs CAR > 0; low
  ###    c) Does the functional groups correlate with CAR expression?
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
  
  ### check there are same number of px
  ### TRUE
  print(identical(unique(JCC_Seurat_Obj$px), unique(Seurat_Obj$px)))
  
  ### check they have both PB and BM
  ### TRUE
  print(identical(unique(JCC_Seurat_Obj$tissue), unique(Seurat_Obj$tissue)))
  
  ### check if CAR+ of the original is the same as the new
  ### FALSE; 185677 vs 184791
  print(identical(length(which(Seurat_Obj$CAR == "CARpos")),
                  nrow(JCC_Seurat_Obj@meta.data)))
  
  ### what are the 185677-184791 = 886 cells?
  ### No specific conditions
  View(Seurat_Obj@meta.data[setdiff(rownames(Seurat_Obj@meta.data)[which(Seurat_Obj$CAR == "CARpos")],
                                    rownames(JCC_Seurat_Obj@meta.data)),])
  
  ### check CAR expression of those cells
  ### No specific conditions
  plot(density(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",
                                            setdiff(rownames(Seurat_Obj@meta.data)[which(Seurat_Obj$CAR == "CARpos")],
                                                    rownames(JCC_Seurat_Obj@meta.data))]))
  
  ### I'll just use all
  
  ### a) DE genes between CAR > 0 VS CAR = 0
  
  ### check current CAR column is CAR > 0
  ### TRUE
  print(identical(length(which(Seurat_Obj$CAR == "CARpos")),
                  length(which(Seurat_Obj@assays$RNA@counts[1,] > 0))))
  
  ### set idents
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj$CAR)
  
  ### CAR > 0 VS CAR = 0
  de_result <- FindMarkers(Seurat_Obj,
                           ident.1 = "CARpos",
                           ident.2 = "CARneg",
                           min.pct = 0.5,
                           logfc.threshold = 0.6,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/2_DE_CARpos_vs_CARneg.xlsx"),
              sheetName = "2_DE_CARpos_vs_CARneg", row.names = FALSE)
  
  ### get entrez ids for the genes
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[which(de_result$p_val_adj < 0.01)],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("Pathways_DE_Genes_CARpos_vs_CARneg"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathways_DE_Genes_CARpos_vs_CARneg"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_Pathways_DE_genes_CARpos_vs_CARneg.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_Pathways_DE_genes_CARpos_vs_CARneg.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  ### up-regulated genes only
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[intersect(which(de_result$p_val_adj < 0.01),
                                                        which(de_result$avg_log2FC > 0))],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("UP_Pathways_DE_Genes_CARpos_vs_CARneg"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("UP_Pathways_DE_Genes_CARpos_vs_CARneg"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_UP_Pathways_DE_genes_CARpos_vs_CARneg.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_UP_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_UP_Pathways_DE_genes_CARpos_vs_CARneg.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_UP_Results"))
  
  ### down-regulated genes only
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[intersect(which(de_result$p_val_adj < 0.01),
                                                        which(de_result$avg_log2FC < 0))],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("DOWN_Pathways_DE_Genes_CARpos_vs_CARneg"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("DOWN_Pathways_DE_Genes_CARpos_vs_CARneg"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_DOWN_Pathways_DE_genes_CARpos_vs_CARneg.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_DOWN_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_DOWN_Pathways_DE_genes_CARpos_vs_CARneg.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_DOWN_Results"))
  
  
  ### see distribution of the CAR expressions
  plot(density(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",]),
       xlim = c(0,50))
  plot(density(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",]),
       xlim = c(0,10))
  
  ### 4 maybe a good threshold
  ### > 4: 126329 cells
  ### > 0 & <= 4: 59348 cells
  print(length(which(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",] > 4)))
  print(length(intersect(which(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",] > 0),
                         which(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",] <= 4))))
  
  ### set new column for High CAR vs Low CAR
  Seurat_Obj$CAR2 <- "Others"
  Seurat_Obj$CAR2[which(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",] > 4)] <- "HighCAR"
  Seurat_Obj$CAR2[intersect(which(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",] > 0),
                            which(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",] <= 4))] <- "LowCAR"
  
  ### set idents
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj$CAR2)
  
  ### High CAR vs Low CAR
  de_result <- FindMarkers(Seurat_Obj,
                           ident.1 = "HighCAR",
                           ident.2 = "LowCAR",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/2_DE_HighCAR_vs_LowCAR.xlsx"),
              sheetName = "2_DE_HighCAR_vs_LowCAR", row.names = FALSE)
  
  ### get entrez ids for the genes
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[which(de_result$p_val_adj < 0.01)],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("Pathways_DE_Genes_HighCAR_vs_LowCAR"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathways_DE_Genes_HighCAR_vs_LowCAR"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_Pathways_DE_genes_HighCAR_vs_LowCAR.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_Pathways_DE_genes_HighCAR_vs_LowCAR.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  
  ### up-regulated genes only
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[intersect(which(de_result$p_val_adj < 0.01),
                                                        which(de_result$avg_log2FC > 0))],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("UP_Pathways_DE_Genes_HighCAR_vs_LowCAR"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("UP_Pathways_DE_Genes_HighCAR_vs_LowCAR"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_UP_Pathways_DE_genes_HighCAR_vs_LowCAR.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_UP_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_UP_Pathways_DE_genes_HighCAR_vs_LowCAR.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_UP_Results"))
  
  ### down-regulated genes only
  de_entrez_ids <- mapIds(org.Hs.eg.db,
                          rownames(de_result)[intersect(which(de_result$p_val_adj < 0.01),
                                                        which(de_result$avg_log2FC < 0))],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  ### GO & KEGG
  pathway_result_GO <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                          org = "human", database = "GO",
                                          title = paste0("DOWN_Pathways_DE_Genes_HighCAR_vs_LowCAR"),
                                          displayNum = 10, imgPrint = TRUE,
                                          dir = outputDir)
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                            org = "human", database = "KEGG",
                                            title = paste0("DOWN_Pathways_DE_Genes_HighCAR_vs_LowCAR"),
                                            displayNum = 10, imgPrint = TRUE,
                                            dir = outputDir)
  write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_DOWN_Pathways_DE_genes_HighCAR_vs_LowCAR.xlsx"),
              row.names = FALSE, sheetName = paste0("GO_DOWN_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_DOWN_Pathways_DE_genes_HighCAR_vs_LowCAR.xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_DOWN_Results"))
  
  ###
  ### c) Does the functional groups correlate with CAR expression?
  ###
  
  ### NEW FUNCTIONAL ANNOTATION BY TAY - 09/27/2021
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters <- "Others"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "1", "5", "7", "9", "10", "11", "12", "15", "19", "21"))] <- "Proliferating"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("2", "4", "6", "17"))] <- "Transitioning"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("3", "8", "14"))] <- "Functional Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("13", "20"))] <- "Dysfunctional"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("16"))] <- "Early Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("18"))] <- "Metabolically Active"
  
  ### color scale
  sjcar19_colors <- c("#640B11", "#AA4C26", "#D39F3A", "#C09969", "#287B66", "#487A8F", "#3B3B53")
  names(sjcar19_colors) <- unique(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters)
  show_col(sjcar19_colors)
  
  ### violin plot of CAR EXP among the functional groups
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters)
  p <- VlnPlot(JCC_Seurat_Obj, features = "JCC-SJCAR19short",
               pt.size = 0, cols = sjcar19_colors)
  p[[1]] <- p[[1]] + geom_boxplot(width=0.1) +
    # stat_compare_means(size = 8) +
    xlab("") +
    # stat_summary(fun=mean, geom="point", size=3, color="black") +
    theme_classic(base_size = 40) +
    theme(legend.key.size = unit(3, 'cm'),
          legend.position = "none",
          legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          legend.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35, color = "black", face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 25, color = "black", face = "bold"),
          axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
  
  ### save the violin plot
  ggsave(file = paste0(outputDir, "2_Violin_CAR_Functional_Groups.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### Now GMP precursor vs other CD8 GMP
  ### are there any CAR EXP differences?
  
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
  
  ### GMP subsisters end up in cluster3 and 8 vs other GMP subsisters
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
  
  ### violin plot of CAR EXP between precursor vs others
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 <- factor(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8,
                                                                    levels = c("GMP_Subsisters_End_Up_In_Cluster_3_And_8",
                                                                               "Other_CD8_GMPs",
                                                                               "Others"))
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8)
  p <- VlnPlot(JCC_Seurat_Obj, features = "JCC-SJCAR19short",
               pt.size = 0, cols = c("#640B11", "#487A8F", "#D39F3A"))
  p[[1]] <- p[[1]] + geom_boxplot(width=0.1) +
    # stat_compare_means(size = 8) +
    xlab("") +
    # stat_summary(fun=mean, geom="point", size=3, color="black") +
    theme_classic(base_size = 40) +
    theme(legend.key.size = unit(3, 'cm'),
          legend.position = "none",
          legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          legend.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35, color = "black", face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 25, color = "black", face = "bold"),
          axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
  
  ### save the violin plot
  ggsave(file = paste0(outputDir, "2_Violin_CAR_Precursor_vs_Others.png"), plot = p, width = 20, height = 10, dpi = 350)
  
  ### DE analysis between precursor vs other CD8 GMPs
  ### See CAR has significant FDR in there
  
  ### set idents
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8)
  
  ### Precursor vs Other CD8 GMPs
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "GMP_Subsisters_End_Up_In_Cluster_3_And_8",
                           ident.2 = "Other_CD8_GMPs",
                           min.pct = 0.2,
                           logfc.threshold = 0,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/2_DE_Precursor_vs_Other_CD8_GMPs.xlsx"),
              sheetName = "2_DE_Precursor_vs_Other_CD8_GMPs", row.names = FALSE)
  
  
  ###
  ### 3. PBMC and BMMC samples are treated identically. It would be interesting to understand how much cells
  ### from these different samples are driving some of the differences in the clusters.
  ###
  ### Show UMAP of PB and BM cells with different color and with labeling cluster/functional annotations
  ### See violin plots of the marker genes in each group - PB vs BM
  ###
  
  ### UMAP
  p <- list()
  p[[1]] <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "tissue",
               pt.size = 3, cols = c("GMP" = "gray", "PB" = "#640B11", "BM" = "#487A8F"),
               order = c("BM", "PB", "GMP")) +
    ggtitle("PB & BM") +
    labs(color="") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40, color = "black", face = "bold"),
          axis.title = element_text(size = 30, color = "black", face = "bold"),
          axis.text = element_text(size = 25, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  p[[1]][[1]]$layers[[1]]$aes_params$alpha <- 0.8
  
  p[[2]] <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap", raster = TRUE,
                    group.by = "tissue",
                    pt.size = 3, cols = c("GMP" = "gray", "PB" = "#640B11", "BM" = "gray"),
                    order = c("PB", "BM", "GMP")) +
    ggtitle("PB ONLY") +
    labs(color="") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40, color = "black", face = "bold"),
          axis.title = element_text(size = 30, color = "black", face = "bold"),
          axis.text = element_text(size = 25, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  p[[2]][[1]]$layers[[1]]$aes_params$alpha <- 0.8
  
  p[[3]] <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap", raster = TRUE,
                    group.by = "tissue",
                    pt.size = 3, cols = c("GMP" = "gray", "PB" = "gray", "BM" = "#487A8F"),
                    order = c("BM", "PB", "GMP")) +
    ggtitle("BM ONLY") +
    labs(color="") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40, color = "black", face = "bold"),
          axis.title = element_text(size = 30, color = "black", face = "bold"),
          axis.text = element_text(size = 25, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  p[[3]][[1]]$layers[[1]]$aes_params$alpha <- 0.8
  
  ### draw the UMAP
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   ncol = 2)
  ggsave(file = paste0(outputDir, "3_UMAP_PB_BM.png"), g, width = 18, height = 10, dpi = 350)
  
  ### BM is only from specific time points
  print(unique(JCC_Seurat_Obj$time2[which(JCC_Seurat_Obj$tissue == "BM")]))
  JCC_Seurat_Obj$time3 <- JCC_Seurat_Obj$time2
  JCC_Seurat_Obj$time3[which(!JCC_Seurat_Obj$time2 %in% c("Wk4", "3mo"))] <- "NA"
  p <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap", raster = TRUE,
          group.by = "time3",
          pt.size = 3, cols = c("NA" = "gray", "Wk4" = "#640B11", "3mo" = "#640B11"),
          order = c("3mo", "Wk4", "NA")) +
  ggtitle("Wk4 & 3mo") +
    labs(color="") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40, color = "black", face = "bold"),
          axis.title = element_text(size = 30, color = "black", face = "bold"),
          axis.text = element_text(size = 25, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.8
  
  ggsave(paste0(outputDir, "3_UMAP_Wk4_3mo.png"), plot = p, width = 13, height = 10, dpi = 350)
  
  ### See violin plots of the marker genes in each group - PB vs BM
  ### Wk4 & 3mo data only
  
  ### what are the clusters in Wk4 & 3mo?
  ### all the clusters except 21
  print(unique(JCC_Seurat_Obj$AllSeuratClusters[which(JCC_Seurat_Obj$time2 %in% c("Wk4", "3mo"))]))
  
  ### dotplot/violin plot similar to Fig1D but PB/BM separated
  ### if there is a lot difference between PB and BM, that would be a problem
  ### meaning one cell type if driving the difference
  
  ### replicate Fig1D
  interesting_genes2 <- c("RPL32", "RPL30", "LAG3", "TOX", "CASP8", "IL7R", "SELL", "BNIP3", "MKI67",
                          "CDC20", "CDK1", "NKG7", "GNLY", "GZMH", "GZMM", "GZMK")
  p <- DotPlot(JCC_Seurat_Obj,
               features = interesting_genes2,
               group.by = "New_Functional_Annotation_Based_On_Clusters") +
    scale_size(range = c(5, 35)) +
    xlab("") +
    ylab("") +
    scale_color_gradientn(colours = c("#487A8F", "#C09969", "#AA4C26")) +
    guides(color = guide_colorbar(title = "Relative Expression")) +
    theme_classic(base_size = 28) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = -45, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'))
  
  ### label the groups with PB/BM
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2 <- paste0(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters,
                                                                        "_", JCC_Seurat_Obj$tissue)
  
  ### choose Wk4/3mo cells only
  temp_obj <- subset(JCC_Seurat_Obj,
                     cells = rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$time2 %in% c("Wk4", "3mo"))])
  
  ### make an order in the column
  temp_obj$New_Functional_Annotation_Based_On_Clusters2 <- factor(temp_obj$New_Functional_Annotation_Based_On_Clusters2)
  
  ### dotplot with PB and BM separately
  p <- DotPlot(temp_obj,
               features = interesting_genes2,
               group.by = "New_Functional_Annotation_Based_On_Clusters2") +
    scale_size(range = c(5, 20)) +
    xlab("") +
    ylab("") +
    scale_color_gradientn(colours = c("#487A8F", "#C09969", "#AA4C26")) +
    guides(color = guide_colorbar(title = "Relative Expression")) +
    theme_classic(base_size = 28) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = -45, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'))
  ggsave(file = paste0(outputDir, "3_Dotplot_PB_BM.png"),
         plot = p, width = 25, height = 10, dpi = 350)
  
  ### violin plot with PB and BM separately
  temp_obj <- SetIdent(object = temp_obj,
                       cells = rownames(temp_obj@meta.data),
                       value = temp_obj$New_Functional_Annotation_Based_On_Clusters2)
  p <- VlnPlot(temp_obj, features = interesting_genes2, pt.size = 0)
  for(i in 1:length(interesting_genes2)) {
    p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
      xlab("") +
      theme_classic() +
      theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
            legend.position = "none")
  }
  
  ### save the violin plot
  ggsave(file = paste0(outputDir, "3_Violin_PB_BM.png"), plot = p, width = 20, height = 15, dpi = 350)
  
  
  ###
  ### 4. expression of CASP8 was highest in state B relative to all other states (Fig 3C)". Is it significant?
  #   higher relative expression of TOX and LAG3. Is it significant? - State C
  ###
  
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
  
  ### the data is too big (can't run for monocle), so we down-sample the cells in each time point
  ### BUT THIS TIME, INCLUDE ALL THE SUBSISTERS (IN ADDITION TO THE ORIGINAL DOWNSAMPLED ONES)
  JCC_Seurat_Obj@meta.data$downsampled2 <- JCC_Seurat_Obj@meta.data$downsampled
  JCC_Seurat_Obj@meta.data$downsampled2[which(JCC_Seurat_Obj@meta.data$ALL_CARpos_Persister == "YES")] <- "YES"
  
  ### Construct a monocle cds
  monocle_metadata <- JCC_Seurat_Obj@meta.data[rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj@meta.data$downsampled2 == "YES")],]
  monocle_metadata$time2 <- factor(monocle_metadata$time2, levels = unique(monocle_metadata$time2))
  monocle_cds2 <- newCellDataSet(as(as.matrix(JCC_Seurat_Obj@assays$RNA@data[,rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj@meta.data$downsampled2 == "YES")]]), 'sparseMatrix'),
                                 phenoData = new('AnnotatedDataFrame', data = monocle_metadata),
                                 featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(JCC_Seurat_Obj@assays$RNA@data),
                                                                                           row.names = row.names(JCC_Seurat_Obj@assays$RNA@data),
                                                                                           stringsAsFactors = FALSE, check.names = FALSE)),
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())
  
  ### run monocle
  monocle_cds2 <- estimateSizeFactors(monocle_cds2)
  monocle_cds2 <- estimateDispersions(monocle_cds2)
  monocle_cds2 <- reduceDimension(monocle_cds2, reduction_method = "DDRTree")
  monocle_cds2 <- orderCells(monocle_cds2)
  
  ### determine the beginning state
  plot_cell_trajectory(monocle_cds2, color_by = "time2")
  plot_cell_trajectory(monocle_cds2, color_by = "State")
  plot_complex_cell_trajectory(monocle_cds2, color_by = "time2")
  plot_complex_cell_trajectory(monocle_cds2, color_by = "State")
  ### this root state should be checked if we want to rerun this -> Sometimes it's 6 and sometimes it's 1
  monocle_cds2 <- orderCells(monocle_cds2, root_state = "6")
  
  ### change the state numbers to alphabets
  monocle_cds2$State <- as.character(monocle_cds2$State)
  monocle_cds2$State[which(monocle_cds2$State == "6")] <- "A"
  monocle_cds2$State[which(monocle_cds2$State == "5")] <- "B"
  monocle_cds2$State[which(monocle_cds2$State == "4")] <- "C"
  monocle_cds2$State[which(monocle_cds2$State == "7")] <- "D"
  monocle_cds2$State[which(monocle_cds2$State == "3")] <- "E"
  monocle_cds2$State[which(monocle_cds2$State == "1")] <- "F"
  monocle_cds2$State[which(monocle_cds2$State == "2")] <- "G"
  monocle_cds2$State <- factor(monocle_cds2$State, levels = c("A", "B", "C", "D", "E", "F", "G"))
  plot_cell_trajectory(monocle_cds2, color_by = "State")
  
  ### make the down-sampled seurat object
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj@meta.data$downsampled2)
  downsampled_seurat_obj <- subset(JCC_Seurat_Obj,
                                   idents = "YES")
  downsampled_seurat_obj$monocle_state <- monocle_cds2@phenoData@data[rownames(downsampled_seurat_obj@meta.data),"State"]
  
  print(identical(rownames(downsampled_seurat_obj@meta.data), colnames(downsampled_seurat_obj@assays$RNA@counts)))
  print(identical(names(Idents(object = downsampled_seurat_obj)), rownames(downsampled_seurat_obj@meta.data)))
  
  ### color scale
  sjcar19_colors <- c("#640B11", "#AA4C26", "#D39F3A", "#C09969", "#287B66", "#487A8F", "#3B3B53")
  names(sjcar19_colors) <- levels(monocle_cds2$State)
  show_col(sjcar19_colors)
  
  ### draw monocle plot for checking
  downsampled_seurat_obj$monocle_state <- factor(downsampled_seurat_obj$monocle_state, levels = rev(levels(downsampled_seurat_obj$monocle_state)))
  DotPlot(downsampled_seurat_obj,
               features = "CASP8",
               group.by = "monocle_state") +
    scale_size(range = c(5, 35)) +
    xlab("") +
    ylab("State") +
    scale_color_gradientn(colours = c("#3B3B53", "#D39F3A", "#640B11"),
                          n.breaks = 3) +
    theme_classic(base_size = 28) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(hjust = 0.5, size = 35, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'))
  
  ### now everything is confirmed so let's begin the analysis
  downsampled_seurat_obj$R4C7 <- "Others"
  downsampled_seurat_obj$R4C7[which(as.character(downsampled_seurat_obj$monocle_state) == "B")] <- "B"
  
  ### set idents
  downsampled_seurat_obj <- SetIdent(object = downsampled_seurat_obj,
                                     cells = rownames(downsampled_seurat_obj@meta.data),
                                     value = downsampled_seurat_obj$R4C7)
  
  ### State B vs Others
  de_result <- FindMarkers(downsampled_seurat_obj,
                           ident.1 = "B",
                           ident.2 = "Others",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/4_DE_StateB_vs_Others.xlsx"),
              sheetName = "4_DE_StateB_vs_Others", row.names = FALSE)
  
  ### now everything is confirmed so let's begin the analysis
  downsampled_seurat_obj$R4C8 <- "Others"
  downsampled_seurat_obj$R4C8[which(as.character(downsampled_seurat_obj$monocle_state) == "C")] <- "C"
  
  ### set idents
  downsampled_seurat_obj <- SetIdent(object = downsampled_seurat_obj,
                                     cells = rownames(downsampled_seurat_obj@meta.data),
                                     value = downsampled_seurat_obj$R4C8)
  
  ### State B vs Others
  de_result <- FindMarkers(downsampled_seurat_obj,
                           ident.1 = "C",
                           ident.2 = "Others",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/4_DE_StateC_vs_Others.xlsx"),
              sheetName = "4_DE_StateC_vs_Others", row.names = FALSE)
  
  
  ###
  ### 5. R4C11
  #      Why were TIGIT, SELL and CD27 selected? According to the supplementary data file, they are not the most
  #      differentially expressed genes between the groups. What's the prediction accuracy of using only these 3 genes
  #      as features in the SVM classifier analysis?
  ###
  
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
  
  ### parameter setting for a classifier
  iteration <- 100
  set.seed(100) # 100-i:26 error, 4321-i:41 error, 1000-i:93 error, 2000-i:95 error, 3000-i:100
  sampleNum <- 100
  genes <- c("TIGIT", "SELL", "CD27")
  methodTypes <- "svmRadial"
  methodNames <- "SVMRadial"
  train_control <- trainControl(method="LOOCV", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
  
  ### iteratively build a classifier
  p <- vector("list", length = iteration)
  acc <- NULL
  auc <- NULL
  for(i in 1:iteration) {
    
    ### progress
    if(i %% 10 == 0) {
      writeLines(paste(i, "/", iteration))
    }
    
    ### normalize the read counts
    ### before the normalization, only keep the samples that will be used in the classifier
    random_samples <- c(sample(which(JCC_Seurat_Obj@meta.data$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8"), sampleNum),
                        sample(which(JCC_Seurat_Obj@meta.data$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 == "Other_CD8_GMPs"), sampleNum))
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(JCC_Seurat_Obj@assays$RNA@counts[genes,random_samples],
                                                                stringsAsFactors = FALSE, check.names = FALSE)+1,
                                         filter_thresh = 1)
    
    # ### reduce the gene size based on variance
    # ### only select high variance genes
    # input_data <- selectTopV(input_data, featureSelectionNum)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(c(rep("GMP_Last", sampleNum),
                                 rep("GMP_Not_Last", sampleNum)),
                               levels = c("GMP_Last", "GMP_Not_Last"))
    
    ### build classifier and test
    ### LOOCV
    model <- train(Class~., data=input_data, trControl=train_control, method=methodTypes)
    roc <- roc(model$pred$obs, model$pred$GMP_Last)
    acc <- c(acc, round(mean(model$results$Accuracy), 3))
    auc <- c(auc, round(as.numeric(roc$auc), 3))
    p[[i]] <- plot.roc(roc, main = paste(methodNames, "Using Gene Expressions\n",
                                         "Accuracy =", acc[i]),
                       legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                       xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    gc()
    
  }
  
  ### save RDS file
  saveRDS(p, file = paste0(outputDir, "5_Classifier_GMP_Precursor_vs_Other_CD8_GMP_3_genes_100_100.RDS"))
  saveRDS(roc, file = paste0(outputDir, "5_Classifier_GMP_Precursor_vs_Other_CD8_GMP_3_genes_100_100_roc.RDS"))
  saveRDS(acc, file = paste0(outputDir, "5_Classifier_GMP_Precursor_vs_Other_CD8_GMP_3_genes_100_100_acc.RDS"))
  saveRDS(auc, file = paste0(outputDir, "5_Classifier_GMP_Precursor_vs_Other_CD8_GMP_3_genes_100_100_auc.RDS"))
  
  ### draw ROC curves
  png(paste0(outputDir, "5_Classifier_GMP_Precursor_vs_Other_CD8_GMP_3_genes_100.png"),
      width = 2000, height = 2000, res = 350)
  par(mfrow=c(4, 3))
  for(i in 1:iteration) {
    plot.roc(p[[i]], main = paste(methodNames, "Using Gene Expressions\n",
                                  "Accuracy =", acc[i]),
             legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
             xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
  }
  dev.off()
  
  ### print average ACC and AUC
  ### "3-Gene-Classifier Avearge ACC: 0.72303 Average AUC: 0.73634"
  print(paste0("3-Gene-Classifier Avearge ACC: ", mean(acc), " Average AUC: ", mean(auc)))
  
  
  ###
  ### 6. R3C4
  #   Findallmarkers among our functional groups to make sure those gene markers represent each group
  ###
  
  ### NEW FUNCTIONAL ANNOTATION BY TAY - 09/27/2021
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters <- "Others"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "1", "5", "7", "9", "10", "11", "12", "15", "19", "21"))] <- "Proliferating"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("2", "4", "6", "17"))] <- "Transitioning"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("3", "8", "14"))] <- "Functional Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("13", "20"))] <- "Dysfunctional"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("16"))] <- "Early Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("18"))] <- "Metabolically Active"
  
  ### set functional groups as idents
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj@meta.data$New_Functional_Annotation_Based_On_Clusters)
  
  ### DE analysis
  de_result <- FindAllMarkers(JCC_Seurat_Obj,
                              min.pct = 0.3,
                              logfc.threshold = 0.3,
                              test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/R3c4_Functional_Groups_AllMarkers.xlsx"),
              sheetName = "R3c4_Functional_Groups_AllMarkers", row.names = FALSE)
  
  
  ###
  ### 7. R2C12
  #      Fig1D - change the x-axis to only the cluster numbers
  ###
  
  ### dot plot2 - like Dave's
  interesting_genes2 <- c("RPL32", "RPL30", "LAG3", "TOX", "CASP8", "IL7R", "SELL", "BNIP3", "MKI67",
                          "CDC20", "CDK1", "NKG7", "GNLY", "GZMH", "GZMM", "GZMK")
  
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2 <- JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Proliferating")] <- "Clusters{0,1,5,7,9,10,11,12,15,19,21}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Transitioning")] <- "Clusters{2,4,6,17}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Metabolically Active")] <- "Clusters{18}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Early Effector")] <- "Clusters{16}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Dysfunctional")] <- "Clusters{13, 20}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Functional Effector")] <- "Clusters{3,8,10}"
  
  p <- DotPlot(JCC_Seurat_Obj,
               features = interesting_genes2,
               group.by = "New_Functional_Annotation_Based_On_Clusters2") +
    scale_size(range = c(5, 45)) +
    xlab("") +
    ylab("") +
    scale_color_gradientn(colours = c("#487A8F", "#C09969", "#AA4C26")) +
    guides(color = guide_colorbar(title = "Relative Expression")) +
    coord_flip() +
    theme_classic(base_size = 28) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = -90, size = 35, vjust = 0.5, hjust = 0, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 25, color = "black", face = "bold"),
          legend.text = element_text(size = 20, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'),
          legend.position = "top")
  ggsave(file = paste0(outputDir, "R2C12_Fig1D_Dotplot_Functional_Group_GEXP.pdf"),
         plot = p, width = 15, height = 35, dpi = 350)
  
  
  ###
  ###  8. R4c1
  #       New Figure - UMAP with cell cycle score
  ###
  
  ### normalization
  JCC_Seurat_Obj <- NormalizeData(JCC_Seurat_Obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### Cell cycle score (will be used later for regression out)
  JCC_Seurat_Obj <- CellCycleScoring(object = JCC_Seurat_Obj,
                                     g2m.features = cc.genes$g2m.genes,
                                     s.features = cc.genes$s.genes)
  
  ### color scale
  sjcar19_colors <- c("#640B11", "#AA4C26", "#D39F3A", "#C09969", "#287B66", "#487A8F", "#3B3B53")
  show_col(sjcar19_colors)
  
  ### draw a UMAP with cell cycle scores
  temp <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
                  group.by = "AllSeuratClusters", label = TRUE,
                  raster = TRUE, pt.size = 3)
  p <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
               group.by = "Phase",
               cols = c("G2M" = "#487A8F", "S" = "#D39F3A", "G1" = "#640B11"),
               order = c("G2M", "S", "G1"),
               raster = TRUE, pt.size = 1) +
    ggtitle("") +
    labs(color="") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  p <- p + geom_shadowtext(data = temp$layers[[2]]$data, aes(x = UMAP_1, y = UMAP_2, label=AllSeuratClusters),
                           size=10, color="cornsilk2", bg.color="black", bg.r=0.2)
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir, "R4C1_UMAP_Cell_Cyel_Phase.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ###
  ### 9. R4C6 - Figure 2D to a table
  ###
  
  ### prepare table for the plot
  plot_df <- data.frame(Cluster=as.character(sapply(levels(JCC_Seurat_Obj$AllSeuratClusters), function(x) rep(x, length(unique(JCC_Seurat_Obj$time2))))),
                        Time=as.character(rep(unique(JCC_Seurat_Obj$time2), length(levels(JCC_Seurat_Obj$AllSeuratClusters)))),
                        Numbers=0,
                        Pcnt=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### calculate numbers
  cnt <- 1
  for(clstr in levels(JCC_Seurat_Obj$AllSeuratClusters)) {
    for(tp in unique(JCC_Seurat_Obj$time2)) {
      plot_df$Numbers[cnt] <- length(intersect(which(JCC_Seurat_Obj$AllSeuratClusters == clstr),
                                               which(JCC_Seurat_Obj$time2 == tp)))
      cnt <- cnt + 1
    }
  }
  
  ### order the plot_df in a GMP/PI - oriented way
  plot_df <- plot_df[order(plot_df$Time),]
  
  ### calculate percentages
  time_sum <- rep(0, length(unique(plot_df$Time)))
  names(time_sum) <- unique(plot_df$Time)
  for(i in 1:length(unique(plot_df$Time))) {
    time_sum[i] <- sum(plot_df[which(plot_df$Time == unique(plot_df$Time)[i]),"Numbers"])
    plot_df$Pcnt[which(plot_df$Time == unique(plot_df$Time)[i])] <- round(plot_df$Numbers[which(plot_df$Time == unique(plot_df$Time)[i])] * 100 / time_sum[i], 1)
  }
  
  ### annotate "%"
  plot_df$Pcnt <- paste0(plot_df$Pcnt, "%")
  
  ### remove the Wk6
  plot_df <- plot_df[which(!plot_df$Time %in% c("Wk6")),]
  
  ### factorize the time point & state
  plot_df$Time <- factor(plot_df$Time, levels = unique(JCC_Seurat_Obj$time2))
  plot_df$Cluster <- factor(plot_df$Cluster, levels = unique(plot_df$Cluster))
  
  ### refine the plot_df
  plot_df2 <- data.frame(Cluster=levels(plot_df$Cluster),
                         matrix(0, length(unique(plot_df$Cluster)), length(unique(plot_df$Time))),
                         stringsAsFactors = FALSE, check.names = FALSE)
  colnames(plot_df2) <- c("Cluster", as.character(unique(plot_df$Time)))
  rownames(plot_df2) <- as.character(unique(plot_df$Cluster))
  for(clstr in as.character(unique(plot_df$Cluster))) {
    for(time in as.character(unique(plot_df$Time))) {
      plot_df2[clstr,time] <- plot_df[intersect(which(plot_df$Cluster == clstr),
                                                which(plot_df$Time == time)),"Pcnt"]
    }
  }
  
  ### reordering
  plot_df2 <- plot_df2[,c("Cluster", "GMP", "Wk1", "Wk2", "Wk3", "Wk4", "Wk8", "3mo", "6mo")]
  
  ### write out the table
  write.xlsx2(plot_df2,
              file = paste0(outputDir, "R4C6_Cluster_Time_Bar_Plot_Table.xlsx"),
              sheetName = "Cluster_Time_Bar_Plot_Table", row.names = FALSE)
  
  
  ###
  ### 10. R4C9 - find out TCRs shared across patients (at least one alpha & beta)
  ###
  
  shared_tcrs <- data.frame(matrix(0, length(unique(JCC_Seurat_Obj$px)), length(unique(JCC_Seurat_Obj$px))),
                            stringsAsFactors = FALSE, check.names = FALSE)
  rownames(shared_tcrs) <- unique(JCC_Seurat_Obj$px)
  colnames(shared_tcrs) <- unique(JCC_Seurat_Obj$px)
  
  for(row_px in unique(JCC_Seurat_Obj$px)) {
    for(col_px in unique(JCC_Seurat_Obj$px)) {
      tcr_intersect <- intersect(JCC_Seurat_Obj$cdr3_one_alpha_beta[which(JCC_Seurat_Obj$px == row_px)],
                                 JCC_Seurat_Obj$cdr3_one_alpha_beta[which(JCC_Seurat_Obj$px == col_px)])
      tcr_intersect <- tcr_intersect[which(!is.na(tcr_intersect))]
      
      single_chain_idx <- sapply(tcr_intersect, function(x) {
        if(is.na(x)) {
          return(0)
        }
        temp <- strsplit(x, ";", fixed = TRUE)[[1]]
        if(length(temp) > 1) {
          return(2)
        } else {
          return(1)
        }
      })
      
      tcr_intersect <- tcr_intersect[which(single_chain_idx == 2)]
      
      if(row_px == col_px) {
        shared_tcrs[row_px,col_px] <- length(tcr_intersect)
      } else {
        shared_tcrs[row_px,col_px] <- paste(c(length(tcr_intersect), tcr_intersect), collapse = ",")
      }
    }
  }
  
  ### write out the result
  write.xlsx2(data.frame(Patient=rownames(shared_tcrs),
                         shared_tcrs,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/R4C9_Shared_TCRs_among_pxs.xlsx"),
              sheetName = "R4C9_R4C9_Shared_TCRs_among_pxs", row.names = FALSE)
  
  ### get indices of cells that have single chain
  single_chain_idx <- sapply(JCC_Seurat_Obj$cdr3_one_alpha_beta, function(x) {
    if(is.na(x)) {
      return(0)
    }
    temp <- strsplit(x, ";", fixed = TRUE)[[1]]
    if(length(temp) > 1) {
      return(2)
    } else {
      return(1)
    }
  })
  
  
  ###
  ### 11. R4MC2 -  add frequency of CD4/CD8 to Fig1C
  ###              how many CD4 and CD8 cells in each cluster?
  ###
  
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
  
  CD4_pcnt <- sapply(levels(JCC_Seurat_Obj$AllSeuratClusters), function(x) {
    return(length(intersect(which(JCC_Seurat_Obj$AllSeuratClusters == x),
                            which(JCC_Seurat_Obj@assays$RNA@counts["CD4",] > 0))) * 100 / length(which(JCC_Seurat_Obj$AllSeuratClusters == x)))
  }, USE.NAMES = TRUE)
  
  CD8_pcnt <- sapply(levels(JCC_Seurat_Obj$AllSeuratClusters), function(x) {
    return(length(intersect(which(JCC_Seurat_Obj$AllSeuratClusters == x),
                            which(JCC_Seurat_Obj@assays$RNA@counts["CD8A",] > 0))) * 100 / length(which(JCC_Seurat_Obj$AllSeuratClusters == x)))
  }, USE.NAMES = TRUE)
  
  CD4_pcnt <- round(CD4_pcnt, 2)
  CD8_pcnt <- round(CD8_pcnt, 2)
  
  CD4_CD8_pcnt <- data.frame(Cluster=names(CD4_pcnt),
                             CD4_Percentage=CD4_pcnt,
                             CD8_Percentage=CD8_pcnt,
                             stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the result
  write.xlsx2(CD4_CD8_pcnt,
              file = paste0(outputDir, "/R4MC2_CD4_CD8_Pcnt.xlsx"),
              sheetName = "R4MC2_CD4_CD8_Pcnt", row.names = FALSE)
  
  ###
  ### 12. R4C3 - Let’s do a giant heatmap of the top 5 (defined by absolute log fold change) DEGs per cluster
  ###
  
  ### load the DE result
  allmarkers_cd4_de_result <- read.table(file = paste0(outputDir, "/JCC212_CARpos_JCC_CD4_allMarkers.tsv"),
                                         header = TRUE,
                                         stringsAsFactors = FALSE, check.names = FALSE)
  allmarkers_cd8_de_result <- read.table(file = paste0(outputDir, "/JCC212_CARpos_JCC_CD8_allMarkers.tsv"),
                                         header = TRUE,
                                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### aggregate the result - top5 from each cluster
  top_5_genes <- vector("list", length(unique(JCC_Seurat_Obj$AllSeuratClusters)))
  names(top_5_genes) <- levels(JCC_Seurat_Obj$AllSeuratClusters)
  for(clstr in names(top_5_genes)) {
    if(clstr %in% c("4", "21")) {
      top_5_genes[[clstr]] <- intersect(allmarkers_cd8_de_result$gene[which(allmarkers_cd8_de_result$cluster == clstr)][1:15],
                                        allmarkers_cd4_de_result$gene[which(allmarkers_cd4_de_result$cluster == clstr)][1:15])[1:5]
    } else if(clstr %in% unique(allmarkers_cd4_de_result$cluster)) {
      top_5_genes[[clstr]] <- allmarkers_cd4_de_result$gene[which(allmarkers_cd4_de_result$cluster == clstr)][1:5]
    } else if(clstr %in% unique(allmarkers_cd8_de_result$cluster)) {
      top_5_genes[[clstr]] <- allmarkers_cd8_de_result$gene[which(allmarkers_cd8_de_result$cluster == clstr)][1:5]
    } else {
      writeLines(paste0("Error"))
    }
  }
  
  ### prepare a heatmap data
  heatdata <- data.frame(matrix(0, nrow = length(unlist(top_5_genes)), ncol = length(levels(JCC_Seurat_Obj$AllSeuratClusters))),
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(heatdata) <- paste0(as.vector(sapply(levels(JCC_Seurat_Obj$AllSeuratClusters), function(x) rep(x, 5))),
                               "_",
                               unlist(top_5_genes))
  colnames(heatdata) <- levels(JCC_Seurat_Obj$AllSeuratClusters)
  for(row in rownames(heatdata)) {
    target_gene <- strsplit(row, split = "_", fixed = TRUE)[[1]][2]
    for(col in colnames(heatdata)) {
      heatdata[row,col] <- mean(JCC_Seurat_Obj@assays$RNA@data[target_gene,
                                                               rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$AllSeuratClusters == col)]])
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
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(levels(JCC_Seurat_Obj$AllSeuratClusters), hue_pal())
  heatclus <- as.vector(sapply(levels(JCC_Seurat_Obj$AllSeuratClusters), function(x) rep(x, 5)))
  
  ### make a heatmap with those genes
  png(paste0(outputDir, "R4C3_Heatmap_Top5_From_Each_Cluster.png"),
      width = 3000, height = 3000, res = 350)
  par(oma=c(0,2,0,3))
  heatmap.2(as.matrix(heatdata), col = colorpanel(36, low = "skyblue", high = "red"),
            scale = "row", dendrogram = "none", trace = "none",
            cexRow = 0.4, cexCol = 2, key.title = "", main = "Top 5 Genes From Each Cluster",
            Rowv = FALSE, labRow = unlist(top_5_genes),
            Colv = FALSE, labCol = levels(JCC_Seurat_Obj$AllSeuratClusters),
            key.xlab = "Norm.Count", key.ylab = "Frequency",
            RowSideColors = cell_colors_clust[heatclus])
  legend("left", inset = -0.1, y.intersp = 0.7,
         box.lty = 0, cex = 1.1,
         title = "Marker Cluster", xpd = TRUE,
         legend=names(cell_colors_clust),
         col=cell_colors_clust,
         pch=15)
  dev.off()
  
  
  ###
  ### 13. R3C5 - UMAP with PB & BM (same time point)
  #              cluster proportions (or functional group proportions) for each BM sample and compare to PB samples from the same time points
  ###
  
  ### which time points are in BM tissues?
  print(unique(JCC_Seurat_Obj$time2[which(JCC_Seurat_Obj$tissue == "BM")]))
  
  ### subset - wk4 and 3mo
  subset_seurat_obj <- subset(JCC_Seurat_Obj,
                              cells = rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$time2 %in% c("Wk4", "3mo"))])
  
  ### UMAP coloring with PB and BM
  p <- DimPlot(object = subset_seurat_obj, reduction = "umap",
               group.by = "tissue",
               pt.size = 1, raster = TRUE) +
    ggtitle("") +
    labs(color="Tissue") +
    theme_classic(base_size = 36) +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.9
  ggsave(paste0(outputDir, "R3C5_UMAP_Wk3_3mo_Tissue.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### cluster proportions
  plot_df <- data.frame(Clusters=levels(subset_seurat_obj$AllSeuratClusters),
                        PB=rep(0, length(levels(subset_seurat_obj$AllSeuratClusters))),
                        BM=rep(0, length(levels(subset_seurat_obj$AllSeuratClusters))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  rownames(plot_df) <- levels(subset_seurat_obj$AllSeuratClusters)
  
  ### fill in the table
  for(clstr in levels(subset_seurat_obj$AllSeuratClusters)) {
    plot_df[clstr, "PB"] <- length(intersect(which(subset_seurat_obj$tissue == "PB"),
                                             which(subset_seurat_obj$AllSeuratClusters == clstr)))
    plot_df[clstr, "BM"] <- length(intersect(which(subset_seurat_obj$tissue == "BM"),
                                             which(subset_seurat_obj$AllSeuratClusters == clstr)))
  }
  
  pb_sum <- sum(plot_df$PB)
  bm_sum <- sum(plot_df$BM)
  
  plot_df$PB_Pcnt <- round(plot_df$PB * 100 / pb_sum, 2)
  plot_df$BM_Pcnt <- round(plot_df$BM * 100 / bm_sum, 2)
  
  ### data table for the plot
  plot_df2 <- data.frame(rbind(as.matrix(plot_df[,c("Clusters", "PB", "PB_Pcnt")]),
                               as.matrix(plot_df[,c("Clusters", "BM", "BM_Pcnt")])),
                         stringsAsFactors = FALSE, check.names = FALSE)
  plot_df2$Tissue <- c(rep("PB", length(levels(subset_seurat_obj$AllSeuratClusters))),
                       rep("BM", length(levels(subset_seurat_obj$AllSeuratClusters))))
  colnames(plot_df2) <- c("Clusters", "Count", "Pcnt", "Tissue")
  
  ### numerize and factorize some columns
  plot_df2$Count <- as.numeric(plot_df2$Count)
  plot_df2$Pcnt <- as.numeric(plot_df2$Pcnt)
  plot_df2$Clusters <- factor(plot_df2$Clusters, levels = levels(subset_seurat_obj$AllSeuratClusters))
  
  ### filter some cluster labels
  plot_df2$Clusters2 <- as.character(plot_df2$Clusters)
  plot_df2$Clusters2[which(as.numeric(plot_df2$Pcnt) < 5)] <- ""
  
  ### proportional bar plots
  p <- ggplot(data=plot_df2, aes_string(x="Tissue", y="Pcnt", fill="Clusters", label="Clusters2")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Cluster Proportions") +
    geom_text(size = 10, position = position_stack(vjust = 0.5)) +
    xlab("") +
    ylab("Percentage") +
    coord_flip() +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(size = 35, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir, "R3C5_PB_BM_Cluster_Proportions.png"), plot = p,
         width = 20, height = 10, dpi = 350)
  
  
  ###
  ### 14. UMAP - GMP effectors - effectors to the cluster3 first vs cluster8 first - coloring differently
  ###
  
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
  ### but no duplicates this time
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8 <- "Others"
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8[setdiff(GMP_Subsisters_PI_Cluster3_idx,
                                                  GMP_Subsisters_PI_Cluster8_idx)] <- "GMP_End_In_Cluster3"
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8[setdiff(GMP_Subsisters_PI_Cluster8_idx,
                                                  GMP_Subsisters_PI_Cluster3_idx)] <- "GMP_End_In_Cluster8"
  
  ### UMAP
  p <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
               group.by = "GMP_End_In_PI_Cluster3_8",
               pt.size = 2, raster = TRUE,
               cols = c("GMP_End_In_Cluster3" = "red", "GMP_End_In_Cluster8" = "blue"),
               order = c("GMP_End_In_Cluster8", "GMP_End_In_Cluster3", "Others")) +
    ggtitle("") +
    labs(color="GMP Effector Precursors") +
    theme_classic(base_size = 36) +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.9
  ggsave(paste0(outputDir, "UMAP_Effector_End_In_Cluster3_8.png"), plot = p, width = 18, height = 10, dpi = 350)
  
  
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
      gsea_result <- data.frame(fgseaMultilevel(pathways = gene_list, stats = signature[[1]], minSize = -Inf, maxSize = Inf))
    }
    
    ### print GSEA plot
    sIdx <- which(gsea_result$padj < fdr_cutoff)
    if(printPlot && length(sIdx) > 0) {
      for(i in sIdx) {
        ### get required values ready
        if(length(signature) > 1) {
          geneset <- gene_list[[i]]
          stats <- signature[[i]]
          stats <- stats[order(-stats)]
          fileName <- names(signature)[i]
        } else {
          geneset <- gene_list[[gsea_result$pathway[i]]]
          stats <- signature[[1]]
          stats <- stats[order(-stats)]
          fileName <- gsea_result$pathway[i]
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
  
  
  ###
  ### 15. GSEA on precursor TF regulons - with signature from cluster3 vs cluster8 (or GMP end up in cluster 3 vs in cluster 8)
  ###
  
  ### get pyscenic regulon file
  f <- "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_regulons.p"
  
  ### load regulon files and load them into R
  pd <- import("pandas")
  
  ### load the *.p pickle data
  pickle_data <- pd$read_pickle(f)
  names(pickle_data) <- sapply(pickle_data, function(x) x$transcription_factor)
  
  ### make an empty regulon list
  regulon_list <- vector("list", length(pickle_data))
  names(regulon_list) <- names(pickle_data)
  
  ### save target genes to the list
  for(tf in names(pickle_data)) {
    regulon_list[[tf]] <- unlist(pickle_data[[tf]]$genes)
  }
  
  ### garbage collection
  gc()
  
  
  ### PI cluster3 vs cluster8
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
  de_result_pi_cluster3_8 <- FindMarkers(JCC_Seurat_Obj,
                                         ident.1 = "PI_Cluster3",
                                         ident.2 = "PI_Cluster8",
                                         min.pct = 0.2,
                                         logfc.threshold = 0.1,
                                         test.use = "wilcox")
  
  ### GMP precursor end in cluster3 vs cluster8
  ### set idents
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8)
  
  ### GMP precursor end in cluster3 vs cluster8
  de_result_precursor_cluster3_8 <- FindMarkers(JCC_Seurat_Obj,
                                                ident.1 = "GMP_End_In_Cluster3",
                                                ident.2 = "GMP_End_In_Cluster8",
                                                min.pct = 0.2,
                                                logfc.threshold = 0.1,
                                                test.use = "wilcox")
  
  ### PI cluster3 vs cluster8
  ### run GSEA
  signat <- de_result_pi_cluster3_8$avg_log2FC
  names(signat) <- rownames(de_result_pi_cluster3_8)
  GSEA_result <- run_gsea(gene_list = regulon_list, signature = list(signat),
                          fdr_cutoff = 0.05,
                          printPlot = TRUE, printPath = paste0(outputDir, "/GSEA/pi_cluster3_8/"))
  GSEA_result <- GSEA_result[order(GSEA_result$padj),]
  
  ### write out the result file
  write.xlsx2(GSEA_result, file = paste0(outputDir, "/GSEA/pi_cluster3_8/GSEA_Precursor_Regulons_pi_cluster3_8.xlsx"),
              sheetName = "GSEA_Precursor_Regulons_pi_cluster3_8", row.names = FALSE)
  
  ### GMP precursor end in cluster3 vs cluster8
  ### run GSEA
  signat <- de_result_precursor_cluster3_8$avg_log2FC
  names(signat) <- rownames(de_result_precursor_cluster3_8)
  GSEA_result <- run_gsea(gene_list = regulon_list, signature = list(signat),
                          fdr_cutoff = 0.05,
                          printPlot = TRUE, printPath = paste0(outputDir, "/GSEA/precursor_cluster3_8/"))
  GSEA_result <- GSEA_result[order(GSEA_result$padj),]
  
  ### write out the result file
  write.xlsx2(GSEA_result, file = paste0(outputDir, "/GSEA/precursor_cluster3_8/GSEA_Precursor_Regulons_precursor_cluster3_8.xlsx"),
              sheetName = "GSEA_Precursor_Regulons_precursor_cluster3_8", row.names = FALSE)
  
  
  ###
  ### 16. Recluster GMP without proliferating clusters - then effector precursors might be clustered together
  ###
  
  ### NEW FUNCTIONAL ANNOTATION BY TAY - 09/27/2021
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters <- "Others"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "1", "5", "7", "9", "10", "11", "12", "15", "19", "21"))] <- "Proliferating"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("2", "4", "6", "17"))] <- "Transitioning"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("3", "8", "14"))] <- "Functional Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("13", "20"))] <- "Dysfunctional"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("16"))] <- "Early Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("18"))] <- "Metabolically Active"
  
  ### get the subset
  GMP_no_prolf_obj <- subset(JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data)[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                                  which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters != "Proliferating"))])
  
  ### reclustering
  ### normalization
  GMP_no_prolf_obj <- NormalizeData(GMP_no_prolf_obj,
                                    normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  GMP_no_prolf_obj <- FindVariableFeatures(GMP_no_prolf_obj,
                                           selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  GMP_no_prolf_obj <- ScaleData(GMP_no_prolf_obj,
                                vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  GMP_no_prolf_obj <- RunPCA(GMP_no_prolf_obj,
                             features = VariableFeatures(object = GMP_no_prolf_obj),
                             npcs = 30)
  ElbowPlot(GMP_no_prolf_obj, ndims = 30, reduction = "pca")
  GMP_no_prolf_obj <- RunUMAP(GMP_no_prolf_obj, dims = 1:30)
  
  ### perform clustering
  GMP_no_prolf_obj <- FindNeighbors(GMP_no_prolf_obj, dims = 1:30)
  GMP_no_prolf_obj <- FindClusters(GMP_no_prolf_obj, resolution = 0.7)
  
  ### UMAP - new clustering
  p <- DimPlot(object = GMP_no_prolf_obj, reduction = "umap",
               group.by = "RNA_snn_res.0.2",
               pt.size = 0.5, raster = TRUE,
               label = TRUE, label.size = 10) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  ggsave(paste0(outputDir, "UMAP_GMP_No_Proliferating_Clusters.png"), plot = p, width = 18, height = 10, dpi = 350)
  
  ### UMAP - effectors
  p <- DimPlot(object = GMP_no_prolf_obj, reduction = "umap",
               group.by = "GMP_End_In_PI_Cluster3_8",
               pt.size = 2, raster = TRUE,
               cols = c("GMP_End_In_Cluster3" = "red", "GMP_End_In_Cluster8" = "blue", "Others" = "lightgray"),
               order = c("GMP_End_In_Cluster8", "GMP_End_In_Cluster3", "Others")) +
    ggtitle("") +
    labs(color="Effector Precursors") +
    theme_classic(base_size = 36) +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  ggsave(paste0(outputDir, "UMAP_GMP_No_Proliferating_Effectors.png"), plot = p, width = 18, height = 10, dpi = 350)
  
  ### UMAP - original cluster
  p <- DimPlot(object = GMP_no_prolf_obj, reduction = "umap",
               group.by = "AllSeuratClusters",
               pt.size = 0.5, raster = TRUE,
               label = TRUE, label.size = 10) +
    ggtitle("") +
    labs(color="Original Clusters") +
    theme_classic(base_size = 36) +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  ggsave(paste0(outputDir, "UMAP_GMP_No_Proliferating_Original_Clusters.png"), plot = p, width = 18, height = 10, dpi = 350)
  
  ### UMAP - functional annotation
  p <- DimPlot(object = GMP_no_prolf_obj, reduction = "umap",
               group.by = "New_Functional_Annotation_Based_On_Clusters",
               pt.size = 1, raster = TRUE) +
    ggtitle("") +
    labs(color="Functions") +
    theme_classic(base_size = 36) +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  ggsave(paste0(outputDir, "UMAP_GMP_No_Proliferating_Functions.png"), plot = p, width = 18, height = 10, dpi = 350)
  
  
  ###
  ### 17. DE analysis for PI CD8s:
  #         - Compare PI effectors vs all other PIs
  #         - Compare PI effectors vs each CD8 dysfunctional/dying cluster
  #         - Compare each dysfunctional/dying cluster to each other
  ###
  
  ### CD4/CD8 annotation
  JCC_Seurat_Obj$CD4_CD8_by_Clusters <- "NA"
  JCC_Seurat_Obj$CD4_CD8_by_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "2", "9", "10", "11", "14", "15", "18"))] <- "CD4"
  JCC_Seurat_Obj$CD4_CD8_by_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("1", "3", "5", "6", "7", "8", "12", "13", "16", "17", "19", "20"))] <- "CD8"
  
  ### - Compare PI_CD8_Effectors vs Other_PI_CD8s
  ### set comparison
  JCC_Seurat_Obj$PI_Effectors <- NA
  JCC_Seurat_Obj$PI_Effectors[intersect(which(JCC_Seurat_Obj$AllSeuratClusters %in% c("3", "8")),
                                        intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                  which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8")))] <- "PI_CD8_Effectors"
  JCC_Seurat_Obj$PI_Effectors[setdiff(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8")),
                                      which(JCC_Seurat_Obj$PI_Effectors == "PI_CD8_Effectors"))] <- "Other_PI_CD8s"
  
  ### set idents
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$PI_Effectors)
  
  ### PI_CD8_Effectors vs Other_PI_CD8s
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "PI_CD8_Effectors",
                           ident.2 = "Other_PI_CD8s",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/DE_PI_CD8_Effectors_vs_Others.xlsx"),
              sheetName = "DE_PI_CD8_Effectors_vs_Others", row.names = FALSE)
  
  ### - Compare PI effectors vs each CD8 dysfunctional/dying cluster
  ### set comparison
  JCC_Seurat_Obj$PI_Effectors2 <- JCC_Seurat_Obj$PI_Effectors
  JCC_Seurat_Obj$PI_Effectors2[intersect(which(JCC_Seurat_Obj$PI_Effectors2 == "Other_PI_CD8s"),
                                         which(JCC_Seurat_Obj$AllSeuratClusters == "13"))] <- "PI_CD8_Exhausted_Effectors"
  JCC_Seurat_Obj$PI_Effectors2[intersect(which(JCC_Seurat_Obj$PI_Effectors2 == "Other_PI_CD8s"),
                                         which(JCC_Seurat_Obj$AllSeuratClusters == "20"))] <- "PI_CD8_Dying"
  
  ### set idents
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$PI_Effectors2)
  
  ### PI_CD8_Effectors vs PI_CD8_Exhausted_Effectors
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "PI_CD8_Effectors",
                           ident.2 = "PI_CD8_Exhausted_Effectors",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/DE_PI_CD8_Effectors_vs_Exhausted_Effectors.xlsx"),
              sheetName = "DE_PI_CD8_Effectors_vs_Exhausted_Effectors", row.names = FALSE)
  
  ### PI_CD8_Effectors vs PI_CD8_Dying
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "PI_CD8_Effectors",
                           ident.2 = "PI_CD8_Dying",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/DE_PI_CD8_Effectors_vs_Dying.xlsx"),
              sheetName = "DE_PI_CD8_Effectors_vs_Dying", row.names = FALSE)
  
  ### - Compare each dysfunctional/dying cluster to each other
  ### PI_CD8_Exhausted_Effectors vs PI_CD8_Dying
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "PI_CD8_Exhausted_Effectors",
                           ident.2 = "PI_CD8_Dying",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/DE_PI_CD8_Exhausted_Effectors_vs_Dying.xlsx"),
              sheetName = "DE_PI_CD8_Exhausted_Effectors_vs_Dying", row.names = FALSE)
  
  
  ### effector cluster3 vs exhausting group cluster 13
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$AllSeuratClusters)
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "3",
                           ident.2 = "13",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/DE_Cluster3_vs_Cluster13.xlsx"),
              sheetName = "DE_Cluster3_vs_Cluster13", row.names = FALSE)
  
  ### effector cluster8 vs exhausting group cluster 13
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj$AllSeuratClusters)
  de_result <- FindMarkers(JCC_Seurat_Obj,
                           ident.1 = "8",
                           ident.2 = "13",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/DE_Cluster8_vs_Cluster13.xlsx"),
              sheetName = "DE_Cluster8_vs_Cluster13", row.names = FALSE)
  
  
  
  ###
  ### 18. Ridge plot of Tay's genes comparing cluster 3, 8, 13, & 20
  ###     CD7, CD2, CD52, CD300A, KLRD1 (CD94), FCGR3A, CX3CR1, KLRC1 (NKG2A)
  ###     KIR2DL3/CD158b2, KLRG1, CD82, CXCR6, CXCR3, CD27, TIGIT
  ###     TNFRSF9 (CD137/4-1BB), CD84, CD82, CD28, ITGA4 (CD49d)
  ###
  
  ### get the target seurat object
  cluster3_8_13_20_seurat <- subset(JCC_Seurat_Obj,
                                    cells = rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("3", "8", "13", "20"))])
  
  ### set target genes
  target_genes <- c("CD2", "CD7", "CD27", "CD52", "CD82", "CD84", "CD300A", "KLRD1", "FCGR3A", "CX3CR1", "KLRC1",
                    "KIR2DL3", "KLRG1", "CXCR3", "CXCR6", "TIGIT", "TNFRSF9", "ITGA4")
  
  ### Ridge plot
  cluster3_8_13_20_seurat <- SetIdent(object = cluster3_8_13_20_seurat,
                                      cells = rownames(cluster3_8_13_20_seurat@meta.data),
                                      value = cluster3_8_13_20_seurat$AllSeuratClusters)

  p <- RidgePlot(cluster3_8_13_20_seurat, features = target_genes)
  for(i in 1:length(target_genes)) {
    p[[i]] <- p[[i]] + labs(y = "Clusters")
  }
  ggsave(paste0(outputDir, "RidgePlot_Tay_Genes_Cluster_3_8_13_20.png"), plot = p, width = 25, height = 15, dpi = 350)
  
  
  
  
  
}