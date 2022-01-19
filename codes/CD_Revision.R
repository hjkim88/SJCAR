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
#               3. R3C5
#                  PBMC and BMMC samples are treated identically. It would be interesting to understand how much cells
#                  from these different samples are driving some of the differences in the clusters.
#               4. R4C7 & R4C8
#                  expression of CASP8 was highest in state B relative to all other states (Fig 3C)". Is it significant?
#                  higher relative expression of TOX and LAG3. Is it significant? - State C
#
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
  
  ### add previous monocle_state to the current one
  monocle_cds2$Original_State <- "NEW"
  monocle_cds2@phenoData@data[rownames(monocle_cds@phenoData@data),"Original_State"] <- monocle_cds$State
  
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
  
  ### color scale
  sjcar19_colors <- c("#640B11", "#AA4C26", "#D39F3A", "#C09969", "#287B66", "#487A8F", "#3B3B53")
  names(sjcar19_colors) <- levels(monocle_cds2$State)
  show_col(sjcar19_colors)
  
  ### draw monocle plot for checking
  plot_cell_trajectory(monocle_cds2, color_by = "State",
                       cell_size = 3, cell_link_size = 3,
                       show_branch_points = FALSE) +
    labs(color="State") +
    scale_color_manual(values = sjcar19_colors, name = "State") +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    theme_classic(base_size = 30) +
    theme(legend.position = "top",
          axis.title = element_text(size = 30, color = "black", face = "bold"),
          axis.text = element_text(size = 25, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  
  
  
}