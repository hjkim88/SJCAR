###
#   File name : CD_Revision.R
#   Author    : Hyunjin Kim
#   Date      : Jan 11, 2022
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. a) DE genes between cluster 3 & 8
#                  b) DE genes between GMP that ends up with cluster 3 vs GMP that ends up with cluster 8
#               2. a) DE genes between CAR > 0 VS CAR = 0
#                  b) DE genes between CAR > 0; HIGH vs CAR > 0; low
#               3. PBMC and BMMC samples are treated identically. It would be interesting to understand how much cells
#                  from these different samples are driving some of the differences in the clusters.
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
  
  ### see distribution of the CAR expressions
  plot(density(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",]),
       xlim = c(0,50))
  plot(density(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",]),
       xlim = c(0,10))
  
  ### 2 maybe a good threshold
  
  
  ### set new column for GMP cells that end up in PI cluster 3 and in PI cluster 8
  ### but no duplicates this time
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8 <- "Others"
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8[setdiff(GMP_Subsisters_PI_Cluster3_idx,
                                                  GMP_Subsisters_PI_Cluster8_idx)] <- "GMP_End_In_Cluster3"
  JCC_Seurat_Obj$GMP_End_In_PI_Cluster3_8[setdiff(GMP_Subsisters_PI_Cluster8_idx,
                                                  GMP_Subsisters_PI_Cluster3_idx)] <- "GMP_End_In_Cluster8"
  
  
  
}