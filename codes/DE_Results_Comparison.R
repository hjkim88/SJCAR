###
#   File name : DE_Results_Comparison.R
#   Author    : Hyunjin Kim
#   Date      : Jun 15, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Compare the two DE results based on logFC and on the FDR
#               1. GMP CAR+ cells last vs do not last
#               2. Responder vs Non-responder
#
#   Instruction
#               1. Source("DE_Results_Comparison.R")
#               2. Run the function "de_results_comparison" - specify the input file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DE_Results_Comparison.R/DE_Results_Comparison.R")
#               > de_results_comparison(de_result1_path="./results/PROTO/DEEP/DE_Results/ALL_Patients_DE_Genes_DESeq2.xlsx",
#                                       de_result2_path="./results/PROTO/DEEP/DE_Results/JCC212_Aggreg21Feb_GMP_CARPOS_RespVSNonMarkers.tsv",
#                                       outputDir="./results/PROTO/DEEP/DE_Results/")
###

de_results_comparison <- function(de_result1_path="./results/PROTO/DEEP/DE_Results/ALL_Patients_DE_Genes_DESeq2.xlsx",
                                  de_result2_path="./results/PROTO/DEEP/DE_Results/JCC212_Aggreg21Feb_GMP_CARPOS_RespVSNonMarkers.tsv",
                                  outputDir="./results/PROTO/DEEP/DE_Results/") {
  
  ### load libraries
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(ggrepel, quietly = TRUE)) {
    install.packages("ggrepel")
    require(ggrepel, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  
  ### load DE results
  de_result1 <- read.xlsx2(file = de_result1_path, sheetIndex = 1,
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(de_result1) <- de_result1$Gene_Symbol
  de_result2 <- read.table(file = de_result2_path, header = TRUE, sep = "\t", row.names = 1,
                           stringsAsFactors = FALSE, check.names = FALSE)
  
  ### get shared gene symbols between the two
  ### In the DE result2, there are identical gene symbols and we use any first ones
  shared_genes <- intersect(rownames(de_result1), rownames(de_result2))
  
  ### combine two DE results
  combined_result <- data.frame(Gene_Symbol=shared_genes,
                                logFC1=de_result1[shared_genes,"avg_logFC"],
                                logFC2=de_result2[shared_genes,"avg_logFC"],
                                adj_p1=de_result1[shared_genes,"FDR"],
                                adj_p2=de_result2[shared_genes,"p_val_adj"],
                                Label="",
                                stringsAsFactors = FALSE, check.names = FALSE)
  combined_result <- combined_result[which(combined_result$adj_p1 != ""),]
  combined_result[2:5] <- sapply(combined_result[2:5], as.numeric)
  
  ### mark 10 interesting genes
  interesting_num <- 20
  # interesting_genes <- c("RPS26", "CD40LG", "TIMP1", "PPDPF", "YBX1", "RPS2", "HLA-DRB1")
  same_sign_gene_idx <- which(sign(combined_result$logFC1) == sign(combined_result$logFC2))
  abs_logFC_sum <- abs(combined_result$logFC1) + abs(combined_result$logFC2)
  top_logFC_idx <- order(-abs_logFC_sum)
  interesting_gene_idx <- intersect(top_logFC_idx, same_sign_gene_idx)[1:interesting_num]
  combined_result$Label[interesting_gene_idx] <- combined_result$Gene_Symbol[interesting_gene_idx]
  
  ### comparison based on logFC
  ggplot(data = combined_result, aes(x=logFC1, y=logFC2)) +
    geom_point(color = "black", size = 1) +
    geom_label_repel(aes(logFC1, logFC2, label = Label), color = "red", box.padding = unit(0.45, "lines")) +
    labs(title="logFC Comparison",
         subtitle=sprintf("S.Cor = %s, p-value = %s",
                          round(cor(combined_result$logFC1, combined_result$logFC2,
                                    use = "pairwise.complete.obs", method = "spearman"), 5),
                          signif(cor.test(combined_result$logFC1, combined_result$logFC2,
                                          method = "spearman")$p.value, 5))
    ) +
    xlab("logFC of (GMP CAR+ Last vs NOT)") +
    ylab("logFC of (Responder vs Non-responder)") +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(file = paste0(outputDir, "DE_logFC_Comparison.png"), width = 10, height = 8, dpi = 300)
  
  ### add new metric
  combined_result$metric1 <- sign(combined_result$logFC1) * log(combined_result$adj_p1)
  combined_result$metric2 <- sign(combined_result$logFC2) * log(combined_result$adj_p2)
  
  ### comparison based on sign(logFC) * log(p-value)
  ggplot(data = combined_result, aes(x=metric1, y=metric2)) +
    geom_point(color = "black", size = 1) +
    geom_label_repel(aes(metric1, metric2, label = Label), color = "red", box.padding = unit(0.45, "lines")) +
    labs(title="sign(logFC) * log(p-value) Comparison",
         subtitle=sprintf("S.Cor = %s, p-value = %s",
                          round(cor(combined_result$metric1, combined_result$metric2,
                                    use = "pairwise.complete.obs", method = "spearman"), 5),
                          signif(cor.test(combined_result$metric1, combined_result$metric2,
                                          method = "spearman")$p.value, 5))
    ) +
    xlab("(GMP CAR+ Last vs NOT)") +
    ylab("(Responder vs Non-responder)") +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(file = paste0(outputDir, "DE_sign(logFC)_x_log(p-value)_Comparison.png"), width = 10, height = 8, dpi = 300)
  
  
  ### pathway analysis
  
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
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title))
              
              png(paste0(dir, "kegg_", title, "_CB.png"), width = 2000, height = 1000)
              print(p[[1]])
              dev.off()
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
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("GO ", title))
              
              png(paste0(dir, "go_", title, "_CB.png"), width = 2000, height = 1000)
              print(p[[2]])
              dev.off()
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
  
  ### get max & min values
  metric1_max <- max(combined_result$metric1[which(combined_result$metric1 != Inf)], na.rm = TRUE)
  metric2_max <- max(combined_result$metric2[which(combined_result$metric2 != Inf)], na.rm = TRUE)
  metric1_min <- min(combined_result$metric1[which(combined_result$metric1 != -Inf)], na.rm = TRUE)
  metric2_min <- min(combined_result$metric2[which(combined_result$metric2 != -Inf)], na.rm = TRUE)
  
  ### Inf -> max+1, -Inf -> min-1
  combined_result$metric1[which(combined_result$metric1 == Inf)] <- metric1_max+1
  combined_result$metric2[which(combined_result$metric2 == Inf)] <- metric2_max+1
  combined_result$metric1[which(combined_result$metric1 == -Inf)] <- metric1_min-1
  combined_result$metric2[which(combined_result$metric2 == -Inf)] <- metric2_min-1
  
  ### scaled new metric
  temp <- scale(combined_result[,c("metric1", "metric2")])
  colnames(temp) <- c("scaled_metric1", "scaled_metric2")
  combined_result <- data.frame(combined_result, temp,
                                stringsAsFactors = FALSE, check.names = FALSE)
  abs_metric_sum <- abs(combined_result$scaled_metric1) + abs(combined_result$scaled_metric2)
  top_metric_idx <- order(-abs_metric_sum)
  input_gene_idx <- intersect(top_metric_idx, same_sign_gene_idx)[1:interesting_num]
  combined_result$Label <- NA
  combined_result$Label[input_gene_idx] <- combined_result$Gene_Symbol[input_gene_idx]
  
  ### comparison based on sign(logFC) * log(p-value)
  ggplot(data = combined_result, aes(x=scaled_metric1, y=scaled_metric2)) +
    geom_point(color = "black", size = 1) +
    geom_label_repel(aes(scaled_metric1, scaled_metric2, label = Label), color = "red", box.padding = unit(0.45, "lines")) +
    labs(title="sign(logFC) * log(p-value) Comparison - Scaled",
         subtitle=sprintf("S.Cor = %s, p-value = %s",
                          round(cor(combined_result$scaled_metric1, combined_result$scaled_metric2,
                                    use = "pairwise.complete.obs", method = "spearman"), 5),
                          signif(cor.test(combined_result$scaled_metric1, combined_result$scaled_metric2,
                                          method = "spearman")$p.value, 5))
    ) +
    xlab("(GMP CAR+ Last vs NOT)") +
    ylab("(Responder vs Non-responder)") +
    geom_smooth(method = lm, color="gray", se=FALSE) +
    theme_classic(base_size = 16)
  ggsave(file = paste0(outputDir, "DE_sign(logFC)_x_log(p-value)_Comparison_Scaled.png"), width = 10, height = 8, dpi = 300)
  
  
  ### get pathway genes
  geneList <- mapIds(org.Hs.eg.db, combined_result$Gene_Symbol[input_gene_idx], "ENTREZID", "SYMBOL")
  
  ### pathway analysis
  pathway_results_GO <- pathwayAnalysis_CP(geneList = geneList,
                                           org = "human", database = "GO",
                                           title = paste0("Pathway_Results_", "Interesting_Genes"),
                                           displayNum = 50, imgPrint = TRUE,
                                           dir = paste0(outputDir))
  pathway_results_KEGG <- pathwayAnalysis_CP(geneList = geneList,
                                             org = "human", database = "KEGG",
                                             title = paste0("Pathway_Results_Interesting_Genes"),
                                             displayNum = 50, imgPrint = TRUE,
                                             dir = paste0(outputDir))
  
  ### save the pathway results in excel file
  write.xlsx2(pathway_results_GO,
              file = paste0(outputDir, "GO_pathway_results_", "Interesting_Genes", ".xlsx"),
              row.names = FALSE, sheetName = "GO_Interesting_Genes")
  write.xlsx2(pathway_results_KEGG,
              file = paste0(outputDir, "KEGG_pathway_results_", "Interesting_Genes", ".xlsx"),
              row.names = FALSE, sheetName = "KEGG_Interesting_Genes")
  
}
