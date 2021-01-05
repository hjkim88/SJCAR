###
#   File name : GMP_Persistency_Marker_Discovery.R
#   Author    : Hyunjin Kim
#   Date      : Dec 29, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. What markers predict whether a lineage (CAR+/CAR-) will persist after infusion?
#               - DE analysis
#               - GO/Pathway analysis
#               - A heatmap to show those genes are actually differentially expressed between the two conditions
#               - Comparison of various factors between persisters & non-persisters (violin plot?)
#               - Using the DE genes to make a classifier
#               - To ensure an individual does not drive the found pattern?
#                 a) if the classifier (LOOCV) worked well, it is a proof
#                 b) multiple 1d scatter plots (colored based on patients) to show that the top DE genes are not one-patient specific
#               - Classifier comparison between 1. using the top DE genes vs 2. using random genes
#               - Revered PCA/UMAP (gene-patient conversion) to show whether those genes are separated by persistency
#
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
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(caret, quietly = TRUE)) {
    install.packages("caret")
    require(caret, quietly = TRUE)
  }
  if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
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
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title)) +
                theme(axis.text = element_text(size = 25))
              
              png(paste0(dir, "kegg_", title, ".png"), width = 2000, height = 1000)
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
                ggtitle(paste0("GO ", title)) +
                theme(axis.text = element_text(size = 25))
              
              png(paste0(dir, "go_", title, ".png"), width = 2000, height = 1000)
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
  
  ### DE analysis
  de_result <- FindMarkers(Seurat_Obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/GMP_CARpos_Persisters_vs_NonPersisters.xlsx"),
              sheetName = "GMP_CARpos_DE_Result", row.names = FALSE)
  
  ### pathway analysis
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                            rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                            "ENTREZID", "SYMBOL"),
                                          org = "human", database = "GO",
                                          title = paste0("Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters"),
                                          displayNum = 30, imgPrint = TRUE,
                                          dir = paste0(outputDir))
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters"),
                                            displayNum = 30, imgPrint = TRUE,
                                            dir = paste0(outputDir))
  if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
    write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters.xlsx"),
                row.names = FALSE, sheetName = paste0("GO_Result"))
  }
  if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
    write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters.xlsx"),
                row.names = FALSE, sheetName = paste0("KEGG_Result"))
  }
  
  #
  ### make a classifier with the DE genes
  #
  ### create outputDir2
  outputDir2 <- paste0(outputDir, "DE_Classifier/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### because of the imbalance of the two cluster sizes, we randomly choose
  ### the same number of samples in each class and iteratively build the classifier
  iteration <- 10
  set.seed(1234)
  featureSelectionNum <- 100
  sampleNum <- 300
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "rf", "LogitBoost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "LogitBoost", "K-NN")
  train_control <- trainControl(method="LOOCV", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  
  ### a small number that will be added to the counts for normalization
  ### log transform should not have 0 values
  log_trans_add <- 1
  
  ### build the classifier 10 times
  for(i in 1:iteration) {
    
    ### randomly choose the "gmp not last" samples
    random_gmp_not_last <- sample(which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO"),
                                  size = length(all_gmp_last))
    
    ### normalize the read counts
    ### before the normalization, only keep the samples that will be used in the classifier
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],
                                                                                             c(sample(all_gmp_last, sampleNum),
                                                                                               sample(random_gmp_not_last, sampleNum))] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE))
    
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
    png(paste0(outputDir2, "Classifier_Using_DEG_GMP_Last_vs_Not_Last_", featureSelectionNum, "_(", i, ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", acc[j]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
    
  }
  
  
  
  
  # ### get GMP CAR+ cells as a SingleCellExperiment object
  # GMP_CARpos_sce <- as.SingleCellExperiment(subset(Seurat_Obj, idents = c("YES", "NO")), assay = "RNA")
  # 
  # ### for memory efficiency, remove the large Seurat object
  # rm(Seurat_Obj)
  # gc()
  # 
  # ### QC
  # ### 26027 genes & 109614 cells
  # row_sum <- rowSums(counts(GMP_CARpos_sce) > 0)
  # col_sum <- colSums(counts(GMP_CARpos_sce) > 0)
  # plot(density(row_sum))
  # plot(density(col_sum))
  # 
  # ### set parameter
  # ### gene.min.pct: in a gene, at least n % of each condition (cells) should be expressed
  # ### cell.min.num: in a cell, at least n genes should be expressed (determine after looking at the density plot above)
  # gene.min.pct <- 0.2
  # cell.min.num <- 2000
  # cell.max.num <- 5000
  # 
  # ### low count filter for genes
  # yes_num <- length(which(GMP_CARpos_sce$ident == "YES"))
  # no_num <- length(which(GMP_CARpos_sce$ident == "NO"))
  # yes_sum <- rowSums(counts(GMP_CARpos_sce)[,which(GMP_CARpos_sce$ident == "YES")] > 0)
  # no_sum <- rowSums(counts(GMP_CARpos_sce)[,which(GMP_CARpos_sce$ident == "NO")] > 0)
  # yes_exp_pct <- yes_sum / yes_num
  # no_exp_pct <- no_sum / no_num
  # keep <- intersect(which(yes_exp_pct >= gene.min.pct),
  #                   which(no_exp_pct >= gene.min.pct))
  # GMP_CARpos_sce <- GMP_CARpos_sce[keep,]
  # 
  # ### low count filter for cells
  # keep <- intersect(which(colSums(counts(GMP_CARpos_sce) > 0) >= cell.min.num),
  #                   which(colSums(counts(GMP_CARpos_sce) > 0) <= cell.max.num))
  # GMP_CARpos_sce <- GMP_CARpos_sce[,keep]
  # 
  # 
  # ### there are too many cells (especially in the second condition), so we will randomly choose some
  # ### to lower the computational complexity
  # set.seed(1234)
  # keep <- sample(which(GMP_CARpos_sce$ident == "NO"), length(which(GMP_CARpos_sce$ident == "YES")))
  # keep <- union(keep, which(GMP_CARpos_sce$ident == "YES"))
  # GMP_CARpos_sce <- GMP_CARpos_sce[,keep]
  # 
  # ### "counts" should be in the first place in the assayNames(GMP_CARpos_sce)
  # nms <- c("counts", setdiff(assayNames(GMP_CARpos_sce), "counts"))
  # assays(GMP_CARpos_sce) <- assays(GMP_CARpos_sce)[nms]
  # 
  # ### epsilon setting as recommended by the ZINB-WaVE integration paper
  # system.time({
  #   zinb <- zinbwave(GMP_CARpos_sce, K=0, observationalWeights=TRUE,
  #                    BPPARAM=SerialParam(), epsilon=1e12)
  # })
  # 
  # ### DE analysis
  # dds <- DESeqDataSet(zinb, design=~condition)
  # dds <- estimateSizeFactors(dds, type="poscounts")
  # library(scran)
  # scr <- computeSumFactors(dds)
  # dat <- data.frame(true=dds$ExpLibSize,
  #                   pos=sizeFactors(dds),
  #                   sum=sizeFactors(scr))
  # dat$true <- dat$true / exp(mean(log(dat$true)))
  # panel.scatter <- function(x,y,...) {
  #   points(x,y,...)
  #   abline(0,1,col="red",lwd=2)
  #   legend("topleft", legend=round(cor(x,y),3))
  # }
  # pairs(dat, panel=panel.scatter)
  
  
  
  
}
