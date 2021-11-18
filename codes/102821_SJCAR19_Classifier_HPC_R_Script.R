### Example script for SJCAR19 classifier
### Oct 28, 2021
### hyunjin.kim@stjude.org

### batch number
batch_num <- 1

### load libraries
if(!require(Seurat, quietly = TRUE)) {
  install.packages("Seurat")
  require(Seurat, quietly = TRUE)
}
if(!require(caret, quietly = TRUE)) {
  install.packages("caret")
  require(caret, quietly = TRUE)
}
if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
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

### load Jeremy's object
JCC_Seurat_Obj <- readRDS(file = "/research/sharedresources/immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds")

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

### change the class names
JCC_Seurat_Obj$Classification_Class <- JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8
JCC_Seurat_Obj$Classification_Class[which(JCC_Seurat_Obj$Classification_Class == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "YES"
JCC_Seurat_Obj$Classification_Class[which(JCC_Seurat_Obj$Classification_Class == "Other_CD8_GMPs")] <- "NO"

### create output path
outputDir3 <- paste0("/research/sharedresources/immunoinformatics/hkim8/SJCAR19_classifier/", batch_num, "/")
dir.create(outputDir3, recursive = TRUE)

### parameter setting for a classifier
iteration <- 10
set.seed(batch_num)
featureSelectionNum <- 100
sampleNum <- 100
methodTypes <- c("svmRadial")
methodNames <- c("SVMRadial")
train_control <- trainControl(method="none", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)

### build the classifier 1000 times with random down-sampling
acc <- vector("list", length = iteration)
roc <- vector("list", length = iteration)
for(i in 1:iteration) {
  
  writeLines(paste(i, "/", iteration))
  
  ### get downsampling indicies
  downsample_idx <- c(sample(which(JCC_Seurat_Obj$Classification_Class == "YES"), sampleNum),
                      sample(which(JCC_Seurat_Obj$Classification_Class == "NO"), sampleNum))
  
  ### get normalized data
  input_data <- normalizeRNASEQwithVST(readCount = data.frame(JCC_Seurat_Obj@assays$RNA@counts[,downsample_idx],
                                                              stringsAsFactors = FALSE, check.names = FALSE))
  
  ### manual LOOCV
  acc[[i]] <- rep(NA, length(downsample_idx))
  roc[[i]] <- rep(NA, length(downsample_idx))
  for(j in 1:length(downsample_idx)) {
    
    ### divide training and test index
    test_idx <- downsample_idx[j]
    training_idx <- downsample_idx[-j]
    
    ### new obj for the training data & set idents with the info
    classifier_seurat_obj <- subset(JCC_Seurat_Obj, cells = rownames(JCC_Seurat_Obj@meta.data)[training_idx])
    classifier_seurat_obj <- SetIdent(object = classifier_seurat_obj,
                                      cells = rownames(classifier_seurat_obj@meta.data),
                                      value = classifier_seurat_obj$Classification_Class)
    
    ### DE analysis
    de_result <- FindMarkers(classifier_seurat_obj,
                             ident.1 = "YES",
                             ident.2 = "NO",
                             min.pct = 0.1,
                             logfc.threshold = 0.1,
                             test.use = "wilcox",
                             verbose = FALSE)
    de_result <- de_result[order(de_result$p_val_adj),]
    de_genes <- rownames(de_result)[1:featureSelectionNum]
    
    ### annotate class for the input data
    input_data2 <- data.frame(t(input_data[intersect(de_genes, rownames(input_data)),rownames(classifier_seurat_obj@meta.data)]), stringsAsFactors = FALSE, check.names = FALSE)
    input_data2$Class <- factor(classifier_seurat_obj$Classification_Class,
                                levels = c("YES", "NO"))
    
    ### annotate class for the test data
    test_data <- data.frame(t(input_data[intersect(de_genes, rownames(input_data)),rownames(JCC_Seurat_Obj@meta.data)[test_idx],drop=FALSE]), stringsAsFactors = FALSE, check.names = FALSE)
    test_data$Class <- factor(JCC_Seurat_Obj$Classification_Class[downsample_idx[j]],
                              levels = c("YES", "NO"))
    
    ### get accuracy and pred prob (not AUC yet)
    model <- train(Class~., data=input_data2, method=methodTypes, trControl = train_control)
    pred_result <- predict(model, newdata = test_data)
    acc[[i]][j] <- round((sum(pred_result == test_data$Class) / nrow(test_data)), 3)
    pred_result <- predict(model, newdata = test_data, type = "prob")
    roc[[i]][j] <- pred_result$YES
    
  }
  
  ### calculate the AUCs
  roc[[i]] <- roc(JCC_Seurat_Obj$Classification_Class[downsample_idx], roc[[i]])
  
  ### draw ROC curve
  pdf(paste0(outputDir3, "Classifier_GMP_Precursor_vs_Other_CD8_GMPs_(", i, ")_Batch", batch_num, ".pdf"),
      width = 5, height = 4)
  plot.roc(roc[[i]], main = paste("Accuracy =", mean(as.numeric(acc[[i]]))),
           legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
           xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
  dev.off()
  
}

### save acc and roc objects
saveRDS(acc, file = paste0(outputDir3, "Classifier_ACC_Batch", batch_num, ".rds"))
saveRDS(roc, file = paste0(outputDir3, "Classifier_ROC_Batch", batch_num, ".rds"))
