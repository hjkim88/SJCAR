###
#   File name : Classifier.R
#   Author    : Hyunjin Kim
#   Date      : Jun 15, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Building some interesting classifier with the SJCAR19 data to classify such as:
#               1. GMP CAR+ cells that will last or not
#               2. Predict time point of how long a given GMP CAR+ cell will last
#               3. Responder vs Non-responder
#               4. Cytokine levels
#
#   Instruction
#               1. Source("Classifier.R")
#               2. Run the function "classifier_run" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Classifier.R/Classifier.R")
#               > classifier_run(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                                outputDir="./results/PROTO/DEEP/")
###

classifier_run <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                           outputDir="./results/PROTO/DEEP/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(caret, quietly = TRUE)) {
    install.packages("caret")
    require(caret, quietly = TRUE)
  }
  if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
  }
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  
  ### set new result directory
  outputDir2 <- paste0(outputDir, "Classifier/")
  dir.create(path = outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### global clonotypes in the Seurat objects
  global_clonotypes <- colnames(Seurat_Obj@meta.data)[grep("global_clonotype", colnames(Seurat_Obj@meta.data), fixed = TRUE)]
  
  ### get patient ids (dir names) from the result directory
  f <- list.dirs(paste0(outputDir, "../"), full.names = FALSE, recursive = FALSE)
  f <- f[grep("SJCAR19", f)]
  
  
  ### 1. GMP CAR+ cells that will last or not
  
  ### set new result directory
  outputDir3 <- paste0(outputDir2, "GMP_Last_vs_GMP_Not_Last/")
  dir.create(path = outputDir3, showWarnings = FALSE, recursive = TRUE)
  
  ### get "GMP CAR+ last" idx
  ### for each patient that has TCR info
  type <- global_clonotypes[1]
  all_gmp_last <- NULL
  for(px in f) {
    ### print progress
    writeLines(paste(px))
    
    ### load the file
    target_file <- read.xlsx2(file = paste0(outputDir, "../", px, "/car_clonotype_frequency_over_time_", px, ".xlsx"),
                              sheetName = type, stringsAsFactors = FALSE, check.names = FALSE,
                              row.names = 1)
    
    ### numerize the table
    for(i in 1:ncol(target_file)) {
      target_file[,i] <- as.numeric(target_file[,i])
    }
    
    ### remove all zero time points
    time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
    
    ### get time points after GMP infusion
    time_points <- setdiff(time_points,
                           c("PreTrans", "Wk-1", "Wk0", "Total"))
    
    if(length(time_points) > 1 && length(which(time_points == "GMP")) > 0) {
      ### select clonotypes
      tidx <- which(time_points == "GMP")
      target_temp <- target_file[which(target_file[,"GMP"] > 0),
                                 time_points[(tidx+1):length(time_points)],drop=FALSE]
      target_clonotypes <- rownames(target_temp)[which(apply(target_temp, 1, sum) > 0)]
      
      if(length(target_clonotypes) > 0) {
        px_time_idx <- intersect(intersect(which(Seurat_Obj@meta.data$Px == px),
                                           which(Seurat_Obj@meta.data$Time == time_points[tidx])),
                                 which(Seurat_Obj@meta.data$CAR == "CARpos"))
        all_gmp_last <- c(all_gmp_last, intersect(px_time_idx,
                                                  which(Seurat_Obj@meta.data[,type] %in% target_clonotypes)))
      }
    }
  }
  
  ### get "GMP CAR+ do not last" idx
  all_gmp_not_last <- setdiff(setdiff(intersect(intersect(which(Seurat_Obj@meta.data$Time == "GMP"),
                                                          which(Seurat_Obj@meta.data$CAR == "CARpos")),
                                                which(!is.na(Seurat_Obj@meta.data$cdr3_nt))),
                                      all_gmp_last),
                              which(Seurat_Obj@meta.data$Px %in% c("SJCAR19-00", "SJCAR19-01")))
  
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
  
  ### a function to select genes based on variance
  selectTopV <- function(x, selectNum) {
    v <- apply(x, 1, var)
    x <- x[order(-v),]
    x <- x[1:selectNum,]
    
    return (x)
  }
  
  ### because of the imbalance of the two cluster sizes, we randomly choose
  ### the same number of samples in each class and iteratively build the classifier
  iteration <- 10
  set.seed(1234)
  featureSelectionNum <- 100
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "rf", "LogitBoost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "LogitBoost", "K-NN")
  train_control <- trainControl(method="LOOCV", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
  for(i in 1:iteration) {
    
    ### randomly choose the "gmp not last" samples
    random_gmp_not_last <- sample(all_gmp_not_last, size = length(all_gmp_last))
    
    ### normalize the read counts
    ### before the normalization, only keep the samples that will be used in the classifier
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[,c(all_gmp_last,
                                                                                                random_gmp_not_last)],
                                                                stringsAsFactors = FALSE, check.names = FALSE))
    
    ### reduce the gene size based on variance
    ### only select high variance genes
    input_data <- selectTopV(input_data, featureSelectionNum)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(c(rep("GMP_Last", length(all_gmp_last)),
                                 rep("GMP_Not_Last", length(random_gmp_not_last))), levels = c("GMP_Last", "GMP_Not_Last"))
    
    ### build classifier and test
    ### LOOCV
    p <- list()
    acc <- NULL
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=input_data, trControl=train_control, method=methodTypes[j])
      roc <- roc(model$pred$obs, model$pred$GMP_Last)
      acc <- c(acc, round(mean(model$results$Accuracy), 3))
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using Gene Expressions\n",
                                           "Accuracy =", acc[j]),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      gc()
    }
    
    ### draw ROC curves
    png(paste0(outputDir3, "Classifier_GMP_Last_vs_Not_Last_", featureSelectionNum, "_(", i, ").png"),
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
  
  
  ### 2. Predict time point of how long a given GMP CAR+ cell will last
  
  ### set new result directory
  outputDir3 <- paste0(outputDir2, "GMP_Time_Predictor/")
  dir.create(path = outputDir3, showWarnings = FALSE, recursive = TRUE)
  
  
  
  
  
  
  
  
}
