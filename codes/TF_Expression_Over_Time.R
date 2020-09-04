###
#   File name : TF_Expression_Over_Time.R
#   Author    : Hyunjin Kim
#   Date      : Sep 2, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Run 'FindAllMarkers()' for every time point and calculate combined p-value of the DE genes
#               using the Fisher's method. Then see if there are TFs of our interest (Steven's list) in the
#               top-ranking DE genes. Draw a plot with the interesting TFs to show how their expression changes
#               over time.
#
#   Instruction
#               1. Source("TF_Expression_Over_Time.R")
#               2. Run the function "tf_expression" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_TF_Expression_Over_Time.R/TF_Expression_Over_Time.R")
#               > tf_expression(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                               outputDir="./results/PROTO/DEEP/TF_Expression/")
###

tf_expression <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                          outputDir="./results/PROTO/DEEP/TF_Expression/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
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
  
  ### active assay = "RNA"
  Seurat_Obj@active.assay <- "RNA"
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  
  ### set Idents of the object before performing DE analysis
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$Time)
  
  ### I tried to run FindAllMarkers() but the sample size is too large,
  ### so I divide the samples for each time point then randomly chose from the others to compare
  
  ### for each time point run DE analysis vs the rest
  tp_of_interest <- setdiff(levels(Seurat_Obj@meta.data$TimeF), c("PreTrans", "Wk-1", "Wk0"))
  for(tp in tp_of_interest) {
    
    ### choose random samples from the rest
    set.seed(1234)
    ctrl_grp_idx <- sample(which(Seurat_Obj@meta.data$Time != tp),
                           length(which(Seurat_Obj@meta.data$Time == tp)))
    
    ###
    de_result <- FindMarkers()
    
    
  }
  
  
  
  
  ### run DE analysis based on all the time points
  all_de_result <- FindAllMarkers(object = Seurat_Obj, test.use = "DESeq2",
                                  logfc.threshold = 0, min.pct = 0.1)
  
  
  
  
  
  
}
