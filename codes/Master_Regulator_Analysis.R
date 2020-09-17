###
#   File name : Master_Regulator_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : Sep 11, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Perform SCENIC to find TF (master regulator) candidates in each time point.
#
#   * SCENIC set-up: https://rawcdn.githack.com/aertslab/SCENIC/701cc7cc4ac762b91479b3bd2eaf5ad5661dd8c2/inst/doc/SCENIC_Setup.html
#   * SCENIC running example: https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/SCENIC_Running.html
#
#   Instruction
#               1. Source("Master_Regulator_Analysis.R")
#               2. Run the function "master_regulator_analysis" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Master_Regulator_Analysis.R/Master_Regulator_Analysis.R")
#               > master_regulator_analysis(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                                           outputDir="./results/PROTO/DEEP/Master_Regulator/")
###

master_regulator_analysis <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                                      outputDir="./results/PROTO/DEEP/Master_Regulator/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(SingleCellExperiment, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("SingleCellExperiment")
    require(SingleCellExperiment, quietly = TRUE)
  }
  if(!require(GENIE3, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("GENIE3")
    require(GENIE3, quietly = TRUE)
  }
  if(!require(RcisTarget, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("RcisTarget")
    require(RcisTarget, quietly = TRUE)
  }
  if(!require(AUCell, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("AUCell")
    require(AUCell, quietly = TRUE)
  }
  if(!require(devtools, quietly = TRUE)) {
    install.packages("devtools")
    require(devtools, quietly = TRUE)
  }
  if(!require(SCENIC, quietly = TRUE)) {
    install_github("aertslab/SCENIC")
    require(SCENIC, quietly = TRUE)
  }
  if(!require(doRNG, quietly = TRUE)) {
    install.packages("doRNG")
    require(doRNG, quietly = TRUE)
  }
  
  ### SCENIC parameter setting
  scenicOptions <- initializeScenic(org="hgnc", dbDir="data", nCores=10)
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
  
  ### active assay = "RNA"
  Seurat_Obj@active.assay <- "RNA"
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  # ### human-specific (hg19) databases for RcisTarget (the motif rankings)
  # dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
  #              "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
  # ### mouse-specific (mm9) databases for RcisTarget (the motif rankings)
  # dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
  #              "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
  
  # for(featherURL in dbFiles) {
  #   download.file(featherURL,
  #                 paste0("./data/", destfile=basename(featherURL)))
  # }
  
  ### set the ident of the object with the Cell type
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$Time)
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  ### run SCENIC for each time point
  time_points <- levels(Seurat_Obj@meta.data$TimeF)
  sampleNum <- 200
  umi_thresh <- 100
  for(tp in time_points) {
    
    ### subset the Seurat object by time
    subset_Seurat_Obj <- subset(Seurat_Obj, idents=tp)
    
    ### filter the cells (a cell should have at least 5 genes expressed larger than 0)
    ### randomly select one and then check if that satisfies the threshold
    set.seed(1234)
    sample_list <- NULL
    pool <- 1:ncol(subset_Seurat_Obj@assays$RNA@counts)
    while(length(sample_list) < sampleNum) {
      random_idx <- sample(pool, 1)
      if(length(which(subset_Seurat_Obj@assays$RNA@counts[,random_idx] > 0)) > umi_thresh) {
        sample_list <- c(sample_list, random_idx)
        pool <- setdiff(pool, sample_list)
      }
    }
    
    ### prepare input data for SCENIC
    exprMat <- as.matrix(subset_Seurat_Obj@assays$RNA@counts[,sample_list])
    cellInfo <- data.frame(seuratCluster=Idents(subset_Seurat_Obj),
                           stringsAsFactors = FALSE, check.names = FALSE)[sample_list,,drop=FALSE]
    
    ### Co-expression network
    genesKept <- geneFiltering(exprMat, scenicOptions)
    exprMat_filtered <- exprMat[genesKept, ]
    runCorrelation(exprMat_filtered, scenicOptions)
    exprMat_filtered_log <- log2(exprMat_filtered+1) 
    runGenie3(exprMat_filtered_log, scenicOptions)
    
    ### Build and score the GRN
    exprMat_log <- log2(exprMat+1)
    scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
    runSCENIC_1_coexNetwork2modules(scenicOptions)
    runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
    runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
    
    
    
    
  }
  
  
  
  
}
