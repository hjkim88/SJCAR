###
#   File name : Apply_Scenic.R
#   Author    : Hyunjin Kim
#   Date      : Feb 2, 2022
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Run Scenic on SJCAR19 data
#
#   Instruction
#               1. Source("Apply_Scenic.R")
#               2. Run the function "run_scenic" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Apply_Scenic.R/Apply_Scenic.R")
#               > run_scenic(Seurat_RObj_path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds",
#                            outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Manuscript/Scenic/")
###

run_scenic <- function(Seurat_RObj_path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds",
                       outputDir="./results/New3/Manuscript/Scenic/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(AUCell, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("AUCell")
    require(AUCell, quietly = TRUE)
  }
  if(!require(RcisTarget, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("RcisTarget")
    require(RcisTarget, quietly = TRUE)
  }
  if(!require(GENIE3, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("GENIE3")
    require(GENIE3, quietly = TRUE)
  }
  if(!require(devtools, quietly = TRUE)) {
    install.packages("devtools")
    require(devtools, quietly = TRUE)
  }
  if(!require(SCENIC, quietly = TRUE)) {
    devtools::install_github("aertslab/SCENIC") 
    require(SCENIC, quietly = TRUE)
  }
  if(!require(doRNG, quietly = TRUE)) {
    install.packages("doRNG")
    require(doRNG, quietly = TRUE)
  }
  if(!require(SeuratDisk, quietly = TRUE)) {
    remotes::install_github("mojaveazure/seurat-disk")
    require(SeuratDisk, quietly = TRUE)
  }
  
  
  ### load Jeremy's object
  JCC_Seurat_Obj <- readRDS(file = Seurat_RObj_path)
  
  ### check whether the orders are the same
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
  print(identical(names(Idents(object = JCC_Seurat_Obj)), rownames(JCC_Seurat_Obj@meta.data)))
  
  ### normalization
  JCC_Seurat_Obj <- NormalizeData(JCC_Seurat_Obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### Cell cycle score (will be used later for regression out)
  JCC_Seurat_Obj <- CellCycleScoring(object = JCC_Seurat_Obj,
                                     g2m.features = cc.genes$g2m.genes,
                                     s.features = cc.genes$s.genes)
  
  ### find variable genes
  JCC_Seurat_Obj <- FindVariableFeatures(JCC_Seurat_Obj,
                                         selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  JCC_Seurat_Obj <- ScaleData(JCC_Seurat_Obj,
                              vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
  
  ### save rds
  # saveRDS(JCC_Seurat_Obj, file = "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC2.rds")
  
  ### download motif ranking db
  ### this is one time thing
  ### but this R download makes error in initializeScenic(), so just download the files directly
  ### mc9nr: Motif collection version 9: 24k motifs
  # dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
  #              "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
  # for(featherURL in dbFiles) {
  #   download.file(featherURL, destfile=paste0("./data/", basename(featherURL)))
  # }
  
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
  
  ### maybe we need to select cells for specific purposes
  precursor_seurat <- subset(JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")])
  non_precursor_seurat <- subset(JCC_Seurat_Obj,
                                 cells = rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 == "Other_CD8_GMPs")])
  
  ### prepare Scenic input
  exprMat <- as.matrix(precursor_seurat@assays$RNA@data)
  cellInfo <- precursor_seurat@meta.data
  
  ### randomly select GMP non-precursor group
  target_px <- unique(JCC_Seurat_Obj@meta.data[colnames(exprMat),"px"])
  cell_num_per_px <- round(ncol(exprMat) / length(target_px))
  set.seed(1234)
  random_idx <- NULL
  for(px in target_px) {
    random_idx <- c(random_idx, sample(which(non_precursor_seurat$px == px), cell_num_per_px))
  }
  exprMat2 <- as.matrix(non_precursor_seurat@assays$RNA@data[,random_idx])
  cellInfo2 <- non_precursor_seurat@meta.data[random_idx,]
  non_precursor_seurat_subset <- subset(non_precursor_seurat,
                                        cells = rownames(non_precursor_seurat@meta.data)[random_idx])
  
  ### combine precursor & non-precursor subset seurat
  combined_seurat <- merge(x = precursor_seurat, y = non_precursor_seurat_subset,
                           add.cell.ids = c("precursor", "non-precursor"))
  
  ### save seurat object as loom
  SaveLoom(precursor_seurat, filename = "./data/precursor_seurat.loom")
  SaveLoom(non_precursor_seurat_subset, filename = "./data/non_precursor_seurat_subset.loom")
  SaveLoom(combined_seurat, filename = "./data/combined_seurat.loom")  
  
  
  ### now we want to see cluster3+8 specifically
  ### but wanna add some other cluster cells as control
  cluster3_8_cells <- NULL
  for(clstr in unique(JCC_Seurat_Obj$AllSeuratClusters)) {
    if(clstr %in% c("3", "8")) {
      cluster3_8_cells <- c(cluster3_8_cells, sample(rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$AllSeuratClusters == clstr)], 1000))
    } else {
      cluster3_8_cells <- c(cluster3_8_cells, sample(rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$AllSeuratClusters == clstr)], 100))
    }
  }
  cluster3_8_seurat <- subset(JCC_Seurat_Obj,
                              cells = cluster3_8_cells)
  
  ### save seurat object as loom
  SaveLoom(cluster3_8_seurat, filename = "./data/cluster3_8_seurat.loom")
  
  
  ### run this at the R terminal
  ### bash: python3: command not found
  ### python --version
  ### Python 3.7.6
  ### pip install loompy
  ### pip install arboreto
  ### pip install pyscenic
  ###
  ### running on my desktop
  # python /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/arboreto_with_multiprocessing.py \
  # /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_seurat.loom \
  # /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/hs_hgnc_tfs.txt \
  # --method grnboost2 \
  # --output /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_adj.tsv \
  # --num_workers 20 \
  # --seed 777
  
  
  
  # ### scenic option
  # ### org <- mgi or hgnc, or dmel
  # ### mgi stands for mouse, hgnc stands for human, dmel stands for fly
  # ### dbDir <- "cisTarget_databases" # RcisTarget databases location
  # scenicOptions <- initializeScenic(org="hgnc", dbDir="./data/cisTarget_databases",
  #                                   dbs = defaultDbNames[["hgnc"]], nCores=1)
  # 
  # ### run Genie3 to make correlations between TF and target genes
  # genesKept <- geneFiltering(exprMat, scenicOptions)
  # exprMat_filtered <- exprMat[genesKept, ]
  # runCorrelation(exprMat_filtered, scenicOptions)
  # exprMat_filtered <- log2(exprMat_filtered+1)
  # scenicOptions@settings$seed <- 123
  # runGenie3(exprMat_filtered, scenicOptions)
  # 
  # ### Build the gene regulatory network: 1. Get co-expression modules 2. Get regulons (with RcisTarget): TF motif analysis)
  # ### Identify cell states: 3. Score GRN (regulons) in the cells (with AUCell) 4. Cluster cells according to the GRN activity
  # scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  # scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod = c("top5perTarget"))
  # scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
  # 
  # saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status
  
}