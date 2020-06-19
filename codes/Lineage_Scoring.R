###
#   File name : Lineage_Scoring.R
#   Author    : Hyunjin Kim
#   Date      : Jun 17, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Give p-values to the GMP CAR+ clones that represent how significant a given CAR+ clone size is
#               at each time point in terms of GMP lasting capability.
#               * A matrix of rows: All GMP CAR+ Clones
#                           & cols: All the GMP+ time points + Total p-value
#               * One-tailed Fisher's exact test + Fisher's method
#
#   Instruction
#               1. Source("Lineage_Scoring.R")
#               2. Run the function "give_lineage_score" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Lineage_Scoring.R/Lineage_Scoring.R")
#               > give_lineage_score(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                                    outputDir="./results/PROTO/")
###

give_lineage_score <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                               outputDir="./results/PROTO/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
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
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
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
  
  ### get patient ids (dir names) from the result directory
  f <- list.dirs(outputDir, full.names = FALSE, recursive = FALSE)
  f <- f[grep("SJCAR19", f)]
  
  ### for each patient, perform the analysis
  for(px in f) {
    ### print progress
    writeLines(paste(px))
    
    ### Patient-CAR+ idx
    target_idx <- intersect(which(Seurat_Obj@meta.data$Px == px),
                            which(Seurat_Obj@meta.data$CAR == "CARpos"))
    
    ### remove NA, "NA", and ""
    target_idx <- setdiff(target_idx,
                          union(union(which(is.na(Seurat_Obj@meta.data$cdr3_nt)),
                                      which(Seurat_Obj@meta.data$cdr3_nt == "NA")),
                                which(Seurat_Obj@meta.data$cdr3_nt == "")))
    
    ### get unique GMP CAR+ clones
    unique_gmp_carpos_clones <- unique(Seurat_Obj@meta.data$global_clonotype_ab_strict0[intersect(target_idx,
                                                                                                  which(Seurat_Obj@meta.data$Time == "GMP"))])
    
    ### get CAR+ time points
    time_points <- unique(Seurat_Obj@meta.data$Time[target_idx])
    
    if(length(unique_gmp_carpos_clones) > 0 && length(time_points) > 1) {
      ### get CDR3 sequences of the unique GMP CAR+ clones
      unique_cdr3_seqs <- unique(Seurat_Obj@meta.data$cdr3_nt[intersect(target_idx,
                                                                        which(Seurat_Obj@meta.data$Time == "GMP"))])
      
      ### make two empty matrices (a: count, b:p-value from FET, c:p-value from OR)
      c_mat <- matrix(0, nrow = length(unique_gmp_carpos_clones), ncol = length(time_points)+1)
      rownames(c_mat) <- unique_gmp_carpos_clones
      colnames(c_mat) <- c(time_points, "Total")
      c_mat <- data.frame(Clone_ID=unique_gmp_carpos_clones,
                          CDR3_NT=unique_cdr3_seqs,
                          c_mat,
                          stringsAsFactors = FALSE, check.names = FALSE)
      p_mat <- matrix(1, nrow = length(unique_gmp_carpos_clones), ncol = length(time_points)+1)
      rownames(p_mat) <- unique_gmp_carpos_clones
      colnames(p_mat) <- c(time_points, "Total")
      p_mat <- data.frame(Clone_ID=unique_gmp_carpos_clones,
                          CDR3_NT=unique_cdr3_seqs,
                          p_mat,
                          stringsAsFactors = FALSE, check.names = FALSE)
      o_mat <- matrix(0, nrow = length(unique_gmp_carpos_clones), ncol = length(time_points)+1)
      rownames(o_mat) <- unique_gmp_carpos_clones
      colnames(o_mat) <- c(time_points, "Total")
      o_mat <- data.frame(Clone_ID=unique_gmp_carpos_clones,
                          CDR3_NT=unique_cdr3_seqs,
                          o_mat,
                          stringsAsFactors = FALSE, check.names = FALSE)
      
      ### get p-values for each clone and for each time point
      for(clone in unique_gmp_carpos_clones) {
        for(tp in time_points) {
          ### get clone size
          c_mat[clone, tp] <- length(intersect(intersect(target_idx,
                                                         which(Seurat_Obj@meta.data$Time == "GMP")),
                                               which(Seurat_Obj@meta.data$global_clonotype_ab_strict0 == clone)))
          
          ### calculate p-value
          ### Fisher's exact test
          ###
          ###           TP CAR+   No TP (GMP) CAR+
          ###          ----------------------------
          ###    Clone |   X            Y
          ### No Clone |   Z            W
          X <- c_mat[clone, tp]
          Y <- length(intersect(intersect(target_idx,
                                          which(Seurat_Obj@meta.data$Time == tp)),
                                which(Seurat_Obj@meta.data$global_clonotype_ab_strict0 == clone)))
          Z <- length(intersect(intersect(target_idx,
                                          which(Seurat_Obj@meta.data$Time == tp)),
                                which(Seurat_Obj@meta.data$global_clonotype_ab_strict0 != clone)))
          W <- length(intersect(intersect(target_idx,
                                          which(Seurat_Obj@meta.data$Time == "GMP")),
                                which(Seurat_Obj@meta.data$global_clonotype_ab_strict0 != clone)))
          p_mat[clone, tp] <- fisher.test(matrix(c(X, Z, Y, W), 2, 2), alternative = "greater")$p.value
          
          ### Odds ratio
          o_mat[clone, tp] <- fisher.test(matrix(c(X, Z, Y, W), 2, 2), alternative = "greater")$estimate
        }
      }
      
      ### fill out the total column
      c_mat$Total <- apply(c_mat[,setdiff(time_points, "GMP")], 1, function(x) sum(x, na.rm = TRUE))
      p_mat$Total <- fisher.method(pvals = p_mat[,setdiff(time_points, "GMP")],
                                   p.corr = "BH",
                                   na.rm = TRUE)
      o_mat$Total <- apply(o_mat[,setdiff(time_points, "GMP")], 1, function(x) mean(x, na.rm = TRUE))
      
      ### order the tables
      c_mat <- c_mat[order(-c_mat$Total),]
      p_mat <- p_mat[order(p_mat$Total),]
      o_mat <- o_mat[order(-o_mat$Total),]
      
      ### write the results
      write.xlsx2(c_mat,
                  file = paste0(outputDir, f, "/GMP_CAR+_Lineage_Statistics_", f, ".xlsx"),
                  sheetName = "Count",
                  append = FALSE)
      write.xlsx2(o_mat,
                  file = paste0(outputDir, f, "/GMP_CAR+_Lineage_Statistics_", f, ".xlsx"),
                  sheetName = "Odds_Ratio",
                  append = TRUE)
      write.xlsx2(p_mat,
                  file = paste0(outputDir, f, "/GMP_CAR+_Lineage_Statistics_", f, ".xlsx"),
                  sheetName = "P-value",
                  append = TRUE)
    }
  }
  
}
