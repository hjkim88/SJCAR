###
#   File name : CD4_CD8_Ratio.R
#   Author    : Hyunjin Kim
#   Date      : Jun 4, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : With CD4, CD8A, and CD8B expression, we determine whether a given cell is a CD4 cell
#               or a CD8 cell. Then we summarize how many CD4 & CD8 cells are in each library.
#
#   Instruction
#               1. Source("CD4_CD8_Ratio.R")
#               2. Run the function "get_cd4_cd8_ratio" - specify the input Seurat object and the output directory
#               3. The ratio summary table will be generated under the output directory
#
#   Example
#               > source("The_directory_of_CD4_CD8_Ratio.R/CD4_CD8_Ratio.R")
#               > get_cd4_cd8_ratio(Seurat_Obj=Seurat_Obj,
#                                   outputDir="./results/")
###

get_cd4_cd8_ratio <- function(Seurat_Obj=Seurat_Obj,
                              outputDir="./results/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  
  ### get gene expressions of CD4/CD8A/CD8B
  cd4_exps <- Seurat_Obj@assays$RNA@counts["CD4",]
  cd8a_exps <- Seurat_Obj@assays$RNA@counts["CD8A",]
  
  ### plot the density
  plot(density(cd4_exps[which(Seurat_Obj@meta.data$CAR == "CARpos")]), main = "CAR+ CD4 Expression Density")
  plot(density(cd8a_exps[which(Seurat_Obj@meta.data$CAR == "CARpos")]), main = "CAR+ CD8A Expression Density")
  plot(cd4_exps[which(Seurat_Obj@meta.data$CAR == "CARpos")],
       cd8a_exps[which(Seurat_Obj@meta.data$CAR == "CARpos")],
       xlab = "CD4 Expression", ylab = "CD8A Expression",
       main = "CAR+ CD4-CD8A Expressions", pch = 16)
  
  ### CD4 & CD8 table
  result_table_p <- data.frame(matrix(0, length(unique(Seurat_Obj@meta.data$Library)), 10),
                             stringsAsFactors = FALSE, check.names = FALSE)
  rownames(result_table_p) <- unique(Seurat_Obj@meta.data$Library)
  colnames(result_table_p) <- c("Patient", "Time",
                              "All_GEX_CD4_Ratio", "All_GEX_CD8_Ratio",
                              "GEX_TCR_CD4_Ratio", "GEX_TCR_CD8_Ratio",
                              "CAR+_All_GEX_CD4_Ratio", "CAR+_All_GEX_CD8_Ratio",
                              "CAR+_GEX_TCR_CD4_Ratio", "CAR+_GEX_TCR_CD8_Ratio")
  
  ### fill out the patient and the time columns
  libs <- strsplit(unique(Seurat_Obj@meta.data$Library), split = "_", fixed = TRUE)
  result_table_p[,"Patient"] <- sapply(libs, function(x) x[2])
  result_table_p[,"Time"] <- sapply(libs, function(x) x[3])
  
  ### we also need the real numbers
  result_table_n <- result_table_p
  colnames(result_table_n) <- c("Patient", "Time",
                                "All_GEX_CD4_Number", "All_GEX_CD8_Number",
                                "GEX_TCR_CD4_Number", "GEX_TCR_CD8_Number",
                                "CAR+_All_GEX_CD4_Number", "CAR+_All_GEX_CD8_Number",
                                "CAR+_GEX_TCR_CD4_Number", "CAR+_GEX_TCR_CD8_Number")
  
  ### fill out the ratios
  for(i in 1:nrow(result_table_p)) {
    all_cell_num <- length(which(Seurat_Obj@meta.data$Library == rownames(result_table_p)[i]))
    gex_tcr_cell_num <- length(intersect(which(!is.na(Seurat_Obj@meta.data$cdr3_nt)),
                                         which(Seurat_Obj@meta.data$Library == rownames(result_table_p)[i])))
    carpos_all_cell_num <- length(intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                            which(Seurat_Obj@meta.data$Library == rownames(result_table_p)[i])))
    carpos_gex_tcr_cell_num <- length(intersect(intersect(which(!is.na(Seurat_Obj@meta.data$cdr3_nt)),
                                                          which(Seurat_Obj@meta.data$CAR == "CARpos")),
                                                which(Seurat_Obj@meta.data$Library == rownames(result_table_p)[i])))
    
    ### CD4
    target_idx <- intersect(which(Seurat_Obj@meta.data$Library == rownames(result_table_p)[i]),
                            which(cd4_exps > 0))
    result_table_p[i,"All_GEX_CD4_Ratio"] <- length(target_idx) / all_cell_num
    result_table_n[i,"All_GEX_CD4_Number"] <- length(target_idx)
    if(gex_tcr_cell_num > 0) {
      result_table_p[i,"GEX_TCR_CD4_Ratio"] <- length(intersect(target_idx,
                                                              which(!is.na(Seurat_Obj@meta.data$cdr3_nt)))) / gex_tcr_cell_num
      result_table_n[i,"GEX_TCR_CD4_Number"] <- length(intersect(target_idx,
                                                                which(!is.na(Seurat_Obj@meta.data$cdr3_nt))))
    }
    target_idx <- intersect(target_idx,
                            which(Seurat_Obj@meta.data$CAR == "CARpos"))
    if(carpos_all_cell_num > 0) {
      result_table_p[i,"CAR+_All_GEX_CD4_Ratio"] <- length(target_idx) / carpos_all_cell_num
      result_table_n[i,"CAR+_All_GEX_CD4_Number"] <- length(target_idx)
    }
    if(carpos_gex_tcr_cell_num > 0) {
      result_table_p[i,"CAR+_GEX_TCR_CD4_Ratio"] <- length(intersect(target_idx,
                                                              which(!is.na(Seurat_Obj@meta.data$cdr3_nt)))) / carpos_gex_tcr_cell_num
      result_table_n[i,"CAR+_GEX_TCR_CD4_Number"] <- length(intersect(target_idx,
                                                                     which(!is.na(Seurat_Obj@meta.data$cdr3_nt))))
    }
    
    ### CD8
    target_idx <- intersect(which(Seurat_Obj@meta.data$Library == rownames(result_table_p)[i]),
                            which(cd8a_exps > 1))
    result_table_p[i,"All_GEX_CD8_Ratio"] <- length(target_idx) / all_cell_num
    result_table_n[i,"All_GEX_CD8_Number"] <- length(target_idx)
    if(gex_tcr_cell_num > 0) {
      result_table_p[i,"GEX_TCR_CD8_Ratio"] <- length(intersect(target_idx,
                                                              which(!is.na(Seurat_Obj@meta.data$cdr3_nt)))) / gex_tcr_cell_num
      result_table_n[i,"GEX_TCR_CD8_Number"] <- length(intersect(target_idx,
                                                                which(!is.na(Seurat_Obj@meta.data$cdr3_nt))))
    }
    target_idx <- intersect(target_idx,
                            which(Seurat_Obj@meta.data$CAR == "CARpos"))
    if(carpos_all_cell_num > 0) {
      result_table_p[i,"CAR+_All_GEX_CD8_Ratio"] <- length(target_idx) / carpos_all_cell_num
      result_table_n[i,"CAR+_All_GEX_CD8_Number"] <- length(target_idx)
    }
    if(carpos_gex_tcr_cell_num > 0) {
      result_table_p[i,"CAR+_GEX_TCR_CD8_Ratio"] <- length(intersect(target_idx,
                                                                   which(!is.na(Seurat_Obj@meta.data$cdr3_nt)))) / carpos_gex_tcr_cell_num
      result_table_n[i,"CAR+_GEX_TCR_CD8_Number"] <- length(intersect(target_idx,
                                                                     which(!is.na(Seurat_Obj@meta.data$cdr3_nt))))
    }
  }
  
  ### x 100 to compute the percentage
  result_table_p[,c("All_GEX_CD4_Ratio", "All_GEX_CD8_Ratio",
                  "GEX_TCR_CD4_Ratio", "GEX_TCR_CD8_Ratio",
                  "CAR+_All_GEX_CD4_Ratio", "CAR+_All_GEX_CD8_Ratio",
                  "CAR+_GEX_TCR_CD4_Ratio", "CAR+_GEX_TCR_CD8_Ratio")] <- result_table_p[,c("All_GEX_CD4_Ratio", "All_GEX_CD8_Ratio",
                                                                                          "GEX_TCR_CD4_Ratio", "GEX_TCR_CD8_Ratio",
                                                                                          "CAR+_All_GEX_CD4_Ratio", "CAR+_All_GEX_CD8_Ratio",
                                                                                          "CAR+_GEX_TCR_CD4_Ratio", "CAR+_GEX_TCR_CD8_Ratio")] * 100
  
  ### round the numbers
  result_table_p[,c("All_GEX_CD4_Ratio", "All_GEX_CD8_Ratio",
                  "GEX_TCR_CD4_Ratio", "GEX_TCR_CD8_Ratio",
                  "CAR+_All_GEX_CD4_Ratio", "CAR+_All_GEX_CD8_Ratio",
                  "CAR+_GEX_TCR_CD4_Ratio", "CAR+_GEX_TCR_CD8_Ratio")] <- signif(result_table_p[,c("All_GEX_CD4_Ratio", "All_GEX_CD8_Ratio",
                                                                                                 "GEX_TCR_CD4_Ratio", "GEX_TCR_CD8_Ratio",
                                                                                                 "CAR+_All_GEX_CD4_Ratio", "CAR+_All_GEX_CD8_Ratio",
                                                                                                 "CAR+_GEX_TCR_CD4_Ratio", "CAR+_GEX_TCR_CD8_Ratio")], digits = 4)
  
  ### only keep GMP info
  # result_table_p <- result_table_p[which(result_table_p$Time == "GMP"),]
  # result_table_n <- result_table_n[which(result_table_n$Time == "GMP"),]
  
  ### only keep GMP+ time points
  result_table_p <- result_table_p[which(result_table_p$Time %in% setdiff(unique(result_table_p$Time),
                                                                          c("PreTrans",
                                                                            "PreTransB",
                                                                            "Wk-1",
                                                                            "Wk-1Run1",
                                                                            "Wk-1Run2",
                                                                            "Wk0"))),]
  result_table_n <- result_table_n[which(result_table_n$Time %in% setdiff(unique(result_table_n$Time),
                                                                          c("PreTrans",
                                                                            "PreTransB",
                                                                            "Wk-1",
                                                                            "Wk-1Run1",
                                                                            "Wk-1Run2",
                                                                            "Wk0"))),]
  
  ### reconstruct the ratio table with the numbers only
  remove_idx <- NULL
  for(i in 1:nrow(result_table_p)) {
    sum_i <- result_table_n[i,3] + result_table_n[i,4]
    if(sum_i == 0) {
      remove_idx <- c(remove_idx, i)
    } else {
      result_table_p[i,3] <- signif(result_table_n[i,3] * 100 / sum_i, digits = 4)
      result_table_p[i,4] <- signif(result_table_n[i,4] * 100 / sum_i, digits = 4)
    }
  }
  if(length(remove_idx) > 0) {
    result_table_p <- result_table_p[-remove_idx,]
    result_table_n <- result_table_n[-remove_idx,]
  }
  
  ### only keep CAR+ info
  result_table_p <- result_table_p[,-c(3,4,5,6,9,10)]
  result_table_n <- result_table_n[,-c(3,4,5,6,9,10)]
  
  ### get average ratios
  new_table <- result_table_n[-which(result_table_n$Time == "GMP"),]
  new_result_table <- data.frame(matrix(0, length(unique(new_table$Patient)), 3),
                                        stringsAsFactors = FALSE, check.names = FALSE)
  rownames(new_result_table) <- unique(new_table$Patient)
  colnames(new_result_table) <- c("Patient", "CAR+_All_GEX_CD4_Ratio", "CAR+_All_GEX_CD8_Ratio")
  new_result_table$Patient <- rownames(new_result_table)
  for(px in unique(new_table$Patient)) {
    cd4_sum <- sum(new_table[which(new_table$Patient == px),3])
    cd8_sum <- sum(new_table[which(new_table$Patient == px),4])
    all_sum <- cd4_sum + cd8_sum
    new_result_table[px,2] <- signif(cd4_sum * 100 / all_sum, digits = 4)
    new_result_table[px,3] <- signif(cd8_sum * 100 / all_sum, digits = 4)
  }
  
  ### write out the result
  write.xlsx2(result_table_p, file = paste0(outputDir, "CD4_CD8_Ratio.xlsx"),
              sheetName = "CD4_CD8_Ratio", row.names = FALSE)
  write.xlsx2(result_table_n, file = paste0(outputDir, "CD4_CD8_Ratio.xlsx"),
              sheetName = "CD4_CD8_Numbers", row.names = FALSE, append = TRUE)
  write.xlsx2(new_result_table, file = paste0(outputDir, "CD4_CD8_Ratio.xlsx"),
              sheetName = "CD4_CD8_Average_Ratio", row.names = FALSE, append = TRUE)
  
}
