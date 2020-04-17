### If we just load tcr data from the gmp and the wk1, then compare the cdr3 nt sequences regardless of a&b chains,
### there are many commonly shared sequences. I would like to make a table of them to show why those do not show
### lineage if we combine into one barcode with considering alpha and beta chains.

### load the tcr data
gmp_tcr <- read.csv(file = "C:/Users/hkim8/SJ/SJCAR19-05/SJCAR19-05_TCR/1757073_JCC212_SJCAR19-05_GMP19016_TCR/filtered_contig_annotations.csv",
                    stringsAsFactors = FALSE)
wk1_tcr <- read.csv(file = "C:/Users/hkim8/SJ/SJCAR19-05/SJCAR19-05_TCR/1757057_JCC212_SJCAR19-05_Wk1_TCR/filtered_contig_annotations.csv",
                    stringsAsFactors = FALSE)

### load the combined data
load("./data/JCC212_Px5_TCR_clonotyped2.Robj")

### get cdr3 nt sequences of the gmp that also appear in the wk1
target_seqs <- unique(gmp_tcr$cdr3_nt[which(gmp_tcr$cdr3_nt %in% wk1_tcr$cdr3_nt)])

### remove the "None"
target_seqs <- target_seqs[which(target_seqs != "None")]

### determine the number of columns in the result matrix
gmp_max <- 0
wk1_max <- 0
for(cdr3 in target_seqs) {
  gmp_idx <- intersect(union(grep(paste0(":",cdr3, ";"), JCC212_Px5@meta.data$cdr3_nt),
                             grep(paste0(":",cdr3, "$"), JCC212_Px5@meta.data$cdr3_nt)),
                       which(JCC212_Px5@meta.data$Time == "GMP"))
  wk1_idx <- intersect(union(grep(paste0(":",cdr3, ";"), JCC212_Px5@meta.data$cdr3_nt),
                             grep(paste0(":",cdr3, "$"), JCC212_Px5@meta.data$cdr3_nt)),
                       which(JCC212_Px5@meta.data$Time == "Wk1"))
  
  if(gmp_max < length(gmp_idx)) {
    gmp_max <- length(gmp_idx)
  }
  if(wk1_max < length(wk1_idx)) {
    wk1_max <- length(wk1_idx)
  }
}

### make an empty matrix for the result
result_mat <- matrix("", length(target_seqs), gmp_max+wk1_max)
rownames(result_mat) <- target_seqs
colnames(result_mat) <- c(paste0(rep("GMP_", gmp_max), 1:gmp_max), paste0(rep("Wk1_", wk1_max), 1:wk1_max))

### fill out the table
for(cdr3 in target_seqs) {
  gmp_idx <- intersect(union(grep(paste0(":",cdr3, ";"), JCC212_Px5@meta.data$cdr3_nt),
                             grep(paste0(":",cdr3, "$"), JCC212_Px5@meta.data$cdr3_nt)),
                       which(JCC212_Px5@meta.data$Time == "GMP"))
  wk1_idx <- intersect(union(grep(paste0(":",cdr3, ";"), JCC212_Px5@meta.data$cdr3_nt),
                             grep(paste0(":",cdr3, "$"), JCC212_Px5@meta.data$cdr3_nt)),
                       which(JCC212_Px5@meta.data$Time == "Wk1"))
  
  ### GMP
  if(length(gmp_idx) > 0) {
    for(i in 1:length(gmp_idx)) {
      result_mat[cdr3,i] <- JCC212_Px5@meta.data$cdr3_nt[gmp_idx[i]]
    }
  }
  ### Wk1
  if(length(wk1_idx) > 0) {
    for(i in 1:length(wk1_idx)) {
      result_mat[cdr3,i+gmp_max] <- JCC212_Px5@meta.data$cdr3_nt[wk1_idx[i]]
    }
  }
}

### save the table
write.xlsx2(data.frame(single_cdr3_nt=rownames(result_mat), result_mat,
                       stringsAsFactors = FALSE, check.names = FALSE),
            file = "./results/tracking_common_cdr3_nts.xlsx",
            row.names = FALSE)
