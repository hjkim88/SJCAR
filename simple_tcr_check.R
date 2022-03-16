### quick check to see shared TCRs between the new px14&15 GMP data
### and the previously created px14&15 PI data

### load library
if(!require(Seurat, quietly = TRUE)) {
  install.packages("Seurat")
  require(Seurat, quietly = TRUE)
}

### get TCR info file paths
TCR_data_dirs <- list.files(path = "Z:/ResearchHome/Groups/thomagrp/home/common/JCC212_SJCAR19/March2022/",
                            pattern = "filtered_contig_annotations.csv$",
                            full.names = TRUE, recursive = TRUE)

tcr <- NULL
for(i in 1:length(TCR_data_dirs)) {
  
  ### load filtered contig annotation data
  tcr_data <- read.csv(file = TCR_data_dirs[i],
                       stringsAsFactors = FALSE, check.names = FALSE)
  
  ### remove the -* at the end of each barcode.
  tcr_data$barcode <- sapply(strsplit(tcr_data$barcode, split = "-", fixed = TRUE), function(x) x[1])
  
  ### remove the rows that do not have CDR3 sequence
  tcr_data <- tcr_data[which(tcr_data$cdr3 != "None"),]
  
  ### remove redundant rows
  tcr_data <- tcr_data[!duplicated(tcr_data[c("barcode", "cdr3")]),]
  
  ### order by "chain" so that TRA rows come first than TRB rows
  ### and secondly, order by "CDR3 Nucleotide" sequence
  tcr_data <- tcr_data[order(as.character(tcr_data$chain), as.character(tcr_data$cdr3_nt)),]
  
  ### annotate TRA & TRB info to the cdr3 sequences
  tcr_data$cdr3 <- paste0(tcr_data$chain, ":", tcr_data$cdr3)
  tcr_data$cdr3_nt <- paste0(tcr_data$chain, ":", tcr_data$cdr3_nt)
  
  ### now merge different TRA & TRB info to one row
  dups <- which(duplicated(tcr_data$barcode))
  if(length(dups) > 0) {
    temp <- tcr_data[dups,]
    tcr_data <- tcr_data[-dups,]
    rownames(tcr_data) <- tcr_data$barcode
    for(barcode in tcr_data$barcode) {
      idx <- which(temp$barcode == barcode)
      tcr_data[barcode,"cdr3"] <- paste(c(tcr_data[barcode,"cdr3"], temp[idx, "cdr3"]), collapse = ";")
      tcr_data[barcode,"cdr3_nt"] <- paste(c(tcr_data[barcode,"cdr3_nt"], temp[idx, "cdr3_nt"]), collapse = ";")
      tcr_data[barcode,"productive"] <- paste(c(tcr_data[barcode,"productive"], temp[idx, "productive"]), collapse = ";")
      tcr_data[barcode,"reads"] <- paste(c(tcr_data[barcode,"reads"], temp[idx, "reads"]), collapse = ";")
      tcr_data[barcode,"umis"] <- paste(c(tcr_data[barcode,"umis"], temp[idx, "umis"]), collapse = ";")
    }
  }
  
  ### only retain informative columns
  tcr_data <- tcr_data[,c("barcode", "raw_clonotype_id", "cdr3", "cdr3_nt", "reads", "umis", "productive")]
  colnames(tcr_data) <- c("barcode", "raw_clonotype_id", "cdr3_aa", "cdr3_nt", "tcr_reads", "tcr_umis", "tcr_productive")
  
  ### add time & patient info
  if(i == 1) {
    tcr_data$time <- "GMP"
    tcr_data$type <- "GMP"
    tcr_data$px <- "SJCAR19-14"
    
  } else if(i == 2) {
    tcr_data$time <- "GMP"
    tcr_data$type <- "GMP"
    tcr_data$px <- "SJCAR19-15"
  }
  
  ### combine the TCR data for all the data
  if(is.null(tcr)) {
    tcr <- tcr_data
  } else {
    tcr <- rbind(tcr, tcr_data)
  }
  
}

### annotate library info
tcr$library <- paste("JCC212", tcr$px, tcr$time, tcr$type, sep = "_")

### combining alpha and beta chains
tcr$cdr3_alpha <- NA
tcr$cdr3_alpha_umis <- NA
tcr$cdr3_beta <- NA
tcr$cdr3_beta_umis <- NA
tcr$cdr3_one_alpha_beta <- NA
for(i in 1:nrow(tcr)) {
  
  ### split CDR3 AA
  aa_split <- strsplit(tcr$cdr3_aa[i], split = ";", fixed = TRUE)[[1]]
  umi_split <- as.numeric(strsplit(tcr$tcr_umis[i], split = ";", fixed = TRUE)[[1]])
  
  ### get alpha & beta indicies
  alpha_idx <- which(startsWith(aa_split, "TRA"))
  beta_idx <- which(startsWith(aa_split, "TRB"))
  
  ### fill out the columns
  tcr$cdr3_alpha[i] <- paste(aa_split[alpha_idx], collapse = ";")
  tcr$cdr3_alpha_umis[i] <- paste(umi_split[alpha_idx], collapse = ";")
  tcr$cdr3_beta[i] <- paste(aa_split[beta_idx], collapse = ";")
  tcr$cdr3_beta_umis[i] <- paste(umi_split[beta_idx], collapse = ";")
  if(length(alpha_idx) > 0 && length(beta_idx) > 0) {
    ### if there are multiple chains that have the same largest number of UMIs, just use the first one
    best_alpha_idx <- alpha_idx[which(umi_split[alpha_idx] == max(umi_split[alpha_idx]))][1]
    best_beta_idx <- beta_idx[which(umi_split[beta_idx] == max(umi_split[beta_idx]))][1]
    tcr$cdr3_one_alpha_beta[i] <- paste0(aa_split[best_alpha_idx], ";", aa_split[best_beta_idx])
  } else if(length(alpha_idx) > 0 && length(beta_idx) == 0) {
    best_alpha_idx <- alpha_idx[which(umi_split[alpha_idx] == max(umi_split[alpha_idx]))][1]
    tcr$cdr3_one_alpha_beta[i] <- aa_split[best_alpha_idx]
  } else if(length(alpha_idx) == 0 && length(beta_idx) > 0) {
    best_beta_idx <- beta_idx[which(umi_split[beta_idx] == max(umi_split[beta_idx]))][1]
    tcr$cdr3_one_alpha_beta[i] <- aa_split[best_beta_idx]
  }
  
  ### progress
  if(i %% 1000 == 0) {
    writeLines(paste(i, "/", nrow(tcr)))
  }
  
}

### change "" to NA
tcr$cdr3_alpha[which(tcr$cdr3_alpha == "")] <- NA
tcr$cdr3_alpha_umis[which(tcr$cdr3_alpha_umis == "")] <- NA
tcr$cdr3_beta[which(tcr$cdr3_beta == "")] <- NA
tcr$cdr3_beta_umis[which(tcr$cdr3_beta_umis == "")] <- NA
tcr$cdr3_one_alpha_beta[which(tcr$cdr3_one_alpha_beta == "")] <- NA


### load Jeremy's object
JCC_Seurat_Obj <- readRDS(file = "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds")

### check whether the orders are the same
print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
print(identical(names(Idents(object = JCC_Seurat_Obj)), rownames(JCC_Seurat_Obj@meta.data)))

### only retain the meta.data
JCC_Seurat_meta <- JCC_Seurat_Obj@meta.data

### remove the large object
rm(JCC_Seurat_Obj)
gc()

### compare px14
px14_gmp_tcrs <- tcr[intersect(which(tcr$time == "GMP"),
                               which(tcr$px == "SJCAR19-14")),]
px14_pi_tcrs <- JCC_Seurat_meta[intersect(which(JCC_Seurat_meta$GMP == "PI"),
                                          which(JCC_Seurat_meta$px == "SJCAR19-14")),]

print(paste("Px14 GMP TCR # = ", length(which(!is.na(px14_gmp_tcrs$cdr3_aa)))))
print(paste("Px14 PI TCR # = ", length(which(!is.na(px14_pi_tcrs$cdr3_aa)))))
print(paste("Px14 Shared TCR # = ", length(intersect(px14_gmp_tcrs$cdr3_one_alpha_beta[which(!is.na(px14_gmp_tcrs$cdr3_one_alpha_beta))],
                                                     px14_pi_tcrs$cdr3_one_alpha_beta[which(!is.na(px14_gmp_tcrs$cdr3_one_alpha_beta))]))))

### compare px15
px15_gmp_tcrs <- tcr[intersect(which(tcr$time == "GMP"),
                               which(tcr$px == "SJCAR19-15")),]
px15_pi_tcrs <- JCC_Seurat_meta[intersect(which(JCC_Seurat_meta$GMP == "PI"),
                                          which(JCC_Seurat_meta$px == "SJCAR19-15")),]

print(paste("Px15 GMP TCR # = ", length(which(!is.na(px15_gmp_tcrs$cdr3_aa)))))
print(paste("Px15 PI TCR # = ", length(which(!is.na(px15_pi_tcrs$cdr3_aa)))))
print(paste("Px15 Shared TCR # = ", length(intersect(px15_gmp_tcrs$cdr3_one_alpha_beta[which(!is.na(px15_gmp_tcrs$cdr3_one_alpha_beta))],
                                                     px15_pi_tcrs$cdr3_one_alpha_beta[which(!is.na(px15_gmp_tcrs$cdr3_one_alpha_beta))]))))




