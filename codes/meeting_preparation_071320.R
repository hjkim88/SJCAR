### This is just to create some useful plots for the immuno-core meeting presentation

Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj"

outputDir="./etc/071320/"
dir.create(path = outputDir, showWarnings = FALSE, recursive = TRUE)

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

### load the Seurat object and save the object name
tmp_env <- new.env()
load(Seurat_RObj_path, tmp_env)
obj_name <- ls(tmp_env)
assign("Seurat_Obj", get(obj_name, envir = tmp_env))
rm(tmp_env)
gc()

### rownames in the meta.data should be in the same order as colnames in the counts
Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]

### set unique patients and time points
px <- unique(Seurat_Obj@meta.data$Px)
tp <- levels(Seurat_Obj@meta.data$TimeF)

### make tables for all the cells (GEX & TCR) - the number of cells
all_gex_mat <- matrix(0, length(px), length(tp))
colnames(all_gex_mat) <- tp
rownames(all_gex_mat) <- px
all_tcr_mat <- matrix(0, length(px), length(tp))
colnames(all_tcr_mat) <- tp
rownames(all_tcr_mat) <- px

### fill out the tables
for(p in px) {
  for(t in tp) {
    all_gex_mat[p, t] <- length(intersect(which(Seurat_Obj@meta.data$Px == p),
                                          which(Seurat_Obj@meta.data$Time == t)))
  }
}
for(p in px) {
  for(t in tp) {
    all_tcr_mat[p, t] <- length(intersect(intersect(which(Seurat_Obj@meta.data$Px == p),
                                                    which(Seurat_Obj@meta.data$Time == t)),
                                          which(!is.na(Seurat_Obj@meta.data$cdr3_aa))))
  }
}

### make tables for CAR+ cells (GEX & TCR) - the number of cells
car_gex_mat <- matrix(0, length(px), length(tp))
colnames(car_gex_mat) <- tp
rownames(car_gex_mat) <- px
car_tcr_mat <- matrix(0, length(px), length(tp))
colnames(car_tcr_mat) <- tp
rownames(car_tcr_mat) <- px

### fill out the tables
for(p in px) {
  for(t in tp) {
    car_gex_mat[p, t] <- length(intersect(intersect(which(Seurat_Obj@meta.data$Px == p),
                                                    which(Seurat_Obj@meta.data$Time == t)),
                                          which(Seurat_Obj@meta.data$CAR == "CARpos")))
  }
}
for(p in px) {
  for(t in tp) {
    car_tcr_mat[p, t] <- length(intersect(intersect(intersect(which(Seurat_Obj@meta.data$Px == p),
                                                              which(Seurat_Obj@meta.data$Time == t)),
                                                    which(!is.na(Seurat_Obj@meta.data$cdr3_aa))),
                                          which(Seurat_Obj@meta.data$CAR == "CARpos")))
  }
}

### write out the tables
write.xlsx2(data.frame(all_gex_mat, check.names = FALSE), file = paste0(outputDir, "all_gex_table.xlsx"))
write.xlsx2(data.frame(all_tcr_mat, check.names = FALSE), file = paste0(outputDir, "all_tcr_table.xlsx"))
write.xlsx2(data.frame(car_gex_mat, check.names = FALSE), file = paste0(outputDir, "car_gex_table.xlsx"))
write.xlsx2(data.frame(car_tcr_mat, check.names = FALSE), file = paste0(outputDir, "car_tcr_table.xlsx"))


### how many cells are in the TCR but not in the GEX?

### get TCR info file paths
TCR_data_dirs="C:/Users/hkim8/SJ/SJCAR19/TCRoutputs/"
TCR_data_dirs <- list.files(path = TCR_data_dirs, pattern = "filtered_contig_annotations.csv$",
                            full.names = TRUE, recursive = TRUE)

### there are two cases of duplicated files that are for exactly the same condition 
### [32] 1757058_JCC212_HealthyDonor32_PreTrans_TCR_filtered_contig_annotations.csv
### [44] 1757070_JCC212_HealthyDonor32_PreTrans_TCR_filtered_contig_annotations.csv
### [33] 1757059_JCC212_HealthyDonor33_PreTrans_TCR_filtered_contig_annotations.csv
### [45] 1757071_JCC212_HealthyDonor33_PreTrans_TCR_filtered_contig_annotations.csv
### I manually checked that the [44] and the [45] have the same libraries (barcodes) as the GEX data (Seurat object)
### so, remove the [32] and the [33] here
TCR_data_dirs <- TCR_data_dirs[-c(32, 33)]

### for each TCR data combine them to the Seurat object
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
  temp <- strsplit(basename(TCR_data_dirs[i]), split = "_", fixed = TRUE)[[1]]
  px <- temp[3]
  temp <- temp[-c(1,2,3)]
  if(temp[1] == "PB" || temp[1] == "BM") {
    tcr_data$time <- temp[2]
    tcr_data$type <- temp[1]
    tcr_data$px <- px
  } else if(grepl("GMP", temp[1]) || grepl("GMP", px)) {
    tcr_data$time <- "GMP"
    tcr_data$type <- "GMP"
    tcr_data$px <- paste0("SJCAR19-Donor", substring(px, 9))
  } else if(temp[2] == "PB" || temp[2] == "BM") {
    tcr_data$time <- temp[1]
    tcr_data$type <- temp[2]
    tcr_data$px <- px
  } else {
    tcr_data$time <- temp[1]
    tcr_data$type <- "PB"
    tcr_data$px <- px
  }
  if(grepl("HealthyDonor", px)) {
    tcr_data$px <- paste0("SJCAR19-Donor", substring(px, 13))
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
idx <- grep("Donor", tcr$px)
tcr$library[idx] <- paste("JCC212", tcr$px[idx], tcr$time[idx], sep = "_")
### there are inconsistent naming between the GEX and the TCR
### there is no rule/consistency for that, so I have to manually change them all
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor02_GMP")] <- "JCC212_SJCAR19-02_GMP"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor00_GMP")] <- "JCC212_SJCAR19-00_GMP"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor01_GMP")] <- "JCC212_SJCAR19-01_GMP"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor03_GMP")] <- "JCC212_SJCAR19-03_GMP"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor04_GMP")] <- "JCC212_SJCAR19-04_GMP"
tcr$library[which(tcr$library == "JCC212_SJCAR19-02_PreTrans_PB")] <- "JCC212_SJCAR19-02_PreTrans"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor32_PreTrans")] <- "JCC212_SJCAR19-Donor32_PreTransB"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor33_PreTrans")] <- "JCC212_SJCAR19-Donor33_PreTransB"
tcr$library[which(tcr$library == "JCC212_SJCAR19-03_PreTrans_PB")] <- "JCC212_SJCAR19-03_PreTransB"
tcr$library[which(tcr$library == "JCC212_SJCAR19-04_PreTrans_PB")] <- "JCC212_SJCAR19-04_PreTransB"
tcr$library[which(tcr$library == "JCC212_SJCAR19-00_PreTrans_PB")] <- "JCC212_SJCAR19-00_PreTransB"
tcr$library[which(tcr$library == "JCC212_SJCAR19-05_PreTrans_PB")] <- "JCC212_SJCAR19-05_PreTransB"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor30_PreTrans")] <- "JCC212_SJCAR19-Donor30_PreTransB"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor05_GMP")] <- "JCC212_SJCAR19-05_GMP"
tcr$library[which(tcr$library == "JCC212_SJCAR19-06_PreTrans_PB")] <- "JCC212_SJCAR19-06_PreTrans"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor07_GMP")] <- "JCC212_SJCAR19-07_GMP19047"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor02_GMP")] <- "JCC212_SJCAR19-02_GMP"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor02_GMP")] <- "JCC212_SJCAR19-02_GMP"
tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor06_GMP")] <- "JCC212_SJCAR19-06_GMP"

tcr$px[which(tcr$px == "SJCAR19-Donor00")] <- "SJCAR19-00"
tcr$px[which(tcr$px == "SJCAR19-Donor01")] <- "SJCAR19-01"
tcr$px[which(tcr$px == "SJCAR19-Donor02")] <- "SJCAR19-02"
tcr$px[which(tcr$px == "SJCAR19-Donor03")] <- "SJCAR19-03"
tcr$px[which(tcr$px == "SJCAR19-Donor04")] <- "SJCAR19-04"
tcr$px[which(tcr$px == "SJCAR19-Donor05")] <- "SJCAR19-05"
tcr$px[which(tcr$px == "SJCAR19-Donor06")] <- "SJCAR19-06"
tcr$px[which(tcr$px == "SJCAR19-Donor07")] <- "SJCAR19-07"


### matrix for scTCR data only
pure_tcr_mat <- matrix(0, length(px), length(tp))
colnames(pure_tcr_mat) <- tp
rownames(pure_tcr_mat) <- px

### fill out the tables
for(p in px) {
  for(t in tp) {
    pure_tcr_mat[p, t] <- length(intersect(which(tcr$px == p),
                                          which(tcr$time == t)))
  }
}

### write out the tables
write.xlsx2(data.frame(pure_tcr_mat, check.names = FALSE), file = paste0(outputDir, "pure_tcr_table.xlsx"))
