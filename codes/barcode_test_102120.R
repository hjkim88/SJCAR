### there are re-sequenced GEX barcodes of the px05 of SJCAR19
### check the number of shared barcodes with the GEX and the TCR

load("~/SJCAR/101520.RData")

new_TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRs_15Oct2020/"
barcode_dir="C:/Users/hkim8/SJ/SJCAR19/GEXbarcodes_test/"
old_barcode_dir="C:/Users/hkim8/SJ/SJCAR19/GEXbarcodes_15Oct2020/"
outputDir="./results/TCR_1021/"

### load libraries
if(!require(xlsx, quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}
if(!require(RColorBrewer, quietly = TRUE)) {
  install.packages("RColorBrewer")
  require(RColorBrewer, quietly = TRUE)
}

### create result dir
dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

### get real tcr names
real_tcr_names <- sapply(names(tcr), function(x) {
  if(startsWith(x, "X_")) {
    return(substring(x, 3))
  } else {
    return(substring(x, 9))
  }
}, USE.NAMES = FALSE)

### match real_tcr_names to the names(real_barcodes)
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-02_Wk1")] <- "JCC212_SJCAR19-02_PB_Wk1"
real_tcr_names[which(real_tcr_names == "JCC212_GMPdonor30")] <- "JCC212_SJCAR19_GMPdonor30"
real_tcr_names[which(real_tcr_names == "JCC212_GMPdonor32")] <- "JCC212_SJCAR19_GMPdonor32"
real_tcr_names[which(real_tcr_names == "JCC212_GMPdonor33")] <- "JCC212_SJCAR19_GMPdonor33"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_Wk-1")] <- "JCC212_SJCAR19-04_PB_Wk-1"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_Wk2")] <- "JCC212_SJCAR19_04_Wk2"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_GMP19003")] <- "JCC212_SJCAR19_04_GMP19003"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_PB_Wk-1")] <- "JCC212_SJCAR19_05_PB_Wk-1"
real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor32_PreTrans")[1]] <- "JCC212_SJCAR19_Donor32_PreTrans"
real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor33_PreTrans")[1]] <- "JCC212_SJCAR19_Donor33_PreTrans"
real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor32_PreTrans")] <- "JCC212_SJCAR19_Donor32_PreTransB"
real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor33_PreTrans")] <- "JCC212_SJCAR19_Donor33_PreTransB"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_Wk8")] <- "JCC212_SJCAR19_04_Wk8"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_Wk2")] <- "JCC212_SJCAR19_05_Wk2"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_Wk3")] <- "JCC212_SJCAR19_05_Wk3"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_PB_Wk4")] <- "JCC212_SJCAR19_05_PB_Wk4"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_BM_Wk4")] <- "JCC212_SJCAR19_05_BM_Wk4"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-03_PreTrans")] <- "JCC212_SJCAR19_03_PreTransB"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_PreTrans")] <- "JCC212_SJCAR19_04_PreTransB"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-00_PreTrans")] <- "JCC212_SJCAR19_00_PreTransB"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_PreTrans")] <- "JCC212_SJCAR19_05_PreTransB"
real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor30_PreTrans")] <- "JCC212_SJCAR19_Donor30_PreTransB"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-06_Wk-1_PB")] <- "JCC212_SJCAR19-06_Wk-1"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR-06_GMP19028")] <- "JCC212_SJCAR19-06_GMP19028"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_3mo_PB")] <- "JCC212_SJCAR19-05_PB_3mo"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_3mo_BM")] <- "JCC212_SJCAR19-05_BM_3mo"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-06_PreTrans")] <- "JCC212_SJCAR19-6_PreTrans"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-06_3_mo_BM")] <- "JCC212_SJCAR19-06_3mo_BM"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-08_GMP")] <- "JCC212_SJCAR19-08_GMP19064"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-13_Wk-1")] <- "JCC212_SJCAR19-13_Wk-1_PB"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_Wk0")] <- "JCC212_SJCAR19_05_PB_Wk0"
real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_Wk1")] <- "JCC212_SJCAR19_05_Wk1"



### get new barcode info file paths
barcode_file_paths <- list.files(path = barcode_dir, pattern = "barcodes.tsv.gz$",
                                 full.names = TRUE, recursive = TRUE)

### load and combine the TCR data
real_barcodes <- vector("list", length = length(barcode_file_paths))
names(real_barcodes) <- sapply(basename(barcode_file_paths), function(x) {
  return(substr(x, 1, nchar(x)-29))
}, USE.NAMES = FALSE)


for(i in 1:length(barcode_file_paths)) {
  ### load barcode data
  real_barcodes[[i]] <- read.table(file = gzfile(barcode_file_paths[i]),
                                   header = FALSE,
                                   stringsAsFactors = FALSE, check.names = FALSE)
  
  ### remove the -* at the end of each barcode.
  real_barcodes[[i]]$V1 <- sapply(strsplit(real_barcodes[[i]]$V1, split = "-", fixed = TRUE), function(x) x[1])
}


### shared barcodes
shared_barcode_mat <- matrix(NA, length(tcr), length(real_barcodes))
rownames(shared_barcode_mat) <- names(tcr)
colnames(shared_barcode_mat) <- names(real_barcodes)
shared_barcode_pct <- matrix(NA, length(tcr), length(real_barcodes))
rownames(shared_barcode_pct) <- names(tcr)
colnames(shared_barcode_pct) <- names(real_barcodes)

for(i in 1:length(tcr)) {
  for(j in 1:length(real_barcodes)) {
    shared_barcode_mat[i,j] <- length(intersect(rownames(tcr[[i]]), real_barcodes[[j]]$V1))
    shared_barcode_pct[i,j] <- 100* shared_barcode_mat[i,j] / nrow(real_barcodes[[j]])
  }
}

### write out the result
write.xlsx2(data.frame(shared_barcode_mat, stringsAsFactors = FALSE, check.names = FALSE),
            file = paste0(outputDir, "Shared_Barcodes_TCR-GEX_px5.xlsx"),
            sheetName = "Shared_Barcodes_Numbers", append = FALSE)
write.xlsx2(data.frame(shared_barcode_pct, stringsAsFactors = FALSE, check.names = FALSE),
            file = paste0(outputDir, "Shared_Barcodes_TCR-GEX_px5.xlsx"),
            sheetName = "Shared_Barcodes_Percentages", append = TRUE)



### get old barcode info file paths
old_barcode_file_paths <- list.files(path = old_barcode_dir, pattern = "barcodes.tsv.gz$",
                                 full.names = TRUE, recursive = TRUE)

### load and combine the TCR data
old_real_barcodes <- vector("list", length = length(old_barcode_file_paths))
names(old_real_barcodes) <- sapply(basename(old_barcode_file_paths), function(x) {
  return(substr(x, 1, nchar(x)-31))
}, USE.NAMES = FALSE)


for(i in 1:length(old_barcode_file_paths)) {
  ### load barcode data
  old_real_barcodes[[i]] <- read.table(file = gzfile(old_barcode_file_paths[i]),
                                   header = FALSE,
                                   stringsAsFactors = FALSE, check.names = FALSE)
  
  ### remove the -* at the end of each barcode.
  old_real_barcodes[[i]]$V1 <- sapply(strsplit(old_real_barcodes[[i]]$V1, split = "-", fixed = TRUE), function(x) x[1])
}


### shared barcodes
shared_barcode_mat <- matrix(NA, length(old_real_barcodes), length(real_barcodes))
rownames(shared_barcode_mat) <- names(old_real_barcodes)
colnames(shared_barcode_mat) <- names(real_barcodes)
shared_barcode_pct <- matrix(NA, length(old_real_barcodes), length(real_barcodes))
rownames(shared_barcode_pct) <- names(old_real_barcodes)
colnames(shared_barcode_pct) <- names(real_barcodes)

for(i in 1:length(old_real_barcodes)) {
  for(j in 1:length(real_barcodes)) {
    shared_barcode_mat[i,j] <- length(intersect(old_real_barcodes[[i]]$V1, real_barcodes[[j]]$V1))
    shared_barcode_pct[i,j] <- 100* shared_barcode_mat[i,j] / nrow(real_barcodes[[j]])
  }
}

### write out the result
write.xlsx2(data.frame(shared_barcode_mat, stringsAsFactors = FALSE, check.names = FALSE),
            file = paste0(outputDir, "Shared_Barcodes_GEX-GEX_px5.xlsx"),
            sheetName = "Shared_Barcodes_Numbers", append = FALSE)
write.xlsx2(data.frame(shared_barcode_pct, stringsAsFactors = FALSE, check.names = FALSE),
            file = paste0(outputDir, "Shared_Barcodes_GEX-GEX_px5.xlsx"),
            sheetName = "Shared_Barcodes_Percentages", append = TRUE)


