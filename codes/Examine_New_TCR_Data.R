###
#   File name : Examine_New_TCR_Data.R
#   Author    : Hyunjin Kim
#   Date      : Sep 22, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Now we have new TCR data. Compare it to the old TCR data and
#               examine how many more lineages can be made with the new data.
#
#   Instruction
#               1. Source("Examine_New_TCR_Data.R")
#               2. Run the function "examine_new_tcr_data" - specify the input directory and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Examine_New_TCR_Data.R/Examine_New_TCR_Data.R")
#               > examine_new_tcr_data(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                                      new_TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRredo/",
#                                      old_TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRoutputs/",
#                                      outputDir="./results/")
###

examine_new_tcr_data <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                                 new_TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRredo/",
                                 old_TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRoutputs/",
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
  
  ### only retain the needed data and remove the rest
  metadata <- Seurat_Obj@meta.data[union(which(Seurat_Obj@meta.data$Px == "SJCAR19-03"),
                                         which(Seurat_Obj@meta.data$Px == "SJCAR19-07")),]
  rm(Seurat_Obj)
  gc()
  
  ### get old TCR info file paths
  old_TCR_data_dirs <- list.files(path = old_TCR_dir, pattern = "filtered_contig_annotations.csv$",
                                  full.names = TRUE, recursive = TRUE)
  old_TCR_data_dirs <- old_TCR_data_dirs[-c(32, 33)]
  
  ### keep only the patient 3 & 7
  old_TCR_data_dirs <- old_TCR_data_dirs[union(grep(pattern = "SJCAR19-03", old_TCR_data_dirs, fixed = TRUE),
                                               grep(pattern = "SJCAR19-07", old_TCR_data_dirs, fixed = TRUE))]
  
  ### get new TCR info file paths
  new_TCR_data_dirs <- list.files(path = new_TCR_dir, pattern = "filtered_contig_annotations.csv$",
                                  full.names = TRUE, recursive = TRUE)
  
  ### combine old & new file paths
  TCR_data_dirs <- c(old_TCR_data_dirs, new_TCR_data_dirs)
  
  ### load and combine the TCR data
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
    } else if(grepl("GMP-redo", temp[1])) {
      tcr_data$time <- "GMP-redo"
      tcr_data$type <- "GMP"
      tcr_data$px <- paste0("SJCAR19-Donor", substring(px, 9))
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
  tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor03_GMP")] <- "JCC212_SJCAR19-03_GMP"
  tcr$library[which(tcr$library == "JCC212_SJCAR19-03_PreTrans_PB")] <- "JCC212_SJCAR19-03_PreTransB"
  tcr$library[which(tcr$library == "JCC212_SJCAR19-06_PreTrans_PB")] <- "JCC212_SJCAR19-06_PreTrans"
  tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor07_GMP")] <- "JCC212_SJCAR19-07_GMP19047"
  tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor03_GMP-redo")] <- "JCC212_SJCAR19-03_GMP-redo"
  tcr$library[which(tcr$library == "JCC212_SJCAR19-Donor07_GMP-redo")] <- "JCC212_SJCAR19-07_GMP-redo"
  
  ### correct some mis-labeled Donor PX
  tcr$px[which(tcr$px == "SJCAR19-Donor03")] <- "SJCAR19-03"
  tcr$px[which(tcr$px == "SJCAR19-Donor07")] <- "SJCAR19-07"
  
  ### we are not interested in PreTrans & Wk-1, so remove those
  tcr <- tcr[-union(which(tcr$time == "PreTrans"),
                    which(tcr$time == "Wk-1")),]
  
  ### factorize the time column
  tcr$timeF <- factor(tcr$time, levels = c("GMP", "GMP-redo", "Wk1", "Wk2", "Wk3", "Wk4", "Wk6", "Wk8"))
  
  ### order the data by timeF
  tcr <- tcr[order(tcr$px, tcr$timeF),]
  
  ### annotate CAR+/-
  tcr$CAR <- NA
  for(i in 1:nrow(tcr)) {
    barcode <- tcr$barcode[i]
    px <- tcr$px[i]
    time <- tcr$time[i]
    if(time == "GMP-redo") {
      time <- "GMP"
    }
    type <- tcr$type[i]
    
    idx <- intersect(intersect(which(metadata$GexCellShort == barcode),
                               which(metadata$Px == px)),
                     intersect(which(metadata$Time == time),
                               which(metadata$Type == type)))
    
    if(length(idx) > 1) {
      writeLines(paste(i, ": There are duplicated barcodes."))
    } else if(length(idx) == 0) {
      # writeLines(paste(i, ": There are no matched barcodes."))
    } else {
      tcr$CAR[i] <- metadata$CAR[idx]
    }
    
    if(i %% 1000 == 0) {
      writeLines(paste(i, "/", nrow(tcr)))
    }
  }
  
  ### initialize the necessary variables
  cloneNums <- vector("list", 2)
  lineageNums_old <- vector("list", 2)
  lineageNums_new <- vector("list", 2)
  lineageNums_combined <- vector("list", 2)
  cloneNums_CAR <- vector("list", 2)
  lineageNums_old_CAR <- vector("list", 2)
  lineageNums_new_CAR <- vector("list", 2)
  lineageNums_combined_CAR <- vector("list", 2)
  names(cloneNums) <- unique(tcr$px)
  names(lineageNums_old) <- unique(tcr$px)
  names(lineageNums_new) <- unique(tcr$px)
  names(lineageNums_combined) <- unique(tcr$px)
  names(cloneNums_CAR) <- unique(tcr$px)
  names(lineageNums_old_CAR) <- unique(tcr$px)
  names(lineageNums_new_CAR) <- unique(tcr$px)
  names(lineageNums_combined_CAR) <- unique(tcr$px)
  
  ### clonotyping
  time_points <- unique(tcr$time)
  tcr$clone_id <- NA
  for(px in unique(tcr$px)) {
    
    ### get patient indices
    px_idx <- which(tcr$px == px)
    
    ### give the strict version of global clonotypes using all alpha & beta chains
    ### since the variables in each of the "cdr3_nt" are ordered,
    ### if the "cdr3_nt" are exactly the same, they are the same clonotype
    
    ### the unique "cdr3_nt" sequences
    unique_seqs <- unique(tcr$cdr3_nt[px_idx])
    
    ### remove NA, "NA", and ""
    unique_seqs <- unique_seqs[intersect(intersect(which(!is.na(unique_seqs)),
                                                   which(unique_seqs != "NA")),
                                         which(unique_seqs != ""))]
    
    ### get unique and duplicated indicies
    unique_idx <- which(!duplicated(tcr$cdr3_nt[px_idx]))
    dup_idx <- which(duplicated(tcr$cdr3_nt[px_idx])) 
    
    ### give clone ids
    tcr$clone_id[px_idx[unique_idx]] <- paste0("clone", 1:length(unique_idx))
    for(i in px_idx[dup_idx]) {
      target_idx <- which(tcr$cdr3_nt == tcr$cdr3_nt[i])
      tcr$clone_id[i] <- tcr$clone_id[target_idx[1]]
    }
    
    ### number of clones and lineages
    
    ### initialize the necessary variables
    # ALL
    cloneNums[[px]] <- rep(NA, length(time_points))
    names(cloneNums[[px]]) <- time_points
    lineageNums_old[[px]] <- rep(NA, length(time_points))
    names(lineageNums_old[[px]]) <- time_points
    lineageNums_new[[px]] <- rep(NA, length(time_points))
    names(lineageNums_new[[px]]) <- time_points
    lineageNums_combined[[px]] <- rep(NA, length(time_points))
    names(lineageNums_combined[[px]]) <- time_points
    # CAR+ only
    cloneNums_CAR[[px]] <- rep(NA, length(time_points))
    names(cloneNums_CAR[[px]]) <- time_points
    lineageNums_old_CAR[[px]] <- rep(NA, length(time_points))
    names(lineageNums_old_CAR[[px]]) <- time_points
    lineageNums_new_CAR[[px]] <- rep(NA, length(time_points))
    names(lineageNums_new_CAR[[px]]) <- time_points
    lineageNums_combined_CAR[[px]] <- rep(NA, length(time_points))
    names(lineageNums_combined_CAR[[px]]) <- time_points
    
    ### get indicies for GMPs
    old_gmp_idx <- intersect(which(tcr$time == "GMP"), px_idx)
    new_gmp_idx <- intersect(which(tcr$time == "GMP-redo"), px_idx)
    old_gmp_CAR_idx <- intersect(which(tcr$CAR == "CARpos"), old_gmp_idx)
    new_gmp_CAR_idx <- intersect(which(tcr$CAR == "CARpos"), new_gmp_idx)
    
    ### get cloneNums and lineageNums
    for(time in time_points) {
      ### get indicies for the given time point
      time_idx <- intersect(which(tcr$time == time), px_idx)
      time_CAR_idx <- intersect(which(tcr$CAR == "CARpos"), time_idx)
      
      ### the number of clones
      cloneNums[[px]][time] <- length(unique(tcr$clone_id[time_idx]))
      cloneNums_CAR[[px]][time] <- length(unique(tcr$clone_id[time_CAR_idx]))
      
      ### the number of lineages between GMP
      # ALL
      lineageNums_old[[px]][time] <- length(unique(intersect(tcr$clone_id[time_idx],
                                                             tcr$clone_id[old_gmp_idx])))
      lineageNums_new[[px]][time] <- length(unique(intersect(tcr$clone_id[time_idx],
                                                             tcr$clone_id[new_gmp_idx])))
      lineageNums_combined[[px]][time] <- length(unique(intersect(tcr$clone_id[time_idx],
                                                                  tcr$clone_id[union(old_gmp_idx,
                                                                                     new_gmp_idx)])))
      # CAR+ only
      lineageNums_old_CAR[[px]][time] <- length(unique(intersect(tcr$clone_id[time_CAR_idx],
                                                                 tcr$clone_id[old_gmp_CAR_idx])))
      lineageNums_new_CAR[[px]][time] <- length(unique(intersect(tcr$clone_id[time_CAR_idx],
                                                                 tcr$clone_id[new_gmp_CAR_idx])))
      lineageNums_combined_CAR[[px]][time] <- length(unique(intersect(tcr$clone_id[time_CAR_idx],
                                                                      tcr$clone_id[union(old_gmp_CAR_idx,
                                                                                          new_gmp_CAR_idx)])))
    }
    
    ### data frame for barplot
    plot_df <- data.frame(cloneNums=cloneNums[[px]],
                          lineageNums_old=lineageNums_old[[px]],
                          lineageNums_new=lineageNums_new[[px]],
                          lineageNums_combined=lineageNums_combined[[px]],
                          cloneNums_CAR=cloneNums_CAR[[px]],
                          lineageNums_old_CAR=lineageNums_old_CAR[[px]],
                          lineageNums_new_CAR=lineageNums_new_CAR[[px]],
                          lineageNums_combined_CAR=lineageNums_combined_CAR[[px]],
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### remove all zero rows
    sums <- apply(plot_df, 1, sum)
    plot_df <- plot_df[-which(sums == 0),]
    
    ### write out the result
    write.xlsx2(data.frame(Time=rownames(plot_df),plot_df),
                file = paste0(outputDir, "Px3_7_TCR_GMP-redo_Summary.xlsx"),
                sheetName = px,
                row.names = FALSE,
                append = TRUE)
    
  }
  
  
  
  
  
  
}
