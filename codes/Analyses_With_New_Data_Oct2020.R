###
#   File name : Analyses_With_New_Data_Oct2020.R
#   Author    : Hyunjin Kim
#   Date      : Oct 19, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Perform all the analyses with the new Seurat Object of SJCAR19
#
#   Instruction
#               1. Source("Analyses_With_New_Data_Oct2020.R")
#               2. Run the function "analyses_with_new_data" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Analyses_With_New_Data_Oct2020.R/Analyses_With_New_Data_Oct2020.R")
#               > analyses_with_new_data(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/SJCAR19_Oct2020_Seurat_Obj2.RDS",
#                                        barcode_dir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/JCC/JCC212_SJCAR19/GEXbarcodes_15Oct2020/",
#                                        TCR_dir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/JCC/JCC212_SJCAR19/TCRs_15Oct2020/",
#                                        clonotype_lineage_info_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/SJCAR19_Clonotype_Lineages.RDS",
#                                        outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/")
###

analyses_with_new_data <- function(Seurat_RObj_path="./data/SJCAR19_Oct2020_Seurat_Obj2.RDS",
                                   barcode_dir="C:/Users/hkim8/SJ/SJCAR19/GEXbarcodes_15Oct2020/",
                                   TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRs_15Oct2020/",
                                   clonotype_lineage_info_path="./data/SJCAR19_Clonotype_Lineages.RDS",
                                   outputDir="./results/New/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(ggalluvial, quietly = TRUE)) {
    install.packages("ggalluvial")
    require(ggalluvial, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  if(!require(viridis, quietly = TRUE)) {
    install.packages("viridis")
    require(viridis, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    library(scales, quietly = TRUE)
  }
  if(!require(slingshot, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("slingshot")
    require(slingshot, quietly = TRUE)
  }
  if(!require(caret, quietly = TRUE)) {
    install.packages("caret")
    require(caret, quietly = TRUE)
  }
  if(!require(e1071, quietly = TRUE)) {
    install.packages("e1071")
    require(e1071, quietly = TRUE)
  }
  if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
  }
  
  ### create outputDir
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### load Seurat object
  Seurat_Obj <- readRDS(Seurat_RObj_path)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
  
  ### active assay = "RNA"
  Seurat_Obj@active.assay <- "RNA"
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  ### set time points
  total_time_points <- c("PreTrans", "PreTransB", "Wk-1", "Wk-1Run1", "Wk-1Run2", "Wk0", "GMP", "GMP-redo",
                         "Wk1", "Wk1b", "Wk2", "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo")
  
  ### assign CAR+ cells
  Seurat_Obj@meta.data$CAR <- sapply(Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short",rownames(Seurat_Obj@meta.data)],
                                     function(x) {
                                       if(x > 0) {
                                         return("CARpos")
                                       } else {
                                         return("CARneg")
                                       }
                                     })
  
  ### assign short version of barcodes
  Seurat_Obj@meta.data$barcode_tag <- sapply(strsplit(rownames(Seurat_Obj@meta.data), split = "-", fixed = TRUE), function(x) x[2])
  Seurat_Obj@meta.data$barcode_short <- sapply(strsplit(rownames(Seurat_Obj@meta.data), split = "-", fixed = TRUE), function(x) x[1])
  
  ### get new barcode info file paths
  barcode_file_paths <- list.files(path = barcode_dir, pattern = "barcodes.tsv.gz$",
                                   full.names = TRUE, recursive = TRUE)
  
  ### load the barcodes
  barcodes <- vector("list", length = length(barcode_file_paths))
  names(barcodes) <- sapply(basename(barcode_file_paths), function(x) {
    return(substr(x, 1, nchar(x)-31))
  }, USE.NAMES = FALSE)
  for(i in 1:length(barcode_file_paths)) {
    ### load barcode data
    barcodes[[i]] <- read.table(file = gzfile(barcode_file_paths[i]),
                                header = FALSE,
                                stringsAsFactors = FALSE, check.names = FALSE)
    
    # ### remove the -* at the end of each barcode.
    barcodes[[i]]$V1 <- sapply(strsplit(barcodes[[i]]$V1, split = "-", fixed = TRUE), function(x) x[1])
  }
  
  ### refine barcode library name
  lib_name <- sapply(names(barcodes), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    px <- temp[2]
    time <- temp[3]
    tissue <- temp[4]
    
    ### exceptions
    # px
    if(px == "SJCAR19-05relapse-2ndInfusionWk1") {
      px <- "SJCAR19-05"
    } else if(px == "SJCAR19-05relapse-2ndInfusionWk2") {
      px <- "SJCAR19-05"
    } else if(px == "SJCAR19-6") {
      px <- "SJCAR19-06"
    } else if(px == "SJCAR19") {
      if(grepl("Donor", temp[3], fixed = TRUE)) {
        px <- paste0(temp[2], "-", strsplit(temp[3], split = "Donor", fixed = TRUE)[[1]][2])
      } else if(grepl("donor", temp[3], fixed = TRUE)) {
        px <- paste0(temp[2], "-", strsplit(temp[3], split = "donor", fixed = TRUE)[[1]][2])
      } else {
        px <- paste0(temp[2], "-", temp[3])
      }
    }
    # time
    if(grepl("GMP-redo", time, fixed = TRUE)) {
      time <- "GMP-redo"
    } else if(grepl("GMP", time, fixed = TRUE)) {
      time <- "GMP"
    }
    if(time == "PB" || time == "BM" || startsWith(time, "0") || startsWith(time, "Donor")) {
      time <- temp[4]
    }
    if(x == "JCC212_SJCAR19-05relapse-2ndInfusionWk1_PB") {
      time <- "Wk1b"
      tissue <- "PB"
    } else if(x == "JCC212_SJCAR19-05relapse-2ndInfusionWk2_PB") {
      time <- "Wk2b"
      tissue <- "PB"
    } else if(x == "JCC212_SJCAR19_04_CD45neg_Wk2") {
      time <- "Wk2b"
      tissue <- "PB"
    } else if (x == "JCC212_SJCAR19_04_GMP19003") {
      time <- "GMP"
      tissue <- "GMP"
    }
    if(time == "PB" || time == "BM") {
      time <- temp[5]
    }
    # tissue
    if(grepl("GMP", x, fixed = TRUE)) {
      tissue <- "GMP"
    } else if(is.na(tissue) || tissue == "PreTrans" || tissue == "PreTransA" || tissue == "PreTransB") {
      tissue <- "PB"
    } else if(startsWith(tissue, "Wk") || endsWith(tissue, "mo")) {
      tissue <- temp[3]
    }
    if(startsWith(tissue, "0")) {
      tissue <- "PB"
    }
    
    return(paste(px, time, tissue, sep = "_"))
    
  }, USE.NAMES = TRUE)
  
  ### px, time, & tissue
  px_name <- sapply(lib_name, function(x) {
    strsplit(x, split = "_", fixed = TRUE)[[1]][1]
  }, USE.NAMES = FALSE)
  names(px_name) <- names(lib_name)
  time_name <- sapply(lib_name, function(x) {
    strsplit(x, split = "_", fixed = TRUE)[[1]][2]
  }, USE.NAMES = FALSE)
  names(time_name) <- names(lib_name)
  tissue_name <- sapply(lib_name, function(x) {
    strsplit(x, split = "_", fixed = TRUE)[[1]][3]
  }, USE.NAMES = FALSE)
  names(tissue_name) <- names(lib_name)
  
  ### assign library
  Seurat_Obj@meta.data$barcode_file_name <- ""
  Seurat_Obj@meta.data$library <- ""
  Seurat_Obj@meta.data$px <- ""
  Seurat_Obj@meta.data$time <- ""
  Seurat_Obj@meta.data$tissue <- ""
  unique_tags <- unique(Seurat_Obj@meta.data$barcode_tag)
  for(i in 1:length(unique_tags)) {
    
    ### find which library-tag shares the barcodes the most
    tag_idx <- which(Seurat_Obj@meta.data$barcode_tag == unique_tags[i])
    max_num <- 0
    max_lib <- ""
    for(j in 1:length(barcodes)) {
      ### get the number of shared barcodes
      current_num <- length(intersect(Seurat_Obj@meta.data$barcode_short[tag_idx],
                                      barcodes[[j]]$V1))
      ### get max_num
      if(current_num > max_num) {
        max_num <- current_num
        max_lib <- names(barcodes)[j]
      }
    }
    
    ### only if the shared barcode percentage is above 80%
    if(max_num > length(tag_idx)*0.8) {
      Seurat_Obj@meta.data$barcode_file_name[tag_idx] <- max_lib
      Seurat_Obj@meta.data$library[tag_idx] <- lib_name[max_lib]
      Seurat_Obj@meta.data$px[tag_idx] <- px_name[max_lib]
      Seurat_Obj@meta.data$time[tag_idx] <- time_name[max_lib]
      Seurat_Obj@meta.data$tissue[tag_idx] <- tissue_name[max_lib]
    }
    
  }
  
  
  ### get new TCR info file paths
  TCR_data_dirs <- list.files(path = TCR_dir, pattern = "filtered_contig_annotations.csv$",
                              full.names = TRUE, recursive = TRUE)
  
  ### load and combine the TCR data
  tcr <- vector("list", length = length(TCR_data_dirs))
  names(tcr) <- sapply(basename(TCR_data_dirs), function(x) {
    return(substr(x, 1, nchar(x)-36))
  }, USE.NAMES = FALSE)
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
    
    ### remove non-productive rows
    tcr_data <- tcr_data[which(tcr_data$productive == "True"),]
    
    ### order by "chain" so that TRA rows come first than TRB rows
    ### and secondly, order by "CDR3 Nucleotide" sequence
    tcr_data <- tcr_data[order(as.character(tcr_data$chain), as.character(tcr_data$cdr3)),]
    
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
    
    ### save the TCR data
    tcr[[i]] <- tcr_data
    
  }
  
  ### there is a typo "DMP". change it to "GMP"
  # grep("DMP", names(tcr), fixed = TRUE)
  names(tcr)[9] <- "1662792_JCC212_GMPdonor33"
  ### there is another typo "SJCAR" should be changed to "SJCAR19" for consistency
  names(tcr)[49] <- "1779067_JCC212_SJCAR19-06_GMP19028"
  ### there is another typo "3_mo". change it to "3mo"
  names(tcr)[71] <- "1879710_JCC212_SJCAR19-06_3mo_BM"
  
  ### refine tcr library name
  tcr_names <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    px <- temp[3]
    time <- temp[4]
    tissue <- temp[5]
    
    ### exceptions
    # px
    if(px == "SJCAR19-05relapse-2ndInfusionWk1") {
      px <- "SJCAR19-05"
    } else if(px == "SJCAR19-05relapse-2ndInfusionWk2") {
      px <- "SJCAR19-05"
    } else if(grepl("Donor", px, fixed = TRUE)) {
      px <- paste0("SJCAR19-", strsplit(px, split = "Donor", fixed = TRUE)[[1]][2])
    } else if(grepl("donor", px, fixed = TRUE)) {
      px <- paste0("SJCAR19-", strsplit(px, split = "donor", fixed = TRUE)[[1]][2])
    }
    # time
    if(grepl("GMP-redo", x, fixed = TRUE)) {
      time <- "GMP-redo"
    } else if(grepl("GMP", x, fixed = TRUE)) {
      time <- "GMP"
    }
    if(time == "PB" || time == "BM") {
      time <- temp[5]
    }
    if(x == "1977283_JCC212_SJCAR19-05relapse-2ndInfusionWk1_PB") {
      time <- "Wk1b"
      tissue <- "PB"
    } else if(x == "1977284_JCC212_SJCAR19-05relapse-2ndInfusionWk2_PB") {
      time <- "Wk2b"
      tissue <- "PB"
    } else if(x == "1757070_JCC212_HealthyDonor32_PreTrans") {
      time <- "PreTransb"
      tissue <- "PB"
    } else if(x == "1757071_JCC212_HealthyDonor33_PreTrans") {
      time <- "PreTransb"
      tissue <- "PB"
    }
    # tissue
    if(grepl("GMP", x, fixed = TRUE)) {
      tissue <- "GMP"
    } else if(is.na(tissue)) {
      tissue <- "PB"
    } else if(startsWith(tissue, "Wk") || endsWith(tissue, "mo")) {
      tissue <- temp[4]
    }
    
    return(paste(px, time, tissue, sep = "_"))
    
  }, USE.NAMES = TRUE)
  
  ### attach TCR info to the meta.data
  Seurat_Obj@meta.data$tcr_file_name <- NA
  Seurat_Obj@meta.data$raw_clonotype_id <- NA
  Seurat_Obj@meta.data$cdr3_aa <- NA
  Seurat_Obj@meta.data$cdr3_nt <- NA
  Seurat_Obj@meta.data$tcr_reads <- NA
  Seurat_Obj@meta.data$tcr_umis <- NA
  unique_libs <- unique(Seurat_Obj@meta.data$library)
  for(i in 1:length(unique_libs)) {
    
    ### find which library shares the barcodes the most
    lib_idx <- which(Seurat_Obj@meta.data$library == unique_libs[i])
    max_num <- 0
    max_lib <- ""
    for(j in 1:length(tcr)) {
      ### get the number of shared barcodes
      current_num <- length(intersect(Seurat_Obj@meta.data$barcode_short[lib_idx],
                                      tcr[[j]]$barcode))
      ### get max_num
      if(current_num > max_num) {
        max_num <- current_num
        max_lib <- names(tcr)[j]
      }
    }
    
    ### attach to TCR info
    shared_barcodes <- intersect(Seurat_Obj@meta.data$barcode_short[lib_idx],
                                 tcr[[max_lib]]$barcode)
    gex_idx <- intersect(which(Seurat_Obj@meta.data$barcode_short %in% shared_barcodes),
                         lib_idx)
    if(identical(Seurat_Obj@meta.data$barcode_short[gex_idx], shared_barcodes)) {
      Seurat_Obj@meta.data$tcr_file_name[gex_idx] <- max_lib
      Seurat_Obj@meta.data$raw_clonotype_id[gex_idx] <- tcr[[max_lib]][shared_barcodes,"raw_clonotype_id"]
      Seurat_Obj@meta.data$cdr3_aa[gex_idx] <- tcr[[max_lib]][shared_barcodes,"cdr3_aa"]
      Seurat_Obj@meta.data$cdr3_nt[gex_idx] <- tcr[[max_lib]][shared_barcodes,"cdr3_nt"]
      Seurat_Obj@meta.data$tcr_reads[gex_idx] <- tcr[[max_lib]][shared_barcodes,"tcr_reads"]
      Seurat_Obj@meta.data$tcr_umis[gex_idx] <- tcr[[max_lib]][shared_barcodes,"tcr_umis"]
    } else {
      stop("ERROR: Barcodes do not match")
    }
    
    ### print some warnings
    if((tcr_names[max_lib] != unique_libs[i])) {
      writeLines(paste0("-\nWARNING: Different libraries matched between GEX and TCR\n",
                        "GEX: ", unique_libs[i], "\n",
                        "TCR: ", tcr_names[max_lib], "\n-"))
    }
    if(max_num < length(lib_idx)*0.5) {
      writeLines(paste("WARNING:", unique_libs[i], "(GEX) -", tcr_names[max_lib], "(TCR) :",
                       max_num, "< 0.5 *", length(lib_idx)))
    }
    
  }
  
  ### global clonotyping for each patient
  Seurat_Obj@meta.data$clonotype_id_by_patient <- NA
  for(px in unique(Seurat_Obj@meta.data$px)) {
    ### print progress
    writeLines(paste(px))
    
    ### get patient indices
    px_idx <- which(Seurat_Obj@meta.data$px == px)
    
    ### patient indices that have TCR info
    px_idx <- intersect(px_idx, which(!is.na(Seurat_Obj@meta.data$cdr3_aa)))
    
    ### give the strict version of global clonotypes using all alpha & beta chains
    ### since the variables in each of the "cdr3_aa" are ordered,
    ### if the "cdr3_aa" are exactly the same, they are the same clonotype
    
    ### get unique and duplicated indicies
    unique_idx <- which(!duplicated(Seurat_Obj@meta.data$cdr3_aa[px_idx]))
    dup_idx <- which(duplicated(Seurat_Obj@meta.data$cdr3_aa[px_idx]))
    
    ### give clone ids
    Seurat_Obj@meta.data$clonotype_id_by_patient[px_idx[unique_idx]] <- paste0(px, "_Clone_", 1:length(unique_idx))
    for(i in px_idx[dup_idx]) {
      target_idx <- which(Seurat_Obj@meta.data$cdr3_aa[px_idx[unique_idx]] == Seurat_Obj@meta.data$cdr3_aa[i])
      Seurat_Obj@meta.data$clonotype_id_by_patient[i] <- Seurat_Obj@meta.data$clonotype_id_by_patient[px_idx[unique_idx][target_idx]]
    }
  }
  
  #
  ### global clonotyping regardless of patients
  #
  ### (across all the patients)
  Seurat_Obj@meta.data$global_clonotype_id <- NA
  
  ### get unique and duplicated indicies
  unique_idx <- intersect(which(!duplicated(Seurat_Obj@meta.data$cdr3_aa)),
                          which(!is.na(Seurat_Obj@meta.data$cdr3_aa)))
  dup_idx <- intersect(which(duplicated(Seurat_Obj@meta.data$cdr3_aa)),
                       which(!is.na(Seurat_Obj@meta.data$cdr3_aa)))
  
  ### give clone ids
  Seurat_Obj@meta.data$global_clonotype_id[unique_idx] <- paste0("Global_Clone_", 1:length(unique_idx))
  for(i in 1:length(dup_idx)) {
    target_idx <- which(Seurat_Obj@meta.data$cdr3_aa[unique_idx] == Seurat_Obj@meta.data$cdr3_aa[dup_idx[i]])
    Seurat_Obj@meta.data$global_clonotype_id[dup_idx[i]] <- Seurat_Obj@meta.data$global_clonotype_id[unique_idx[target_idx]]
    
    if(i %% 1000 == 0) {
      writeLines(paste(i, "/", length(dup_idx)))
    }
  }
  
  
  ### UMAP by patients
  p <- DimPlot(Seurat_Obj, reduction = "umap", group.by = "time", split.by = "px", pt.size = 1, ncol = 3) +
    ggtitle("UMAP of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir, "/", "UMAP_Plot.png"), plot = p, width = 20, height = 12, dpi = 300)
  
  
  ###
  ### temporarily exclude Px5 - Wk1 cells for all the analyses
  ### since we noticed there are weirdly many CAR+ cells in there
  ###
  # Seurat_Obj <- SetIdent(object = Seurat_Obj,
  #                        cells = rownames(Seurat_Obj@meta.data),
  #                        value = Seurat_Obj@meta.data$library)
  # Seurat_Obj <- subset(Seurat_Obj, idents = unique(Seurat_Obj@meta.data$library)[-which(Seurat_Obj@meta.data$library == "SJCAR19-05_Wk1_PB")])
  # print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
  
  
  ### examine the proportion of CAR+ cells in each patient - barplot 
  plot_df <- data.frame(LIB=rep(unique(Seurat_Obj@meta.data$library), 2),
                        CAR=rep(c("CARpos", "CARneg"), length(unique(Seurat_Obj@meta.data$library))),
                        NUM=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:nrow(plot_df)) {
    plot_df$NUM[i] <- length(intersect(which(Seurat_Obj@meta.data$library == plot_df$LIB[i]),
                                       which(Seurat_Obj@meta.data$CAR == plot_df$CAR[i])))
  }
  ggplot(data=plot_df, aes_string(x="LIB", y="NUM", fill="CAR", label="NUM")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Cell Numbers - CARpos vs CARneg") +
    geom_text(size = 3, position = position_stack(vjust = 1)) +
    coord_flip() +
    scale_y_continuous(expand = c(0,300)) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 10),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir, "/", "CAR_Distribution_In_Each_Library.png"), width = 20, height = 10, dpi = 300)
  
  
  ### Alluvial plot - visualization of the lineage tracing
  
  ### theme that draws dotted lines for each y-axis ticks
  ### this function is from "immunarch" package
  theme_cleveland2 <- function(rotate = TRUE) {
    if (rotate) {
      theme(
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(
          colour = "grey70",
          linetype = "dashed"
        )
      )
    }
    else {
      theme(
        panel.grid.major.x = element_line(
          colour = "grey70",
          linetype = "dashed"
        ), panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )
    }
  }
  
  ### for each patient, make a statistics table
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    
    ### print progress
    writeLines(paste(patient))
    
    ### indicies that assign to the given patient
    px_idx <- which(Seurat_Obj@meta.data$px == patient)
    
    ### remove NA, "NA", and ""
    px_idx <- intersect(px_idx, intersect(intersect(which(!is.na(Seurat_Obj@meta.data$cdr3_aa)),
                                                    which(Seurat_Obj@meta.data$cdr3_aa != "NA")),
                                          which(Seurat_Obj@meta.data$cdr3_aa != "")))
    
    if(length(px_idx) > 0) {
      
      ### output directory for each patient
      outputDir2 <- paste0(outputDir, "/", patient, "/")
      dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
      
      ### set target indicies
      target_idx <- intersect(px_idx, which(!is.na(Seurat_Obj@meta.data$clonotype_id_by_patient)))
      
      ### get unique clonotypes
      dups <- Seurat_Obj@meta.data$clonotype_id_by_patient[target_idx][which(duplicated(Seurat_Obj@meta.data$clonotype_id_by_patient[target_idx]))]
      unique_clonotypes <- unique(dups)
      unique_clonotypes <- unique_clonotypes[intersect(intersect(which(!is.na(unique_clonotypes)),
                                                                 which(unique_clonotypes != "NA")),
                                                       which(unique_clonotypes != ""))]
      
      ### get unique time points
      unique_time_points <- unique(intersect(total_time_points, Seurat_Obj@meta.data$time[target_idx]))
      
      ### empty frequency table
      frequency_over_time <- matrix(0, length(unique_clonotypes), length(unique_time_points))
      rownames(frequency_over_time) <- unique_clonotypes
      colnames(frequency_over_time) <- unique_time_points
      
      ### empty CAR+ frequency table
      car_frequency_over_time <- matrix(0, length(unique_clonotypes), length(unique_time_points))
      rownames(car_frequency_over_time) <- unique_clonotypes
      colnames(car_frequency_over_time) <- unique_time_points
      
      ### empty frequency proportion table
      proportion_over_time <- matrix(0, length(unique_clonotypes), length(unique_time_points))
      rownames(proportion_over_time) <- unique_clonotypes
      colnames(proportion_over_time) <- unique_time_points
      
      ### empty CAR+ frequency proportion table
      car_proportion_over_time <- matrix(0, length(unique_clonotypes), length(unique_time_points))
      rownames(car_proportion_over_time) <- unique_clonotypes
      colnames(car_proportion_over_time) <- unique_time_points
      
      ### start time
      start_time <- Sys.time()
      
      ### set progress bar
      pb <- txtProgressBar(min = 0, max = length(unique_clonotypes)*length(unique_time_points), style = 3)
      
      ### fill out the frequency tables
      cnt <- 0
      for(time in unique_time_points) {
        for(clonotype in unique_clonotypes) {
          ### indicies for the specific patient, time, and clonotype
          target_idx <- intersect(intersect(px_idx,
                                            which(Seurat_Obj@meta.data$time == time)),
                                  which(Seurat_Obj@meta.data$clonotype_id_by_patient == clonotype))
          
          frequency_over_time[clonotype,time] <- length(target_idx)
          car_frequency_over_time[clonotype,time] <- length(intersect(target_idx,
                                                                      which(Seurat_Obj@meta.data$CAR == "CARpos")))
          
          ### update the progress bar
          cnt <- cnt + 1
          setTxtProgressBar(pb, cnt)
        }
        
        ### fill out the proportion tables
        cells_time_num <- length(intersect(px_idx,
                                           which(Seurat_Obj@meta.data$time == time)))
        car_cells_time_num <- length(intersect(intersect(px_idx,
                                                         which(Seurat_Obj@meta.data$time == time)),
                                               which(Seurat_Obj@meta.data$CAR == "CARpos")))
        if(car_cells_time_num > 0) {
          proportion_over_time[,time] <- signif(100*frequency_over_time[,time]/cells_time_num, digits = 4)
          car_proportion_over_time[,time] <- signif(100*car_frequency_over_time[,time]/car_cells_time_num, digits = 4)
        } else if(cells_time_num > 0) {
          proportion_over_time[,time] <- signif(100*frequency_over_time[,time]/cells_time_num, digits = 4)
        } else {
          proportion_over_time[,time] <- 0
          car_proportion_over_time[,time] <- 0
        }
      }
      close(pb)
      
      ### end time
      end_time <- Sys.time()
      
      ### print out the running time
      cat(paste("Running Time:",
                signif(as.numeric(difftime(end_time, start_time, units = "mins")), digits = 3),
                "mins"))
      
      ### add total column
      frequency_over_time <- data.frame(Clonotype=unique_clonotypes,
                                        frequency_over_time,
                                        Total=apply(frequency_over_time, 1, sum),
                                        stringsAsFactors = FALSE, check.names = FALSE)
      car_frequency_over_time <- data.frame(Clonotype=unique_clonotypes,
                                            car_frequency_over_time,
                                            Total=apply(car_frequency_over_time, 1, sum),
                                            stringsAsFactors = FALSE, check.names = FALSE)
      proportion_over_time <- data.frame(Clonotype=unique_clonotypes,
                                         proportion_over_time,
                                         Total=apply(proportion_over_time, 1, sum),
                                         stringsAsFactors = FALSE, check.names = FALSE)
      car_proportion_over_time <- data.frame(Clonotype=unique_clonotypes,
                                             car_proportion_over_time,
                                             Total=apply(car_proportion_over_time, 1, sum),
                                             stringsAsFactors = FALSE, check.names = FALSE)
      
      ### order the data frames by the Total column
      frequency_over_time <- frequency_over_time[order(-frequency_over_time[,"Total"]),]
      proportion_over_time <- proportion_over_time[order(-proportion_over_time[,"Total"]),]
      car_frequency_over_time <- car_frequency_over_time[order(-car_frequency_over_time[,"Total"]),]
      car_proportion_over_time <- car_proportion_over_time[order(-car_proportion_over_time[,"Total"]),]
      
      ### trim the car tables
      car_proportion_over_time <- car_proportion_over_time[which(car_frequency_over_time[,"Total"] > 1),]
      car_frequency_over_time <- car_frequency_over_time[which(car_frequency_over_time[,"Total"] > 1),]
      
      ### save the tables in Excel format
      write.xlsx2(frequency_over_time, file = paste0(outputDir2, "clonotype_frequency_over_time_", patient, ".xlsx"),
                  sheetName = "Clonotype_Frequency", row.names = FALSE)
      gc()
      write.xlsx2(proportion_over_time, file = paste0(outputDir2, "clonotype_proportion_over_time_", patient, ".xlsx"),
                  sheetName = "Clonotype_Proportion", row.names = FALSE)
      gc()
      write.xlsx2(car_frequency_over_time, file = paste0(outputDir2, "car_clonotype_frequency_over_time_", patient, ".xlsx"),
                  sheetName = "CARpos_Clonotype_Frequency", row.names = FALSE)
      gc()
      write.xlsx2(car_proportion_over_time, file = paste0(outputDir2, "car_clonotype_proportion_over_time_", patient, ".xlsx"),
                  sheetName = "CARpos_Clonotype_Proportion", row.names = FALSE)
      gc()
      
    }
    
  }
  
  ### draw an alluvial plot for each patient
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    
    ### print progress
    writeLines(paste(patient))
    
    ### load the file
    target_file <- read.xlsx2(file = paste0(outputDir, patient, "/car_clonotype_frequency_over_time_", patient, ".xlsx"),
                              sheetIndex = 1, stringsAsFactors = FALSE, check.names = FALSE,
                              row.names = 1)
    
    ### numerize the table
    for(i in 1:ncol(target_file)) {
      target_file[,i] <- as.numeric(target_file[,i])
    }
    
    ### remove all zero time points
    time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
    
    ### get time points except the Total
    time_points <- setdiff(time_points, c("Total"))
    
    if(length(time_points) > 1) {
      ###  get lineages
      lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
        return(length(which(x > 0)) > 1)  
      })),time_points]
      
      ### get an input data frame for the alluvial plot
      total_rows <- length(which(lineage_table[,time_points] > 0))
      if(total_rows > 0) {
        plot_df <- data.frame(Time=rep("", total_rows),
                              Clone_Size=rep(0, total_rows),
                              Clone=rep("", total_rows),
                              CDR3=rep("", total_rows))
        cnt <- 1
        for(i in 1:nrow(lineage_table)) {
          for(tp in time_points) {
            if(lineage_table[i,tp] > 0) {
              plot_df[cnt,] <- c(tp,
                                 lineage_table[i,tp],
                                 rownames(lineage_table)[i],
                                 "CDR3")
              cnt <- cnt + 1
            }
          }
        }
        plot_df$Time <- factor(plot_df$Time, levels = time_points)
        
        ### numerize the clone_size column
        plot_df$Clone_Size <- as.numeric(plot_df$Clone_Size)
        
        ### draw the alluvial plot (GMP - GMP-redo)
        ggplot(plot_df,
               aes(x = Time, stratum = Clone, alluvium = Clone,
                   y = Clone_Size,
                   fill = Clone, label = Clone)) +
          ggtitle(paste("Clonal Tracing in the CAR+ cells of", patient)) +
          geom_flow() +
          geom_stratum(alpha = 1) +
          geom_text(stat = "stratum", size = 2) +
          rotate_x_text(90) +
          theme_pubr(legend = "none") +
          theme(axis.title.x = element_blank()) +
          theme_cleveland2() +
          scale_fill_viridis(discrete = T) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
        ggsave(file = paste0(outputDir, patient, "/Car+_Clonal_Tracing_", patient, ".png"), width = 20, height = 10, dpi = 300)
        
        ### 1-on-1 alluvial plot with GMP
        if(length(grep("^GMP$", time_points)) == 1) {
          p <- vector("list", length = length(time_points)-1)
          names(p) <- time_points[-grep("^GMP$", time_points)]
          for(tp in names(p)) {
            p[[tp]] <- ggplot(plot_df[union(which(plot_df$Time == "GMP"),
                                            which(plot_df$Time == tp)),],
                              aes(x = Time, stratum = Clone, alluvium = Clone,
                                  y = Clone_Size,
                                  fill = Clone, label = Clone)) +
              geom_flow() +
              geom_stratum(alpha = 1) +
              geom_text(stat = "stratum", size = 2) +
              rotate_x_text(90) +
              theme_pubr(legend = "none") +
              theme(axis.title.x = element_blank()) +
              theme_cleveland2() +
              scale_fill_viridis(discrete = T) +
              scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          }
          
          ### arrange the plots and save
          fName <- paste0("Car+_Clonal_Tracing_", patient, "_1-on-1")
          if(length(p) == 1) {
            rowNum <- 1
            colNum <- 1
          } else if(length(p) == 2) {
            rowNum <- 1
            colNum <- 2
          } else if(length(p) < 5) {
            rowNum <- 2
            colNum <- 2
          } else if(length(p) < 7) {
            rowNum <- 3
            colNum <- 2
          } else {
            rowNum <- ceiling(sqrt(length(p)))
            colNum <- ceiling(sqrt(length(p)))
          }
          g <- arrangeGrob(grobs = p,
                           nrow = rowNum,
                           ncol = colNum,
                           top = fName)
          ggsave(file = paste0(outputDir, patient, "/", fName, ".png"), g, width = 20, height = 12, dpi = 300)
        }
        
        ### 1-on-1 alluvial plot with GMP-redo
        if(length(grep("^GMP-redo$", time_points)) == 1) {
          p <- vector("list", length = length(time_points)-1)
          names(p) <- time_points[-grep("^GMP-redo$", time_points)]
          for(tp in names(p)) {
            p[[tp]] <- ggplot(plot_df[union(which(plot_df$Time == "GMP-redo"),
                                            which(plot_df$Time == tp)),],
                              aes(x = Time, stratum = Clone, alluvium = Clone,
                                  y = Clone_Size,
                                  fill = Clone, label = Clone)) +
              geom_flow() +
              geom_stratum(alpha = 1) +
              geom_text(stat = "stratum", size = 2) +
              rotate_x_text(90) +
              theme_pubr(legend = "none") +
              theme(axis.title.x = element_blank()) +
              theme_cleveland2() +
              scale_fill_viridis(discrete = T) +
              scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          }
          
          ### arrange the plots and save
          fName <- paste0("Car+_Clonal_Tracing_", patient, "_1-on-1(2)")
          if(length(p) == 1) {
            rowNum <- 1
            colNum <- 1
          } else if(length(p) == 2) {
            rowNum <- 1
            colNum <- 2
          } else if(length(p) < 5) {
            rowNum <- 2
            colNum <- 2
          } else if(length(p) < 7) {
            rowNum <- 3
            colNum <- 2
          } else {
            rowNum <- ceiling(sqrt(length(p)))
            colNum <- ceiling(sqrt(length(p)))
          }
          g <- arrangeGrob(grobs = p,
                           nrow = rowNum,
                           ncol = colNum,
                           top = fName)
          ggsave(file = paste0(outputDir, patient, "/", fName, ".png"), g, width = 20, height = 12, dpi = 300)
        }
        
        ### draw the alluvial plot (GMP-redo - GMP)
        if(length(grep("GMP", time_points)) == 2) {
          gmp_idx <- which(time_points == "GMP")
          gmp_redo_idx <- which(time_points == "GMP-redo")
          time_points[gmp_idx] <- "GMP-redo"
          time_points[gmp_redo_idx] <- "GMP"
          plot_df$Time <- factor(plot_df$Time, levels = time_points)
          ggplot(plot_df,
                 aes(x = Time, stratum = Clone, alluvium = Clone,
                     y = Clone_Size,
                     fill = Clone, label = Clone)) +
            ggtitle(paste("Clonal Tracing in the CAR+ cells of", patient)) +
            geom_flow() +
            geom_stratum(alpha = 1) +
            geom_text(stat = "stratum", size = 2) +
            rotate_x_text(90) +
            theme_pubr(legend = "none") +
            theme(axis.title.x = element_blank()) +
            theme_cleveland2() +
            scale_fill_viridis(discrete = T) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          ggsave(file = paste0(outputDir, patient, "/Car+_Clonal_Tracing_", patient, "(2).png"), width = 20, height = 10, dpi = 300)
        }
      }
      
      gc()
    }
    
    ### NOW IT'S NOT CAR+ ONLY BUT WITH ALL THE TCR CELLS!
    ### load the file
    target_file <- read.xlsx2(file = paste0(outputDir, patient, "/clonotype_frequency_over_time_", patient, ".xlsx"),
                              sheetIndex = 1, stringsAsFactors = FALSE, check.names = FALSE,
                              row.names = 1)
    
    ### numerize the table
    for(i in 1:ncol(target_file)) {
      target_file[,i] <- as.numeric(target_file[,i])
    }
    
    ### remove all zero time points
    time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
    
    ### get time points except the Total
    time_points <- setdiff(time_points, c("Total"))
    
    if(length(time_points) > 1) {
      ###  get lineages
      lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
        return(length(which(x > 0)) > 1)  
      })),time_points]
      
      ### get an input data frame for the alluvial plot
      total_rows <- length(which(lineage_table[,time_points] > 0))
      if(total_rows > 0) {
        plot_df <- data.frame(Time=rep("", total_rows),
                              Clone_Size=rep(0, total_rows),
                              Clone=rep("", total_rows),
                              CDR3=rep("", total_rows))
        cnt <- 1
        for(i in 1:nrow(lineage_table)) {
          for(tp in time_points) {
            if(lineage_table[i,tp] > 0) {
              plot_df[cnt,] <- c(tp,
                                 lineage_table[i,tp],
                                 rownames(lineage_table)[i],
                                 "CDR3")
              cnt <- cnt + 1
            }
          }
        }
        plot_df$Time <- factor(plot_df$Time, levels = time_points)
        
        ### numerize the clone_size column
        plot_df$Clone_Size <- as.numeric(plot_df$Clone_Size)
        
        ### draw the alluvial plot
        ggplot(plot_df,
               aes(x = Time, stratum = Clone, alluvium = Clone,
                   y = Clone_Size,
                   fill = Clone, label = Clone)) +
          ggtitle(paste("Clonal Tracing with the TCR of", patient)) +
          geom_flow() +
          geom_stratum(alpha = 1) +
          rotate_x_text(90) +
          theme_pubr(legend = "none") +
          theme(axis.title.x = element_blank()) +
          theme_cleveland2() +
          scale_fill_viridis(discrete = T) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
        ggsave(file = paste0(outputDir, patient, "/Clonal_Tracing_", patient, ".png"), width = 20, height = 10, dpi = 300)
        
        ### 1-on-1 alluvial plot with GMP
        if(length(grep("^GMP$", time_points)) == 1) {
          p <- vector("list", length = length(time_points)-1)
          names(p) <- time_points[-grep("^GMP$", time_points)]
          for(tp in names(p)) {
            p[[tp]] <- ggplot(plot_df[union(which(plot_df$Time == "GMP"),
                                            which(plot_df$Time == tp)),],
                              aes(x = Time, stratum = Clone, alluvium = Clone,
                                  y = Clone_Size,
                                  fill = Clone, label = Clone)) +
              geom_flow() +
              geom_stratum(alpha = 1) +
              geom_text(stat = "stratum", size = 2) +
              rotate_x_text(90) +
              theme_pubr(legend = "none") +
              theme(axis.title.x = element_blank()) +
              theme_cleveland2() +
              scale_fill_viridis(discrete = T) +
              scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          }
          
          ### arrange the plots and save
          fName <- paste0("Clonal_Tracing_", patient, "_1-on-1")
          if(length(p) == 1) {
            rowNum <- 1
            colNum <- 1
          } else if(length(p) == 2) {
            rowNum <- 1
            colNum <- 2
          } else if(length(p) < 5) {
            rowNum <- 2
            colNum <- 2
          } else if(length(p) < 7) {
            rowNum <- 3
            colNum <- 2
          } else {
            rowNum <- ceiling(sqrt(length(p)))
            colNum <- ceiling(sqrt(length(p)))
          }
          g <- arrangeGrob(grobs = p,
                           nrow = rowNum,
                           ncol = colNum,
                           top = fName)
          ggsave(file = paste0(outputDir, patient, "/", fName, ".png"), g, width = 20, height = 12, dpi = 300)
        }
        
        ### 1-on-1 alluvial plot with GMP-redo
        if(length(grep("^GMP-redo$", time_points)) == 1) {
          p <- vector("list", length = length(time_points)-1)
          names(p) <- time_points[-grep("^GMP-redo$", time_points)]
          for(tp in names(p)) {
            p[[tp]] <- ggplot(plot_df[union(which(plot_df$Time == "GMP-redo"),
                                            which(plot_df$Time == tp)),],
                              aes(x = Time, stratum = Clone, alluvium = Clone,
                                  y = Clone_Size,
                                  fill = Clone, label = Clone)) +
              geom_flow() +
              geom_stratum(alpha = 1) +
              geom_text(stat = "stratum", size = 2) +
              rotate_x_text(90) +
              theme_pubr(legend = "none") +
              theme(axis.title.x = element_blank()) +
              theme_cleveland2() +
              scale_fill_viridis(discrete = T) +
              scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          }
          
          ### arrange the plots and save
          fName <- paste0("Clonal_Tracing_", patient, "_1-on-1(2)")
          if(length(p) == 1) {
            rowNum <- 1
            colNum <- 1
          } else if(length(p) == 2) {
            rowNum <- 1
            colNum <- 2
          } else if(length(p) < 5) {
            rowNum <- 2
            colNum <- 2
          } else if(length(p) < 7) {
            rowNum <- 3
            colNum <- 2
          } else {
            rowNum <- ceiling(sqrt(length(p)))
            colNum <- ceiling(sqrt(length(p)))
          }
          g <- arrangeGrob(grobs = p,
                           nrow = rowNum,
                           ncol = colNum,
                           top = fName)
          ggsave(file = paste0(outputDir, patient, "/", fName, ".png"), g, width = 20, height = 12, dpi = 300)
        }
        
        ### draw the alluvial plot (GMP-redo - GMP)
        if(length(grep("GMP", time_points)) == 2) {
          gmp_idx <- which(time_points == "GMP")
          gmp_redo_idx <- which(time_points == "GMP-redo")
          time_points[gmp_idx] <- "GMP-redo"
          time_points[gmp_redo_idx] <- "GMP"
          plot_df$Time <- factor(plot_df$Time, levels = time_points)
          ggplot(plot_df,
                 aes(x = Time, stratum = Clone, alluvium = Clone,
                     y = Clone_Size,
                     fill = Clone, label = Clone)) +
            ggtitle(paste("Clonal Tracing with the TCR of", patient)) +
            geom_flow() +
            geom_stratum(alpha = 1) +
            rotate_x_text(90) +
            theme_pubr(legend = "none") +
            theme(axis.title.x = element_blank()) +
            theme_cleveland2() +
            scale_fill_viridis(discrete = T) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          ggsave(file = paste0(outputDir, patient, "/Clonal_Tracing_", patient, "(2).png"), width = 20, height = 10, dpi = 300)
        }
      }
      
      gc()
    }
    
  }
  
  
  ### because we think there are too many CAR+ cells in patient 5 - Wk1,
  ### we want to see how the CAR counts of the cells are
  print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
  target_idx <- intersect(intersect(which(Seurat_Obj@meta.data$px == "SJCAR19-05"),
                                    which(Seurat_Obj@meta.data$time == "Wk1")),
                          which(Seurat_Obj@meta.data$CAR == "CARpos"))
  
  ### density plot
  ggplot(data.frame(CARcount=Seurat_Obj@assays$RNA@counts["JCC-SJCAR19short", target_idx]),
         aes(x = CARcount)) +
    geom_histogram(binwidth = 1) +
    ggtitle(paste("CAR counts of the Px5 Wk1 CAR+ Cells")) +
    theme_classic(base_size = 16)
  ggsave(file = paste0(outputDir, "/CAR_Counts_Px5_Wk1.png"), width = 8, height = 6, dpi = 300)
  
  
  ### make a clonotype lineage table RDS file
  SJCAR19_Clonotype_Frequency <- vector("list", length = 2)
  names(SJCAR19_Clonotype_Frequency) <- c("ALL", "CARPOSONLY")
  SJCAR19_Clonotype_Frequency[[1]] <- vector("list", length = length(unique(Seurat_Obj@meta.data$px)))
  names(SJCAR19_Clonotype_Frequency[[1]]) <- unique(Seurat_Obj@meta.data$px)
  SJCAR19_Clonotype_Frequency[[2]] <- vector("list", length = length(unique(Seurat_Obj@meta.data$px)))
  names(SJCAR19_Clonotype_Frequency[[2]]) <- unique(Seurat_Obj@meta.data$px)
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    ### print progress
    writeLines(paste(patient))
    
    ### load the file
    SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[patient]] <- read.xlsx2(file = paste0(outputDir, patient, "/car_clonotype_frequency_over_time_", patient, ".xlsx"),
                                                                         sheetIndex = 1, stringsAsFactors = FALSE, check.names = FALSE,
                                                                         row.names = 1)
    
    ### numerize the table
    for(i in 1:ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[patient]])) {
      SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[patient]][,i] <- as.numeric(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[patient]][,i])
    }
    
    ### NOW IT'S NOT CAR+ ONLY BUT WITH ALL THE TCR CELLS!
    ### load the file
    SJCAR19_Clonotype_Frequency[["ALL"]][[patient]] <- read.xlsx2(file = paste0(outputDir, patient, "/clonotype_frequency_over_time_", patient, ".xlsx"),
                                                                  sheetIndex = 1, stringsAsFactors = FALSE, check.names = FALSE,
                                                                  row.names = 1)
    
    ### numerize the table
    for(i in 1:ncol(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])) {
      SJCAR19_Clonotype_Frequency[["ALL"]][[patient]][,i] <- as.numeric(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]][,i])
    }
  }
  
  ### Save the clonotype lineage info
  saveRDS(SJCAR19_Clonotype_Frequency, file = paste0(dirname(Seurat_RObj_path), "/SJCAR19_Clonotype_Lineages.RDS" ))
  
  
  # ******************************************************************************************
  # Pathway Analysis with clusterProfiler package
  # Input: geneList     = a vector of gene Entrez IDs for pathway analysis [numeric or character]
  #        org          = organism that will be used in the analysis ["human" or "mouse"]
  #                       should be either "human" or "mouse"
  #        database     = pathway analysis database (KEGG or GO) ["KEGG" or "GO"]
  #        title        = title of the pathway figure [character]
  #        pv_threshold = pathway analysis p-value threshold (not DE analysis threshold) [numeric]
  #        displayNum   = the number of pathways that will be displayed [numeric]
  #                       (If there are many significant pathways show the few top pathways)
  #        imgPrint     = print a plot of pathway analysis [TRUE/FALSE]
  #        dir          = file directory path of the output pathway figure [character]
  #
  # Output: Pathway analysis results in figure - using KEGG and GO pathways
  #         The x-axis represents the number of DE genes in the pathway
  #         The y-axis represents pathway names
  #         The color of a bar indicates adjusted p-value from the pathway analysis
  #         For Pathview Result, all colored genes are found DE genes in the pathway,
  #         and the color indicates log2(fold change) of the DE gene from DE analysis
  # ******************************************************************************************
  pathwayAnalysis_CP <- function(geneList,
                                 org,
                                 database,
                                 title="Pathway_Results",
                                 pv_threshold=0.05,
                                 displayNum=Inf,
                                 imgPrint=TRUE,
                                 dir="./") {
    
    ### load library
    if(!require(clusterProfiler, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("clusterProfiler")
      require(clusterProfiler, quietly = TRUE)
    }
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    
    ### collect gene list (Entrez IDs)
    geneList <- geneList[which(!is.na(geneList))]
    
    if(!is.null(geneList)) {
      ### make an empty list
      p <- list()
      
      if(database == "KEGG") {
        ### KEGG Pathway
        kegg_enrich <- enrichKEGG(gene = geneList, organism = org, pvalueCutoff = pv_threshold)
        
        if(is.null(kegg_enrich)) {
          writeLines("KEGG Result does not exist")
          return(NULL)
        } else {
          kegg_enrich@result <- kegg_enrich@result[which(kegg_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(kegg_enrich@result) <= displayNum)) {
              result <- kegg_enrich@result
              description <- kegg_enrich@result$Description
            } else {
              result <- kegg_enrich@result[1:displayNum,]
              description <- kegg_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(kegg_enrich) > 0) {
              p[[1]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title)) +
                theme(axis.text = element_text(size = 25))
              
              png(paste0(dir, "kegg_", title, ".png"), width = 2000, height = 1000)
              print(p[[1]])
              dev.off()
            } else {
              writeLines("KEGG Result does not exist")
            }
          }
          
          return(kegg_enrich@result)
        }
      } else if(database == "GO") {
        ### GO Pathway
        if(org == "human") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else if(org == "mouse") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else {
          go_enrich <- NULL
          writeLines(paste("Unknown org variable:", org))
        }
        
        if(is.null(go_enrich)) {
          writeLines("GO Result does not exist")
          return(NULL)
        } else {
          go_enrich@result <- go_enrich@result[which(go_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(go_enrich@result) <= displayNum)) {
              result <- go_enrich@result
              description <- go_enrich@result$Description
            } else {
              result <- go_enrich@result[1:displayNum,]
              description <- go_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(go_enrich) > 0) {
              p[[2]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("GO ", title)) +
                theme(axis.text = element_text(size = 25))
              
              png(paste0(dir, "go_", title, ".png"), width = 2000, height = 1000)
              print(p[[2]])
              dev.off()
            } else {
              writeLines("GO Result does not exist")
            }
          }
          
          return(go_enrich@result)
        }
      } else {
        stop("database prameter should be \"GO\" or \"KEGG\"")
      }
    } else {
      writeLines("geneList = NULL")
    }
  }
  
  
  ### persisters vs non-persisters
  ### DE + pathway analysis
  ### With all the patients
  #
  # load clonotype lineages info
  SJCAR19_Clonotype_Frequency <- readRDS(clonotype_lineage_info_path)
  
  ### GMP CAR+ persistent clones
  pClones <- NULL
  for(i in 1:length(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]])) {
    gmp_idx <- which(colnames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) == "GMP")
    gmp_redo_idx <- which(colnames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) == "GMP-redo")
    last_gmp_idx <- max(gmp_idx, gmp_redo_idx)
    
    ### if at least GMP or GMP-redo exist and there are at least one afterward-time point
    if((last_gmp_idx != -Inf) && (ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) - last_gmp_idx > 1)) {
      ### collect persistent clones that appeared in GMP and persist afterwards
      if(nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) > 0) {
        for(j in 1:nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])) {
          for(k in last_gmp_idx:(ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])-1)) {
            if((SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,"GMP"] > 0 ||
                SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,"GMP-redo"] > 0) &&
               SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,k] > 0) {
              pClones <- c(pClones, rownames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])[j])
              break;
            }
          }
        }
      }
    }
  }
  
  ### GMP CAR+ persistent cells
  pIdx <- intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient %in% pClones),
                    intersect(union(which(Seurat_Obj@meta.data$time == "GMP"),
                                    which(Seurat_Obj@meta.data$time == "GMP-redo")),
                              which(Seurat_Obj@meta.data$CAR == "CARpos")))
  
  ### GMP CAR+ non-persistent cells
  npIdx <- setdiff(intersect(union(which(Seurat_Obj@meta.data$time == "GMP"),
                                   which(Seurat_Obj@meta.data$time == "GMP-redo")),
                             which(Seurat_Obj@meta.data$CAR == "CARpos")),
                   which(Seurat_Obj@meta.data$clonotype_id_by_patient %in% pClones))
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  ### annotate GMP CAR+ persisters
  Seurat_Obj@meta.data$GMP_CARpos_Persister <- NA
  Seurat_Obj@meta.data$GMP_CARpos_Persister[pIdx] <- "YES"
  Seurat_Obj@meta.data$GMP_CARpos_Persister[npIdx] <- "NO"
  
  ### because there are too many cells in non-persisters
  ### randomly select some from those and perform DE analysis
  set.seed(1234)
  Seurat_Obj@meta.data$GMP_CARpos_Persister[sample(npIdx, length(npIdx) - length(pIdx))] <- NA
  
  ### set ident of the Seurat object with persister info
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$GMP_CARpos_Persister)
  
  ### DE analysis between persisters vs non-persisters
  de_result <- FindMarkers(Seurat_Obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/GMP_CARpos_Persisters_vs_NonPersisters.xlsx"),
              sheetName = "GMP_CARpos_DE_Result", row.names = FALSE)
  
  ### pathway analysis
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                            rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                            "ENTREZID", "SYMBOL"),
                                          org = "human", database = "GO",
                                          title = paste0("Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters"),
                                          displayNum = 30, imgPrint = TRUE,
                                          dir = paste0(outputDir))
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters"),
                                            displayNum = 30, imgPrint = TRUE,
                                            dir = paste0(outputDir))
  if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
    write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters.xlsx"),
                row.names = FALSE, sheetName = paste0("GO_Result"))
  }
  if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
    write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters.xlsx"),
                row.names = FALSE, sheetName = paste0("KEGG_Result"))
  }
  
  
  ### a function for color brewer
  cell_pal <- function(cell_vars, pal_fun) {
    if (is.numeric(cell_vars)) {
      pal <- pal_fun(100)
      return(pal[cut(cell_vars, breaks = 100)])
    } else {
      categories <- sort(unique(cell_vars))
      pal <- setNames(pal_fun(length(categories)), categories)
      return(pal[cell_vars])
    }
  }
  
  #' @title Plot Slingshot output
  #' @name plot-SlingshotDataSet
  #' @aliases plot-SlingshotDataSet plot,SlingshotDataSet,ANY-method
  #'
  #' @description Tools for visualizing lineages inferred by \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{"curves"}, or \code{"both"} (by partial matching),
  #'   see Details for more.
  #' @param linInd integer, an index indicating which lineages should be plotted
  #'   (default is to plot all lineages). If \code{col} is a vector, it will be
  #'   subsetted by \code{linInd}.
  #' @param show.constraints logical, whether or not the user-specified initial
  #'   and terminal clusters should be specially denoted by green and red dots,
  #'   respectively.
  #' @param add logical, indicates whether the output should be added to an
  #'   existing plot.
  #' @param dims numeric, which dimensions to plot (default is \code{1:2}).
  #' @param asp numeric, the y/x aspect ratio, see \code{\link{plot.window}}.
  #' @param cex numeric, amount by which points should be magnified, see
  #'   \code{\link{par}}.
  #' @param lwd numeric, the line width, see \code{\link{par}}.
  #' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
  #' @param ... additional parameters to be passed to \code{\link{lines}}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' plot(sds, type = 'b')
  #'
  #' # add to existing plot
  #' plot(rd, col = 'grey50')
  #' lines(sds, lwd = 3)
  #'
  #' @import graphics
  #' @import grDevices
  #' @export
  setMethod(
    f = "plot",
    signature = signature(x = "SlingshotDataSet"),
    definition = function(x, type = NULL,
                          linInd = NULL,
                          show.constraints = FALSE,
                          constraints.col = NULL,
                          add = FALSE,
                          dims = seq_len(2),
                          asp = 1,
                          cex = 2,
                          lwd = 2,
                          col = 1,
                          ...) {
      col <- rep(col, length(slingLineages(x)))
      curves <- FALSE
      lineages <- FALSE
      if(is.null(type)){
        if(length(slingCurves(x)) > 0){
          type <- 'curves'
        }else if(length(slingLineages(x)) > 0){
          type <- 'lineages'
        }else{
          stop('No lineages or curves detected.')
        }
      }else{
        type <- c('curves','lineages','both')[pmatch(type,
                                                     c('curves','lineages','both'))]
        if(is.na(type)){
          stop('Unrecognized type argument.')
        }
      }
      
      if(type %in% c('lineages','both')){
        lineages <- TRUE
      }
      if(type %in% c('curves','both')){
        curves <- TRUE
      }
      
      if(lineages & (length(slingLineages(x))==0)){
        stop('No lineages detected.')
      }
      if(curves & (length(slingCurves(x))==0)){
        stop('No curves detected.')
      }
      
      if(is.null(linInd)){
        linInd <- seq_along(slingLineages(x))
      }else{
        linInd <- as.integer(linInd)
        if(!all(linInd %in% seq_along(slingLineages(x)))){
          if(any(linInd %in% seq_along(slingLineages(x)))){
            linInd.removed <-
              linInd[! linInd %in% seq_along(slingLineages(x))]
            linInd <-
              linInd[linInd %in% seq_along(slingLineages(x))]
            message('Unrecognized lineage indices (linInd): ',
                    paste(linInd.removed, collapse = ", "))
          }else{
            stop('None of the provided lineage indices',
                 ' (linInd) were found.')
          }
        }
      }
      
      if(lineages){
        X <- reducedDim(x)
        clusterLabels <- slingClusterLabels(x)
        connectivity <- slingAdjacency(x)
        clusters <- rownames(connectivity)
        nclus <- nrow(connectivity)
        centers <- t(vapply(clusters,function(clID){
          w <- clusterLabels[,clID]
          return(apply(X, 2, weighted.mean, w = w))
        }, rep(0,ncol(X))))
        rownames(centers) <- clusters
        X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
        clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                       drop = FALSE]
        linC <- slingParams(x)
        clus2include <- unique(unlist(slingLineages(x)[linInd]))
      }
      
      if(!add){
        xs <- NULL
        ys <- NULL
        if(lineages){
          xs <- c(xs, centers[,dims[1]])
          ys <- c(ys, centers[,dims[2]])
        }
        if(curves){
          npoints <- nrow(slingCurves(x)[[1]]$s)
          xs <- c(xs, as.numeric(vapply(slingCurves(x),
                                        function(c){ c$s[,dims[1]] }, rep(0,npoints))))
          ys <- c(ys, as.numeric(vapply(slingCurves(x),
                                        function(c){ c$s[,dims[2]] }, rep(0,npoints))))
        }
        plot(x = NULL, y = NULL, asp = asp,
             xlim = range(xs), ylim = range(ys),
             xlab = colnames(reducedDim(x))[dims[1]],
             ylab = colnames(reducedDim(x))[dims[2]])
      }
      
      if(lineages){
        for(i in seq_len(nclus-1)){
          for(j in seq(i+1,nclus)){
            if(connectivity[i,j]==1 &
               all(clusters[c(i,j)] %in% clus2include)){
              lines(centers[c(i,j), dims],
                    lwd = lwd, col = col[1], ...)
            }
          }
        }
        points(centers[clusters %in% clus2include, dims],
               cex = cex, pch = 16, col = col[1])
        if(show.constraints && !is.null(constraints.col)){
          for(const in names(constraints.col)) {
            points(centers[clusters %in% const, dims,
                           drop=FALSE], cex = cex / 2,
                   col = constraints.col[const], pch = 16)
            text(x = centers[clusters %in% const, dims[1]]+0.3,
                 y = centers[clusters %in% const, dims[2]]+0.8,
                 labels = const,
                 cex = cex / 3,
                 col = "black")
          }
        }
      }
      if(curves){
        for(ii in seq_along(slingCurves(x))[linInd]){
          c <- slingCurves(x)[[ii]]
          lines(c$s[c$ord, dims], lwd = lwd, col = col[ii], ...)
        }
      }
      invisible(NULL)
    }
  )
  
  #' @title Pairs plot of Slingshot output
  #' @name pairs-SlingshotDataSet
  #'
  #' @description A tool for quickly visualizing lineages inferred by
  #'   \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
  #'   Details for more.
  #' @param show.constraints logical, whether or not the user-specified initial
  #'   and terminal clusters should be specially denoted by green and red dots,
  #'   respectively.
  #' @param col character, color vector for points.
  #' @param pch integer or character specifying the plotting symbol, see
  #'   \code{\link{par}}.
  #' @param cex numeric, amount by which points should be magnified, see
  #'   \code{\link{par}}.
  #' @param lwd numeric, the line width, see \code{\link{par}}.
  #' @param ... additional parameters for \code{plot} or \code{axis}, see
  #'   \code{\link{pairs}}.
  #' @param labels character, the names of the variables, see \code{\link{pairs}}.
  #' @param horInd see \code{\link{pairs}}.
  #' @param verInd see \code{\link{pairs}}.
  #' @param lower.panel see \code{\link{pairs}}.
  #' @param upper.panel see \code{\link{pairs}}.
  #' @param diag.panel see \code{\link{pairs}}.
  #' @param text.panel see \code{\link{pairs}}.
  #' @param label.pos see \code{\link{pairs}}.
  #' @param line.main see \code{\link{pairs}}.
  #' @param cex.labels see \code{\link{pairs}}.
  #' @param font.labels see \code{\link{pairs}}.
  #' @param row1attop see \code{\link{pairs}}.
  #' @param gap see \code{\link{pairs}}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' pairs(sds, type = 'curves')
  #'
  #' @export
  pairs.SlingshotDataSet <-
    function (x, type = NULL, show.constraints = FALSE, col = NULL,
              constraints.col = NULL,
              pch = 16, cex=1, lwd=2, ...,
              labels, horInd = seq_len(nc), verInd = seq_len(nc),
              lower.panel = FALSE, upper.panel = TRUE,
              diag.panel = NULL, text.panel = textPanel,
              label.pos = 0.5 + has.diag/3, line.main = 3,
              cex.labels = NULL, font.labels = 1,
              row1attop = TRUE, gap = 1,
              xlim=NULL, ylim=NULL) {
      #####
      lp.sling <- lower.panel
      up.sling <- upper.panel
      panel <- points
      if(!up.sling){
        upper.panel <- NULL
      }else{
        upper.panel <- panel
      }
      if(!lower.panel){
        lower.panel <- NULL
      }else{
        lower.panel <- panel
      }
      log = ""
      sds <- x
      x <- reducedDim(sds)
      curves <- FALSE
      lineages <- FALSE
      if(is.null(type)){
        if(length(slingCurves(sds)) > 0){
          type <- 'curves'
        }else if(length(slingLineages(sds)) > 0){
          type <- 'lineages'
        }else{
          stop('No lineages or curves detected.')
        }
      }else{
        type <- c('curves','lineages','both')[pmatch(type,
                                                     c('curves','lineages',
                                                       'both'))]
        if(is.na(type)){
          stop('Unrecognized type argument.')
        }
      }
      if(type %in% c('lineages','both')){
        lineages <- TRUE
      }
      if(type %in% c('curves','both')){
        curves <- TRUE
      }
      if(lineages & (length(slingLineages(sds))==0)){
        stop('No lineages detected.')
      }
      if(curves & (length(slingCurves(sds))==0)){
        stop('No curves detected.')
      }
      if(lineages){
        forest <- slingAdjacency(sds)
        clusters <- rownames(forest)
        nclus <- nrow(forest)
        centers <- t(vapply(clusters,function(clID){
          w <- slingClusterLabels(sds)[,clID]
          return(apply(x, 2, weighted.mean, w = w))
        }, rep(0,ncol(reducedDim(sds)))))
        rownames(centers) <- clusters
        linC <- slingParams(sds)
      }
      range.max <- max(apply(x,2,function(xi){
        r <- range(xi, na.rm = TRUE)
        return(abs(r[2] - r[1]))
      }))
      plot.ranges <- apply(x,2,function(xi){
        mid <- (max(xi,na.rm = TRUE) + min(xi,na.rm = TRUE))/2
        return(c(mid - range.max/2, mid + range.max/2))
      })
      if(is.null(col)){
        if(requireNamespace("RColorBrewer", quietly = TRUE)) {
          cc <- c(RColorBrewer::brewer.pal(9, "Set1")[-c(1,3,6)],
                  RColorBrewer::brewer.pal(7, "Set2")[-2],
                  RColorBrewer::brewer.pal(6, "Dark2")[-5],
                  RColorBrewer::brewer.pal(8, "Set3")[-c(1,2)])
        } else {
          cc <- seq_len(100)
        }
        col <- cc[apply(slingClusterLabels(sds),1,which.max)]
      }
      #####
      if(doText <- missing(text.panel) || is.function(text.panel))
        textPanel <-
        function(x = 0.5, y = 0.5, txt, cex, font)
          text(x, y, txt, cex = cex, font = font)
      
      localAxis <- function(side, x, y, xpd, bg, col=NULL, lwd=NULL, main,
                            oma, ...) {
        ## Explicitly ignore any color argument passed in as
        ## it was most likely meant for the data points and
        ## not for the axis.
        xpd <- NA
        if(side %% 2L == 1L && xl[j]) xpd <- FALSE
        if(side %% 2L == 0L && yl[i]) xpd <- FALSE
        if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
      }
      
      localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
      localLowerPanel <- function(..., main, oma, font.main, cex.main)
        lower.panel(...)
      localUpperPanel <- function(..., main, oma, font.main, cex.main)
        upper.panel(...)
      localDiagPanel <- function(..., main, oma, font.main, cex.main)
        diag.panel(...)
      
      dots <- list(...); nmdots <- names(dots)
      if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for(i in seq_along(names(x))) {
          if(is.factor(x[[i]]) || is.logical(x[[i]]))
            x[[i]] <- as.numeric(x[[i]])
          if(!is.numeric(unclass(x[[i]])))
            stop("non-numeric argument to 'pairs'")
        }
      } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")
      panel <- match.fun(panel)
      if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel <- match.fun(lower.panel)
      if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel <- match.fun(upper.panel)
      if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))
        diag.panel <- match.fun( diag.panel)
      
      if(row1attop) {
        tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
        tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
      }
      
      nc <- ncol(x)
      if (nc < 2L) stop("only one column in the argument to 'pairs'")
      if(!all(horInd >= 1L & horInd <= nc))
        stop("invalid argument 'horInd'")
      if(!all(verInd >= 1L & verInd <= nc))
        stop("invalid argument 'verInd'")
      if(doText) {
        if (missing(labels)) {
          labels <- colnames(x)
          if (is.null(labels)) labels <- paste("var", 1L:nc)
        }
        else if(is.null(labels)) doText <- FALSE
      }
      oma <- if("oma" %in% nmdots) dots$oma
      main <- if("main" %in% nmdots) dots$main
      if (is.null(oma))
        oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)
      opar <- par(mfrow = c(length(horInd), length(verInd)),
                  mar = rep.int(gap/2, 4), oma = oma)
      on.exit(par(opar))
      dev.hold(); on.exit(dev.flush(), add = TRUE)
      
      xl <- yl <- logical(nc)
      if (is.numeric(log)) xl[log] <- yl[log] <- TRUE
      else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}
      for (i in if(row1attop) verInd else rev(verInd))
        for (j in horInd) {
          l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
          
          if(is.null(xlim) & !is.null(ylim))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = plot.ranges[,j], ylim=ylim)
          else if(!is.null(xlim) & is.null(ylim))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = xlim, ylim = plot.ranges[,i])
          else if(!is.null(xlim) & !is.null(ylim))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = xlim, ylim=ylim)
          else
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = plot.ranges[,j], ylim = plot.ranges[,i])
          
          if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
            box()
            if(i == 1  && (!(j %% 2L) || !has.upper || !has.lower ))
              localAxis(1L + 2L*row1attop, x[, j], x[, i], ...)
            if(i == nc && (  j %% 2L  || !has.upper || !has.lower ))
              localAxis(3L - 2L*row1attop, x[, j], x[, i], ...)
            if(j == 1  && (!(i %% 2L) || !has.upper || !has.lower ))
              localAxis(2L, x[, j], x[, i], ...)
            if(j == nc && (  i %% 2L  || !has.upper || !has.lower ))
              localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if(i == j) {
              if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
              if (doText) {
                par(usr = c(0, 1, 0, 1))
                if(is.null(cex.labels)) {
                  l.wid <- strwidth(labels, "user")
                  cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                }
                xlp <- if(xl[i]) 10^0.5 else 0.5
                ylp <- if(yl[j]) 10^label.pos else label.pos
                text.panel(xlp, ylp, labels[i],
                           cex = cex.labels, font = font.labels)
              }
            } else if(i < j){
              if(up.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        seg.col <- 1
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col, ...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex=2*cex)
                  if(show.constraints && is.null(constraints.col)){
                    if(any(linC$start.given)){
                      st.ind <- clusters %in%
                        linC$start.clus[linC$start.given]
                      points(centers[st.ind,j],
                             centers[st.ind,i], cex = cex,
                             col = 'green3',
                             pch = pch)
                    }
                    if(any(linC$end.given)){
                      en.ind <- clusters %in%
                        linC$end.clus[linC$end.given]
                      points(centers[en.ind,j],
                             centers[en.ind,i], cex = cex,
                             col = 'red2', pch = pch)
                    }
                  } else if(show.constraints && !is.null(constraints.col)){
                    for(const in names(constraints.col)) {
                      points(centers[clusters %in% const, j, drop=FALSE],
                             centers[clusters %in% const, i, drop=FALSE],
                             cex = cex, pch = 16,
                             col = constraints.col[const])
                    }
                  }
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)], lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            else{
              if(lp.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        if(clusters[ii] %in%
                           linC$start.clus |
                           clusters[jj] %in%
                           linC$start.clus){
                          seg.col <- 'green3'
                        }else if(clusters[ii] %in%
                                 linC$end.clus[
                                   linC$end.given] |
                                 clusters[jj] %in%
                                 linC$end.clus[
                                   linC$end.given]){
                          seg.col <- 'red2'
                        }else{
                          seg.col <- 1
                        }
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col,...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex = 2*cex)
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)],lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            if (any(par("mfg") != mfg))
              stop("the 'panel' function made a new plot")
          } else par(new = FALSE)
          
        }
      if (!is.null(main)) {
        font.main <- if("font.main" %in% nmdots){
          dots$font.main
        }else par("font.main")
        cex.main <- if("cex.main" %in%
                       nmdots) dots$cex.main else par("cex.main")
        mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main,
              font = font.main)
      }
      invisible(NULL)
    }
  
  #' @name plot3d-SlingshotDataSet
  #' @title Plot Slingshot output in 3D
  #'
  #' @description Tools for visualizing lineages inferred by \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
  #'   Details for more.
  #' @param linInd integer, an index indicating which lineages should be plotted
  #'   (default is to plot all lineages). If \code{col} is a vector, it will be
  #'   subsetted by \code{linInd}.
  #' @param add logical, indicates whether the output should be added to an
  #'   existing plot.
  #' @param dims numeric, which dimensions to plot (default is \code{1:3}).
  #' @param aspect either a logical indicating whether to adjust the aspect ratio
  #'   or a new ratio, see \code{\link[rgl:plot3d]{plot3d}}.
  #' @param size numeric, size of points for MST (default is \code{10}), see
  #'   \code{\link[rgl:plot3d]{plot3d}}.
  #' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
  #' @param ... additional parameters to be passed to \code{lines3d}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' \dontrun{
  #' library(rgl)
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' rd <- cbind(rd, rnorm(nrow(rd)))
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' plot3d(sds, type = 'b')
  #'
  #' # add to existing plot
  #' plot3d(rd, col = 'grey50', aspect = 'iso')
  #' plot3d(sds, lwd = 3, add = TRUE)
  #' }
  # #' @importFrom rgl plot3d
  #' @export
  plot3d.SlingshotDataSet <- function(x,
                                      type = NULL,
                                      linInd = NULL,
                                      add = FALSE,
                                      dims = seq_len(3),
                                      aspect = 'iso',
                                      size = 10,
                                      col = 1,
                                      col2 = NULL,
                                      ...){
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("Package 'rgl' is required for 3D plotting.",
           call. = FALSE)
    }
    col <- rep(col, length(slingLineages(x)))
    n <- nrow(reducedDim(x))
    curves <- FALSE
    lineages <- FALSE
    if(is.null(type)){
      if(length(slingCurves(x)) > 0){
        type <- 'curves'
      }else if(length(slingLineages(x)) > 0){
        type <- 'lineages'
      }else{
        stop('No lineages or curves detected.')
      }
    }else{
      type <- c('curves','lineages','both')[pmatch(type,c('curves','lineages',
                                                          'both'))]
      if(is.na(type)){
        stop('Unrecognized type argument.')
      }
    }
    
    if(type %in% c('lineages','both')){
      lineages <- TRUE
    }
    if(type %in% c('curves','both')){
      curves <- TRUE
    }
    
    if(lineages & (length(slingLineages(x))==0)){
      stop('No lineages detected.')
    }
    if(curves & (length(slingCurves(x))==0)){
      stop('No curves detected.')
    }
    
    if(is.null(linInd)){
      linInd <- seq_along(slingLineages(x))
    }else{
      linInd <- as.integer(linInd)
      if(!all(linInd %in% seq_along(slingLineages(x)))){
        if(any(linInd %in% seq_along(slingLineages(x)))){
          linInd.removed <-
            linInd[! linInd %in% seq_along(slingLineages(x))]
          linInd <-
            linInd[linInd %in% seq_along(slingLineages(x))]
          message('Unrecognized lineage indices (linInd): ',
                  paste(linInd.removed, collapse = ", "))
        }else{
          stop('None of the provided lineage indices',
               ' (linInd) were found.')
        }
      }
    }
    
    if(lineages){
      X <- reducedDim(x)
      clusterLabels <- slingClusterLabels(x)
      connectivity <- slingAdjacency(x)
      clusters <- rownames(connectivity)
      nclus <- nrow(connectivity)
      centers <- t(vapply(clusters,function(clID){
        w <- clusterLabels[,clID]
        return(apply(X, 2, weighted.mean, w = w))
      }, rep(0,ncol(X))))
      rownames(centers) <- clusters
      X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
      clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                     drop = FALSE]
      clus2include <- unique(unlist(slingLineages(x)[linInd]))
    }
    
    if(!add){
      xs <- NULL
      ys <- NULL
      zs <- NULL
      if(lineages){
        xs <- c(xs, centers[,dims[1]])
        ys <- c(ys, centers[,dims[2]])
        zs <- c(zs, centers[,dims[3]])
      }
      if(curves){
        npoints <- nrow(slingCurves(x)[[1]]$s)
        xs <- c(xs, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[1]] }, rep(0,npoints))))
        ys <- c(ys, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[2]] }, rep(0,npoints))))
        zs <- c(zs, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[3]] }, rep(0,npoints))))
      }
      rgl::plot3d(x = NULL, y = NULL, z = NULL, aspect = aspect,
                  xlim = range(xs), ylim = range(ys), zlim = range(zs),
                  xlab = colnames(reducedDim(x))[dims[1]],
                  ylab = colnames(reducedDim(x))[dims[2]],
                  zlab = colnames(reducedDim(x))[dims[3]])
    }
    
    if(lineages){
      for(i in seq_len(nclus-1)){
        for(j in seq(i+1,nclus)){
          if(connectivity[i,j]==1 &
             all(clusters[c(i,j)] %in% clus2include)){
            rgl::lines3d(x = centers[c(i,j),dims[1]],
                         y = centers[c(i,j),dims[2]],
                         z = centers[c(i,j),dims[3]],
                         col = col[1], ...)
          }
        }
      }
      rgl::points3d(centers[clusters %in% clus2include, dims],
                    size = size/2, col = col2[clusters[clusters %in% clus2include]])
      rgl::points3d(centers[clusters %in% clus2include, dims],
                    size = size, col = col[1])
    }
    if(curves){
      for(ii in seq_along(slingCurves(x))[linInd]){
        c <- slingCurves(x)[[ii]]
        rgl::lines3d(c$s[c$ord,dims], col = col[ii], ...)
      }
    }
    invisible(NULL)
  }
  
  ### the plot3d.SlingshotDataSet of the slingshot package is incomplete and too simple,
  ### so, i'm implementing a 3d plot function myself
  slingshot_3d_lineages <- function(slingshot_obj, color, title,
                                    print=FALSE, outputDir=NULL,
                                    width=1200, height=800) {
    
    ### load libraries
    if(!require(Seurat, quietly = TRUE)) {
      install.packages("Seurat")
      require(Seurat, quietly = TRUE)
    }
    if(!require(slingshot, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("slingshot")
      require(slingshot, quietly = TRUE)
    }
    if(!require(rgl, quietly = TRUE)) {
      install.packages("rgl")
      require(rgl, quietly = TRUE)
    }
    if(!require(rmarkdown, quietly = TRUE)) {
      install.packages("rmarkdown")
      require(rmarkdown, quietly = TRUE)
    }
    
    #
    ### 3D Slingshot
    #
    
    ### draw 3D PCA
    par3d(windowRect = c(50, 50, width+50, height+50))
    plot3d.SlingshotDataSet(slingshot_obj, dims = 1:3, col = "black", col2 = color, type = "lineages", add = TRUE)
    plot3d(slingshot_obj@reducedDim, col = apply(slingshot_obj@clusterLabels, 1, function(x) color[names(x)[which(x == 1)]]),
           size = 5, alpha = 0.5, aspect = FALSE, add = TRUE)
    axes3d(edges=c("x+-", "y+-", "z++"), lwd = 2,
           labels=TRUE, tick = FALSE, nticks = 3, box = TRUE, expand = 1.05)
    mtext3d(text = expression(bold("PC1")), edge="x+-", line = -2, at = min(slingshot_obj@reducedDim[,1]), pos = NA)
    mtext3d(text = expression(bold("PC2")), edge="y+-", line = -2, at = min(slingshot_obj@reducedDim[,2]), pos = NA)
    mtext3d(text = expression(bold("PC3")), edge="z++", line = -2, at = max(slingshot_obj@reducedDim[,3]), pos = NA)
    decorate3d(xlim = NULL, ylim = NULL, zlim = NULL, 
               xlab = "", ylab = "", zlab = "", 
               box = FALSE, axes = FALSE, main = title, sub = NULL,
               top = TRUE, aspect = FALSE, expand = 1.05)
    legend3d("topright", legend = names(color), title = "Clusters",
             col = color, pch = 19, cex=2)
    if(print) {
      writeWebGL(dir=outputDir, filename = paste0(outputDir, title, ".html"),
                 width=width, height = height)
    }
    
  }
  
  
  ### get pca
  pca_map <- Embeddings(Seurat_Obj, reduction = "pca")[rownames(Seurat_Obj@meta.data), 1:10]
  
  ### get colors for time points
  cell_colors_clust <- cell_pal(total_time_points, hue_pal())
  
  ### for each patient, run pseudotime analysis
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    
    ### print progress
    writeLines(paste(patient))
    
    ### patient index
    px_idx <- which(Seurat_Obj@meta.data$px == patient)
    
    ### pca with the given patient only
    sub_pca_map <- pca_map[px_idx,]
    
    ### patient specific colors
    cell_colors_px <- cell_colors_clust[which(names(cell_colors_clust) %in% unique(Seurat_Obj@meta.data$time[px_idx]))]
    
    ### get slingshot object
    slingshot_obj <- slingshot(sub_pca_map,
                               clusterLabels = factor(Seurat_Obj@meta.data$time[px_idx],
                                                      levels = names(cell_colors_px)), 
                               reducedDim = "PCA")
    
    ### Trajectory inference
    png(paste0(outputDir, "/", patient, "/Trajectory_Inference_PCA_", patient, ".png"), width = 2500, height = 1500, res = 200)
    plot(sub_pca_map,
         main = paste(patient, "Trajectory Inference Based On Time Points (PCA)"),
         col = cell_colors_clust[as.character(Seurat_Obj@meta.data$time[px_idx])],
         pch = 19, cex = 1.5)
    lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
          show.constraints = TRUE, constraints.col = cell_colors_px)
    legend("bottomleft", legend = names(cell_colors_px), col = cell_colors_px,
           pch = 19)
    dev.off()
    
    ### Trajectory inference on multi dimentional PCA
    png(paste0(outputDir, "/", patient, "/Trajectory_Inference_Multi-PCA_", patient, ".png"), width = 2500, height = 1500, res = 200)
    pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
          show.constraints = TRUE, constraints.col = cell_colors_px, cex = 0.8,
          horInd = 1:5, verInd = 1:5, main = paste(patient, "Trajectory Inference Based On Time Points (PCA)"))
    par(xpd = TRUE)
    legend("bottomleft", legend = names(cell_colors_px), col = cell_colors_px,
           pch = 19, title = "Time")
    dev.off()
    
    # ### 3D Slingshot
    # slingshot_3d_lineages(slingshot_obj = slingshot_obj,
    #                       color = cell_colors_px,
    #                       title = paste("Trajectory_Inference_3D-PCA_", patient),
    #                       print = TRUE,
    #                       outputDir = paste0(outputDir, "/", patient, "/"),
    #                       width = 1200,
    #                       height = 800)
    # rgl.close()
    
    
    gc()
    
  }
  
  #
  ### build a classifier of predicting GMP CAR+ persisters among GMP CAR+ cells
  #
  # load clonotype lineages info
  SJCAR19_Clonotype_Frequency <- readRDS(clonotype_lineage_info_path)
  
  ### GMP CAR+ persistent clones
  pClones <- NULL
  for(i in 1:length(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]])) {
    gmp_idx <- which(colnames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) == "GMP")
    gmp_redo_idx <- which(colnames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) == "GMP-redo")
    last_gmp_idx <- max(gmp_idx, gmp_redo_idx)
    
    ### if at least GMP or GMP-redo exist and there are at least one afterward-time point
    if((last_gmp_idx != -Inf) && (ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) - last_gmp_idx > 1)) {
      ### collect persistent clones that appeared in GMP and persist afterwards
      if(nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) > 0) {
        for(j in 1:nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])) {
          for(k in last_gmp_idx:(ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])-1)) {
            if((SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,"GMP"] > 0 ||
                SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,"GMP-redo"] > 0) &&
               SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,k] > 0) {
              pClones <- c(pClones, rownames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])[j])
              break;
            }
          }
        }
      }
    }
  }
  
  ### GMP CAR+ persistent cells
  pIdx <- intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient %in% pClones),
                    intersect(union(which(Seurat_Obj@meta.data$time == "GMP"),
                                    which(Seurat_Obj@meta.data$time == "GMP-redo")),
                              which(Seurat_Obj@meta.data$CAR == "CARpos")))
  
  ### GMP CAR+ non-persistent cells
  npIdx <- setdiff(intersect(union(which(Seurat_Obj@meta.data$time == "GMP"),
                                   which(Seurat_Obj@meta.data$time == "GMP-redo")),
                             which(Seurat_Obj@meta.data$CAR == "CARpos")),
                   which(Seurat_Obj@meta.data$clonotype_id_by_patient %in% pClones))
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  ### annotate GMP CAR+ persisters
  Seurat_Obj@meta.data$GMP_CARpos_Persister <- NA
  Seurat_Obj@meta.data$GMP_CARpos_Persister[pIdx] <- "YES"
  Seurat_Obj@meta.data$GMP_CARpos_Persister[npIdx] <- "NO"
  
  ### because there are too many cells in non-persisters
  ### randomly select some from those and perform DE analysis
  set.seed(1234)
  Seurat_Obj@meta.data$GMP_CARpos_Persister[sample(npIdx, length(npIdx) - length(pIdx))] <- NA
  
  #'******************************************************************************
  #' A function to transform RNA-Seq data with VST in DESeq2 package
  #' readCount: RNA-Seq rawcounts in a matrix or in a data frame form
  #'            Rows are genes and columns are samples
  #' filter_thresh: The function filters out genes that have at least one sample
  #'                with counts larger than the 'filter_thresh' value
  #'                e.g., if the 'filter_thresh' = 1, then it removes genes
  #'                that have counts <= 1 across all the samples
  #'                if 0, then there will be no filtering
  #'******************************************************************************
  normalizeRNASEQwithVST <- function(readCount, filter_thresh=1) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filter_thresh > 0) {
      ### Remove rubbish rows - this will decrease the number of rows
      keep = apply(counts(deSeqData), 1, function(r){
        return(sum(r > filter_thresh) > 0)
      })
      deSeqData <- deSeqData[keep,]
    }
    
    ### VST
    vsd <- varianceStabilizingTransformation(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  ### a function to select genes based on variance
  selectTopV <- function(x, selectNum) {
    v <- apply(x, 1, var)
    x <- x[order(-v),]
    x <- x[1:selectNum,]
    
    return (x)
  }
  
  ### parameter setting for a classifier
  iteration <- 10
  set.seed(2990)
  featureSelectionNum <- 200
  sampleNum <- 100
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "rf", "LogitBoost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "LogitBoost", "K-NN")
  train_control <- trainControl(method="LOOCV", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  print(identical(rownames(Seurat_Obj@meta.data), colnames(Seurat_Obj@assays$RNA@counts)))
  
  ### iteratively build a classifier
  for(i in 1:iteration) {
    
    ### normalize the read counts
    ### before the normalization, only keep the samples that will be used in the classifier
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[,c(sample(which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES"), sampleNum),
                                                                                                sample(which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO"), sampleNum))],
                                                                stringsAsFactors = FALSE, check.names = FALSE))
    
    ### reduce the gene size based on variance
    ### only select high variance genes
    input_data <- selectTopV(input_data, featureSelectionNum)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(c(rep("GMP_Last", sampleNum),
                                 rep("GMP_Not_Last", sampleNum)),
                               levels = c("GMP_Last", "GMP_Not_Last"))
    
    ### build classifier and test
    ### LOOCV
    p <- list()
    acc <- NULL
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=input_data, trControl=train_control, method=methodTypes[j])
      roc <- roc(model$pred$obs, model$pred$GMP_Last)
      acc <- c(acc, round(mean(model$results$Accuracy), 3))
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using Gene Expressions\n",
                                           "Accuracy =", acc[j]),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      gc()
    }
    
    ### draw ROC curves
    png(paste0(outputDir, "Classifier_GMP_Last_vs_Not_Last_", featureSelectionNum, "_NEW_(", i, ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", acc[j]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
    
    gc()
    
  }
  
}
