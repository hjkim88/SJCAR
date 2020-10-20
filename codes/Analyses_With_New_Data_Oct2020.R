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
#               > analyses_with_new_data(Seurat_RObj_path="./data/SJCAR19_Oct2020_Seurat_Obj2.RDS",
#                                        barcode_dir="C:/Users/hkim8/SJ/SJCAR19/GEXbarcodes_15Oct2020/",
#                                        TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRs_15Oct2020/",
#                                        outputDir="./results/New/")
###

analyses_with_new_data <- function(Seurat_RObj_path="./data/SJCAR19_Oct2020_Seurat_Obj2.RDS",
                                   barcode_dir="C:/Users/hkim8/SJ/SJCAR19/GEXbarcodes_15Oct2020/",
                                   TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRs_15Oct2020/",
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
  
  
  
  
}
