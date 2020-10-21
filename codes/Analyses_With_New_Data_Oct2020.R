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
