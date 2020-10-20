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
#               > analyses_with_new_data(Seurat_RObj_path="./data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                        barcode_dir="C:/Users/hkim8/SJ/SJCAR19/GEXbarcodes_15Oct2020/",
#                                        TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRs_15Oct2020/",
#                                        outputDir="./results/New/")
###

analyses_with_new_data <- function(Seurat_RObj_path="./data/SJCAR19_Oct2020_Seurat_Obj.RDS",
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
  
  ### refine library name
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
  Seurat_Obj@meta.data$file_name <- ""
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
      Seurat_Obj@meta.data$file_name[tag_idx] <- max_lib
      Seurat_Obj@meta.data$library[tag_idx] <- lib_name[max_lib]
      Seurat_Obj@meta.data$px[tag_idx] <- px_name[max_lib]
      Seurat_Obj@meta.data$time[tag_idx] <- time_name[max_lib]
      Seurat_Obj@meta.data$tissue[tag_idx] <- tissue_name[max_lib]
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
