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
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
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
  
  ### assign library
  lib_idx <- which(Seurat_Obj@meta.data$barcode_short %in% barcodes[[30]]$V1)
  length(intersect(barcodes[[30]]$V1, Seurat_Obj@meta.data$barcode_short))
  
  
  
  
  
  
  
  
}
