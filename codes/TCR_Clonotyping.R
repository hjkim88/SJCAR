###
#   File name : TCR_Clonotyping.R
#   Author    : Hyunjin Kim
#   Date      : Apr 3, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : After combining the GEX and TCR data, the "clonotype_id" represents
#               clonotypes in each time point. Therefore, I would like to do
#               a global clnotypying with all the time points.
#   Instruction
#               1. Source("TCR_Clonotyping.R")
#               2. Run the function "clonotyping" - specify the input file path and the output file path
#               3. The result Robj file will be generated in the output file path
#
#   Example
#               > source("The_directory_of_TCR_Clonotyping.R/TCR_Clonotyping.R")
#               > clonotyping(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_combined.Robj",
#                             outRobjPath="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped.Robj")
###

clonotyping <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_combined.Robj",
                        outRobjPath="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped.Robj") {
  
  ### load library
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### get alpha only and beta only sequences
  cdr3 <- strsplit(Seurat_Obj@meta.data[,"cdr3_nt"], split = ";", fixed = TRUE)
  ### only retain alpha chains
  tcr_a <- sapply(1:length(cdr3), function(x) {
    return(paste(cdr3[[x]][grep("TRA", cdr3[[x]])], collapse = ";"))
  })
  ### only retain beta chains
  tcr_b <- sapply(1:length(cdr3), function(x) {
    return(paste(cdr3[[x]][grep("TRB", cdr3[[x]])], collapse = ";"))
  })
  
  
  ### clonotyping types
  ### 1. global_clonotype_ab_strict0: if the "cdr3_nt" are exactly the same, they are the same clonotype
  ### 2. global_clonotype_a_strict0: if the alpha chains in the "cdr3_nt" are exactly the same, they are the same clonotype
  ### 3. global_clonotype_b_strict0: if the beta chains in the "cdr3_nt" are exactly the same, they are the same clonotype
  ### 4. global_clonotype_ab_lenient0: if only one alpha chain and only one beta chain are the same, they are the same clonotype
  ### the suffix 0 can be changed to 1 or 2, which means the number of different NT sequences that can be accepted
  ### this can affect the chain comparison when clonotyping
  ### e.g., if the suffix = 1, ATGCATGC and ATGCATGG can be the same chain
  ### e.g., if the suffix = 2, ATGCATGC and ATGCATTT can be the same chain
  ### therefore, there would be 4 (described above) x 3 (suffix 0/1/2) = 12 clonotyping types in total
  
  ### a function to determine whether two sequences are in the same clone or not
  is_same_clone <- function(seq1, seq2,
                            option = c("ab_strict", "a_strict", "b_strict", "ab_lenient"),
                            gap = c(0, 1, 2)) {
    
  }
  
  
  ### create columns for the clonotypes
  Seurat_Obj@meta.data$global_clonotype_ab_strict0 <- NA
  Seurat_Obj@meta.data$global_clonotype_a_strict0 <- NA
  Seurat_Obj@meta.data$global_clonotype_b_strict0 <- NA
  Seurat_Obj@meta.data$global_clonotype_ab_lenient0 <- NA
  Seurat_Obj@meta.data$global_clonotype_ab_lenient0[which(!is.na(Seurat_Obj@meta.data$cdr3_nt))] <- ""
  
  ### each library has its own barcode system
  ### so run them separately
  for(lib in unique(Seurat_Obj@meta.data$Library)) {
    ### indicies that assign to the given library
    lib_idx <- which(Seurat_Obj@meta.data$Library == lib)
    
    ### 1 
    ### give the strict version of global clonotypes using all alpha & beta chains
    ### since the variables in each of the "cdr3_nt" are ordered,
    ### if the "cdr3_nt" are exactly the same, they are the same clonotype
    
    ### the unique "cdr3_nt" sequences
    unique_seqs <- unique(Seurat_Obj@meta.data$cdr3_nt[lib_idx])
    
    ### remove NA, "NA", and ""
    unique_seqs <- unique_seqs[-intersect(intersect(which(is.na(unique_seqs)),
                                                    which(unique_seqs == "NA")),
                                          which(unique_seqs == ""))]
    
    ### give clonotypes
    for(i in 1:length(unique_seqs)) {
      idx <- intersect(which(Seurat_Obj@meta.data$cdr3_nt == unique_seqs[i]), lib_idx)
      Seurat_Obj@meta.data$global_clonotype_ab_strict0[idx] <- paste0("clonotype", i)
    }
    
    
    ### 2
    ### give the strict version of global clonotypes using alpha chains only
    
    ### the unique alpha chain sequences
    unique_seqs <- unique(tcr_a[lib_idx])
    
    ### remove NA, "NA", and ""
    unique_seqs <- unique_seqs[-intersect(intersect(which(is.na(unique_seqs)),
                                                    which(unique_seqs == "NA")),
                                          which(unique_seqs == ""))]
    
    ### give clonotypes
    for(i in 1:length(unique_seqs)) {
      idx <- intersect(which(tcr_a == unique_seqs[i]), lib_idx)
      Seurat_Obj@meta.data$global_clonotype_a_strict0[idx] <- paste0("clonotype", i)
    }
    
    
    ### 3
    ### give the strict version of global clonotypes using beta chains only
    
    ### the unique alpha chain sequences
    unique_seqs <- unique(tcr_b[lib_idx])
    
    ### remove NA, "NA", and ""
    unique_seqs <- unique_seqs[-intersect(intersect(which(is.na(unique_seqs)),
                                                    which(unique_seqs == "NA")),
                                          which(unique_seqs == ""))]
    
    ### give clonotypes
    for(i in 1:length(unique_seqs)) {
      idx <- intersect(which(tcr_b == unique_seqs[i]), lib_idx)
      Seurat_Obj@meta.data$global_clonotype_b_strict0[idx] <- paste0("clonotype", i)
    }
    
    
    ### 4
    ### give the lenient version of global clonotypes - at least one alpha and one beta
    ### the unique "cdr3_nt" sequences
    unique_seqs <- unique(Seurat_Obj@meta.data$cdr3_nt[lib_idx])
    
    ### remove NA, "NA", and ""
    unique_seqs <- unique_seqs[-intersect(intersect(which(is.na(unique_seqs)),
                                                    which(unique_seqs == "NA")),
                                          which(unique_seqs == ""))]
    
    ### give clonotypes
    for(i in 1:length(unique_seqs)) {
      idx <- NULL
      for(j in lib_idx) {
        if(is_same_clone(unique_seqs[i], Seurat_Obj@meta.data$cdr3_nt[j],
                         option = "ab_lenient", gap = 0)) {
          idx <- c(idx, j)
        }
      }
      
      Seurat_Obj@meta.data$global_clonotype_ab_lenient0[idx] <- paste(Seurat_Obj@meta.data$global_clonotype_ab_lenient0[idx],
                                                                      paste0("clonotype", i), sep = ";")
    }
    
    ### 8 cases - 4 types x GAP 1 & 2
    
    
  }
  
  
  ### change the object name to the original one
  assign(obj_name, Seurat_Obj)
  rm(Seurat_Obj)
  gc()
  
  ### save the combined Seurat object
  save(list = c(obj_name), file = outRobjPath)
  
}
