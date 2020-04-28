###
#   File name : TCR_Clonotyping.R
#   Author    : Hyunjin Kim
#   Date      : Apr 3, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : After combining the GEX and TCR data, the "clonotype_id" represents
#               clonotypes in each time point. Therefore, I would like to do
#               a global clnotypying with all the time points.
#
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
  if(!require(Biostrings, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Biostrings")
    require(Biostrings, quietly = TRUE)
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
  
  ### if two sequences are given, identify they can be in the same group
  is_same_seq <- function(seq1, seq2, gap = c(0, 1, 2)) {
    ### transform the sequences to DNAStrings
    seq1 <- DNAString(seq1)
    seq2 <- DNAString(seq2)
    
    ### get length of the sequences
    len1 <- nchar(seq1)
    len2 <- nchar(seq2)
    
    ### calculate total gap between the two sequences
    result <- FALSE
    if(abs(len1 - len2) <= gap[1]) {
      min_gap <- gap[1] + 1
      
      for(i in -gap[1]:gap[1]) {
        ### current_gap: initial gap
        ### comparison_bp_len: bps to be compared
        if(i < 0) {
          current_gap <- -i + abs(len2 - len1 - i)
          comparison_bp_len <- min((len1 + i), len2)
          
          ### compare the two matching sequences
          for(j in 1:comparison_bp_len) {
            if(seq1[-i+j] != seq2[j]) {
              current_gap = current_gap + 1
            }
            if(current_gap >= min_gap) {
              break
            }
          }
        } else if(i > 0) {
          current_gap <- i + abs(len1 - len2 + i)
          comparison_bp_len <- min((len2 - i), len1)
          
          ### compare the two matching sequences
          for(j in 1:comparison_bp_len) {
            if(seq1[j] != seq2[i+j]) {
              current_gap = current_gap + 1
            }
            if(current_gap >= min_gap) {
              break
            }
          }
        } else {
          current_gap <- abs(len1 - len2)
          comparison_bp_len <- min(len1, len2)
          
          ### compare the two matching sequences
          for(j in 1:comparison_bp_len) {
            if(seq1[j] != seq2[j]) {
              current_gap = current_gap + 1
            }
            if(current_gap >= min_gap) {
              break
            }
          }
        }
        
        ### update min_gap
        min_gap <- min(min_gap, current_gap)
      }
      
      if(min_gap <= gap[1]) {
        result <- TRUE
      }
    }
    
    return(result)
  }
  
  ### a function to determine whether two sequence sets are in the same clone or not
  is_same_clone <- function(seq1, seq2,
                            option = c("ab_strict", "a_strict", "b_strict", "ab_lenient"),
                            gap = c(0, 1, 2)) {
    
    ### split sequences into chains
    temp1 <- strsplit(seq1, split = ";", fixed = TRUE)[[1]]
    temp2 <- strsplit(seq2, split = ";", fixed = TRUE)[[1]]
    chain_a1 <- temp1[grep("TRA", temp1)]
    chain_b1 <- temp1[grep("TRB", temp1)]
    chain_a2 <- temp2[grep("TRA", temp2)]
    chain_b2 <- temp2[grep("TRB", temp2)]
    chain_a1 <- sapply(chain_a1, function(x) {
      return(strsplit(x, split = ":", fixed = TRUE)[[1]][2])
    })
    chain_b1 <- sapply(chain_b1, function(x) {
      return(strsplit(x, split = ":", fixed = TRUE)[[1]][2])
    })
    chain_a2 <- sapply(chain_a2, function(x) {
      return(strsplit(x, split = ":", fixed = TRUE)[[1]][2])
    })
    chain_b2 <- sapply(chain_b2, function(x) {
      return(strsplit(x, split = ":", fixed = TRUE)[[1]][2])
    })
    
    ### are they in the same clone?
    result <- TRUE
    if(option[1] == "ab_strict" &&
       length(chain_a1) > 0 &&
       length(chain_a2) > 0 &&
       length(chain_b1) > 0 &&
       length(chain_b2) > 0) {
      if((length(chain_a1) == length(chain_a2)) &&
         (length(chain_b1) == length(chain_b2))) {
        all_pairs <- 0
        for(i in 1:length(chain_a1)) {
          for(j in 1:length(chain_a2)) {
            if(is_same_seq(chain_a1[i], chain_a2[j], gap = gap[1])) {
              all_pairs <- all_pairs + 1
              break
            }
          }
        }
        if(all_pairs != length(chain_a1)) {
          result <- FALSE
        } else {
          all_pairs <- 0
          for(i in 1:length(chain_a2)) {
            for(j in 1:length(chain_a1)) {
              if(is_same_seq(chain_a2[i], chain_a1[j], gap = gap[1])) {
                all_pairs <- all_pairs + 1
                break
              }
            }
          }
          if(all_pairs != length(chain_a2)) {
            result <- FALSE
          } else {
            all_pairs <- 0
            for(i in 1:length(chain_b1)) {
              for(j in 1:length(chain_b2)) {
                if(is_same_seq(chain_b1[i], chain_b2[j], gap = gap[1])) {
                  all_pairs <- all_pairs + 1
                  break
                }
              }
            }
            if(all_pairs != length(chain_b1)) {
              result <- FALSE
            } else {
              all_pairs <- 0
              for(i in 1:length(chain_b2)) {
                for(j in 1:length(chain_b1)) {
                  if(is_same_seq(chain_b2[i], chain_b1[j], gap = gap[1])) {
                    all_pairs <- all_pairs + 1
                    break
                  }
                }
              }
              if(all_pairs != length(chain_b1)) {
                result <- FALSE
              }
            }
          }
        }
      } else {
        result <- FALSE
      }
    } else if(option[1] == "a_strict" && length(chain_a1) > 0 && length(chain_a2) > 0) {
      if(length(chain_a1) == length(chain_a2)) {
        all_pairs <- 0
        for(i in 1:length(chain_a1)) {
          for(j in 1:length(chain_a2)) {
            if(is_same_seq(chain_a1[i], chain_a2[j], gap = gap[1])) {
              all_pairs <- all_pairs + 1
              break
            }
          }
        }
        if(all_pairs != length(chain_a1)) {
          result <- FALSE
        } else {
          all_pairs <- 0
          for(i in 1:length(chain_a2)) {
            for(j in 1:length(chain_a1)) {
              if(is_same_seq(chain_a2[i], chain_a1[j], gap = gap[1])) {
                all_pairs <- all_pairs + 1
                break
              }
            }
          }
          if(all_pairs != length(chain_a2)) {
            result <- FALSE
          }
        }
      } else {
        result <- FALSE
      }
    } else if(option[1] == "b_strict" && length(chain_b1) > 0 && length(chain_b2) > 0) {
      if(length(chain_b1) == length(chain_b2)) {
        all_pairs <- 0
        for(i in 1:length(chain_b1)) {
          for(j in 1:length(chain_b2)) {
            if(is_same_seq(chain_b1[i], chain_b2[j], gap = gap[1])) {
              all_pairs <- all_pairs + 1
              break
            }
          }
        }
        if(all_pairs != length(chain_b1)) {
          result <- FALSE
        } else {
          all_pairs <- 0
          for(i in 1:length(chain_b2)) {
            for(j in 1:length(chain_b1)) {
              if(is_same_seq(chain_b2[i], chain_b1[j], gap = gap[1])) {
                all_pairs <- all_pairs + 1
                break
              }
            }
          }
          if(all_pairs != length(chain_b2)) {
            result <- FALSE
          }
        }
      } else {
        result <- FALSE
      }
    } else if(option[1] == "ab_lenient" &&
              length(chain_a1) > 0 &&
              length(chain_a2) > 0 &&
              length(chain_b1) > 0 &&
              length(chain_b2) > 0) {
      any_pairs <- 0
      for(i in 1:length(chain_a1)) {
        for(j in 1:length(chain_a2)) {
          if(is_same_seq(chain_a1[i], chain_a2[j], gap = gap[1])) {
            any_pairs <- 1
            break
          }
        }
        if(any_pairs > 0) {
          break
        }
      }
      if(any_pairs == 0) {
        result <- FALSE
      } else {
        any_pairs <- 0
        for(i in 1:length(chain_b1)) {
          for(j in 1:length(chain_b2)) {
            if(is_same_seq(chain_b1[i], chain_b2[j], gap = gap[1])) {
              any_pairs <- 1
              break
            }
          }
          if(any_pairs > 0) {
            break
          }
        }
        if(any_pairs == 0) {
          result <- FALSE
        }
      }
    } else {
      result <- FALSE
    }
    
    return(result)
    
  }
  
  
  ### create columns for the clonotypes
  Seurat_Obj@meta.data$global_clonotype_ab_strict0 <- NA
  Seurat_Obj@meta.data$global_clonotype_a_strict0 <- NA
  Seurat_Obj@meta.data$global_clonotype_b_strict0 <- NA
  Seurat_Obj@meta.data$global_clonotype_ab_lenient0 <- NA
  Seurat_Obj@meta.data$global_clonotype_ab_lenient0[which(!is.na(Seurat_Obj@meta.data$cdr3_nt))] <- ""
  result_mat <- matrix(NA, nrow(Seurat_Obj@meta.data), 8)
  rownames(result_mat) <- rownames(Seurat_Obj@meta.data)
  colnames(result_mat) <- c("global_clonotype_ab_strict1",
                            "global_clonotype_a_strict1",
                            "global_clonotype_b_strict1",
                            "global_clonotype_ab_lenient1",
                            "global_clonotype_ab_strict2",
                            "global_clonotype_a_strict2",
                            "global_clonotype_b_strict2",
                            "global_clonotype_ab_lenient2")
  result_mat[which(!is.na(Seurat_Obj@meta.data$cdr3_nt)),] <- ""
  
  ### set progress bar
  cnt <- 0
  pb <- txtProgressBar(min = 0, max = length(unique(Seurat_Obj@meta.data$Library)), style = 3)
  
  ### start time
  start_time <- Sys.time()
  
  ### each library has its own barcode system
  ### so run them separately
  for(lib in unique(Seurat_Obj@meta.data$Library)) {
    ### indicies that assign to the given library
    lib_idx <- which(Seurat_Obj@meta.data$Library == lib)
    
    ### remove NA, "NA", and ""
    lib_idx <- intersect(lib_idx, intersect(intersect(which(!is.na(Seurat_Obj@meta.data$cdr3_nt)),
                                                      which(Seurat_Obj@meta.data$cdr3_nt != "NA")),
                                            which(Seurat_Obj@meta.data$cdr3_nt != "")))
    
    ### 1 
    ### give the strict version of global clonotypes using all alpha & beta chains
    ### since the variables in each of the "cdr3_nt" are ordered,
    ### if the "cdr3_nt" are exactly the same, they are the same clonotype
    
    ### the unique "cdr3_nt" sequences
    unique_seqs <- unique(Seurat_Obj@meta.data$cdr3_nt[lib_idx])
    
    ### remove NA, "NA", and ""
    unique_seqs <- unique_seqs[intersect(intersect(which(!is.na(unique_seqs)),
                                                    which(unique_seqs != "NA")),
                                          which(unique_seqs != ""))]
    
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
    unique_seqs <- unique_seqs[intersect(intersect(which(!is.na(unique_seqs)),
                                                   which(unique_seqs != "NA")),
                                         which(unique_seqs != ""))]
    
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
    unique_seqs <- unique_seqs[intersect(intersect(which(!is.na(unique_seqs)),
                                                   which(unique_seqs != "NA")),
                                         which(unique_seqs != ""))]
    
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
    unique_seqs <- unique_seqs[intersect(intersect(which(!is.na(unique_seqs)),
                                                   which(unique_seqs != "NA")),
                                         which(unique_seqs != ""))]
    
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
    for(col in colnames(result_mat)) {
      ### parse the options
      type <- paste(strsplit(col, split = "_", fixed = TRUE)[[1]][3:4], collapse = "_")
      gap <- as.integer(substr(type, nchar(type), nchar(type)))
      type <- substr(type, 1, nchar(type)-1)
      
      ### give clonotypes
      for(i in 1:length(unique_seqs)) {
        idx <- NULL
        for(j in lib_idx) {
          if(is_same_clone(unique_seqs[i], Seurat_Obj@meta.data$cdr3_nt[j],
                           option = type, gap = gap)) {
            idx <- c(idx, j)
          }
        }
        result_mat[idx,col] <- paste(result_mat[idx, col], paste0("clonotype", i), sep = ";")
      }
    }
    
    ### save the partial data
    assign(paste0("vec", cnt), Seurat_Obj@meta.data$global_clonotype_ab_lenient0[lib_idx,])
    assign(paste0("mat", cnt), result_mat[lib_idx,])
    save(list = c("lib_idx", paste0("vec", cnt), paste0("mat", cnt)),
         file = paste0(dirname(outRobjPath), "/partial_result", cnt, ".rda"))
    
    ### update the progress bar
    cnt <- cnt + 1
    setTxtProgressBar(pb, cnt)
  }
  
  ### close the connections
  close(pb)
  
  ### merge result_mat and the Seurat object
  Seurat_Obj@meta.data <- cbind(Seurat_Obj@meta.data, result_mat)
  
  ### end time
  end_time <- Sys.time()
  
  ### print out the running time
  cat(paste("Running Time:",
            signif(as.numeric(difftime(end_time, start_time, units = "mins")), digits = 3),
            "mins"))
  
  ### change the object name to the original one
  assign(obj_name, Seurat_Obj)
  rm(Seurat_Obj)
  gc()
  
  ### save the combined Seurat object
  save(list = c(obj_name), file = outRobjPath)
  
}
