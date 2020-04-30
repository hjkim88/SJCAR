###
#   File name : TCR_Clonotyping_HPC.R
#   Author    : Hyunjin Kim
#   Date      : Apr 27, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : This is similar to TCR_Clonotyping.R but for the HPC cluster run.
#               The function takes additional options like "lib", "option", and "gap"
#               so that it can be separately run with the specific options.
#
#   Instruction
#               1. Source("TCR_Clonotyping_HPC.R")
#               2. Run the function "clonotyping_hpc" - specify the input file path and the output directory
#               3. The result RDA file will be generated under the output directory
#
#   Example
#               > source("The_directory_of_TCR_Clonotyping_HPC.R/TCR_Clonotyping_HPC.R")
#               > clonotyping_hpc(metadata_path="./data/metadata_hpc.rda",
#                                 lib="JCC212_SJCAR19-05_Wk8_PB",
#                                 option="ab_lenient",
#                                 gap=0,
#                                 outputDir="./data/")
###

clonotyping_hpc <- function(metadata_path="./data/metadata_hpc.rda",
                            lib="JCC212_SJCAR19-05_Wk8_PB",
                            option="ab_lenient",
                            gap=0,
                            outputDir="./data/") {
  
  ### load library
  if(!require(Biostrings, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Biostrings")
    require(Biostrings, quietly = TRUE)
  }
  
  ### load the data
  load(metadata_path)
  
  ### make new output directory
  outputDir <- paste0(outputDir, "/", lib, "/")
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ### get alpha only and beta only sequences
  cdr3 <- strsplit(metadata[,"cdr3_nt"], split = ";", fixed = TRUE)
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
    
    ### chain existence
    is_chain_a <- length(chain_a1) > 0 && length(chain_a2) > 0
    is_chain_b <- length(chain_b1) > 0 && length(chain_b2) > 0
    is_at_least_one <- is_chain_a || is_chain_b
    is_both <- is_chain_a && is_chain_b
    
    ### are they in the same clone?
    result <- TRUE
    if(option[1] == "ab_strict" && is_at_least_one) {
      if((length(chain_a1) == length(chain_a2)) &&
         (length(chain_b1) == length(chain_b2))) {
        all_pairs <- 0
        if(is_chain_a) {
          for(i in 1:length(chain_a1)) {
            for(j in 1:length(chain_a2)) {
              if(is_same_seq(chain_a1[i], chain_a2[j], gap = gap[1])) {
                all_pairs <- all_pairs + 1
                break
              }
            }
          }
        }
        if(all_pairs != length(chain_a1)) {
          result <- FALSE
        } else {
          all_pairs <- 0
          if(is_chain_a) {
            for(i in 1:length(chain_a2)) {
              for(j in 1:length(chain_a1)) {
                if(is_same_seq(chain_a2[i], chain_a1[j], gap = gap[1])) {
                  all_pairs <- all_pairs + 1
                  break
                }
              }
            }
          }
          if(all_pairs != length(chain_a2)) {
            result <- FALSE
          } else {
            all_pairs <- 0
            if(is_chain_b) {
              for(i in 1:length(chain_b1)) {
                for(j in 1:length(chain_b2)) {
                  if(is_same_seq(chain_b1[i], chain_b2[j], gap = gap[1])) {
                    all_pairs <- all_pairs + 1
                    break
                  }
                }
              }
            }
            if(all_pairs != length(chain_b1)) {
              result <- FALSE
            } else {
              all_pairs <- 0
              if(is_chain_b) {
                for(i in 1:length(chain_b2)) {
                  for(j in 1:length(chain_b1)) {
                    if(is_same_seq(chain_b2[i], chain_b1[j], gap = gap[1])) {
                      all_pairs <- all_pairs + 1
                      break
                    }
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
    } else if(option[1] == "a_strict" && is_chain_a) {
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
    } else if(option[1] == "b_strict" && is_chain_b) {
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
    } else if(option[1] == "ab_lenient" && is_at_least_one) {
      if(is_both) {
        any_pairs <- 0
        if(is_chain_a) {
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
        }
        if(any_pairs == 0) {
          result <- FALSE
        } else {
          any_pairs <- 0
          if(is_chain_b) {
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
          }
          if(any_pairs == 0) {
            result <- FALSE
          }
        }
      } else {
        if(is_chain_a) {
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
          }
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
      }
    } else {
      result <- FALSE
    }
    
    return(result)
    
  }
  
  ### create columns for clonotype
  metadata$clonotype <- NA
  metadata$clonotype[which(!is.na(metadata$cdr3_nt))] <- ""
  
  ### indicies that assign to the given library
  lib_idx <- which(metadata$Library == lib)
  
  ### remove NA, "NA", and ""
  lib_idx <- intersect(lib_idx, intersect(intersect(which(!is.na(metadata$cdr3_nt)),
                                                    which(metadata$cdr3_nt != "NA")),
                                          which(metadata$cdr3_nt != "")))
  
  ### get patient id
  patient <- unique(metadata$Px[lib_idx])
  
  ### get patient indices
  px_idx <- which(metadata$Px == patient)
  
  ### remove NA, "NA", and ""
  px_idx <- intersect(px_idx, intersect(intersect(which(!is.na(metadata$cdr3_nt)),
                                                    which(metadata$cdr3_nt != "NA")),
                                          which(metadata$cdr3_nt != "")))
  
  ### the unique "cdr3_nt" sequences
  if(option == "a_strict") {
    unique_seqs <- unique(tcr_a[px_idx])
  } else if(option == "b_strict") {
    unique_seqs <- unique(tcr_b[px_idx])
  } else {
    unique_seqs <- unique(metadata$cdr3_nt[px_idx])
  }
  
  ### remove NA, "NA", and ""
  unique_seqs <- unique_seqs[intersect(intersect(which(!is.na(unique_seqs)),
                                                 which(unique_seqs != "NA")),
                                       which(unique_seqs != ""))]
  
  ### partial file name
  pfName <- paste("partial", lib, option, gap, sep = "_")
  
  ### partial directory
  partialDir <- paste0(outputDir, paste("partial", lib, option, gap, sep = "_"), "/")
  dir.create(partialDir, showWarnings = FALSE, recursive = TRUE)
  
  ### start time
  start_time <- Sys.time()
  
  ### give clonotypes
  if(length(unique_seqs) > 0) {
    for(i in 1:length(unique_seqs)) {
      idx <- NULL
      for(j in lib_idx) {
        if(is_same_clone(unique_seqs[i], metadata$cdr3_nt[j],
                         option = option, gap = gap)) {
          idx <- c(idx, j)
        }
      }
      metadata$clonotype[idx] <- paste(metadata$clonotype[idx], paste0("clonotype", i), sep = ";")
      
      ### save partial results
      partialDir <- paste0(outputDir, paste("partial", lib, option, gap, sep = "_"), "/")
      dir.create(partialDir, showWarnings = FALSE, recursive = TRUE)
      write.table(t(c(paste0("clonotype", i), idx)), file = paste0(partialDir, pfName, ".txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    }
  }
  
  ### end time
  end_time <- Sys.time()
  
  ### print out the running time
  cat(paste("Running Time:",
            signif(as.numeric(difftime(end_time, start_time, units = "mins")), digits = 3),
            "mins"))
  
  ### result file name
  fName <- paste(lib, option, gap, sep = "_")
  
  ### only the library data
  metadata <- metadata[which(metadata$Library == lib),]
  
  ### remove the very first ";" in the clonotype column
  target_idx <- intersect(which(!is.na(metadata$clonotype)),
                          which(metadata$clonotype != ""))
  metadata$clonotype[target_idx] <- substring(metadata$clonotype[target_idx], 2)
  
  ### change the column name
  colnames(metadata)[4] <- paste0("global_clonotype_", option, gap)
    
  ### save the result
  save(list = c("metadata"), file = paste0(outputDir, fName, ".rda"))
  
}
