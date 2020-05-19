###
#   File name : Clonotype_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : May 5, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Now we have combined GEX and TCR info (GE, clonotypes, and lineages) of the SJCAR19 data
#               I would like to see if there is any clonotype that appeared frequently in the CAR+ cells
#               across PATIENTS, and also want to see gene expression patterns of the CAR+ cells that
#               have a lineage, and their pathways (biological functions).
#
#   Instruction
#               1. Source("Clonotype_Analysis.R")
#               2. Run the function "clonotype_analysis" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Clonotype_Analysis.R/Clonotype_Analysis.R")
#               > clonotype_analysis(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                                    outputDir="./results/PROTO/")
###

clonotype_analysis <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                               outputDir="./results/PROTO/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
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
  
  ### get patient ids (dir names) from the result directory
  f <- list.dirs(outputDir, full.names = FALSE, recursive = FALSE)
  f <- f[grep("SJCAR19", f)]
  
  ### set new result directory
  outputDir2 <- paste0(outputDir, "DEEP/")
  dir.create(path = outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### global clonotypes in the Seurat objects
  global_clonotypes <- colnames(Seurat_Obj@meta.data)[grep("global_clonotype", colnames(Seurat_Obj@meta.data), fixed = TRUE)]
  
  ### clonotypes in CAR+ cells
  clonotypes_car <- vector("list", length = length(f))
  names(clonotypes_car) <- f
  
  ### all nt sequences of interest
  nt_seqs <- vector("list", length = length(global_clonotypes))
  names(nt_seqs) <- global_clonotypes
  
  ### get the clonotypes from the files
  for(i in 1:length(f)) {
    
    ### different types
    clonotypes_car[[i]] <- vector("list", length = 3)
    names(clonotypes_car[[i]]) <- global_clonotypes
    
    ### for each type
    for(j in 1:length(global_clonotypes)) {
      ### load the file
      target_file <- read.xlsx2(file = paste0(outputDir, f[i], "/car_clonotype_frequency_over_time_", f[i], ".xlsx"),
                                sheetName = global_clonotypes[j], stringsAsFactors = FALSE, check.names = FALSE)
      rownames(target_file) <- target_file[,1]
      
      ### numerize the table
      for(k in 2:ncol(target_file)) {
        target_file[,k] <- as.numeric(target_file[,k])
      }
      
      ### save the clonotypes
      clonotypes_car[[i]][[j]] <- vector("list", length = nrow(target_file))
      names(clonotypes_car[[i]][[j]]) <- rownames(target_file)
      
      ### save the cell barcodes
      for(k in 1:nrow(target_file)) {
        ### indicies that have the clonotype in the patient
        target_idx <- intersect(which(Seurat_Obj@meta.data[,global_clonotypes[j]] == rownames(target_file)[k]),
                                which(Seurat_Obj@meta.data$Px == f[i]))
        
        ### save the info
        clonotypes_car[[i]][[j]][[k]] <- Seurat_Obj@meta.data[target_idx,
                                                              c("GexCellFull", "Library", "Px", "Time", "Type", "CAR", "cdr3_nt", global_clonotypes[j])]
        
        ### save nt sequences
        if(grepl("_ab_", global_clonotypes[j])) {
          nt_seqs[[j]] <- c(nt_seqs[[j]], unique(Seurat_Obj@meta.data$cdr3_nt[target_idx]))
        } else if(grepl("_a_", global_clonotypes[j])) {
          nt_seqs[[j]] <- c(nt_seqs[[j]], unique(tcr_a[target_idx]))
        } else if(grepl("_b_", global_clonotypes[j])) {
          nt_seqs[[j]] <- c(nt_seqs[[j]], unique(tcr_b[target_idx]))
        } else {
          stop("ERROR: Check the existence of the global_clonotypes")
        }
      }
      
    }
    
  }
  
  ### get unique nt sequences
  for(i in 1:length(nt_seqs)) {
    nt_seqs[[i]] <- unique(nt_seqs[[i]])
  }
  
  ### for each clonotyping type
  shared_nt_sequences <- vector("list", length = length(nt_seqs))
  names(shared_nt_sequences) <- names(nt_seqs)
  for(i in 1:length(nt_seqs)) {
    ### make an empty matrix for shared interesting nt sequences across patients
    shared_nt_sequences[[i]] <- matrix(0, nrow = length(nt_seqs[[i]]), ncol = length(f))
    rownames(shared_nt_sequences[[i]]) <- nt_seqs[[i]]
    colnames(shared_nt_sequences[[i]]) <- f
    
    ### start time
    start_time <- Sys.time()
    
    ### set progress bar
    pb <- txtProgressBar(min = 0, max = length(nt_seqs[[i]])*length(f), style = 3)
    
    ### get the numbers for the matrix
    cnt <- 0
    for(sq in nt_seqs[[i]]) {
      for(px in f) {
        if(grepl("_ab_", names(nt_seqs)[i])) {
          shared_nt_sequences[[i]][sq,px] <- length(intersect(which(Seurat_Obj@meta.data$Px == px),
                                                         which(Seurat_Obj@meta.data$cdr3_nt == sq)))
        } else if(grepl("_a_", names(nt_seqs)[i])) {
          shared_nt_sequences[[i]][sq,px] <- length(intersect(which(Seurat_Obj@meta.data$Px == px),
                                                         which(tcr_a == sq)))
        } else if(grepl("_b_", names(nt_seqs)[i])) {
          shared_nt_sequences[[i]][sq,px] <- length(intersect(which(Seurat_Obj@meta.data$Px == px),
                                                         which(tcr_b == sq)))
        } else {
          stop("ERROR: Check the existence of the global_clonotypes")
        }
        
        cnt <- cnt + 1
        setTxtProgressBar(pb, cnt)
      }
    }
    close(pb)
    
    ### end time
    end_time <- Sys.time()
    
    ### print out the running time
    cat(paste("Running Time:",
              signif(as.numeric(difftime(end_time, start_time, units = "mins")), digits = 3),
              "mins"))
    
    ### filter the result
    keep_idx <- NULL
    for(j in 1:nrow(shared_nt_sequences[[i]])) {
      if(length(which(shared_nt_sequences[[i]][j,] > 0)) > 1) {
        keep_idx <- c(keep_idx, j)
      }
    }
    shared_nt_sequences[[i]] <- shared_nt_sequences[[i]][keep_idx,,drop=FALSE]
  }
  
  ### save the results in Excel file
  for(i in 1:length(shared_nt_sequences)) {
    ### order the results based on:
    ### 1. the number of different time points
    ### 2. the total number of appearance
    time_points <- apply(shared_nt_sequences[[i]], 1, function(x) {
      return(length(which(x > 0)))
    })
    total_appear <- apply(shared_nt_sequences[[i]], 1, sum)
    shared_nt_sequences[[i]] <- shared_nt_sequences[[i]][order(-time_points, -total_appear),,drop=FALSE]
    
    write.xlsx2(data.frame(shared_nt_sequences[[i]],
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir2, "/", fName, ".xlsx"),
                sheetName = names(shared_nt_sequences)[i],
                append = TRUE)
  }
  
  ### draw a plot with the results
  top_threshold <- 5
  p <- vector("list", length = length(nt_seqs))
  names(p) <- names(nt_seqs)
  for(i in 1:length(shared_nt_sequences)) {
    ### select top 5 nt sequences
    target_idx <- 1:min(nrow(shared_nt_sequences[[i]]), top_threshold)
    wrapped_row_names <- sapply(rownames(shared_nt_sequences[[i]]), function(x) {
      paste(strsplit(x, split = ";", fixed = TRUE)[[1]], collapse = "\n")
    })
    
    ### prepare a dataframe for the plot
    plot_df <- data.frame(Number=as.vector(shared_nt_sequences[[i]][target_idx,,drop=FALSE]),
                          Patient=as.vector(sapply(colnames(shared_nt_sequences[[i]][target_idx,,drop=FALSE]), function(x) rep(x, nrow(shared_nt_sequences[[i]][target_idx,,drop=FALSE])))),
                          NT_AA=rep(wrapped_row_names[target_idx], ncol(shared_nt_sequences[[i]][target_idx,,drop=FALSE])),
                          stringsAsFactors = FALSE, check.names = FALSE)
    plot_df$Number[which(plot_df$Number == 0)] <- NA
    
    ### draw a bar plot with the result
    p[[i]] <- ggplot(plot_df, aes_string(x="Patient", y="Number", fill="NT_AA", group="NT_AA")) +
      labs(x="", y="The Clone Size in the Patients") +
      geom_bar(position = "dodge", stat = "identity") +
      geom_text(aes_string(label="Number", color="NT_AA", group="NT_AA"),
                position=position_dodge(width=1), size=3.5, hjust=0.5, vjust=-0.25,
                show.legend = FALSE) +
      ggtitle(paste0("CAR+ Clonotypes Across Patients - ", names(nt_seqs)[i])) +
      guides(fill=guide_legend(title=NULL)) +
      ylim(0, max(plot_df$Number, na.rm = TRUE)*1.1) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 14),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14),
            axis.title.y = element_text(size = 10))
    
  }
  
  ### arrange the plots and save
  fName <- paste0("Top_", top_threshold, "_CARpos_Clonotypes_Across_Patients")
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   ncol = 1,
                   top = fName)
  ggsave(file = paste0(outputDir2, "/", fName, ".png"), g, width = 20, height = 12, dpi = 300)
  
  
  ###
  ### Gene expression profiling
  ###
  
  ### set new group for DE
  new.ident <- paste0(Seurat_Obj@meta.data$Library, "_", Seurat_Obj@meta.data$CAR)
  Idents(object = Seurat_Obj) <- new.ident
  
  ### start time
  start_time <- Sys.time()
  
  ### set progress bar
  pb <- txtProgressBar(min = 0, max = length(unique(Seurat_Obj@meta.data$Library)), style = 3)
  
  ### in every library, compare CAR+ and CAR-
  compare_car <- list()
  cnt <- 0
  for(lib in unique(Seurat_Obj@meta.data$Library)) {
    ### get indicies for the library
    lib_idx <- which(Seurat_Obj@meta.data$Library == lib)
    
    ### run only if there are both CARpos and CARneg in the given library
    car <- unique(Seurat_Obj@meta.data$CAR[lib_idx])
    if(length(car) == 2) {
      compare_car <- c(compare_car, list(FindMarkers(Seurat_Obj,
                                                     ident.1 = paste0(lib, "_", car[1]),
                                                     ident.2 = paste0(lib, "_", car[2]),
                                                     logfc.threshold = 0)))
      names(compare_car)[length(compare_car)] <- lib
    }
    
    cnt <- cnt + 1
    setTxtProgressBar(pb, cnt)
  }
  close(pb)
  
  ### end time
  end_time <- Sys.time()
  
  ### print out the running time
  cat(paste("Running Time:",
            signif(as.numeric(difftime(end_time, start_time, units = "mins")), digits = 3),
            "mins"))
  
  ### get significant results
  sig_results <- data.frame()
  for(lib in names(compare_car)) {
    significant_idx <- which(compare_car[[lib]][,"p_val_adj"] < 0.05)
    if(length(significant_idx) > 0) {
      sig_results <- rbind(sig_results, data.frame(compare_car[[lib]], lib=lib, stringsAsFactors = FALSE, check.names = FALSE))
    }
  }
  
  ### write the significant result
  write.xlsx2(sig_results,
              file = paste0(outputDir2, "/DE_CAR_pos_vs_neg_0.05.xlsx"))
  
  
  ###
  ### find markers for each clonotype in each patient
  ###
  
  ### set frequency cut-off
  ### we are finding markers for clonotypes that appeared equal to or more than the cut-off in the patient
  # freq_cut_off <- 10
  ### we are finding markers with the top clonotypes in each patient
  top_freq <- 5
  
  ### create an empty maker list
  markers <- vector("list", length = length(global_clonotypes))
  names(markers) <- global_clonotypes
  
  ### for each clonotyping type
  for(type in global_clonotypes) {
    
    ### show progress
    writeLines(paste("###", type))
    
    ### new idents for the analysis
    new.ident <- Seurat_Obj@meta.data[,type]
    
    ### create an empty maker list for the type
    markers[[type]] <- vector("list", length = length(unique(Seurat_Obj@meta.data$Px)))
    names(markers[[type]]) <- unique(Seurat_Obj@meta.data$Px)
    
    ### for each patient update idents
    ### we will find markers for clonotypes that have more than one frequencies
    for(px in unique(Seurat_Obj@meta.data$Px)) {
      ### indicies for the given patient
      px_idx <- which(Seurat_Obj@meta.data$Px == px)
      
      ### load the car clonotype frequency data
      freq_d <- read.xlsx2(file = paste0(outputDir, px, "/car_clonotype_frequency_over_time_", px, ".xlsx"),
                           sheetName = type, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
      
      ### numerize the table
      for(j in 1:ncol(freq_d)) {
        freq_d[,j] <- as.numeric(freq_d[,j])
      }
      
      ### select clonotypes that appeared in GMP and also in any of the later time points
      freq_d <- freq_d[which(freq_d$GMP > 0),
                             setdiff(colnames(freq_d),
                                     c("Clonotype", "PreTrans", "GMP", "Wk-1", "Wk0", "Total"))]
      target_clonotypes <- rownames(freq_d)[which(apply(freq_d, 1, sum) > 0)]
      
      ###
      
      # ### duplicated clonotypes in the patient
      # dups <- unique(Seurat_Obj@meta.data[px_idx[which(duplicated(Seurat_Obj@meta.data[px_idx,type]))],type])
      # 
      # ### remove NA, "", and "NA"
      # dups <- dups[intersect(intersect(which(!is.na(dups)),
      #                                  which(dups != "")),
      #                        which(dups != "NA"))]
      
      # ### remove dups that appeared less than freq_cut_off in the list
      # dup_remove_idx <- NULL
      # for(i in 1:length(dups)) {
      #   freq <- length(intersect(px_idx, which(Seurat_Obj@meta.data[,type] == dups[i])))
      #   if(freq < freq_cut_off) {
      #     dup_remove_idx <- c(dup_remove_idx, i)
      #   }
      # }
      # if(length(dup_remove_idx) > 0) {
      #   dups <- dups[-dup_remove_idx]
      # }
      
      ### create an empty maker list for the type
      markers[[type]][[px]] <- vector("list", length = length(dups))
      names(markers[[type]][[px]]) <- dups
      
      ### get unique clonotype indicies in the patient (to be modified to NA) 
      remove_idx <- setdiff(px_idx, intersect(px_idx, which(Seurat_Obj@meta.data[,type] %in% dups)))
      
      new.ident[remove_idx] <- NA
    }
    
    ### set new idents
    Seurat_Obj@meta.data$Marker_Cluster <- new.ident
    Idents(object = Seurat_Obj) <- Seurat_Obj@meta.data$Px
    
    ### start time
    start_time <- Sys.time()
    
    ### find markers in each patient
    cnt <- 0
    for(px in unique(Seurat_Obj@meta.data$Px)) {
      ### find markers for each clonotype
      for(clono in names(markers[[type]][[px]])) {
        ### find differentially expressed genes
        markers[[type]][[px]][[clono]] <- FindMarkers(Seurat_Obj,
                                                      ident.1 = clono,
                                                      group.by = "Marker_Cluster",
                                                      subset.ident = px)
        
        ### update the progress
        cnt <- cnt + 1
        writeLines(paste(cnt))
      }
    }
    
    ### end time
    end_time <- Sys.time()
    
    ### print out the running time
    cat(paste("Running Time:",
              signif(as.numeric(difftime(end_time, start_time, units = "mins")), digits = 3),
              "mins"))
    
  }
  
  ### save the marker discovery results
  save(list = c("markers"), file = paste0(outputDir2, "markers.RDATA"))
  
  
  
  ###
  ### with the shared NT sequences across all the patients,
  ### we would like to see conserved genes among the cells that have the same NT sequences (same clonotype)
  ###
  
  ### only use the clonotypes that were appeared in different time points more than the cut-off
  time_point_cut_off <- 2
  target_clonotypes <- data.frame()
  for(i in 1:length(shared_nt_sequences)) {
    for(j in 1:nrow(shared_nt_sequences[[i]])) {
      if(length(which(shared_nt_sequences[[i]][j,] > 0)) > time_point_cut_off) {
        target_clonotypes <- rbind(target_clonotypes, shared_nt_sequences[[i]][j,,drop=FALSE])
      }
    }
  }
  
  ### give new idents for the gene conservation analysis
  new.ident <- rep(NA, nrow(Seurat_Obj@meta.data))
  for(nt_seq in rownames(target_clonotypes)) {
    new.ident[which(Seurat_Obj@meta.data$)]
  }
  Idents(object = Seurat_Obj) <- new.ident
  
  ### compute the conserved genes for each of the selected clonotype
  for(nt_seq in rownames(target_clonotypes)) {
    
    
    
  }
  
  
}
