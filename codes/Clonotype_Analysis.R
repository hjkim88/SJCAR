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
  ### GMP cells that also appear later vs not
  ### Wk1 cells that also appear later vs not
  ### ...
  ### 3mo cells that also appear later vs not
  ### to find genes that makes the persistence
  ### what are the GE patterns of the persistence?
  ###
  
  ### an empty list for saving the marker results
  markers <- vector("list", length = length(f))
  names(markers) <- f
  
  ### for each patient that has TCR info
  for(px in f) {
    ### print progress
    writeLines(paste(px))
    
    ### an empty list for saving the marker results
    markers[[px]] <- vector("list", length = length(global_clonotypes))
    names(markers[[px]]) <- global_clonotypes
    
    ### for each type
    for(type in global_clonotypes) {
      ### load the file
      target_file <- read.xlsx2(file = paste0(outputDir, px, "/car_clonotype_frequency_over_time_", px, ".xlsx"),
                                sheetName = type, stringsAsFactors = FALSE, check.names = FALSE,
                                row.names = 1)
      
      ### numerize the table
      for(i in 1:ncol(target_file)) {
        target_file[,i] <- as.numeric(target_file[,i])
      }
      
      ### remove all zero time points
      time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
      
      ### get time points after GMP infusion
      time_points <- setdiff(time_points,
                             c("PreTrans", "Wk-1", "Wk0", "Total"))
      
      ### an empty list for saving the marker results
      markers[[px]][[type]] <- vector("list", length = length(time_points)-1)
      names(markers[[px]][[type]]) <- time_points[-length(time_points)]
      
      if(length(time_points) > 1) {
        ### for each time point - current time point cells that appears later vs not
        for(i in 1:(length(time_points)-1)) {
          ### select clonotypes
          target_file <- target_file[which(target_file[,time_points[i]] > 0),
                                    time_points[(i+1):length(time_points)]]
          target_clonotypes <- rownames(target_file)[which(apply(target_file, 1, sum) > 0)]
          
          if(length(target_clonotypes) > 0) {
            ### set ident.1 and ident.2 for DE analysis
            ### ident.1 - current time point cells that appears later
            ### ident.2 - current time point cells that never appears later
            new.ident <- rep(NA, nrow(Seurat_Obj@meta.data))
            px_time_idx <- intersect(which(Seurat_Obj@meta.data$Px == px),
                                     which(Seurat_Obj@meta.data$Time == time_points[i]))
            new.ident[intersect(px_time_idx,
                                which(Seurat_Obj@meta.data[,type] %in% target_clonotypes))] <- "ident1"
            new.ident[intersect(px_time_idx,
                                which(is.na(Seurat_Obj@meta.data[,type])))] <- "ident2"
            Idents(object = Seurat_Obj) <- new.ident
            
            ### perform DE analysis
            markers[[px]][[type]][[i]] <- FindMarkers(Seurat_Obj,
                                                      ident.1 = "ident1",
                                                      ident.2 = "ident2")
          }
        }
      }
    }
  }
  
  ### save the marker discovery results
  save(list = c("markers"), file = paste0(outputDir2, "markers.RDATA"))
  
  ### find common genes among time points and also among patients
  
  
}
