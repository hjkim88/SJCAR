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
#                                    geneRIFPath="./data/generifs_basic.txt",
#                                    outputDir="./results/PROTO/")
###

clonotype_analysis <- function(Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                               geneRIFPath="./data/generifs_basic.txt",
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
  if(!require(DESeq2, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("DESeq2")
    require(DESeq2, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(metaseqR, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("metaseqR")
    require(metaseqR, quietly = TRUE)
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
  if(!require(OmicCircos, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("OmicCircos")
    require(OmicCircos, quietly = TRUE)
  }
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  
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
  ### CAR+ vs CAR- for every library
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
          target_temp <- target_file[which(target_file[,time_points[i]] > 0),
                                     time_points[(i+1):length(time_points)],drop=FALSE]
          target_clonotypes <- rownames(target_temp)[which(apply(target_temp, 1, sum) > 0)]
          
          if(length(target_clonotypes) > 0) {
            ### set ident.1 and ident.2 for DE analysis
            ### ident.1 - current time point cells that appears later
            ### ident.2 - current time point cells that never appears later
            new.ident <- rep(NA, nrow(Seurat_Obj@meta.data))
            px_time_idx <- intersect(intersect(which(Seurat_Obj@meta.data$Px == px),
                                               which(Seurat_Obj@meta.data$Time == time_points[i])),
                                     which(Seurat_Obj@meta.data$CAR == "CARpos"))
            new.ident[intersect(px_time_idx,
                                which(Seurat_Obj@meta.data[,type] %in% target_clonotypes))] <- "ident1"
            new.ident[setdiff(px_time_idx,
                              which(Seurat_Obj@meta.data[,type] %in% target_clonotypes))] <- "ident2"
            
            ### there should be at least 3 samples in each class for DE analysis
            if((length(which(new.ident == "ident1")) >= 3) && (length(which(new.ident == "ident2")) >= 3)) {
              Idents(object = Seurat_Obj) <- new.ident
              
              ### perform DE analysis
              markers[[px]][[type]][[i]] <- FindMarkers(Seurat_Obj,
                                                        ident.1 = "ident1",
                                                        ident.2 = "ident2",
                                                        logfc.threshold = 0,
                                                        test.use = "DESeq2")
              
              ### garbage collection
              gc()
            }
          }
        }
      }
    }
  }
  
  ### save the marker discovery results
  save(list = c("markers"), file = paste0(outputDir2, "markers.RDATA"))
  
  ###
  ### find common genes among time points and also among patients
  ###
  
  ### set p-value threshold for determining DE genes
  pv_threshold <- 0.05
  
  ### get existing time points
  exs_time_points <- setdiff(levels(Seurat_Obj@meta.data$TimeF), c("PreTrans", "Wk-1", "Wk0"))
  
  ### for each time point
  de_genes <- vector("list", length = length(global_clonotypes))
  names(de_genes) <- global_clonotypes
  for(type in global_clonotypes) {
    de_genes[[type]] <- vector("list", length = length(exs_time_points))
    names(de_genes[[type]]) <- exs_time_points
    for(tp in exs_time_points) {
      de_genes[[type]][[tp]] <- data.frame()
      
      ### combine all the rows
      for(px in names(markers)) {
        if(!is.null(markers[[px]][[type]][[tp]])) {
          thresh_idx <- which(markers[[px]][[type]][[tp]][,"p_val_adj"] < pv_threshold)
          if(length(thresh_idx) > 0) {
            de_genes[[type]][[tp]] <- rbind(de_genes[[type]][[tp]],
                                            cbind(paste0(px, "-", tp),
                                                  rownames(markers[[px]][[type]][[tp]][thresh_idx,]),
                                                  markers[[px]][[type]][[tp]][thresh_idx,]))
          }
        }
      }
      if(length(de_genes[[type]][[tp]]) > 0) {
        colnames(de_genes[[type]][[tp]]) <- c("Px-Time", "Gene_Symbol", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
        
        ### unfactorize the columns
        de_genes[[type]][[tp]][,"Px-Time"] <- as.character(de_genes[[type]][[tp]][,"Px-Time"])
        de_genes[[type]][[tp]][,"Gene_Symbol"] <- as.character(de_genes[[type]][[tp]][,"Gene_Symbol"])
        
        ### add gene name column
        de_genes[[type]][[tp]]$Gene_Name <- mapIds(org.Hs.eg.db, keys=rownames(de_genes[[type]][[tp]]),
                                                   column="GENENAME", keytype="SYMBOL")
        
        ### add count column
        de_genes[[type]][[tp]]$Count <- 1
        
        ### handle duplicated genes
        dup_idx <- which(duplicated(de_genes[[type]][[tp]][,"Gene_Symbol"]))
        if(length(dup_idx) > 0) {
          dups <- de_genes[[type]][[tp]][dup_idx,,drop=FALSE]
          de_genes[[type]][[tp]] <- de_genes[[type]][[tp]][-dup_idx,,drop=FALSE]
          for(i in 1:nrow(dups)) {
            match_idx <- which(de_genes[[type]][[tp]][,"Gene_Symbol"] == dups[i,"Gene_Symbol"])
            
            ### combine multiple rows into one
            combine_cols <- c("Px-Time", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
            de_genes[[type]][[tp]][match_idx,combine_cols] <- paste(de_genes[[type]][[tp]][match_idx,combine_cols],
                                                        dups[i,combine_cols],
                                                        sep = ";")
            de_genes[[type]][[tp]][match_idx,"Count"] <- de_genes[[type]][[tp]][match_idx,"Count"] + 1
          }
        }
        
        ### add combined p_val
        de_genes[[type]][[tp]]$Combined_p_val <- de_genes[[type]][[tp]][,"p_val"]
        multi_idx <- which(de_genes[[type]][[tp]][,"Count"] > 1)
        de_genes[[type]][[tp]][multi_idx,"Combined_p_val"] <- sapply(de_genes[[type]][[tp]][multi_idx,"Combined_p_val"], function(x) {
          temp <- as.numeric(strsplit(x, split = ";", fixed = TRUE)[[1]])
          return(fisher.method(matrix(temp, nrow = 1))[,"p.value"])
        })
        
        ### add combined adj.p_val
        de_genes[[type]][[tp]]$Combined_p_val_adj <- p.adjust(de_genes[[type]][[tp]]$Combined_p_val, method = "BH")
        
        ### order the data frame by 1. combined adj.p and 2. Count 
        de_genes[[type]][[tp]] <- de_genes[[type]][[tp]][order(de_genes[[type]][[tp]][,"Combined_p_val_adj"],
                                                               -de_genes[[type]][[tp]][,"Count"]),]
      }
    }
  }
  
  ### remove NULL results
  for(name1 in names(de_genes)) {
    for(name2 in names(de_genes[[name1]])) {
      if(length(de_genes[[name1]][[name2]]) == 0) {
        de_genes[[name1]][[name2]] <- NULL
      }
    }
  }
  
  ### save the combined results
  save(list = c("de_genes"), file = paste0(outputDir2, "combined_markers.RDATA"))
  
  ### save the combined results in Excel files
  for(type in names(de_genes)) {
    for(tp in names(de_genes[[type]])) {
      write.xlsx2(de_genes[[type]][[tp]], file = paste0(outputDir2, "combined_de_genes_", type, ".xlsx"),
                  row.names = FALSE, sheetName = tp, append = TRUE)
    }
  }
  
  
  ### pathway analysis
  
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
                ggtitle(paste0("KEGG ", title))
              
              png(paste0(dir, "kegg_", title, "_CB.png"), width = 2000, height = 1000)
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
                ggtitle(paste0("GO ", title))
              
              png(paste0(dir, "go_", title, "_CB.png"), width = 2000, height = 1000)
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
  
  ### create empty pathway result list
  pathway_results_GO <- vector("list", length = length(global_clonotypes))
  names(pathway_results_GO) <- global_clonotypes
  pathway_results_KEGG <- vector("list", length = length(global_clonotypes))
  names(pathway_results_KEGG) <- global_clonotypes
  
  ### pathway analysis
  for(type in global_clonotypes) {
    pathway_results_GO[[type]] <- vector("list", length = length(de_genes[[type]]))
    names(pathway_results_GO[[type]]) <- names(de_genes[[type]])
    pathway_results_KEGG[[type]] <- vector("list", length = length(de_genes[[type]]))
    names(pathway_results_KEGG[[type]]) <- names(de_genes[[type]])
    
    for(tp in names(pathway_results_GO[[type]])) {
      pathway_results_GO[[type]][[tp]] <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, de_genes[[type]][[tp]][,"Gene_Symbol"], "ENTREZID", "SYMBOL"),
                                                             org = "human", database = "GO",
                                                             title = paste0("Pathway_Results_", type, "_", tp),
                                                             displayNum = 50, imgPrint = TRUE,
                                                             dir = paste0(outputDir2))
      pathway_results_KEGG[[type]][[tp]] <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, de_genes[[type]][[tp]][,"Gene_Symbol"], "ENTREZID", "SYMBOL"),
                                                               org = "human", database = "KEGG",
                                                               title = paste0("Pathway_Results_", type, "_", tp),
                                                               displayNum = 50, imgPrint = TRUE,
                                                               dir = paste0(outputDir2))
    }
  }
  ### save the pathway results in Excel files
  for(type in global_clonotypes) {
    for(tp in names(pathway_results_GO[[type]])) {
      write.xlsx2(pathway_results_GO[[type]][[tp]], file = paste0(outputDir2, "GO_pathway_results_", type, ".xlsx"),
                  row.names = FALSE, sheetName = tp, append = TRUE)
    }
    for(tp in names(pathway_results_KEGG[[type]])) {
      write.xlsx2(pathway_results_KEGG[[type]][[tp]], file = paste0(outputDir2, "KEGG_pathway_results_", type, ".xlsx"),
                  row.names = FALSE, sheetName = tp, append = TRUE)
    }
  }
  
  ### merge the DE results in each patient
  merged_de_genes <- vector("list", length = length(global_clonotypes))
  names(merged_de_genes) <- global_clonotypes
  for(type in global_clonotypes) {
    
    merged_de_genes[[type]] <- data.frame()
    for(tp in names(de_genes[[type]])) {
      merged_de_genes[[type]] <- rbind(merged_de_genes[[type]],
                                       de_genes[[type]][[tp]])  
    }
    
    write.xlsx2(merged_de_genes[[type]], file = paste0(outputDir2, "Merged_DE_Genes.xlsx"),
                sheetName = type, append = TRUE, row.names = FALSE)
  }
  
  
  #
  ### GeneRIF
  #
  
  ### load geneRIF
  geneRIF <- read.table(file = geneRIFPath, header = FALSE, sep = "\t", check.names = FALSE)
  
  ### get merged & unique DE genes
  input_de_genes <- lapply(merged_de_genes, function(x) {
    y <- unique(x[,"Gene_Symbol"])
    return(mapIds(org.Hs.eg.db, y, "ENTREZID", "SYMBOL"))
  })
  
  ### create an empty list for the results
  geneRIF_results <- vector("list", length = length(global_clonotypes))
  names(geneRIF_results) <- global_clonotypes
  
  ### GeneRIF annotation
  for(type in global_clonotypes) {
    ### extract geneRIF for the specific genes
    add_info <- geneRIF[which(geneRIF$V2 %in% input_de_genes[[type]]),c(2,3,5)]
    colnames(add_info) <- c("Gene_Name", "PubMed_ID", "Gene_RIF")
    
    ### change Entrez IDs to gene symbols
    add_info[,1] <- mapIds(org.Hs.eg.db, as.character(add_info[,1]), "SYMBOL", "ENTREZID")
    
    ### remove duplicates
    rIdx <- duplicated.data.frame(add_info)
    add_info <- add_info[!rIdx,]
    
    ### write out the result
    write.xlsx2(add_info, file = paste0(outputDir2, "Merged_DE_Genes_GeneRIF_Annotated.xlsx"),
                sheetName = type, append = TRUE, row.names = FALSE)
  }
  
  ### just use all the DE genes for pathway analysis
  merged_pathway_results_GO <- vector("list", length = length(global_clonotypes))
  names(merged_pathway_results_GO) <- global_clonotypes
  merged_pathway_results_KEGG <- vector("list", length = length(global_clonotypes))
  names(merged_pathway_results_KEGG) <- global_clonotypes
  for(type in names(merged_de_genes)) {
    target_genes <- unique(merged_de_genes[[type]][,"Gene_Symbol"])
    
    merged_pathway_results_GO[[type]] <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, target_genes, "ENTREZID", "SYMBOL"),
                                                            org = "human", database = "GO",
                                                            title = paste0("Merged_Pathway_Results_", type),
                                                            displayNum = 50, imgPrint = TRUE,
                                                            dir = paste0(outputDir2))
    merged_pathway_results_KEGG[[type]] <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, target_genes, "ENTREZID", "SYMBOL"),
                                                              org = "human", database = "KEGG",
                                                              title = paste0("Merged_Pathway_Results_", type),
                                                              displayNum = 50, imgPrint = TRUE,
                                                              dir = paste0(outputDir2))
  }
  ### save the merged pathway results in Excel files
  for(type in global_clonotypes) {
    write.xlsx2(merged_pathway_results_GO[[type]], file = paste0(outputDir2, "GO_merged_pathway_results.xlsx"),
                row.names = FALSE, sheetName = type, append = TRUE)
    write.xlsx2(merged_pathway_results_KEGG[[type]], file = paste0(outputDir2, "KEGG_merged_pathway_results.xlsx"),
                row.names = FALSE, sheetName = type, append = TRUE)
  }
  
  ### UMAP plot with gene expression for some interesting genes
  
  
  
  
  #
  ### merged lineage abstract table
  #
  for(px in f) {
    ### load the file
    target_file <- read.xlsx2(file = paste0(outputDir, px, "/lineage_abstract_table_", px, ".xlsx"),
                              sheetName = px, row.names = 1,
                              stringsAsFactors = FALSE, check.names = FALSE)
    
    ### numerize the matrix
    rowN <- rownames(target_file)
    target_file <- sapply(target_file, as.numeric)
    rownames(target_file) <- rowN
    
    ### merge
    if(px == f[1]) {
      result_file <- target_file
    } else {
      for(i in 1:nrow(result_file)) {
        for(j in 1:ncol(result_file)) {
          result_file[i,j] <- result_file[i,j] + target_file[i,j]
        }
      }
    }
    
    ### progress
    writeLines(paste(px))
    
    ### garbage collection
    gc()
  }
  
  ### save the result table
  write.xlsx2(result_file, file = paste0(outputDir2, "Merged_lineage_abstract_table.xlsx"),
              sheetName = "Merged_Lineage_Abstract")
  
  
  #
  ### DE analysis - GMP_Last_All vs Control
  ### the above ones are done separately for each patient
  ### here, only consider ab_strict0 and GMP CAR+ cells
  #
  
  ### for each patient that has TCR info
  type <- global_clonotypes[1]
  all_gmp_last <- NULL
  for(px in f) {
    ### print progress
    writeLines(paste(px))
    
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
    
    if(length(time_points) > 1 && length(which(time_points == "GMP")) > 0) {
      ### select clonotypes
      tidx <- which(time_points == "GMP")
      target_temp <- target_file[which(target_file[,"GMP"] > 0),
                                 time_points[(tidx+1):length(time_points)],drop=FALSE]
      target_clonotypes <- rownames(target_temp)[which(apply(target_temp, 1, sum) > 0)]
      
      if(length(target_clonotypes) > 0) {
        px_time_idx <- intersect(intersect(which(Seurat_Obj@meta.data$Px == px),
                                           which(Seurat_Obj@meta.data$Time == time_points[tidx])),
                                 which(Seurat_Obj@meta.data$CAR == "CARpos"))
        all_gmp_last <- c(all_gmp_last, intersect(px_time_idx,
                                                  which(Seurat_Obj@meta.data[,type] %in% target_clonotypes)))
      }
    }
  }
  
  ### Ident configure
  new.ident <- rep(NA, nrow(Seurat_Obj@meta.data))
  new.ident[all_gmp_last] <- "ident1"
  new.ident[setdiff(setdiff(intersect(intersect(which(Seurat_Obj@meta.data$Time == "GMP"),
                                        which(Seurat_Obj@meta.data$CAR == "CARpos")),
                              which(!is.na(Seurat_Obj@meta.data$cdr3_nt))),
                            all_gmp_last),
                    which(Seurat_Obj@meta.data$Px %in% c("SJCAR19-00", "SJCAR19-01")))] <- "ident2"
  Idents(object = Seurat_Obj) <- new.ident
  
  ### DE analysis with Wilcoxon test
  all_de_result <- FindMarkers(Seurat_Obj,
                               ident.1 = "ident1",
                               ident.2 = "ident2",
                               logfc.threshold = 0,
                               min.pct = 0.1)
  
  ### rearange the columns
  all_de_result <- data.frame(Gene_Symbol=rownames(all_de_result),
                              all_de_result[,-which(colnames(all_de_result) == "p_val_adj")],
                              FDR=p.adjust(all_de_result$p_val, method = "BH"),
                              stringsAsFactors = FALSE, check.names = FALSE)
  
  ### save in Excel
  write.xlsx2(all_de_result, file = paste0(outputDir2, "ALL_Patients_DE_Genes.xlsx"),
              sheetName = "ALL", row.names = FALSE)
  
  ### pathway analysis
  target_genes <- all_de_result$Gene_Symbol[which(all_de_result$FDR < 0.05)]
  all_pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, target_genes, "ENTREZID", "SYMBOL"),
                                              org = "human", database = "GO",
                                              title = paste0("All_Pathway_Results"),
                                              displayNum = 50, imgPrint = TRUE,
                                              dir = paste0(outputDir2))
  write.xlsx2(all_pathway_result_GO, file = paste0(outputDir2, "GO_all_pathway_results.xlsx"),
              row.names = FALSE, sheetName = "GO")
  all_pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, target_genes, "ENTREZID", "SYMBOL"),
                                                org = "human", database = "KEGG",
                                                title = paste0("All_Pathway_Results"),
                                                displayNum = 50, imgPrint = TRUE,
                                                dir = paste0(outputDir2))
  write.xlsx2(all_pathway_result_KEGG, file = paste0(outputDir2, "KEGG_all_pathway_results.xlsx"),
              row.names = FALSE, sheetName = "KEGG")
  
  ### DE analysis with DESeq2
  
  ### first becuase we are in short of memory for DESeq2,
  ### we need to extract the exact data from the object and use it only
  
  ### extract the data
  subset_Seurat_Obj <- subset(Seurat_Obj, idents=c("ident1", "ident2"))
  
  ### remove the original Seurat object
  rm(Seurat_Obj)
  gc()
  
  ### perform DESeq2
  all_de_result <- FindMarkers(subset_Seurat_Obj,
                               ident.1 = "ident1",
                               ident.2 = "ident2",
                               logfc.threshold = 0,
                               min.pct = 0.1,
                               test.use = "DESeq2")
  
  ### rearange the columns
  all_de_result <- data.frame(Gene_Symbol=rownames(all_de_result),
                              all_de_result[,-which(colnames(all_de_result) == "p_val_adj")],
                              FDR=p.adjust(all_de_result$p_val, method = "BH"),
                              stringsAsFactors = FALSE, check.names = FALSE)
  
  ### save in Excel
  write.xlsx2(all_de_result, file = paste0(outputDir2, "ALL_Patients_DE_Genes_DESeq2.xlsx"),
              sheetName = "ALL_DESeq2", row.names = FALSE)

  ### pathway analysis
  target_genes <- all_de_result$Gene_Symbol[which(all_de_result$FDR < 0.05)]
  all_pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, target_genes, "ENTREZID", "SYMBOL"),
                                              org = "human", database = "GO",
                                              title = paste0("All_Pathway_Results_DESeq2"),
                                              displayNum = 50, imgPrint = TRUE,
                                              dir = paste0(outputDir2))
  write.xlsx2(all_pathway_result_GO, file = paste0(outputDir2, "GO_all_pathway_results_DESeq2.xlsx"),
              row.names = FALSE, sheetName = "GO")
  all_pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, target_genes, "ENTREZID", "SYMBOL"),
                                                org = "human", database = "KEGG",
                                                title = paste0("All_Pathway_Results_DESeq2"),
                                                displayNum = 50, imgPrint = TRUE,
                                                dir = paste0(outputDir2))
  write.xlsx2(all_pathway_result_KEGG, file = paste0(outputDir2, "KEGG_all_pathway_results_DESeq2.xlsx"),
              row.names = FALSE, sheetName = "KEGG")
  
  
  ### load geneRIF
  geneRIF <- read.table(file = geneRIFPath, header = FALSE, sep = "\t", check.names = FALSE)
  
  ### get merged & unique DE genes
  input_de_genes <- mapIds(org.Hs.eg.db, target_genes, "ENTREZID", "SYMBOL")
  
  ### extract geneRIF for the specific genes
  add_info <- geneRIF[which(geneRIF$V2 %in% input_de_genes),c(2,3,5)]
  colnames(add_info) <- c("Gene_Name", "PubMed_ID", "Gene_RIF")
  
  ### change Entrez IDs to gene symbols
  add_info[,1] <- mapIds(org.Hs.eg.db, as.character(add_info[,1]), "SYMBOL", "ENTREZID")
  
  ### remove duplicates
  rIdx <- duplicated.data.frame(add_info)
  add_info <- add_info[!rIdx,]
  
  ### write out the result
  write.xlsx2(add_info, file = paste0(outputDir2, "DE_Genes_DESeq2_GeneRIF_Annotated.xlsx"),
              sheetName = "DESeq2_GeneRIF", row.names = FALSE)
  
  
  ### Circular plot - visualization of the lineage tracing in each patient
  for(px in f) {
    
    ### print progress
    writeLines(paste(px))
    
    ### load the file
    target_file <- read.xlsx2(file = paste0(outputDir, px, "/car_clonotype_frequency_over_time_", px, ".xlsx"),
                              sheetIndex = 1, stringsAsFactors = FALSE, check.names = FALSE,
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
    
    if(length(time_points) > 1 && length(which(time_points == "GMP")) > 0) {
      ### select clonotypes
      tidx <- which(time_points == "GMP")
      target_temp <- target_file[which(target_file[,"GMP"] > 0),
                                 time_points[(tidx+1):length(time_points)],drop=FALSE]
      target_clonotypes <- rownames(target_temp)[which(apply(target_temp, 1, sum) > 0)]
      
      if(length(target_clonotypes) > 0) {
        target_file <- target_file[target_clonotypes,time_points]
        
      }
    }
    
  }
  
  
}
