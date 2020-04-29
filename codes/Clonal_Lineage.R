###
#   File name : Clonal_Lineage.R
#   Author    : Hyunjin Kim
#   Date      : Apr 12, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : For every clonotype make a table that how proportion of clonotypes change
#               over time - the absolute number and the percentage. Also draw
#               a UMAP for visualization of Clone size among CAR+ cells with TCR info.
#               Track the clonotypes in GMP - CAR+ cells to distinguish which are
#               never going to be seen again and which appeared in ealier time again.
#
#   Instruction
#               1. Source("Clonal_Lineage.R")
#               2. Run the function "lineage_analysis" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Clonal_Lineage.R/Clonal_Lineage.R")
#               > lineage_analysis(Seurat_RObj_path="./data/JCC212_Px5_TCR_clonotyped2.Robj",
#                                  outputDir="./results/")
###

lineage_analysis <- function(Seurat_RObj_path="./data/JCC212_Px5_TCR_clonotyped2.Robj",
                             outputDir="./results/") {
  
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
  
  ### global clonotypes in the Seurat objects
  global_clonotypes <- colnames(Seurat_Obj@meta.data)[grep("global_clonotype", colnames(Seurat_Obj@meta.data), fixed = TRUE)]
  
  for(patient in unique(Seurat_Obj@meta.data$Px)) {
    
    ### indicies that assign to the given patient
    px_idx <- which(Seurat_Obj@meta.data$Px == patient)
    
    ### remove NA, "NA", and ""
    px_idx <- intersect(px_idx, intersect(intersect(which(!is.na(Seurat_Obj@meta.data$cdr3_nt)),
                                                    which(Seurat_Obj@meta.data$cdr3_nt != "NA")),
                                          which(Seurat_Obj@meta.data$cdr3_nt != "")))
    
    ### make a table that shows how clonotype frequency (and proportion) changes over time
    unique_clonotypes <- vector("list", length(global_clonotypes))
    names(unique_clonotypes) <- global_clonotypes
    
    unique_time_points <- vector("list", length(global_clonotypes))
    names(unique_time_points) <- global_clonotypes
    
    frequency_over_time <- vector("list", length(global_clonotypes))
    names(frequency_over_time) <- global_clonotypes
    
    car_frequency_over_time <- vector("list", length(global_clonotypes))
    names(car_frequency_over_time) <- global_clonotypes
    
    proportion_over_time <- vector("list", length(global_clonotypes))
    names(proportion_over_time) <- global_clonotypes
    
    car_proportion_over_time <- vector("list", length(global_clonotypes))
    names(car_proportion_over_time) <- global_clonotypes
    
    for(type in global_clonotypes) {
      target_idx <- intersect(px_idx, which(!is.na(Seurat_Obj@meta.data[,type])))
      
      unique_clonotypes[[type]] <- unique(Seurat_Obj@meta.data[target_idx,type])
      
      unique_time_points[[type]] <- levels(Seurat_Obj@meta.data$TimeF[target_idx])
      
      frequency_over_time[[type]] <- matrix(0, length(unique_clonotypes[[type]]), length(unique_time_points[[type]]))
      rownames(frequency_over_time[[type]]) <- unique_clonotypes[[type]]
      colnames(frequency_over_time[[type]]) <- unique_time_points[[type]]
      
      car_frequency_over_time[[type]] <- matrix(0, length(unique_clonotypes[[type]]), length(unique_time_points[[type]]))
      rownames(car_frequency_over_time[[type]]) <- unique_clonotypes[[type]]
      colnames(car_frequency_over_time[[type]]) <- unique_time_points[[type]]
      
      proportion_over_time[[type]] <- matrix(0, length(unique_clonotypes[[type]]), length(unique_time_points[[type]]))
      rownames(proportion_over_time[[type]]) <- unique_clonotypes[[type]]
      colnames(proportion_over_time[[type]]) <- unique_time_points[[type]]
      
      car_proportion_over_time[[type]] <- matrix(0, length(unique_clonotypes[[type]]), length(unique_time_points[[type]]))
      rownames(car_proportion_over_time[[type]]) <- unique_clonotypes[[type]]
      colnames(car_proportion_over_time[[type]]) <- unique_time_points[[type]]
      
      ### start time
      start_time <- Sys.time()
      
      ### set progress bar
      pb <- txtProgressBar(min = 0, max = length(unique_clonotypes_strict)*length(unique_time_points), style = 3)
      
      ### fill out the frequency tables
      cnt <- 0
      for(time in unique_time_points[[type]]) {
        for(clonotype in unique_clonotypes[[type]]) {
          frequency_over_time[[type]][clonotype,time] <- length(intersect(which(Seurat_Obj@meta.data[,type] == clonotype),
                                                                  which(Seurat_Obj@meta.data$Time == time)))
          car_frequency_over_time[[type]][clonotype,time] <- length(intersect(intersect(which(Seurat_Obj@meta.data[,type] == clonotype),
                                                                                which(Seurat_Obj@meta.data$Time == time)),
                                                                      which(Seurat_Obj@meta.data$CAR == "CARpos")))
          
          ### update the progress bar
          cnt <- cnt + 1
          setTxtProgressBar(pb, cnt)
        }
        
        ### fill out the proportion tables
        cells_time_num <- length(which(Seurat_Obj@meta.data$Time == time))
        car_cells_time_num <- length(intersect(which(Seurat_Obj@meta.data$Time == time),
                                               which(Seurat_Obj@meta.data$CAR == "CARpos")))
        if(car_cells_time_num > 0) {
          proportion_over_time[[type]][,time] <- signif(100*frequency_over_time[[type]][,time]/cells_time_num, digits = 4)
          car_proportion_over_time[[type]][,time] <- signif(100*car_frequency_over_time[[type]][,time]/car_cells_time_num, digits = 4)
        } else if(cells_time_num > 0) {
          proportion_over_time[[type]][,time] <- signif(100*frequency_over_time[[type]][,time]/cells_time_num, digits = 4)
        } else {
          proportion_over_time[[type]][,time] <- 0
          car_proportion_over_time[[type]][,time] <- 0
        }
      }
      close(pb)
      
      ### end time
      end_time <- Sys.time()
      
      ### print out the running time
      cat(paste("Running Time:",
                signif(as.numeric(difftime(end_time, start_time, units = "mins")), digits = 3),
                "mins"))
      
      
      
    }
    
    
    
    
    
  }
  
  
  
  unique_clonotypes_strict <- unique(Seurat_Obj@meta.data$global_clonotype_strict[which(!is.na(Seurat_Obj@meta.data$global_clonotype_strict))])
  unique_time_points <- levels(Seurat_Obj@meta.data$TimeF)
  frequency_over_time <- matrix(0, length(unique_clonotypes_strict), length(unique_time_points))
  rownames(frequency_over_time) <- unique_clonotypes_strict
  colnames(frequency_over_time) <- unique_time_points
  car_frequency_over_time <- matrix(0, length(unique_clonotypes_strict), length(unique_time_points))
  rownames(car_frequency_over_time) <- unique_clonotypes_strict
  colnames(car_frequency_over_time) <- unique_time_points
  proportion_over_time <- matrix(0, length(unique_clonotypes_strict), length(unique_time_points))
  rownames(proportion_over_time) <- unique_clonotypes_strict
  colnames(proportion_over_time) <- unique_time_points
  car_proportion_over_time <- matrix(0, length(unique_clonotypes_strict), length(unique_time_points))
  rownames(car_proportion_over_time) <- unique_clonotypes_strict
  colnames(car_proportion_over_time) <- unique_time_points
  
  ### start time
  start_time <- Sys.time()
  
  ### set progress bar
  pb <- txtProgressBar(min = 0, max = length(unique_clonotypes_strict)*length(unique_time_points), style = 3)
  
  ### fill out the frequency tables
  cnt <- 0
  for(time in unique_time_points) {
    for(clonotype in unique_clonotypes_strict) {
      frequency_over_time[clonotype,time] <- length(intersect(which(Seurat_Obj@meta.data$global_clonotype_strict == clonotype),
                                                              which(Seurat_Obj@meta.data$Time == time)))
      car_frequency_over_time[clonotype,time] <- length(intersect(intersect(which(Seurat_Obj@meta.data$global_clonotype_strict == clonotype),
                                                                            which(Seurat_Obj@meta.data$Time == time)),
                                                        which(Seurat_Obj@meta.data$CAR == "CARpos")))
      
      ### update the progress bar
      cnt <- cnt + 1
      setTxtProgressBar(pb, cnt)
    }
    
    ### fill out the proportion tables
    cells_time_num <- length(which(Seurat_Obj@meta.data$Time == time))
    car_cells_time_num <- length(intersect(which(Seurat_Obj@meta.data$Time == time),
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
  frequency_over_time <- data.frame(Clonotype=unique_clonotypes_strict,
                                    frequency_over_time, Total=apply(frequency_over_time, 1, sum),
                                    stringsAsFactors = FALSE, check.names = FALSE)
  car_frequency_over_time <- data.frame(Clonotype=unique_clonotypes_strict,
                                        car_frequency_over_time, Total=apply(car_frequency_over_time, 1, sum),
                                        stringsAsFactors = FALSE, check.names = FALSE)
  proportion_over_time <- data.frame(Clonotype=unique_clonotypes_strict,
                                     proportion_over_time, Total=apply(proportion_over_time, 1, sum),
                                     stringsAsFactors = FALSE, check.names = FALSE)
  car_proportion_over_time <- data.frame(Clonotype=unique_clonotypes_strict,
                                         car_proportion_over_time, Total=apply(car_proportion_over_time, 1, sum),
                                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### order the data frames by the Total column
  frequency_over_time <- frequency_over_time[order(-frequency_over_time[,"Total"]),]
  proportion_over_time <- proportion_over_time[order(-proportion_over_time[,"Total"]),]
  car_frequency_over_time <- car_frequency_over_time[order(-car_frequency_over_time[,"Total"]),]
  car_proportion_over_time <- car_proportion_over_time[order(-car_proportion_over_time[,"Total"]),]
  
  ### save the tables in Excel format
  write.xlsx2(frequency_over_time, file = paste0(outputDir, "clonotype_frequency_over_time.xlsx"),
              sheetName = "clonotype_frequency_over_time", row.names = FALSE)
  write.xlsx2(proportion_over_time, file = paste0(outputDir, "clonotype_frequency_over_time.xlsx"),
              sheetName = "proportion_frequency_over_time", row.names = FALSE, append = TRUE)
  write.xlsx2(car_frequency_over_time, file = paste0(outputDir, "car_clonotype_frequency_over_time.xlsx"),
              sheetName = "car_clonotype_frequency_over_time", row.names = FALSE)
  write.xlsx2(car_proportion_over_time, file = paste0(outputDir, "car_clonotype_frequency_over_time.xlsx"),
              sheetName = "car_proportion_frequency_over_time", row.names = FALSE, append = TRUE)
  
  ### visualize the frequency and proportion of some interesting clonotypes with bar plots
  top10_carts <- t(car_frequency_over_time[1:10,5:14])
  top10_carts <- data.frame(top10_carts, Time=rownames(top10_carts),
                            stringsAsFactors = FALSE, check.names = FALSE)
  top10_carts$Time <- factor(top10_carts$Time, levels = levels(Seurat_Obj@meta.data$TimeF))
  p <- vector("list", length = ncol(top10_carts)-1)
  names(p) <- colnames(top10_carts)[1:length(p)]
  for(clonotype in names(p)) {
    p[[clonotype]] <- ggplot(top10_carts, aes_string(x="Time", y=clonotype)) +
      labs(x="", y="Clone Size") +
      geom_bar(fill = "grey50", position = "dodge", stat = "identity") +
      ggtitle(paste0(clonotype)) +
      guides(fill=guide_legend(title=NULL)) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 14),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14),
            axis.title.y = element_text(size = 14))
  }
  
  ### arrange the plots and save
  fName <- paste0("Top10_CAR+_Tcells_Lineage_Tracing")
  g <- arrangeGrob(grobs = p,
                   nrow = 5,
                   ncol = 2,
                   top = fName)
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 12, height = 12, dpi = 300)
  
  
  ### clone size of all the CAR+ T cells on UMAP
  ### CAR+ T cells which have TCR info
  
  ### set min-max values for x-y axes
  x_range <- c(min(Seurat_Obj@meta.data$Full_UMAP1), max(Seurat_Obj@meta.data$Full_UMAP1))
  y_range <- c(min(Seurat_Obj@meta.data$Full_UMAP2), max(Seurat_Obj@meta.data$Full_UMAP2))
  
  ### One plot for each group
  unique_time <- levels(Seurat_Obj@meta.data$TimeF)
  p <- vector("list", length(unique_time))
  names(p) <- unique_time
  for(time in unique_time) {
    ### data creation
    temp <- Seurat_Obj@meta.data[which(Seurat_Obj@meta.data$Time == time),]
    temp <- merge(temp, car_frequency_over_time[,c("Clonotype", time),drop=FALSE],
                  by.x = "global_clonotype_strict", by.y = "Clonotype",
                  all.x = TRUE)
    plot_df <- data.frame(X=temp$Full_UMAP1,
                          Y=temp$Full_UMAP2,
                          Time=temp$TimeF,
                          Clonotype=temp$global_clonotype_strict,
                          Clone_Size=temp[,time])
    plot_df$Clone_Size[intersect(which(is.na(plot_df$Clone_Size)),
                                 intersect(which(!is.na(temp$cdr3_aa)),
                                           which(temp$CAR == "CARpos")))] <- 1
    plot_df <- plot_df[which(!is.na(plot_df$Clone_Size)),]
    
    ### UMAP for each time point
    if(sum(plot_df$Clone_Size) == 0) {
      p[[time]] <- NULL
    } else {
      p[[time]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="Clone_Size"), size=1.5, alpha=0.6) +
        xlab("UMAP1") + ylab("UMAP2") +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle(time) +
        theme_classic(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_gradient(low="blue", high="red")
    }
  }
  
  ### arrange the plots and save
  fName <- paste0("UMAP_CAR+_Tcells_Clone_Size")
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   top = fName)
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 20, height = 12, dpi = 300)
  
  
  ### the number of lineages - considering 1. alpha & beta, 2. beta only, 3. alpha only
  ### Lineages: the number of lineages that appeared more than one time points in the Seurat object
  ### Lineages_CARpos: the number of lineages in CAR+ cells that appeared more than one time points in the Seurat object
  ### Clonotypes: the number of clonotypes in the Seurat object
  ### Clonotypes_CARpos: the number of clonotypes in CAR+ cells in the Seurat object
  ### Filtered_Clonotypes: This is similar to the "Clonotypes" but unique (that exist only one in the data)
  ###                    clonotypes were filtered out 
  lineage_abstract_table <- matrix(0, 5, 4)
  rownames(lineage_abstract_table) <- c("Lineages", "Lineages_CARpos",
                                        "Clonotypes", "Clonotypes_CARpos", "Filtered_Clonotypes")
  colnames(lineage_abstract_table) <- c("Alpha&Beta", "Alpha_Only", "Beta_Only", "Only_one_Aplha&Beta")
  
  only_one_strict_clonotypes_idx <- intersect(which(!is.na(Seurat_Obj@meta.data$cdr3_aa)),
                                              which(is.na(Seurat_Obj@meta.data$global_clonotype_strict)))
  only_one_lenient_clonotypes_idx <- intersect(which(!is.na(Seurat_Obj@meta.data$cdr3_aa)),
                                              which(is.na(Seurat_Obj@meta.data$global_clonotype_lenient)))
  only_one_strict_clonotypes_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                     only_one_strict_clonotypes_idx)
  only_one_lenient_clonotypes_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                      only_one_lenient_clonotypes_idx)
  
  ### Lineages
  ab <- unique(Seurat_Obj@meta.data$global_clonotype_strict)
  ab <- ab[!is.na(ab)]
  cnt <- 0
  for(i in 1:length(ab)) {
    if(length(unique(Seurat_Obj@meta.data$Time[which(Seurat_Obj@meta.data$global_clonotype_strict == ab[i])])) > 1) {
      cnt <- cnt+1
    }
  }
  lineage_abstract_table["Lineages","Alpha&Beta"] <- cnt
  
  b <- unique(Seurat_Obj@meta.data$global_clonotype_lenient)
  b <- b[!is.na(b)]
  cnt <- 0
  for(i in 1:length(b)) {
    if(length(unique(Seurat_Obj@meta.data$Time[which(Seurat_Obj@meta.data$global_clonotype_lenient == b[i])])) > 1) {
      cnt <- cnt+1
    }
  }
  lineage_abstract_table["Lineages","Beta_Only"] <- cnt
  
  ### Lineages - CAR+
  cnt <- 0
  for(i in 1:length(ab)) {
    if(length(unique(Seurat_Obj@meta.data$Time[intersect(which(Seurat_Obj@meta.data$global_clonotype_strict == ab[i]),
                                                         which(Seurat_Obj@meta.data$CAR == "CARpos"))])) > 1) {
      cnt <- cnt+1
    }
  }
  lineage_abstract_table["Lineages_CARpos","Alpha&Beta"] <- cnt
  
  cnt <- 0
  for(i in 1:length(b)) {
    if(length(unique(Seurat_Obj@meta.data$Time[intersect(which(Seurat_Obj@meta.data$global_clonotype_lenient == b[i]),
                                                         which(Seurat_Obj@meta.data$CAR == "CARpos"))])) > 1) {
      cnt <- cnt+1
    }
  }
  lineage_abstract_table["Lineages_CARpos","Beta_Only"] <- cnt
  
  ### Clonotypes
  lineage_abstract_table["Clonotypes","Alpha&Beta"] <- length(ab) +
    length(only_one_strict_clonotypes_idx)
    
  lineage_abstract_table["Clonotypes","Beta_Only"] <- length(b) +
    length(only_one_lenient_clonotypes_idx)
  
  ### Clonotypes - CAR+
  lineage_abstract_table["Clonotypes_CARpos","Alpha&Beta"] <- length(unique(Seurat_Obj@meta.data$global_clonotype_strict[which(Seurat_Obj@meta.data$CAR == "CARpos")])) +
    length(only_one_strict_clonotypes_carpos_idx)
  
  lineage_abstract_table["Clonotypes_CARpos","Beta_Only"] <- length(unique(Seurat_Obj@meta.data$global_clonotype_lenient[which(Seurat_Obj@meta.data$CAR == "CARpos")])) +
    length(only_one_lenient_clonotypes_carpos_idx)
  
  ### Filtered_Clonotypes
  lineage_abstract_table["Filtered_Clonotypes","Alpha&Beta"] <- length(ab)
  lineage_abstract_table["Filtered_Clonotypes","Beta_Only"] <- length(b)
  
  
  ### New task after the meeting at Apr 20
  
  
  
  
  
  
}
