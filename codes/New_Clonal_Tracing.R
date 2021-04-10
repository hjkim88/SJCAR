###
#   File name : New_Clonal_Tracing.R
#   Author    : Hyunjin Kim
#   Date      : Apr 5, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Define various clonotypes and evaluate the results
#
#   Instruction
#               1. Source("New_Clonal_Tracing.R")
#               2. Run the function "new_clonal_tracing" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_New_Clonal_Tracing.R/New_Clonal_Tracing.R")
#               > new_clonal_tracing(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                    outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/New_Tracing/")
###

new_clonal_tracing <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
                               outputDir="./results/New3/New_Tracing/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
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
  
  ### various clonotyping approaches
  Seurat_Obj@meta.data$cdr3_alpha <- NA
  Seurat_Obj@meta.data$cdr3_alpha_umis <- NA
  Seurat_Obj@meta.data$cdr3_beta <- NA
  Seurat_Obj@meta.data$cdr3_beta_umis <- NA
  Seurat_Obj@meta.data$cdr3_one_alpha_beta <- NA
  for(i in 1:nrow(Seurat_Obj@meta.data)) {
    
    ### split CDR3 AA
    aa_split <- strsplit(Seurat_Obj@meta.data$cdr3_aa[i], split = ";", fixed = TRUE)[[1]]
    umi_split <- as.numeric(strsplit(Seurat_Obj@meta.data$tcr_umis[i], split = ";", fixed = TRUE)[[1]])
    
    ### get alpha & beta indicies
    alpha_idx <- which(startsWith(aa_split, "TRA"))
    beta_idx <- which(startsWith(aa_split, "TRB"))
    
    ### fill out the columns
    Seurat_Obj@meta.data$cdr3_alpha[i] <- paste(aa_split[alpha_idx], collapse = ";")
    Seurat_Obj@meta.data$cdr3_alpha_umis[i] <- paste(umi_split[alpha_idx], collapse = ";")
    Seurat_Obj@meta.data$cdr3_beta[i] <- paste(aa_split[beta_idx], collapse = ";")
    Seurat_Obj@meta.data$cdr3_beta_umis[i] <- paste(umi_split[beta_idx], collapse = ";")
    if(length(alpha_idx) > 0 && length(beta_idx) > 0) {
      best_alpha_idx <- alpha_idx[which(umi_split[alpha_idx] == max(umi_split[alpha_idx]))]
      best_beta_idx <- beta_idx[which(umi_split[beta_idx] == max(umi_split[beta_idx]))]
      Seurat_Obj@meta.data$cdr3_one_alpha_beta[i] <- paste0(aa_split[best_alpha_idx], ";", aa_split[best_beta_idx])
    } else if(length(alpha_idx) > 0 && length(beta_idx) == 0) {
      best_alpha_idx <- alpha_idx[which(umi_split[alpha_idx] == max(umi_split[alpha_idx]))]
      Seurat_Obj@meta.data$cdr3_one_alpha_beta[i] <- aa_split[best_alpha_idx]
    } else if(length(alpha_idx) == 0 && length(beta_idx) > 0) {
      best_beta_idx <- beta_idx[which(umi_split[beta_idx] == max(umi_split[beta_idx]))]
      Seurat_Obj@meta.data$cdr3_one_alpha_beta[i] <- aa_split[best_beta_idx]
    }
    
    ### progress
    if(i %% 10000 == 0) {
      writeLines(paste(i, "/", nrow(Seurat_Obj@meta.data)))
    }
    
  }
  
  ### change "" to NA
  Seurat_Obj@meta.data$cdr3_alpha[which(Seurat_Obj@meta.data$cdr3_alpha == "")] <- NA
  Seurat_Obj@meta.data$cdr3_alpha_umis[which(Seurat_Obj@meta.data$cdr3_alpha_umis == "")] <- NA
  Seurat_Obj@meta.data$cdr3_beta[which(Seurat_Obj@meta.data$cdr3_beta == "")] <- NA
  Seurat_Obj@meta.data$cdr3_beta_umis[which(Seurat_Obj@meta.data$cdr3_beta_umis == "")] <- NA
  Seurat_Obj@meta.data$cdr3_one_alpha_beta[which(Seurat_Obj@meta.data$cdr3_one_alpha_beta == "")] <- NA
  
  ### global clonotyping for each patient
  Seurat_Obj@meta.data$clonotype_id_by_patient_alpha <- NA
  Seurat_Obj@meta.data$clonotype_id_by_patient_beta <- NA
  Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta <- NA
  for(px in unique(Seurat_Obj@meta.data$px)) {
    ### print progress
    writeLines(paste(px))
    
    #
    ### alpha only
    #
    ### patient indices that have TCR info
    px_idx <- intersect(which(Seurat_Obj@meta.data$px == px),
                        which(!is.na(Seurat_Obj@meta.data$cdr3_alpha)))
    
    ### get unique and duplicated indicies
    unique_idx <- which(!duplicated(Seurat_Obj@meta.data$cdr3_alpha[px_idx]))
    dup_idx <- which(duplicated(Seurat_Obj@meta.data$cdr3_alpha[px_idx]))
    
    ### give clone ids
    Seurat_Obj@meta.data$clonotype_id_by_patient_alpha[px_idx[unique_idx]] <- paste0(px, "_Clone_", 1:length(unique_idx))
    for(i in px_idx[dup_idx]) {
      target_idx <- which(Seurat_Obj@meta.data$cdr3_alpha[px_idx[unique_idx]] == Seurat_Obj@meta.data$cdr3_alpha[i])
      Seurat_Obj@meta.data$clonotype_id_by_patient_alpha[i] <- Seurat_Obj@meta.data$clonotype_id_by_patient_alpha[px_idx[unique_idx][target_idx]]
    }
    
    #
    ### beta only
    #
    ### patient indices that have TCR info
    px_idx <- intersect(which(Seurat_Obj@meta.data$px == px),
                        which(!is.na(Seurat_Obj@meta.data$cdr3_beta)))
    
    ### get unique and duplicated indicies
    unique_idx <- which(!duplicated(Seurat_Obj@meta.data$cdr3_beta[px_idx]))
    dup_idx <- which(duplicated(Seurat_Obj@meta.data$cdr3_beta[px_idx]))
    
    ### give clone ids
    Seurat_Obj@meta.data$clonotype_id_by_patient_beta[px_idx[unique_idx]] <- paste0(px, "_Clone_", 1:length(unique_idx))
    for(i in px_idx[dup_idx]) {
      target_idx <- which(Seurat_Obj@meta.data$cdr3_beta[px_idx[unique_idx]] == Seurat_Obj@meta.data$cdr3_beta[i])
      Seurat_Obj@meta.data$clonotype_id_by_patient_beta[i] <- Seurat_Obj@meta.data$clonotype_id_by_patient_beta[px_idx[unique_idx][target_idx]]
    }
    
    #
    ### one alpha & beta per cell
    #
    ### patient indices that have TCR info
    px_idx <- intersect(which(Seurat_Obj@meta.data$px == px),
                        which(!is.na(Seurat_Obj@meta.data$cdr3_one_alpha_beta)))
    
    ### get unique and duplicated indicies
    unique_idx <- which(!duplicated(Seurat_Obj@meta.data$cdr3_one_alpha_beta[px_idx]))
    dup_idx <- which(duplicated(Seurat_Obj@meta.data$cdr3_one_alpha_beta[px_idx]))
    
    ### give clone ids
    Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[px_idx[unique_idx]] <- paste0(px, "_Clone_", 1:length(unique_idx))
    for(i in px_idx[dup_idx]) {
      target_idx <- which(Seurat_Obj@meta.data$cdr3_one_alpha_beta[px_idx[unique_idx]] == Seurat_Obj@meta.data$cdr3_one_alpha_beta[i])
      Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[i] <- Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[px_idx[unique_idx][target_idx]]
    }
    
  }
  
  ### lineage data
  lineages <- vector("list", length = 4)
  names(lineages) <- c("Strict", "Alpha", "Beta", "One_From_Each")
  
  ### making lineage table
  for(type in names(lineages)) {
    
    ### progress
    writeLines(paste(type))
    
    ### set clonotype column based on the given type
    if(type == "Strict") {
      clonotype_col_name <- "clonotype_id_by_patient"
    } else if(type == "Alpha") {
      clonotype_col_name <- "clonotype_id_by_patient_alpha"
    } else if(type == "Beta") {
      clonotype_col_name <- "clonotype_id_by_patient_beta"
    } else {
      clonotype_col_name <- "clonotype_id_by_patient_one_alpha_beta"
    }
    
    ### with all the cells and with car-pos cells only
    lineages[[type]] <- vector("list", length = 2)
    names(lineages[[type]]) <- c("All", "CARpos")
    
    ### for each cell type
    for(cell_type in names(lineages[[type]])) {
      
      ### empty list for each patient
      lineages[[type]][[cell_type]] <- vector("list", length = length(unique(Seurat_Obj@meta.data$px)))
      names(lineages[[type]][[cell_type]]) <- unique(Seurat_Obj@meta.data$px)
      
      ### for each patient
      for(px in unique(Seurat_Obj@meta.data$px)) {
        
        ### patient indicies
        px_idx <- which(Seurat_Obj@meta.data$px == px)
        
        ### get time points for the given patient
        px_time_points <- unique(Seurat_Obj@meta.data$time[px_idx])
        px_time_points <- px_time_points[order(factor(px_time_points, levels = total_time_points))]
        
        ### clonotype table and lineage table
        lineages[[type]][[cell_type]][[px]] <- vector("list", length = 3)
        names(lineages[[type]][[cell_type]][[px]]) <- c("Cell_#", "Clonotype_#", "Lineage_#")
        
        ### number of cells
        lineages[[type]][[cell_type]][[px]][["Cell_#"]] <- rep(0, length(px_time_points))
        names(lineages[[type]][[cell_type]][[px]][["Cell_#"]]) <- px_time_points
        for(tp in px_time_points) {
          if(cell_type == "All") {
            lineages[[type]][[cell_type]][[px]][["Cell_#"]][tp] <- length(intersect(px_idx,
                                                                                    which(Seurat_Obj@meta.data$time == tp)))
          } else {
            lineages[[type]][[cell_type]][[px]][["Cell_#"]][tp] <- length(intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                                                    intersect(px_idx,
                                                                                              which(Seurat_Obj@meta.data$time == tp))))
          }
        }
        
        ### number of unique clonotypes
        lineages[[type]][[cell_type]][[px]][["Clonotype_#"]] <- rep(0, length(px_time_points))
        names(lineages[[type]][[cell_type]][[px]][["Clonotype_#"]]) <- px_time_points
        for(tp in px_time_points) {
          if(cell_type == "All") {
            lineages[[type]][[cell_type]][[px]][["Clonotype_#"]][tp] <- length(unique(Seurat_Obj@meta.data[,clonotype_col_name][intersect(px_idx,
                                                                                                                                          which(Seurat_Obj@meta.data$time == tp))]))
          } else {
            lineages[[type]][[cell_type]][[px]][["Clonotype_#"]][tp] <- length(unique(Seurat_Obj@meta.data[,clonotype_col_name][intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                                                                                                          intersect(px_idx,
                                                                                                                                                    which(Seurat_Obj@meta.data$time == tp)))]))
          }
        }
        
        ### get index after GMP
        gmp_idx <- which(px_time_points == "GMP")
        gmp_redo_idx <- which(px_time_points == "GMP-redo")
        after_gmp_idx <- max(gmp_idx, gmp_redo_idx) + 1
        
        ### an empty vector
        lineages[[type]][[cell_type]][[px]][["Lineage_#"]] <- rep(0, length(px_time_points))
        names(lineages[[type]][[cell_type]][[px]][["Lineage_#"]]) <- px_time_points
        
        if((after_gmp_idx != -Inf) && (after_gmp_idx <= length(px_time_points))) {
          ### number of unique lineages
          for(tp in px_time_points[after_gmp_idx:length(px_time_points)]) {
            if(cell_type == "All") {
              lineages[[type]][[cell_type]][[px]][["Lineage_#"]][tp] <- length(unique(intersect(Seurat_Obj@meta.data[,clonotype_col_name][intersect(px_idx,
                                                                                                                                                    which(Seurat_Obj@meta.data$time == tp))],
                                                                                                Seurat_Obj@meta.data[,clonotype_col_name][intersect(px_idx,
                                                                                                                                                    which(Seurat_Obj@meta.data$time %in% c("GMP", "GMP-redo")))])))
            } else {
              lineages[[type]][[cell_type]][[px]][["Lineage_#"]][tp] <- length(unique(intersect(Seurat_Obj@meta.data[,clonotype_col_name][intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                                                                                                                    intersect(px_idx,
                                                                                                                                                              which(Seurat_Obj@meta.data$time == tp)))],
                                                                                                Seurat_Obj@meta.data[,clonotype_col_name][intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                                                                                                                    intersect(px_idx,
                                                                                                                                                              which(Seurat_Obj@meta.data$time %in% c("GMP", "GMP-redo"))))])))
            }
          }
        }
        
      }
      
    }
    
  }
  
  ### data frame for bar plot - all the cells
  plot_df <- data.frame(matrix(0, length(lineages)*length(unique(Seurat_Obj@meta.data$px)), 5))
  colnames(plot_df) <- c("Cell_Num", "Clonotype_Num", "Lineage_Num", "Type", "Patient")
  plot_df[,"Cell_Num"] <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) {
    temp <- NULL
    for(type in names(lineages)) {
      temp <- c(temp, sum(lineages[[type]][["All"]][[x]][["Cell_#"]]))
    }
    return(temp)
  }))
  plot_df[,"Clonotype_Num"] <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) {
    temp <- NULL
    for(type in names(lineages)) {
      temp <- c(temp, sum(lineages[[type]][["All"]][[x]][["Clonotype_#"]]))
    }
    return(temp)
  }))
  plot_df[,"Lineage_Num"] <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) {
    temp <- NULL
    for(type in names(lineages)) {
      temp <- c(temp, sum(lineages[[type]][["All"]][[x]][["Lineage_#"]]))
    }
    return(temp)
  }))
  plot_df[,"Type"] <- rep(names(lineages), length(unique(Seurat_Obj@meta.data$px)))
  plot_df[,"Patient"] <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) rep(x, length(lineages))))
  
  ### bar plot for each patient
  p <- vector("list", 2)
  p[[1]] <- ggplot(plot_df, aes_string(x="Patient", y="Clonotype_Num", fill="Type")) +
    labs(x="", y="Unique Clonotype # (All)") +
    geom_bar(position = "dodge", stat = "identity") +
    ggtitle("The Number of Unique Clonotypes from All the Cells") +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#EC6258")) +
    guides(fill=guide_legend(title=NULL)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size = 25),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 35))
  p[[2]] <- ggplot(plot_df, aes_string(x="Patient", y="Lineage_Num", fill="Type")) +
    labs(x="", y="Unique Lineage # (All)") +
    geom_bar(position = "dodge", stat = "identity") +
    ggtitle("The Number of Unique Lineages from All the Cells") +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#EC6258")) +
    guides(fill=guide_legend(title=NULL)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size = 25),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 35))
  
  ### arrange the plots and save
  fName <- paste0("Clonotypes_and_Lineages_in_All_Cells")
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   top = "")
  ggsave(file = paste0(outputDir, fName, ".png"), g, width = 20, height = 10, dpi = 300)
  
  
  ### data frame for bar plot - CAR+ cells only
  plot_df <- data.frame(matrix(0, length(lineages)*length(unique(Seurat_Obj@meta.data$px)), 5))
  colnames(plot_df) <- c("Cell_Num", "Clonotype_Num", "Lineage_Num", "Type", "Patient")
  plot_df[,"Cell_Num"] <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) {
    temp <- NULL
    for(type in names(lineages)) {
      temp <- c(temp, sum(lineages[[type]][["CARpos"]][[x]][["Cell_#"]]))
    }
    return(temp)
  }))
  plot_df[,"Clonotype_Num"] <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) {
    temp <- NULL
    for(type in names(lineages)) {
      temp <- c(temp, sum(lineages[[type]][["CARpos"]][[x]][["Clonotype_#"]]))
    }
    return(temp)
  }))
  plot_df[,"Lineage_Num"] <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) {
    temp <- NULL
    for(type in names(lineages)) {
      temp <- c(temp, sum(lineages[[type]][["CARpos"]][[x]][["Lineage_#"]]))
    }
    return(temp)
  }))
  plot_df[,"Type"] <- rep(names(lineages), length(unique(Seurat_Obj@meta.data$px)))
  plot_df[,"Patient"] <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) rep(x, length(lineages))))
  
  ### bar plot for each patient
  p <- vector("list", 2)
  p[[1]] <- ggplot(plot_df, aes_string(x="Patient", y="Clonotype_Num", fill="Type")) +
    labs(x="", y="Unique Clonotype # (CAR+)") +
    geom_bar(position = "dodge", stat = "identity") +
    ggtitle("The Number of Unique Clonotypes from CAR+ Cells") +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#EC6258")) +
    guides(fill=guide_legend(title=NULL)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size = 25),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 35))
  p[[2]] <- ggplot(plot_df, aes_string(x="Patient", y="Lineage_Num", fill="Type")) +
    labs(x="", y="Unique Lineage # (CAR+)") +
    geom_bar(position = "dodge", stat = "identity") +
    ggtitle("The Number of Unique Lineages from CAR+ Cells") +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#EC6258")) +
    guides(fill=guide_legend(title=NULL)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size = 25),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 35))
  
  ### arrange the plots and save
  fName <- paste0("Clonotypes_and_Lineages_in_CARpos_Cells")
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   top = "")
  ggsave(file = paste0(outputDir, fName, ".png"), g, width = 20, height = 10, dpi = 300)
  
  
  ### for each patient and for each time point - clonotype # for CAR+ cells
  p <- vector("list", length = length(unique(Seurat_Obj@meta.data$px)))
  names(p) <- unique(Seurat_Obj@meta.data$px)
  
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    ### patient indicies
    px_idx <- which(Seurat_Obj@meta.data$px == px)
    
    ### get time points for the given patient
    px_time_points <- unique(Seurat_Obj@meta.data$time[px_idx])
    px_time_points <- px_time_points[order(factor(px_time_points, levels = total_time_points))]
    
    ### make an empty matrix for plot
    plot_df <- data.frame(matrix(NA, length(px_time_points)*length(lineages), 3),
                          stringsAsFactors = FALSE, check.names = FALSE)
    colnames(plot_df) <- c("Number", "Time", "Type")
    
    ### for each time point, fill out the matrix
    for(i in 1:length(px_time_points)) {
      ### for each clonotyping type
      for(j in 1:length(lineages)) {
        ### the number of clonotypes
        plot_df[(i-1)*length(lineages)+j,"Number"] <- lineages[[j]][["CARpos"]][[patient]][["Clonotype_#"]][px_time_points[i]]
      }
    }
    plot_df[,"Time"] <- as.vector(sapply(px_time_points, function(x) rep(x, length(lineages))))
    plot_df[,"Type"] <- rep(names(lineages), length(px_time_points))
    
    ### make 0 values to NA
    plot_df$Number[which(plot_df$Number == 0)] <- NA
    
    ### make "Time" factorized
    plot_df$Time <- factor(plot_df$Time, levels = total_time_points)
    
    ### save the plot to the list
    p[[patient]] <- ggplot(plot_df, aes_string(x="Time", y="Number", fill="Type", group="Type")) +
      labs(x="", y="The Number of Clonotypes (CAR+)") +
      geom_bar(position = "dodge", stat = "identity") +
      geom_text(aes_string(label="Number", color="Type", group="Type"),
                position=position_dodge(width=1), size=5, angle = 45, hjust=0, vjust=0,
                show.legend = FALSE) +
      ggtitle(patient) +
      scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#EC6258")) +
      scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#EC6258")) +
      guides(fill=guide_legend(title=NULL)) +
      ylim(0, max(plot_df$Number, na.rm = TRUE)*1.15) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size = 25),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35),
            axis.title.y = element_text(size = 25),
            legend.position="top",
            legend.text = element_text(size = 20))
  }
  
  ### arrange the plots and save
  fName <- "Clonotypes_per_TP_CARpos1"
  g <- arrangeGrob(grobs = p[1:6],
                   ncol = 3,
                   top = "")
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 24, height = 12, dpi = 300)
  
  ### arrange the plots and save
  fName <- "Clonotypes_per_TP_CARpos2"
  g <- arrangeGrob(grobs = p[7:12],
                   ncol = 3,
                   top = "")
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 24, height = 12, dpi = 300)
  
  
  ### for each patient and for each time point - lineage # for CAR+ cells
  p <- vector("list", length = length(unique(Seurat_Obj@meta.data$px)))
  names(p) <- unique(Seurat_Obj@meta.data$px)
  
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    ### patient indicies
    px_idx <- which(Seurat_Obj@meta.data$px == px)
    
    ### get time points for the given patient
    px_time_points <- unique(Seurat_Obj@meta.data$time[px_idx])
    px_time_points <- px_time_points[order(factor(px_time_points, levels = total_time_points))]
    
    ### make an empty matrix for plot
    plot_df <- data.frame(matrix(NA, length(px_time_points)*length(lineages), 3),
                          stringsAsFactors = FALSE, check.names = FALSE)
    colnames(plot_df) <- c("Number", "Time", "Type")
    
    ### for each time point, fill out the matrix
    for(i in 1:length(px_time_points)) {
      ### for each clonotyping type
      for(j in 1:length(lineages)) {
        ### the number of lineages
        plot_df[(i-1)*length(lineages)+j,"Number"] <- lineages[[j]][["CARpos"]][[patient]][["Lineage_#"]][px_time_points[i]]
      }
    }
    plot_df[,"Time"] <- as.vector(sapply(px_time_points, function(x) rep(x, length(lineages))))
    plot_df[,"Type"] <- rep(names(lineages), length(px_time_points))
    
    ### make 0 values to NA
    plot_df$Number[which(plot_df$Number == 0)] <- NA
    
    ### make "Time" factorized
    plot_df$Time <- factor(plot_df$Time, levels = total_time_points)
    
    ### save the plot to the list
    p[[patient]] <- ggplot(plot_df, aes_string(x="Time", y="Number", fill="Type", group="Type")) +
      labs(x="", y="The Number of Lineages (CAR+)") +
      geom_bar(position = "dodge", stat = "identity") +
      geom_text(aes_string(label="Number", color="Type", group="Type"),
                position=position_dodge(width=1), size=5, angle = 45, hjust=0, vjust=0,
                show.legend = FALSE) +
      ggtitle(patient) +
      scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#EC6258")) +
      scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#EC6258")) +
      guides(fill=guide_legend(title=NULL)) +
      ylim(0, max(plot_df$Number, na.rm = TRUE)*1.15) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1, size = 25),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35),
            axis.title.y = element_text(size = 25),
            legend.position="top",
            legend.text = element_text(size = 20))
  }
  
  ### arrange the plots and save
  fName <- "Lineages_per_TP_CARpos1"
  g <- arrangeGrob(grobs = p[1:6],
                   ncol = 3,
                   top = "")
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 24, height = 12, dpi = 300)
  
  ### arrange the plots and save
  fName <- "Lineages_per_TP_CARpos2"
  g <- arrangeGrob(grobs = p[7:12],
                   ncol = 3,
                   top = "")
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 24, height = 12, dpi = 300)
  
  
  ### save the Seurat object
  saveRDS(Seurat_Obj, file = Seurat_RObj_path)
  
}
