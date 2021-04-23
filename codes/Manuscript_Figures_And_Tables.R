###
#   File name : Manuscript_Figures_And_Tables.R
#   Author    : Hyunjin Kim
#   Date      : Apr 15, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : 1. Alluvial plots of lineages - CAR+ a) all patients into one b) separated by patient
#               2. Number of lineages & clonotypes table CD4+, CD8+, CAR+, CAR-, per patient
#                  (alpha only, beta only, one from each, strict)
#               3. CAR+ % for each library - proportional bar graph & table - CAR+, CAR-, CD4+, CD8+
#               4. Clone size between CAR+ lineages vs non-lineage CAR+ after infusion
#               5. a) PCA/UMAP plot of CAR+s over time (coloring based on time)
#                  b) Clustering + UMAP 
#                  c) Pseudotime analysis on PCA
#                  d) PCA/UMAP plot of lineages with size=1 vs lineages with size > 1 (coloring differently)
#               6. Time series DE analysis & pathway analysis 
#
#   Instruction
#               1. Source("Manuscript_Figures_And_Tables.R")
#               2. Run the function "manuscript_prep" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Manuscript_Figures_And_Tables.R/Manuscript_Figures_And_Tables.R")
#               > manuscript_prep(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                 px_result_dir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/",
#                                 outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Manuscript/")
###

manuscript_prep <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj_Total2.RDS",
                            px_result_dir="./results/New3/",
                            outputDir="./results/New3/Manuscript/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
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
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
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
  
  ### combine some seprated time points into one
  Seurat_Obj@meta.data$time2 <- Seurat_Obj@meta.data$time
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "GMP-redo")] <- "GMP"
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "PreTransB")] <- "PreTrans"
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "Wk1b")] <- "Wk1"
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "Wk-1Run1")] <- "Wk-1"
  Seurat_Obj@meta.data$time2[which(Seurat_Obj@meta.data$time2 == "Wk-1Run2")] <- "Wk-1"
  
  ### set time points
  total_time_points <- c("PreTrans", "Wk-1", "Wk0", "GMP", "Wk1", "Wk2", 
                         "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo")
  gmp_after_time_points <- c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", 
                             "Wk6", "Wk8", "3mo", "6mo", "9mo")
  
  ### theme that draws dotted lines for each y-axis ticks
  ### this function is from "immunarch" package
  theme_cleveland2 <- function(rotate = TRUE) {
    if (rotate) {
      theme(
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(
          colour = "grey70",
          linetype = "dashed"
        )
      )
    }
    else {
      theme(
        panel.grid.major.x = element_line(
          colour = "grey70",
          linetype = "dashed"
        ), panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )
    }
  }
  
  #
  ### 1. Alluvial plots of lineages - CAR+ a) all patients into one b) separated by patient
  #
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/1/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### run for each patient
  total_plot_df <- NULL
  p <- vector("list", length= length(unique(Seurat_Obj@meta.data$px)))
  names(p) <- unique(Seurat_Obj@meta.data$px)
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    
    ### print progress
    writeLines(paste(patient))
    
    ### load the file
    target_file <- read.xlsx2(file = paste0(px_result_dir, patient, "/car_clonotype_frequency_over_time_", patient, ".xlsx"),
                              sheetName = paste0("CARpos_Clonotype_Frequency_One_"), stringsAsFactors = FALSE, check.names = FALSE,
                              row.names = 1)
    
    ### numerize the table
    for(i in 1:ncol(target_file)) {
      target_file[,i] <- as.numeric(target_file[,i])
    }
    
    ### combine some redundant time points to one
    if(length(which(colnames(target_file) == "GMP-redo")) > 0) {
      target_file[,"GMP"] <- target_file[,"GMP"] + target_file[,"GMP-redo"]
      target_file <- target_file[,-which(colnames(target_file) == "GMP-redo")]
    }
    if(length(which(colnames(target_file) == "PreTransB")) > 0) {
      if(length(which(colnames(target_file) == "PreTrans")) > 0) {
        target_file[,"PreTrans"] <- target_file[,"PreTrans"] + target_file[,"PreTransB"]
        target_file <- target_file[,-which(colnames(target_file) == "PreTransB")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "PreTransB")] <- "PreTrans"
      }
    }
    if(length(which(colnames(target_file) == "Wk1b")) > 0) {
      if(length(which(colnames(target_file) == "Wk1")) > 0) {
        target_file[,"Wk1"] <- target_file[,"Wk1"] + target_file[,"Wk1b"]
        target_file <- target_file[,-which(colnames(target_file) == "Wk1b")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "Wk1b")] <- "Wk1"
      }
    }
    if(length(which(colnames(target_file) == "Wk-1Run1")) > 0) {
      if(length(which(colnames(target_file) == "Wk-1")) > 0) {
        target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run1"]
        target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run1")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "Wk-1Run1")] <- "Wk-1"
      }
    }
    if(length(which(colnames(target_file) == "Wk-1Run2")) > 0) {
      if(length(which(colnames(target_file) == "Wk-1")) > 0) {
        target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run2"]
        target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run2")]
      } else {
        colnames(target_file)[which(colnames(target_file) == "Wk-1Run2")] <- "Wk-1"
      }
    }
    
    ### remove all zero time points
    time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
    
    ### get time points except the Total
    time_points <- setdiff(time_points, c("Total"))
    
    ### remove before GMP time points
    time_points <- intersect(time_points, gmp_after_time_points)
    
    ### draw when there are at least two time points
    if(length(time_points) > 1) {
      ###  get lineages
      lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
        return(length(which(x > 0)) > 1)  
      })),time_points]
      
      ### get an input data frame for the alluvial plot
      total_rows <- length(which(lineage_table[,time_points] > 0))
      if(total_rows > 0) {
        plot_df <- data.frame(Time=rep("", total_rows),
                              Clone_Size=rep(0, total_rows),
                              Clone=rep("", total_rows),
                              CDR3=rep("", total_rows))
        cnt <- 1
        for(i in 1:nrow(lineage_table)) {
          for(tp in time_points) {
            if(lineage_table[i,tp] > 0) {
              plot_df[cnt,] <- c(tp,
                                 lineage_table[i,tp],
                                 rownames(lineage_table)[i],
                                 "CDR3")
              cnt <- cnt + 1
            }
          }
        }
        plot_df$Time <- factor(plot_df$Time, levels = intersect(time_points, unique(plot_df$Time)))
        
        ### numerize the clone_size column
        plot_df$Clone_Size <- as.numeric(plot_df$Clone_Size)
        
        ### draw an alluvial plot
        p[[patient]] <- ggplot(plot_df,
               aes(x = Time, stratum = Clone, alluvium = Clone,
                   y = Clone_Size,
                   fill = Clone, label = Clone)) +
          ggtitle(paste(patient)) +
          geom_flow() +
          geom_stratum(alpha = 1) +
          # geom_text(stat = "stratum", size = 2) +
          rotate_x_text(90) +
          theme_pubr(legend = "none") +
          # theme_cleveland2() +
          scale_fill_viridis(discrete = T) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
          theme_classic(base_size = 36) +
          theme(axis.text.x = element_text(size = 30),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 30),
                legend.position = "none")
        ggsave(file = paste0(outputDir2, "Car+_Clonal_Tracing_", patient, ".png"),
               plot = p[[patient]],
               width = 20, height = 10, dpi = 350)
        
        ### combine the plot data into one
        if(is.null(total_plot_df)) {
          total_plot_df <- plot_df
        } else {
          total_plot_df <- rbind(total_plot_df, plot_df)
        }
      }
    }
    
    gc()
  }
  
  ### draw an alluvial plot with all-patient-combined plot data
  total_plot_df$Time <- factor(total_plot_df$Time, levels = intersect(total_time_points, unique(total_plot_df$Time)))
  ggplot(total_plot_df,
         aes(x = Time, stratum = Clone, alluvium = Clone,
             y = Clone_Size,
             fill = Clone, label = Clone)) +
    ggtitle(paste("All patients - Px00-Px15")) +
    geom_flow() +
    geom_stratum(alpha = 1) +
    # geom_text(stat = "stratum", size = 2) +
    rotate_x_text(90) +
    theme_pubr(legend = "none") +
    # theme_cleveland2() +
    scale_fill_viridis(discrete = T) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30),
          legend.position = "none")
  ggsave(file = paste0(outputDir2, "All_Car+_Clonal_Tracing.png"), width = 20, height = 10, dpi = 350)
  
  ### combine plos of the selected 9 patients into one
  g <- arrangeGrob(grobs = p[c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                               "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                               "SJCAR19-09", "SJCAR19-10", "SJCAR19-11")],
                   nrow = 3,
                   ncol = 3,
                   top = "")
  ggsave(file = paste0(outputDir2, "Car+_Clonal_Tracing_9pxs_in_one.png"), g, width = 40, height = 15, dpi = 350)
  
  
  ### 2. Number of lineages & clonotypes table CD4+, CD8+, CAR+, CAR-, per patient
  ### (alpha only, beta only, one from each, strict)
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/2/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### set table colnames & rownames
  table_colnames <- c("Alpha", "Beta", "One_From_Each", "Strict")
  clonotype_types <- c("clonotype_id_by_patient_alpha", "clonotype_id_by_patient_beta",
                       "clonotype_id_by_patient_one_alpha_beta", "clonotype_id_by_patient")
  names(table_colnames) <- clonotype_types
  names(clonotype_types) <- table_colnames
  table_rownames <- c("CD4+/CAR+ Cell #", "CD4+/CAR- Cell #",
                      "CD8+/CAR+ Cell #", "CD8+/CAR- Cell #",
                      "CD4+/CAR+ Clonotype #", "CD4+/CAR- Clonotype #",
                      "CD8+/CAR+ Clonotype #", "CD8+/CAR- Clonotype #",
                      "CD4+/CAR+ Lineage #", "CD4+/CAR- Lineage #",
                      "CD8+/CAR+ Lineage #", "CD8+/CAR- Lineage #",
                      "CD4+/CAR+ GMP Lineage #", "CD8+/CAR+ GMP Lineage #")
  px_table_rownames <- as.vector(sapply(unique(Seurat_Obj@meta.data$px), function(x) paste(x, table_rownames)))
  
  ### make an empty table
  result_table <- matrix(0, length(px_table_rownames), length(table_colnames))
  colnames(result_table) <- table_colnames
  rownames(result_table) <- px_table_rownames
  
  ### for each clonotype type to get lineages
  total_lineages <- vector("list", length = length(table_colnames))
  names(total_lineages) <- table_colnames
  total_gmp_lineages <- vector("list", length = length(table_colnames))
  names(total_gmp_lineages) <- table_colnames
  for(type in table_colnames) {
    total_lineages[[type]] <- NULL
    total_gmp_lineages[[type]] <- NULL
    ### run for each patient
    for(patient in unique(Seurat_Obj@meta.data$px)) {
      
      ### print progress
      writeLines(paste(patient))
      
      ### load the file
      target_file <- read.xlsx2(file = paste0(px_result_dir, patient, "/car_clonotype_frequency_over_time_", patient, ".xlsx"),
                                sheetName = paste0("CARpos_Clonotype_Frequency_", substr(type, 1, 4)), stringsAsFactors = FALSE, check.names = FALSE,
                                row.names = 1)
      
      ### numerize the table
      for(i in 1:ncol(target_file)) {
        target_file[,i] <- as.numeric(target_file[,i])
      }
      
      ### combine some redundant time points to one
      if(length(which(colnames(target_file) == "GMP-redo")) > 0) {
        target_file[,"GMP"] <- target_file[,"GMP"] + target_file[,"GMP-redo"]
        target_file <- target_file[,-which(colnames(target_file) == "GMP-redo")]
      }
      if(length(which(colnames(target_file) == "PreTransB")) > 0) {
        if(length(which(colnames(target_file) == "PreTrans")) > 0) {
          target_file[,"PreTrans"] <- target_file[,"PreTrans"] + target_file[,"PreTransB"]
          target_file <- target_file[,-which(colnames(target_file) == "PreTransB")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "PreTransB")] <- "PreTrans"
        }
      }
      if(length(which(colnames(target_file) == "Wk1b")) > 0) {
        if(length(which(colnames(target_file) == "Wk1")) > 0) {
          target_file[,"Wk1"] <- target_file[,"Wk1"] + target_file[,"Wk1b"]
          target_file <- target_file[,-which(colnames(target_file) == "Wk1b")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "Wk1b")] <- "Wk1"
        }
      }
      if(length(which(colnames(target_file) == "Wk-1Run1")) > 0) {
        if(length(which(colnames(target_file) == "Wk-1")) > 0) {
          target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run1"]
          target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run1")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "Wk-1Run1")] <- "Wk-1"
        }
      }
      if(length(which(colnames(target_file) == "Wk-1Run2")) > 0) {
        if(length(which(colnames(target_file) == "Wk-1")) > 0) {
          target_file[,"Wk-1"] <- target_file[,"Wk-1"] + target_file[,"Wk-1Run2"]
          target_file <- target_file[,-which(colnames(target_file) == "Wk-1Run2")]
        } else {
          colnames(target_file)[which(colnames(target_file) == "Wk-1Run2")] <- "Wk-1"
        }
      }
      
      ### remove all zero time points
      time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
      
      ### get time points except the Total
      time_points <- setdiff(time_points, c("Total"))
      
      ### draw when there are at least two time points
      if(length(time_points) > 1) {
        ###  get lineages
        lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
          return(length(which(x > 0)) > 1)  
        })),time_points]
        lineages <- rownames(lineage_table)
        
        ### get gmp lineages
        gmp_lineage_table <- lineage_table[which(lineage_table$GMP > 0),]
        gmp_lineages <- rownames(gmp_lineage_table)
        
        ### combine lineages
        if(length(lineages) > 0) {
          ### combine
          if(is.null(total_lineages[[type]])) {
            total_lineages[[type]] <- lineages
          } else {
            total_lineages[[type]] <- c(total_lineages[[type]], lineages)
          }
        }
        
        ### combine gmp lineages
        if(length(gmp_lineages) > 0) {
          ### combine
          if(is.null(total_gmp_lineages[[type]])) {
            total_gmp_lineages[[type]] <- gmp_lineages
          } else {
            total_gmp_lineages[[type]] <- c(total_gmp_lineages[[type]], gmp_lineages)
          }
        }
      }
      
      gc()
    }
  }
  
  ### some pre-calculated indices
  cd4_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4"),
                              which(Seurat_Obj@meta.data$CAR == "CARpos"))
  cd4_carneg_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4"),
                              which(Seurat_Obj@meta.data$CAR == "CARneg"))
  cd8_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"),
                              which(Seurat_Obj@meta.data$CAR == "CARpos"))
  cd8_carneg_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"),
                              which(Seurat_Obj@meta.data$CAR == "CARneg"))
  lineages_idx <- vector("list", length = length(total_lineages))
  names(lineages_idx) <- names(total_lineages)
  for(type in names(lineages_idx)) {
    lineages_idx[[type]] <- which(Seurat_Obj@meta.data[,clonotype_types[type]] %in% total_lineages[[type]])
  }
  gmp_lineages_idx <- vector("list", length = length(total_gmp_lineages))
  names(gmp_lineages_idx) <- names(total_gmp_lineages)
  for(type in names(gmp_lineages_idx)) {
    gmp_lineages_idx[[type]] <- which(Seurat_Obj@meta.data[,clonotype_types[type]] %in% total_gmp_lineages[[type]])
  }
  
  ### fill out the table
  for(i in 1:length(table_colnames)) {
    for(j in 1:length(unique(Seurat_Obj@meta.data$px))) {
      
      ### fill out for each column and for each patient
      px <- unique(Seurat_Obj@meta.data$px)[j]
      result_table[paste(px, "CD4+/CAR+ Cell #"),
                   table_colnames[i]] <- length(intersect(which(Seurat_Obj@meta.data$px == px),
                                                          cd4_carpos_idx))
      result_table[paste(px, "CD4+/CAR- Cell #"),
                   table_colnames[i]] <- length(intersect(which(Seurat_Obj@meta.data$px == px),
                                                          cd4_carneg_idx))
      result_table[paste(px, "CD8+/CAR+ Cell #"),
                   table_colnames[i]] <- length(intersect(which(Seurat_Obj@meta.data$px == px),
                                                          cd8_carpos_idx))
      result_table[paste(px, "CD8+/CAR- Cell #"),
                   table_colnames[i]] <- length(intersect(which(Seurat_Obj@meta.data$px == px),
                                                          cd8_carneg_idx))
      result_table[paste(px, "CD4+/CAR+ Clonotype #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                      cd4_carpos_idx),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD4+/CAR- Clonotype #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                      cd4_carneg_idx),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR+ Clonotype #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                      cd8_carpos_idx),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR- Clonotype #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                      cd8_carneg_idx),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD4+/CAR+ Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd4_carpos_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD4+/CAR- Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd4_carneg_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR+ Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd8_carpos_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR- Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd8_carneg_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD4+/CAR+ GMP Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(gmp_lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd4_carpos_idx)),
                                                                            clonotype_types[i]]))
      result_table[paste(px, "CD8+/CAR+ GMP Lineage #"),
                   table_colnames[i]] <- length(unique(Seurat_Obj@meta.data[intersect(gmp_lineages_idx[[i]],
                                                                                      intersect(which(Seurat_Obj@meta.data$px == px),
                                                                                                cd8_carpos_idx)),
                                                                            clonotype_types[i]]))
      
    }
  }
  
  ### save the result table
  write.xlsx2(result_table,
              sheetName = "Lineage_Statistics_Table",
              file = paste0(outputDir2, "Lineage_Statistics_Table_Per_Px.xlsx"))
  
  ### make a px-combined table
  specific_pxs <- c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                    "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                    "SJCAR19-09", "SJCAR19-10", "SJCAR19-11")
  specifix_pxs_idx <- as.vector(sapply(specific_pxs, function(x) grep(x, rownames(result_table), fixed = TRUE)))
  total_result_table <- matrix(0, length(table_rownames), length(table_colnames))
  rownames(total_result_table) <- table_rownames
  colnames(total_result_table) <- table_colnames
  for(r in table_rownames) {
    for(c in table_colnames) {
      total_result_table[r,c] <- sum(result_table[intersect(specifix_pxs_idx,
                                                            grep(r, rownames(result_table), fixed = TRUE)),
                                                  c])
    }
  }
  
  ### save the total result table
  write.xlsx2(total_result_table,
              sheetName = "Total_Lineage_Statistics_Table",
              file = paste0(outputDir2, "Lineage_Statistics_Table_Total.xlsx"))
  
  
  #
  ### 3. CAR+ % for each library - proportional bar graph & table - CAR+, CAR-, CD4+, CD8+
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/3/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### some pre-calculated indices
  cd4_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4"),
                              which(Seurat_Obj@meta.data$CAR == "CARpos"))
  cd4_carneg_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4"),
                              which(Seurat_Obj@meta.data$CAR == "CARneg"))
  cd8_carpos_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"),
                              which(Seurat_Obj@meta.data$CAR == "CARpos"))
  cd8_carneg_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"),
                              which(Seurat_Obj@meta.data$CAR == "CARneg"))
  
  ### examine the proportion of cells in each patient - barplot 
  plot_df <- data.frame(LIB=as.vector(sapply(unique(Seurat_Obj@meta.data$library), function(x) rep(x,4))),
                        TYPE=rep(c("CD4+ CAR+", "CD4+ CAR-", "CD8+ CAR+", "CD8+ CAR-"),
                                length(unique(Seurat_Obj@meta.data$library))),
                        NUM=0,
                        PCNT=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:nrow(plot_df)) {
    type <- strsplit(plot_df$TYPE[i], split = " ", fixed = TRUE)[[1]]
    if(type[1] == "CD4+" && type[2] == "CAR+") {
      target_idx <- cd4_carpos_idx
    } else if(type[1] == "CD4+" && type[2] == "CAR-") {
      target_idx <- cd4_carneg_idx
    } else if(type[1] == "CD8+" && type[2] == "CAR+") {
      target_idx <- cd8_carpos_idx
    } else {
      target_idx <- cd8_carneg_idx
    }
    plot_df$NUM[i] <- length(intersect(which(Seurat_Obj@meta.data$library == plot_df$LIB[i]),
                                       target_idx))
    plot_df$PCNT[i] <- round(plot_df$NUM[i]*100/length(which(Seurat_Obj@meta.data$library == plot_df$LIB[i])), digits = 0)
  }
  ### pcnt < 5 -> ""
  plot_df$PCNT[which(as.numeric(plot_df$PCNT) < 5)] <- ""
  plot_df$PCNT <- as.character(plot_df$PCNT)
  ggplot(data=plot_df, aes_string(x="LIB", y="NUM", fill="TYPE", label="PCNT")) +
    geom_bar(position = "stack", stat = "identity") +
    ggtitle("Cell %") +
    geom_text(size = 3, position = position_stack(vjust = 1)) +
    coord_flip() +
    scale_y_continuous(expand = c(0,300)) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 10),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir2, "Cell_Distribution_In_Each_Library.png"), width = 20, height = 10, dpi = 400)
  
  
  #
  ### 4. Clone size between CAR+ lineages vs non-lineage CAR+ after infusion
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/4/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### AFTER-INFUSION - GMP SUBSISTER CAR+ VS THE REST CAR+
  
  ### get AI indicies
  ai_indicies <- which(Seurat_Obj@meta.data$time2 %in% c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo"))
  
  ### target lineages
  subsister_clones_ai <- unique(intersect(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")],
                                          Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[ai_indicies]))
  rest_carpos_clones_ai <- unique(intersect(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[intersect(Seurat_Obj@meta.data$CAR == "CARpos",
                                                                                                                  which(is.na(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister)))],
                                            Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[ai_indicies]))
  
  ### get persister vs rest carpos clone sizes in AI time points
  subsister_clone_sizes_ai <- rep(0, length(subsister_clones_ai))
  names(subsister_clone_sizes_ai) <- subsister_clones_ai
  for(clone in subsister_clones_ai) {
    subsister_clone_sizes_ai[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[ai_indicies] == clone))
  }
  
  rest_carpos_clone_sizes_ai <- rep(0, length(rest_carpos_clones_ai))
  names(rest_carpos_clone_sizes_ai) <- rest_carpos_clones_ai
  for(clone in rest_carpos_clones_ai) {
    rest_carpos_clone_sizes_ai[clone] <- length(which(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[ai_indicies] == clone))
  }
  
  ### prepare a data frame for the plot
  plot_df <- data.frame(Clone_Size=c(subsister_clone_sizes_ai,
                                     rest_carpos_clone_sizes_ai),
                        Group=c(rep("Subsisters_After_Infusion", length(subsister_clone_sizes_ai)),
                                rep("Rest_CARpos_After_Infusion", length(rest_carpos_clone_sizes_ai))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df$Group <- factor(plot_df$Group, levels = c("Subsisters_After_Infusion", "Rest_CARpos_After_Infusion"))
  
  ### draw a violin plot
  ggplot(plot_df, aes_string(x="Group", y="Clone_Size", fill = "Group")) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill="white") +
    ylim(c(-0.5, 2)) +
    geom_text(data = data.frame(Group=c("Subsisters_After_Infusion", "Rest_CARpos_After_Infusion"),
                                Median=c(median(subsister_clone_sizes_ai),
                                         median(rest_carpos_clone_sizes_ai)),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Median"),
              size = 10, hjust = -2, vjust = -2) +
    geom_text(data = data.frame(Group=c("Subsisters_After_Infusion", "Rest_CARpos_After_Infusion"),
                                Median=c(median(subsister_clone_sizes_ai),
                                         median(rest_carpos_clone_sizes_ai)),
                                Length=c(paste("n =", length(subsister_clone_sizes_ai)),
                                         paste("n =", length(rest_carpos_clone_sizes_ai))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "Group", y = "Median", label = "Length"),
              size = 5, hjust = 0.5, vjust = 35) +
    labs(title="", x="", y = "Clone Size", fill = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_classic(base_size = 36) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 30))
  ggsave(paste0(outputDir, "/", "Violin_After_Infusion_Clone_Size_S_vs_R.png"), width = 20, height = 12, dpi = 500)
  
  ### density plot
  ggplot(plot_df, aes_string(x="Clone_Size", col="Group")) +
    geom_density(size = 2) +
    xlim(c(0,20)) +
    geom_text(data = data.frame(Group=c("Subsisters_After_Infusion", "Rest_CARpos_After_Infusion"),
                                x = c(18, 18),
                                y = c(0.3, 0.2),
                                Length=c(paste("n =", length(subsister_clone_sizes_ai)),
                                         paste("n =", length(rest_carpos_clone_sizes_ai))),
                                stringsAsFactors = FALSE, check.names = FALSE),
              aes_string(x = "x", y = "y", label = "Length"),
              size = 5, hjust = 0.5, vjust = 0.5) +
    labs(title="", x="Clone Size", y = "Density", col = "") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_classic(base_size = 36) +
    theme(legend.text = element_text(size = 35))
  ggsave(paste0(outputDir, "/", "Density_After_Infusion_Clone_Size_S_vs_R.png"), width = 15, height = 12, dpi = 500)
  
  
  #
  ### 5. a) PCA/UMAP plot of CAR+s over time (coloring based on time)
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/5/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get CARpos-only seurat object
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$CAR)
  sub_seurat_obj <- subset(Seurat_Obj, idents = c("CARpos"))
  
  ### after gmp time points only
  after_gmp_time_points <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk6",
                             "Wk8", "3mo", "6mo", "9mo")
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$time2)
  sub_seurat_obj <- subset(sub_seurat_obj, idents = intersect(after_gmp_time_points,
                                                              unique(sub_seurat_obj@meta.data$time2)))
  
  #
  ### run pca on the sub_seurat_obj
  #
  ### normalization
  sub_seurat_obj <- NormalizeData(sub_seurat_obj,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  ### find variable genes
  sub_seurat_obj <- FindVariableFeatures(sub_seurat_obj,
                                         selection.method = "vst", nfeatures = 2000)
  ### scaling
  sub_seurat_obj <- ScaleData(sub_seurat_obj,
                              vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  ### PCA
  sub_seurat_obj <- RunPCA(sub_seurat_obj,
                           features = VariableFeatures(object = sub_seurat_obj),
                           npcs = 15)
  ### UMAP
  sub_seurat_obj <- RunUMAP(sub_seurat_obj, dims = 1:15)
  
  ### get seurat object for some specific patients
  sub_seurat_obj <- SetIdent(object = sub_seurat_obj,
                             cells = rownames(sub_seurat_obj@meta.data),
                             value = sub_seurat_obj@meta.data$px)
  sub_seurat_obj2 <- subset(sub_seurat_obj, idents = c("SJCAR19-02", "SJCAR19-04", "SJCAR19-05",
                                                       "SJCAR19-06", "SJCAR19-07", "SJCAR19-08",
                                                       "SJCAR19-09", "SJCAR19-10", "SJCAR19-11"))
  
  ### factorize the time2 column
  sub_seurat_obj2@meta.data$time2 <- factor(sub_seurat_obj2@meta.data$time2,
                                            levels = intersect(total_time_points, sub_seurat_obj2@meta.data$time2))
  
  ### PCA with time by each patient
  p <- DimPlot(object = sub_seurat_obj2, reduction = "pca",
               group.by = "time2", split.by = "px",
               pt.size = 3, ncol = 3) +
    ggtitle("") +
    labs(color="Time") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30)) +
    scale_colour_brewer(palette = "OrRd")
  # p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "PCA_CARpos_Time.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### PCA with time with all patients
  p <- DimPlot(object = sub_seurat_obj2, reduction = "pca",
               group.by = "time2",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Time") +
    scale_colour_brewer(palette = "OrRd") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "PCA_CARpos_Time_ALL.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  ### UMAP with time by each patient
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "time2", split.by = "px",
               pt.size = 3, ncol = 3) +
    ggtitle("") +
    labs(color="Time") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30)) +
    scale_colour_brewer(palette = "OrRd")
  # p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir2, "UMAP_CARpos_Time.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### UMAP with time with all patients
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "time2",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Time") +
    scale_colour_brewer(palette = "OrRd") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "UMAP_CARpos_Time_ALL.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  
  
  
  
  
  
  
  
  
  
  
    
}
