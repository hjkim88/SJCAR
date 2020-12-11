###
#   File name : CD4_CD8_Lineage_Investigation.R
#   Author    : Hyunjin Kim
#   Date      : Nov 8, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Annotate CD4/CD8 to the cells and identify CD4/CD8 - associated lineages.
#
#   Instruction
#               1. Source("CD4_CD8_Lineage_Investigation.R")
#               2. Run the function "cd4_cd8_investigation" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_CD4_CD8_Lineage_Investigation.R/CD4_CD8_Lineage_Investigation.R")
#               > cd4_cd8_investigation(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                       clonotype_lineage_info_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Clonotype_Lineages.RDS",
#                                       outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results2/")
###

cd4_cd8_investigation <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
                                  clonotype_lineage_info_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Clonotype_Lineages.RDS",
                                  outputDir="./results/New2/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(SingleR, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("SingleR")
    require(SingleR, quietly = TRUE)
  }
  if(!require(SingleCellExperiment, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("SingleCellExperiment")
    require(SingleCellExperiment, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
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
  
  ### CD4/CD8 by expression
  exp_quotient <- Seurat_Obj@assays$RNA@counts["CD4",] / Seurat_Obj@assays$RNA@counts["CD8A",]
  cd4_idx <- which(exp_quotient > 1)
  cd8_idx <- which(exp_quotient < 1)
  both_idx <- which(exp_quotient == 1)
  Seurat_Obj@meta.data$CD4_CD8_by_Exp <- NA
  Seurat_Obj@meta.data$CD4_CD8_by_Exp[cd4_idx] <- "CD4"
  Seurat_Obj@meta.data$CD4_CD8_by_Exp[cd8_idx] <- "CD8"
  Seurat_Obj@meta.data$CD4_CD8_by_Exp[both_idx] <- "BOTH"
  rm(exp_quotient)
  gc()
  
  ### download the hematopoietic cell population data
  NHD <- NovershternHematopoieticData()
  
  ### set idents with the libary
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$library)
  
  ### for each library add annotations of CD4/CD8 using singleR
  Seurat_Obj@meta.data$CD4_CD8_by_NHD_main <- NA
  Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine <- NA
  for(lib in unique(Seurat_Obj@meta.data$library)) {
    ### print progress
    writeLines(paste(lib))
    
    ### get gene expressions from the subset of the Seurat object
    target_mat <- as.SingleCellExperiment(subset(Seurat_Obj, idents = lib), assay = "RNA")
    
    ### get CD4/CD8 annotation using the main labels
    NHD.main <- SingleR(test = target_mat, ref = list(NHD), labels = list(NHD$label.main))$labels
    
    ### get CD4/CD8 annotation using the fine labels
    NHD.fine <- SingleR(test = target_mat, ref = list(NHD), labels = list(NHD$label.fine))$labels
    
    ### add annotions to the meta data
    Seurat_Obj@meta.data[colnames(target_mat),"CD4_CD8_by_NHD_main"] <- NHD.main
    Seurat_Obj@meta.data[colnames(target_mat),"CD4_CD8_by_NHD_fine"] <- NHD.fine
    
    gc()
  }
  
  ### change NA to "NA" in CD4/CD8 columns
  Seurat_Obj@meta.data$CD4_CD8_by_Exp[which(is.na(Seurat_Obj@meta.data$CD4_CD8_by_Exp))] <- "NA"
  Seurat_Obj@meta.data$CD4_CD8_by_NHD_main[which(is.na(Seurat_Obj@meta.data$CD4_CD8_by_NHD_main))] <- "NA"
  Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine[which(is.na(Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine))] <- "NA"
  
  ### CD4/CD8 Consensus
  Seurat_Obj@meta.data$CD4_CD8_by_Consensus <- NA
  for(i in 1:nrow(Seurat_Obj@meta.data)) {
    cd4_cnt <- 0
    cd8_cnt <- 0
    if(Seurat_Obj@meta.data$CD4_CD8_by_Exp[i] == "CD4") {
      cd4_cnt <- cd4_cnt + 1
    } else if (Seurat_Obj@meta.data$CD4_CD8_by_Exp[i] == "CD8") {
      cd8_cnt <- cd8_cnt + 1
    }
    if(grepl("CD4", Seurat_Obj@meta.data$CD4_CD8_by_NHD_main[i], fixed = TRUE)) {
      cd4_cnt <- cd4_cnt + 1
    } else if(grepl("CD8", Seurat_Obj@meta.data$CD4_CD8_by_NHD_main[i], fixed = TRUE)) {
      cd8_cnt <- cd8_cnt + 1
    }
    if(grepl("CD4", Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine[i], fixed = TRUE)) {
      cd4_cnt <- cd4_cnt + 1
    } else if(grepl("CD8", Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine[i], fixed = TRUE)) {
      cd8_cnt <- cd8_cnt + 1
    }
    
    ### decide with the consensus
    if(cd4_cnt > cd8_cnt) {
      Seurat_Obj@meta.data$CD4_CD8_by_Consensus[i] <- "CD4"
    } else if(cd4_cnt < cd8_cnt) {
      Seurat_Obj@meta.data$CD4_CD8_by_Consensus[i] <- "CD8"
    } else if((cd4_cnt != 0) && (cd4_cnt == cd8_cnt)) {
      Seurat_Obj@meta.data$CD4_CD8_by_Consensus[i] <- "BOTH"
    }
    
    ### print progress
    if(i %% 10000 == 0) {
      writeLines(paste(i, "/", nrow(Seurat_Obj@meta.data)))
    }
  }
  
  ### change "NA" to NA in CD4/CD8 columns
  Seurat_Obj@meta.data$CD4_CD8_by_Exp[which(Seurat_Obj@meta.data$CD4_CD8_by_Exp == "NA")] <- NA
  Seurat_Obj@meta.data$CD4_CD8_by_NHD_main[which(Seurat_Obj@meta.data$CD4_CD8_by_NHD_main == "NA")] <- NA
  Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine[which(Seurat_Obj@meta.data$CD4_CD8_by_NHD_fine == "NA")] <- NA
  
  ### load clonotype lineages info
  SJCAR19_Clonotype_Frequency <- readRDS(clonotype_lineage_info_path)
  
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
  
  ### draw an alluvial plot for each patient
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    
    ### print progress
    writeLines(paste(patient))
    
    ### create a dir
    outputDir2 <- paste0(outputDir, patient, "/CD4_CD8/")
    dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
    
    ### get clonotype frequency data of the patient - CARpos only
    target_file <- SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[patient]]
    
    ### remove all zero time points
    time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
    
    ### get time points except the Total
    time_points <- setdiff(time_points, c("Total"))
    
    if(length(time_points) > 1) {
      ###  get lineages
      lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
        return(length(which(x > 0)) > 1)  
      })),time_points]
      
      if(nrow(lineage_table) > 0) {
        ### separate the CD4 and CD8 cells
        lineage_table <- data.frame(clone_id=rownames(lineage_table),
                                    lineage_table,
                                    stringsAsFactors = FALSE, check.names = FALSE)
        lineage_table_CD <- data.frame(sapply(lineage_table, function(x) c(rbind(x, x, x))),
                                              stringsAsFactors = FALSE, check.names = FALSE)
        lineage_table_CD <- data.frame(clone_id=lineage_table_CD$clone_id,
                                       cell_type=rep(c("CD4", "CD8", "ALL"), nrow(lineage_table)),
                                       sapply(lineage_table_CD[,time_points],
                                              as.numeric),
                                       stringsAsFactors = FALSE, check.names = FALSE)
        rownames(lineage_table_CD) <- paste(c(rbind(rownames(lineage_table),
                                                    rownames(lineage_table),
                                                    rownames(lineage_table))),
                                            c("CD4", "CD8", "ALL"), sep = "_")
        
        ### specific time point indicies
        tp_indicies <- lapply(time_points, function(x) which(Seurat_Obj@meta.data$time == x))
        names(tp_indicies) <- time_points
        
        ### fill out the new table
        lineage_table_CD$total_count <- 0
        for(i in 1:nrow(lineage_table_CD)) {
          if(lineage_table_CD$cell_type[i] != "ALL") {
            car_clone_idx <- intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                       intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient == lineage_table_CD$clone_id[i]),
                                                 which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == lineage_table_CD$cell_type[i])))
            for(tp in time_points) {
              lineage_table_CD[i,tp] <- length(intersect(car_clone_idx, tp_indicies[[tp]]))
            }
          }
          lineage_table_CD$total_count[i] <- sum(lineage_table_CD[i,time_points])
        }
        
        ### save the table as Excel file
        write.xlsx2(lineage_table_CD, file = paste0(outputDir2, patient, "_CARpos_CD4_CD8_Lineage_Table.xlsx"),
                    sheetName = "CARpos_CD4_CD8_Lineage_Table", row.names = FALSE)
        
        ### get an input data frame for the alluvial plot
        total_rows <- length(which(lineage_table[,time_points] > 0))
        if(total_rows > 0) {
          plot_df <- data.frame(Time=rep("", total_rows),
                                Clone_Size=rep(0, total_rows),
                                Clone=rep("", total_rows),
                                CD4=rep("", total_rows))
          cnt <- 1
          for(i in 1:nrow(lineage_table)) {
            for(tp in time_points) {
              if(lineage_table[i,tp] > 0) {
                if(lineage_table_CD[paste0(lineage_table$clone_id[i], "_CD4"),tp] > 0) {
                  plot_df[cnt,] <- c(tp,
                                     lineage_table[i,tp],
                                     lineage_table$clone_id[i],
                                     paste0("CD4:", lineage_table_CD[paste0(lineage_table$clone_id[i], "_CD4"),tp]))
                } else {
                  plot_df[cnt,] <- c(tp,
                                     lineage_table[i,tp],
                                     lineage_table$clone_id[i],
                                     "")
                }
                cnt <- cnt + 1
              }
            }
          }
          plot_df$Time <- factor(plot_df$Time, levels = time_points)
          
          ### numerize the clone_size column
          plot_df$Clone_Size <- as.numeric(plot_df$Clone_Size)
          
          ### draw the alluvial plot (GMP - GMP-redo)
          ggplot(plot_df,
                 aes(x = Time, stratum = Clone, alluvium = Clone,
                     y = Clone_Size,
                     fill = Clone, label = CD4)) +
            ggtitle(paste("Clonal Tracing in the CAR+ cells of", patient)) +
            geom_flow() +
            geom_stratum(alpha = 1) +
            geom_text(stat = "stratum", size = 2, col = "red") +
            rotate_x_text(90) +
            theme_pubr(legend = "none") +
            theme(axis.title.x = element_blank()) +
            theme_cleveland2() +
            scale_fill_viridis(discrete = T) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          ggsave(file = paste0(outputDir2, "CARpos_Clonal_Tracing_", patient, ".png"), width = 20, height = 10, dpi = 300)
          
          ### 1-on-1 alluvial plot with GMP
          if(length(grep("^GMP$", time_points)) == 1) {
            p <- vector("list", length = length(time_points)-1)
            names(p) <- time_points[-grep("^GMP$", time_points)]
            for(tp in names(p)) {
              p[[tp]] <- ggplot(plot_df[union(which(plot_df$Time == "GMP"),
                                              which(plot_df$Time == tp)),],
                                aes(x = Time, stratum = Clone, alluvium = Clone,
                                    y = Clone_Size,
                                    fill = Clone, label = CD4)) +
                geom_flow() +
                geom_stratum(alpha = 1) +
                geom_text(stat = "stratum", size = 2, col = "red") +
                rotate_x_text(90) +
                theme_pubr(legend = "none") +
                theme(axis.title.x = element_blank()) +
                theme_cleveland2() +
                scale_fill_viridis(discrete = T) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
            }
            
            ### arrange the plots and save
            fName <- paste0("CARpos_Clonal_Tracing_", patient, "_1-on-1")
            if(length(p) == 1) {
              rowNum <- 1
              colNum <- 1
            } else if(length(p) == 2) {
              rowNum <- 1
              colNum <- 2
            } else if(length(p) < 5) {
              rowNum <- 2
              colNum <- 2
            } else if(length(p) < 7) {
              rowNum <- 3
              colNum <- 2
            } else {
              rowNum <- ceiling(sqrt(length(p)))
              colNum <- ceiling(sqrt(length(p)))
            }
            g <- arrangeGrob(grobs = p,
                             nrow = rowNum,
                             ncol = colNum,
                             top = fName)
            ggsave(file = paste0(outputDir2, fName, ".png"), g, width = 20, height = 12, dpi = 300)
          }
          
          ### 1-on-1 alluvial plot with GMP-redo
          if(length(grep("^GMP-redo$", time_points)) == 1) {
            p <- vector("list", length = length(time_points)-1)
            names(p) <- time_points[-grep("^GMP-redo$", time_points)]
            for(tp in names(p)) {
              p[[tp]] <- ggplot(plot_df[union(which(plot_df$Time == "GMP-redo"),
                                              which(plot_df$Time == tp)),],
                                aes(x = Time, stratum = Clone, alluvium = Clone,
                                    y = Clone_Size,
                                    fill = Clone, label = CD4)) +
                geom_flow() +
                geom_stratum(alpha = 1) +
                geom_text(stat = "stratum", size = 2, col = "red") +
                rotate_x_text(90) +
                theme_pubr(legend = "none") +
                theme(axis.title.x = element_blank()) +
                theme_cleveland2() +
                scale_fill_viridis(discrete = T) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
            }
            
            ### arrange the plots and save
            fName <- paste0("CARpos_Clonal_Tracing_", patient, "_1-on-1(2)")
            if(length(p) == 1) {
              rowNum <- 1
              colNum <- 1
            } else if(length(p) == 2) {
              rowNum <- 1
              colNum <- 2
            } else if(length(p) < 5) {
              rowNum <- 2
              colNum <- 2
            } else if(length(p) < 7) {
              rowNum <- 3
              colNum <- 2
            } else {
              rowNum <- ceiling(sqrt(length(p)))
              colNum <- ceiling(sqrt(length(p)))
            }
            g <- arrangeGrob(grobs = p,
                             nrow = rowNum,
                             ncol = colNum,
                             top = fName)
            ggsave(file = paste0(outputDir2, fName, ".png"), g, width = 20, height = 12, dpi = 300)
          }
          
          ### draw the alluvial plot (GMP-redo - GMP)
          if(length(grep("GMP", time_points)) == 2) {
            gmp_idx <- which(time_points == "GMP")
            gmp_redo_idx <- which(time_points == "GMP-redo")
            time_points[gmp_idx] <- "GMP-redo"
            time_points[gmp_redo_idx] <- "GMP"
            plot_df$Time <- factor(plot_df$Time, levels = time_points)
            ggplot(plot_df,
                   aes(x = Time, stratum = Clone, alluvium = Clone,
                       y = Clone_Size,
                       fill = Clone, label = CD4)) +
              ggtitle(paste("Clonal Tracing in the CAR+ cells of", patient)) +
              geom_flow() +
              geom_stratum(alpha = 1) +
              geom_text(stat = "stratum", size = 2, col = "red") +
              rotate_x_text(90) +
              theme_pubr(legend = "none") +
              theme(axis.title.x = element_blank()) +
              theme_cleveland2() +
              scale_fill_viridis(discrete = T) +
              scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
            ggsave(file = paste0(outputDir2, "CARpos_Clonal_Tracing_", patient, "(2).png"), width = 20, height = 10, dpi = 300)
          }
        }
        
        gc()
      }
    }
    
    
    ### get clonotype frequency data of the patient - ALL cells
    target_file <- SJCAR19_Clonotype_Frequency[["ALL"]][[patient]]
    
    ### remove all zero time points
    time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
    
    ### get time points except the Total
    time_points <- setdiff(time_points, c("Total"))
    
    if(length(time_points) > 1) {
      ###  get lineages
      lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
        return(length(which(x > 0)) > 1)  
      })),time_points]
      
      ### separate the CD4 and CD8 cells
      lineage_table <- data.frame(clone_id=rownames(lineage_table),
                                  lineage_table,
                                  stringsAsFactors = FALSE, check.names = FALSE)
      lineage_table_CD <- data.frame(sapply(lineage_table, function(x) c(rbind(x, x, x))),
                                     stringsAsFactors = FALSE, check.names = FALSE)
      lineage_table_CD <- data.frame(clone_id=lineage_table_CD$clone_id,
                                     cell_type=rep(c("CD4", "CD8", "ALL"), nrow(lineage_table)),
                                     sapply(lineage_table_CD[,time_points],
                                            as.numeric),
                                     stringsAsFactors = FALSE, check.names = FALSE)
      rownames(lineage_table_CD) <- paste(c(rbind(rownames(lineage_table),
                                                  rownames(lineage_table),
                                                  rownames(lineage_table))),
                                          c("CD4", "CD8", "ALL"), sep = "_")
      
      ### specific time point indicies
      tp_indicies <- lapply(time_points, function(x) which(Seurat_Obj@meta.data$time == x))
      names(tp_indicies) <- time_points
      
      ### fill out the new table
      lineage_table_CD$total_count <- 0
      for(i in 1:nrow(lineage_table_CD)) {
        if(lineage_table_CD$cell_type[i] != "ALL") {
          clone_idx <- intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient == lineage_table_CD$clone_id[i]),
                                 which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == lineage_table_CD$cell_type[i]))
          for(tp in time_points) {
            lineage_table_CD[i,tp] <- length(intersect(clone_idx, tp_indicies[[tp]]))
          }
        }
        lineage_table_CD$total_count[i] <- sum(lineage_table_CD[i,time_points])
      }
      
      ### save the table as Excel file
      write.xlsx2(lineage_table_CD, file = paste0(outputDir2, patient, "_ALL_CD4_CD8_Lineage_Table.xlsx"),
                  sheetName = "ALL_CD4_CD8_Lineage_Table", row.names = FALSE)
      
      ### get an input data frame for the alluvial plot
      total_rows <- length(which(lineage_table[,time_points] > 0))
      if(total_rows > 0) {
        plot_df <- data.frame(Time=rep("", total_rows),
                              Clone_Size=rep(0, total_rows),
                              Clone=rep("", total_rows),
                              CD4=rep("", total_rows))
        cnt <- 1
        for(i in 1:nrow(lineage_table)) {
          for(tp in time_points) {
            if(lineage_table[i,tp] > 0) {
              if(lineage_table_CD[paste0(lineage_table$clone_id[i], "_CD4"),tp] > 0) {
                plot_df[cnt,] <- c(tp,
                                   lineage_table[i,tp],
                                   lineage_table$clone_id[i],
                                   paste0("CD4:", lineage_table_CD[paste0(lineage_table$clone_id[i], "_CD4"),tp]))
              } else {
                plot_df[cnt,] <- c(tp,
                                   lineage_table[i,tp],
                                   lineage_table$clone_id[i],
                                   "")
              }
              cnt <- cnt + 1
            }
          }
        }
        plot_df$Time <- factor(plot_df$Time, levels = time_points)
        
        ### numerize the clone_size column
        plot_df$Clone_Size <- as.numeric(plot_df$Clone_Size)
        
        ### draw the alluvial plot (GMP - GMP-redo)
        ggplot(plot_df,
               aes(x = Time, stratum = Clone, alluvium = Clone,
                   y = Clone_Size,
                   fill = Clone, label = CD4)) +
          ggtitle(paste("Clonal Tracing in ALL cells of", patient)) +
          geom_flow() +
          geom_stratum(alpha = 1) +
          geom_text(stat = "stratum", size = 2, col = "red") +
          rotate_x_text(90) +
          theme_pubr(legend = "none") +
          theme(axis.title.x = element_blank()) +
          theme_cleveland2() +
          scale_fill_viridis(discrete = T) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
        ggsave(file = paste0(outputDir2, "ALL_Clonal_Tracing_", patient, ".png"), width = 20, height = 10, dpi = 300)
        
        ### 1-on-1 alluvial plot with GMP
        if(length(grep("^GMP$", time_points)) == 1) {
          p <- vector("list", length = length(time_points)-1)
          names(p) <- time_points[-grep("^GMP$", time_points)]
          for(tp in names(p)) {
            p[[tp]] <- ggplot(plot_df[union(which(plot_df$Time == "GMP"),
                                            which(plot_df$Time == tp)),],
                              aes(x = Time, stratum = Clone, alluvium = Clone,
                                  y = Clone_Size,
                                  fill = Clone, label = CD4)) +
              geom_flow() +
              geom_stratum(alpha = 1) +
              geom_text(stat = "stratum", size = 2, col = "red") +
              rotate_x_text(90) +
              theme_pubr(legend = "none") +
              theme(axis.title.x = element_blank()) +
              theme_cleveland2() +
              scale_fill_viridis(discrete = T) +
              scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          }
          
          ### arrange the plots and save
          fName <- paste0("ALL_Clonal_Tracing_", patient, "_1-on-1")
          if(length(p) == 1) {
            rowNum <- 1
            colNum <- 1
          } else if(length(p) == 2) {
            rowNum <- 1
            colNum <- 2
          } else if(length(p) < 5) {
            rowNum <- 2
            colNum <- 2
          } else if(length(p) < 7) {
            rowNum <- 3
            colNum <- 2
          } else {
            rowNum <- ceiling(sqrt(length(p)))
            colNum <- ceiling(sqrt(length(p)))
          }
          g <- arrangeGrob(grobs = p,
                           nrow = rowNum,
                           ncol = colNum,
                           top = fName)
          ggsave(file = paste0(outputDir2, fName, ".png"), g, width = 20, height = 12, dpi = 300)
        }
        
        ### 1-on-1 alluvial plot with GMP-redo
        if(length(grep("^GMP-redo$", time_points)) == 1) {
          p <- vector("list", length = length(time_points)-1)
          names(p) <- time_points[-grep("^GMP-redo$", time_points)]
          for(tp in names(p)) {
            p[[tp]] <- ggplot(plot_df[union(which(plot_df$Time == "GMP-redo"),
                                            which(plot_df$Time == tp)),],
                              aes(x = Time, stratum = Clone, alluvium = Clone,
                                  y = Clone_Size,
                                  fill = Clone, label = CD4)) +
              geom_flow() +
              geom_stratum(alpha = 1) +
              geom_text(stat = "stratum", size = 2, col = "red") +
              rotate_x_text(90) +
              theme_pubr(legend = "none") +
              theme(axis.title.x = element_blank()) +
              theme_cleveland2() +
              scale_fill_viridis(discrete = T) +
              scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          }
          
          ### arrange the plots and save
          fName <- paste0("ALL_Clonal_Tracing_", patient, "_1-on-1(2)")
          if(length(p) == 1) {
            rowNum <- 1
            colNum <- 1
          } else if(length(p) == 2) {
            rowNum <- 1
            colNum <- 2
          } else if(length(p) < 5) {
            rowNum <- 2
            colNum <- 2
          } else if(length(p) < 7) {
            rowNum <- 3
            colNum <- 2
          } else {
            rowNum <- ceiling(sqrt(length(p)))
            colNum <- ceiling(sqrt(length(p)))
          }
          g <- arrangeGrob(grobs = p,
                           nrow = rowNum,
                           ncol = colNum,
                           top = fName)
          ggsave(file = paste0(outputDir2, fName, ".png"), g, width = 20, height = 12, dpi = 300)
        }
        
        ### draw the alluvial plot (GMP-redo - GMP)
        if(length(grep("GMP", time_points)) == 2) {
          gmp_idx <- which(time_points == "GMP")
          gmp_redo_idx <- which(time_points == "GMP-redo")
          time_points[gmp_idx] <- "GMP-redo"
          time_points[gmp_redo_idx] <- "GMP"
          plot_df$Time <- factor(plot_df$Time, levels = time_points)
          ggplot(plot_df,
                 aes(x = Time, stratum = Clone, alluvium = Clone,
                     y = Clone_Size,
                     fill = Clone, label = CD4)) +
            ggtitle(paste("Clonal Tracing in the CAR+ cells of", patient)) +
            geom_flow() +
            geom_stratum(alpha = 1) +
            geom_text(stat = "stratum", size = 2, col = "red") +
            rotate_x_text(90) +
            theme_pubr(legend = "none") +
            theme(axis.title.x = element_blank()) +
            theme_cleveland2() +
            scale_fill_viridis(discrete = T) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
          ggsave(file = paste0(outputDir2, "ALL_Clonal_Tracing_", patient, "(2).png"), width = 20, height = 10, dpi = 300)
        }
      }
      
      gc()
    }
    
  }
  
  
  ### check if there are lineages that have mixed CD4/CD8 cells
  ### that means something is wrong with CD4 or CD8 cells because they should all be the same
  
  ### load clonotype lineages info
  SJCAR19_Clonotype_Frequency <- readRDS(clonotype_lineage_info_path)
  
  ### for each patient produce how many lineages (clones) have the mixed cells
  for(patient in unique(Seurat_Obj@meta.data$px)) {
    
    ### print progress
    writeLines(paste(patient))
    
    ### create a dir
    outputDir2 <- paste0(outputDir, patient, "/CD4_CD8/")
    dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
    
    ### get clonotype frequency data of the patient - CARpos only
    target_file <- SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[patient]]
    
    ### remove all zero time points
    time_points <- colnames(target_file)[which(apply(target_file, 2, sum) != 0)]
    
    ### get time points except the Total
    time_points <- setdiff(time_points, c("Total"))
    
    if(length(time_points) > 1) {
      ###  get lineages
      lineage_table <- target_file[which(apply(target_file[,time_points], 1, function(x) {
        return(length(which(x > 0)) > 1)  
      })),time_points]
      
      if(nrow(lineage_table) > 0) {
        ### separate the CD4 and CD8 cells
        lineage_table <- data.frame(clone_id=rownames(lineage_table),
                                    lineage_table,
                                    stringsAsFactors = FALSE, check.names = FALSE)
        lineage_table_CD <- data.frame(sapply(lineage_table, function(x) c(rbind(x, x, x))),
                                       stringsAsFactors = FALSE, check.names = FALSE)
        lineage_table_CD <- data.frame(clone_id=lineage_table_CD$clone_id,
                                       cell_type=rep(c("CD4", "CD8", "ALL"), nrow(lineage_table)),
                                       sapply(lineage_table_CD[,time_points],
                                              as.numeric),
                                       stringsAsFactors = FALSE, check.names = FALSE)
        rownames(lineage_table_CD) <- paste(c(rbind(rownames(lineage_table),
                                                    rownames(lineage_table),
                                                    rownames(lineage_table))),
                                            c("CD4", "CD8", "ALL"), sep = "_")
        
        ### specific time point indicies
        tp_indicies <- lapply(time_points, function(x) which(Seurat_Obj@meta.data$time == x))
        names(tp_indicies) <- time_points
        
        ### fill out the new table
        lineage_table_CD$total_count <- 0
        for(i in 1:nrow(lineage_table_CD)) {
          if(lineage_table_CD$cell_type[i] != "ALL") {
            car_clone_idx <- intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                       intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient == lineage_table_CD$clone_id[i]),
                                                 which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == lineage_table_CD$cell_type[i])))
            for(tp in time_points) {
              lineage_table_CD[i,tp] <- length(intersect(car_clone_idx, tp_indicies[[tp]]))
            }
          }
          lineage_table_CD$total_count[i] <- sum(lineage_table_CD[i,time_points])
        }
        
        ### check whether a clone has mixed cells
        total_cellNum <- 0
        mixed_cellNum <- 0
        total_cellNum_specific <- NULL
        mixed_cellNum_specific <- NULL
        mixed_clones <- NULL
        for(clone in unique(lineage_table_CD$clone_id)) {
          if((lineage_table_CD[paste0(clone, "_CD4"),"total_count"] > 0) && (lineage_table_CD[paste0(clone, "_CD8"),"total_count"] > 0)) {
            temp <- ifelse(lineage_table_CD[paste0(clone, "_CD4"),"total_count"] > lineage_table_CD[paste0(clone, "_CD8"),"total_count"],
                           as.numeric(lineage_table_CD[paste0(clone, "_CD8"),"total_count"]),
                           as.numeric(lineage_table_CD[paste0(clone, "_CD4"),"total_count"]))
            mixed_cellNum <- mixed_cellNum + temp
            mixed_cellNum_specific <- c(mixed_cellNum_specific, temp)
            mixed_clones <- c(mixed_clones, clone)
            total_cellNum_specific <- c(total_cellNum_specific, lineage_table_CD[paste0(clone, "_ALL"),"total_count"])
          }
          total_cellNum <- total_cellNum + lineage_table_CD[paste0(clone, "_ALL"),"total_count"]
        }
        
        ### save result if exists
        if(length(mixed_clones) > 0) {
          result_table <- data.frame(Clone_ID=mixed_clones,
                                     Mixed_CellNum=mixed_cellNum_specific,
                                     Total_CellNum=total_cellNum_specific)
          
          write.xlsx2(result_table, file = paste0(outputDir2, patient, "_CD4_CD8_Mixed_Clones_And_Cells.xlsx"),
                      sheetName = "CD4_CD8_Mixed_Clones_And_Cells", row.names = FALSE)
        }
        
        ### save total result
        writeLines(paste(patient, "Total cell #:", total_cellNum, ",",
                         "Mixed cell #:", mixed_cellNum,
                         "Error %:", round(100*mixed_cellNum/total_cellNum, digits = 3), "%"))
      }
    }
    
  }
  
  #
  ### CD4/CD8 ratio in CAR+ and ALL
  #
  
  ### make an empty data frame
  cd4_cd8_ratio <- data.frame(Library=unique(Seurat_Obj@meta.data$library),
                              CD4_CD8_ALL_Num="",
                              CD4_CD8_CARpos_Num="",
                              CD4_CD8_ALL_Ratio="",
                              CD4_CD8_CARpos_Ratio="",
                              stringsAsFactors = FALSE, check.names = FALSE)
  rownames(cd4_cd8_ratio) <- cd4_cd8_ratio$Library
  
  ### fill the data frame
  for(lib in cd4_cd8_ratio$Library) {
    
    ### prepare the numbers
    total_num_all <- length(which(Seurat_Obj@meta.data$library == lib))
    total_num_car <- length(intersect(which(Seurat_Obj@meta.data$library == lib),
                                      which(Seurat_Obj@meta.data$CAR == "CARpos")))
    cd4_num_all <- length(intersect(which(Seurat_Obj@meta.data$library == lib),
                                    which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4")))
    cd8_num_all <- length(intersect(which(Seurat_Obj@meta.data$library == lib),
                                    which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8")))
    cd4_num_car <- length(intersect(intersect(which(Seurat_Obj@meta.data$library == lib),
                                              which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD4")),
                                    which(Seurat_Obj@meta.data$CAR == "CARpos")))
    cd8_num_car <- length(intersect(intersect(which(Seurat_Obj@meta.data$library == lib),
                                              which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8")),
                                    which(Seurat_Obj@meta.data$CAR == "CARpos")))
    cd4_ratio_all <- round(100 * cd4_num_all / total_num_all, digits = 3)
    cd8_ratio_all <- round(100 * cd8_num_all / total_num_all, digits = 3)
    cd4_ratio_car <- round(100 * cd4_num_car / total_num_car, digits = 3)
    cd8_ratio_car <- round(100 * cd8_num_car / total_num_car, digits = 3)
    
    ### fill the data frame
    cd4_cd8_ratio[lib,"CD4_CD8_ALL_Num"] <- paste("CD4:", cd4_num_all, ",", "CD8:", cd8_num_all)
    cd4_cd8_ratio[lib,"CD4_CD8_CARpos_Num"] <- paste("CD4:", cd4_num_car, ",", "CD8:", cd8_num_car)
    cd4_cd8_ratio[lib,"CD4_CD8_ALL_Ratio"] <- paste("CD4:", cd4_ratio_all, "% ,", "CD8:", cd8_ratio_all, "%")
    cd4_cd8_ratio[lib,"CD4_CD8_CARpos_Ratio"] <- paste("CD4:", cd4_ratio_car, "% ,", "CD8:", cd8_ratio_car, "%")
    
  }
  
  ### save the result as an Excel file
  write.xlsx2(cd4_cd8_ratio, file = paste0(outputDir, "CD4_CD8_Ratio.xlsx"),
              sheetName = "CD4_CD8_Ratio", row.names = FALSE)
  
}
