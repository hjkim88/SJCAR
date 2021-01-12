###
#   File name : Epitope_Spreading.R
#   Author    : Hyunjin Kim
#   Date      : Dec 9, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Analyses about epitope spreading in SJCAR19 data.
#
#   Instruction
#               1. Source("Epitope_Spreading.R")
#               2. Run the function "epitope_spreading_investigation" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Epitope_Spreading.R/Epitope_Spreading.R")
#               > epitope_spreading_investigation(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                                 clonotype_lineage_info_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Clonotype_Lineages.RDS",
#                                                 vdjdb_file_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/VDJDB_Homo_Sapiens_121320.xlsx",
#                                                 outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Lineages_by_CAR/")
###

epitope_spreading_investigation <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
                                            clonotype_lineage_info_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Clonotype_Lineages.RDS",
                                            vdjdb_file_path="./data/VDJDB_Homo_Sapiens_121320.xlsx",
                                            outputDir="./results/New2/Lineages_by_CAR/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  options(java.parameters = "-Xmx10240m")
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(immunarch, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("immunarch")
    require(immunarch, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
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
  
  ### load clonotype lineages info
  SJCAR19_Clonotype_Frequency <- readRDS(clonotype_lineage_info_path)
  
  ### load VDJDB
  # vdjdb2 <- dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.6.4/private/vdjdb.slim.txt.gz", "vdjdb")
  vdjdb <- read.xlsx2(file = vdjdb_file_path,
                      sheetIndex = 1,
                      stringsAsFactors = FALSE, check.names = FALSE)
  
  ### get combined CDR3 (e.g., TRA:CSVFK) and remove duplicates
  vdjdb$Combined_CDR3 <- paste(vdjdb$Gene, vdjdb$CDR3, sep = ":")
  vdjdb <- vdjdb[!duplicated(vdjdb$Combined_CDR3),]
  
  ### epitope spreading finding in each patient
  SJCAR19_Lineages_by_CAR <- vector("list", length = length(SJCAR19_Clonotype_Frequency[["ALL"]]))
  names(SJCAR19_Lineages_by_CAR) <- names(SJCAR19_Clonotype_Frequency[["ALL"]])
  for(patient in names(SJCAR19_Clonotype_Frequency[["ALL"]])) {
    
    ### separate the CAR+ and CAR- cells
    SJCAR19_Lineages_by_CAR[[patient]] <- SJCAR19_Clonotype_Frequency[["ALL"]][[patient]]
    SJCAR19_Lineages_by_CAR[[patient]] <- data.frame(Clone_ID=rownames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]]),
                                                     V_Gene=sapply(rownames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]]), function(x) {
                                                       return(paste(unique(Seurat_Obj@meta.data$v_gene[which(Seurat_Obj@meta.data$clonotype_id_by_patient == x)]), collapse = "-"))
                                                     }),
                                                     J_Gene=sapply(rownames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]]), function(x) {
                                                       return(paste(unique(Seurat_Obj@meta.data$j_gene[which(Seurat_Obj@meta.data$clonotype_id_by_patient == x)]), collapse = "-"))
                                                     }),
                                                     C_Gene=sapply(rownames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]]), function(x) {
                                                       return(paste(unique(Seurat_Obj@meta.data$c_gene[which(Seurat_Obj@meta.data$clonotype_id_by_patient == x)]), collapse = "-"))
                                                     }),
                                                     CDR3_AA=sapply(rownames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]]), function(x) {
                                                       return(paste(unique(Seurat_Obj@meta.data$cdr3_aa[which(Seurat_Obj@meta.data$clonotype_id_by_patient == x)]), collapse = "-"))
                                                     }),
                                                     CDR3_NT=sapply(rownames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]]), function(x) {
                                                       return(paste(unique(Seurat_Obj@meta.data$cdr3_nt[which(Seurat_Obj@meta.data$clonotype_id_by_patient == x)]), collapse = "-"))
                                                     }),
                                                     sapply(SJCAR19_Lineages_by_CAR[[patient]], as.numeric),
                                                     stringsAsFactors = FALSE, check.names = FALSE)
    SJCAR19_Lineages_by_CAR[[patient]] <- data.frame(sapply(SJCAR19_Lineages_by_CAR[[patient]],
                                                            function(x) c(rbind(x, x, x))),
                                                     stringsAsFactors = FALSE, check.names = FALSE)
    SJCAR19_Lineages_by_CAR[[patient]] <- data.frame(Cell_Type=rep(c("CARneg", "CARpos", "ALL"), nrow(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])),
                                                     SJCAR19_Lineages_by_CAR[[patient]],
                                                     stringsAsFactors = FALSE, check.names = FALSE)
    rownames(SJCAR19_Lineages_by_CAR[[patient]]) <- paste(SJCAR19_Lineages_by_CAR[[patient]]$Clone_ID,
                                                          c("CARneg", "CARpos", "ALL"),
                                                          sep = "_")
        
    ### numerize the numeric columns
    SJCAR19_Lineages_by_CAR[[patient]][,colnames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])] <- sapply(SJCAR19_Lineages_by_CAR[[patient]][,colnames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])], as.numeric)
    
    ### get time points
    time_points <- colnames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])[1:(ncol(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])-1)]
    
    ### fill out the new table
    for(i in 1:nrow(SJCAR19_Lineages_by_CAR[[patient]])) {
      if(SJCAR19_Lineages_by_CAR[[patient]]$Cell_Type[i] != "ALL") {
        clone_idx <- intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient == SJCAR19_Lineages_by_CAR[[patient]]$Clone_ID[i]),
                               which(Seurat_Obj@meta.data$CAR == SJCAR19_Lineages_by_CAR[[patient]]$Cell_Type[i]))
        for(tp in time_points) {
          SJCAR19_Lineages_by_CAR[[patient]][i,tp] <- length(intersect(clone_idx,
                                                                       which(Seurat_Obj@meta.data$time == tp)))
        }
        SJCAR19_Lineages_by_CAR[[patient]]$Total[i] <- sum(SJCAR19_Lineages_by_CAR[[patient]][i,time_points])
      }
    }
    
  }
  
  ### save the result
  saveRDS(SJCAR19_Lineages_by_CAR, file = paste0(outputDir, "SJCAR19_Lineages_by_CAR.RDS"))
  for(patient in names(SJCAR19_Lineages_by_CAR)) {
    write.xlsx2(SJCAR19_Lineages_by_CAR[[patient]],
                file = paste0(outputDir, "SJCAR19_Lineages_by_CAR.xlsx"),
                sheetName = patient,
                row.names = FALSE,
                append = TRUE)
  }
  
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
  
  ### CARneg threshold
  ### Any lineages that have AFTER-GMP CAR- cells larger than the threshold will be considered interesting
  carneg_threshold <- 10
  
  ### find epitope spreading candidates
  interesting_clones <- vector("list", length = length(SJCAR19_Lineages_by_CAR))
  names(interesting_clones) <- names(SJCAR19_Lineages_by_CAR)
  for(patient in names(SJCAR19_Lineages_by_CAR)) {
    
    ### get index after GMP
    gmp_idx <- which(colnames(SJCAR19_Lineages_by_CAR[[patient]]) == "GMP")
    gmp_redo_idx <- which(colnames(SJCAR19_Lineages_by_CAR[[patient]]) == "GMP-redo")
    after_gmp_idx <- max(gmp_idx, gmp_redo_idx)
    
    if(after_gmp_idx != -Inf && (after_gmp_idx+1) < ncol(SJCAR19_Lineages_by_CAR[[patient]])) {
      ### output directory for the given patient
      outputDir2 <- paste0(outputDir, patient, "/")
      dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
      
      ### interesting indicies
      interesting_idx <- NULL
      for(i in 1:nrow(SJCAR19_Lineages_by_CAR[[patient]])) {
        if(SJCAR19_Lineages_by_CAR[[patient]]$Cell_Type[i] == "CARneg") {
          after_gmp_CARneg_sum <- sum(SJCAR19_Lineages_by_CAR[[patient]][i,after_gmp_idx:(ncol(SJCAR19_Lineages_by_CAR[[patient]])-1)])
          if(after_gmp_CARneg_sum > carneg_threshold) {
            interesting_idx <- c(interesting_idx, i)
          }
        }
      }
      
      ### extract interesting things only
      all_idx <- c(interesting_idx, (interesting_idx+1), (interesting_idx+2))
      all_idx <- all_idx[order(all_idx)]
      interesting_clones[[patient]] <- SJCAR19_Lineages_by_CAR[[patient]][all_idx,]
      
      #
      ### Alluvial plot
      #
      ### filtering for alluvial plot (lineages only & ALL only)
      plot_data_table <- interesting_clones[[patient]][which(interesting_clones[[patient]]$Cell_Type == "ALL"),]
      
      ### get an input data frame for the alluvial plot
      time_points <- colnames(interesting_clones[[patient]])[5:(ncol(interesting_clones[[patient]])-1)]
      total_rows <- length(which(plot_data_table[,time_points] > 0))
      plot_df <- data.frame(Time=rep("", total_rows),
                            Clone_Size=rep(0, total_rows),
                            Clone=rep("", total_rows),
                            CDR3=rep("", total_rows))
      cnt <- 1
      for(i in 1:nrow(plot_data_table)) {
        for(tp in time_points) {
          if(plot_data_table[i,tp] > 0) {
            plot_df[cnt,] <- c(tp,
                               plot_data_table[i,tp],
                               plot_data_table$Clone_ID[i],
                               plot_data_table$CDR3_AA[i])
            cnt <- cnt + 1
          }
        }
      }
      plot_df$Time <- factor(plot_df$Time, levels = time_points)
      
      ### numerize the clone_size column
      plot_df$Clone_Size <- as.numeric(plot_df$Clone_Size)
      
      ### add the number of CARpos cells
      plot_df$CARpos_Num <- ""
      for(i in 1:nrow(plot_df)) {
        num <- interesting_clones[[patient]][paste0(plot_df$Clone[i], "_CARpos"),as.character(plot_df$Time[i])]
        if(num > 0) {
          plot_df$CARpos_Num[i] <- num
        }
      }
      
      ### draw the alluvial plot
      ggplot(plot_df,
             aes(x = Time, stratum = Clone, alluvium = Clone,
                 y = Clone_Size,
                 fill = CDR3, label = CARpos_Num)) +
        ggtitle(paste(patient, "Clonal Tracing (Epitope Spreading - Related)")) +
        geom_stratum(alpha = 1) +
        geom_text(stat = "stratum", size = 3, col = "black") +
        geom_flow() +
        rotate_x_text(90) +
        theme_pubr(legend = "none") +
        theme(axis.title.x = element_blank()) +
        theme_cleveland2() +
        scale_fill_viridis(discrete = T) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
      ggsave(file = paste0(outputDir2, patient, "_Interesting_Clonal_Tracing_CAR.png"), width = 18, height = 9, dpi = 300)
      
      #
      ### annotate VDJDB
      #
      interesting_clones[[patient]]$MHC_A <- ""
      interesting_clones[[patient]]$MHC_B <- ""
      interesting_clones[[patient]]$MHC_Class <- ""
      interesting_clones[[patient]]$Epitope <- ""
      interesting_clones[[patient]]$Epitope_Gene <- ""
      interesting_clones[[patient]]$Epitope_Species <- ""
      interesting_clones[[patient]]$Reference <- ""
      for(i in seq(1, nrow(interesting_clones[[patient]]), 3)) {
        ### get each CDR_AA sequences for each clone
        cdr3_list <- strsplit(interesting_clones[[patient]]$CDR3_AA[i], split = ";", fixed = TRUE)[[1]]
        
        ### get annotation
        for(cdr3 in cdr3_list) {
          target_idx <- which(vdjdb$Combined_CDR3 == cdr3)
          if(length(target_idx) > 0) {
            if(interesting_clones[[patient]]$MHC_A[i] == "") {
              interesting_clones[[patient]]$MHC_A[i] <- vdjdb$`MHC A`[target_idx]
              interesting_clones[[patient]]$MHC_A[i+1] <- vdjdb$`MHC A`[target_idx]
              interesting_clones[[patient]]$MHC_A[i+2] <- vdjdb$`MHC A`[target_idx]
              interesting_clones[[patient]]$MHC_B[i] <- vdjdb$`MHC B`[target_idx]
              interesting_clones[[patient]]$MHC_B[i+1] <- vdjdb$`MHC B`[target_idx]
              interesting_clones[[patient]]$MHC_B[i+2] <- vdjdb$`MHC B`[target_idx]
              interesting_clones[[patient]]$MHC_Class[i] <- vdjdb$`MHC class`[target_idx]
              interesting_clones[[patient]]$MHC_Class[i+1] <- vdjdb$`MHC class`[target_idx]
              interesting_clones[[patient]]$MHC_Class[i+2] <- vdjdb$`MHC class`[target_idx]
              interesting_clones[[patient]]$Epitope[i] <- vdjdb$Epitope[target_idx]
              interesting_clones[[patient]]$Epitope[i+1] <- vdjdb$Epitope[target_idx]
              interesting_clones[[patient]]$Epitope[i+2] <- vdjdb$Epitope[target_idx]
              interesting_clones[[patient]]$Epitope_Gene[i] <- vdjdb$`Epitope gene`[target_idx]
              interesting_clones[[patient]]$Epitope_Gene[i+1] <- vdjdb$`Epitope gene`[target_idx]
              interesting_clones[[patient]]$Epitope_Gene[i+2] <- vdjdb$`Epitope gene`[target_idx]
              interesting_clones[[patient]]$Epitope_Species[i] <- vdjdb$`Epitope species`[target_idx]
              interesting_clones[[patient]]$Epitope_Species[i+1] <- vdjdb$`Epitope species`[target_idx]
              interesting_clones[[patient]]$Epitope_Species[i+2] <- vdjdb$`Epitope species`[target_idx]
              interesting_clones[[patient]]$Reference[i] <- vdjdb$Reference[target_idx]
              interesting_clones[[patient]]$Reference[i+1] <- vdjdb$Reference[target_idx]
              interesting_clones[[patient]]$Reference[i+2] <- vdjdb$Reference[target_idx]
            } else {
              interesting_clones[[patient]]$MHC_A[i] <- paste0(interesting_clones[[patient]]$MHC_A[i], ";",
                                                               vdjdb$`MHC A`[target_idx])
              interesting_clones[[patient]]$MHC_A[i+1] <- paste0(interesting_clones[[patient]]$MHC_A[i+1], ";",
                                                                 vdjdb$`MHC A`[target_idx])
              interesting_clones[[patient]]$MHC_A[i+2] <- paste0(interesting_clones[[patient]]$MHC_A[i+2], ";",
                                                                 vdjdb$`MHC A`[target_idx])
              interesting_clones[[patient]]$MHC_B[i] <- paste0(interesting_clones[[patient]]$MHC_B[i], ";",
                                                               vdjdb$`MHC B`[target_idx])
              interesting_clones[[patient]]$MHC_B[i+1] <- paste0(interesting_clones[[patient]]$MHC_B[i+1], ";",
                                                                 vdjdb$`MHC B`[target_idx])
              interesting_clones[[patient]]$MHC_B[i+2] <- paste0(interesting_clones[[patient]]$MHC_B[i+2], ";",
                                                                 vdjdb$`MHC B`[target_idx])
              interesting_clones[[patient]]$MHC_Class[i] <- paste0(interesting_clones[[patient]]$MHC_Class[i], ";",
                                                                   vdjdb$`MHC class`[target_idx])
              interesting_clones[[patient]]$MHC_Class[i+1] <- paste0(interesting_clones[[patient]]$MHC_Class[i+1], ";",
                                                                     vdjdb$`MHC class`[target_idx])
              interesting_clones[[patient]]$MHC_Class[i+2] <- paste0(interesting_clones[[patient]]$MHC_Class[i+2], ";",
                                                                     vdjdb$`MHC class`[target_idx])
              interesting_clones[[patient]]$Epitope[i] <- paste0(interesting_clones[[patient]]$Epitope[i], ";",
                                                                 vdjdb$Epitope[target_idx])
              interesting_clones[[patient]]$Epitope[i+1] <- paste0(interesting_clones[[patient]]$Epitope[i+1], ";",
                                                                   vdjdb$Epitope[target_idx])
              interesting_clones[[patient]]$Epitope[i+2] <- paste0(interesting_clones[[patient]]$Epitope[i+2], ";",
                                                                   vdjdb$Epitope[target_idx])
              interesting_clones[[patient]]$Epitope_Gene[i] <- paste0(interesting_clones[[patient]]$Epitope_Gene[i], ";",
                                                                      vdjdb$`Epitope gene`[target_idx])
              interesting_clones[[patient]]$Epitope_Gene[i+1] <- paste0(interesting_clones[[patient]]$Epitope_Gene[i+1], ";",
                                                                        vdjdb$`Epitope gene`[target_idx])
              interesting_clones[[patient]]$Epitope_Gene[i+2] <- paste0(interesting_clones[[patient]]$Epitope_Gene[i+2], ";",
                                                                        vdjdb$`Epitope gene`[target_idx])
              interesting_clones[[patient]]$Epitope_Species[i] <- paste0(interesting_clones[[patient]]$Epitope_Species[i], ";",
                                                                         vdjdb$`Epitope species`[target_idx])
              interesting_clones[[patient]]$Epitope_Species[i+1] <- paste0(interesting_clones[[patient]]$Epitope_Species[i+1], ";",
                                                                           vdjdb$`Epitope species`[target_idx])
              interesting_clones[[patient]]$Epitope_Species[i+2] <- paste0(interesting_clones[[patient]]$Epitope_Species[i+2], ";",
                                                                           vdjdb$`Epitope species`[target_idx])
              interesting_clones[[patient]]$Reference[i] <- paste0(interesting_clones[[patient]]$Reference[i], ";",
                                                                   vdjdb$Reference[target_idx])
              interesting_clones[[patient]]$Reference[i+1] <- paste0(interesting_clones[[patient]]$Reference[i+1], ";",
                                                                     vdjdb$Reference[target_idx])
              interesting_clones[[patient]]$Reference[i+2] <- paste0(interesting_clones[[patient]]$Reference[i+2], ";",
                                                                     vdjdb$Reference[target_idx])
            }
          }
        }
      }
      
      ### save the result
      write.xlsx2(interesting_clones[[patient]],
                  file = paste0(outputDir, "Interesting_Lineages_by_CAR.xlsx"),
                  sheetName = patient,
                  row.names = FALSE,
                  append = TRUE)
      gc()
      
    }
  }
  
  ### split the CAR+/CAR- lineage table to convey whether a cell is a CD4/CD8
  
  ### make a new object for that
  SJCAR19_Lineages_in_Full <- vector("list", length = length(SJCAR19_Lineages_by_CAR))
  names(SJCAR19_Lineages_in_Full) <- names(SJCAR19_Lineages_by_CAR)
  for(patient in names(SJCAR19_Lineages_in_Full)) {
    ### add CD4/CD8 sections
    ### CAR+/CAR-/ALL x CD4/CD8/ALL - 9 rows per clone
    SJCAR19_Lineages_in_Full[[patient]] <- apply(SJCAR19_Lineages_by_CAR[[patient]], 2, function(x) rbind(x, x, x))
    
    ### column name change 'Cell_Type' -> 'CAR_Type'
    colnames(SJCAR19_Lineages_in_Full[[patient]])[1] <- "CAR_Type"
    
    ### add a column to describe CD4/CD8 info
    SJCAR19_Lineages_in_Full[[patient]] <- data.frame(CD_Type=rep(c("CD4", "CD8", "ALL"), nrow(SJCAR19_Lineages_by_CAR[[patient]])),
                                                      SJCAR19_Lineages_in_Full[[patient]],
                                                      stringsAsFactors = FALSE, check.names = FALSE)
    
    ### numerize the numeric columns
    SJCAR19_Lineages_in_Full[[patient]][,colnames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])] <- sapply(SJCAR19_Lineages_in_Full[[patient]][,colnames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])], as.numeric)
    
    ### get time points
    time_points <- colnames(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])[1:(ncol(SJCAR19_Clonotype_Frequency[["ALL"]][[patient]])-1)]
    
    ### fill out the table
    for(i in 1:nrow(SJCAR19_Lineages_in_Full[[patient]])) {
      if(SJCAR19_Lineages_in_Full[[patient]]$CD_Type[i] != "ALL") {
        if(SJCAR19_Lineages_in_Full[[patient]]$CAR_Type[i] == "ALL") {
          clone_idx <- intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient == SJCAR19_Lineages_in_Full[[patient]]$Clone_ID[i]),
                                 which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == SJCAR19_Lineages_in_Full[[patient]]$CD_Type[i]))
        } else {
          clone_idx <- intersect(which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == SJCAR19_Lineages_in_Full[[patient]]$CD_Type[i]),
                                 intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient == SJCAR19_Lineages_in_Full[[patient]]$Clone_ID[i]),
                                           which(Seurat_Obj@meta.data$CAR == SJCAR19_Lineages_in_Full[[patient]]$CAR_Type[i])))
        }
        
        for(tp in time_points) {
          SJCAR19_Lineages_in_Full[[patient]][i,tp] <- length(intersect(clone_idx,
                                                                        which(Seurat_Obj@meta.data$time == tp)))
        }
        SJCAR19_Lineages_in_Full[[patient]]$Total[i] <- sum(SJCAR19_Lineages_in_Full[[patient]][i,time_points])
      }
    }
    
    ### annotate VDJDB
    SJCAR19_Lineages_in_Full[[patient]]$MHC_A <- ""
    SJCAR19_Lineages_in_Full[[patient]]$MHC_B <- ""
    SJCAR19_Lineages_in_Full[[patient]]$MHC_Class <- ""
    SJCAR19_Lineages_in_Full[[patient]]$Epitope <- ""
    SJCAR19_Lineages_in_Full[[patient]]$Epitope_Gene <- ""
    SJCAR19_Lineages_in_Full[[patient]]$Epitope_Species <- ""
    SJCAR19_Lineages_in_Full[[patient]]$Reference <- ""
    CDR3_AA <- unique(SJCAR19_Lineages_in_Full[[patient]]$CDR3_AA)
    for(i in 1:length(CDR3_AA)) {
      ### get each CDR_AA sequences for each clone
      cdr3_list <- strsplit(CDR3_AA[i], split = ";", fixed = TRUE)[[1]]
      
      ### get index of the clone in the table
      table_idx <- which(SJCAR19_Lineages_in_Full[[patient]]$CDR3_AA == CDR3_AA[i])
      
      ### get annotation
      for(cdr3 in cdr3_list) {
        ### get index of the CDR3 in the vdj
        target_idx <- which(vdjdb$Combined_CDR3 == cdr3)
        if(length(target_idx) > 0) {
          if(SJCAR19_Lineages_in_Full[[patient]]$MHC_A[table_idx[1]] == "") {
            SJCAR19_Lineages_in_Full[[patient]]$MHC_A[table_idx] <- vdjdb$`MHC A`[target_idx]
            SJCAR19_Lineages_in_Full[[patient]]$MHC_B[table_idx] <- vdjdb$`MHC B`[target_idx]
            SJCAR19_Lineages_in_Full[[patient]]$MHC_Class[table_idx] <- vdjdb$`MHC class`[target_idx]
            SJCAR19_Lineages_in_Full[[patient]]$Epitope[table_idx] <- vdjdb$Epitope[target_idx]
            SJCAR19_Lineages_in_Full[[patient]]$Epitope_Gene[table_idx] <- vdjdb$`Epitope gene`[target_idx]
            SJCAR19_Lineages_in_Full[[patient]]$Epitope_Species[table_idx] <- vdjdb$`Epitope species`[target_idx]
            SJCAR19_Lineages_in_Full[[patient]]$Reference[table_idx] <- vdjdb$Reference[target_idx]
          } else {
            SJCAR19_Lineages_in_Full[[patient]]$MHC_A[table_idx] <- paste0(SJCAR19_Lineages_in_Full[[patient]]$MHC_A[table_idx], ";",
                                                                           vdjdb$`MHC A`[target_idx])
            SJCAR19_Lineages_in_Full[[patient]]$MHC_B[table_idx] <- paste0(SJCAR19_Lineages_in_Full[[patient]]$MHC_B[table_idx], ";",
                                                                           vdjdb$`MHC B`[target_idx])
            SJCAR19_Lineages_in_Full[[patient]]$MHC_Class[table_idx] <- paste0(SJCAR19_Lineages_in_Full[[patient]]$MHC_Class[table_idx], ";",
                                                                               vdjdb$`MHC class`[target_idx])
            SJCAR19_Lineages_in_Full[[patient]]$Epitope[table_idx] <- paste0(SJCAR19_Lineages_in_Full[[patient]]$Epitope[table_idx], ";",
                                                                             vdjdb$Epitope[target_idx])
            SJCAR19_Lineages_in_Full[[patient]]$Epitope_Gene[table_idx] <- paste0(SJCAR19_Lineages_in_Full[[patient]]$Epitope_Gene[table_idx], ";",
                                                                                  vdjdb$`Epitope gene`[target_idx])
            SJCAR19_Lineages_in_Full[[patient]]$Epitope_Species[table_idx] <- paste0(SJCAR19_Lineages_in_Full[[patient]]$Epitope_Species[table_idx], ";",
                                                                                     vdjdb$`Epitope species`[target_idx])
            SJCAR19_Lineages_in_Full[[patient]]$Reference[table_idx] <- paste0(SJCAR19_Lineages_in_Full[[patient]]$Reference[table_idx], ";",
                                                                               vdjdb$Reference[target_idx])
          }
        }
      }
    }
    
    #
    ### lineages in CAR- cells only, but annotate CD4 cells in the alluvial plot
    #
    
    ### get index after GMP
    gmp_idx <- which(colnames(SJCAR19_Lineages_in_Full[[patient]]) == "GMP")
    gmp_redo_idx <- which(colnames(SJCAR19_Lineages_in_Full[[patient]]) == "GMP-redo")
    after_gmp_idx <- max(gmp_idx, gmp_redo_idx)
    total_idx <- which(colnames(SJCAR19_Lineages_in_Full[[patient]]) == "Total")
    
    if(after_gmp_idx != -Inf && (after_gmp_idx+1) < total_idx) {
      ### get the interesting clones only & CAR- cells only
      plot_data_table <- SJCAR19_Lineages_in_Full[[patient]][intersect(intersect(which(SJCAR19_Lineages_in_Full[[patient]]$CAR_Type == "CARneg"),
                                                                                 which(SJCAR19_Lineages_in_Full[[patient]]$CD_Type == "ALL")),
                                                                       which(SJCAR19_Lineages_in_Full[[patient]]$Clone_ID %in% unique(interesting_clones[[patient]]$Clone_ID))),]
      
      ### get an input data frame for the alluvial plot
      time_points <- colnames(SJCAR19_Lineages_in_Full[[patient]])[6:(ncol(SJCAR19_Lineages_in_Full[[patient]])-8)]
      total_rows <- length(which(plot_data_table[,time_points] > 0))
      plot_df <- data.frame(Time=rep("", total_rows),
                            Clone_Size=rep(0, total_rows),
                            Clone=rep("", total_rows),
                            CDR3=rep("", total_rows))
      cnt <- 1
      for(i in 1:nrow(plot_data_table)) {
        for(tp in time_points) {
          if(plot_data_table[i,tp] > 0) {
            plot_df[cnt,] <- c(tp,
                               plot_data_table[i,tp],
                               plot_data_table$Clone_ID[i],
                               plot_data_table$CDR3_AA[i])
            cnt <- cnt + 1
          }
        }
      }
      plot_df$Time <- factor(plot_df$Time, levels = time_points)
      
      ### numerize the clone_size column
      plot_df$Clone_Size <- as.numeric(plot_df$Clone_Size)
      
      ### add the number of CD4 cells
      plot_df$CD4_Num <- ""
      for(i in 1:nrow(plot_df)) {
        t_idx <- intersect(intersect(which(SJCAR19_Lineages_in_Full[[patient]]$Clone_ID == plot_df$Clone[i]),
                                     which(SJCAR19_Lineages_in_Full[[patient]]$CAR_Type == "CARneg")),
                           which(SJCAR19_Lineages_in_Full[[patient]]$CD_Type == "CD4"))
        num <- SJCAR19_Lineages_in_Full[[patient]][t_idx,as.character(plot_df$Time[i])]
        if(num > 0) {
          plot_df$CD4_Num[i] <- num
        }
      }
      
      ### draw the alluvial plot
      ggplot(plot_df,
             aes(x = Time, stratum = Clone, alluvium = Clone,
                 y = Clone_Size,
                 fill = CDR3, label = CD4_Num)) +
        ggtitle(paste(patient, "Clonal Tracing (Epitope Spreading) CAR- Only")) +
        geom_stratum(alpha = 1) +
        geom_text(stat = "stratum", size = 3, col = "black") +
        geom_flow() +
        rotate_x_text(90) +
        theme_pubr(legend = "none") +
        theme(axis.title.x = element_blank()) +
        theme_cleveland2() +
        scale_fill_viridis(discrete = T) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
      ggsave(file = paste0(outputDir, patient, "/", patient, "_Interesting_Clonal_Tracing_CARneg_Only.png"), width = 18, height = 9, dpi = 300)
    }
  }
  
  ### save the result
  saveRDS(SJCAR19_Lineages_in_Full, file = paste0(outputDir, "SJCAR19_Lineages_in_Full.RDS"))
  for(patient in names(SJCAR19_Lineages_in_Full)) {
    write.xlsx2(SJCAR19_Lineages_in_Full[[patient]],
                file = paste0(outputDir, "SJCAR19_Lineages_in_Full.xlsx"),
                sheetName = patient,
                row.names = FALSE,
                append = TRUE)
    gc()
  }
  
}
