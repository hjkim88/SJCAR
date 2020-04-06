###
#   File name : UMAP_with_GEX_TCR.R
#   Author    : Hyunjin Kim
#   Date      : Apr 3, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Draw a UMAP plot that shows how TCR clonality changes over time
#               and besides that, how GEX clusters change over time.
#
#   * There are some important questions here:
#     1. Which clonotypes appeared most across all the time points?
#     2. Clonotype frequency. This is different from the (1) because
#        the (1) examines the numer of "time points" that a clonotype appeared in,
#        while this examines the number of total appearance of a "clonotype".
#     3. When we show the (1) with a UMAP plot, do they show clonal lineage?
#        A UMAP split by all the time points
#     4. Similar to the (3), but using the (2) this time.
#     5. Can we show clonal lineage before and after GMP based on CAR+ cells?
#
#   Instruction
#               1. Source("UMAP_with_GEX_TCR.R")
#               2. Run the function "umap_gex_tcr" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_UMAP_with_GEX_TCR.R/UMAP_with_GEX_TCR.R")
#               > umap_gex_tcr(Seurat_RObj_path="./data/JCC212_Px5_TCR_clonotyped.Robj",
#                              outputDir="./results/")
###

umap_gex_tcr <- function(Seurat_RObj_path="./data/JCC212_Px5_TCR_clonotyped.Robj",
                         outputDir="./results/") {
  
  ### load library
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
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    require(scales, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### STRICT VERSION
  ### calculate the followings:
  ### appearance - compute which clonotype appears most "across different time points"
  ### frequency - compute freqency of the clonotypes regardless of the time points
  unique_clonotypes_strict <- unique(Seurat_Obj@meta.data$global_clonotype_strict[which(!is.na(Seurat_Obj@meta.data$global_clonotype_strict))])
  clone_info_strict <- data.frame(clonotype=unique_clonotypes_strict,
                                  appearance=NA,
                                  appearance_details=NA,
                                  appearance_car_pos=NA,
                                  appearance_car_pos_details=NA,
                                  frequency=NA,
                                  frequency_car_pos=NA,
                                  frequency_gmp=NA,
                                  tcr_seq=NA,
                                  stringsAsFactors = FALSE, check.names = FALSE)
  rownames(clone_info_strict) <- clone_info_strict$clonotype
  for(clone in rownames(clone_info_strict)) {
    clone_info_strict[clone,"appearance"] <- length(unique(Seurat_Obj@meta.data$Time[which(Seurat_Obj@meta.data$global_clonotype_strict == clone)]))
  }
  for(clone in rownames(clone_info_strict)) {
    temp <- unique(Seurat_Obj@meta.data$TimeF[which(Seurat_Obj@meta.data$global_clonotype_strict == clone)])
    clone_info_strict[clone,"appearance_details"] <- paste(temp[order(temp)], collapse = ";")
  }
  for(clone in rownames(clone_info_strict)) {
    clone_info_strict[clone,"appearance_car_pos"] <- length(unique(Seurat_Obj@meta.data$Time[intersect(which(Seurat_Obj@meta.data$global_clonotype_strict == clone),
                                                                                                       which(Seurat_Obj@meta.data$CAR == "CARpos"))]))
  }
  for(clone in rownames(clone_info_strict)) {
    temp <- unique(Seurat_Obj@meta.data$TimeF[intersect(which(Seurat_Obj@meta.data$global_clonotype_strict == clone),
                                                        which(Seurat_Obj@meta.data$CAR == "CARpos"))])
    clone_info_strict[clone,"appearance_car_pos_details"] <- paste(temp[order(temp)], collapse = ";")
  }
  for(clone in rownames(clone_info_strict)) {
    clone_info_strict[clone,"frequency"] <- length(which(Seurat_Obj@meta.data$global_clonotype_strict == clone))
  }
  for(clone in rownames(clone_info_strict)) {
    clone_info_strict[clone,"frequency_car_pos"] <- length(intersect(which(Seurat_Obj@meta.data$global_clonotype_strict == clone),
                                                                      which(Seurat_Obj@meta.data$CAR == "CARpos")))
  }
  for(clone in rownames(clone_info_strict)) {
    clone_info_strict[clone,"frequency_gmp"] <- length(intersect(which(Seurat_Obj@meta.data$global_clonotype_strict == clone),
                                                                      which(Seurat_Obj@meta.data$CAR == "GMP")))
  }
  for(clone in rownames(clone_info_strict)) {
    clone_info_strict[clone,"tcr_seq"] <- Seurat_Obj@meta.data$global_cdr3_aa_strict[which(Seurat_Obj@meta.data$global_clonotype_strict == clone)[1]]
  }
  ### order the data frame by frequency_car_pos and save it
  clone_info_strict <- clone_info_strict[order(-clone_info_strict$frequency_car_pos,
                                               -clone_info_strict$appearance_car_pos),]
  write.xlsx2(clone_info_strict, file = paste0(outputDir, "clone_info_strict.xlsx"),
              sheetName = "clone_info_strict", row.names = FALSE)
  
  ### LENIENT VERSION
  ### calculate the followings:
  ### appearance - compute which clonotype appears most "across different time points"
  ### frequency - compute freqency of the clonotypes regardless of the time points
  unique_clonotypes_lenient <- unique(Seurat_Obj@meta.data$global_clonotype_lenient[which(!is.na(Seurat_Obj@meta.data$global_clonotype_lenient))])
  clone_info_lenient <- data.frame(clonotype=unique_clonotypes_lenient,
                                   appearance=NA,
                                   appearance_details=NA,
                                   appearance_car_pos=NA,
                                   appearance_car_pos_details=NA,
                                   frequency=NA,
                                   frequency_car_pos=NA,
                                   frequency_gmp=NA,
                                   tcr_seq=NA,
                                   stringsAsFactors = FALSE, check.names = FALSE)
  rownames(clone_info_lenient) <- clone_info_lenient$clonotype
  for(clone in rownames(clone_info_lenient)) {
    clone_info_lenient[clone,"appearance"] <- length(unique(Seurat_Obj@meta.data$Time[which(Seurat_Obj@meta.data$global_clonotype_lenient == clone)]))
  }
  for(clone in rownames(clone_info_lenient)) {
    temp <- unique(Seurat_Obj@meta.data$TimeF[which(Seurat_Obj@meta.data$global_clonotype_lenient == clone)])
    clone_info_lenient[clone,"appearance_details"] <- paste(temp[order(temp)], collapse = ";")
  }
  for(clone in rownames(clone_info_lenient)) {
    clone_info_lenient[clone,"appearance_car_pos"] <- length(unique(Seurat_Obj@meta.data$Time[intersect(which(Seurat_Obj@meta.data$global_clonotype_lenient == clone),
                                                                                                       which(Seurat_Obj@meta.data$CAR == "CARpos"))]))
  }
  for(clone in rownames(clone_info_lenient)) {
    temp <- unique(Seurat_Obj@meta.data$TimeF[intersect(which(Seurat_Obj@meta.data$global_clonotype_lenient == clone),
                                                        which(Seurat_Obj@meta.data$CAR == "CARpos"))])
    clone_info_lenient[clone,"appearance_car_pos_details"] <- paste(temp[order(temp)], collapse = ";")
  }
  for(clone in rownames(clone_info_lenient)) {
    clone_info_lenient[clone,"frequency"] <- length(which(Seurat_Obj@meta.data$global_clonotype_lenient == clone))
  }
  for(clone in rownames(clone_info_lenient)) {
    clone_info_lenient[clone,"frequency_car_pos"] <- length(intersect(which(Seurat_Obj@meta.data$global_clonotype_lenient == clone),
                                                                      which(Seurat_Obj@meta.data$CAR == "CARpos")))
  }
  for(clone in rownames(clone_info_lenient)) {
    clone_info_lenient[clone,"frequency_gmp"] <- length(intersect(which(Seurat_Obj@meta.data$global_clonotype_lenient == clone),
                                                                       which(Seurat_Obj@meta.data$CAR == "GMP")))
  }
  for(clone in rownames(clone_info_lenient)) {
    clone_info_lenient[clone,"tcr_seq"] <- Seurat_Obj@meta.data$global_cdr3_aa_lenient[which(Seurat_Obj@meta.data$global_clonotype_lenient == clone)[1]]
  }
  ### order the data frame by frequency_car_pos and save it
  clone_info_lenient <- clone_info_lenient[order(-clone_info_lenient$frequency_car_pos,
                                                 -clone_info_lenient$appearance_car_pos),]
  write.xlsx2(clone_info_lenient, file = paste0(outputDir, "clone_info_lenient.xlsx"),
              sheetName = "clone_info_lenient", row.names = FALSE)
  
  
  ### create a data frame for ggplot
  time <- factor(as.character(Seurat_Obj@meta.data$TimeF),
                  levels = intersect(levels(Seurat_Obj@meta.data$TimeF), unique(Seurat_Obj@meta.data$TimeF)))
  cluster <- factor(as.character(Seurat_Obj@meta.data$seurat_clusters),
                    levels = intersect(levels(Seurat_Obj@meta.data$seurat_clusters), unique(Seurat_Obj@meta.data$seurat_clusters)))
  car <- factor(Seurat_Obj@meta.data$CAR,
                levels = unique(Seurat_Obj@meta.data$CAR)[order(unique(Seurat_Obj@meta.data$CAR))])
  clonotype_strict <- factor(Seurat_Obj@meta.data$global_clonotype_strict,
                             levels = unique(Seurat_Obj@meta.data$global_clonotype_strict)[order(unique(Seurat_Obj@meta.data$global_clonotype_strict))])
  clonotype_lenient <- factor(Seurat_Obj@meta.data$global_clonotype_lenient,
                              levels = unique(Seurat_Obj@meta.data$global_clonotype_lenient)[order(unique(Seurat_Obj@meta.data$global_clonotype_lenient))])
  temp <- merge(Seurat_Obj@meta.data[,c("GexCellFull", "global_clonotype_strict")], clone_info_strict,
                by.x = "global_clonotype_strict", by.y = "clonotype",
                all.x = TRUE)
  rownames(temp) <- temp$GexCellFull
  clone_size_strict <- as.numeric(temp[rownames(Seurat_Obj@meta.data),"frequency"])
  clone_size_car_pos_strict <- as.numeric(temp[rownames(Seurat_Obj@meta.data),"frequency_car_pos"])
  clone_size_gmp_strict <- as.numeric(temp[rownames(Seurat_Obj@meta.data),"frequency_gmp"])
  temp <- merge(Seurat_Obj@meta.data[,c("GexCellFull", "global_clonotype_lenient")], clone_info_lenient,
                by.x = "global_clonotype_lenient", by.y = "clonotype",
                all.x = TRUE)
  rownames(temp) <- temp$GexCellFull
  clone_size_lenient <- as.numeric(temp[rownames(Seurat_Obj@meta.data),"frequency"])
  clone_size_car_pos_lenient <- as.numeric(temp[rownames(Seurat_Obj@meta.data),"frequency_car_pos"])
  clone_size_gmp_lenient <- as.numeric(temp[rownames(Seurat_Obj@meta.data),"frequency_gmp"])
  rm(temp)
  gc()
  plot_df <- data.frame(X=Seurat_Obj@meta.data$Full_UMAP1,
                        Y=Seurat_Obj@meta.data$Full_UMAP2,
                        Time=time,
                        Cluster=cluster,
                        Clonotype_strict=clonotype_strict,
                        Clonotype_lenient=clonotype_lenient,
                        Clone_size_strict=clone_size_strict,
                        Clone_size_lenient=clone_size_lenient,
                        Clone_size_car_pos_strict=clone_size_car_pos_strict,
                        Clone_size_car_pos_lenient=clone_size_car_pos_lenient,
                        Clone_size_gmp_strict=clone_size_gmp_strict,
                        Clone_size_gmp_lenient=clone_size_gmp_lenient,
                        CAR=car)
  
  ### only retain the cells with TCR info
  plot_df <- plot_df[which(!is.na(Seurat_Obj@meta.data$cdr3_aa)),]
  
  ### set min-max values for x-y axes
  x_range <- c(min(plot_df$X), max(plot_df$X))
  y_range <- c(min(plot_df$Y), max(plot_df$Y))
  
  ### set unique time points
  unique_time <- unique(plot_df$Time)
  unique_time <- unique_time[order(unique_time)]
  
  ### for each attribute, print a plot
  for(attr in colnames(plot_df)[-c(1,2,3,4)]) {
    
    ### One plot with all the groups + each group
    p <- vector("list", length(unique_time)+1)
    
    ### one plot for all the time points colored by the given attribute
    if(is.numeric(plot_df[,attr]) || attr == "CAR") {
      ### plot for all the time points
      p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col=attr), size=1.5, alpha=0.5) +
        xlab("UMAP1") + ylab("UMAP2") +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle("All") +
        theme_classic(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5))
    } else if(is.factor(plot_df[,attr])) {
      ### set colors for each group
      col_palette <- hue_pal()(length(unique(plot_df[,attr])))
      names(col_palette) <- unique(plot_df[,attr])
      
      ### plot for all the time points
      p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col=attr), size=1.5, alpha=0.5) +
        xlab("UMAP1") + ylab("UMAP2") +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle("All") +
        theme_classic(base_size = 16) +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    } else {
      stop("ERROR: the attribute is not a factor nor a numeric value")
    }
    
    ### separated plot for each time point
    for(i in 1:length(unique_time)) {
      if(is.numeric(plot_df[,attr]) || attr == "CAR") {
        temp <- unique(plot_df[which(plot_df$Time == unique_time[i]),attr])
        if((length(temp) == 1) && is.na(temp)) {
          p[[i+1]] <- ggplot(plot_df[which(plot_df$Time == unique_time[i]),], aes_string(x="X", y="Y")) +
            geom_point(col="grey50", size=1.5, alpha=0.5) +
            xlab("UMAP1") + ylab("UMAP2") +
            xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
            ggtitle(unique_time[i]) +
            theme_classic(base_size = 16) +
            theme(plot.title = element_text(hjust = 0.5))
        } else {
          p[[i+1]] <- ggplot(plot_df[which(plot_df$Time == unique_time[i]),], aes_string(x="X", y="Y")) +
            geom_point(aes_string(col=attr), size=1.5, alpha=0.5) +
            xlab("UMAP1") + ylab("UMAP2") +
            xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
            ggtitle(unique_time[i]) +
            theme_classic(base_size = 16) +
            theme(plot.title = element_text(hjust = 0.5))
        }
      } else if(is.factor(plot_df[,attr])) {
        temp <- unique(plot_df[which(plot_df$Time == unique_time[i]),attr])
        if((length(temp) == 1) && is.na(temp)) {
          p[[i+1]] <- ggplot(plot_df[which(plot_df$Time == unique_time[i]),], aes_string(x="X", y="Y")) +
            geom_point(col="grey50", size=1.5, alpha=0.5) +
            xlab("UMAP1") + ylab("UMAP2") +
            xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
            ggtitle(unique_time[i]) +
            theme_classic(base_size = 16) +
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
        } else {
          p[[i+1]] <- ggplot(plot_df[which(plot_df$Time == unique_time[i]),], aes_string(x="X", y="Y")) +
            geom_point(aes_string(col=attr), size=1.5, alpha=0.5) +
            xlab("UMAP1") + ylab("UMAP2") +
            xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
            ggtitle(unique_time[i]) +
            theme_classic(base_size = 16) +
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
        }
      } else {
        stop("ERROR: the attribute is not a factor nor a numeric value")
      }
    }
    
    ### arrange the plots and save
    fName <- paste0("UMAP_", attr)
    g <- arrangeGrob(grobs = p,
                     nrow = 3,
                     ncol = 4,
                     top = fName)
    ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 22, height = 12, dpi = 300)
    
  }
  
}
