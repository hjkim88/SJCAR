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
  if(!require(ggsci, quietly = TRUE)) {
    install.packages("ggsci")
    require(ggsci, quietly = TRUE)
  }
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### global clonotypes in the Seurat objects
  global_clonotypes <- colnames(Seurat_Obj@meta.data)[grep("global_clonotype", colnames(Seurat_Obj@meta.data), fixed = TRUE)]
  
  ### for each patient, make a statistics table
  for(patient in unique(Seurat_Obj@meta.data$Px)) {
    
    ### print progress
    writeLines(paste(patient))
    
    ### indicies that assign to the given patient
    px_idx <- which(Seurat_Obj@meta.data$Px == patient)
    
    ### remove NA, "NA", and ""
    px_idx <- intersect(px_idx, intersect(intersect(which(!is.na(Seurat_Obj@meta.data$cdr3_nt)),
                                                    which(Seurat_Obj@meta.data$cdr3_nt != "NA")),
                                          which(Seurat_Obj@meta.data$cdr3_nt != "")))
    
    if(length(px_idx) > 0) {
      
      ### output directory for each patient
      outputDir2 <- paste0(outputDir, "/", patient, "/")
      dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
      
      ### set min-max values for x-y axes (will be used later)
      x_range <- c(min(Seurat_Obj@meta.data$Full_UMAP1[px_idx]), max(Seurat_Obj@meta.data$Full_UMAP1[px_idx]))
      y_range <- c(min(Seurat_Obj@meta.data$Full_UMAP2[px_idx]), max(Seurat_Obj@meta.data$Full_UMAP2[px_idx]))
      
      ### make an empty lineage table (will be used later)
      lineage_abstract_table <- matrix(0, 4, length(global_clonotypes))
      rownames(lineage_abstract_table) <- c("Lineages", "Lineages_CARpos",
                                            "Clonotypes", "Clonotypes_CARpos")
      colnames(lineage_abstract_table) <- global_clonotypes
      
      ### make a table that shows how clonotype frequency (and proportion) changes over time
      for(type in global_clonotypes) {
        target_idx <- intersect(px_idx, which(!is.na(Seurat_Obj@meta.data[,type])))
        
        dups <- Seurat_Obj@meta.data[target_idx, type][which(duplicated(Seurat_Obj@meta.data[target_idx, type]))]
        unique_clonotypes <- unique(dups)
        unique_clonotypes <- unique_clonotypes[intersect(intersect(which(!is.na(unique_clonotypes)),
                                                                   which(unique_clonotypes != "NA")),
                                                         which(unique_clonotypes != ""))]
        
        unique_time_points <- levels(Seurat_Obj@meta.data$TimeF[target_idx])
        
        frequency_over_time <- matrix(0, length(unique_clonotypes), length(unique_time_points))
        rownames(frequency_over_time) <- unique_clonotypes
        colnames(frequency_over_time) <- unique_time_points
        
        car_frequency_over_time <- matrix(0, length(unique_clonotypes), length(unique_time_points))
        rownames(car_frequency_over_time) <- unique_clonotypes
        colnames(car_frequency_over_time) <- unique_time_points
        
        proportion_over_time <- matrix(0, length(unique_clonotypes), length(unique_time_points))
        rownames(proportion_over_time) <- unique_clonotypes
        colnames(proportion_over_time) <- unique_time_points
        
        car_proportion_over_time <- matrix(0, length(unique_clonotypes), length(unique_time_points))
        rownames(car_proportion_over_time) <- unique_clonotypes
        colnames(car_proportion_over_time) <- unique_time_points
        
        ### start time
        start_time <- Sys.time()
        
        ### set progress bar
        pb <- txtProgressBar(min = 0, max = length(unique_clonotypes)*length(unique_time_points), style = 3)
        
        ### fill out the frequency tables
        cnt <- 0
        for(time in unique_time_points) {
          for(clonotype in unique_clonotypes) {
            ### indicies for the specific patient, time, and clonotype
            target_idx <- intersect(intersect(px_idx,
                                              which(Seurat_Obj@meta.data$Time == time)),
                                    which(Seurat_Obj@meta.data[,type] == clonotype))
            
            frequency_over_time[clonotype,time] <- length(target_idx)
            car_frequency_over_time[clonotype,time] <- length(intersect(target_idx,
                                                                        which(Seurat_Obj@meta.data$CAR == "CARpos")))
            
            ### update the progress bar
            cnt <- cnt + 1
            setTxtProgressBar(pb, cnt)
          }
          
          ### fill out the proportion tables
          cells_time_num <- length(intersect(px_idx,
                                             which(Seurat_Obj@meta.data$Time == time)))
          car_cells_time_num <- length(intersect(intersect(px_idx,
                                                           which(Seurat_Obj@meta.data$Time == time)),
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
        frequency_over_time <- data.frame(Clonotype=unique_clonotypes,
                                          frequency_over_time,
                                          Total=apply(frequency_over_time, 1, sum),
                                          stringsAsFactors = FALSE, check.names = FALSE)
        car_frequency_over_time <- data.frame(Clonotype=unique_clonotypes,
                                              car_frequency_over_time,
                                              Total=apply(car_frequency_over_time, 1, sum),
                                              stringsAsFactors = FALSE, check.names = FALSE)
        proportion_over_time <- data.frame(Clonotype=unique_clonotypes,
                                           proportion_over_time,
                                           Total=apply(proportion_over_time, 1, sum),
                                           stringsAsFactors = FALSE, check.names = FALSE)
        car_proportion_over_time <- data.frame(Clonotype=unique_clonotypes,
                                               car_proportion_over_time,
                                               Total=apply(car_proportion_over_time, 1, sum),
                                               stringsAsFactors = FALSE, check.names = FALSE)
        
        ### order the data frames by the Total column
        frequency_over_time <- frequency_over_time[order(-frequency_over_time[,"Total"]),]
        proportion_over_time <- proportion_over_time[order(-proportion_over_time[,"Total"]),]
        car_frequency_over_time <- car_frequency_over_time[order(-car_frequency_over_time[,"Total"]),]
        car_proportion_over_time <- car_proportion_over_time[order(-car_proportion_over_time[,"Total"]),]
        
        ### trim the car tables
        car_proportion_over_time <- car_proportion_over_time[which(car_frequency_over_time[,"Total"] > 1),]
        car_frequency_over_time <- car_frequency_over_time[which(car_frequency_over_time[,"Total"] > 1),]
        
        ### save the tables in Excel format
        write.xlsx2(frequency_over_time, file = paste0(outputDir2, "clonotype_frequency_over_time_", patient, ".xlsx"),
                    sheetName = type, row.names = FALSE, append = TRUE)
        gc()
        write.xlsx2(proportion_over_time, file = paste0(outputDir2, "clonotype_proportion_over_time_", patient, ".xlsx"),
                    sheetName = type, row.names = FALSE, append = TRUE)
        gc()
        write.xlsx2(car_frequency_over_time, file = paste0(outputDir2, "car_clonotype_frequency_over_time_", patient, ".xlsx"),
                    sheetName = type, row.names = FALSE, append = TRUE)
        gc()
        write.xlsx2(car_proportion_over_time, file = paste0(outputDir2, "car_clonotype_proportion_over_time_", patient, ".xlsx"),
                    sheetName = type, row.names = FALSE, append = TRUE)
        gc()
        
        ###
        ### visualize the frequency and proportion of some interesting clonotypes with bar plots
        ###
        top10_carts <- t(car_frequency_over_time[1:10,2:(ncol(car_frequency_over_time)-1)])
        top10_carts <- data.frame(top10_carts,
                                  Time=rownames(top10_carts),
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
        fName <- paste0("Top10_CAR+_Tcells_Lineage_Tracing_", patient, "_", type)
        g <- arrangeGrob(grobs = p,
                         nrow = 5,
                         ncol = 2,
                         top = fName)
        ggsave(file = paste0(outputDir2, "/", fName, ".png"), g, width = 12, height = 12, dpi = 300)
        
        
        ###
        ### similar to the bar plot above but with the top GMP clonotypes
        ###
        frequency_over_time <- frequency_over_time[order(-as.numeric(frequency_over_time$GMP)),]
        top10_gmps <- t(frequency_over_time[1:10,2:(ncol(frequency_over_time)-1)])
        top10_gmps <- data.frame(top10_gmps,
                                 Time=rownames(top10_gmps),
                                 stringsAsFactors = FALSE, check.names = FALSE)
        top10_gmps$Time <- factor(top10_gmps$Time, levels = levels(Seurat_Obj@meta.data$TimeF))
        p <- vector("list", length = ncol(top10_gmps)-1)
        names(p) <- colnames(top10_gmps)[1:length(p)]
        for(clonotype in names(p)) {
          p[[clonotype]] <- ggplot(top10_gmps, aes_string(x="Time", y=clonotype)) +
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
        fName <- paste0("Top10_GMP_Clone_Size_Lineage_Tracing_", patient, "_", type)
        g <- arrangeGrob(grobs = p,
                         nrow = 5,
                         ncol = 2,
                         top = fName)
        ggsave(file = paste0(outputDir2, "/", fName, ".png"), g, width = 12, height = 12, dpi = 300)
        
        
        ###
        ### clone size of all the CAR+ T cells on UMAP
        ### CAR+ T cells which have TCR info
        ###
        
        ### One plot for each group
        p <- vector("list", length(unique_time_points))
        names(p) <- unique_time_points
        for(time in unique_time_points) {
          ### data creation
          temp <- Seurat_Obj@meta.data[intersect(which(Seurat_Obj@meta.data$Px == patient),
                                                 which(Seurat_Obj@meta.data$Time == time)),]
          temp <- merge(temp, car_frequency_over_time[,c("Clonotype", time),drop=FALSE],
                        by.x = type, by.y = "Clonotype",
                        all.x = TRUE)
          plot_df <- data.frame(X=temp$Full_UMAP1,
                                Y=temp$Full_UMAP2,
                                Time=temp$TimeF,
                                Clonotype=temp[,type],
                                Clone_Size=temp[,time])
          ### we only used dups, so
          ### give 0 value to the only-one-clonotypes across time points
          ### give 1 value if it is a CAR+ cell
          plot_df$Clone_Size[intersect(intersect(intersect(which(!is.na(plot_df$Clonotype)),
                                                           which(plot_df$Clonotype != "")),
                                                 which(is.na(plot_df$Clone_Size))),
                                       which(temp$CAR == "CARneg"))] <- 0
          plot_df$Clone_Size[intersect(intersect(intersect(which(!is.na(plot_df$Clonotype)),
                                                           which(plot_df$Clonotype != "")),
                                                 which(is.na(plot_df$Clone_Size))),
                                       which(temp$CAR == "CARpos"))] <- 1
          ### remove cells that do not have TCR info
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
        fName <- paste0("UMAP_CAR+_Tcells_Clone_Size_", patient, "_", type)
        g <- arrangeGrob(grobs = p,
                         nrow = 3,
                         top = fName)
        ggsave(file = paste0(outputDir2, "/", fName, ".png"), g, width = 20, height = 12, dpi = 300)
        
        
        ###
        ### the number of lineages
        ### Lineages: the number of lineages that appeared more than one time points in the Seurat object
        ### Lineages_CARpos: the number of lineages in CAR+ cells that appeared more than one time points in the Seurat object
        ### Clonotypes: the number of clonotypes in the Seurat object
        ### Clonotypes_CARpos: the number of clonotypes in CAR+ cells in the Seurat object
        ###
        
        ### Lineages
        chains <- unique(Seurat_Obj@meta.data[px_idx,type])
        chains <- chains[!is.na(chains)]
        chains <- chains[which(chains != "")]
        cnt <- 0
        for(i in 1:length(chains)) {
          if(length(unique(Seurat_Obj@meta.data$Time[intersect(px_idx,
                                                               which(Seurat_Obj@meta.data[,type] == chains[i]))])) > 1) {
            cnt <- cnt+1
          }
        }
        lineage_abstract_table["Lineages",type] <- cnt
        
        ### Lineages - CAR+
        cnt <- 0
        for(i in 1:length(chains)) {
          if(length(unique(Seurat_Obj@meta.data$Time[intersect(px_idx,
                                                               intersect(which(Seurat_Obj@meta.data[,type] == chains[i]),
                                                                         which(Seurat_Obj@meta.data$CAR == "CARpos")))])) > 1) {
            cnt <- cnt+1
          }
        }
        lineage_abstract_table["Lineages_CARpos",type] <- cnt
        
        ### Clonotypes
        lineage_abstract_table["Clonotypes",type] <- length(chains)
        
        ### Clonotypes - CAR+
        car_chains <- unique(Seurat_Obj@meta.data[intersect(px_idx,
                                                            which(Seurat_Obj@meta.data$CAR == "CARpos")),type])
        car_chains <- car_chains[!is.na(car_chains)]
        lineage_abstract_table["Clonotypes_CARpos",type] <- length(car_chains)
      }
      
      ### save the lineage table
      write.xlsx2(lineage_abstract_table,
                  file = paste0(outputDir2, "lineage_abstract_table_", patient, ".xlsx"),
                  sheetName = patient)
    
    }
    
  }
  
  ###
  ### the number of clonotypes that appeared in GMP and also in the later time points in CAR+ cells
  ###
  
  ### get directories
  f <- list.dirs(outputDir, full.names = FALSE, recursive = FALSE)
  f <- f[startsWith(f, "SJCAR19")]
  
  ### make an empty matrix for the result
  clonotype_gmp_all_patients <- matrix(0, length(global_clonotypes), length(f))
  rownames(clonotype_gmp_all_patients) <- global_clonotypes
  colnames(clonotype_gmp_all_patients) <- f
  
  ### for each patient
  for(i in 1:length(f)) {
    ### for each type of clonotyping approach
    for(type in global_clonotypes) {
      
      ### load the ingridient table
      car_frequency_over_time <- read.xlsx2(file = paste0(outputDir, f[i], "/car_clonotype_frequency_over_time_", f[i], ".xlsx"),
                                            sheetName = type, stringsAsFactors = FALSE, check.names = FALSE)
      rownames(car_frequency_over_time) <- car_frequency_over_time$Clonotype
      
      ### numerize the table
      for(j in 2:ncol(car_frequency_over_time)) {
        car_frequency_over_time[,j] <- as.numeric(car_frequency_over_time[,j])
      }
      
      ### get the number
      temp <- car_frequency_over_time[which(car_frequency_over_time$GMP > 0),
                                      setdiff(colnames(car_frequency_over_time),
                                              c("Clonotype", "PreTrans", "GMP", "Wk-1", "Wk0", "Total"))]
      temp <- apply(temp, 1, sum)
      clonotype_gmp_all_patients[type,i] <- length(which(temp > 0))
      
      ### garbage collection
      gc()
      
    }
  }
  
  ### save the table
  write.xlsx2(data.frame(clonotype_gmp_all_patients, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "clonotype_gmp_later_carpos.xlsx"),
              sheetName = "clonotype_gmp_all_patients",
              check.names = FALSE)
  
  ### make a graph
  p <- vector("list", length = length(global_clonotypes))
  names(p) <- global_clonotypes
  for(type in names(p)) {
    plot_df <- data.frame(t(clonotype_gmp_all_patients[type,,drop=FALSE]), check.names = FALSE)
    plot_df$Px <- rownames(plot_df)
    p[[type]] <- ggplot(plot_df, aes_string(x="Px", y=type)) +
      labs(x="", y="The Number of GMP-Post GMP CAR+ Cases") +
      geom_bar(fill = "grey50", position = "dodge", stat = "identity") +
      ggtitle(paste0(type)) +
      guides(fill=guide_legend(title=NULL)) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 14),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14),
            axis.title.y = element_text(size = 10))
  }
  
  ### arrange the plots and save
  fName <- paste0("clonotype_gmp_later_carpos")
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   top = fName)
  ggsave(file = paste0(outputDir, fName, ".png"), g, width = 15, height = 10, dpi = 300)
  
  
  ###
  ### the number of clonotypes and lineages in CAR+ cells of all the patients
  ###
  
  ### get directories
  f <- list.dirs(outputDir, full.names = FALSE, recursive = FALSE)
  
  ### make empty matrices for the results
  clonotype_all_patients_car <- matrix(0, length(global_clonotypes), length(f))
  rownames(clonotype_all_patients_car) <- global_clonotypes
  colnames(clonotype_all_patients_car) <- f
  lineage_all_patients_car <- matrix(0, length(global_clonotypes), length(f))
  rownames(lineage_all_patients_car) <- global_clonotypes
  colnames(lineage_all_patients_car) <- f
  
  ### for each patient
  for(i in 1:length(f)) {
    ### load the ingridient table
    lineage_abstract_table <- read.xlsx2(file = paste0(outputDir, f[i], "/lineage_abstract_table_", f[i], ".xlsx"),
                                         sheetIndex = 1, row.names = 1,
                                         stringsAsFactors = FALSE, check.names = FALSE)
    
    ### fill out the matrices
    for(type in global_clonotypes) {
      clonotype_all_patients_car[type, f[i]] <- as.numeric(lineage_abstract_table["Clonotypes_CARpos",type])
      lineage_all_patients_car[type,f[i]] <- as.numeric(lineage_abstract_table["Lineages_CARpos",type])
    }
  }
  
  ### save the tables
  write.xlsx2(data.frame(clonotype_all_patients_car, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "clonotype_carpos.xlsx"),
              sheetName = "clonotype_carpos",
              check.names = FALSE)
  write.xlsx2(data.frame(lineage_all_patients_car, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "lineage_carpos.xlsx"),
              sheetName = "lineage_carpos",
              check.names = FALSE)
  
  ### bar plot
  p <- vector("list", 2)
  plot_df <- data.frame(matrix(0, length(global_clonotypes)*length(f), 3))
  colnames(plot_df) <- c("Number", "Type", "Patient")
  plot_df[,"Number"] <- as.vector(clonotype_all_patients_car)
  plot_df[,"Type"] <- rep(global_clonotypes, length(f))
  plot_df[,"Patient"] <- as.vector(sapply(f, function(x) rep(x, 3)))
  p[[1]] <- ggplot(plot_df, aes_string(x="Patient", y="Number", fill="Type")) +
    labs(x="", y="The Number of Clonotypes in CAR+ Cells") +
    geom_bar(position = "dodge", stat = "identity") +
    ggtitle("The Number of Clonotypes in CAR+ Cells") +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    guides(fill=guide_legend(title=NULL)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 14),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14),
          axis.title.y = element_text(size = 10))
  plot_df <- data.frame(matrix(0, length(global_clonotypes)*length(f), 3))
  colnames(plot_df) <- c("Number", "Type", "Patient")
  plot_df[,"Number"] <- as.vector(lineage_all_patients_car)
  plot_df[,"Type"] <- rep(global_clonotypes, length(f))
  plot_df[,"Patient"] <- as.vector(sapply(f, function(x) rep(x, 3)))
  p[[2]] <- ggplot(plot_df, aes_string(x="Patient", y="Number", fill="Type")) +
    labs(x="", y="The Number of Lineages across CAR+ Cells") +
    geom_bar(position = "dodge", stat = "identity") +
    ggtitle("The Number of Lineages across CAR+ Cells") +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    guides(fill=guide_legend(title=NULL)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 14),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14),
          axis.title.y = element_text(size = 10))
  
  ### arrange the plots and save
  fName <- paste0("Clonotypes_and_Lineages_in_CARpos_Cells")
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   top = "")
  ggsave(file = paste0(outputDir, fName, ".png"), g, width = 20, height = 10, dpi = 300)
  
  
  ###
  ### clonotype diversity after GMP
  ###
  
  ### get target patients
  target_px <- unique(Seurat_Obj@meta.data$Px[intersect(intersect(intersect(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                                                                            which(Seurat_Obj@meta.data$Time != "GMP")),
                                                                  which(Seurat_Obj@meta.data$Time != "PreTrans")),
                                                        which(!is.na(Seurat_Obj@meta.data$cdr3_nt)))])
  target_px <- target_px[order(target_px)]
  
  ### get target time points
  target_tp <- unique(Seurat_Obj@meta.data$Time[which(Seurat_Obj@meta.data$Px %in% target_px)])
  target_tp <- factor(target_tp, levels = levels(Seurat_Obj@meta.data$TimeF))
  target_tp <- as.character(target_tp[order(target_tp)])
  target_tp <- setdiff(target_tp, c("PreTrans", "Wk-1", "Wk0"))
  
  ### for each patient, make a clonotype diversity table and draw a plot
  p <- vector("list", length = length(target_px))
  names(p) <- target_px
  car_p <- vector("list", length = length(target_px))
  names(car_p) <- target_px
  norm_p <- vector("list", length = length(target_px))
  names(norm_p) <- target_px
  car_norm_p <- vector("list", length = length(target_px))
  names(car_norm_p) <- target_px
  for(patient in target_px) {
    ### make an empty matrix for plot
    plot_df <- data.frame(matrix(NA, length(target_tp)*length(global_clonotypes), 3),
                          stringsAsFactors = FALSE, check.names = FALSE)
    colnames(plot_df) <- c("Number", "Time", "Type")
    car_plot_df <- data.frame(matrix(NA, length(target_tp)*length(global_clonotypes), 3),
                              stringsAsFactors = FALSE, check.names = FALSE)
    colnames(car_plot_df) <- c("Number", "Time", "Type")
    norm_plot_df <- data.frame(matrix(NA, length(target_tp)*length(global_clonotypes), 3),
                               stringsAsFactors = FALSE, check.names = FALSE)
    colnames(norm_plot_df) <- c("Number", "Time", "Type")
    car_norm_plot_df <- data.frame(matrix(NA, length(target_tp)*length(global_clonotypes), 3),
                                   stringsAsFactors = FALSE, check.names = FALSE)
    colnames(car_norm_plot_df) <- c("Number", "Time", "Type")
    
    ### for each time point, fill out the matrix
    for(i in 1:length(target_tp)) {
      ### for each clonotyping type
      for(j in 1:length(global_clonotypes)) {
        ### the number of cells in the patient that have TCR info
        nc_idx <- intersect(intersect(which(Seurat_Obj@meta.data$Px == patient),
                                      which(Seurat_Obj@meta.data$Time == target_tp[i])),
                            intersect(which(Seurat_Obj@meta.data[,global_clonotypes[j]] != ""),
                                      which(!is.na(Seurat_Obj@meta.data[,global_clonotypes[j]]))))
        nc <- length(nc_idx)
        
        ### the number of CAR+ cells in the patient that have TCR info
        ncc_idx <- intersect(nc_idx, which(Seurat_Obj@meta.data$CAR == "CARpos"))
        ncc <- length(ncc_idx)
        
        ### the number of clonotypes
        plot_df[(i-1)*3+j,"Number"] <- length(unique(Seurat_Obj@meta.data[nc_idx,
                                                                          global_clonotypes[j]]))
        car_plot_df[(i-1)*3+j,"Number"] <- length(unique(Seurat_Obj@meta.data[ncc_idx,
                                                                              global_clonotypes[j]]))
        norm_plot_df[(i-1)*3+j,"Number"] <- signif(plot_df[(i-1)*3+j,"Number"] / nc, 3)
        car_norm_plot_df[(i-1)*3+j,"Number"] <- signif(car_plot_df[(i-1)*3+j,"Number"] / ncc, 3)
      }
    }
    plot_df[,"Time"] <- as.vector(sapply(target_tp, function(x) rep(x, length(global_clonotypes))))
    plot_df[,"Type"] <- rep(global_clonotypes, length(target_tp))
    car_plot_df[,"Time"] <- as.vector(sapply(target_tp, function(x) rep(x, length(global_clonotypes))))
    car_plot_df[,"Type"] <- rep(global_clonotypes, length(target_tp))
    norm_plot_df[,"Time"] <- as.vector(sapply(target_tp, function(x) rep(x, length(global_clonotypes))))
    norm_plot_df[,"Type"] <- rep(global_clonotypes, length(target_tp))
    car_norm_plot_df[,"Time"] <- as.vector(sapply(target_tp, function(x) rep(x, length(global_clonotypes))))
    car_norm_plot_df[,"Type"] <- rep(global_clonotypes, length(target_tp))
    
    ### make 0 values to NA
    plot_df$Number[which(plot_df$Number == 0)] <- NA
    car_plot_df$Number[which(car_plot_df$Number == 0)] <- NA
    norm_plot_df$Number[which(norm_plot_df$Number == 0)] <- NA
    car_norm_plot_df$Number[which(car_norm_plot_df$Number == 0)] <- NA
    
    ### make "Time" factorized
    plot_df$Time <- factor(plot_df$Time, levels = levels(Seurat_Obj@meta.data$TimeF))
    car_plot_df$Time <- factor(car_plot_df$Time, levels = levels(Seurat_Obj@meta.data$TimeF))
    norm_plot_df$Time <- factor(norm_plot_df$Time, levels = levels(Seurat_Obj@meta.data$TimeF))
    car_norm_plot_df$Time <- factor(car_norm_plot_df$Time, levels = levels(Seurat_Obj@meta.data$TimeF))
    
    ### save the plot to the list
    p[[patient]] <- ggplot(plot_df, aes_string(x="Time", y="Number", fill="Type", group="Type")) +
      labs(x="", y="The Number of Clonotypes") +
      geom_bar(position = "dodge", stat = "identity") +
      geom_text(aes_string(label="Number", color="Type", group="Type"),
                position=position_dodge(width=1), size=3.5, hjust=0.5, vjust=-0.25,
                show.legend = FALSE) +
      ggtitle(patient) +
      scale_fill_npg() +
      scale_color_npg() +
      guides(fill=guide_legend(title=NULL)) +
      ylim(0, max(plot_df$Number, na.rm = TRUE)*1.1) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 14),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14),
            axis.title.y = element_text(size = 10), legend.position="top")
    car_p[[patient]] <- ggplot(car_plot_df, aes_string(x="Time", y="Number", fill="Type", group="Type")) +
      labs(x="", y="The Number of Clonotypes in CAR+") +
      geom_bar(position = "dodge", stat = "identity") +
      geom_text(aes_string(label="Number", color="Type", group="Type"),
                position=position_dodge(width=1), size=3.5, hjust=0.5, vjust=-0.25,
                show.legend = FALSE) +
      ggtitle(patient) +
      scale_fill_npg() +
      scale_color_npg() +
      guides(fill=guide_legend(title=NULL)) +
      ylim(0, max(car_plot_df$Number, na.rm = TRUE)*1.1) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 14),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14),
            axis.title.y = element_text(size = 10), legend.position="top")
    norm_p[[patient]] <- ggplot(norm_plot_df, aes_string(x="Time", y="Number", fill="Type", group="Type")) +
      labs(x="", y="Clonotypes Norm.#") +
      geom_bar(position = "dodge", stat = "identity") +
      geom_text(aes_string(label="Number", color="Type", group="Type"),
                position=position_dodge(width=1), size=3.5, hjust=0.5, vjust=-0.25,
                show.legend = FALSE) +
      ggtitle(patient) +
      scale_fill_npg() +
      scale_color_npg() +
      guides(fill=guide_legend(title=NULL)) +
      ylim(0, max(norm_plot_df$Number, na.rm = TRUE)*1.1) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 14),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14),
            axis.title.y = element_text(size = 10), legend.position="top")
    car_norm_p[[patient]] <- ggplot(car_norm_plot_df, aes_string(x="Time", y="Number", fill="Type", group="Type")) +
      labs(x="", y="Clonotypes Norm.# in CAR+") +
      geom_bar(position = "dodge", stat = "identity") +
      geom_text(aes_string(label="Number", color="Type", group="Type"),
                position=position_dodge(width=1), size=3.5, hjust=0.5, vjust=-0.25,
                show.legend = FALSE) +
      ggtitle(patient) +
      scale_fill_npg() +
      scale_color_npg() +
      guides(fill=guide_legend(title=NULL)) +
      ylim(0, max(car_norm_plot_df$Number, na.rm = TRUE)*1.1) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 14),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 14),
            axis.title.y = element_text(size = 10), legend.position="top")
  }
  
  ### arrange the plots and save
  fName <- "Clonotype_Diversity_After_GMP"
  g <- arrangeGrob(grobs = p,
                   nrow = 4,
                   ncol = 2,
                   top = fName)
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 24, height = 12, dpi = 300)
  fName <- "Clonotype_Diversity_After_GMP_In_CAR+"
  g <- arrangeGrob(grobs = car_p,
                   nrow = 4,
                   ncol = 2,
                   top = fName)
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 24, height = 12, dpi = 300)
  fName <- "Clonotype_Diversity_After_GMP_Norm"
  g <- arrangeGrob(grobs = norm_p,
                   nrow = 4,
                   ncol = 2,
                   top = fName)
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 24, height = 12, dpi = 300)
  fName <- "Clonotype_Diversity_After_GMP_In_CAR+_Norm"
  g <- arrangeGrob(grobs = car_norm_p,
                   nrow = 4,
                   ncol = 2,
                   top = fName)
  ggsave(file = paste0(outputDir, "/", fName, ".png"), g, width = 24, height = 12, dpi = 300)
  
}
