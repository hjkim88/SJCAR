###
#   File name : Final_Figures_And_Tables.R
#   Author    : Hyunjin Kim
#   Date      : Apr 25, 2022
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : We submitted a manuscript to Cancer Discovery.
#               This code generates all the figures and tables of the manuscript.
#
#   Instruction
#               1. Source("Final_Figures_And_Tables.R")
#               2. Run the function "generate_final" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Final_Figures_And_Tables.R/Final_Figures_And_Tables.R")
#               > generate_final(Seurat_RObj_path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds",
#                                outputDir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19/Final/")
###

generate_final <- function(Seurat_RObj_path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR19_data/data/NEW_SJCAR_SEURAT_OBJ/CARpos_JCC.rds",
                           outputDir="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/SJCAR/Final/") {
  
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
  if(!require(ggforce, quietly = TRUE)) {
    install.packages("ggforce")
    require(ggforce, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  if(!require(viridis, quietly = TRUE)) {
    install.packages("viridis")
    require(viridis, quietly = TRUE)
  }
  if(!require(grid, quietly = TRUE)) {
    install.packages("grid")
    require(grid, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    require(scales, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(ggrepel, quietly = TRUE)) {
    install.packages("ggrepel")
    require(ggrepel, quietly = TRUE)
  }
  if(!require(monocle, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("monocle")
    require(monocle, quietly = TRUE)
  }
  if(!require(ggthemes, quietly = TRUE)) {
    install.packages("ggthemes")
    require(ggthemes, quietly = TRUE)
  }
  if(!require(caret, quietly = TRUE)) {
    install.packages("caret")
    require(caret, quietly = TRUE)
  }
  if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
  }
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    require(gplots, quietly = TRUE)
  }
  if(!require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    require(RColorBrewer, quietly = TRUE)
  }
  if(!require(wesanderson, quietly = TRUE)) {
    install.packages("wesanderson")
    require(wesanderson, quietly = TRUE)
  }
  if(!require(dplyr, quietly = TRUE)) {
    install.packages("dplyr")
    require(dplyr, quietly = TRUE)
  }
  if(!require(OmicCircos, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("OmicCircos")
    require(OmicCircos, quietly = TRUE)
  }
  if(!require(igraph, quietly = TRUE)) {
    install.packages("igraph")
    require(igraph, quietly = TRUE)
  }
  if(!require(RedeR, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("RedeR")
    require(RedeR, quietly = TRUE)
  }
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(shadowtext, quietly = TRUE)) {
    install.packages("shadowtext")
    require(shadowtext, quietly = TRUE)
  }
  
  ### load Jeremy's object
  JCC_Seurat_Obj <- readRDS(file = Seurat_RObj_path)
  
  ### check whether the orders are the same
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
  print(identical(names(Idents(object = JCC_Seurat_Obj)), rownames(JCC_Seurat_Obj@meta.data)))
  
  ### create outputDir
  dir.create(outputDir, showWarnings = TRUE, recursive = TRUE)
  
  
  ### Fig1B
  sjcar19_colors <- c("#16101E", "#D0B78F", "#8C8781", "#C7904F", "#133D31", "#82A5B8", "#3B3B53", "#4C8493", "#C31517",
                      "#D94C21", "#3E89A8", "#AA4C26",  "#CAA638", "#640B11", "#629488", "#BD7897", "#3C2D16", "#25245D",
                      "#E64E46", "#73BCAA", "#7047C1", "#286278")
  names(sjcar19_colors) <- levels(JCC_Seurat_Obj$AllSeuratClusters)
  show_col(sjcar19_colors)
  
  p <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "AllSeuratClusters",
               pt.size = 1,
               label = TRUE,
               label.color = "cornsilk2",
               label.size = 10,
               cols = c(sjcar19_colors)) +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 40) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40, color = "black", face = "bold"),
          axis.title = element_text(size = 25, color = "black", face = "bold"),
          axis.text = element_text(size = 20, color = "black", face = "bold"),
          legend.position = "none")
  p <- p + geom_shadowtext(data = p$layers[[2]]$data, aes(x = UMAP_1, y = UMAP_2, label=AllSeuratClusters),
                           size=10, color="cornsilk2", bg.color="black", bg.r=0.2)
  p$layers[[2]] <- NULL
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir, "Fig1B.pdf"), plot = p, width = 12, height = 10, dpi = 350)
  
  
  ### Fig1D
  interesting_genes <- c("RPL32", "RPL30", "LAG3", "TOX", "CASP8", "IL7R", "SELL", "BNIP3", "MKI67",
                         "CDC20", "CDK1", "NKG7", "GNLY", "GZMH", "GZMM", "GZMK")
  
  ### NEW FUNCTIONAL ANNOTATION BY TAY - 09/27/2021
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters <- "Others"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "1", "5", "7", "9", "10", "11", "12", "15", "19", "21"))] <- "Proliferating"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("2", "4", "6", "17"))] <- "Transitioning"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("3", "8", "14"))] <- "Functional Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("13", "20"))] <- "Dysfunctional"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("16"))] <- "Early Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("18"))] <- "Metabolically Active"
  
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2 <- JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Proliferating")] <- "Clusters{0,1,5,7,9,10,11,12,15,19,21}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Transitioning")] <- "Clusters{2,4,6,17}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Metabolically Active")] <- "Clusters{18}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Early Effector")] <- "Clusters{16}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Dysfunctional")] <- "Clusters{13, 20}"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters2[which(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters == "Functional Effector")] <- "Clusters{3,8,10}"
  
  ### https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=11
  p <- DotPlot(JCC_Seurat_Obj,
               features = interesting_genes,
               group.by = "New_Functional_Annotation_Based_On_Clusters2") +
    scale_size(range = c(5, 45)) +
    xlab("") +
    ylab("") +
    scale_color_gradientn(colours = c("#313695", "#ffffbf", "#a50026")) +
    guides(color = guide_colorbar(title = "Relative Expression")) +
    coord_flip() +
    theme_classic(base_size = 28) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = -90, size = 35, vjust = 0.5, hjust = 0, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 25, color = "black", face = "bold"),
          legend.text = element_text(size = 20, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'),
          legend.position = "top")
  ggsave(file = paste0(outputDir, "Fig1D.pdf"),
         plot = p, width = 15, height = 35, dpi = 350)
  
  
  ### Fig2B
  ### gray to sjcar19 blue color
  blue_color_scale <- colorRampPalette(c("gray85", "#08519c"))(length(unique(JCC_Seurat_Obj$time2)))
  names(blue_color_scale) <- unique(JCC_Seurat_Obj$time2)
  blue_color_scale["6mo"] <- "#08306b"
  show_col(blue_color_scale)
  
  ### Make GMP cells gray; color PI cells by their time points
  p <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
               group.by = "time2",
               pt.size = 2.5, raster = TRUE,
               cols = blue_color_scale[JCC_Seurat_Obj$time2],
               order = rev(unique(JCC_Seurat_Obj$time2))) +
    ggtitle("") +
    labs(color="Time") +
    theme_classic(base_size = 36) +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.9
  ggsave(paste0(outputDir, "Fig2B.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### Fig2C
  ### make a plot matrix
  target_genes <- c("PRF1", "GZMB", "GZMM", "GZMH", "GZMK", "GZMA", "TBX21", "KLRD1",
                    "KLRG1", "GNLY", "EOMES", "TOX", "PDCD1", "LAG3", "TIGIT", "CASP8")
  target_tp <- unique(JCC_Seurat_Obj$time2)
  target_tp <- target_tp[-which(target_tp == "Wk6")]
  heatmap_mat <- matrix(0, length(target_genes), length(target_tp))
  rownames(heatmap_mat) <- target_genes
  colnames(heatmap_mat) <- target_tp
  
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@data)))
  
  for(g in target_genes) {
    for(tp in target_tp) {
      heatmap_mat[g,tp] <- mean(as.numeric(JCC_Seurat_Obj@assays$RNA@data[g,which(JCC_Seurat_Obj$time2 == tp)]))
    }
  }
  
  ### set the custom distance and clustering functions
  hclustfunc <- function(x) hclust(x, method="complete")
  distfunc <- function(x) dist(x, method="euclidean")
  
  ### add cell # in the col names
  adv_col_names <- colnames(heatmap_mat)
  cell_num <- c("(118,749)", "(21,085)", "(17,612)", "(19,872)", "(5,562)", "(1,535)", "(368)", "(7)")
  adv_col_names <- paste(adv_col_names, cell_num, sep = "\n")
  
  ### draw a heatmap
  pdf(paste0(outputDir, "Fig2C.pdf"),
      width = 30, height = 25)
  par(oma=c(10,0,3,5), xpd = TRUE, cex.main = 4, font = 2, font.lab = 2, font.axis = 2)
  heatmap.2(heatmap_mat, col = colorpanel(24, low = "#313695", mid = "#ffffbf", high = "#a50026"),
            scale = "row", dendrogram = "none", trace = "none",
            Rowv = FALSE, Colv = FALSE,
            cexRow = 6, cexCol = 5,
            key.title = "", keysize = 2, density.info = "none",
            key.xlab = "Z-Score", key.ylab = "", key.par = list(mar=c(6,0,0,0), xpd=TRUE, cex.lab=2, cex.axis=1.5, cex=2, font=1, font.axis=1),
            # main = "Average Gene Expression Across Post-Infusion",
            hclustfun = hclustfunc, distfun = distfunc,
            lmat = rbind(3:4,2:1),
            labRow = rownames(heatmap_mat), labCol = adv_col_names,
            offsetRow=-148, offsetCol = 10, lwid = c(1.5, 8), lhei = c(1.5, 8),
            adjRow = c(1, 0.5), srtCol = 0, adjCol = c(0.5, 0.5))
  dev.off()
  
  
  ### Fig2D
  ### prepare table for the plot
  plot_df <- data.frame(Cluster=as.character(sapply(levels(JCC_Seurat_Obj$AllSeuratClusters), function(x) rep(x, length(unique(JCC_Seurat_Obj$time2))))),
                        Time=as.character(rep(unique(JCC_Seurat_Obj$time2), length(levels(JCC_Seurat_Obj$AllSeuratClusters)))),
                        Numbers=0,
                        Pcnt=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### calculate numbers
  cnt <- 1
  for(clstr in levels(JCC_Seurat_Obj$AllSeuratClusters)) {
    for(tp in unique(JCC_Seurat_Obj$time2)) {
      plot_df$Numbers[cnt] <- length(intersect(which(JCC_Seurat_Obj$AllSeuratClusters == clstr),
                                               which(JCC_Seurat_Obj$time2 == tp)))
      cnt <- cnt + 1
    }
  }
  
  ### order the plot_df in a GMP/PI - oriented way
  plot_df <- plot_df[order(plot_df$Time),]
  
  ### calculate percentages
  time_sum <- rep(0, length(unique(plot_df$Time)))
  names(time_sum) <- unique(plot_df$Time)
  for(i in 1:length(unique(plot_df$Time))) {
    time_sum[i] <- sum(plot_df[which(plot_df$Time == unique(plot_df$Time)[i]),"Numbers"])
    plot_df$Pcnt[which(plot_df$Time == unique(plot_df$Time)[i])] <- round(plot_df$Numbers[which(plot_df$Time == unique(plot_df$Time)[i])] * 100 / time_sum[i], 1)
  }
  
  ### remove NaN rows
  plot_df <- plot_df[which(!is.nan(plot_df$Pcnt)),]
  
  ### preserve the pcnt
  plot_df$Pcnt2 <- plot_df$Pcnt
  
  ### pcnt < 10 -> ""
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) < 10)] <- ""
  plot_df$Pcnt <- as.character(plot_df$Pcnt)
  
  ### add cluster2 column
  plot_df$Cluster2 <- ""
  plot_df$Cluster2[which(as.numeric(plot_df$Numbers) != 0)] <- as.character(plot_df$Cluster[which(as.numeric(plot_df$Numbers) != 0)])
  plot_df$Cluster2[which(plot_df$Pcnt == "")] <- ""
  
  ### annotate "%"
  plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) != 0)] <- paste0(plot_df$Pcnt[which(as.numeric(plot_df$Pcnt) != 0)], "%")
  plot_df$Pcnt[which(plot_df$Pcnt == 0)] <- ""
  
  ### factorize the time point & state
  plot_df$Time <- factor(plot_df$Time, levels = unique(JCC_Seurat_Obj$time2))
  plot_df$Cluster <- factor(plot_df$Cluster, levels = unique(plot_df$Cluster))
  
  ### remove the Wk6
  plot_df <- plot_df[which(!plot_df$Time %in% c("Wk6")),]
  
  ### color scale
  sjcar19_colors <- c("#16101E", "#D0B78F", "#8C8781", "#C7904F", "#133D31", "#82A5B8", "#3B3B53", "#4C8493", "#C31517",
                      "#D94C21", "#3E89A8", "#AA4C26",  "#CAA638", "#640B11", "#629488", "#BD7897", "#3C2D16", "#25245D",
                      "#E64E46", "#73BCAA", "#7047C1", "#286278")
  names(sjcar19_colors) <- unique(plot_df$Cluster)
  show_col(sjcar19_colors)
  
  ### draw a proportional bar plot
  p <- ggplot(data=plot_df, aes_string(x="Time", y="Pcnt2", fill="Cluster", label="Pcnt")) +
    geom_bar(position = "stack", stat = "identity") +
    # ggtitle("Proportion of Cells") +
    xlab("") + ylab("Percentage") +
    # geom_text(size = 3, position = position_stack(vjust = 0.5), vjust = 3, color = "blue") +
    geom_text(aes_string(x="Time", y="Pcnt2", label = "Cluster2"),
              position = position_stack(vjust = 0.5),
              size = 15, color = "white") +
    coord_flip() +
    scale_fill_manual(values = sjcar19_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(size = 35, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir, "Fig2D.pdf"), plot = p,
         width = 20, height = 10, dpi = 350)
  
  
  ### Fig3A
  ### cell # for each time point
  cell_num <- sapply(unique(JCC_Seurat_Obj@meta.data$time2), function(x) {
    return(length(which(JCC_Seurat_Obj@meta.data$time2 == x)))
  })
  
  ### remove wk6 & 6mo data since there are only 1 & 7 cells
  cell_num <- cell_num[-c(6,9)]
  
  ### the data is too big (can't run for monocle), so we down-sample the cells in each time point
  set.seed(1234)
  fixed_min_cell_num <- min(cell_num)
  JCC_Seurat_Obj@meta.data$downsampled <- "NO"
  for(tp in names(cell_num)) {
    JCC_Seurat_Obj@meta.data$downsampled[sample(which(JCC_Seurat_Obj@meta.data$time2 == tp), fixed_min_cell_num)] <- "YES"
  }
  
  ### the data is too big (can't run for monocle), so we down-sample the cells in each time point
  ### BUT THIS TIME, INCLUDE ALL THE SUBSISTERS (IN ADDITION TO THE ORIGINAL DOWNSAMPLED ONES)
  JCC_Seurat_Obj@meta.data$downsampled2 <- JCC_Seurat_Obj@meta.data$downsampled
  JCC_Seurat_Obj@meta.data$downsampled2[which(JCC_Seurat_Obj@meta.data$ALL_CARpos_Persister == "YES")] <- "YES"
  
  ### Construct a monocle cds
  monocle_metadata <- JCC_Seurat_Obj@meta.data[rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj@meta.data$downsampled2 == "YES")],]
  monocle_metadata$time2 <- factor(monocle_metadata$time2, levels = unique(monocle_metadata$time2))
  monocle_cds2 <- newCellDataSet(as(as.matrix(JCC_Seurat_Obj@assays$RNA@data[,rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj@meta.data$downsampled2 == "YES")]]), 'sparseMatrix'),
                                 phenoData = new('AnnotatedDataFrame', data = monocle_metadata),
                                 featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(JCC_Seurat_Obj@assays$RNA@data),
                                                                                           row.names = row.names(JCC_Seurat_Obj@assays$RNA@data),
                                                                                           stringsAsFactors = FALSE, check.names = FALSE)),
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())
  
  ### run monocle
  monocle_cds2 <- estimateSizeFactors(monocle_cds2)
  monocle_cds2 <- estimateDispersions(monocle_cds2)
  monocle_cds2 <- reduceDimension(monocle_cds2, reduction_method = "DDRTree")
  monocle_cds2 <- orderCells(monocle_cds2)
  
  ### determine the beginning state
  plot_cell_trajectory(monocle_cds2, color_by = "State")
  ### this root state should be checked if we want to rerun this -> Sometimes it's 6 and sometimes it's 1
  monocle_cds2 <- orderCells(monocle_cds2, root_state = "6")
  
  ### add previous monocle_state to the current one
  monocle_cds2$Original_State <- "NEW"
  monocle_cds2@phenoData@data[rownames(monocle_cds@phenoData@data),"Original_State"] <- monocle_cds$State
  
  ### change the state numbers to alphabets
  monocle_cds2$State <- as.character(monocle_cds2$State)
  monocle_cds2$State[which(monocle_cds2$State == "6")] <- "A"
  monocle_cds2$State[which(monocle_cds2$State == "5")] <- "B"
  monocle_cds2$State[which(monocle_cds2$State == "4")] <- "C"
  monocle_cds2$State[which(monocle_cds2$State == "7")] <- "D"
  monocle_cds2$State[which(monocle_cds2$State == "3")] <- "E"
  monocle_cds2$State[which(monocle_cds2$State == "1")] <- "F"
  monocle_cds2$State[which(monocle_cds2$State == "2")] <- "G"
  monocle_cds2$State <- factor(monocle_cds2$State, levels = c("A", "B", "C", "D", "E", "F", "G"))
  
  ### color scale
  sjcar19_colors <- c("#d73027", "#f46d43", "#fdae61", "#fee090", "#abd9e9", "#74add1", "#4575b4")
  names(sjcar19_colors) <- levels(monocle_cds2$State)
  
  p <- plot_cell_trajectory(monocle_cds2, color_by = "State",
                            cell_size = 3, cell_link_size = 2,
                            show_branch_points = FALSE) +
    labs(color="State") +
    scale_color_manual(values = sjcar19_colors, name = "State") +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    theme_classic(base_size = 30) +
    theme(legend.position = "top",
          axis.title = element_text(size = 30, color = "black", face = "bold"),
          axis.text = element_text(size = 25, color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"))
  ggsave(file = paste0(outputDir, "Fig3A.pdf"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  
  ### Fig3B
  ### represent GMP and PI on the pseudotime plot w/o individual time points
  two_color_scale <- c("#d7191c", "#2c7bb6")
  names(two_color_scale) <- c("GMP", "PI")
  
  p <- plot_cell_trajectory(monocle_cds2, color_by = "GMP", cell_size = 3, cell_link_size = 2, show_branch_points = FALSE) +
    labs(color="") +
    scale_color_manual(labels = c("GMP", "Post-Infusion"), values = two_color_scale, name = "") +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    theme_classic(base_size = 36) +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0, vjust = 0.5, size = 30, color = "black", face = "bold"),
          plot.subtitle = element_text(hjust = 0, vjust = 0.5, size = 25, color = "black", face = "bold"),
          axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
  ggsave(file = paste0(outputDir, "Fig3B.pdf"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  
  ### Fig3C
  ### make the down-sampled seurat object
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj@meta.data$downsampled2)
  downsampled_seurat_obj <- subset(JCC_Seurat_Obj,
                                   idents = "YES")
  downsampled_seurat_obj$monocle_state <- monocle_cds2@phenoData@data[rownames(downsampled_seurat_obj@meta.data),"State"]
  
  print(identical(rownames(downsampled_seurat_obj@meta.data), colnames(downsampled_seurat_obj@assays$RNA@counts)))
  print(identical(names(Idents(object = downsampled_seurat_obj)), rownames(downsampled_seurat_obj@meta.data)))
  
  ### set genes of interest
  interesting_genes <- c("CASP8", "LAG3", "TOX")
  
  downsampled_seurat_obj$monocle_state <- factor(downsampled_seurat_obj$monocle_state, levels = rev(levels(downsampled_seurat_obj$monocle_state)))
  p <- DotPlot(downsampled_seurat_obj,
               features = interesting_genes,
               group.by = "monocle_state") +
    scale_size(range = c(5, 35)) +
    xlab("") +
    ylab("State") +
    scale_color_gradientn(colours = c("#313695", "#ffffbf", "#a50026"),
                          n.breaks = 3) +
    theme_classic(base_size = 28) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(hjust = 0.5, size = 35, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'))
  ggsave(file = paste0(outputDir, "Fig3C.pdf"),
         plot = p, width = 18, height = 10, dpi = 350)
  
  
  ### Fig3D
  ### make the down-sampled seurat object
  JCC_Seurat_Obj <- SetIdent(object = JCC_Seurat_Obj,
                             cells = rownames(JCC_Seurat_Obj@meta.data),
                             value = JCC_Seurat_Obj@meta.data$downsampled2)
  downsampled_seurat_obj <- subset(JCC_Seurat_Obj,
                                   idents = "YES")
  downsampled_seurat_obj$monocle_state <- monocle_cds2@phenoData@data[rownames(downsampled_seurat_obj@meta.data),"State"]
  
  print(identical(rownames(downsampled_seurat_obj@meta.data), colnames(downsampled_seurat_obj@assays$RNA@counts)))
  print(identical(names(Idents(object = downsampled_seurat_obj)), rownames(downsampled_seurat_obj@meta.data)))
  
  ### set genes of interest
  interesting_genes <- c("NKG7", "GNLY", "GZMB", "GZMK")
  
  p <- DotPlot(downsampled_seurat_obj,
               features = interesting_genes,
               group.by = "monocle_state") +
    scale_size(range = c(5, 35)) +
    xlab("") +
    ylab("State") +
    scale_color_gradientn(colours = c("#313695", "#ffffbf", "#a50026"),
                          n.breaks = 5) +
    theme_classic(base_size = 28) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(hjust = 0.5, size = 35, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'))
  ggsave(file = paste0(outputDir, "Fig3D.pdf"),
         plot = p, width = 18, height = 10, dpi = 350)
  
  
  ### Fig4A
  ### load clonotype frequency data
  px_result_dir="./results/New3/"
  
  ### set time points
  total_time_points <- c("PreTrans", "Wk-1", "Wk0", "GMP", "Wk1", "Wk2", 
                         "Wk3", "Wk4", "Wk6", "Wk8", "3mo", "6mo", "9mo")
  gmp_after_time_points <- c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", 
                             "Wk6", "Wk8", "3mo", "6mo", "9mo")
  
  ### run for each patient
  total_plot_df <- NULL
  p <- vector("list", length= length(unique(JCC_Seurat_Obj@meta.data$px)))
  names(p) <- unique(JCC_Seurat_Obj@meta.data$px)
  for(patient in unique(JCC_Seurat_Obj@meta.data$px)) {
    
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
  sjcar19_color_scale <- colorRampPalette(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#abd9e9", "#74add1", "#4575b4", "#313695"))(length(unique(total_plot_df$Clone)))
  total_plot_df$Time <- factor(total_plot_df$Time, levels = intersect(total_time_points, unique(total_plot_df$Time)))
  p <- ggplot(total_plot_df,
              aes(x = Time, stratum = Clone, alluvium = Clone,
                  y = Clone_Size,
                  fill = Clone, label = Clone)) +
    ggtitle(paste("")) +
    geom_flow() +
    geom_stratum(alpha = 1) +
    # geom_text(stat = "stratum", size = 2) +
    rotate_x_text(90) +
    theme_pubr(legend = "none") +
    # theme_cleveland2() +
    scale_fill_manual(values = sjcar19_color_scale) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(size = 35, color = "black", face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'),
          legend.position = "none")
  ggsave(file = paste0(outputDir, "Fig4A.pdf"), plot = p, width = 20, height = 10, dpi = 350)
  
  
  ### Fig4B_1
  ### set some annotations
  JCC_Seurat_Obj$ALL_CARpos_Persister2 <- JCC_Seurat_Obj$ALL_CARpos_Persister
  JCC_Seurat_Obj$ALL_CARpos_Persister2[which(is.na(JCC_Seurat_Obj$ALL_CARpos_Persister2))] <- "NO"
  
  JCC_Seurat_Obj$Seurat_Clusters_Subsisters2 <- JCC_Seurat_Obj$ALL_CARpos_Persister2
  JCC_Seurat_Obj$Seurat_Clusters_Subsisters2[which(JCC_Seurat_Obj$Seurat_Clusters_Subsisters2 == "YES")] <- JCC_Seurat_Obj$time2[which(JCC_Seurat_Obj$Seurat_Clusters_Subsisters2 == "YES")]
  JCC_Seurat_Obj$Seurat_Clusters_Subsisters2[which(JCC_Seurat_Obj$Seurat_Clusters_Subsisters2 == "NO")] <- "Non-Subsisters"
  
  JCC_Seurat_Obj$Seurat_Clusters_Subsisters3 <- JCC_Seurat_Obj$Seurat_Clusters_Subsisters2
  JCC_Seurat_Obj$Seurat_Clusters_Subsisters3[which(JCC_Seurat_Obj$Seurat_Clusters_Subsisters3 %in% c("Non-Subsisters", "GMP"))] <- "NA"
  
  ### contruct arrow plot df - remove 'FROM GMP' arrows
  gmp_subsisters_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                                                  which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  pi_subsister_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                                which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  print(identical(gmp_subsisters_clones[order(gmp_subsisters_clones)], pi_subsister_clones[order(pi_subsister_clones)]))
  
  umap_map <- Embeddings(JCC_Seurat_Obj, reduction = "umap")[rownames(JCC_Seurat_Obj@meta.data), 1:2]
  print(identical(rownames(umap_map), rownames(JCC_Seurat_Obj@meta.data)))
  
  arrow_df <- data.frame(x1=0,
                         y1=0,
                         x2=0,
                         y2=0,
                         stringsAsFactors = FALSE, check.names = FALSE)
  cnt <- 1
  for(i in 1:length(gmp_subsisters_clones)) {
    target_indicies <- which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta == gmp_subsisters_clones[i])
    existing_time <- unique(JCC_Seurat_Obj$time2[target_indicies])
    for(j in 1:(length(existing_time)-1)) {
      if(existing_time[j] != "GMP") {
        target_indicies2 <- intersect(target_indicies,
                                      which(JCC_Seurat_Obj$time2 == existing_time[j]))
        arrow_df <- rbind(arrow_df, c(0, 0, 0, 0))
        arrow_df$x1[cnt] <- mean(umap_map[target_indicies2,
                                          "UMAP_1"])
        arrow_df$y1[cnt] <- mean(umap_map[target_indicies2,
                                          "UMAP_2"])
        target_indicies2 <- intersect(target_indicies,
                                      which(JCC_Seurat_Obj$time2 == existing_time[j+1]))
        arrow_df$x2[cnt] <- mean(umap_map[target_indicies2,
                                          "UMAP_1"])
        arrow_df$y2[cnt] <- mean(umap_map[target_indicies2,
                                          "UMAP_2"])
        cnt <- cnt + 1
      }
    }
  }
  arrow_df <- arrow_df[-nrow(arrow_df),]
  
  ### color scale
  sjcar19_colors <- c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6", "lightgray")
  names(sjcar19_colors) <- c("Wk1", "Wk2", "Wk3", "Wk4", "Wk8", "NA")
  
  ### add arrows to the previous UMAP - no GMP
  p <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "Seurat_Clusters_Subsisters3",
               pt.size = 5,
               cols = sjcar19_colors,
               order = rev(unique(JCC_Seurat_Obj$Seurat_Clusters_Subsisters3))) +
    ggtitle("") +
    labs(color="") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_text(size = 48),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 24)) +
    geom_curve(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = arrow_df,
      arrow = arrow(length = unit(0.03, "npc"))
    )
  ggsave(paste0(outputDir, "Fig4B_1.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### Fig4B_2
  ### GMP clusters vs PI clusters - how subsister lineages are moving
  gmp_subsister_idx <- intersect(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% gmp_subsisters_clones),
                                 which(JCC_Seurat_Obj$GMP == "GMP"))
  pi_subsister_idx <- intersect(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% pi_subsister_clones),
                                which(JCC_Seurat_Obj$GMP == "PI"))
  plot_df <- data.frame(GMP_PI=c(rep("GMP", length(gmp_subsister_idx)), rep("PI", length(pi_subsister_idx))),
                        Cluster=c(JCC_Seurat_Obj@meta.data[gmp_subsister_idx, "AllSeuratClusters"],
                                  JCC_Seurat_Obj@meta.data[pi_subsister_idx, "AllSeuratClusters"]),
                        Clone=c(JCC_Seurat_Obj@meta.data[gmp_subsister_idx, "clonotype_id_by_patient_one_alpha_beta"],
                                JCC_Seurat_Obj@meta.data[pi_subsister_idx, "clonotype_id_by_patient_one_alpha_beta"]),
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df$Clone_Cluster <- paste0(plot_df$Clone, "_", plot_df$Cluster)
  plot_df$Clone_Cluster_Size <- 1
  ### GMP
  dup_idx <- which(plot_df$GMP_PI == "GMP")[which(duplicated(plot_df$Clone_Cluster[which(plot_df$GMP_PI == "GMP")]))]
  for(i in dup_idx) {
    target_idx <- intersect(which(plot_df$Clone_Cluster == plot_df$Clone_Cluster[i]),
                            which(plot_df$GMP_PI == "GMP"))[1]
    plot_df$Clone_Cluster_Size[target_idx] <- plot_df$Clone_Cluster_Size[target_idx] + 1
  }
  plot_df <- plot_df[-dup_idx,]
  ### PI
  dup_idx <- which(plot_df$GMP_PI == "PI")[which(duplicated(plot_df$Clone_Cluster[which(plot_df$GMP_PI == "PI")]))]
  for(i in dup_idx) {
    target_idx <- intersect(which(plot_df$Clone_Cluster == plot_df$Clone_Cluster[i]),
                            which(plot_df$GMP_PI == "PI"))[1]
    plot_df$Clone_Cluster_Size[target_idx] <- plot_df$Clone_Cluster_Size[target_idx] + 1
  }
  plot_df <- plot_df[-dup_idx,]
  
  ### set parameters for circos
  plot_df$GMP_PI_Cluster <- paste0(plot_df$GMP_PI, "_", plot_df$Cluster)
  seg.name <- unique(plot_df$GMP_PI_Cluster)
  seg.name <- factor(seg.name, levels = c(paste0("GMP_", levels(JCC_Seurat_Obj$AllSeuratClusters)),
                                          paste0("PI_", levels(JCC_Seurat_Obj$AllSeuratClusters))))
  seg.name <- as.character(seg.name[order(seg.name)])
  sample.num <- sum(plot_df$Clone_Cluster_Size)
  
  ### set seg.f
  seg.f <- matrix("NA", sample.num, 9)
  colnames(seg.f) <- c("seg.name", "seg.start", "seg.end", "optional_col1", "optional_col2", "from_to", "clone", "time", "barcode")
  
  seg.f[,"seg.name"] <- unlist(sapply(seg.name, function(x) rep(x, sum(plot_df$Clone_Cluster_Size[which(plot_df$GMP_PI_Cluster == x)]))))
  
  for(sn in unique(seg.f[,"seg.name"])) {
    target_idx <- which(seg.f[,"seg.name"] == sn)
    for(i in 1:length(target_idx)) {
      seg.f[target_idx[i],"seg.start"] <- i-1
      seg.f[target_idx[i],"seg.end"] <- i
    }
  }
  seg.f <- data.frame(seg.f)
  
  ### set seg.v
  seg.v <- matrix(0, sample.num, 5)
  colnames(seg.v) <- c("seg.name", "sample", "exp_TIGIT", "exp_SELL", "exp_CD27")
  
  seg.v[,"seg.name"] <- seg.f[,"seg.name"]
  seg.v[,"sample"] <- seg.f[,"seg.end"]
  
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@data)))
  
  total_link_num <- 0
  for(sn in unique(seg.v[,"seg.name"])) {
    temp <- strsplit(sn, split = "_", fixed = TRUE)[[1]]
    target_tp <- temp[1]
    target_cluster <- temp[2]
    target_clones <- plot_df$Clone[intersect(which(plot_df$GMP_PI == temp[1]),
                                             which(plot_df$Cluster == temp[2]))]
    target_index <- intersect(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% target_clones),
                              intersect(which(JCC_Seurat_Obj$GMP == target_tp),
                                        which(JCC_Seurat_Obj$AllSeuratClusters == target_cluster)))
    seg.v_index <- which(seg.v[,"seg.name"] == sn)
    
    if(length(target_index) != length(seg.v_index)) {
      writeLines(paste("ERROR: target_index - ", sn))
    }
    
    seg.v[seg.v_index,"exp_TIGIT"] <- JCC_Seurat_Obj@assays$RNA@data["TIGIT",target_index]
    seg.v[seg.v_index,"exp_SELL"] <- JCC_Seurat_Obj@assays$RNA@data["SELL",target_index]
    seg.v[seg.v_index,"exp_CD27"] <- JCC_Seurat_Obj@assays$RNA@data["CD27",target_index]
    
    seg.f[seg.v_index,"clone"] <- JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[target_index]
    seg.f[seg.v_index,"time"] <- JCC_Seurat_Obj$time2[target_index]
    seg.f[seg.v_index,"barcode"] <- rownames(JCC_Seurat_Obj@meta.data)[target_index]
    
    opposite_tp <- setdiff(unique(JCC_Seurat_Obj$GMP), target_tp)
    
    for(i in seg.v_index) {
      target_index2 <- intersect(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta == seg.f$clone[i]),
                                 which(JCC_Seurat_Obj$GMP == opposite_tp))
      target_clusters2 <- JCC_Seurat_Obj$AllSeuratClusters[target_index2]
      target_clusters2 <- as.character(target_clusters2[order(target_clusters2)])
      
      seg.f$from_to[i] <- paste(target_clusters2, collapse = ";")
      
      total_link_num <- total_link_num + length(target_clusters2)
    }
  }
  seg.v <- data.frame(seg.v)
  
  ### reorder seg.f and seg.v
  unique_seg.names <- unique(seg.f$seg.name)
  for(seg in unique_seg.names) {
    ### get new ordered indicies for the given segment
    total_idx <- which(seg.f$seg.name == seg)
    single_idx <- intersect(total_idx,
                            which(!grepl(pattern = ";", x = seg.f$from_to, fixed = TRUE)))
    multiple_idx <- intersect(total_idx,
                              grep(pattern = ";", x = seg.f$from_to, fixed = TRUE))
    new_order_idx <- c(single_idx[order(as.numeric(seg.f$from_to[single_idx]), seg.f$clone[single_idx])], multiple_idx)
    
    ### reorder
    seg.f$from_to[total_idx] <- seg.f$from_to[new_order_idx]
    seg.f$clone[total_idx] <- seg.f$clone[new_order_idx]
    seg.f$time[total_idx] <- seg.f$time[new_order_idx]
    seg.f$barcode[total_idx] <- seg.f$barcode[new_order_idx]
    
    seg.v$exp_TIGIT[total_idx] <- seg.v$exp_TIGIT[new_order_idx]
    seg.v$exp_SELL[total_idx] <- seg.v$exp_SELL[new_order_idx]
    seg.v$exp_CD27[total_idx] <- seg.v$exp_CD27[new_order_idx]
  }
  
  ### which time points are subsist to PI cluster3 & 8?
  result_table <- data.frame(Barcode=seg.f$barcode,
                             Segment=seg.f$seg.name,
                             Time=seg.f$time,
                             GMP_PI=sapply(seg.f$seg.name, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1]),
                             Cluster=sapply(seg.f$seg.name, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2]),
                             From_To_Clusters=seg.f$from_to,
                             Clone=seg.f$clone,
                             stringsAsFactors = FALSE, check.names = FALSE)
  
  ### only select PI clusters
  result_table <- result_table[grep(pattern = "PI_", result_table$Segment, fixed = TRUE),]
  colnames(result_table)[which(colnames(result_table) == "From_To_Clusters")] <- "Lineage_From"
  
  ### expanded result_table
  expanded_result_table <- result_table
  multiple_idx <- grep(pattern = ";", expanded_result_table$Lineage_From, fixed = TRUE)
  for(idx in multiple_idx) {
    sub_clusters <- strsplit(expanded_result_table$Lineage_From[idx], split = ";", fixed = TRUE)[[1]]
    for(cluster in sub_clusters) {
      temp_row <- expanded_result_table[idx,]
      temp_row["Lineage_From"] <- cluster
      expanded_result_table <- rbind(expanded_result_table, temp_row)
    }
  }
  expanded_result_table <- expanded_result_table[-multiple_idx,]
  
  ### remove duplicates
  expanded_result_table$Size <- 1
  nodup_idx <- which(!duplicated(expanded_result_table))
  for(idx in nodup_idx) {
    expanded_result_table$Size[idx] <- length(intersect(intersect(intersect(which(expanded_result_table$Barcode == expanded_result_table$Barcode[idx]),
                                                                            which(expanded_result_table$Time == expanded_result_table$Time[idx])),
                                                                  intersect(which(expanded_result_table$Cluster == expanded_result_table$Cluster[idx]),
                                                                            which(expanded_result_table$Lineage_From == expanded_result_table$Lineage_From[idx]))),
                                                        which(expanded_result_table$Clone == expanded_result_table$Clone[idx])))
  }
  expanded_result_table <- expanded_result_table[nodup_idx,]
  expanded_result_table$Lineage_From <- paste0("GMP_", expanded_result_table$Lineage_From)
  
  ### PI_Cluster -> Time_Cluster
  expanded_result_table <- expanded_result_table[,-which(colnames(expanded_result_table) == "Barcode")]
  expanded_result_table$Segment2 <- paste0(expanded_result_table$Time, "_", expanded_result_table$Cluster)
  ### size adjustment - clone & lineage_from & segment2
  expanded_result_table$Segment3 <- paste0(expanded_result_table$Clone, "_", expanded_result_table$Lineage_From, "_", expanded_result_table$Segment2)
  nodup_idx <- which(!duplicated(expanded_result_table$Segment3))
  for(idx in nodup_idx) {
    expanded_result_table$Size[idx] <- sum(expanded_result_table$Size[which(expanded_result_table$Segment3 == expanded_result_table$Segment3[idx])])
  }
  expanded_result_table <- expanded_result_table[nodup_idx,]
  
  ### add PI_TO_PI column
  expanded_result_table$PI_TO_PI <- ""
  for(clone in unique(expanded_result_table$Clone)) {
    clone_idx <- which(expanded_result_table$Clone == clone)
    
    unique_tps <- unique(expanded_result_table$Time[clone_idx])
    
    if(length(unique_tps) > 1) {
      expanded_result_table$PI_TO_PI[clone_idx] <- "PI_TO_PI"
    }
  }
  
  ### add PI to PI connections
  interesting_table <- expanded_result_table[which(expanded_result_table$PI_TO_PI == "PI_TO_PI"),]
  ### remove unnecessary columns
  interesting_table <- interesting_table[,-which(colnames(interesting_table) %in% c("Segment3", "PI_TO_PI", ""))]
  ### if clone & segment2 are the same just combine them
  interesting_table$Temp_Col <- paste0(interesting_table$Clone, "_", interesting_table$Segment2)
  nodup_idx <- which(!duplicated(interesting_table$Temp_Col))
  for(idx in nodup_idx) {
    interesting_table$Size[idx] <- sum(interesting_table$Size[which(interesting_table$Temp_Col == interesting_table$Temp_Col[idx])])
  }
  interesting_table <- interesting_table[nodup_idx,]
  
  unique_seg_clones <- unique(interesting_table$Clone)
  for(clone in unique_seg_clones) {
    unique_tps <- unique(interesting_table$Time[which(interesting_table$Clone == clone)])
    unique_tps <- factor(unique_tps, levels = unique(JCC_Seurat_Obj$time2))
    unique_tps <- as.character(unique_tps[order(unique_tps)])
    
    for(i in 1:(length(unique_tps)-1)) {
      interesting_idx1 <- intersect(which(interesting_table$Clone == clone),
                                    which(interesting_table$Time == unique_tps[i]))
      interesting_idx2 <- intersect(which(interesting_table$Clone == clone),
                                    which(interesting_table$Time == unique_tps[i+1]))
      
      unique_segs1 <- unique(interesting_table$Segment2[interesting_idx1])
      unique_segs2 <- unique(interesting_table$Segment2[interesting_idx2])
      
      for(seg1 in unique_segs1) {
        for(seg2 in unique_segs2) {
          idx1 <- intersect(interesting_idx1,
                            which(interesting_table$Segment2 == seg1))
          idx2 <- intersect(interesting_idx2,
                            which(interesting_table$Segment2 == seg2))
          
          expanded_result_table <- rbind(expanded_result_table, c(interesting_table$Segment[idx2],
                                                                  interesting_table$Time[idx2],
                                                                  interesting_table$GMP_PI[idx2],
                                                                  interesting_table$Cluster[idx2],
                                                                  interesting_table$Segment2[idx1],
                                                                  interesting_table$Clone[idx2],
                                                                  as.integer(as.numeric(interesting_table$Size[idx1])*as.numeric(interesting_table$Size[idx2])),
                                                                  interesting_table$Segment2[idx2],
                                                                  paste0(interesting_table$Clone[idx2], "_", interesting_table$Segment2[idx1], "_", interesting_table$Segment2[idx2]),
                                                                  "PI_TO_PI_ADDED"))
        }
      }
    }
  }
  expanded_result_table$Size <- as.numeric(expanded_result_table$Size)
  
  ### there should not be duplicates in Segment3 column
  nodup_idx <- which(!duplicated(expanded_result_table$Segment3))
  for(idx in nodup_idx) {
    expanded_result_table$Size[idx] <- sum(expanded_result_table$Size[which(expanded_result_table$Segment3 == expanded_result_table$Segment3[idx])])
  }
  expanded_result_table <- expanded_result_table[nodup_idx,]
  
  ### save only the PI_TO_PI for future use
  pi_to_pi_result_table <- expanded_result_table[which(expanded_result_table$PI_TO_PI == "PI_TO_PI_ADDED"),]
  
  ### limit the connections that finally ended up with PI_3 or PI_8
  remove_idx <- NULL
  for(clone in unique(expanded_result_table$Clone)) {
    clone_idx <- which(expanded_result_table$Clone == clone)
    unique_seg <- unique(expanded_result_table$Segment[clone_idx])
    
    if(length(which(unique_seg %in% c("PI_3", "PI_8"))) == 0) {
      remove_idx <- c(remove_idx, clone_idx)
    }
  }
  if(length(remove_idx) > 0) {
    expanded_result_table <- expanded_result_table[-remove_idx,]
  }
  
  ### because there are too many nodes, just combine the GMPs
  expanded_result_table$Lineage_From2 <- sapply(expanded_result_table$Lineage_From, function(x) {
    if(grepl("GMP", x, fixed = TRUE)) {
      return("GMP")
    } else {
      return(x)
    }
  })
  expanded_result_table$Segment4 <- paste0(expanded_result_table$Lineage_From2, "_", expanded_result_table$Segment2)
  nodup_idx <- which(!duplicated(expanded_result_table$Segment4))
  for(idx in nodup_idx) {
    expanded_result_table$Size[idx] <- sum(expanded_result_table$Size[which(expanded_result_table$Segment4 == expanded_result_table$Segment4[idx])])
  }
  expanded_result_table <- expanded_result_table[nodup_idx,]
  
  ### remove the clone column because it's a wrong column now
  expanded_result_table <- expanded_result_table[,-which(colnames(expanded_result_table) %in% c("Clone", "Segment3"))]
  
  ### connection size recalculation - just all the involved cells based on the Segment4
  expanded_result_table$Size <- 0
  subsister_idx <- which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES")
  for(i in 1:nrow(expanded_result_table)) {
    from_time <- strsplit(expanded_result_table$Lineage_From2[i], "_", fixed = TRUE)[[1]][1]
    from_cluster <- strsplit(expanded_result_table$Lineage_From2[i], "_", fixed = TRUE)[[1]][2]
    to_time <- expanded_result_table$Time[i]
    to_cluster <- expanded_result_table$Cluster[i]
    
    if(is.na(from_cluster)) {
      clones1 <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(subsister_idx,
                                                                                        which(JCC_Seurat_Obj$time2 == from_time))])
      clones2 <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(subsister_idx,
                                                                                        intersect(which(JCC_Seurat_Obj$time2 == to_time),
                                                                                                  which(JCC_Seurat_Obj$AllSeuratClusters == to_cluster)))])
    } else {
      clones1 <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(subsister_idx,
                                                                                        intersect(which(JCC_Seurat_Obj$time2 == from_time),
                                                                                                  which(JCC_Seurat_Obj$AllSeuratClusters == from_cluster)))])
      clones2 <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(subsister_idx,
                                                                                        intersect(which(JCC_Seurat_Obj$time2 == to_time),
                                                                                                  which(JCC_Seurat_Obj$AllSeuratClusters == to_cluster)))])
    }
    intersected_clones <- intersect(clones1, clones2)
    
    if(length(intersected_clones) > 0) {
      if(is.na(from_cluster)) {
        expanded_result_table$Size[i] <- length(intersect(intersect(subsister_idx,
                                                                    which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% intersected_clones)),
                                                          which(JCC_Seurat_Obj$time2 == from_time))) + length(intersect(intersect(subsister_idx,
                                                                                                                                  which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% intersected_clones)),
                                                                                                                        intersect(which(JCC_Seurat_Obj$time2 == to_time),
                                                                                                                                  which(JCC_Seurat_Obj$AllSeuratClusters == to_cluster))))
      } else {
        expanded_result_table$Size[i] <- length(intersect(intersect(subsister_idx,
                                                                    which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% intersected_clones)),
                                                          intersect(which(JCC_Seurat_Obj$time2 == from_time),
                                                                    which(JCC_Seurat_Obj$AllSeuratClusters == from_cluster)))) + length(intersect(intersect(subsister_idx,
                                                                                                                                                            which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% intersected_clones)),
                                                                                                                                                  intersect(which(JCC_Seurat_Obj$time2 == to_time),
                                                                                                                                                            which(JCC_Seurat_Obj$AllSeuratClusters == to_cluster))))
      }
    }
  }
  
  ### pi_to_pi - remove duplicates and size adjustment - just count the unique single lineages
  pi_to_pi_result_table$Size <- 1
  pi_to_pi_result_table$Segment4 <- paste0(pi_to_pi_result_table$Lineage_From, "-", pi_to_pi_result_table$Segment2)
  nodup_idx <- which(!duplicated(pi_to_pi_result_table$Segment4))
  for(idx in nodup_idx) {
    pi_to_pi_result_table$Size[idx] <- sum(pi_to_pi_result_table$Size[which(pi_to_pi_result_table$Segment4 == pi_to_pi_result_table$Segment4[idx])])
  }
  pi_to_pi_result_table <- pi_to_pi_result_table[nodup_idx,]
  
  ### remove the 'Clone' and 'Segment3' columns because they are wrong now
  pi_to_pi_result_table <- pi_to_pi_result_table[,-which(colnames(pi_to_pi_result_table) %in% c("Clone", "Segment3"))]
  
  ### add 'from_cluster' column
  pi_to_pi_result_table$From_Cluster <- sapply(pi_to_pi_result_table$Lineage_From, function(x) strsplit(x, "_", TRUE)[[1]][2])
  
  ### make a data frame for alluvial plot
  plot_df3 <- data.frame(PI_Level=c(rep("Earliest Post-Infusion", nrow(pi_to_pi_result_table)), rep("Final Post-Infusion", nrow(pi_to_pi_result_table))),
                         Cluster=c(pi_to_pi_result_table$From_Cluster, pi_to_pi_result_table$Cluster),
                         Connection_Identifier=c(pi_to_pi_result_table$Segment4, pi_to_pi_result_table$Segment4),
                         Size=c(pi_to_pi_result_table$Size, pi_to_pi_result_table$Size),
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### Wes Anderson color palette
  sjcar19_color_scale <- colorRampPalette(c("#d73027", "#f46d43", "#fdae61", "#fee090", "#abd9e9", "#74add1", "#4575b4"))(length(unique(plot_df3$Cluster)))
  
  ### draw an alluvial plot
  plot_df3$Cluster <- factor(plot_df3$Cluster, levels = levels(JCC_Seurat_Obj$AllSeuratClusters))
  p <- ggplot(plot_df3,
              aes(x = PI_Level, stratum = Cluster, alluvium = Connection_Identifier,
                  y = Size,
                  fill = Cluster, label = Cluster)) +
    ggtitle("") +
    ylab("Unique Lineages") +
    geom_flow() +
    geom_stratum(alpha = 1) +
    geom_label_repel(stat = "stratum", size = 12, show.legend = FALSE, col = "black") +
    rotate_x_text(90) +
    scale_fill_manual(values = sjcar19_color_scale) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(size = 25, color = "black", face = "bold"),
          axis.title = element_text(size = 30, color = "black", face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'),
          legend.position = "right")
  ggsave(file = paste0(outputDir, "Fig4B_2.pdf"), plot = p,
         width = 15, height = 8, dpi = 350)
  
  
  ### Fig4C_1
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
