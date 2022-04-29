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
#               > generate_final(Seurat_RObj_path="Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Inputs/CARpos_JCC2.rds",
#                                outputDir="Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Outputs/")
###

generate_final <- function(Seurat_RObj_path="Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Inputs/CARpos_JCC2.rds",
                           outputDir="Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Outputs/") {
  
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
  if(!require(vegan, quietly = TRUE)) {
    install.packages("vegan")
    require(vegan, quietly = TRUE)
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
  interesting_genes <- c("RPL32", "RPL30", "LAG3", "TOX", "CASP8", "IL7R", "SELL", "MKI67",
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
  
  
  ### Fig2A
  ### NEW FUNCTIONAL ANNOTATION BY TAY - 09/27/2021
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters <- "Others"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "1", "5", "7", "9", "10", "11", "12", "15", "19", "21"))] <- "Proliferating"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("2", "4", "6", "17"))] <- "Transitioning"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("3", "8", "14"))] <- "Functional Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("13", "20"))] <- "Dysfunctional"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("16"))] <- "Early Effector"
  JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("18"))] <- "Metabolically Active"
  
  propsdf <- data.frame(time = character(), annotation = character(), numcells = numeric(), numcellstot = numeric(), stringsAsFactors = FALSE)
  for (c in unique(JCC_Seurat_Obj$New_Functional_Annotation_Based_On_Clusters)){
    for (t in unique(JCC_Seurat_Obj$time2)){
      
      temp <- JCC_Seurat_Obj@meta.data[which(JCC_Seurat_Obj@meta.data$time2 == t),]
      temp2 <- temp[which(temp$New_Functional_Annotation_Based_On_Clusters == c),]
      propsdf[nrow(propsdf) + 1,] <- c(t, c, nrow(temp2), nrow(temp))
    }
  }
  propsdf$numcells <- as.numeric(propsdf$numcells)
  propsdf$numcellstot <- as.numeric(propsdf$numcellstot)
  
  propsdf$proportion <- round(propsdf$numcells/propsdf$numcellstot, 3)
  
  ##removing Wk6 because we only have 1 car positive cell
  propsdf <- propsdf[which(propsdf$time != "Wk6" & propsdf$time != "6mo"),]
  
  propsdf$time <- factor(propsdf$time, levels = c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", "Wk8", "3mo"))
  
  ### add cell # in the col names
  adv_col_names <- levels(propsdf$time)
  cell_num <- c("(118,749)", "(21,085)", "(17,612)", "(19,872)", "(5,562)", "(1,535)", "(368)")
  adv_col_names <- paste(adv_col_names, cell_num, sep = "\n")
  
  sjcar19_colors <- c("#fee090", "#abd9e9", "#f46d43", "#74add1", "#d73027", "#4575b4")
  names(sjcar19_colors) <- unique(propsdf$annotation)
  
  p <- ggplot(data = propsdf,
         aes(x = time,
             y = proportion,
             alluvium = (annotation))) +
    geom_alluvium(aes(fill = (annotation)),
                  alpha = 1,
                  width = 1/2,
                  decreasing = NA, #this part stops it from changing the order in the bars, try T to see the difference
                  size = 0.7) +
    theme_classic(base_size = 36) +
    scale_fill_manual(values=sjcar19_colors)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = expand_scale(mult = c(0.01, 0.01))) + 
    ylab("Proportion of CAR+ Cells") +
    theme(axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=35, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0, vjust = 0.5, size = 30, color = "black", face = "bold"),
          plot.subtitle = element_text(hjust = 0, vjust = 0.5, size = 25, color = "black", face = "bold"),
          axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold")) +
    geom_vline(xintercept = "Wk2", linetype = "dashed", size=.8)
  
  ggsave(file = paste0(outputDir, "Fig2A.pdf"),
         plot = p, width = 22, height = 15, dpi = 350)
  
  
  ### Fig2B
  ### gray to sjcar19 blue color
  blue_color_scale <- colorRampPalette(c("gray85", "#08519c"))(length(unique(JCC_Seurat_Obj$time2)))
  names(blue_color_scale) <- unique(JCC_Seurat_Obj$time2)
  blue_color_scale["6mo"] <- "#08306b"
  
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
  target_tp <- target_tp[-which(target_tp %in% c("Wk6", "6mo"))]
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
  cell_num <- c("(118,749)", "(21,085)", "(17,612)", "(19,872)", "(5,562)", "(1,535)", "(368)")
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
  plot_df$Time <- factor(plot_df$Time, levels = rev(unique(JCC_Seurat_Obj$time2)))
  plot_df$Cluster <- factor(plot_df$Cluster, levels = unique(plot_df$Cluster))
  
  ### remove the Wk6 and 6mo
  plot_df <- plot_df[which(!plot_df$Time %in% c("Wk6", "6mo")),]
  
  ### color scale
  sjcar19_colors <- c("#16101E", "#D0B78F", "#8C8781", "#C7904F", "#133D31", "#82A5B8", "#3B3B53", "#4C8493", "#C31517",
                      "#D94C21", "#3E89A8", "#AA4C26",  "#CAA638", "#640B11", "#629488", "#BD7897", "#3C2D16", "#25245D",
                      "#E64E46", "#73BCAA", "#7047C1", "#286278")
  names(sjcar19_colors) <- unique(plot_df$Cluster)
  
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
  px_result_dir="Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Inputs/"
  
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
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40, color = "black", face = "bold"),
          axis.title = element_text(size = 25, color = "black", face = "bold"),
          axis.text = element_text(size = 20, color = "black", face = "bold")) +
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
  ### get some useful info
  JCC_Seurat_Obj$ALL_CARpos_Persister2 <- JCC_Seurat_Obj$ALL_CARpos_Persister
  JCC_Seurat_Obj$ALL_CARpos_Persister2[which(is.na(JCC_Seurat_Obj$ALL_CARpos_Persister2))] <- "NO"
  
  JCC_Seurat_Obj$Seurat_Clusters_Subsisters2 <- JCC_Seurat_Obj$ALL_CARpos_Persister2
  JCC_Seurat_Obj$Seurat_Clusters_Subsisters2[which(JCC_Seurat_Obj$Seurat_Clusters_Subsisters2 == "YES")] <- JCC_Seurat_Obj$time2[which(JCC_Seurat_Obj$Seurat_Clusters_Subsisters2 == "YES")]
  JCC_Seurat_Obj$Seurat_Clusters_Subsisters2[which(JCC_Seurat_Obj$Seurat_Clusters_Subsisters2 == "NO")] <- "Non-Subsisters"
  
  gmp_subsisters_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                                                  which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  pi_subsister_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                                which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  print(identical(gmp_subsisters_clones[order(gmp_subsisters_clones)], pi_subsister_clones[order(pi_subsister_clones)]))
  
  umap_map <- Embeddings(JCC_Seurat_Obj, reduction = "umap")[rownames(JCC_Seurat_Obj@meta.data), 1:2]
  
  ### only show GMP-> PI lineages
  ### should show GMP time points as well
  arrow_df <- data.frame(x1=0,
                         y1=0,
                         x2=0,
                         y2=0,
                         stringsAsFactors = FALSE, check.names = FALSE)
  cnt <- 1
  for(i in 1:length(gmp_subsisters_clones)) {
    target_indicies <- which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta == gmp_subsisters_clones[i])
    existing_time <- unique(JCC_Seurat_Obj$time2[target_indicies])
    if(length(existing_time) > 2 && existing_time[1] == "GMP") {
      target_indicies2 <- intersect(target_indicies,
                                    which(JCC_Seurat_Obj$time2 == existing_time[1]))
      arrow_df <- rbind(arrow_df, c(0, 0, 0, 0))
      arrow_df$x1[cnt] <- mean(umap_map[target_indicies2,
                                        "UMAP_1"])
      arrow_df$y1[cnt] <- mean(umap_map[target_indicies2,
                                        "UMAP_2"])
      target_indicies2 <- intersect(target_indicies,
                                    which(JCC_Seurat_Obj$time2 == existing_time[2]))
      arrow_df$x2[cnt] <- mean(umap_map[target_indicies2,
                                        "UMAP_1"])
      arrow_df$y2[cnt] <- mean(umap_map[target_indicies2,
                                        "UMAP_2"])
      cnt <- cnt + 1
    }
  }
  arrow_df <- arrow_df[-nrow(arrow_df),]
  
  ### color scale
  sjcar19_colors <- c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#91bfdb", "#4575b4", "lightgray")
  names(sjcar19_colors) <- c("GMP", "Wk1", "Wk2", "Wk3", "Wk4", "Wk8", "NA")
  
  ### add arrows to the previous UMAP
  JCC_Seurat_Obj$Seurat_Clusters_Subsisters2[which(JCC_Seurat_Obj$Seurat_Clusters_Subsisters2 == "Non-Subsisters")] <- "NA"
  p <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap", raster = TRUE,
               group.by = "Seurat_Clusters_Subsisters2",
               pt.size = 5,
               cols = sjcar19_colors,
               order = rev(unique(JCC_Seurat_Obj$Seurat_Clusters_Subsisters2))) +
    ggtitle("") +
    labs(color="") +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_classic(base_size = 48) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 40, color = "black", face = "bold"),
          axis.title = element_text(size = 25, color = "black", face = "bold"),
          axis.text = element_text(size = 20, color = "black", face = "bold")) +
    geom_curve(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = arrow_df,
      arrow = arrow(length = unit(0.03, "npc"))
    )
  ggsave(paste0(outputDir, "Fig4C_1.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### Fig4C_2
  ### set some useful info
  gmp_subsisters_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                                                  which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  pi_subsister_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                                which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  print(identical(gmp_subsisters_clones[order(gmp_subsisters_clones)], pi_subsister_clones[order(pi_subsister_clones)]))
  
  ### CD4/CD8 annotation
  JCC_Seurat_Obj$CD4_CD8_by_Clusters <- "NA"
  JCC_Seurat_Obj$CD4_CD8_by_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "2", "9", "10", "11", "14", "15", "18"))] <- "CD4"
  JCC_Seurat_Obj$CD4_CD8_by_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("1", "3", "5", "6", "7", "8", "12", "13", "16", "17", "19", "20"))] <- "CD8"
  
  plot_df3 <- data.frame(GMP_PI="",
                         Cluster="",
                         Connection_Identifier="",
                         Size=1,
                         stringsAsFactors = FALSE, check.names = FALSE)
  for(clone in gmp_subsisters_clones) {
    gmp_target_idx <- intersect(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta == clone),
                                intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                          which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8")))
    gmp_target_clusters <- unique(as.character(JCC_Seurat_Obj$AllSeuratClusters[gmp_target_idx]))
    pi_target_idx <- intersect(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta == clone),
                               intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                         which(JCC_Seurat_Obj$CD4_CD8_by_Clusters == "CD8")))
    pi_target_clusters <- unique(as.character(JCC_Seurat_Obj$AllSeuratClusters[pi_target_idx]))
    
    for(clstr1 in gmp_target_clusters) {
      for(clstr2 in pi_target_clusters) {
        plot_df3 <- rbind(plot_df3,
                          c("GMP Product", clstr1, paste0(clstr1, "_", clstr2), 1))
        plot_df3 <- rbind(plot_df3,
                          c("Post-Infusion", clstr2, paste0(clstr1, "_", clstr2), 1))
      }
    }
  }
  plot_df3 <- plot_df3[-1,]
  plot_df3$Size <- as.numeric(plot_df3$Size)
  
  ### sum up the duplicates
  nodup_idx <- which(!duplicated(plot_df3))
  for(i in nodup_idx) {
    target_idx <- intersect(intersect(which(plot_df3$GMP_PI == plot_df3$GMP_PI[i]),
                                      which(plot_df3$Cluster == plot_df3$Cluster[i])),
                            which(plot_df3$Connection_Identifier == plot_df3$Connection_Identifier[i]))
    plot_df3$Size[i] <- length(target_idx)
  }
  plot_df3 <- plot_df3[nodup_idx,]
  
  ### draw an alluvial plot
  plot_df3$Cluster <- factor(plot_df3$Cluster, levels = levels(JCC_Seurat_Obj$AllSeuratClusters))
  sjcar19_color_scale <- colorRampPalette(c("#d73027", "#f46d43", "#fdae61", "#fee090", "#abd9e9", "#74add1", "#4575b4"))(length(unique(plot_df3$Cluster)))
  p <- ggplot(plot_df3,
              aes(x = GMP_PI, stratum = Cluster, alluvium = Connection_Identifier,
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
  ggsave(file = paste0(outputDir, "Fig4C_2.pdf"), plot = p,
         width = 15, height = 8, dpi = 350)
  
  
  ### Fig4D
  ### setting some info
  gmp_subsisters_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                                                  which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  pi_subsister_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                                which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  print(identical(gmp_subsisters_clones[order(gmp_subsisters_clones)], pi_subsister_clones[order(pi_subsister_clones)]))
  
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
  
  ### get all the nodes
  node_names <- unique(c(expanded_result_table$Lineage_From2, expanded_result_table$Segment2))
  node_names <- as.character(node_names[order(node_names)])
  
  ### make adjacency matrix for network
  ### rows: outbound
  ### columns: inbound
  adj_mat <- matrix(0, length(node_names), length(node_names))
  rownames(adj_mat) <- node_names
  colnames(adj_mat) <- node_names
  
  ### fill out the adj matrix
  out_node_names <- unique(expanded_result_table$Lineage_From2)
  # out_node_names <- factor(out_node_names, levels = as.vector(sapply(unique(JCC_Seurat_Obj$time2), function(x) paste0(x, "_", levels(JCC_Seurat_Obj$AllSeuratClusters)))))
  out_node_names <- as.character(out_node_names[order(out_node_names)])
  in_node_names <- unique(expanded_result_table$Segment2)
  # in_node_names <- factor(in_node_names, levels = as.vector(sapply(unique(JCC_Seurat_Obj$time2), function(x) paste0(x, "_", levels(JCC_Seurat_Obj$AllSeuratClusters)))))
  in_node_names <- as.character(in_node_names[order(in_node_names)])
  for(r in out_node_names) {
    for(c in in_node_names) {
      target_idx <- intersect(which(expanded_result_table$Lineage_From2 == r),
                              which(expanded_result_table$Segment2 == c))
      if(length(target_idx) > 1) {
        writeLines(paste("ERROR: duplicated rows exist:", r, "\t", c, "-", length(target_idx)))
      } else if(length(target_idx) > 0) {
        adj_mat[r,c] <- expanded_result_table$Size[intersect(which(expanded_result_table$Lineage_From2 == r),
                                                             which(expanded_result_table$Segment2 == c))]
      }
    }
  }
  
  ### diagonal <- 0
  diag(adj_mat) <- 0
  
  ### remove nodes that do not have any connections (edges)
  remove_idx <- NULL
  for(i in 1:nrow(adj_mat)) {
    if((sum(adj_mat[i,]) == 0) && (sum(adj_mat[,i]) == 0)) {
      remove_idx <- c(remove_idx, i)
    }
  }
  if(length(remove_idx) > 0) {
    adj_mat <- adj_mat[-remove_idx, -remove_idx]
  }
  
  ### make an igraph
  g <- graph_from_adjacency_matrix(adj_mat, mode = "directed", weighted = TRUE)
  # coords <- layout_(g, as_tree())
  # plot(g, layout = coords)
  
  ### node color set
  total_clusters <- unique(c(sapply(expanded_result_table$Lineage_From2, function(x) {
    return(strsplit(x, split = "_", fixed = TRUE)[[1]][2])
  }), sapply(expanded_result_table$Segment2, function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]][2]
    return(strsplit(x, split = "_", fixed = TRUE)[[1]][2])
  })))
  total_clusters <- total_clusters[-which(is.na(total_clusters))]
  total_clusters <- total_clusters[order(as.numeric(total_clusters))]
  total_clusters <- c("GMP", total_clusters)
  sjcar19_color_scale <- c("#e0f3f8", "#d73027", "#f46d43", "#fdae61", "#fee090", "#abd9e9", "#74add1", "#4575b4")
  names(sjcar19_color_scale) <- total_clusters
  
  ### node and edge width + colors
  E(g)$width <- E(g)$weight + 5
  E(g)$edgeColor <- c(rep("#de77ae", length(which(adj_mat["GMP",] > 0))), rep("#542788", length(E(g))-length(which(adj_mat["GMP",] > 0))))
  V(g)$nodeSize <- sapply(V(g)$name, function(x) {
    if(x == "GMP") {
      r <- length(intersect(subsister_idx,
                            which(JCC_Seurat_Obj$time2 == x)))
    } else {
      temp <- strsplit(x, "_", fixed = TRUE)[[1]]
      r <- length(intersect(subsister_idx,
                            intersect(which(JCC_Seurat_Obj$time2 == temp[1]),
                                      which(JCC_Seurat_Obj$AllSeuratClusters == temp[2]))))
    }
    return(r)
  })
  V(g)$nodeSize <- (V(g)$nodeSize / max(V(g)$nodeSize, na.rm = TRUE)) * 100 + 20
  # V(g)$label.color <- sapply(V(g)$name, function(x) {
  #   tp <- strsplit(x, split = "_", fixed = TRUE)[[1]][1]
  #   return(wa_color_scale2[tp])
  # })
  V(g)$label.color <- "black"
  V(g)$color <- sapply(V(g)$name, function(x) {
    cls <- strsplit(x, split = "_", fixed = TRUE)[[1]][2]
    if(is.na(cls)) {
      cls <- "GMP"
    }
    return(sjcar19_color_scale[cls])
  })
  
  ### edge sizes should be smaller than node sizes
  if(max(E(g)$width) > (0.3 * max(V(g)$nodeSize))) {
    E(g)$width <- (E(g)$width / max(E(g)$width, na.rm = TRUE)) * (0.3 * max(V(g)$nodeSize)) 
  }
  
  ### increase the vertex label size
  V(g)$label.cex <- 3
  
  ### load RedeR screen and plot the graph
  rdp<-RedPort()
  calld(rdp)
  addGraph(rdp, g, layout.kamada.kawai(g))
  
  ### add legends
  # color
  addLegend.color(rdp, colvec=sjcar19_color_scale[-1], labvec=names(sjcar19_color_scale)[-1], title="Cluster",
                  vertical=FALSE, position="bottomleft", dxborder=10, dyborder=550, size=80, ftsize=50)
  
  # size
  circleLabel <- floor(seq(min(V(g)$nodeSize),max(V(g)$nodeSize),(max(V(g)$nodeSize) - min(V(g)$nodeSize))/4))
  circleSize <- (circleLabel / max(circleLabel)) * 100
  diag(adj_mat) <- NA
  circleLabel <- floor(seq(min(adj_mat, na.rm = TRUE), max(adj_mat, na.rm = TRUE),
                           ((max(adj_mat, na.rm = TRUE) - min(adj_mat, na.rm = TRUE))/4)))
  ### circle size in the legend should be at least 1
  if(circleSize[1] == 0) {
    circleLabel[1] <- 1
    circleSize[1] <- (1 / max(circleLabel)) * 100
  }
  if(circleLabel[1] == 0) {
    circleLabel[1] <- 1
  }
  addLegend.size(rdp,sizevec=circleSize,labvec=circleLabel,title="Cell Size",
                 position="bottomleft", dxborder=10, dyborder=10, ftsize=50)
  
  ### AFTER THIS, YOU HAVE TO MANUALLY LAYOUT THE NODES
  
  ### I ALSO HAVE A SAVED REDER FILE, SO YOU COULD LOAD IT
  ### run the following code and then in the popped up window, load the REDER file:
  ### Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Inputs/network_042822.reder
  rdp<-RedPort()
  calld(rdp)
  
  
  ### Fig4E
  ### setting some info
  gmp_subsisters_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                                                  which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  pi_subsister_clones <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                                which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  print(identical(gmp_subsisters_clones[order(gmp_subsisters_clones)], pi_subsister_clones[order(pi_subsister_clones)]))
  
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
  
  ### set link.pg.v
  link.pg.v <- data.frame(matrix(0, total_link_num/2, 6))
  colnames(link.pg.v) <- c("seg1", "start1", "end1", "seg2", "start2", "end2")
  
  cnt <- 1
  seg.f_gmp_idx <- which(startsWith(seg.f$seg.name, "GMP_"))
  for(i in seg.f_gmp_idx) {
    from_seg <- seg.f$seg.name[i]
    from_start <- seg.f$seg.start[i]
    from_end <- seg.f$seg.end[i]
    
    target_tp <- strsplit(from_seg, split = "_", fixed = TRUE)[[1]][1]
    opposite_tp <- setdiff(unique(JCC_Seurat_Obj$GMP), target_tp)
    
    to_clusters <- strsplit(seg.f$from_to[i], split = ";", fixed = TRUE)[[1]]
    unique_to_clusters <- unique(to_clusters)
    
    for(to_cluster in unique_to_clusters) {
      target_idx <- intersect(which(seg.f$seg.name == paste0(opposite_tp, "_", to_cluster)),
                              which(seg.f$clone == seg.f$clone[i]))
      for(idx in target_idx) {
        link.pg.v$seg1[cnt] <- from_seg
        link.pg.v$start1[cnt] <- from_start
        link.pg.v$end1[cnt] <- from_end
        link.pg.v$seg2[cnt] <- seg.f$seg.name[idx]
        link.pg.v$start2[cnt] <- seg.f$seg.start[idx]
        link.pg.v$end2[cnt] <- seg.f$seg.end[idx]
        cnt <- cnt + 1
      }
    }
  }
  
  ### there are too many links, so reduce them with threshold
  link_thresh <- 50
  unique_link_seg <- unique(link.pg.v$seg1)
  retain_idx <- NULL
  for(seg in unique_link_seg) {
    seg2_list <- sapply(unique(link.pg.v$seg2[which(link.pg.v$seg1 == seg)]), function(x) {
      return(length(intersect(which(link.pg.v$seg1 == seg),
                              which(link.pg.v$seg2 == x))))
    })
    
    target_segs <- names(seg2_list)[which(seg2_list > link_thresh)]
    if(length(target_segs) > 0) {
      writeLines(paste(seg, "-", target_segs))
    }
    
    retain_idx <- c(retain_idx, intersect(which(link.pg.v$seg1 == seg),
                                          which(link.pg.v$seg2 %in% target_segs)))
  }
  link.pg.v <- link.pg.v[retain_idx,]
  
  ## Add an alpha value to a colour
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  ### plot color
  connection <- paste0(link.pg.v$seg1, "_", link.pg.v$seg2)
  unique_connection <- unique(connection)
  line_colors3 <- add.alpha(c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#91bfdb", "#4575b4"), alpha = 0.7)
  names(line_colors3) <- unique_connection
  
  unique_seg.names <- unique(seg.f$seg.name)
  unique_seg_GMP_PI <- sapply(unique_seg.names, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1])
  node_colors <- colorRampPalette(colors=c("#de77ae","#542788"))(length(unique(unique_seg_GMP_PI)))
  names(node_colors) <- unique(unique_seg_GMP_PI)
  colorType <- node_colors[unique_seg_GMP_PI]
  cell_colors <- viridis(length(unique(seg.f$time)))
  names(cell_colors) <- unique(seg.f$time)
  colorType2 <- cell_colors[seg.f$time]
  
  ### set db
  db <- segAnglePo(seg.f, seg=seg.name, angle.start = 0, angle.end = 360)
  
  ### draw circular plot and save as pdf (without gene expression version)
  pdf(paste0(outputDir, "Fig4E.pdf"), width = 5, height = 5)
  par(mar = c(4, 2, 2, 4), xpd = TRUE)
  plot(c(1,1000), c(1,1000), type="n", axes=FALSE, xlab="", ylab="", main="")
  circos(xc=500, yc=440, R=440, cir=db, W=30, type="chr", col=colorType, print.chr.lab=TRUE, scale=FALSE, cex = 3)
  # circos(xc=500, yc=440, R=440, cir=db, W=30, mapping=seg.f, type="arc2", B=FALSE, col=colorType2, lwd=5, cutoff=0, cex = 3)
  # circos(xc=500, yc=440, R=400, cir=db, W=50, mapping=seg.v, col.v=3, type="l", B=TRUE, col="deepskyblue", lwd=2, scale=FALSE)
  # circos(xc=500, yc=440, R=340, cir=db, W=50, mapping=seg.v, col.v=4, type="l", B=TRUE, col="blue", lwd=2, scale=FALSE)
  # circos(xc=500, yc=440, R=280, cir=db, W=50, mapping=seg.v, col.v=5, type="l", B=TRUE, col="darkblue", lwd=2, scale=FALSE)
  circos(xc=500, yc=440, R=410, cir=db, W=50, mapping=link.pg.v, type="link.pg", lwd=2, col=line_colors3[connection])
  # legend("bottomright", 
  #        legend=c("TIGIT EXP", "SELL EXP", "CD27 EXP"),
  #        col=c("deepskyblue", "blue","darkblue"),
  #        pch=15, cex = 0.8, xpd = TRUE, inset = c(-0.12, -0.12))
  # legend("bottomright", 
  #        legend=names(cell_colors),
  #        col=cell_colors,
  #        pch=15, cex = 0.8, xpd = TRUE, inset = c(-0.2, -0.2))
  dev.off()
  
  
  ### Fig5A
  ### PI subsisters in the cluster3&8
  cluster3_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                   which(JCC_Seurat_Obj$AllSeuratClusters == "3")),
                                                                         which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
  cluster8_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                   which(JCC_Seurat_Obj$AllSeuratClusters == "8")),
                                                                         which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
  cluster38_pi_subsisters <- c(cluster3_pi_subsisters, cluster8_pi_subsisters)
  
  ### PI subsister clones in the cluster3&8
  cluster3_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster3_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
  cluster8_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster8_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
  cluster38_pi_subsister_clones <- unique(c(cluster3_pi_subsister_clones, cluster8_pi_subsister_clones))
  
  ### remove NA clones
  cluster3_pi_subsister_clones <- cluster3_pi_subsister_clones[which(!is.na(cluster3_pi_subsister_clones))]
  cluster8_pi_subsister_clones <- cluster8_pi_subsister_clones[which(!is.na(cluster8_pi_subsister_clones))]
  cluster38_pi_subsister_clones <- cluster38_pi_subsister_clones[which(!is.na(cluster38_pi_subsister_clones))]
  
  ### get the GMP subsister indicies that end up in PI cluster 3 & 8
  GMP_Subsisters_PI_Cluster38_idx <- intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                               which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% cluster38_pi_subsister_clones))
  
  ### add a column for the info
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 <- "Others"
  JCC_Seurat_Obj@meta.data[cluster38_pi_subsisters, "GMP_Subsisters_End_Up_In_Cluster38"] <- "PI_Subsisters_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38[GMP_Subsisters_PI_Cluster38_idx] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38[intersect(which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "Others"),
                                                              which(JCC_Seurat_Obj$GMP_CARpos_Persister == "YES"))] <- "Other_GMP_Subsisters"
  
  ### GMP subsisters end up in cluster3 and 8 vs all other GMPs
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 <- "Others"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2[setdiff(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                              which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8"))] <- "Other_GMPs"
  
  ### GMP subsisters end up in cluster3 and 8 vs all other CD8 GMPs
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 <- "Others"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[intersect(which(JCC_Seurat_Obj$AllSeuratClusters %in% c("1", "3", "5", "6", "7", "8", "12", "13", "16", "17", "19", "20")),
                                                                    which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 == "Other_GMPs"))] <- "Other_CD8_GMPs"
  
  ### filter out unwanted cells
  temp_seurat_obj <- subset(JCC_Seurat_Obj, cells = rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 %in% c("GMP_Subsisters_End_Up_In_Cluster_3_And_8", "Other_CD8_GMPs"))])
  
  ### check whether the orders are the same
  print(identical(rownames(temp_seurat_obj@meta.data), colnames(temp_seurat_obj@assays$RNA@counts)))
  print(identical(names(Idents(object = temp_seurat_obj)), rownames(temp_seurat_obj@meta.data)))
  
  ### factorize the column
  temp_seurat_obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 <- factor(temp_seurat_obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8,
                                                                     levels = c("GMP_Subsisters_End_Up_In_Cluster_3_And_8",
                                                                                "Other_CD8_GMPs"))
  
  ### dot plot - heatmap
  temp_seurat_obj <- SetIdent(object = temp_seurat_obj,
                              cells = rownames(temp_seurat_obj@meta.data),
                              value = temp_seurat_obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8)
  p <- DotPlot(temp_seurat_obj,
               features = c("EOMES", "TIGIT", "IFITM1", "IFITM2", "SELL", "CD27", "GNLY",
                            "GZMH", "KLRD1", "GZMK", "IFNG", "LAG3", "LEF1", "IL7R"),
               cols = c("#4575b4", "#d73027"),
               group.by = "GMP_Subsisters_End_Up_In_Cluster38_2_CD8") +
    scale_size(range = c(2, 18)) +
    coord_flip() +
    xlab("") +
    ylab("") +
    scale_y_discrete(labels = c("GMP Effector Precursors", "GMP Controls")) +
    theme_classic(base_size = 30) +
    theme(legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          legend.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0, vjust = 0.5, size = 30, color = "black", face = "bold"),
          plot.subtitle = element_text(hjust = 0, vjust = 0.5, size = 25, color = "black", face = "bold"),
          axis.title = element_blank(),
          axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 90, size = 25, vjust = 0.5, hjust = 1, color = "black", face = "bold"))
  ggsave(file = paste0(outputDir, "Fig5A.pdf"),
         plot = p, width = 8, height = 13, dpi = 350)
  
  ### second version - rotated
  p <- DotPlot(temp_seurat_obj,
               features = rev(c("EOMES", "TIGIT", "IFITM1", "IFITM2", "SELL", "CD27", "GNLY",
                                "GZMH", "KLRD1", "GZMK", "IFNG", "LAG3", "LEF1", "IL7R")),
               cols = c("#4575b4", "#d73027"),
               group.by = "GMP_Subsisters_End_Up_In_Cluster38_2_CD8") +
    scale_size(range = c(2, 18)) +
    xlab("") +
    ylab("") +
    scale_y_discrete(labels = c("GMP Effector Precursors", "GMP Controls")) +
    theme_classic(base_size = 30) +
    theme(legend.position = "top",
          legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          legend.text = element_text(angle = 0, size = 15, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0, vjust = 0.5, size = 30, color = "black", face = "bold"),
          plot.subtitle = element_text(hjust = 0, vjust = 0.5, size = 25, color = "black", face = "bold"),
          axis.title = element_blank(),
          axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text.x = element_text(angle = -45, size = 25, vjust = 0, hjust = 0.5, color = "black", face = "bold"))
  ggsave(file = paste0(outputDir, "Fig5A(2).pdf"),
         plot = p, width = 22, height = 5, dpi = 350)
  
  
  ### Fig5B
  ### IF YOU RUN THIS ON A PC, IT WILL TAKE VERY LONG TIME
  ### SO I RECOMMEND RUNNING IT ON HPC IN PARALLEL
  ### THE RESULT FIGURE Fig5B IS DERIVED FROM THE 541TH ITERATION
  
  ### new output directory
  outputDir2 <- paste0(outputDir, "/Classifier/")
  dir.create(outputDir2)
  
  ### change the class names
  JCC_Seurat_Obj$Classification_Class <- JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8
  JCC_Seurat_Obj$Classification_Class[which(JCC_Seurat_Obj$Classification_Class == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "YES"
  JCC_Seurat_Obj$Classification_Class[which(JCC_Seurat_Obj$Classification_Class == "Other_CD8_GMPs")] <- "NO"
  
  ### parameter setting for a classifier
  iteration <- 1000
  set.seed(1234)
  featureSelectionNum <- 100
  sampleNum <- 100
  methodTypes <- c("svmRadial")
  methodNames <- c("SVMRadial")
  train_control <- trainControl(method="none", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  print(identical(rownames(JCC_Seurat_Obj@meta.data), colnames(JCC_Seurat_Obj@assays$RNA@counts)))
  
  ### build the classifier 1000 times with random down-sampling
  acc <- vector("list", length = iteration)
  roc <- vector("list", length = iteration)
  for(i in 1:iteration) {
    
    writeLines(paste(i, "/", iteration))
    
    ### get downsampling indicies
    downsample_idx <- c(sample(which(JCC_Seurat_Obj$Classification_Class == "YES"), sampleNum),
                        sample(which(JCC_Seurat_Obj$Classification_Class == "NO"), sampleNum))
    
    ### get normalized data
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(JCC_Seurat_Obj@assays$RNA@counts[,downsample_idx],
                                                                stringsAsFactors = FALSE, check.names = FALSE))
    
    ### manual LOOCV
    acc[[i]] <- rep(NA, length(downsample_idx))
    roc[[i]] <- rep(NA, length(downsample_idx))
    for(j in 1:length(downsample_idx)) {
      
      ### divide training and test index
      test_idx <- downsample_idx[j]
      training_idx <- downsample_idx[-j]
      
      ### new obj for the training data & set idents with the info
      classifier_seurat_obj <- subset(JCC_Seurat_Obj, cells = rownames(JCC_Seurat_Obj@meta.data)[training_idx])
      classifier_seurat_obj <- SetIdent(object = classifier_seurat_obj,
                                        cells = rownames(classifier_seurat_obj@meta.data),
                                        value = classifier_seurat_obj$Classification_Class)
      
      ### DE analysis
      de_result <- FindMarkers(classifier_seurat_obj,
                               ident.1 = "YES",
                               ident.2 = "NO",
                               min.pct = 0.1,
                               logfc.threshold = 0.1,
                               test.use = "wilcox",
                               verbose = FALSE)
      de_result <- de_result[order(de_result$p_val_adj),]
      de_genes <- rownames(de_result)[1:featureSelectionNum]
      
      ### annotate class for the input data
      input_data2 <- data.frame(t(input_data[intersect(de_genes, rownames(input_data)),rownames(classifier_seurat_obj@meta.data)]), stringsAsFactors = FALSE, check.names = FALSE)
      input_data2$Class <- factor(classifier_seurat_obj$Classification_Class,
                                  levels = c("YES", "NO"))
      
      ### annotate class for the test data
      test_data <- data.frame(t(input_data[intersect(de_genes, rownames(input_data)),rownames(JCC_Seurat_Obj@meta.data)[test_idx],drop=FALSE]), stringsAsFactors = FALSE, check.names = FALSE)
      test_data$Class <- factor(JCC_Seurat_Obj$Classification_Class[downsample_idx[j]],
                                levels = c("YES", "NO"))
      
      ### get accuracy and pred prob (not AUC yet)
      model <- train(Class~., data=input_data2, method=methodTypes, trControl = train_control)
      pred_result <- predict(model, newdata = test_data)
      acc[[i]][j] <- round((sum(pred_result == test_data$Class) / nrow(test_data)), 3)
      pred_result <- predict(model, newdata = test_data, type = "prob")
      roc[[i]][j] <- pred_result$YES
      
    }
    
    ### calculate the AUCs
    roc[[i]] <- roc(JCC_Seurat_Obj$Classification_Class[downsample_idx], roc[[i]])
    
    ### draw ROC curve
    pdf(paste0(outputDir2, "Classifier_GMP_Precursor_vs_Other_CD8_GMPs_(", i, ").pdf"),
        width = 5, height = 4)
    plot.roc(roc[[i]], main = paste("Accuracy =", mean(as.numeric(acc[[i]]))),
             legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
             xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    dev.off()
    
  }
  
  
  ### FigS1A
  ### GMP/PI UMAP
  p <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
               group.by = "GMP",
               cols = c("GMP" = "#d73027", "PI" = "#4575b4"),
               pt.size = 2, raster = TRUE, label = FALSE) +
    ggtitle("") +
    labs(color="GMP/PI") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.9
  ggsave(paste0(outputDir, "FigS1A.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### FigS1B
  ### CD4/CD8 annotation
  JCC_Seurat_Obj$CD4_CD8_by_Clusters <- "NA"
  JCC_Seurat_Obj$CD4_CD8_by_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("0", "2", "9", "10", "11", "14", "15", "18"))] <- "CD4"
  JCC_Seurat_Obj$CD4_CD8_by_Clusters[which(JCC_Seurat_Obj$AllSeuratClusters %in% c("1", "3", "5", "6", "7", "8", "12", "13", "16", "17", "19", "20"))] <- "CD8"
  
  ### UMAP
  p <- DimPlot(object = JCC_Seurat_Obj, reduction = "umap",
               group.by = "CD4_CD8_by_Clusters",
               cols = c("CD8" = "#de77ae", "CD4" = "#542788", "NA" = "gray"),
               order = c("CD8", "CD4", "NA"),
               pt.size = 2, raster = TRUE) +
    ggtitle("") +
    labs(color="CD4/CD8") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.9
  ggsave(paste0(outputDir, "FigS1B.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### FigS1C
  ### violin plot of cluster3 vs cluster 8
  JCC_Seurat_Obj_3_8 <- subset(JCC_Seurat_Obj,
                               cells = rownames(JCC_Seurat_Obj@meta.data)[union(which(JCC_Seurat_Obj$AllSeuratClusters == "3"),
                                                                                which(JCC_Seurat_Obj$AllSeuratClusters == "8"))])
  
  ### color scale
  sjcar19_colors <- c("#4575b4", "#d73027")
  names(sjcar19_colors) <- c("3", "8")
  
  ### violin plot
  JCC_Seurat_Obj_3_8 <- SetIdent(object = JCC_Seurat_Obj_3_8,
                                 cells = rownames(JCC_Seurat_Obj_3_8@meta.data),
                                 value = JCC_Seurat_Obj_3_8@meta.data$AllSeuratClusters)
  p <- VlnPlot(JCC_Seurat_Obj_3_8, features = c("GZMK"),
               pt.size = 0, cols = sjcar19_colors) +
    ggtitle("GZMK Expression") +
    xlab("Cluster") +
    # stat_summary(fun=mean, geom="point", size=3, color="black") +
    theme_classic(base_size = 40) +
    theme(legend.key.size = unit(3, 'cm'),
          legend.position = "none",
          legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          legend.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35, color = "black", face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 25, color = "black", face = "bold"),
          axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
  
  ### save the violin plot
  ggsave(file = paste0(outputDir, "FigS1C.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### FigS1D
  ### UMAP with cell density
  sjcar19_colors <- c("white", "#313695", "#4575b4", "#74add1", "#abd9e9", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026")
  p <- ggplot(JCC_Seurat_Obj@meta.data, aes(x = JCC_Seurat_Obj@reductions$umap@cell.embeddings[,1],
                                            y = JCC_Seurat_Obj@reductions$umap@cell.embeddings[,2])) +
    geom_density_2d_filled(contour_var = "ndensity") +
    scale_fill_manual(values = sjcar19_colors) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(fill = "Level") +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key = element_rect()) +
    scale_x_continuous(expand = c(0, 0), limits = c(-11, 11)) + scale_y_continuous(expand = c(0, 0), limits = c(-8, 9))
  ggsave(paste0(outputDir, "FigS1D.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### FigS1E
  ### dot plot3
  interesting_genes3 <- c("MKI67", "CDK1", "TYMS", "TUBB", "LTB", "IL7R", "TCF7", "IL10RA", "GNLY",
                          "GZMK", "PRF1", "NKG7", "GZMB", "MCM5", "MCM7", "CDCA7", "GZMA",
                          "PLK1", "CDC20", "HLA-DPA1", "CD74", "IL2RA", "HILPDA", "BNIP3",
                          "TOX", "LAG3", "CD27", "TIGIT", "ENO1", "PGK1", "TOP2A", "CASP8",
                          "HIST1H1E", "HIST1H4F", "HIST1H3C", "HIST1H2AL", "HIST1H2AI",
                          "HIST1H3H", "HIST1H1C", "HIST1H1D", "HIST1H2AH")
  p <- DotPlot(JCC_Seurat_Obj,
               features = interesting_genes3,
               group.by = "AllSeuratClusters") +
    scale_size(range = c(5, 15)) +
    xlab("") +
    ylab("Clusters") +
    scale_color_gradientn(colours = c("#313695", "#ffffbf", "#a50026")) +
    theme_classic(base_size = 35) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.title = element_text(size = 35, hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 90, size = 30, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'))
  ggsave(file = paste0(outputDir, "FigS1E.pdf"),
         plot = p, width = 36, height = 14, dpi = 350)
  
  
  ### FigS2A
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
  monocle_cds2 <- orderCells(monocle_cds2, root_state = "1")
  
  ### by pseudotime
  p <- plot_cell_trajectory(monocle_cds2, color_by = "Pseudotime", cell_size = 3, cell_link_size = 3, show_branch_points = FALSE) +
    labs(color="Pseudotime") +
    theme_classic(base_size = 36) +
    scale_color_gradientn(colours = c("#313695", "#ffffbf", "#a50026"),
                          n.breaks = 3) +
    theme(legend.position = "top",
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 30, color = "black", face = "bold"),
          legend.title = element_text(size = 36, color = "black", face = "bold"),
          legend.text = element_text(size = 30, color = "black", face = "bold"))
  ggsave(file = paste0(outputDir, "FigS2A.pdf"),
         plot = p,
         width = 15, height = 10, dpi = 350)
  
  
  ### FigS2B
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
  
  ### by each cluster
  sjcar19_colors <- c("#16101E", "#D0B78F", "#8C8781", "#C7904F", "#133D31", "#82A5B8", "#3B3B53", "#4C8493", "#C31517",
                      "#D94C21", "#3E89A8", "#AA4C26",  "#CAA638", "#640B11", "#629488", "#BD7897", "#3C2D16", "#25245D",
                      "#E64E46", "#73BCAA", "#7047C1", "#286278")
  names(sjcar19_colors) <- levels(JCC_Seurat_Obj$AllSeuratClusters)
  p <- plot_cell_trajectory(monocle_cds2, color_by = "AllSeuratClusters", cell_size = 3, cell_link_size = 2, show_branch_points = FALSE) +
    labs(color="Clusters") +
    scale_color_manual(values = sjcar19_colors) +
    theme_classic(base_size = 36) +
    theme(legend.position = "none",
          text = element_text(size = 36, color = "black", face = "bold"),
          axis.title = element_text(size = 36, color = "black", face = "bold"),
          axis.text = element_text(size = 30, color = "black", face = "bold"),
          legend.title = element_text(size = 36, color = "black", face = "bold"),
          legend.text = element_text(size = 30, color = "black", face = "bold")) +
    facet_wrap(~AllSeuratClusters, ncol = 7)
  ggsave(file = paste0(outputDir, "FigS2B.pdf"),
         plot = p,
         width = 20, height = 15, dpi = 350)
  
  
  ### FigS3B
  ### target lineages
  subsister_clones_gmp <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                                                 which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  rest_carpos_clones_gmp <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                                                                   which(JCC_Seurat_Obj$ALL_CARpos_Persister == "NO"))])
  rest_carpos_clones_gmp <- rest_carpos_clones_gmp[which(!is.na(rest_carpos_clones_gmp))]
  
  ### gmp indicies
  gmp_carpos_indicies <- intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                   which(JCC_Seurat_Obj$CAR == "CARpos"))
  
  ### get persister vs rest carpos clone sizes in GMP
  subsister_clone_sizes_gmp <- rep(1, length(subsister_clones_gmp))
  names(subsister_clone_sizes_gmp) <- subsister_clones_gmp
  for(clone in subsister_clones_gmp) {
    subsister_clone_sizes_gmp[clone] <- length(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[gmp_carpos_indicies] == clone))
  }
  
  rest_carpos_clone_sizes_gmp <- rep(0, length(rest_carpos_clones_gmp))
  names(rest_carpos_clone_sizes_gmp) <- rest_carpos_clones_gmp
  for(clone in rest_carpos_clones_gmp) {
    rest_carpos_clone_sizes_gmp[clone] <- length(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[gmp_carpos_indicies] == clone))
  }
  
  
  ### target lineages
  subsister_clones_pi <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                                which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))])
  rest_carpos_clones_pi <- unique(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                                  which(JCC_Seurat_Obj$ALL_CARpos_Persister == "NO"))])
  rest_carpos_clones_pi <- rest_carpos_clones_pi[which(!is.na(rest_carpos_clones_pi))]
  
  ### pi indicies
  pi_carpos_indicies <- intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                  which(JCC_Seurat_Obj$CAR == "CARpos"))
  
  ### get persister vs rest carpos clone sizes in PI time points
  subsister_clone_sizes_pi <- rep(1, length(subsister_clones_pi))
  names(subsister_clone_sizes_pi) <- subsister_clones_pi
  for(clone in subsister_clones_pi) {
    subsister_clone_sizes_pi[clone] <- length(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[pi_carpos_indicies] == clone))
  }
  
  rest_carpos_clone_sizes_pi <- rep(0, length(rest_carpos_clones_pi))
  names(rest_carpos_clone_sizes_pi) <- rest_carpos_clones_pi
  for(clone in rest_carpos_clones_pi) {
    rest_carpos_clone_sizes_pi[clone] <- length(which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta[pi_carpos_indicies] == clone))
  }
  
  ### draw a cumulative plot
  max_clone_size <- Reduce(max, c(subsister_clone_sizes_gmp, rest_carpos_clone_sizes_gmp,
                                  subsister_clone_sizes_pi, rest_carpos_clone_sizes_pi))
  pdf(file = paste0(outputDir, "FigS3B.pdf"), width = 5, height = 4)
  plot(ecdf(subsister_clone_sizes_gmp), col="#d7191c", lwd = 1,
       main = "Cumulative Clone Sizes",
       sub = paste0("KS test p-value: ",
                    "GMP=", formatC(ks.test(subsister_clone_sizes_gmp, rest_carpos_clone_sizes_gmp)$p.value, format = "e", digits = 2), ", ",
                    "PI=", formatC(ks.test(subsister_clone_sizes_pi, rest_carpos_clone_sizes_pi)$p.value, format = "e", digits = 2)),
       xlab = "Clone Size",
       ylab = "Proportion",
       xaxs="i", xlim=c(0,max_clone_size), ylim=c(0.75,1))
  lines(ecdf(rest_carpos_clone_sizes_gmp), col="#fdae61", lwd = 1)
  lines(ecdf(subsister_clone_sizes_pi), col="#abd9e9", lwd = 1)
  lines(ecdf(rest_carpos_clone_sizes_pi), col="#2c7bb6", lwd = 1)
  legend("bottomright", 
         legend=c("GMP Lineages", "GMP Controls",
                  "Post-Infusion Lineages","Post-Infusion Controls"),
         col=c("#d7191c", "#fdae61","#abd9e9","#2c7bb6"),
         pch=15)
  dev.off()
  
  
  ### FigS4A
  ### PI subsisters in the cluster3&8
  cluster3_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                   which(JCC_Seurat_Obj$AllSeuratClusters == "3")),
                                                                         which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
  cluster8_pi_subsisters <- rownames(JCC_Seurat_Obj@meta.data)[intersect(intersect(which(JCC_Seurat_Obj$GMP == "PI"),
                                                                                   which(JCC_Seurat_Obj$AllSeuratClusters == "8")),
                                                                         which(JCC_Seurat_Obj$ALL_CARpos_Persister == "YES"))]
  cluster38_pi_subsisters <- c(cluster3_pi_subsisters, cluster8_pi_subsisters)
  
  ### PI subsister clones in the cluster3&8
  cluster3_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster3_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
  cluster8_pi_subsister_clones <- unique(JCC_Seurat_Obj@meta.data[cluster8_pi_subsisters,"clonotype_id_by_patient_one_alpha_beta"])
  cluster38_pi_subsister_clones <- unique(c(cluster3_pi_subsister_clones, cluster8_pi_subsister_clones))
  
  ### remove NA clones
  cluster3_pi_subsister_clones <- cluster3_pi_subsister_clones[which(!is.na(cluster3_pi_subsister_clones))]
  cluster8_pi_subsister_clones <- cluster8_pi_subsister_clones[which(!is.na(cluster8_pi_subsister_clones))]
  cluster38_pi_subsister_clones <- cluster38_pi_subsister_clones[which(!is.na(cluster38_pi_subsister_clones))]
  
  ### get the GMP subsister indicies that end up in PI cluster 3 & 8
  GMP_Subsisters_PI_Cluster38_idx <- intersect(which(JCC_Seurat_Obj$GMP == "GMP"),
                                               which(JCC_Seurat_Obj$clonotype_id_by_patient_one_alpha_beta %in% cluster38_pi_subsister_clones))
  
  ### add a column for the info
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 <- "Others"
  JCC_Seurat_Obj@meta.data[cluster38_pi_subsisters, "GMP_Subsisters_End_Up_In_Cluster38"] <- "PI_Subsisters_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38[GMP_Subsisters_PI_Cluster38_idx] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38[intersect(which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "Others"),
                                                              which(JCC_Seurat_Obj$GMP_CARpos_Persister == "YES"))] <- "Other_GMP_Subsisters"
  
  ### GMP subsisters end up in cluster3 and 8 vs other CD8 GMP subsisters
  JCC_Seurat_Obj@meta.data$GMP_Subsisters_End_Up_In_Cluster38_CD8 <- "Others"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_CD8[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_CD8[intersect(which(JCC_Seurat_Obj$AllSeuratClusters %in% c("1", "3", "5", "6", "7", "8", "12", "13", "16", "17", "19", "20")),
                                                                  which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "Other_GMP_Subsisters"))] <- "Other_CD8_GMP_Subsisters"
  
  ### GMP subsisters end up in cluster3 and 8 vs all other GMPs
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 <- "Others"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2[setdiff(which(JCC_Seurat_Obj$GMP == "GMP"),
                                                              which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8"))] <- "Other_GMPs"
  
  ### GMP subsisters end up in cluster3 and 8 vs all other CD8 GMPs
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 <- "Others"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP_Subsisters_End_Up_In_Cluster_3_And_8"
  JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[intersect(which(JCC_Seurat_Obj$AllSeuratClusters %in% c("1", "3", "5", "6", "7", "8", "12", "13", "16", "17", "19", "20")),
                                                                    which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2 == "Other_GMPs"))] <- "Other_CD8_GMPs"
  
  ### genes of interest
  genes_of_interest <- c("TIGIT", "KLRD1", "CD86", "IL2RA", "CD70",
                         "LAG3", "CD7", "SELL", "CD27", "IL7R")
  genes_of_interest <- intersect(genes_of_interest,
                                 rownames(JCC_Seurat_Obj@assays$RNA@counts))
  
  ### Ridge plot
  temp_obj <- subset(JCC_Seurat_Obj,
                     cells = rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 %in% c("GMP_Subsisters_End_Up_In_Cluster_3_And_8",
                                                                                                                                     "Other_CD8_GMPs"))])
  temp_obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[which(temp_obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 == "GMP_Subsisters_End_Up_In_Cluster_3_And_8")] <- "GMP Precursors"
  temp_obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8[which(temp_obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8 == "Other_CD8_GMPs")] <- "GMP Controls"
  temp_obj <- SetIdent(object = temp_obj,
                       cells = rownames(temp_obj@meta.data),
                       value = temp_obj$GMP_Subsisters_End_Up_In_Cluster38_2_CD8)
  p <- RidgePlot(temp_obj, features = genes_of_interest,
                 ncol = 2,
                 cols = c("#d73027", "#4575b4"))
  for(i in 1:length(genes_of_interest)) {
    p[[i]] <- p[[i]] +
      labs(y = "") +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35, color = "black", face = "bold"),
            axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
            axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
  }
  ggsave(paste0(outputDir, "FigS4A.pdf"), plot = p, width = 15, height = 18, dpi = 350)
  
  
  ### FigS6
  ### load turtle data first
  GSE125881_path <- "Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Inputs/GSE125881_RAW/"
  f <- list.files(path = GSE125881_path)
  f_name <- sapply(f, function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    temp2 <- strsplit(temp[3], split = "-", fixed = TRUE)[[1]]
    return(paste0(temp[2], "_", temp2[1]))
  })
  
  turtle_tcrs <- vector("list", length = length(f))
  names(turtle_tcrs) <- f_name
  
  for(i in 1:length(f)) {
    temp <- read.csv(file = paste0(GSE125881_path, f[i]),
                     stringsAsFactors = FALSE, check.names = FALSE)
    temp$CDR3_AA <- paste0(temp$chain, ":", temp$cdr3s)
    
    tcr_frequency <- temp$frequency
    names(tcr_frequency) <- temp$CDR3_AA
    
    turtle_tcrs[[i]] <- tcr_frequency
  }
  
  ### load our tcr data
  TCR_contig_path <- "Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Inputs/ContigAnnots/"
  f2 <- list.files(path = TCR_contig_path)
  f_name2 <- sapply(f2, function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[3], "_", temp[4]))
  })
  f_name2[1] <- "SJCAR19-02_Wk3"
  f_name2[2] <- "SJCAR19-02_GMP"
  f_name2[3] <- "SJCAR19-00_GMP"
  f_name2[4] <- "SJCAR19-01_GMP"
  f_name2[5] <- "SJCAR19-03_Wk3"
  f_name2[6] <- "SJCAR19-03_GMP"
  f_name2[7] <- "SJCAR19-04_GMP"
  f_name2[10] <- "SJCAR19-05_GMP"
  f_name2[11] <- "SJCAR19-06_GMP"
  f_name2[14] <- "SJCAR19-07_GMP"
  
  our_tcrs <- vector("list", length = length(f2))
  names(our_tcrs) <- f_name2
  
  for(i in 1:length(f2)) {
    our_tcrs[[i]] <- read.csv(file = paste0(TCR_contig_path, f2[i]),
                              stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  ### filter and calculate the frequency
  for(i in 1:length(our_tcrs)) {
    temp <- our_tcrs[[i]]
    temp <- temp[which(temp$cdr3 != "None"),]
    
    temp$CDR3_AA <- paste0(temp$chain, ":", temp$cdr3)
    unique_tcrs <- unique(temp$CDR3_AA)
    
    tcr_frequency <- rep(1, length(unique_tcrs))
    names(tcr_frequency) <- unique_tcrs
    
    dup_tcrs <- unique(temp$CDR3_AA[which(duplicated(temp$CDR3_AA))])
    
    for(dup in dup_tcrs) {
      tcr_frequency[dup] <- length(which(temp$CDR3_AA == dup))
    }
    
    our_tcrs[[i]] <- tcr_frequency
  }
  
  ### divide the TCRs into alpha and beta
  our_tcrs_alpha <- our_tcrs
  our_tcrs_beta <- our_tcrs
  for(i in 1:length(our_tcrs)) {
    alpha_idx <- grep("TRA:", names(our_tcrs_alpha[[i]]))
    beta_idx <- grep("TRB:", names(our_tcrs_beta[[i]]))
    
    our_tcrs_alpha[[i]] <- our_tcrs_alpha[[i]][alpha_idx]
    our_tcrs_beta[[i]] <- our_tcrs_beta[[i]][beta_idx]
  }
  
  turtle_tcrs_alpha <- turtle_tcrs
  turtle_tcrs_beta <- turtle_tcrs
  for(i in 1:length(turtle_tcrs)) {
    alpha_idx <- grep("TRA:", names(turtle_tcrs_alpha[[i]]))
    beta_idx <- grep("TRB:", names(turtle_tcrs_beta[[i]]))
    
    turtle_tcrs_alpha[[i]] <- turtle_tcrs_alpha[[i]][alpha_idx]
    turtle_tcrs_beta[[i]] <- turtle_tcrs_beta[[i]][beta_idx]
  }
  
  ### we have GMP & GMP-redo, so combine them
  info_table <- data.frame(library=names(our_tcrs),
                           px=sapply(names(our_tcrs), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1]),
                           time=sapply(names(our_tcrs), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][2]),
                           stringsAsFactors = FALSE, check.names = FALSE)
  info_table$time2 <- info_table$time
  info_table$time2[which(info_table$time2 == "GMP-redo")] <- "GMP"
  info_table$new_lib <- paste0(info_table$px, "_", info_table$time2)
  
  dup_idx <- which(duplicated(info_table$new_lib))
  for(idx in dup_idx) {
    original_idx <- which(info_table$new_lib == info_table$new_lib[idx])[1]
    
    ### alpha
    temp1 <- our_tcrs_alpha[[original_idx]]
    temp2 <- our_tcrs_alpha[[idx]]
    
    common_clones <- intersect(names(temp1), names(temp2))
    
    temp3 <- c(temp1, temp2[setdiff(names(temp2), names(temp1))])
    for(clone in common_clones) {
      temp3[clone] <- temp3[clone] + temp2[clone]
    }
    
    our_tcrs_alpha[[original_idx]] <- temp3
    
    ### beta
    temp1 <- our_tcrs_beta[[original_idx]]
    temp2 <- our_tcrs_beta[[idx]]
    
    common_clones <- intersect(names(temp1), names(temp2))
    
    temp3 <- c(temp1, temp2[setdiff(names(temp2), names(temp1))])
    for(clone in common_clones) {
      temp3[clone] <- temp3[clone] + temp2[clone]
    }
    
    our_tcrs_beta[[original_idx]] <- temp3
  }
  our_tcrs_alpha <- our_tcrs_alpha[-dup_idx]
  our_tcrs_beta <- our_tcrs_beta[-dup_idx]
  
  ### rename
  names(our_tcrs_alpha) <- info_table$new_lib[-dup_idx]
  names(our_tcrs_beta) <- info_table$new_lib[-dup_idx]
  
  ### the number of unique TCRs - normalized by total TCR #
  normalized_unique_tcr_num_turtle <- rep(0, length(turtle_tcrs))
  names(normalized_unique_tcr_num_turtle) <- names(turtle_tcrs)
  normalized_unique_tcr_num_ours <- rep(0, length(our_tcrs))
  names(normalized_unique_tcr_num_ours) <- names(our_tcrs)
  for(i in 1:length(turtle_tcrs)) {
    normalized_unique_tcr_num_turtle[i] <- length(turtle_tcrs[[i]]) / sum(turtle_tcrs[[i]])
  }
  for(i in 1:length(our_tcrs)) {
    normalized_unique_tcr_num_ours[i] <- length(our_tcrs[[i]]) / sum(our_tcrs[[i]])
  }
  
  ### the number of unique TCRs - normalized by total TCR # (beta only)
  normalized_unique_tcr_num_turtle_beta <- rep(0, length(turtle_tcrs_beta))
  names(normalized_unique_tcr_num_turtle_beta) <- names(turtle_tcrs_beta)
  normalized_unique_tcr_num_ours_beta <- rep(0, length(our_tcrs_beta))
  names(normalized_unique_tcr_num_ours_beta) <- names(our_tcrs_beta)
  len_turtle_beta <- rep(0, length(turtle_tcrs_beta))
  names(len_turtle_beta) <- names(turtle_tcrs_beta)
  len_ours_beta <- rep(0, length(our_tcrs_beta))
  names(len_ours_beta) <- names(our_tcrs_beta)
  for(i in 1:length(turtle_tcrs_beta)) {
    normalized_unique_tcr_num_turtle_beta[i] <- length(turtle_tcrs_beta[[i]]) / sum(turtle_tcrs_beta[[i]])
    len_turtle_beta[i] <- sum(turtle_tcrs_beta[[i]])
  }
  for(i in 1:length(our_tcrs_beta)) {
    normalized_unique_tcr_num_ours_beta[i] <- length(our_tcrs_beta[[i]]) / sum(our_tcrs_beta[[i]])
    len_ours_beta[i] <- sum(our_tcrs_beta[[i]])
  }
  
  ### various diversity indicies
  shannon_turtle_beta <- rep(0, length(turtle_tcrs_beta))
  names(shannon_turtle_beta) <- names(turtle_tcrs_beta)
  shannon_our_beta <- rep(0, length(our_tcrs_beta))
  names(shannon_our_beta) <- names(our_tcrs_beta)
  
  simpson_turtle_beta <- rep(0, length(turtle_tcrs_beta))
  names(simpson_turtle_beta) <- names(turtle_tcrs_beta)
  simpson_our_beta <- rep(0, length(our_tcrs_beta))
  names(simpson_our_beta) <- names(our_tcrs_beta)
  
  fisher_turtle_beta <- rep(0, length(turtle_tcrs_beta))
  names(fisher_turtle_beta) <- names(turtle_tcrs_beta)
  fisher_our_beta <- rep(0, length(our_tcrs_beta))
  names(fisher_our_beta) <- names(our_tcrs_beta)
  
  richness_turtle_beta <- rep(0, length(turtle_tcrs_beta))
  names(richness_turtle_beta) <- names(turtle_tcrs_beta)
  richness_our_beta <- rep(0, length(our_tcrs_beta))
  names(richness_our_beta) <- names(our_tcrs_beta)
  
  ### fill in the indicies
  for(i in 1:length(turtle_tcrs_beta)) {
    shannon_turtle_beta[i] <- diversity(turtle_tcrs_beta[[i]], index = "shannon")
    simpson_turtle_beta[i] <- diversity(turtle_tcrs_beta[[i]], index = "simpson")
    fisher_turtle_beta[i] <- fisher.alpha(turtle_tcrs_beta[[i]])
    richness_turtle_beta[i] <- specnumber(turtle_tcrs_beta[[i]])
  }
  for(i in 1:length(our_tcrs_beta)) {
    shannon_our_beta[i] <- diversity(our_tcrs_beta[[i]], index = "shannon")
    simpson_our_beta[i] <- diversity(our_tcrs_beta[[i]], index = "simpson")
    fisher_our_beta[i] <- fisher.alpha(our_tcrs_beta[[i]])
    richness_our_beta[i] <- specnumber(our_tcrs_beta[[i]])
  }
  
  ### make a plot data frame
  plot_df <- data.frame(Name=c(names(normalized_unique_tcr_num_turtle_beta),
                               names(normalized_unique_tcr_num_ours_beta)),
                        Source=c(rep("Shieh et al.", length(normalized_unique_tcr_num_turtle_beta)),
                                 rep("SJCAR19 Cohort", length(normalized_unique_tcr_num_ours_beta))),
                        Time=c(c("IP", "Wk3", "Wk5", "4mo", "IP", "Wk2", "Wk4", "3mo", "IP", "Wk2", "Wk4", "3mo", "IP", "Wk2", "Wk4", "3mo"),
                               sapply(names(normalized_unique_tcr_num_ours_beta), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][2])),
                        Patient=c(sapply(names(normalized_unique_tcr_num_turtle_beta), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1]),
                                  sapply(names(normalized_unique_tcr_num_ours_beta), function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])),
                        NormRich=c(normalized_unique_tcr_num_turtle_beta,
                                   normalized_unique_tcr_num_ours_beta),
                        ShannonIdx=c(shannon_turtle_beta,
                                     shannon_our_beta),
                        SimpsonIdx=c(simpson_turtle_beta,
                                     simpson_our_beta),
                        FisherIdx=c(fisher_turtle_beta,
                                    fisher_our_beta),
                        Richness=c(richness_turtle_beta,
                                   richness_our_beta),
                        TotalNum=c(len_turtle_beta,
                                   len_ours_beta),
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### only IP/GMP
  plot_df2 <- plot_df[which(plot_df$Time %in% c("IP", "GMP", "GMP-redo")),]
  
  ### beeswarm plot
  plot_df2$Source <- factor(plot_df2$Source, levels = c("SJCAR19 Cohort", "Shieh et al."))
  p <- ggplot(plot_df2, aes_string(x="Source", y="ShannonIdx")) +
    geom_boxplot() +
    ggtitle("Diversity - Shannon Index") +
    geom_beeswarm(aes_string(col="Source"), na.rm = TRUE, show.legend = FALSE, size = 5, cex = 3) +
    stat_compare_means(size = 5, hjust = -0.5) +
    xlab("") + ylab("") +
    labs(col="") +
    scale_color_manual(values = c("#4575b4", "#d73027")) +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30, color = "black", face = "bold"),
          plot.subtitle = element_text(hjust = 0, vjust = 0.5, size = 25, color = "black", face = "bold"),
          legend.text = element_text(hjust = 0, vjust = 0.5, size = 30, color = "black", face = "bold"),
          axis.title = element_blank(),
          axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir, "FigS6.pdf"),
         plot = p, width = 12, height = 8, dpi = 350)
  
  
  ### TableS1
  
  
  
  
  
  ### R4C3
  ### load the DE result
  allmarkers_cd4_de_result <- read.table(file = paste0("Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Inputs/JCC212_CARpos_JCC_CD4_allMarkers.tsv"),
                                         header = TRUE,
                                         stringsAsFactors = FALSE, check.names = FALSE)
  allmarkers_cd8_de_result <- read.table(file = paste0("Z:/ResearchHome/Groups/thomagrp/home/common/Hyunjin/JCC212_SJCAR19/Final/Inputs/JCC212_CARpos_JCC_CD8_allMarkers.tsv"),
                                         header = TRUE,
                                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### aggregate the result - top5 from each cluster
  top_5_genes <- vector("list", length(unique(JCC_Seurat_Obj$AllSeuratClusters)))
  names(top_5_genes) <- levels(JCC_Seurat_Obj$AllSeuratClusters)
  for(clstr in names(top_5_genes)) {
    if(clstr %in% c("4", "21")) {
      top_5_genes[[clstr]] <- intersect(allmarkers_cd8_de_result$gene[which(allmarkers_cd8_de_result$cluster == clstr)][1:15],
                                        allmarkers_cd4_de_result$gene[which(allmarkers_cd4_de_result$cluster == clstr)][1:15])[1:5]
    } else if(clstr %in% unique(allmarkers_cd4_de_result$cluster)) {
      top_5_genes[[clstr]] <- allmarkers_cd4_de_result$gene[which(allmarkers_cd4_de_result$cluster == clstr)][1:5]
    } else if(clstr %in% unique(allmarkers_cd8_de_result$cluster)) {
      top_5_genes[[clstr]] <- allmarkers_cd8_de_result$gene[which(allmarkers_cd8_de_result$cluster == clstr)][1:5]
    } else {
      writeLines(paste0("Error"))
    }
  }
  
  ### prepare a heatmap data
  heatdata <- data.frame(matrix(0, nrow = length(unlist(top_5_genes)), ncol = length(levels(JCC_Seurat_Obj$AllSeuratClusters))),
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(heatdata) <- paste0(as.vector(sapply(levels(JCC_Seurat_Obj$AllSeuratClusters), function(x) rep(x, 5))),
                               "_",
                               unlist(top_5_genes))
  colnames(heatdata) <- levels(JCC_Seurat_Obj$AllSeuratClusters)
  for(row in rownames(heatdata)) {
    target_gene <- strsplit(row, split = "_", fixed = TRUE)[[1]][2]
    for(col in colnames(heatdata)) {
      heatdata[row,col] <- mean(JCC_Seurat_Obj@assays$RNA@data[target_gene,
                                                               rownames(JCC_Seurat_Obj@meta.data)[which(JCC_Seurat_Obj$AllSeuratClusters == col)]])
    }
  }
  
  ### a function for color brewer
  cell_pal <- function(cell_vars, pal_fun) {
    if (is.numeric(cell_vars)) {
      pal <- pal_fun(100)
      return(pal[cut(cell_vars, breaks = 100)])
    } else {
      categories <- sort(unique(cell_vars))
      pal <- setNames(pal_fun(length(categories)), categories)
      return(pal[cell_vars])
    }
  }
  
  ### get colors for the clustering result
  cell_colors_clust <- c("#16101E", "#D0B78F", "#8C8781", "#C7904F", "#133D31", "#82A5B8", "#3B3B53", "#4C8493", "#C31517",
                      "#D94C21", "#3E89A8", "#AA4C26",  "#CAA638", "#640B11", "#629488", "#BD7897", "#3C2D16", "#25245D",
                      "#E64E46", "#73BCAA", "#7047C1", "#286278")
  names(cell_colors_clust) <- levels(JCC_Seurat_Obj$AllSeuratClusters)
  heatclus <- as.vector(sapply(levels(JCC_Seurat_Obj$AllSeuratClusters), function(x) rep(x, 5)))
  
  ### make a heatmap with those genes
  png(paste0(outputDir, "R4C3_Heatmap_Top5_From_Each_Cluster.png"),
      width = 3000, height = 3000, res = 350)
  par(oma=c(0,2,0,3))
  heatmap.2(as.matrix(heatdata), col = colorpanel(36, low = "#2c7bb6", high = "#d7191c"),
            scale = "row", dendrogram = "none", trace = "none",
            cexRow = 0.4, cexCol = 2, key.title = "", main = "Top 5 Genes From Each Cluster",
            Rowv = FALSE, labRow = unlist(top_5_genes),
            Colv = FALSE, labCol = levels(JCC_Seurat_Obj$AllSeuratClusters),
            key.xlab = "Norm.Count", key.ylab = "Frequency",
            RowSideColors = cell_colors_clust[heatclus])
  legend("left", inset = -0.1, y.intersp = 0.7,
         box.lty = 0, cex = 1.1,
         title = "Marker Cluster", xpd = TRUE,
         legend=names(cell_colors_clust),
         col=cell_colors_clust,
         pch=15)
  dev.off()
  
  ### change the pink purple to
  #762a83
  #1b7837
  
  
  
  
  
  
  
  
  
}
