###
#   File name : PCA_with_all_time.R
#   Author    : Hyunjin Kim
#   Date      : Apr 1, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Put all the cells in one and draw a PCA to find if there is
#               any difference among different time points
#
#   Instruction
#               1. Source("PCA_with_all_time.R")
#               2. Run the function "draw_umap" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PCA_with_all_time.R/PCA_with_all_time.R")
#               > draw_umap(Seurat_RObj_path="./data/JCC212_Px5.Robj",
#                           outputDir="./results/")
###

draw_pca <- function(Seurat_RObj_path="./data/JCC212_Px5.Robj",
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
  
  ### load the object
  load(Seurat_RObj_path)
  
  # ### see columns of the meta data
  # colnames(JCC212_Px5@meta.data)
  # 
  # ### draw an UMAP plot split by different time points
  # DimPlot(JCC212_Px5, reduction = "umap",
  #         split.by = "Time",
  #         cols = rep("black", nrow(JCC212_Px5@meta.data))) +
  #   theme(legend.position = "none")
  # 
  # ### draw one UMAP plot with all the time points
  # DimPlot(JCC212_Px5, reduction = "umap",
  #         cols = rep("black", nrow(JCC212_Px5@meta.data))) +
  #   theme(legend.position = "none")
  # 
  # ### draw one UMAP plot with all the time points and color each time point
  # DimPlot(JCC212_Px5, reduction = "umap",
  #         group.by = "Time")
  # 
  # ### draw one PCA plot with all the time points and color each time point
  # DimPlot(JCC212_Px5, reduction = "pca",
  #         group.by = "Time")
  
  ### put PCA results in the slot
  JCC212_Px5$Full_PC1 <- JCC212_Px5@reductions$pca@cell.embeddings[,"PC_1"]
  JCC212_Px5$Full_PC2 <- JCC212_Px5@reductions$pca@cell.embeddings[,"PC_2"]
  
  # ### draw one PCA plot with all the time points and color each time point - transparent dots
  # ggplot(JCC212_Px5@meta.data, aes_string(x="Full_PC1", y="Full_PC2", group="TimeF")) +
  #   geom_point(aes_string(color="TimeF"), size=1.5, alpha=0.5) +
  #   xlab("PC1") + ylab("PC2") +
  #   theme_classic(base_size = 16)
  
  
  ### a function that returns multi-figures from PCA, TSNE, AND UMAP
  ### One plot with all the groups + each group
  ### x: a numeric vector of first component
  ### y: a numeric vector of second component
  ### group: a factor vector of group info
  ### type: type of the dimensionality reduction method
  multiReducPlot <- function(x, y, group, type=c("PCA", "TSNE", "UMAP"), isPrint=FALSE) {
    
    ### load library
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
    
    ### create a data frame for ggplot
    group <- factor(as.character(group),
                    levels = intersect(levels(group), unique(group)))
    plot_df <- data.frame(X=x, Y=y, Group=group)
    
    ### set x & y axes labels
    if(type[1] == "PCA") {
      x_label <- "PC1"
      y_label <- "PC2"
    } else if(type[1] == "TSNE") {
      x_label <- "TSNE1"
      y_label <- "TSNE2"
    } else if(type[1] == "UMAP") {
      x_label <- "UMAP1"
      y_label <- "UMAP2"
    } else {
      stop("ERROR: type parameter should be either \"PCA\", \"TSNE\", or \"UMAP\".")
    }
    
    ### set colors for each group
    col_palette <- hue_pal()(length(levels(group)))
    names(col_palette) <- levels(group)
    
    ### 1. One plot with all the groups + each group
    p <- vector("list", length(levels(group))+1)
    x_range <- c(min(plot_df$X), max(plot_df$X))
    y_range <- c(min(plot_df$Y), max(plot_df$Y))
    p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
      geom_point(aes_string(col="Group"), size=1.5, alpha=0.5) +
      xlab(x_label) + ylab(y_label) +
      xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
      ggtitle("All") +
      theme_classic(base_size = 16) +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    for(i in 2:(length(levels(group))+1)) {
      p[[i]] <- ggplot(plot_df[which(plot_df$Group == levels(group)[i-1]),], aes_string(x="X", y="Y")) +
        geom_point(col=col_palette[levels(group)[i-1]], size=1.5, alpha=0.5) +
        xlab(x_label) + ylab(y_label) +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle(levels(group)[i-1]) +
        theme_classic(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, 
                                        color = col_palette[levels(group)[i-1]]))
    }
    ### arrange the plots
    fName <- paste0("Whole Group & Each")
    g <- arrangeGrob(grobs = p,
                     nrow = ceiling(sqrt(length(levels(group))+1)),
                     ncol = ceiling(sqrt(length(levels(group))+1)),
                     top = fName)
    
    if(isPrint) {
      plot(g)
    }
    
    return(g)
    
  }
  
  ### Put all the cells in one and draw a UMAP to find if there is any difference among different time points
  g <- multiReducPlot(x = JCC212_Px5$Full_PC1, y = JCC212_Px5$Full_PC2,
                      group = JCC212_Px5$TimeF, type = "PCA")
  ggsave(file = paste0(outputDir, "Px5_PCA_dissected.png"), g, width = 22, height = 12, dpi = 300)
  
}
