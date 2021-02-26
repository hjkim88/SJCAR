###
#   File name : GMP_Persistency_Marker_Discovery.R
#   Author    : Hyunjin Kim
#   Date      : Dec 29, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : What markers predict whether a lineage (CAR+/CAR-) will persist after infusion?
#               - DE analysis
#               - GO/Pathway analysis
#               - A heatmap to show those genes are actually differentially expressed between the two conditions
#               - Comparison of various factors between persisters & non-persisters (violin plot?)
#               - Using the DE genes to make a classifier
#               - To ensure an individual does not drive the found pattern?
#                 a) if the classifier (LOOCV) worked well, it is a proof
#                 b) multiple 1d scatter plots (colored based on patients) to show that the top DE genes are not one-patient specific
#               - Classifier comparison between 1. using the top DE genes vs 2. using random genes
#               - if there are too many persisters from one patient, then build a classifier with the patient only
#                 and test with other cells
#               - Reversed PCA/UMAP (gene-patient conversion) to show whether those genes are separated by persistency
#
#               * After meeting with Jeremy at Jan 11,
#               === TO-DO List ===
#               1. Make statistics
#               - The number of lineages
#               - The number of persister cells
#               - Statistics info with lineages specifically from GMP
#               - CD4/CD8 - associated lineage info
#               2. Classifier (Remove CD4 cells)
#               - Use only one random cell per each lineage and permute to see if there are many changes
#               - Use all the cells of all the lineages but LOPOCV (Leave One Patient Out -> Patient-based)
#               3. Attach V & J info, and also the responder info to the Seurat object
#
#   Instruction
#               1. Source("GMP_Persistency_Marker_Discovery.R")
#               2. Run the function "persistency_study" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_GMP_Persistency_Marker_Discovery.R/GMP_Persistency_Marker_Discovery.R")
#               > persistency_study(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                   clonotype_lineage_info_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Clonotype_Lineages.RDS",
#                                   lineage_in_full_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Lineages_by_CAR/SJCAR19_Lineages_in_Full.RDS",
#                                   outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Persistency/")
###

persistency_study <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj.RDS",
                              clonotype_lineage_info_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Clonotype_Lineages.RDS",
                              lineage_in_full_path="./results/New2/Lineages_by_CAR/SJCAR19_Lineages_in_Full.RDS",
                              outputDir="./results/New2/Persistency/") {
  
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
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(grid, quietly = TRUE)) {
    install.packages("grid")
    require(grid, quietly = TRUE)
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
  if(!require(ggsci, quietly = TRUE)) {
    install.packages("ggsci")
    require(ggsci, quietly = TRUE)
  }
  if(!require(viridis, quietly = TRUE)) {
    install.packages("viridis")
    require(viridis, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(caret, quietly = TRUE)) {
    install.packages("caret")
    require(caret, quietly = TRUE)
  }
  if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
  }
  if(!require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    require(RColorBrewer, quietly = TRUE)
  }
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    library(gplots, quietly = TRUE)
  }
  if(!require(uwot, quietly = TRUE)) {
    install.packages("uwot")
    require(uwot, quietly = TRUE)
  }
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(ggthemes, quietly = TRUE)) {
    install.packages("ggthemes")
    require(ggthemes, quietly = TRUE)
  }
  if(!require(msigdbr, quietly = TRUE)) {
    install.packages("msigdbr")
    library(msigdbr, quietly = TRUE)
  }
  
  # ******************************************************************************************
  # Pathway Analysis with clusterProfiler package
  # Input: geneList     = a vector of gene Entrez IDs for pathway analysis [numeric or character]
  #        org          = organism that will be used in the analysis ["human" or "mouse"]
  #                       should be either "human" or "mouse"
  #        database     = pathway analysis database (KEGG or GO) ["KEGG" or "GO"]
  #        title        = title of the pathway figure [character]
  #        pv_threshold = pathway analysis p-value threshold (not DE analysis threshold) [numeric]
  #        displayNum   = the number of pathways that will be displayed [numeric]
  #                       (If there are many significant pathways show the few top pathways)
  #        imgPrint     = print a plot of pathway analysis [TRUE/FALSE]
  #        dir          = file directory path of the output pathway figure [character]
  #
  # Output: Pathway analysis results in figure - using KEGG and GO pathways
  #         The x-axis represents the number of DE genes in the pathway
  #         The y-axis represents pathway names
  #         The color of a bar indicates adjusted p-value from the pathway analysis
  #         For Pathview Result, all colored genes are found DE genes in the pathway,
  #         and the color indicates log2(fold change) of the DE gene from DE analysis
  # ******************************************************************************************
  pathwayAnalysis_CP <- function(geneList,
                                 org,
                                 database,
                                 title="Pathway_Results",
                                 pv_threshold=0.05,
                                 displayNum=Inf,
                                 imgPrint=TRUE,
                                 dir="./") {
    
    ### load library
    if(!require(clusterProfiler, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("clusterProfiler")
      require(clusterProfiler, quietly = TRUE)
    }
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    
    ### collect gene list (Entrez IDs)
    geneList <- geneList[which(!is.na(geneList))]
    
    if(!is.null(geneList)) {
      ### make an empty list
      p <- list()
      
      if(database == "KEGG") {
        ### KEGG Pathway
        kegg_enrich <- enrichKEGG(gene = geneList, organism = org, pvalueCutoff = pv_threshold)
        
        if(is.null(kegg_enrich)) {
          writeLines("KEGG Result does not exist")
          return(NULL)
        } else {
          kegg_enrich@result <- kegg_enrich@result[which(kegg_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(kegg_enrich@result) <= displayNum)) {
              result <- kegg_enrich@result
              description <- kegg_enrich@result$Description
            } else {
              result <- kegg_enrich@result[1:displayNum,]
              description <- kegg_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(kegg_enrich) > 0) {
              p[[1]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title)) +
                theme(axis.text = element_text(size = 25))
              
              png(paste0(dir, "kegg_", title, ".png"), width = 2000, height = 1000)
              print(p[[1]])
              dev.off()
            } else {
              writeLines("KEGG Result does not exist")
            }
          }
          
          return(kegg_enrich@result)
        }
      } else if(database == "GO") {
        ### GO Pathway
        if(org == "human") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else if(org == "mouse") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else {
          go_enrich <- NULL
          writeLines(paste("Unknown org variable:", org))
        }
        
        if(is.null(go_enrich)) {
          writeLines("GO Result does not exist")
          return(NULL)
        } else {
          go_enrich@result <- go_enrich@result[which(go_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(go_enrich@result) <= displayNum)) {
              result <- go_enrich@result
              description <- go_enrich@result$Description
            } else {
              result <- go_enrich@result[1:displayNum,]
              description <- go_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(go_enrich) > 0) {
              p[[2]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("GO ", title)) +
                theme(axis.text = element_text(size = 25))
              
              png(paste0(dir, "go_", title, ".png"), width = 2000, height = 1000)
              print(p[[2]])
              dev.off()
            } else {
              writeLines("GO Result does not exist")
            }
          }
          
          return(go_enrich@result)
        }
      } else {
        stop("database prameter should be \"GO\" or \"KEGG\"")
      }
    } else {
      writeLines("geneList = NULL")
    }
  }
  
  #'******************************************************************************
  #' A function to transform RNA-Seq data with VST in DESeq2 package
  #' readCount: RNA-Seq rawcounts in a matrix or in a data frame form
  #'            Rows are genes and columns are samples
  #' filter_thresh: The function filters out genes that have at least one sample
  #'                with counts larger than the 'filter_thresh' value
  #'                e.g., if the 'filter_thresh' = 1, then it removes genes
  #'                that have counts <= 1 across all the samples
  #'                if 0, then there will be no filtering
  #'******************************************************************************
  normalizeRNASEQwithVST <- function(readCount, filter_thresh=1) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filter_thresh > 0) {
      ### Remove rubbish rows - this will decrease the number of rows
      keep = apply(counts(deSeqData), 1, function(r){
        return(sum(r > filter_thresh) > 0)
      })
      deSeqData <- deSeqData[keep,]
    }
    
    ### VST
    vsd <- varianceStabilizingTransformation(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  ###
  #   This function was downloaded from: https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
  ###
  heatmap.3 <- function(x,
                        Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                        distfun = dist,
                        hclustfun = hclust,
                        dendrogram = c("both","row", "column", "none"),
                        symm = FALSE,
                        scale = c("none","row", "column"),
                        na.rm = TRUE,
                        revC = identical(Colv,"Rowv"),
                        add.expr,
                        breaks,
                        symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                        col = "heat.colors",
                        colsep,
                        rowsep,
                        sepcolor = "white",
                        sepwidth = c(0.05, 0.05),
                        cellnote,
                        notecex = 1,
                        notecol = "cyan",
                        na.color = par("bg"),
                        trace = c("none", "column","row", "both"),
                        tracecol = "cyan",
                        hline = median(breaks),
                        vline = median(breaks),
                        linecol = tracecol,
                        margins = c(5,5),
                        ColSideColors,
                        RowSideColors,
                        side.height.fraction=0.3,
                        cexRow = 0.2 + 1/log10(nr),
                        cexCol = 0.2 + 1/log10(nc),
                        labRow = NULL,
                        labCol = NULL,
                        key = TRUE,
                        keysize = 1.5,
                        density.info = c("none", "histogram", "density"),
                        denscol = tracecol,
                        symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                        densadj = 0.25,
                        main = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        lmat = NULL,
                        lhei = NULL,
                        lwid = NULL,
                        ColSideColorsSize = 1,
                        RowSideColorsSize = 1,
                        KeyValueName="Value",...){
    
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
      if (is.list(x))
        return(all(sapply(x, invalid)))
      else if (is.vector(x))
        return(all(is.na(x)))
      else return(FALSE)
    }
    
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
      "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale=\"row\" or scale=\"column\" when breaks are",
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
      Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
      Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
        if (is.logical(Colv) && (Colv))
          dendrogram <- "column"
        else dedrogram <- "none"
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
        if (is.logical(Rowv) && (Rowv))
          dendrogram <- "row"
        else dendrogram <- "none"
        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
      rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
      colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x)))
        (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x)))
        (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      
      if (!missing(ColSideColors)) {
        #if (!is.matrix(ColSideColors))
        #stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
          stop("'ColSideColors' must be a matrix of nrow(x) rows")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
        lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
      }
      
      if (!missing(RowSideColors)) {
        #if (!is.matrix(RowSideColors))
        #stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
          stop("'RowSideColors' must be a matrix of ncol(x) columns")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }
    
    if (length(lhei) != nrow(lmat))
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors)){
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      } else {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = t(RowSideColors[,rowInd, drop=F])
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
          rsc.colors[rsc.i] = rsc.name
          rsc[rsc == rsc.name] = rsc.i
          rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(rownames(RowSideColors)) > 0) {
          axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE, cex.axis = cexRow/2)
        }
      }
    }
    
    if (!missing(ColSideColors)) {
      
      if (!is.matrix(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      } else {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop=F]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
          csc.colors[csc.i] = csc.name
          csc[csc == csc.name] = csc.i
          csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
          axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE, cex.axis = cexCol/2)
        }
      }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals, col = linecol,
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
        if (!is.null(hline)) {
          abline(h = i + hline, col = linecol, lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote))
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
           col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      par(mar = c(5, 4, 2, 1), cex = 0.75)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
      }
      
      z <- seq(min.raw, max.raw, length = length(col))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(1, at = xv, labels = lv)
      if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
      else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
      else mtext(side = 1, KeyValueName, line = 2)
      if (density.info == "density") {
        dens <- density(x, adjust = densadj, na.rm = TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
              lwd = 1)
        axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
        title("Color Key\nand Density Plot")
        par(cex = 0.5)
        mtext(side = 2, "Density", line = 2)
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
              col = denscol)
        axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
        title("Color Key\nand Histogram")
        par(cex = 0.5)
        mtext(side = 2, "Count", line = 2)
      }
      else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
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
  
  pClones <- NULL
  indeterminate_patient_pool <- NULL
  for(i in 1:length(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]])) {
    gmp_idx <- which(colnames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) == "GMP")
    gmp_redo_idx <- which(colnames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) == "GMP-redo")
    last_gmp_idx <- max(gmp_idx, gmp_redo_idx)
    
    if(length(gmp_redo_idx) > 0) {
      ### if at least GMP or GMP-redo exist and there are at least one afterward-time point
      if((last_gmp_idx != -Inf) && (ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) - last_gmp_idx > 1)) {
        ### collect persistent clones that appeared in GMP and persist afterwards
        if(nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) > 0) {
          for(j in 1:nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])) {
            for(k in last_gmp_idx:(ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])-1)) {
              if((SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,"GMP"] > 0 ||
                  SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,"GMP-redo"] > 0) &&
                 SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,k] > 0) {
                pClones <- c(pClones, rownames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])[j])
                break;
              }
            }
          }
        }
      }
    } else {
      ### if at least GMP or GMP-redo exist and there are at least one afterward-time point
      if((last_gmp_idx != -Inf) && (ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) - last_gmp_idx > 1)) {
        ### collect persistent clones that appeared in GMP and persist afterwards
        if(nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) > 0) {
          for(j in 1:nrow(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])) {
            for(k in last_gmp_idx:(ncol(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])-1)) {
              if(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,"GMP"] > 0 &&
                 SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]][j,k] > 0) {
                pClones <- c(pClones, rownames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]])[j])
                break;
              }
            }
          }
        }
      }
    }
    
    ### what are the patients that have less than 3 after-infusion time points?
    total_idx <- which(colnames(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]][[i]]) == "Total")
    if(total_idx - last_gmp_idx <= 3) {
      indeterminate_patient_pool <- c(indeterminate_patient_pool, names(SJCAR19_Clonotype_Frequency[["CARPOSONLY"]])[i])
    }
  }
  
  ### All the CAR+ persistent cells
  pIdx <- intersect(which(Seurat_Obj@meta.data$clonotype_id_by_patient %in% pClones),
                    which(Seurat_Obj@meta.data$CAR == "CARpos"))
  
  ### All the CAR+ non-persistent cells
  ### some patients have no later-time point data, so we cannot say they are non-persisters
  ### even if there is no lineage found - e.g., Px00 & Px01
  ### there should be at least 3 after-infusion time points to determine non-persisters
  npIdx <- setdiff(which(Seurat_Obj@meta.data$CAR == "CARpos"),
                   which(Seurat_Obj@meta.data$clonotype_id_by_patient %in% pClones))
  
  ### what are the patients that have less than 3 after-infusion time points?
  npIdx <- setdiff(npIdx, which(Seurat_Obj@meta.data$px %in% indeterminate_patient_pool))
  
  ### remove cells which do not have TCR info (we don't know about those cells yet)
  npIdx <- setdiff(npIdx, which(is.na(Seurat_Obj@meta.data$cdr3_aa)))
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
  
  ### annotate GMP CAR+ persisters
  Seurat_Obj@meta.data$ALL_CARpos_Persister <- NA
  Seurat_Obj@meta.data$ALL_CARpos_Persister[pIdx] <- "YES"
  Seurat_Obj@meta.data$ALL_CARpos_Persister[npIdx] <- "NO"
  
  
  ### set idents with the libary
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$ALL_CARpos_Persister)
  
  ### temporary seurat object for umap
  umap_seurat_obj <- subset(Seurat_Obj, idents = c("YES", "NO"))
  
  ### UMAP with GMP
  p <- DimPlot(object = umap_seurat_obj, reduction = "umap",
               group.by = "ALL_CARpos_Persister", split.by = NULL,
               pt.size = 2) +
    ggtitle("UMAP of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir, "/", "UMAP_Plot_Persistency1.png"), plot = p, width = 20, height = 12, dpi = 300)
  
  ### UMAP with GMP by each patient
  p <- DimPlot(object = umap_seurat_obj, reduction = "umap",
               group.by = "ALL_CARpos_Persister", split.by = "px",
               pt.size = 2, ncol = 3) +
    ggtitle("UMAP of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir, "/", "UMAP_Plot_Persistency2.png"), plot = p, width = 20, height = 12, dpi = 300)
  
  ### remove the temporary object
  rm(umap_seurat_obj)
  gc()
  
  ### DE analysis
  de_result <- FindMarkers(Seurat_Obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/GMP_CARpos_Persisters_vs_NonPersisters.xlsx"),
              sheetName = "GMP_CARpos_DE_Result", row.names = FALSE)
  
  ### pathway analysis
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                            rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                            "ENTREZID", "SYMBOL"),
                                          org = "human", database = "GO",
                                          title = paste0("Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters"),
                                          displayNum = 30, imgPrint = TRUE,
                                          dir = paste0(outputDir))
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters"),
                                            displayNum = 30, imgPrint = TRUE,
                                            dir = paste0(outputDir))
  if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
    write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters.xlsx"),
                row.names = FALSE, sheetName = paste0("GO_Result"))
  }
  if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
    write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_Pathway_Result_GMP_CARpos_Persisters_vs_NonPersisters.xlsx"),
                row.names = FALSE, sheetName = paste0("KEGG_Result"))
  }
  
  ### A function for scaling for heatmap
  scale_h <- function(data, type, na.rm=TRUE) {
    
    if(type == "row") {
      scaled <- t(scale(t(data)))
    } else if(type == "col") {
      scaled <- scale(data)
    } else {
      stop("Type is required: row or col")
    }
    
    if(na.rm == TRUE && (length(which(is.na(scaled))) > 0))  {
      scaled <- scaled[-unique(which(is.na(scaled), arr.ind = TRUE)[,1]),]
    }
    
    return(scaled)
  }
  
  ### hierarchical clustering functions
  dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
  hclust.ave <- function(x) hclust(x, method="average")
  
  #
  ### make a classifier with the DE genes
  #
  ### create outputDir2
  outputDir2 <- paste0(outputDir, "DE_Classifier/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### because of the imbalance of the two cluster sizes, we randomly choose
  ### the same number of samples in each class and iteratively build the classifier
  iteration <- 10
  set.seed(2990)
  featureSelectionNum <- 100
  sampleNum <- 100
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "parRF", "glmboost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "Linear_Model", "K-NN")
  train_control <- trainControl(method="LOOCV", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO")
  
  ### a small number that will be added to the counts for normalization
  ### log transform should not have 0 values
  log_trans_add <- 1
  
  ### build the classifier 10 times
  for(i in 1:iteration) {
    
    ### get random samples for each condition
    cond1_samps <- rownames(Seurat_Obj@meta.data)[sample(all_gmp_last, sampleNum)]
    cond2_samps <- rownames(Seurat_Obj@meta.data)[sample(all_gmp_not_last, sampleNum)]
    
    ### normalize the read counts
    ### before the normalization, only keep the samples that will be used in the classifier
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],
                                                                                             c(cond1_samps,
                                                                                               cond2_samps)] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE),
                                         filter_thresh = 0)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(c(rep("GMP_Last", sampleNum),
                                 rep("GMP_Not_Last", sampleNum)),
                               levels = c("GMP_Last", "GMP_Not_Last"))
    
    ### build classifier and test
    ### LOOCV
    p <- list()
    acc <- NULL
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=input_data, trControl=train_control, method=methodTypes[j])
      roc <- roc(model$pred$obs, model$pred$GMP_Last)
      acc <- c(acc, round(mean(model$results$Accuracy), 3))
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using DE Genes\n",
                                           "Accuracy =", acc[j]),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      gc()
    }
    
    ### draw ROC curves
    png(paste0(outputDir2, "Classifier_Using_DEG_GMP_Last_vs_Not_Last_", featureSelectionNum, "_(", i, ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", acc[j]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
    
    
    #
    ### Classifier using random genes
    #
    ### normalize the read counts
    ### before the normalization, only keep the samples that will be used in the classifier
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[sample(rownames(Seurat_Obj@assays$RNA@counts), featureSelectionNum),
                                                                                             c(cond1_samps,
                                                                                               cond2_samps)] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE),
                                         filter_thresh = 0)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(c(rep("GMP_Last", sampleNum),
                                 rep("GMP_Not_Last", sampleNum)),
                               levels = c("GMP_Last", "GMP_Not_Last"))
    
    ### build classifier and test
    ### LOOCV
    p <- list()
    acc <- NULL
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=input_data, trControl=train_control, method=methodTypes[j])
      roc <- roc(model$pred$obs, model$pred$GMP_Last)
      acc <- c(acc, round(mean(model$results$Accuracy), 3))
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using DE Genes\n",
                                           "Accuracy =", acc[j]),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      gc()
    }
    
    ### draw ROC curves
    png(paste0(outputDir2, "Classifier_Using_RG_GMP_Last_vs_Not_Last_", featureSelectionNum, "_(", i, ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", acc[j]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
    
    #
    ### A heatmap to show those genes are actually differentially expressed between the two conditions
    #
    # ### make a matrix for the heatmap
    # heatmap_mat <- as.matrix(normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],
    #                                                                                                     c(cond1_samps,
    #                                                                                                       cond2_samps)] + log_trans_add,
    #                                                                        stringsAsFactors = FALSE, check.names = FALSE)),
    #                                                 filter_trash = 0)
    
    ### make a matrix for the heatmap
    heatmap_mat <- as.matrix(scale_h(data.frame(Seurat_Obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],
                                                                             c(cond1_samps,
                                                                               cond2_samps)],
                                                                           stringsAsFactors = FALSE, check.names = FALSE),
                                     type = "row"))
    
    ### set rowside colors
    uniqueV <- c("YES", "NO")
    colors1 <- c("gray", "black")
    names(colors1) <- uniqueV
    uniqueV <- unique(Seurat_Obj@meta.data[colnames(heatmap_mat),"px"])
    colors2 <- colorRampPalette(brewer.pal(length(uniqueV), "Spectral"))(length(uniqueV))
    names(colors2) <- uniqueV
    
    ### order the colors2
    colors2 <- colors2[order(names(colors2))]
    
    ### because there are some outliers in positive values
    ### we scale more for the heatmap
    median_v <- median(heatmap_mat)
    min_v <- min(heatmap_mat)
    max_v <- max(heatmap_mat)
    if((median_v - min_v) > (max_v - median_v)) {
      imit_v <- median_v - (max_v - median_v)
      heatmap_mat[which(heatmap_mat < limit_v)] <- limit_v
    } else {
      limit_v <- median_v + (median_v - min_v)
      heatmap_mat[which(heatmap_mat > limit_v)] <- limit_v
    }
    
    ### heatmap
    png(paste0(outputDir2, "Heatmap_DEG_GMP_Last_vs_Not_Last_", featureSelectionNum, "_(", i, ").png"),
        width = 3000, height = 2200, res = 300)
    par(oma=c(0,0,3,5))
    heatmap.3(as.matrix(heatmap_mat),
              xlab = "", ylab = "", col=colorpanel(100, low = "green", high = "red"),
              scale="none", key=TRUE, keysize=1.5, density.info="density",
              dendrogram = "none", trace = "none",
              labRow = rownames(heatmap_mat), labCol = FALSE,
              Rowv = FALSE, Colv = FALSE,
              distfun=dist.spear, hclustfun=hclust.ave,
              ColSideColors = cbind(colors1[as.character(Seurat_Obj@meta.data[colnames(heatmap_mat),"GMP_CARpos_Persister"])],
                                    colors2[as.character(Seurat_Obj@meta.data[colnames(heatmap_mat),"px"])]),
              cexRow = 0.3, cexCol = 1.5, na.rm = TRUE,
              main = paste0(nrow(heatmap_mat), " DE Genes x ",
                            ncol(heatmap_mat), " Cells"))
    legend("topright", inset = 0, xpd = TRUE, title = "Persistency",
           legend = c("Persisters", "Non-Persisters"),
           fill = colors1, cex = 0.8, box.lty = 1)
    legend("left", inset = 0, xpd = TRUE, title = "Patients",
           legend = names(colors2),
           fill = colors2, cex = 0.8, box.lty = 1)
    dev.off()
    
    
    gc()
  }
  
  ### select 10 from each patient and build the classifier & the heatmap
  sample_num_from_each_px <- 10
  set.seed(2990)
  cond1_samps <- NULL
  cond2_samps <- NULL
  for(px in unique(Seurat_Obj@meta.data$px[all_gmp_last])) {
    px_idx <- which(Seurat_Obj@meta.data$px == px)
    cond1_samps <- c(cond1_samps, rownames(Seurat_Obj@meta.data)[sample(intersect(all_gmp_last,
                                                                                  px_idx),
                                                                        sample_num_from_each_px)])
    cond2_samps <- c(cond2_samps, rownames(Seurat_Obj@meta.data)[sample(intersect(all_gmp_not_last,
                                                                                  px_idx),
                                                                        sample_num_from_each_px)])
  }
  
  ### normalize the read counts
  ### before the normalization, only keep the samples that will be used in the classifier
  input_data <- normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],
                                                                                           c(cond1_samps,
                                                                                             cond2_samps)] + log_trans_add,
                                                              stringsAsFactors = FALSE, check.names = FALSE),
                                       filter_thresh = 0)
  
  ### annotate class for the input data
  input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
  input_data$Class <- factor(c(rep("GMP_Last", length(cond1_samps)),
                               rep("GMP_Not_Last", length(cond2_samps))),
                             levels = c("GMP_Last", "GMP_Not_Last"))
  
  ### build classifier and test
  ### LOOCV
  p <- list()
  acc <- NULL
  for(j in 1:length(methodTypes)) {
    writeLines(paste(methodTypes[j]))
    model <- train(Class~., data=input_data, trControl=train_control, method=methodTypes[j])
    roc <- roc(model$pred$obs, model$pred$GMP_Last)
    acc <- c(acc, round(mean(model$results$Accuracy), 3))
    p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "With DEGs & Selected Samples\n",
                                         "Accuracy =", acc[j]),
                       legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                       xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    gc()
  }
  
  ### draw ROC curves
  png(paste0(outputDir2, "Classifier_Using_DEG_GMP_Last_vs_Not_Last_With_Selected_Samples_", featureSelectionNum, ".png"),
      width = 2000, height = 2000, res = 350)
  par(mfrow=c(3, 2))
  for(j in 1:length(methodTypes)) {
    plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                  "Accuracy =", acc[j]),
             legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
             xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
  }
  dev.off()
  
  
  ### make a matrix for the heatmap
  heatmap_mat <- as.matrix(scale_h(data.frame(Seurat_Obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],
                                                                           c(cond1_samps,
                                                                             cond2_samps)],
                                              stringsAsFactors = FALSE, check.names = FALSE),
                                   type = "row"))
  
  ### set rowside colors
  uniqueV <- c("YES", "NO")
  colors1 <- c("gray", "black")
  names(colors1) <- uniqueV
  uniqueV <- unique(Seurat_Obj@meta.data[colnames(heatmap_mat),"px"])
  colors2 <- colorRampPalette(brewer.pal(length(uniqueV), "Spectral"))(length(uniqueV))
  names(colors2) <- uniqueV
  
  ### order the colors2
  colors2 <- colors2[order(names(colors2))]
  
  ### because there are some outliers in positive values
  ### we scale more for the heatmap
  median_v <- median(heatmap_mat)
  min_v <- min(heatmap_mat)
  max_v <- max(heatmap_mat)
  if((median_v - min_v) > (max_v - median_v)) {
    imit_v <- median_v - (max_v - median_v)
    heatmap_mat[which(heatmap_mat < limit_v)] <- limit_v
  } else {
    limit_v <- median_v + (median_v - min_v)
    heatmap_mat[which(heatmap_mat > limit_v)] <- limit_v
  }
  
  ### heatmap
  png(paste0(outputDir2, "Heatmap_DEG_GMP_Last_vs_Not_Last_With_Selected_Samples_", featureSelectionNum, ".png"),
      width = 3000, height = 2200, res = 300)
  par(oma=c(0,0,3,5))
  heatmap.3(as.matrix(heatmap_mat),
            xlab = "", ylab = "", col=colorpanel(100, low = "green", high = "red"),
            scale="none", key=TRUE, keysize=1.5, density.info="density",
            dendrogram = "none", trace = "none",
            labRow = rownames(heatmap_mat), labCol = FALSE,
            Rowv = FALSE, Colv = FALSE,
            distfun=dist.spear, hclustfun=hclust.ave,
            ColSideColors = cbind(colors1[as.character(Seurat_Obj@meta.data[colnames(heatmap_mat),"GMP_CARpos_Persister"])],
                                  colors2[as.character(Seurat_Obj@meta.data[colnames(heatmap_mat),"px"])]),
            cexRow = 0.3, cexCol = 1.5, na.rm = TRUE,
            main = paste0(nrow(heatmap_mat), " DE Genes x ",
                          ncol(heatmap_mat), " Cells"))
  legend("topright", inset = 0, xpd = TRUE, title = "Persistency",
         legend = c("Persisters", "Non-Persisters"),
         fill = colors1, cex = 0.8, box.lty = 1)
  legend("left", inset = 0, xpd = TRUE, title = "Patients",
         legend = names(colors2),
         fill = colors2, cex = 0.8, box.lty = 1)
  dev.off()
  
  #
  ### Reversed UMAP (gene-patient conversion) to show whether those genes are separated by persistency
  #
  
  ### get random samples for each condition
  set.seed(1234)
  persister_samps <- rownames(Seurat_Obj@meta.data)[sample(all_gmp_last, 1000)]
  DEGs <- rownames(de_result)[1:featureSelectionNum]
  RGs <- sample(rownames(Seurat_Obj@assays$RNA@counts), featureSelectionNum)
  
  ### Reversed UMAP
  input_data <- normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[c(DEGs, RGs),
                                                                                           persister_samps] + log_trans_add,
                                                              stringsAsFactors = FALSE, check.names = FALSE),
                                       filter_thresh = 0)
  reversed_umap <- umap(input_data)
  rownames(reversed_umap) <- c(DEGs, RGs)
  colnames(reversed_umap) <- c("UMAP1", "UMAP2")
  reversed_umap <- data.frame(reversed_umap,
                              Cluster=c(rep("DEGs", length(DEGs)), rep("RGs", length(RGs))),
                              Gene=rownames(reversed_umap),
                              stringsAsFactors = FALSE, check.names = FALSE)
  
  ### draw UMAP
  ggplot(reversed_umap, aes_string(x="UMAP1", y="UMAP2", label="Gene")) +
    geom_point(aes_string(col="Cluster"), size=2, alpha=0.8) +
    xlab("UMAP1") + ylab("UMAP2") +
    labs(col="DEGs/Random Genes") +
    ggtitle(paste0("Reversed UMAP with 100 DEGs & 100 RGs")) +
    scale_color_brewer(palette="Dark2") +
    theme_classic(base_size = 16)
  ggsave(file = paste0(outputDir2, "Reversed_UMAP_DEG_vs_RG.png"),
         width = 8, height = 6, dpi = 500)
  
  
  # ### get GMP CAR+ cells as a SingleCellExperiment object
  # GMP_CARpos_sce <- as.SingleCellExperiment(subset(Seurat_Obj, idents = c("YES", "NO")), assay = "RNA")
  # 
  # ### for memory efficiency, remove the large Seurat object
  # rm(Seurat_Obj)
  # gc()
  # 
  # ### QC
  # ### 26027 genes & 109614 cells
  # row_sum <- rowSums(counts(GMP_CARpos_sce) > 0)
  # col_sum <- colSums(counts(GMP_CARpos_sce) > 0)
  # plot(density(row_sum))
  # plot(density(col_sum))
  # 
  # ### set parameter
  # ### gene.min.pct: in a gene, at least n % of each condition (cells) should be expressed
  # ### cell.min.num: in a cell, at least n genes should be expressed (determine after looking at the density plot above)
  # gene.min.pct <- 0.2
  # cell.min.num <- 2000
  # cell.max.num <- 5000
  # 
  # ### low count filter for genes
  # yes_num <- length(which(GMP_CARpos_sce$ident == "YES"))
  # no_num <- length(which(GMP_CARpos_sce$ident == "NO"))
  # yes_sum <- rowSums(counts(GMP_CARpos_sce)[,which(GMP_CARpos_sce$ident == "YES")] > 0)
  # no_sum <- rowSums(counts(GMP_CARpos_sce)[,which(GMP_CARpos_sce$ident == "NO")] > 0)
  # yes_exp_pct <- yes_sum / yes_num
  # no_exp_pct <- no_sum / no_num
  # keep <- intersect(which(yes_exp_pct >= gene.min.pct),
  #                   which(no_exp_pct >= gene.min.pct))
  # GMP_CARpos_sce <- GMP_CARpos_sce[keep,]
  # 
  # ### low count filter for cells
  # keep <- intersect(which(colSums(counts(GMP_CARpos_sce) > 0) >= cell.min.num),
  #                   which(colSums(counts(GMP_CARpos_sce) > 0) <= cell.max.num))
  # GMP_CARpos_sce <- GMP_CARpos_sce[,keep]
  # 
  # 
  # ### there are too many cells (especially in the second condition), so we will randomly choose some
  # ### to lower the computational complexity
  # set.seed(1234)
  # keep <- sample(which(GMP_CARpos_sce$ident == "NO"), length(which(GMP_CARpos_sce$ident == "YES")))
  # keep <- union(keep, which(GMP_CARpos_sce$ident == "YES"))
  # GMP_CARpos_sce <- GMP_CARpos_sce[,keep]
  # 
  # ### "counts" should be in the first place in the assayNames(GMP_CARpos_sce)
  # nms <- c("counts", setdiff(assayNames(GMP_CARpos_sce), "counts"))
  # assays(GMP_CARpos_sce) <- assays(GMP_CARpos_sce)[nms]
  # 
  # ### epsilon setting as recommended by the ZINB-WaVE integration paper
  # system.time({
  #   zinb <- zinbwave(GMP_CARpos_sce, K=0, observationalWeights=TRUE,
  #                    BPPARAM=SerialParam(), epsilon=1e12)
  # })
  # 
  # ### DE analysis
  # dds <- DESeqDataSet(zinb, design=~condition)
  # dds <- estimateSizeFactors(dds, type="poscounts")
  # library(scran)
  # scr <- computeSumFactors(dds)
  # dat <- data.frame(true=dds$ExpLibSize,
  #                   pos=sizeFactors(dds),
  #                   sum=sizeFactors(scr))
  # dat$true <- dat$true / exp(mean(log(dat$true)))
  # panel.scatter <- function(x,y,...) {
  #   points(x,y,...)
  #   abline(0,1,col="red",lwd=2)
  #   legend("topleft", legend=round(cor(x,y),3))
  # }
  # pairs(dat, panel=panel.scatter)
  
  
  ### Make statistics
  ### 1.  The number of the lineages
  ### 2.  The number of CD4 lineages
  ### 3.  The number of CD8 lineages
  ### 4.  The number of the CAR+ lineages
  ### 5. The number of the CAR- lineages
  ### 6.  The number of CD4 cells in the CD4 lineages
  ### 7.  The number of CD8 cells in the CD8 lineages
  ### 8.  The number of CAR+ cells in the CAR+ lineages
  ### 9.  The number of CAR- cells in the CAR- lineages
  ### 10. The number of total cells in the lineages
  ### 11. The number of lineages from GMP
  ### 12. The number of cells from GMP
  ### 13. The number of CD4 lineages from GMP
  ### 14. The number of CD8 lineages from GMP
  ### 15. The number of CAR+ lineages from GMP 
  ### 16. The number of CAR- lineages from GMP
  ### 17. The number of CD4 cells in the CD4 lineages from GMP
  ### 18. The number of CD8 cells in the CD8 lineages from GMP
  ### 19. The number of CAR+ cells in the CAR+ lineages from GMP
  ### 20. The number of CAR- cells in the CAR- lineages from GMP
  ###
  ### rows are patients and the 12 columns are the above info
  
  ### load full lineage info
  SJCAR19_Lineages_in_Full <- readRDS(lineage_in_full_path)
  
  ### make an empty data frame
  SJCAR19_Stats_Table <- data.frame(matrix(NA,
                                           length(SJCAR19_Clonotype_Frequency[["ALL"]]),
                                           20),
                                    stringsAsFactors = FALSE, check.names = FALSE)
  rownames(SJCAR19_Stats_Table) <- names(SJCAR19_Clonotype_Frequency[["ALL"]])
  colnames(SJCAR19_Stats_Table) <- c("Lineage_#",
                                     "CD4_Lineage_#",
                                     "CD8_Lineage_#",
                                     "CARpos_Lineage_#",
                                     "CARneg_Lineage_#",
                                     "Lineage_CD4_Cell_#",
                                     "Lineage_CD8_Cell_#",
                                     "Lineage_CARpos_Cell_#",
                                     "Lineage_CARneg_Cell_#",
                                     "Lineage_Total_Cell_#",
                                     "GMP_Persister_Lineage_#",
                                     "GMP_Persister_Cell_#",
                                     "GMP_Persister_CD4_Lineage_#",
                                     "GMP_Persister_CD8_Lineage_#",
                                     "GMP_Persister_CARpos_Lineage_#",
                                     "GMP_Persister_CARneg_Lineage_#",
                                     "GMP_Persister_CD4_Cell_#",
                                     "GMP_Persister_CD8_Cell_#",
                                     "GMP_Persister_CARpos_Cell_#",
                                     "GMP_Persister_CARneg_Cell_#")
  
  ### for each patient fill out the stats table
  for(px in rownames(SJCAR19_Stats_Table)) {
    ### print progress
    writeLines(paste(px))
    
    ### Lineage_#
    target_table <- SJCAR19_Clonotype_Frequency[["ALL"]][[px]]
    if(ncol(target_table) > 2) {
      lineage_idx <- which(apply(target_table[,1:(ncol(target_table)-1)], 1, function(x) length(which(x > 0))) > 1)
      SJCAR19_Stats_Table[px,"Lineage_#"] <- length(lineage_idx)
      
      ### Lineage_Total_Cell_#
      SJCAR19_Stats_Table[px,"Lineage_Total_Cell_#"] <- sum(as.numeric(unlist(target_table[lineage_idx,"Total"])))
      
      ### get index after GMP
      gmp_idx <- which(colnames(target_table) == "GMP")
      gmp_redo_idx <- which(colnames(target_table) == "GMP-redo")
      after_gmp_idx <- max(gmp_idx, gmp_redo_idx) + 1
      
      ### GMP_Persister_Lineage_#
      ### GMP_Persister_Cell_#
      if((after_gmp_idx != -Inf) && (after_gmp_idx < ncol(target_table))) {
        if(length(gmp_idx) > 0 && length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- (target_table[,gmp_idx] > 0) | (target_table[,gmp_redo_idx] > 0)
        } else if(length(gmp_idx) > 0) {
          is_gmp_cell_exist <- target_table[,gmp_idx] > 0
        } else if(length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- target_table[,gmp_redo_idx] > 0
        } else {
          stop("ERROR")
        }
        is_after_gmp_exist <- apply(target_table[,after_gmp_idx:(ncol(target_table)-1),drop=FALSE], 1, function(x) sum(x) > 0)
        is_gmp_lineage_exist <- is_gmp_cell_exist & is_after_gmp_exist
        
        SJCAR19_Stats_Table[px,"GMP_Persister_Lineage_#"] <- length(which(is_gmp_lineage_exist == TRUE))
        SJCAR19_Stats_Table[px,"GMP_Persister_Cell_#"] <- sum(as.numeric(unlist(target_table[is_gmp_lineage_exist,after_gmp_idx:(ncol(target_table)-1)])))
      } else {
        SJCAR19_Stats_Table[px,"GMP_Persister_Lineage_#"] <- 0
        SJCAR19_Stats_Table[px,"GMP_Persister_Cell_#"] <- 0
      }
      
      ### CD4_Lineage_#
      target_table_cd4 <- SJCAR19_Lineages_in_Full[[px]][intersect(which(SJCAR19_Lineages_in_Full[[px]]$CD_Type == "CD4"),
                                                                   which(SJCAR19_Lineages_in_Full[[px]]$CAR_Type == "ALL")),
                                                         (which(colnames(SJCAR19_Lineages_in_Full[[px]]) == "CDR3_NT")+1):which(colnames(SJCAR19_Lineages_in_Full[[px]]) == "Total")]
      CD4_lineage_idx <- which(apply(target_table_cd4[,1:(ncol(target_table_cd4)-1)], 1, function(x) length(which(x > 0))) > 1)
      SJCAR19_Stats_Table[px,"CD4_Lineage_#"] <- length(CD4_lineage_idx)
      
      ### CD8_Lineage_#
      target_table_cd8 <- SJCAR19_Lineages_in_Full[[px]][intersect(which(SJCAR19_Lineages_in_Full[[px]]$CD_Type == "CD8"),
                                                                   which(SJCAR19_Lineages_in_Full[[px]]$CAR_Type == "ALL")),
                                                         (which(colnames(SJCAR19_Lineages_in_Full[[px]]) == "CDR3_NT")+1):which(colnames(SJCAR19_Lineages_in_Full[[px]]) == "Total")]
      CD8_lineage_idx <- which(apply(target_table_cd8[,1:(ncol(target_table_cd8)-1)], 1, function(x) length(which(x > 0))) > 1)
      SJCAR19_Stats_Table[px,"CD8_Lineage_#"] <- length(CD8_lineage_idx)
      
      ### Lineage_CD4_Cell_#
      SJCAR19_Stats_Table[px,"Lineage_CD4_Cell_#"] <- sum(as.numeric(unlist(target_table_cd4[CD4_lineage_idx,"Total"])))
      
      ### Lineage_CD8_Cell_#
      SJCAR19_Stats_Table[px,"Lineage_CD8_Cell_#"] <- sum(as.numeric(unlist(target_table_cd8[CD8_lineage_idx,"Total"])))
      
      ### CARpos_Lineage_#
      target_table_carpos <- SJCAR19_Lineages_in_Full[[px]][intersect(which(SJCAR19_Lineages_in_Full[[px]]$CAR_Type == "CARpos"),
                                                                      which(SJCAR19_Lineages_in_Full[[px]]$CD_Type == "ALL")),
                                                            (which(colnames(SJCAR19_Lineages_in_Full[[px]]) == "CDR3_NT")+1):which(colnames(SJCAR19_Lineages_in_Full[[px]]) == "Total")]
      carpos_lineage_idx <- which(apply(target_table_carpos[,1:(ncol(target_table_carpos)-1)], 1, function(x) length(which(x > 0))) > 1)
      SJCAR19_Stats_Table[px,"CARpos_Lineage_#"] <- length(carpos_lineage_idx)
      
      ### CARneg_Lineage_#
      target_table_carneg <- SJCAR19_Lineages_in_Full[[px]][intersect(which(SJCAR19_Lineages_in_Full[[px]]$CAR_Type == "CARneg"),
                                                                      which(SJCAR19_Lineages_in_Full[[px]]$CD_Type == "ALL")),
                                                            (which(colnames(SJCAR19_Lineages_in_Full[[px]]) == "CDR3_NT")+1):which(colnames(SJCAR19_Lineages_in_Full[[px]]) == "Total")]
      carneg_lineage_idx <- which(apply(target_table_carneg[,1:(ncol(target_table_carneg)-1)], 1, function(x) length(which(x > 0))) > 1)
      SJCAR19_Stats_Table[px,"CARneg_Lineage_#"] <- length(carneg_lineage_idx)
      
      ### Lineage_CARpos_Cell_#
      SJCAR19_Stats_Table[px,"Lineage_CARpos_Cell_#"] <- sum(as.numeric(unlist(target_table_carpos[carpos_lineage_idx,"Total"])))
      
      ### Lineage_CARneg_Cell_#
      SJCAR19_Stats_Table[px,"Lineage_CARneg_Cell_#"] <- sum(as.numeric(unlist(target_table_carneg[carneg_lineage_idx,"Total"])))
      
      #
      ### GMP_Persister_CD4_Lineage_#
      ### GMP_Persister_CD4_Cell_#
      #
      ### get index after GMP
      gmp_idx <- which(colnames(target_table_cd4) == "GMP")
      gmp_redo_idx <- which(colnames(target_table_cd4) == "GMP-redo")
      after_gmp_idx <- max(gmp_idx, gmp_redo_idx) + 1
      
      if((after_gmp_idx != -Inf) && (after_gmp_idx < ncol(target_table_cd4))) {
        if(length(gmp_idx) > 0 && length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- (target_table_cd4[,gmp_idx] > 0) | (target_table_cd4[,gmp_redo_idx] > 0)
        } else if(length(gmp_idx) > 0) {
          is_gmp_cell_exist <- target_table_cd4[,gmp_idx] > 0
        } else if(length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- target_table_cd4[,gmp_redo_idx] > 0
        } else {
          stop("ERROR")
        }
        is_after_gmp_exist <- apply(target_table_cd4[,after_gmp_idx:(ncol(target_table_cd4)-1),drop=FALSE], 1, function(x) sum(x) > 0)
        is_gmp_lineage_exist <- is_gmp_cell_exist & is_after_gmp_exist
        
        SJCAR19_Stats_Table[px,"GMP_Persister_CD4_Lineage_#"] <- length(which(is_gmp_lineage_exist == TRUE))
        SJCAR19_Stats_Table[px,"GMP_Persister_CD4_Cell_#"] <- sum(as.numeric(unlist(target_table_cd4[is_gmp_lineage_exist,after_gmp_idx:(ncol(target_table_cd4)-1)])))
      } else {
        SJCAR19_Stats_Table[px,"GMP_Persister_CD4_Lineage_#"] <- 0
        SJCAR19_Stats_Table[px,"GMP_Persister_CD4_Cell_#"] <- 0
      }
      
      #
      ### GMP_Persister_CD8_Lineage_#
      ### GMP_Persister_CD8_Cell_#
      #
      ### get index after GMP
      gmp_idx <- which(colnames(target_table_cd8) == "GMP")
      gmp_redo_idx <- which(colnames(target_table_cd8) == "GMP-redo")
      after_gmp_idx <- max(gmp_idx, gmp_redo_idx) + 1
      
      if((after_gmp_idx != -Inf) && (after_gmp_idx < ncol(target_table_cd8))) {
        if(length(gmp_idx) > 0 && length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- (target_table_cd8[,gmp_idx] > 0) | (target_table_cd8[,gmp_redo_idx] > 0)
        } else if(length(gmp_idx) > 0) {
          is_gmp_cell_exist <- target_table_cd8[,gmp_idx] > 0
        } else if(length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- target_table_cd8[,gmp_redo_idx] > 0
        } else {
          stop("ERROR")
        }
        is_after_gmp_exist <- apply(target_table_cd8[,after_gmp_idx:(ncol(target_table_cd8)-1),drop=FALSE], 1, function(x) sum(x) > 0)
        is_gmp_lineage_exist <- is_gmp_cell_exist & is_after_gmp_exist
        
        SJCAR19_Stats_Table[px,"GMP_Persister_CD8_Lineage_#"] <- length(which(is_gmp_lineage_exist == TRUE))
        SJCAR19_Stats_Table[px,"GMP_Persister_CD8_Cell_#"] <- sum(as.numeric(unlist(target_table_cd8[is_gmp_lineage_exist,after_gmp_idx:(ncol(target_table_cd8)-1)])))
      } else {
        SJCAR19_Stats_Table[px,"GMP_Persister_CD8_Lineage_#"] <- 0
        SJCAR19_Stats_Table[px,"GMP_Persister_CD8_Cell_#"] <- 0
      }
      
      #
      ### GMP_Persister_CARpos_Lineage_#
      ### GMP_Persister_CARpos_Cell_#
      #
      ### get index after GMP
      gmp_idx <- which(colnames(target_table_carpos) == "GMP")
      gmp_redo_idx <- which(colnames(target_table_carpos) == "GMP-redo")
      after_gmp_idx <- max(gmp_idx, gmp_redo_idx) + 1
      
      if((after_gmp_idx != -Inf) && (after_gmp_idx < ncol(target_table_carpos))) {
        if(length(gmp_idx) > 0 && length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- (target_table_carpos[,gmp_idx] > 0) | (target_table_carpos[,gmp_redo_idx] > 0)
        } else if(length(gmp_idx) > 0) {
          is_gmp_cell_exist <- target_table_carpos[,gmp_idx] > 0
        } else if(length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- target_table_carpos[,gmp_redo_idx] > 0
        } else {
          stop("ERROR")
        }
        is_after_gmp_exist <- apply(target_table_carpos[,after_gmp_idx:(ncol(target_table_carpos)-1),drop=FALSE], 1, function(x) sum(x) > 0)
        is_gmp_lineage_exist <- is_gmp_cell_exist & is_after_gmp_exist
        
        SJCAR19_Stats_Table[px,"GMP_Persister_CARpos_Lineage_#"] <- length(which(is_gmp_lineage_exist == TRUE))
        SJCAR19_Stats_Table[px,"GMP_Persister_CARpos_Cell_#"] <- sum(as.numeric(unlist(target_table_carpos[is_gmp_lineage_exist,after_gmp_idx:(ncol(target_table_carpos)-1)])))
      } else {
        SJCAR19_Stats_Table[px,"GMP_Persister_CARpos_Lineage_#"] <- 0
        SJCAR19_Stats_Table[px,"GMP_Persister_CARpos_Cell_#"] <- 0
      }
      
      #
      ### GMP_Persister_CARneg_Lineage_#
      ### GMP_Persister_CARneg_Cell_#
      #
      ### get index after GMP
      gmp_idx <- which(colnames(target_table_carneg) == "GMP")
      gmp_redo_idx <- which(colnames(target_table_carneg) == "GMP-redo")
      after_gmp_idx <- max(gmp_idx, gmp_redo_idx) + 1
      
      if((after_gmp_idx != -Inf) && (after_gmp_idx < ncol(target_table_carneg))) {
        if(length(gmp_idx) > 0 && length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- (target_table_carneg[,gmp_idx] > 0) | (target_table_carneg[,gmp_redo_idx] > 0)
        } else if(length(gmp_idx) > 0) {
          is_gmp_cell_exist <- target_table_carneg[,gmp_idx] > 0
        } else if(length(gmp_redo_idx) > 0) {
          is_gmp_cell_exist <- target_table_carneg[,gmp_redo_idx] > 0
        } else {
          stop("ERROR")
        }
        is_after_gmp_exist <- apply(target_table_carneg[,after_gmp_idx:(ncol(target_table_carneg)-1),drop=FALSE], 1, function(x) sum(x) > 0)
        is_gmp_lineage_exist <- is_gmp_cell_exist & is_after_gmp_exist
        
        SJCAR19_Stats_Table[px,"GMP_Persister_CARneg_Lineage_#"] <- length(which(is_gmp_lineage_exist == TRUE))
        SJCAR19_Stats_Table[px,"GMP_Persister_CARneg_Cell_#"] <- sum(as.numeric(unlist(target_table_carneg[is_gmp_lineage_exist,after_gmp_idx:(ncol(target_table_carneg)-1)])))
      } else {
        SJCAR19_Stats_Table[px,"GMP_Persister_CARneg_Lineage_#"] <- 0
        SJCAR19_Stats_Table[px,"GMP_Persister_CARneg_Cell_#"] <- 0
      }
    }
  }
  
  ### save the stats table
  output <- data.frame(Patient=rownames(SJCAR19_Stats_Table), SJCAR19_Stats_Table,
                       stringsAsFactors = FALSE, check.names = FALSE)
  output[which(is.na(output), arr.ind = TRUE)] <- 0
  write.xlsx2(output,
              file = paste0(outputDir, "SJCAR19_Statistics.xlsx"),
              sheetName = "SJCAR19_Statistics", row.names = FALSE)
  
  
  ### Classifier (Remove CD4 cells)
  ### Use only one random cell per each lineage and permute to see if there are many changes
  #
  ### gmp_last vs gmp_not_last
  ### 3428 vs 102958 cells
  ### 2126 vs 58114 CD8 cells
  ### 1223 vs 94240 lineages
  ### 837 vs 51671 CD8 lineages
  
  ### create outputDir2
  outputDir2 <- paste0(outputDir, "DE_Classifier_One_Cell_Per_Lineage/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### set parameters
  set.seed(1234)
  featureSelectionNum <- 100
  cv_k <- 10
  resampling_p <- 0.8
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "parRF", "glmboost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "Linear_Model", "KNN")
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  all_gmp_last <- intersect(all_gmp_last,
                            which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  all_gmp_not_last <- intersect(all_gmp_not_last,
                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### a table that indicates which cell comes from which lineage (persister vs non-persister)
  persister_cell_table <- data.frame(Persistance=c(rep("YES", length(all_gmp_last)),
                                                   rep("NO", length(all_gmp_not_last))),
                                     Cell_Name=rownames(Seurat_Obj@meta.data)[c(all_gmp_last, all_gmp_not_last)],
                                     Index=c(all_gmp_last, all_gmp_not_last),
                                     Clonotype=Seurat_Obj@meta.data$clonotype_id_by_patient[c(all_gmp_last, all_gmp_not_last)],
                                     stringsAsFactors = FALSE, check.names = FALSE)
  
  ### unique lineages
  all_gmp_last_lineages <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient[all_gmp_last])
  all_gmp_not_last_lineages <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient[all_gmp_not_last])
  
  ### a small number that will be added to the counts for normalization
  ### log transform should not have 0 values
  log_trans_add <- 1
  
  ### resampling numbers
  trainingNum <- round(length(all_gmp_last_lineages)*resampling_p)
  testNum <- length(all_gmp_last_lineages) - trainingNum
  
  ### performance evaluation
  eval_acc <- vector("list", length(methodTypes))
  names(eval_acc) <- methodNames
  eval_auc <- vector("list", length(methodTypes))
  names(eval_auc) <- methodNames
  
  ### resampling
  training_clonotypes <- vector("list", length = cv_k)
  test_clonotypes <- vector("list", length = cv_k)
  for(i in 1:length(training_clonotypes)) {
    training_clonotypes[[i]] <- c(sample(all_gmp_last_lineages, trainingNum), sample(all_gmp_not_last_lineages, trainingNum))
    test_clonotypes[[i]] <- c(sample(setdiff(all_gmp_last_lineages, training_clonotypes[[i]]), testNum),
                              sample(setdiff(all_gmp_not_last_lineages, training_clonotypes[[i]]), testNum))
  }
  
  ### get random samples from the sample lineages (one cell per lineage)
  training_samps <- vector("list", length = cv_k)
  test_samps <- vector("list", length = cv_k)
  for(i in 1:cv_k) {
    temp <- NULL
    for(lineage in training_clonotypes[[i]]) {
      lineage_pool <- which(persister_cell_table$Clonotype == lineage)
      ### if there one integer vector, sample(I, 1) == sample(1:I, 1) == random # between 1 & I, so be careful
      ### if one character vector, sample(C, 1) == C
      if(length(lineage_pool) > 1) {
        temp <- c(temp, persister_cell_table$Cell_Name[sample(lineage_pool, 1)])
      } else {
        temp <- c(temp, persister_cell_table$Cell_Name[lineage_pool])
      }
    }
    training_samps[[i]] <- temp
    temp <- NULL
    for(lineage in test_clonotypes[[i]]) {
      lineage_pool <- which(persister_cell_table$Clonotype == lineage)
      ### if there one integer vector, sample(I, 1) == sample(1:I, 1) == random # between 1 & I, so be careful
      ### if one character vector, sample(C, 1) == C
      if(length(lineage_pool) > 1) {
        temp <- c(temp, persister_cell_table$Cell_Name[sample(lineage_pool, 1)])
      } else {
        temp <- c(temp, persister_cell_table$Cell_Name[lineage_pool])
      }
    }
    test_samps[[i]] <- temp
  }
  
  ### build the classifier 10 times
  for(i in 1:cv_k) {
    
    ### wirte progress
    writeLines(paste(i))
    
    ### new obj for the training data & set idents with the info
    classifier_seurat_obj <- subset(Seurat_Obj, cells = training_samps[[i]])
    classifier_seurat_obj <- SetIdent(object = classifier_seurat_obj,
                                      cells = rownames(classifier_seurat_obj@meta.data),
                                      value = classifier_seurat_obj@meta.data$GMP_CARpos_Persister)
    
    ### DE analysis
    de_result <- FindMarkers(classifier_seurat_obj,
                             ident.1 = "YES",
                             ident.2 = "NO",
                             min.pct = 0.1,
                             logfc.threshold = 0.1,
                             test.use = "wilcox")
    
    ### normalize the read counts
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE),
                                         filter_thresh = 0)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(classifier_seurat_obj@meta.data$GMP_CARpos_Persister,
                               levels = c("YES", "NO"))
    
    ### new obj for the training data & set idents with the info
    classifier_seurat_obj <- subset(Seurat_Obj, cells = test_samps[[i]])
    
    ### normalize the read counts
    test_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],] + log_trans_add,
                                                               stringsAsFactors = FALSE, check.names = FALSE),
                                        filter_thresh = 0)
    
    ### annotate class for the test data
    test_data <- data.frame(t(test_data), stringsAsFactors = FALSE, check.names = FALSE)
    test_data$Class <- factor(classifier_seurat_obj@meta.data$GMP_CARpos_Persister,
                              levels = c("YES", "NO"))
    
    ### train control options
    train_control <- trainControl(method="none", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
    
    ### build classifier and test
    p <- list()
    acc <- NA
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=input_data, method=methodTypes[j], trControl = train_control)
      pred_result <- predict(model, newdata = test_data)
      acc <- round((sum(pred_result == test_data$Class) / nrow(test_data)), 3)
      pred_result <- predict(model, newdata = test_data, type = "prob")
      roc <- roc(test_data$Class, pred_result$YES)
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using DE Genes\n",
                                           "Accuracy =", acc),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      eval_acc[[methodNames[j]]] <- c(eval_acc[[methodNames[j]]], acc)
      eval_auc[[methodNames[j]]] <- c(eval_auc[[methodNames[j]]], as.numeric(roc$auc))
      gc()
    }
    
    ### draw ROC curves
    png(paste0(outputDir2, "Classifier_Using_DEG_GMP_Last_vs_Not_Last_", featureSelectionNum, "_One_Cell_Per_Lineage_(", i, ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", eval_acc[[methodNames[j]]][i]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
  
    gc()
  }
  
  ### save the ACC & AUC values
  saveRDS(eval_acc, file = paste0(outputDir2, "eval_acc.RDS"))
  saveRDS(eval_auc, file = paste0(outputDir2, "eval_auc.RDS"))
  
  ### Draw line graphs with iteration ACC & AUC
  plot_df <- matrix(0, cv_k*2, length(methodTypes)+2)
  colnames(plot_df) <- c("Iteration", methodNames, "Measure")
  plot_df <- data.frame(plot_df, stringsAsFactors = FALSE, check.names = FALSE)
  
  ### fill out the table
  plot_df$Iteration <- c(1:cv_k, 1:cv_k)
  for(mname in methodNames) {
    plot_df[1:cv_k,mname] <- eval_acc[[mname]]
    plot_df[(cv_k+1):(cv_k*2),mname] <- eval_auc[[mname]]
    plot_df[,"Measure"] <- c(rep("ACC", cv_k), rep("AUC", cv_k))
  }
  
  ### line graph generation
  p <- vector("list", length = length(methodTypes))
  names(p) <- methodNames
  for(mname in methodNames) {
    p[[mname]] <- ggplot(plot_df, aes_string(x= "Iteration", y=mname, group="Measure")) +
      geom_line(aes_string(color="Measure", linetype="Measure"), size=2) +
      geom_point() +
      geom_text(label = as.character(round(plot_df[,mname], digits = 2)),
                vjust = "inward", hjust = "inward") +
      ylim(c(0, 1)) +
      theme_classic(base_size = 16) +
      scale_color_npg() +
      ggtitle(mname) +
      theme(axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(breaks = seq(0, cv_k, by = 1))
  }
  
  ### arrange the plots and save
  fName <- paste0("Classifier_Result_Using_DEG_GMP_Last_vs_Not_Last_", featureSelectionNum, "_One_Cell_Per_Lineage")
  rowNum <- 3
  colNum <- 2
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir2, fName, ".png"), g, width = 17, height = 11, dpi = 300)
  
  
  #
  ### Classifier (Remove CD4 cells)
  ### Use all the cells of all the lineages but LOPOCV (Leave One Patient Out -> Patient-based)
  #
  ### create outputDir2
  outputDir2 <- paste0(outputDir, "DE_Classifier_LOPOCV/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  all_gmp_last <- intersect(all_gmp_last,
                            which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  all_gmp_not_last <- intersect(all_gmp_not_last,
                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### a table that indicates which cell comes from which lineage (persister vs non-persister)
  persister_cell_table <- data.frame(Persistance=c(rep("YES", length(all_gmp_last)),
                                                   rep("NO", length(all_gmp_not_last))),
                                     Cell_Name=rownames(Seurat_Obj@meta.data)[c(all_gmp_last, all_gmp_not_last)],
                                     Index=c(all_gmp_last, all_gmp_not_last),
                                     Clonotype=Seurat_Obj@meta.data$clonotype_id_by_patient[c(all_gmp_last, all_gmp_not_last)],
                                     stringsAsFactors = FALSE, check.names = FALSE)
  
  ### unique lineages
  all_gmp_last_lineages <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient[all_gmp_last])
  all_gmp_not_last_lineages <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient[all_gmp_not_last])
  
  ### because of the imbalance of the two cluster sizes, we randomly choose
  ### the same number of samples in each class and iteratively build the classifier
  iteration <- 10
  set.seed(2990)
  featureSelectionNum <- 100
  testSampleNum <- 1000
  target_px <- unique(intersect(Seurat_Obj@meta.data$px[all_gmp_last],
                                Seurat_Obj@meta.data$px[all_gmp_not_last]))
  cv_k <- length(target_px)
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "parRF", "glmboost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "Linear_Model", "KNN")
  log_trans_add <- 1
  
  ### performance evaluation
  eval_acc <- vector("list", length(methodTypes))
  names(eval_acc) <- methodNames
  eval_auc <- vector("list", length(methodTypes))
  names(eval_auc) <- methodNames
  
  ### prepare training & test samples
  training_samps <- vector("list", length = cv_k)
  test_samps <- vector("list", length = cv_k)
  for(i in 1:cv_k) {
    test_samps[[i]] <- rownames(Seurat_Obj@meta.data)[intersect(c(all_gmp_last, sample(all_gmp_not_last, length(all_gmp_last))),
                                                                which(Seurat_Obj@meta.data$px == target_px[i]))]
    if(length(test_samps[[i]]) > testSampleNum) {
      test_samps[[i]] <- sample(test_samps[[i]], testSampleNum)
    }
    ### because there are too many cells in Non-persisters, we select (5 x # test samples) random samples
    training_samps[[i]] <- sample(rownames(Seurat_Obj@meta.data)[setdiff(c(all_gmp_last, all_gmp_not_last),
                                                                         which(Seurat_Obj@meta.data$px == target_px[i]))],
                                  5*length(test_samps[[i]]))
  }
  
  ### build the classifier 10 times
  for(i in 1:cv_k) {
    
    ### wirte progress
    writeLines(paste(i))
    
    ### new obj for the training data & set idents with the info
    classifier_seurat_obj <- subset(Seurat_Obj, cells = training_samps[[i]])
    classifier_seurat_obj <- SetIdent(object = classifier_seurat_obj,
                                      cells = rownames(classifier_seurat_obj@meta.data),
                                      value = classifier_seurat_obj@meta.data$GMP_CARpos_Persister)
    
    ### DE analysis
    de_result <- FindMarkers(classifier_seurat_obj,
                             ident.1 = "YES",
                             ident.2 = "NO",
                             min.pct = 0.1,
                             logfc.threshold = 0.1,
                             test.use = "wilcox")
    
    ### normalize the read counts
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE),
                                         filter_thresh = 0)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(classifier_seurat_obj@meta.data$GMP_CARpos_Persister,
                               levels = c("YES", "NO"))
    
    ### new obj for the training data & set idents with the info
    classifier_seurat_obj <- subset(Seurat_Obj, cells = test_samps[[i]])
    
    ### normalize the read counts
    test_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],] + log_trans_add,
                                                               stringsAsFactors = FALSE, check.names = FALSE),
                                        filter_thresh = 0)
    
    ### annotate class for the test data
    test_data <- data.frame(t(test_data), stringsAsFactors = FALSE, check.names = FALSE)
    test_data$Class <- factor(classifier_seurat_obj@meta.data$GMP_CARpos_Persister,
                              levels = c("YES", "NO"))
    
    ### train control options
    train_control <- trainControl(method="none", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
    
    ### build classifier and test
    p <- list()
    acc <- NA
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=input_data, method=methodTypes[j], trControl = train_control)
      pred_result <- predict(model, newdata = test_data)
      acc <- round((sum(pred_result == test_data$Class) / nrow(test_data)), 3)
      pred_result <- predict(model, newdata = test_data, type = "prob")
      roc <- roc(test_data$Class, pred_result$YES)
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using DE Genes\n",
                                           "Accuracy =", acc),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      eval_acc[[methodNames[j]]] <- c(eval_acc[[methodNames[j]]], acc)
      eval_auc[[methodNames[j]]] <- c(eval_auc[[methodNames[j]]], as.numeric(roc$auc))
      gc()
    }
    
    ### draw ROC curves
    png(paste0(outputDir2, "Classifier_Using_DEG_GMP_Last_vs_Not_Last_", featureSelectionNum, "_LOPOCV_(", i, ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", eval_acc[[methodNames[j]]][i]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
    
    gc()
  }
  
  
  ### A function to perform 2D PCA and save a plot
  ### normalizedMat: rows are genes and columns are samples
  ### grp: group information of the samples - color
  ### grp2: group2 information of the samples - shape
  ### num: the number of top genes to be used based on variance (-1 [default]: use all the genes)
  ### component: to draw a plot with PC1 & PC2 or PC2 & PC3
  ### title: title of the plot
  ### suppliment: if it is TRUE, this function additionally generates figures and tables
  ###             related to contributions. For example, which genes are highly contributed to
  ###             the PC1, PC2, etc.
  ### outDir: output directory for the plot
  pca_plot <- function(normalizedMat, grp, grp2=NULL,
                       num = -1, component=c("PC1&PC2", "PC2&PC3"),
                       title="PCA_Plot",
                       suppliment=FALSE,
                       outDir="./") {
    ### load library
    if(!require(ggfortify, quietly = TRUE)) {
      install.packages("ggfortify")
      library(ggfortify, quietly = TRUE)
    }
    if(!require(FactoMineR, quietly = TRUE)) {
      install.packages("FactoMineR")
      library(FactoMineR, quietly = TRUE)
    }
    if(!require(factoextra, quietly = TRUE)) {
      install.packages("factoextra")
      library(factoextra, quietly = TRUE)
    }
    if(!require(xlsx, quietly = TRUE)) {
      install.packages("xlsx")
      require(xlsx, quietly = TRUE)
    }
    
    ### select the top genes based on variance
    if(num >= 0 && num <= nrow(normalizedMat)) {
      v <- apply(normalizedMat, 1, var)
      v <- v[order(-v)]
      top_genes <- names(v)[1:num]
    } else {
      top_genes <- rownames(normalizedMat)
    }
    
    ### PCA
    pca_result <- PCA(t(normalizedMat[top_genes,]), graph = FALSE)
    colnames(pca_result$ind$coord) <- paste0("PC", 1:ncol(pca_result$ind$coord))
    colnames(pca_result$var$contrib) <- paste0("PC", 1:ncol(pca_result$var$contrib))
    if(is.null(grp2)) {
      pca_group <- data.frame(pca_result$ind$coord, group=grp)
    } else {
      pca_group <- data.frame(pca_result$ind$coord, group=grp, group2=grp2)
    }
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save as png
    if(component[1] == "PC1&PC2") {
      if(is.null(grp2)) {
        ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
          labs(title=paste0(title, "_PC1-2")) +
          geom_point() +
          # geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 16)
      } else {
        ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
          labs(title=paste0(title, "_PC1-2")) +
          geom_point(aes(shape=group2)) +
          # geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 16)  
      }
      ggsave(filename = paste0(outDir, title, "_PC1-2", ".png"), width = 10, height = 8)
      
      if(suppliment) {
        fviz_contrib(pca_result, choice = "var", axes = 1, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC1_contribution.png"), width = 12, height = 8)
        
        fviz_contrib(pca_result, choice = "var", axes = 2, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC2_contribution.png"), width = 12, height = 8)
      }
    } else if(component[1] == "PC2&PC3") {
      if(is.null(grp2)) {
        ggplot(pca_group,aes(x=PC2,y=PC3,col=group)) +
          labs(title=paste0(title, "_PC2-3")) +
          geom_point() +
          # geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 16)
      } else {
        ggplot(pca_group,aes(x=PC2,y=PC3,col=group)) +
          labs(title=paste0(title, "_PC2-3")) +
          geom_point(aes(shape=group2)) +
          # geom_text(aes(label=colnames(normalizedMat)),hjust="inward", vjust="inward") +
          scale_color_manual(values = colors) +
          theme_classic(base_size = 16)  
      }
      ggsave(filename = paste0(outDir, title, "_PC2-3", ".png"), width = 10, height = 8)
      
      if(suppliment) {
        fviz_contrib(pca_result, choice = "var", axes = 2, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC2_contribution.png"), width = 12, height = 8)
        
        fviz_contrib(pca_result, choice = "var", axes = 3, top = 30)
        ggsave(filename = paste0(outDir, title, "_PC3_contribution.png"), width = 12, height = 8)
      }
    } else {
      stop("\"component\" parameter should be \"PC1&PC2\" or \"PC2&PC3\"")
    }
    
    if(suppliment) {
      write.xlsx2(data.frame(Gene_Symbol=rownames(pca_result$var$contrib), pca_result$var$contrib,
                             stringsAsFactors = FALSE, check.names = FALSE),
                  file = paste0(outDir, title, "_PC_contribution.xlsx"),
                  sheetName = "PCA_contribution", row.names = FALSE)
    }
  }
  
  
  #
  ### Split patients into two groups (weighting for # of CAR+ lineages or cells at GMP),
  ### use one group as training and test on the other, then flip them and do the same
  #
  ### first we need to combine the two Seurat objects from different patients
  Seurat_Obj_C2 <- readRDS(file = "./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj_C2.RDS")
  Seurat_Obj_Total <- merge(x = Seurat_Obj, y = Seurat_Obj_C2,
                            add.cell.ids = c("Set1", "Set2"),
                            project = "SJCAR19_Oct2020")
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj_Total@meta.data <- Seurat_Obj_Total@meta.data[colnames(Seurat_Obj_Total@assays$RNA@counts),]
  print(identical(rownames(Seurat_Obj_Total@meta.data), colnames(Seurat_Obj_Total@assays$RNA@counts)))
  
  ### active assay = "RNA"
  Seurat_Obj_Total@active.assay <- "RNA"
  
  ### check whether the orders are the same
  print(identical(names(Idents(object = Seurat_Obj_Total)), rownames(Seurat_Obj_Total@meta.data)))
  
  ### set batch info to the meta data
  temp <- sapply(rownames(Seurat_Obj_Total@meta.data), function(x) substr(x, 1, 4))
  Seurat_Obj_Total@meta.data$Batch <- NA
  Seurat_Obj_Total@meta.data$Batch[which(temp == "Set1")] <- "Batch1"
  Seurat_Obj_Total@meta.data$Batch[which(temp == "Set2")] <- "Batch2"
  
  ### remove the source seurat objects
  rm(Seurat_Obj)
  rm(Seurat_Obj_C2)
  gc()
  
  ### save the new total seurat object
  saveRDS(Seurat_Obj_Total, file = "./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj_Total.RDS")
  
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj_Total@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj_Total@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  all_gmp_last <- intersect(all_gmp_last,
                            which(Seurat_Obj_Total@meta.data$CD4_CD8_by_Consensus == "CD8"))
  all_gmp_not_last <- intersect(all_gmp_not_last,
                                which(Seurat_Obj_Total@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### print out the numbers for each patient
  for(px in unique(Seurat_Obj_Total@meta.data$px)) {
    writeLines(paste(px, length(intersect(all_gmp_last,
                                          which(Seurat_Obj_Total@meta.data$px == px)))))
  }
  
  ### create outputDir2
  outputDir2 <- paste0(outputDir, "DE_Classifier_Two_Group/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### all-in-one function for the two group classifier
  two_group_classifier <- function(Given_Seurat_Obj,
                                   Given_gmp_last_idx,
                                   Given_gmp_not_last_idx,
                                   Group_Info,
                                   seed.k=1234,
                                   featureSelectionNum=100,
                                   second_condition_multiplier=1,
                                   methodTypes = c("svmLinear", "svmRadial", "gbm", "parRF", "glmboost", "knn"),
                                   methodNames = NULL,
                                   file_nums = c(1,2),
                                   output_dir) {
    
    ### set parameters for the classifier
    set.seed(seed.k)
    target_px <- unique(intersect(Given_Seurat_Obj@meta.data$px[Given_gmp_last_idx],
                                  Given_Seurat_Obj@meta.data$px[Given_gmp_not_last_idx]))
    if(is.null(methodNames)) {
      methodNames <- methodTypes
    }
    log_trans_add <- 1
    
    ### performance evaluation
    eval_acc <- vector("list", length(methodTypes))
    names(eval_acc) <- methodNames
    eval_auc <- vector("list", length(methodTypes))
    names(eval_auc) <- methodNames
    
    ### prepare training & test samples
    samps1 <- rownames(Given_Seurat_Obj@meta.data)[intersect(which(Group_Info == "G1"),
                                                             Given_gmp_last_idx)]
    not_last_pool <- intersect(which(Group_Info == "G1"),
                               Given_gmp_not_last_idx)
    pxs <- unique(Given_Seurat_Obj@meta.data$px[which(Group_Info == "G1")])
    samps1_last_num <- length(samps1) * second_condition_multiplier
    sampleNum_per_px <- round(samps1_last_num / length(pxs))
    for(px in pxs) {
      target_pool <- intersect(not_last_pool,
                               which(Given_Seurat_Obj@meta.data$px == px))
      if(length(target_pool) > sampleNum_per_px) {
        samps1 <- c(samps1, rownames(Given_Seurat_Obj@meta.data)[sample(target_pool,
                                                                        sampleNum_per_px)])
      } else {
        samps1 <- c(samps1, rownames(Given_Seurat_Obj@meta.data)[target_pool])
      }
    }
    set_num <- samps1_last_num*(second_condition_multiplier+1)/second_condition_multiplier
    if(samps1_last_num < length(not_last_pool) && length(samps1) < set_num) {
      samps1 <- c(samps1, sample(setdiff(rownames(Given_Seurat_Obj@meta.data)[not_last_pool], samps1),
                                 (set_num-length(samps1))))
    }
    
    samps2 <- rownames(Given_Seurat_Obj@meta.data)[intersect(which(Group_Info == "G2"),
                                                             Given_gmp_last_idx)]
    not_last_pool <- intersect(which(Group_Info == "G2"),
                               Given_gmp_not_last_idx)
    pxs <- unique(Given_Seurat_Obj@meta.data$px[which(Group_Info == "G2")])
    samps2_last_num <- length(samps2) * second_condition_multiplier
    sampleNum_per_px <- round(samps2_last_num / length(pxs))
    for(px in pxs) {
      target_pool <- intersect(not_last_pool,
                               which(Given_Seurat_Obj@meta.data$px == px))
      if(length(target_pool) > sampleNum_per_px) {
        samps2 <- c(samps2, rownames(Given_Seurat_Obj@meta.data)[sample(target_pool,
                                                                        sampleNum_per_px)])
      } else {
        samps2 <- c(samps2, rownames(Given_Seurat_Obj@meta.data)[target_pool])
      }
    }
    set_num <- samps2_last_num*(second_condition_multiplier+1)/second_condition_multiplier
    if(samps2_last_num < length(not_last_pool) && length(samps2) < set_num) {
      samps2 <- c(samps2, sample(setdiff(rownames(Given_Seurat_Obj@meta.data)[not_last_pool], samps2),
                                 (set_num-length(samps2))))
    }
    
    #
    ### perform classifier a) training: samps1 test: samps2, b) training: samps2 test: samps1
    #
    ### new obj for the training data & set idents with the info
    classifier_seurat_obj <- subset(Given_Seurat_Obj, cells = samps1)
    classifier_seurat_obj <- SetIdent(object = classifier_seurat_obj,
                                      cells = rownames(classifier_seurat_obj@meta.data),
                                      value = classifier_seurat_obj@meta.data$GMP_CARpos_Persister)
    
    ### DE analysis
    de_result <- FindMarkers(classifier_seurat_obj,
                             ident.1 = "YES",
                             ident.2 = "NO",
                             min.pct = 0.1,
                             logfc.threshold = 0.1,
                             test.use = "wilcox")
    
    ### write out the DE result
    write.xlsx2(data.frame(Gene=rownames(de_result),
                           de_result,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(output_dir, "/GMP_CD8_CARpos_Persisters_vs_NonPersisters_", file_nums[1], ".xlsx"),
                sheetName = "GMP_CD8_CARpos_DE_Result", row.names = FALSE)
    
    ### normalize the read counts
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE),
                                         filter_thresh = 0)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(classifier_seurat_obj@meta.data[rownames(input_data),"GMP_CARpos_Persister"],
                               levels = c("YES", "NO"))
    
    ### new obj for the test data & set idents with the info
    classifier_seurat_obj <- subset(Given_Seurat_Obj, cells = samps2)
    classifier_seurat_obj <- SetIdent(object = classifier_seurat_obj,
                                      cells = rownames(classifier_seurat_obj@meta.data),
                                      value = classifier_seurat_obj@meta.data$GMP_CARpos_Persister)
    
    ### DE analysis for the test data
    de_result2 <- FindMarkers(classifier_seurat_obj,
                              ident.1 = "YES",
                              ident.2 = "NO",
                              min.pct = 0.1,
                              logfc.threshold = 0.1,
                              test.use = "wilcox")
    
    ### write out the DE result
    write.xlsx2(data.frame(Gene=rownames(de_result2),
                           de_result2,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(output_dir, "/GMP_CD8_CARpos_Persisters_vs_NonPersisters_", file_nums[2], ".xlsx"),
                sheetName = "GMP_CD8_CARpos_DE_Result", row.names = FALSE)
    
    ### normalize the read counts
    test_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],] + log_trans_add,
                                                               stringsAsFactors = FALSE, check.names = FALSE),
                                        filter_thresh = 0)
    
    ### annotate class for the test data
    test_data <- data.frame(t(test_data), stringsAsFactors = FALSE, check.names = FALSE)
    test_data$Class <- factor(classifier_seurat_obj@meta.data[rownames(test_data),"GMP_CARpos_Persister"],
                              levels = c("YES", "NO"))
    
    ### train control options
    train_control <- trainControl(method="none", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
    
    # 1
    ### build classifier and test
    p <- list()
    acc <- NA
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=input_data, method=methodTypes[j], trControl = train_control)
      pred_result <- predict(model, newdata = test_data)
      acc <- round((sum(pred_result == test_data$Class) / nrow(test_data)), 3)
      pred_result <- predict(model, newdata = test_data, type = "prob")
      roc <- roc(test_data$Class, pred_result$YES)
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using DE Genes\n",
                                           "Accuracy =", acc),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      eval_acc[[methodNames[j]]] <- c(eval_acc[[methodNames[j]]], acc)
      eval_auc[[methodNames[j]]] <- c(eval_auc[[methodNames[j]]], as.numeric(roc$auc))
      gc()
    }
    
    ### draw ROC curves
    png(paste0(output_dir, "/Classifier_Using_DEG_GMP_Last_vs_Not_Last_", featureSelectionNum, "_Two_Group_(", file_nums[1], ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", eval_acc[[methodNames[j]]][1]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
    
    ### PCA
    pca_plot(normalizedMat = data.frame(t(test_data[,-ncol(test_data)]),
                                        stringsAsFactors = FALSE, check.names = FALSE),
             grp = as.character(test_data$Class),
             title = paste0("PCA_Classifier_", featureSelectionNum, "_Two_Group_(", file_nums[1], ")"),
             outDir = output_dir)
    
    # 2
    ### new obj for the training data
    classifier_seurat_obj <- subset(Given_Seurat_Obj, cells = samps1)
    
    ### normalize the read counts
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result2)[1:featureSelectionNum],] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE),
                                         filter_thresh = 0)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(classifier_seurat_obj@meta.data[rownames(input_data),"GMP_CARpos_Persister"],
                               levels = c("YES", "NO"))
    
    ### new obj for the test data
    classifier_seurat_obj <- subset(Given_Seurat_Obj, cells = samps2)
    
    ### normalize the read counts
    test_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result2)[1:featureSelectionNum],] + log_trans_add,
                                                               stringsAsFactors = FALSE, check.names = FALSE),
                                        filter_thresh = 0)
    
    ### annotate class for the test data
    test_data <- data.frame(t(test_data), stringsAsFactors = FALSE, check.names = FALSE)
    test_data$Class <- factor(classifier_seurat_obj@meta.data[rownames(test_data),"GMP_CARpos_Persister"],
                              levels = c("YES", "NO"))
    
    ### build classifier and test
    p <- list()
    acc <- NA
    for(j in 1:length(methodTypes)) {
      writeLines(paste(methodTypes[j]))
      model <- train(Class~., data=test_data, method=methodTypes[j], trControl = train_control)
      pred_result <- predict(model, newdata = input_data)
      acc <- round((sum(pred_result == input_data$Class) / nrow(input_data)), 3)
      pred_result <- predict(model, newdata = input_data, type = "prob")
      roc <- roc(input_data$Class, pred_result$YES)
      p[[j]] <- plot.roc(roc, main = paste(methodNames[j], "Using DE Genes\n",
                                           "Accuracy =", acc),
                         legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                         xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
      eval_acc[[methodNames[j]]] <- c(eval_acc[[methodNames[j]]], acc)
      eval_auc[[methodNames[j]]] <- c(eval_auc[[methodNames[j]]], as.numeric(roc$auc))
      gc()
    }
    
    ### draw ROC curves
    png(paste0(output_dir, "/Classifier_Using_DEG_GMP_Last_vs_Not_Last_", featureSelectionNum, "_Two_Group_(", file_nums[2], ").png"),
        width = 2000, height = 2000, res = 350)
    par(mfrow=c(3, 2))
    for(j in 1:length(methodTypes)) {
      plot.roc(p[[j]], main = paste(methodNames[j], "Using Gene Expressions\n",
                                    "Accuracy =", eval_acc[[methodNames[j]]][2]),
               legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
               xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    }
    dev.off()
    
    ### PCA
    pca_plot(normalizedMat = data.frame(t(input_data[,-ncol(input_data)]),
                                        stringsAsFactors = FALSE, check.names = FALSE),
             grp = as.character(input_data$Class),
             title = paste0("PCA_Classifier_", featureSelectionNum, "_Two_Group_(", file_nums[2], ")"),
             outDir = output_dir)
    
    ### print out the sample ratio
    writeLines(paste("#_Set1_GMP_Last:", length(which(input_data$Class == "YES")),
                     "#_Set1:GMP_NOT_Last:", length(which(input_data$Class == "NO"))))
    writeLines(paste("#_Set2_GMP_Last:", length(which(test_data$Class == "YES")),
                     "#_Set2:GMP_NOT_Last:", length(which(test_data$Class == "NO"))))
    
  }
  
  ### set classification methods
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "parRF", "glmboost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "Linear_Model", "KNN")
  
  ### divide the patients into two groups
  #
  ### GMP_Persister_CD8_Cell_#
  ### 8941 (11) + 1347 (3) + 118 (12) + 162 (5) + 72 (9) + 70 (13) = 10710
  ### 4013 (6) + 2189 (2) + 1954 (4) + 496 (10) + 282 (7) + 256 (8) = 9190
  #
  ### GMP_Persister_CARpos_Lineage_#
  ### 108 (11) + 0 (3) + 0 (12) + 18 (5) + 4 (9) + 1 (13) = 131
  ### 20 (6) + 4 (2) + 6 (4) + 34 (10) + 9 (7) + 25 (8) = 98
  
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-03",
                                                                     "SJCAR19-05",
                                                                     "SJCAR19-09",
                                                                     "SJCAR19-11",
                                                                     "SJCAR19-12",
                                                                     "SJCAR19-13"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-02",
                                                                     "SJCAR19-04",
                                                                     "SJCAR19-06",
                                                                     "SJCAR19-07",
                                                                     "SJCAR19-08",
                                                                     "SJCAR19-10"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(1,2),
                       output_dir = outputDir2)
  
  
  ### divide the patients into two groups
  #
  ### GMP Persister CD8 Cell Num
  ### 966 (6) + 89 (10) + 63 (4) + 19 (12) + 16 (2) + 8 (13)
  ### 523 (7) + 254 (11) + 118 (8) + 57 (5) + 32 (3) + 8 (9)
  
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06",
                                                                     "SJCAR19-10",
                                                                     "SJCAR19-04",
                                                                     "SJCAR19-12",
                                                                     "SJCAR19-02",
                                                                     "SJCAR19-13"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-07",
                                                                     "SJCAR19-11",
                                                                     "SJCAR19-08",
                                                                     "SJCAR19-05",
                                                                     "SJCAR19-03",
                                                                     "SJCAR19-09"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(3,4),
                       output_dir = outputDir2)
  
  
  ### divide the patients into two groups
  #
  ### GMP Persister CD8 Cell Num
  ### 966 (6) + 118 (8) + 57 (5) + 19 (12) + 16 (2) + 8 (9)
  ### 523 (7) + 254 (11) + 89 (10) + 63 (4) + 32 (3) + 8 (13)
  
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06",
                                                                     "SJCAR19-08",
                                                                     "SJCAR19-05",
                                                                     "SJCAR19-12",
                                                                     "SJCAR19-02",
                                                                     "SJCAR19-09"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-07",
                                                                     "SJCAR19-11",
                                                                     "SJCAR19-10",
                                                                     "SJCAR19-04",
                                                                     "SJCAR19-03",
                                                                     "SJCAR19-13"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(5,6),
                       output_dir = outputDir2)
  
  
  ### divide the patients into two groups
  #
  ### GMP Persister CD8 Cell Num
  ### 966 (6) + 254 (11) + 32 (3) + 19 (12) + 16 (2) + 8 (9)
  ### 523 (7) + 118 (8) + 89 (10) + 63 (4) + 57 (5) + 8 (13)
  
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06",
                                                                     "SJCAR19-11",
                                                                     "SJCAR19-03",
                                                                     "SJCAR19-12",
                                                                     "SJCAR19-02",
                                                                     "SJCAR19-09"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-07",
                                                                     "SJCAR19-08",
                                                                     "SJCAR19-10",
                                                                     "SJCAR19-04",
                                                                     "SJCAR19-05",
                                                                     "SJCAR19-13"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(7,8),
                       output_dir = outputDir2)
  
  
  ### divide the patients into two groups
  #
  ### GMP Persister CD8 Cell Num
  ### 966 (6) + 254 (11) + 32 (3) + 16 (2) + 8 (9)
  ### 523 (7) + 118 (8) + 89 (10) + 63 (4) + 57 (5)
  
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06",
                                                                     "SJCAR19-11",
                                                                     "SJCAR19-03",
                                                                     "SJCAR19-02",
                                                                     "SJCAR19-09"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-07",
                                                                     "SJCAR19-08",
                                                                     "SJCAR19-10",
                                                                     "SJCAR19-04",
                                                                     "SJCAR19-05"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(9,10),
                       output_dir = outputDir2)
  
  
  ### divide the patients into two groups
  #
  ### GMP Persister CD8 Cell Num
  ### 966 (6) + 89 (10) + 32 (3) + 16 (2) + 8 (9)
  ### 523 (7) + 118 (8) + 254 (11) + 63 (4) + 57 (5)
  
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06",
                                                                     "SJCAR19-10",
                                                                     "SJCAR19-03",
                                                                     "SJCAR19-02",
                                                                     "SJCAR19-09"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-07",
                                                                     "SJCAR19-08",
                                                                     "SJCAR19-11",
                                                                     "SJCAR19-04",
                                                                     "SJCAR19-05"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(11,12),
                       output_dir = outputDir2)
  
  
  ### we want to perform the classifier with random cells
  ### but this time, regardless of patient info,
  ### and use the same number for training & testing
  #
  ### divide the patients into two groups randomly
  set.seed(1234)
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[sample(length(Seurat_Obj_Total$Classifier_Group), round(length(Seurat_Obj_Total$Classifier_Group)/2))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(is.na(Seurat_Obj_Total$Classifier_Group))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(13,14),
                       output_dir = outputDir2)
  
  #
  ### PCA and UMAP comparing cells in a persister lineage VS not but with ALL the genes
  #
  
  ### set idents with the library
  Seurat_Obj_Total <- SetIdent(object = Seurat_Obj_Total,
                               cells = rownames(Seurat_Obj_Total@meta.data),
                               value = Seurat_Obj_Total@meta.data$GMP_CARpos_Persister)
  
  ### extract GMP cells only
  Seurat_Obj_GMP <- subset(Seurat_Obj_Total, idents = c("YES", "NO"))
  
  ### normalization
  Seurat_Obj_GMP <- NormalizeData(Seurat_Obj_GMP,
                                  normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  Seurat_Obj_GMP <- FindVariableFeatures(Seurat_Obj_GMP,
                                         selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  Seurat_Obj_GMP <- ScaleData(Seurat_Obj_GMP,
                              vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### PCA
  Seurat_Obj_GMP <- RunPCA(Seurat_Obj_GMP,
                           features = VariableFeatures(object = Seurat_Obj_GMP),
                           npcs = 15)
  
  ### UMAP
  Seurat_Obj_GMP <- RunUMAP(Seurat_Obj_GMP, dims = 1:15)
  
  ### UMAP with GMP
  p <- DimPlot(object = Seurat_Obj_GMP, reduction = "umap",
               group.by = "GMP_CARpos_Persister", split.by = NULL,
               pt.size = 2) +
       ggtitle("UMAP of SJCAR19 Data") +
       theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
       labs(color="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir, "/", "UMAP_Plot_GMP1.png"), plot = p, width = 20, height = 12, dpi = 300)
  
  ### UMAP with GMP by each patient
  p <- DimPlot(object = Seurat_Obj_GMP, reduction = "umap",
               group.by = "GMP_CARpos_Persister", split.by = "px",
               pt.size = 2, ncol = 3) +
       ggtitle("UMAP of SJCAR19 Data") +
       theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
       labs(color="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir, "/", "UMAP_Plot_GMP2.png"), plot = p, width = 20, height = 12, dpi = 300)
  
  
  ### divide the patients into two groups BUT REMOVE NON-RESPONDERS THIS TIME
  #
  ### GMP Persister CD8 Cell Num
  ### 966 (6) + 63 (4) + 19 (12) + 16 (2)
  ### 254 (11) + 118 (8) + 89 (10) + 57 (5) + 32 (3) + 8 (13)
  
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06",
                                                                     "SJCAR19-04",
                                                                     "SJCAR19-12",
                                                                     "SJCAR19-02"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-11",
                                                                     "SJCAR19-08",
                                                                     "SJCAR19-10",
                                                                     "SJCAR19-05",
                                                                     "SJCAR19-03",
                                                                     "SJCAR19-13"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(15,16),
                       output_dir = outputDir2)
  
  ### divide the patients into two groups BUT REMOVE NON-RESPONDERS THIS TIME
  #
  ### GMP Persister CD8 Cell Num
  ### 966 (6) + 63 (4) + 19 (12) + 16 (2)
  ### 254 (11) + 118 (8) + 89 (10) + 57 (5) + 32 (3) + 8 (13)
  
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06",
                                                                     "SJCAR19-04",
                                                                     "SJCAR19-12",
                                                                     "SJCAR19-02"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-11",
                                                                     "SJCAR19-08",
                                                                     "SJCAR19-10",
                                                                     "SJCAR19-05",
                                                                     "SJCAR19-03",
                                                                     "SJCAR19-13"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 4321,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(17,18),
                       output_dir = outputDir2)
  
  
  ### set classification methods
  methodTypes <- c("svmLinear", "svmRadial", "svmLinearWeights", "parRF", "glmboost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "SVMLinearWeights", "RandomForest", "Linear_Model", "KNN")
  
  ### If you take each patient’s data and split it 50/50 test and training, does the accuracy/AUC look as bad
  set.seed(1234)
  Seurat_Obj_Total <- SetIdent(object = Seurat_Obj_Total,
                         cells = rownames(Seurat_Obj_Total@meta.data),
                         value = Seurat_Obj_Total@meta.data$px)
  outputDir3 <- paste0(outputDir2, "50_50_by_each_patient/")
  dir.create(outputDir3, showWarnings = FALSE, recursive = TRUE)
  for(px in unique(Seurat_Obj_Total@meta.data$px[which(Seurat_Obj_Total@meta.data$GMP_CARpos_Persister == "YES")])) {
    ### get patient's seurat object
    Seurat_Obj_px <- subset(Seurat_Obj_Total, idents = px)
    
    ### the indicies of the persisters
    px_gmp_last <- which(Seurat_Obj_px@meta.data$GMP_CARpos_Persister == "YES")
    px_gmp_not_last <- which(Seurat_Obj_px@meta.data$GMP_CARpos_Persister == "NO")
    
    ### only use the CD8 cells
    px_gmp_last <- intersect(px_gmp_last,
                             which(Seurat_Obj_px@meta.data$CD4_CD8_by_Consensus == "CD8"))
    px_gmp_not_last <- intersect(px_gmp_not_last,
                                 which(Seurat_Obj_px@meta.data$CD4_CD8_by_Consensus == "CD8"))
    
    #
    ### DE genes
    #
    ### set idents
    Seurat_Obj_px@meta.data$GMP_CD8_CARpos_Persister <- NA
    Seurat_Obj_px@meta.data$GMP_CD8_CARpos_Persister[px_gmp_last] <- "YES"
    Seurat_Obj_px@meta.data$GMP_CD8_CARpos_Persister[px_gmp_not_last] <- "NO"
    Seurat_Obj_px <- SetIdent(object = Seurat_Obj_px,
                              cells = rownames(Seurat_Obj_px@meta.data),
                              value = Seurat_Obj_px@meta.data$GMP_CD8_CARpos_Persister)
    
    ### DE analysis
    de_result <- FindMarkers(Seurat_Obj_px,
                             ident.1 = "YES",
                             ident.2 = "NO",
                             min.pct = 0.1,
                             logfc.threshold = 0.2,
                             test.use = "wilcox")
    
    ### write out the DE result
    write.xlsx2(data.frame(Gene=rownames(de_result),
                           de_result,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir3, px, "_GMP_CD8_CARpos_Persisters_vs_NonPersisters.xlsx"),
                sheetName = "GMP_CD8_CARpos_DE_Result", row.names = FALSE)
    
    ### calculate the number of one condition in one group
    half_gmp_last_num <- floor(length(px_gmp_last) / 2)
    
    ### set groups
    Seurat_Obj_px$Classifier_Group <- NA
    Seurat_Obj_px$Classifier_Group[sample(px_gmp_last, half_gmp_last_num)] <- "G1"
    Seurat_Obj_px$Classifier_Group[sample(setdiff(px_gmp_last,
                                                  which(Seurat_Obj_px$Classifier_Group == "G1")), half_gmp_last_num)] <- "G2"
    Seurat_Obj_px$Classifier_Group[sample(px_gmp_not_last, half_gmp_last_num)] <- "G1"
    Seurat_Obj_px$Classifier_Group[sample(setdiff(px_gmp_not_last,
                                                  which(Seurat_Obj_px$Classifier_Group == "G1")), half_gmp_last_num)] <- "G2"
    
    ### perform classification
    two_group_classifier(Given_Seurat_Obj = Seurat_Obj_px,
                         Given_gmp_last_idx = px_gmp_last,
                         Given_gmp_not_last_idx = px_gmp_not_last,
                         Group_Info = Seurat_Obj_px$Classifier_Group,
                         seed.k = 1234,
                         featureSelectionNum = 100,
                         methodTypes = methodTypes,
                         methodNames = methodNames,
                         file_nums = c(paste0(px, "-1"), paste0(px, "-2")),
                         output_dir = outputDir3)
  }
  
  
  #
  ### Kaity's request - dot plot
  #
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  all_gmp_last <- intersect(all_gmp_last,
                            which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  all_gmp_not_last <- intersect(all_gmp_not_last,
                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### GMP CD8 CARpos persisters
  Seurat_Obj@meta.data$GMP_CD8_CARpos_Persister <- NA
  Seurat_Obj@meta.data$GMP_CD8_CARpos_Persister[all_gmp_last] <- "YES"
  Seurat_Obj@meta.data$GMP_CD8_CARpos_Persister[all_gmp_not_last] <- "NO"
  
  ### set idents
  Seurat_Obj <- SetIdent(object = Seurat_Obj,
                         cells = rownames(Seurat_Obj@meta.data),
                         value = Seurat_Obj@meta.data$GMP_CD8_CARpos_Persister)
  
  ### DE analysis
  de_result <- FindMarkers(Seurat_Obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/GMP_CD8_CARpos_Persisters_vs_NonPersisters.xlsx"),
              sheetName = "GMP_CD8_CARpos_DE_Result", row.names = FALSE)
  
  ### pathway analysis
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                            rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                            "ENTREZID", "SYMBOL"),
                                          org = "human", database = "GO",
                                          title = paste0("Pathway_Result_GMP_CD8_CARpos_Persisters_vs_NonPersisters"),
                                          displayNum = 30, imgPrint = TRUE,
                                          dir = paste0(outputDir))
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathway_Result_GMP_CD8_CARpos_Persisters_vs_NonPersisters"),
                                            displayNum = 30, imgPrint = TRUE,
                                            dir = paste0(outputDir))
  if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
    write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_Pathway_Result_GMP_CD8_CARpos_Persisters_vs_NonPersisters.xlsx"),
                row.names = FALSE, sheetName = paste0("GO_Result"))
  }
  if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
    write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_Pathway_Result_GMP_CD8_CARpos_Persisters_vs_NonPersisters.xlsx"),
                row.names = FALSE, sheetName = paste0("KEGG_Result"))
  }
  
  ### beeswarm to compare gene expressions between p & np
  set.seed(1234)
  iter <- 10
  for(i in 1:iter) {
    target_data <- data.frame(Seurat_Obj@assays$RNA@counts[rownames(de_result)[1:25],
                                                           c(sample(which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES"), 50),
                                                             sample(which(Seurat_Obj@meta.data$GMP_CARpos_Persister== "NO"), 50))],
                              stringsAsFactors = FALSE, check.names = FALSE)
    plot_df <- data.frame(matrix(NA, 2500, 4),
                          stringsAsFactors = FALSE, check.names = FALSE)
    colnames(plot_df) <- c("Gene_Name", "Cell_Name", "GEX", "Persistency")
    for(j in 1:nrow(target_data)) {
      for(k in 1:ncol(target_data)) {
        plot_df[100*(j-1)+k,"Gene_Name"] <- rownames(target_data)[j]
        plot_df[100*(j-1)+k,"Cell_Name"] <- colnames(target_data)[k]
        plot_df[100*(j-1)+k,"GEX"] <- as.numeric(target_data[j,k])
      }
      plot_df[(100*(j-1)+1):(100*j),"Persistency"] <- c(rep("YES", 50), rep("NO", 50))
    }
    
    p <- vector("list", length = 25)
    names(p) <- rownames(de_result)[1:25]
    for(j in 1:length(p)) {
      temp_plot_df <- plot_df[which(plot_df$Gene_Name == names(p)[j]),]
      p[[j]] <- ggplot(temp_plot_df, aes_string(x="Persistency", y="GEX")) +
        geom_boxplot() +
        geom_beeswarm(aes_string(col="Persistency"), na.rm = TRUE) +
        stat_compare_means() +
        xlab("") + ylab("Gene Expression") +
        labs(col="Is_Persistent") +
        ggtitle(paste0(names(p)[j])) +
        theme_classic(base_size = 16)
    }
    
    g <- arrangeGrob(grobs = p,
                     nrow = 5,
                     ncol = 5,
                     top = paste0("GMP_CD8_CARpos_Persisters_vs_NonPersisters"))
    ggsave(file = paste0(outputDir, "/GMP_CD8_CARpos_Persisters_vs_NonPersisters(", i, ").png"), g,
           width = 30, height = 20, dpi = 300)
  }
  
  ### only get the cells of interests
  target_Seurat_Obj <- subset(Seurat_Obj, idents = c("YES", "NO"))
  
  ### Dotplot2
  p <- NULL
  p[[1]] <- DotPlot(target_Seurat_Obj,
          features = rownames(de_result)[30:1],
          cols = c("skyblue", "pink"),
          group.by = "GMP_CD8_CARpos_Persister") +
    coord_flip() +
    # ggtitle("Persisters vs Non-Persisters Avg Gene Exp") +
    xlab("Top 1:30 DE Genes") +
    ylab("Is_Persistent") +
    theme_calc() +
    theme(plot.title = element_text(hjust = 0.5))
  p[[2]] <- DotPlot(target_Seurat_Obj,
                    features = rownames(de_result)[60:31],
                    cols = c("skyblue", "pink"),
                    group.by = "GMP_CD8_CARpos_Persister") +
    coord_flip() +
    # ggtitle("Persisters vs Non-Persisters Avg Gene Exp") +
    xlab("Top 31:60 DE Genes") +
    ylab("Is_Persistent") +
    theme_calc() +
    theme(plot.title = element_text(hjust = 0.5))
  g <- arrangeGrob(grobs = p,
                   nrow = 1,
                   ncol = 2,
                   top = paste0("GMP CD8 CARpos Persisters vs Non-Persisters Avg Gene Exp"))
  ggsave(file = paste0(outputDir, "GMP_CD8_CARpos_Persisters_vs_NonPersisters.png"),
         g, width = 20, height = 12, dpi = 300)
  
  ### Genes of interest selected by Kaity
  genes_of_interest <- c("CD52", "SELL", "GNLY", "GZMK", "CCL5", "ID2", "GZMB", "TNFRSF18")
  DotPlot(target_Seurat_Obj,
          features = genes_of_interest,
          cols = c("skyblue", "pink"),
          group.by = "GMP_CD8_CARpos_Persister") +
    coord_flip() +
    # ggtitle("Persisters vs Non-Persisters Avg Gene Exp") +
    xlab("Genes of Interest") +
    ylab("Is_Persistent") +
    theme_calc(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir, "GMP_CD8_CARpos_Persisters_vs_NonPersisters2.png"),
         width = 12, height = 12, dpi = 300)
  
  
  #
  ### build a classifier with Px6 and use other patients' samples for testing
  #
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-02",
                                                                     "SJCAR19-03",
                                                                     "SJCAR19-04",
                                                                     "SJCAR19-05",
                                                                     "SJCAR19-08",
                                                                     "SJCAR19-10",
                                                                     "SJCAR19-11",
                                                                     "SJCAR19-12",
                                                                     "SJCAR19-13"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(19,20),
                       output_dir = outputDir2)
  
  #
  ### build a classifier with Px6 and use other patients' samples for testing
  #
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-11"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(21,22),
                       output_dir = outputDir2)
  
  #
  ### build a classifier with Px6 and use other patients' samples for testing
  #
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-08"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(23,24),
                       output_dir = outputDir2)
  
  #
  ### build a classifier with Px6 and use other patients' samples for testing
  #
  Seurat_Obj_Total$Classifier_Group <- NA
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-06"))] <- "G1"
  Seurat_Obj_Total$Classifier_Group[which(Seurat_Obj_Total$px %in% c("SJCAR19-10"))] <- "G2"
  
  ### perform classification
  two_group_classifier(Given_Seurat_Obj = Seurat_Obj_Total,
                       Given_gmp_last_idx = all_gmp_last,
                       Given_gmp_not_last_idx = all_gmp_not_last,
                       Group_Info = Seurat_Obj_Total$Classifier_Group,
                       seed.k = 1234,
                       featureSelectionNum = 100,
                       methodTypes = methodTypes,
                       methodNames = methodNames,
                       file_nums = c(25,26),
                       output_dir = outputDir2)
  
  #
  ### What if you double the number of non-persisters sampled?
  #
  
  ### set classification methods
  methodTypes <- c("svmLinear", "svmRadial", "svmLinearWeights", "parRF", "glmboost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "SVMLinearWeights", "RandomForest", "Linear_Model", "KNN")
  
  ### If you take each patient’s data and split it 50/50 test and training, does the accuracy/AUC look as bad
  set.seed(1234)
  Seurat_Obj_Total <- SetIdent(object = Seurat_Obj_Total,
                               cells = rownames(Seurat_Obj_Total@meta.data),
                               value = Seurat_Obj_Total@meta.data$px)
  outputDir3 <- paste0(outputDir2, "50_100_by_each_patient/")
  dir.create(outputDir3, showWarnings = FALSE, recursive = TRUE)
  for(px in unique(Seurat_Obj_Total@meta.data$px[which(Seurat_Obj_Total@meta.data$GMP_CARpos_Persister == "YES")])) {
    ### get patient's seurat object
    Seurat_Obj_px <- subset(Seurat_Obj_Total, idents = px)
    
    ### the indicies of the persisters
    px_gmp_last <- which(Seurat_Obj_px@meta.data$GMP_CARpos_Persister == "YES")
    px_gmp_not_last <- which(Seurat_Obj_px@meta.data$GMP_CARpos_Persister == "NO")
    
    ### only use the CD8 cells
    px_gmp_last <- intersect(px_gmp_last,
                             which(Seurat_Obj_px@meta.data$CD4_CD8_by_Consensus == "CD8"))
    px_gmp_not_last <- intersect(px_gmp_not_last,
                                 which(Seurat_Obj_px@meta.data$CD4_CD8_by_Consensus == "CD8"))
    
    #
    ### DE genes
    #
    ### set idents
    Seurat_Obj_px@meta.data$GMP_CD8_CARpos_Persister <- NA
    Seurat_Obj_px@meta.data$GMP_CD8_CARpos_Persister[px_gmp_last] <- "YES"
    Seurat_Obj_px@meta.data$GMP_CD8_CARpos_Persister[px_gmp_not_last] <- "NO"
    Seurat_Obj_px <- SetIdent(object = Seurat_Obj_px,
                              cells = rownames(Seurat_Obj_px@meta.data),
                              value = Seurat_Obj_px@meta.data$GMP_CD8_CARpos_Persister)
    
    ### DE analysis
    de_result <- FindMarkers(Seurat_Obj_px,
                             ident.1 = "YES",
                             ident.2 = "NO",
                             min.pct = 0.1,
                             logfc.threshold = 0.2,
                             test.use = "wilcox")
    
    ### write out the DE result
    write.xlsx2(data.frame(Gene=rownames(de_result),
                           de_result,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir3, px, "_GMP_CD8_CARpos_Persisters_vs_NonPersisters.xlsx"),
                sheetName = "GMP_CD8_CARpos_DE_Result", row.names = FALSE)
    
    ### calculate the number of one condition in one group
    half_gmp_last_num <- floor(length(px_gmp_last) / 2)
    
    ### set groups
    Seurat_Obj_px$Classifier_Group <- NA
    Seurat_Obj_px$Classifier_Group[sample(px_gmp_last, half_gmp_last_num)] <- "G1"
    Seurat_Obj_px$Classifier_Group[sample(setdiff(px_gmp_last,
                                                  which(Seurat_Obj_px$Classifier_Group == "G1")), half_gmp_last_num)] <- "G2"
    Seurat_Obj_px$Classifier_Group[sample(px_gmp_not_last, floor(length(px_gmp_not_last)/2))] <- "G1"
    Seurat_Obj_px$Classifier_Group[sample(setdiff(px_gmp_not_last,
                                                  which(Seurat_Obj_px$Classifier_Group == "G1")), floor(length(px_gmp_not_last)/2))] <- "G2"
    
    ### perform classification
    two_group_classifier(Given_Seurat_Obj = Seurat_Obj_px,
                         Given_gmp_last_idx = px_gmp_last,
                         Given_gmp_not_last_idx = px_gmp_not_last,
                         Group_Info = Seurat_Obj_px$Classifier_Group,
                         seed.k = 1234,
                         featureSelectionNum = 100,
                         methodTypes = methodTypes,
                         methodNames = methodNames,
                         second_condition_multiplier = 2,
                         file_nums = c(paste0(px, "-1"), paste0(px, "-2")),
                         output_dir = outputDir3)
  }
  
  #
  ### If you change SJCAR19-08 to k=5 cross validation, does anything change?
  #
  
  ### If you take each patient’s data and split it 50/50 test and training, does the accuracy/AUC look as bad
  Seurat_Obj_Total <- SetIdent(object = Seurat_Obj_Total,
                               cells = rownames(Seurat_Obj_Total@meta.data),
                               value = Seurat_Obj_Total@meta.data$px)
  outputDir3 <- paste0(outputDir2, "5k_CV_Px08/")
  dir.create(outputDir3, showWarnings = FALSE, recursive = TRUE)
  px <- "SJCAR19-08"
  
  ### get patient's seurat object
  Seurat_Obj_px <- subset(Seurat_Obj_Total, idents = px)
  
  ### the indicies of the persisters
  px_gmp_last <- which(Seurat_Obj_px@meta.data$GMP_CARpos_Persister == "YES")
  px_gmp_not_last <- which(Seurat_Obj_px@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  px_gmp_last <- intersect(px_gmp_last,
                           which(Seurat_Obj_px@meta.data$CD4_CD8_by_Consensus == "CD8"))
  px_gmp_not_last <- intersect(px_gmp_not_last,
                               which(Seurat_Obj_px@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### set options
  set.seed(1234)
  featureSelectionNum = 100
  log_trans_add <- 1
  methodTypes <- c("svmLinear", "svmRadial", "svmLinearWeights", "parRF", "glmboost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "SVMLinearWeights", "RandomForest", "Linear_Model", "KNN")
  
  ### set 5k CV cell list
  cv_sample_num <- round(length(px_gmp_last)/5)
  cv_cell_list <- vector("list", 5)
  gmp_last_so_far <- NULL
  gmp_not_last_so_far <- NULL
  for(i in 1:5) {
    ### get random cells
    cv_cell_list[[i]] <- vector("list", 2)
    names(cv_cell_list[[i]]) <- c("GMP_Last", "GMP_Not_Last")
    
    if(i != 5) {
      cv_cell_list[[i]][["GMP_Last"]] <- sample(setdiff(px_gmp_last, gmp_last_so_far), cv_sample_num)
      cv_cell_list[[i]][["GMP_Not_Last"]] <- sample(setdiff(px_gmp_not_last, gmp_not_last_so_far), cv_sample_num)
    } else {
      cv_cell_list[[i]][["GMP_Last"]] <- setdiff(px_gmp_last, gmp_last_so_far)
      cv_cell_list[[i]][["GMP_Not_Last"]] <- sample(setdiff(px_gmp_not_last, gmp_not_last_so_far), length(cv_cell_list[[i]][["GMP_Last"]]))
    }
    gmp_last_so_far <- c(gmp_last_so_far, cv_cell_list[[i]][["GMP_Last"]])
    gmp_not_last_so_far <- c(gmp_not_last_so_far, cv_cell_list[[i]][["GMP_Not_Last"]])
  }
  
  ### perform classification
  all_samps <- rownames(Seurat_Obj_px@meta.data)[c(gmp_last_so_far, gmp_not_last_so_far)]
  eval_acc <- vector("list", 5)
  eval_roc <- vector("list", 5)
  for(i in 1:5) {
    ### set samples for cv
    ts_samps <- rownames(Seurat_Obj_px@meta.data)[c(cv_cell_list[[i]][[1]], cv_cell_list[[i]][[2]])]
    tr_samps <- setdiff(all_samps, ts_samps)
    
    ### new obj for the training data & set idents with the info - TR SAMPS
    classifier_seurat_obj <- subset(Seurat_Obj_px, cells = tr_samps)
    classifier_seurat_obj <- SetIdent(object = classifier_seurat_obj,
                                      cells = rownames(classifier_seurat_obj@meta.data),
                                      value = classifier_seurat_obj@meta.data$GMP_CARpos_Persister)
    
    ### DE analysis
    de_result <- FindMarkers(classifier_seurat_obj,
                             ident.1 = "YES",
                             ident.2 = "NO",
                             min.pct = 0.1,
                             logfc.threshold = 0.1,
                             test.use = "wilcox")
    
    ### normalize the read counts
    input_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],] + log_trans_add,
                                                                stringsAsFactors = FALSE, check.names = FALSE),
                                         filter_thresh = 0)
    
    ### annotate class for the input data
    input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
    input_data$Class <- factor(classifier_seurat_obj@meta.data[rownames(input_data),"GMP_CARpos_Persister"],
                               levels = c("YES", "NO"))
    
    ### new obj for the test data & set idents with the info - TS SAMPS
    classifier_seurat_obj <- subset(Seurat_Obj_px, cells = ts_samps)
    
    ### normalize the read counts
    test_data <- normalizeRNASEQwithVST(readCount = data.frame(classifier_seurat_obj@assays$RNA@counts[rownames(de_result)[1:featureSelectionNum],] + log_trans_add,
                                                               stringsAsFactors = FALSE, check.names = FALSE),
                                        filter_thresh = 0)
    
    ### annotate class for the test data
    test_data <- data.frame(t(test_data), stringsAsFactors = FALSE, check.names = FALSE)
    test_data$Class <- factor(classifier_seurat_obj@meta.data[rownames(test_data),"GMP_CARpos_Persister"],
                              levels = c("YES", "NO"))
    
    ### train control options
    train_control <- trainControl(method="none", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
    
    ### lists for acc & auc
    eval_acc[[i]] <- vector("list", length(methodTypes))
    eval_roc[[i]] <- vector("list", length(methodTypes))
    names(eval_acc[[i]]) <- methodTypes
    names(eval_roc[[i]]) <- methodTypes
    
    for(j in 1:length(methodTypes)) {
      model <- train(Class~., data=input_data, method=methodTypes[j], trControl = train_control)
      pred_result <- predict(model, newdata = test_data)
      eval_acc[[i]][[j]] <- round((sum(pred_result == test_data$Class) / nrow(test_data)), 3)
      pred_result <- predict(model, newdata = test_data, type = "prob")
      eval_roc[[i]][[j]] <- roc(test_data$Class, pred_result$YES)
      gc()
    }
  }
  
  ### create a ROC Curve from the 5k CV
  p <- NULL
  for(i in 1:length(methodTypes)) {
    temp_roc_list <- vector("list", 5)
    temp_acc <- NULL
    temp_auc <- NULL
    for(j in 1:5) {
      temp_roc_list[[j]] <- eval_roc[[j]][[i]]
      temp_acc <- c(temp_acc, eval_acc[[j]][[i]])
      temp_auc <- c(temp_auc, temp_roc_list[[j]]$auc)
      names(temp_roc_list)[j] <- paste0(j, " (ACC = ", round(temp_acc[j], digits = 3),
                                        ", AUC = ", round(temp_auc[j], digits = 3), ")")
    }
    p[[i]] <- ggroc(temp_roc_list, legacy.axes = TRUE, aes=c("linetype", "color"), size = 2) +
      ggtitle(paste0(methodNames[i], " (Avg ACC = ", round(mean(temp_acc), digits = 3),
                     ", Avg AUC = ", round(mean(temp_auc), digits = 3), ")")) +
      geom_abline(intercept=0, slope=1, color="black", size=1) +
      scale_color_npg() +
      labs(color="Iteration (k=5)", linetype="Iteration (k=5)") +
      theme_classic(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  ### arrange the plots and save
  fName <- paste0("Px08_Classifier_5k_Cross_Validation")
  rowNum <- 3
  colNum <- 2
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir3, fName, ".png"), g, width = 20, height = 12, dpi = 300)
  
  
  #
  ### GSEA with the important genes of the PC1
  #
  
  #'****************************************************************************************
  #' Gene Set Enrichment Analysis function
  #' 
  #' It receives gene list (character vector) and signature profiles (named numeric vector)
  #' as inputs, performs GSEA and returns a table of GSEA result table and draws
  #' a GSEA plot. It is basically a statistical significance test to check how the
  #' given given genes are biased on the signature profiles.
  #' 
  #' Whether there are multiple gene sets or multiple signatures,
  #' multiple testing (FDR computation) is performed.
  #' But if the input gene set and the input signature are both lists with multiple
  #' items (The length of the two are both more than 1) then we return an error message.
  #' 
  #' The plot file names will be determined by names(gene_list) or names(signature)
  #' If length(gene_list) > 1, then names(gene_list) will be used and
  #' if length(signature) > 1, then names(signature) will be used as file names.
  #' If there is no list names, then file names will be "GSEA_Plot_i.png".
  #' Here, i indicates that the plot is from the i-th row of the GSEA result table.
  #' 
  #' * Some plot drawing codes were from Rtoolbox/R/ReplotGSEA.R written by Thomas Kuilman. 
  #'****************************************************************************************
  #' @title	run_gsea
  #' 
  #' @param gene_list   A list of character vectors containing gene names to be tested
  #' @param signature   A list of named numeric vectors of signature values for GSEA. The gene_list
  #'                    should be included in the names(signature)
  #' @param printPlot   If TRUE, it also generates GSEA plot of the results
  #'                    (Default = FALSE)
  #' @param fdr_cutoff  When printing GSEA plots, print them with the FDR < fdr_cutoff only
  #'                    (Default = 0.05)
  #' @param heatmap_color_type  when 'relative', the heatmap of the GSEA colors the bottom half of the
  #'                            absolute range of the signature as blue and the upper half as red
  #'                            when 'absolute', the heatmap of GSEA colors the negative signature as blue
  #'                            and the positives as red
  #' @param printPath   When printing GSEA plots, print them in the designated path
  #'                    (Default = "./")
  #' @param width       The width of the plot file
  #'                    (Default = 2000)
  #' @param height      The height of the plot file
  #'                    (Default = 1200)
  #' @param res         The resolution of the plot file
  #'                    (Default = 130)
  #' 
  #' @return 	          It tests bias of the "gene_list" on the "signature" range and
  #'                    returns a table including p-values and FDRs (adjusted p-values)
  #'                    If fdr_cutoff == TRUE, it also generates a GSEA plot with the result
  #' 
  run_gsea <- function(gene_list,
                       signature,
                       printPlot = FALSE,
                       fdr_cutoff = 0.05,
                       heatmap_color_type = c("relative", "absolute"),
                       width = 2000,
                       height = 1200,
                       res = 130,
                       printPath = "./") {
    
    ### load required libraries
    if(!require("fgsea", quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("fgsea")
      require("fgsea", quietly = TRUE)
    }
    if(!require("limma", quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("limma")
      require("limma", quietly = TRUE)
    } 
    if(!require("checkmate", quietly = TRUE)) {
      install.packages("checkmate")
      require("checkmate", quietly = TRUE)
    }
    
    ### argument checking
    assertList(gene_list)
    assertList(signature)
    assertLogical(printPlot)
    assertNumeric(fdr_cutoff)
    assertIntegerish(width)
    assertIntegerish(height)
    assertIntegerish(res)
    assertString(printPath)
    if(length(gene_list) > 1 && length(signature) > 1) {
      stop("ERROR: \"gene_list\" and \"signature\" cannot be both \"list\"")
    }
    
    ### set random seed
    set.seed(1234)
    
    ### run GSEA
    ### if there are more than one signatures
    if(length(signature) > 1) {
      ### combine GSEA results of every signature inputs
      for(i in 1:length(signature)) {
        temp <- data.frame(fgsea(pathways = gene_list, stats = signature[[i]], nperm = 1000))
        if(i == 1) {
          gsea_result <- temp
        } else {
          gsea_result <- rbind(gsea_result, temp)
        }
      }
      
      ### compute FDRs
      corrected_gsea_result <- gsea_result[order(gsea_result$pval),]
      corrected_gsea_result$padj <- p.adjust(corrected_gsea_result$pval, method = "BH")
      gsea_result <- corrected_gsea_result[rownames(gsea_result),]
    } else {
      ### if there are more than one gene sets
      gsea_result <- data.frame(fgsea(pathways = gene_list, stats = signature[[1]], nperm = 1000))
    }
    
    ### print GSEA plot
    sIdx <- which(gsea_result$padj < fdr_cutoff)
    if(printPlot && length(sIdx) > 0) {
      for(i in sIdx) {
        ### get required values ready
        if(length(signature) > 1) {
          geneset <- gene_list[[1]]
          stats <- signature[[i]]
          stats <- stats[order(-stats)]
          fileName <- names(signature)[i]
        } else {
          geneset <- gene_list[[i]]
          stats <- signature[[1]]
          stats <- stats[order(-stats)]
          fileName <- names(gene_list)[i]
        }
        if(is.null(fileName)) {
          fileName <- paste0("GSEA_Plot_", i)
        }
        stats <- stats[!is.na(stats)]
        gsea.hit.indices <- which(names(stats) %in% geneset)
        es.temp <- calcGseaStat(stats, gsea.hit.indices, returnAllExtremes = TRUE)
        if(es.temp$res >= 0) {
          gsea.es.profile <- es.temp$tops
        } else {
          gsea.es.profile <- es.temp$bottoms
        }
        enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
        metric.range <- c(min(stats), max(stats))
        gsea.p.value <- round(gsea_result$pval[i] ,5)
        gsea.fdr <- round(gsea_result$padj[i] ,5)
        gsea.enrichment.score <- round(gsea_result$ES[i], 5)
        gsea.normalized.enrichment.score <- round(gsea_result$NES[i], 5)
        
        ### print GSEA result plot
        png(paste0(printPath, fileName, ".png"), width = width, height = height, res = res)
        
        ### set layout
        layout.show(layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2)))
        
        ### draw the GSEA plot
        par(mar = c(0, 5, 2, 2))
        plot(c(1, gsea.hit.indices, length(stats)),
             c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
             xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
             ylim = enrichment.score.range,
             main = list(fileName, font = 1, cex = 1),
             panel.first = {
               abline(h = seq(round(enrichment.score.range[1], digits = 1),
                              enrichment.score.range[2], 0.1),
                      col = "gray95", lty = 2)
               abline(h = 0, col = "gray50", lty = 2)
             }
        )
        
        ### add informative text to the GSEA plot
        plot.coordinates <- par("usr")
        if(es.temp$res < 0) {
          text(length(stats) * 0.01, plot.coordinates[3] * 0.98,
               paste("P-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
                     gsea.enrichment.score, "\nNormalized ES:",
                     gsea.normalized.enrichment.score), adj = c(0, 0))
        } else {
          text(length(stats) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
               paste("P-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
                     gsea.enrichment.score, "\nNormalized ES:",
                     gsea.normalized.enrichment.score), adj = c(1, 1))
        }
        
        ### draw hit indices
        par(mar = c(0, 5, 0, 2))
        plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
             ylab = "", xlim = c(1, length(stats)))
        abline(v = gsea.hit.indices, lwd = 0.75)
        
        ### create color palette for the heatmap
        par(mar = c(0, 5, 0, 2))
        if(heatmap_color_type[1] == "relative") {
          rank.colors <- stats - metric.range[1]
          rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
          rank.colors <- ceiling(rank.colors * 511 + 1)
          rank.colors <- colorRampPalette(c("blue", "white", "red"))(512)[rank.colors]
        } else {
          rank.colors1 <- stats[which(stats >= 0)]
          rank.colors1 <- rank.colors1 - min(rank.colors1)
          rank.colors1 <- rank.colors1 / (max(rank.colors1) - min(rank.colors1))
          rank.colors1 <- ceiling(rank.colors1 * 255 + 1)
          rank.colors1 <- colorRampPalette(c("white", "red"))(256)[rank.colors1]
          rank.colors2 <- stats[which(stats < 0)]
          rank.colors2 <- rank.colors2 - min(rank.colors2)
          rank.colors2 <- rank.colors2 / (max(rank.colors2) - min(rank.colors2))
          rank.colors2 <- ceiling(rank.colors2 * 255 + 1)
          rank.colors2 <- colorRampPalette(c("blue", "white"))(256)[rank.colors2]
          rank.colors <- c(rank.colors1, rank.colors2)
        }
        
        ### draw the heatmap
        rank.colors <- rle(rank.colors)
        barplot(matrix(rank.colors$lengths), col = rank.colors$values,
                border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(stats)))
        box()
        text(length(stats) / 2, 0.7,
             labels = "Signature")
        text(length(stats) * 0.01, 0.7, "Largest", adj = c(0, NA))
        text(length(stats) * 0.99, 0.7, "Smallest", adj = c(1, NA))
        
        ### draw signature values
        par(mar = c(5, 5, 0, 2))
        rank.metric <- rle(round(stats, digits = 2))
        plot(stats, type = "n", xaxs = "i",
             xlab = "Rank in ordered gene list", xlim = c(0, length(stats)),
             ylim = metric.range, yaxs = "i",
             ylab = "Signature values",
             panel.first = abline(h = seq(metric.range[1] / 2,
                                          metric.range[2] - metric.range[1] / 4,
                                          metric.range[2] / 2), col = "gray95", lty = 2))
        
        barplot(rank.metric$values, col = "lightgrey", lwd = 0.1,
                xlim = c(0, length(stats)), ylim = c(-1, 1),
                width = rank.metric$lengths, border = NA,
                space = 0, add = TRUE, xaxt = "n")
        box()
        
        ### print out the file
        dev.off()
      }
    }
    
    return(gsea_result)
    
  }
  
  ### db preparation
  # MSIGDB
  m_df <- msigdbr(species = "Homo sapiens")
  m_list <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  #
  ### Best Predictor Patient's GMP CAR+ CD8 cells vs those of every others
  #
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  all_gmp_last <- intersect(all_gmp_last,
                            which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  all_gmp_not_last <- intersect(all_gmp_not_last,
                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### only get the cells of interests
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[c(all_gmp_last, all_gmp_not_last)])
  
  ### set new column for DE analysis comparison
  target_Seurat_Obj@meta.data$New_Group <- paste0(target_Seurat_Obj@meta.data$ALL_CARpos_Persister,
                                                  "_",
                                                  target_Seurat_Obj@meta.data$px)
  
  ### set idents with the new info
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$New_Group)
  
  ### DE analysis for all the comparisons
  for(px in setdiff(unique(target_Seurat_Obj@meta.data$px), "SJCAR19-06")) {
    
    ### DE analysis for persisters
    de_result <- FindMarkers(target_Seurat_Obj,
                             ident.1 = "YES_SJCAR19-06",
                             ident.2 = paste0("YES_", px),
                             min.pct = 0.1,
                             logfc.threshold = 0.1,
                             test.use = "wilcox")
    
    ### write out the DE result
    write.xlsx2(data.frame(Gene=rownames(de_result),
                           de_result,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir, "/GMP_CARpos_CD8_Persister_Best_vs_Others.xlsx"),
                sheetName = paste0("SJCAR19-06_vs_", px), row.names = FALSE, append = TRUE)
    
    ### pathway analysis
    pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "GO",
                                            title = paste0("Persisters_SJCAR19-06_vs_", px),
                                            displayNum = 30, imgPrint = TRUE,
                                            dir = paste0(outputDir))
    pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                                "ENTREZID", "SYMBOL"),
                                              org = "human", database = "KEGG",
                                              title = paste0("Persisters_SJCAR19-06_vs_", px),
                                              displayNum = 30, imgPrint = TRUE,
                                              dir = paste0(outputDir))
    if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
      write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_Persisters_SJCAR19-06_vs_", px, ".xlsx"),
                  row.names = FALSE, sheetName = paste0("GO_Result"))
    }
    if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
      write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_Persisters_SJCAR19-06_vs_", px, ".xlsx"),
                  row.names = FALSE, sheetName = paste0("KEGG_Result"))
    }
    
    #
    ### GSEA
    #
    output_dir <- paste0(outputDir, "Persister_SJCAR19-06_vs_", px, "/")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    ### signature preparation
    signat <- de_result$avg_logFC
    names(signat) <- rownames(de_result)
    
    ### run GSEA
    GSEA_result <- run_gsea(gene_list = m_list, signature = list(signat), printPlot = FALSE)
    GSEA_result <- GSEA_result[order(GSEA_result$pval),]
    
    ### only get pathways that have pval < 0.01 & size > 30 & up-regulating results (enriched with important) only
    pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                        which(GSEA_result$size > 30)),
                                              which(GSEA_result$NES > 2))]
    if(length(pathways) < 5) {
      pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                          which(GSEA_result$size > 30)),
                                                which(GSEA_result$NES > 1.8))]
    }
    if(length(pathways) < 5) {
      pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                          which(GSEA_result$size > 30)),
                                                which(GSEA_result$NES > 1.6))]
    }
    if(length(pathways) > 30) {
      pathways <- GSEA_result$pathway[1:30]
    }
    
    ### run GSEA again with the significant result - plot printing
    GSEA_result2 <- run_gsea(gene_list = m_list[pathways], signature = list(signat),
                             printPlot = TRUE, printPath = output_dir)
    
    ### write out the result file
    write.xlsx2(GSEA_result, file = paste0(output_dir, "/GSEA_GMP_CARpos_CD8_Persister_Best_vs_Others.xlsx"),
                sheetName = paste0("SJCAR19-06_vs_", px), row.names = FALSE, append = TRUE)
    
    
    ### DE analysis for non-persisters
    de_result <- FindMarkers(target_Seurat_Obj,
                             ident.1 = "NO_SJCAR19-06",
                             ident.2 = paste0("NO_", px),
                             min.pct = 0.1,
                             logfc.threshold = 0.1,
                             test.use = "wilcox")
    
    ### write out the DE result
    write.xlsx2(data.frame(Gene=rownames(de_result),
                           de_result,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir, "/GMP_CARpos_CD8_Non-Persister_Best_vs_Others.xlsx"),
                sheetName = paste0("SJCAR19-06_vs_", px), row.names = FALSE, append = TRUE)
    
    ### pathway analysis
    pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "GO",
                                            title = paste0("Non-Persisters_SJCAR19-06_vs_", px),
                                            displayNum = 30, imgPrint = TRUE,
                                            dir = paste0(outputDir))
    pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                rownames(de_result)[which(de_result$p_val_adj < 0.05)],
                                                                "ENTREZID", "SYMBOL"),
                                              org = "human", database = "KEGG",
                                              title = paste0("Non-Persisters_SJCAR19-06_vs_", px),
                                              displayNum = 30, imgPrint = TRUE,
                                              dir = paste0(outputDir))
    if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
      write.xlsx2(pathway_result_GO, file = paste0(outputDir, "GO_Non-Persisters_SJCAR19-06_vs_", px, ".xlsx"),
                  row.names = FALSE, sheetName = paste0("GO_Result"))
    }
    if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
      write.xlsx2(pathway_result_KEGG, file = paste0(outputDir, "KEGG_Non-Persisters_SJCAR19-06_vs_", px, ".xlsx"),
                  row.names = FALSE, sheetName = paste0("KEGG_Result"))
    }
    
    #
    ### GSEA
    #
    output_dir <- paste0(outputDir, "Non-Persister_SJCAR19-06_vs_", px, "/")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    ### signature preparation
    signat <- de_result$avg_logFC
    names(signat) <- rownames(de_result)
    
    ### run GSEA
    GSEA_result <- run_gsea(gene_list = m_list, signature = list(signat), printPlot = FALSE)
    GSEA_result <- GSEA_result[order(GSEA_result$pval),]
    
    ### only get pathways that have pval < 0.01 & size > 30 & up-regulating results (enriched with important) only
    pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                        which(GSEA_result$size > 30)),
                                              which(GSEA_result$NES > 2))]
    if(length(pathways) < 5) {
      pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                          which(GSEA_result$size > 30)),
                                                which(GSEA_result$NES > 1.8))]
    }
    if(length(pathways) < 5) {
      pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                          which(GSEA_result$size > 30)),
                                                which(GSEA_result$NES > 1.6))]
    }
    if(length(pathways) > 30) {
      pathways <- GSEA_result$pathway[1:30]
    }
    
    ### run GSEA again with the significant result - plot printing
    GSEA_result2 <- run_gsea(gene_list = m_list[pathways], signature = list(signat),
                             printPlot = TRUE, printPath = output_dir)
    
    ### write out the result file
    write.xlsx2(GSEA_result, file = paste0(output_dir, "/GSEA_GMP_CARpos_CD8_Non-Persister_Best_vs_Others.xlsx"),
                sheetName = paste0("SJCAR19-06_vs_", px), row.names = FALSE, append = TRUE)
    
  }
  
  ### write out the number of cells in each comparison
  for(px in unique(target_Seurat_Obj@meta.data$px)) {
    writeLines(paste(px, "\n# Persisters: ", length(which(target_Seurat_Obj@meta.data$New_Group == paste0("YES_", px))),
                     ", # Non-Persisters:", length(which(target_Seurat_Obj@meta.data$New_Group == paste0("NO_", px))), "\n"))
  }
  
  
  #'****************************************************************************************
  #' Compare different ranks of the same elements
  #' 
  #' It receives (1) two named numeric vectors that contain two different rank metrics
  #' OR (2) receives two character vectors of the ordered/ranked element names
  #' 
  #' In case of (1), if the two given vectors are A and B,
  #' names(A) and names(B) should have the same elements.
  #' length(intersect(names(A), names(B))) = length(A) = length(B).
  #' In case of (2), if the two given vectors are A and B,
  #' A and B should have the same elements.
  #' length(intersect(A, B)) = length(A) = length(B).
  #' 
  #' The function will generate 4 figures in one plot.
  #' If the two given vectors are A and B,
  #' 1. a barplot of ranked metric of A
  #' 2. a line plot that shows lines connect top same elements between ranked A and ranked B
  #' 3. a barplot of ranked metric of B
  #' 4. a scatter plot of A and B
  #'****************************************************************************************
  #' @title	compare_two_different_ranks
  #' 
  #' @param A           A vector of named numeric vector or a character vector
  #' @param B           A vector of named numeric vector or a character vector
  #' @param A_name      The title of the given rank A
  #'                    (Default = "A")      
  #' @param B_name      The title of the given rank B
  #'                    (Default = "B")
  #' @param ordering    Order of the rank metric ["decreasing" or "increasing"]
  #'                    (Default = "decreasing")
  #' @param top         The number of top hubs from each of A and B that will be
  #'                    presented in the line graph
  #'                    (Default = 100)
  #' @param alternative The alternative hypothesis for enrichment permutation ["greater" or "less"]
  #'                    (Default = "greater")
  #' @param permutation The number of permutation test
  #'                    (Default = 10000)
  #' @param fileName    The name of the plot file
  #'                    (Default = "Comparison_of_Ranks_Between_A_and_B")
  #' @param printPath   Print the plot in the designated path
  #'                    (Default = "./")
  #' @param width       The width of the plot file
  #'                    (Default = 24)
  #' @param height      The height of the plot file
  #'                    (Default = 12)
  #' 
  #' @return 	          It does not return any objects
  #' 
  compare_two_different_ranks <- function(A,
                                          B,
                                          A_name = "A",
                                          B_name = "B",
                                          ordering = "decreasing",
                                          top = 100,
                                          alternative = "greater",
                                          permutation = 10000,
                                          fileName = "Comparison_of_Ranks_Between_A_and_B",
                                          printPath = "./",
                                          width = 24,
                                          height = 10) {
    
    ### load required libraries
    if(!require(ggplot2, quietly = TRUE)) {
      install.packages("ggplot2")
      require(ggplot2, quietly = TRUE)
    }
    if(!require(gridExtra, quietly = TRUE)) {
      install.packages("gridExtra")
      require(gridExtra, quietly = TRUE)
    }
    if(!require("checkmate", quietly = TRUE)) {
      install.packages("checkmate")
      require("checkmate", quietly = TRUE)
    }
    
    ### argument checking
    assertVector(A)
    assertVector(B)
    assertString(A_name)
    assertString(B_name)
    assertChoice(ordering, c("decreasing", "increasing"))
    assertIntegerish(top)
    assertChoice(alternative, c("greater", "less"))
    assertIntegerish(permutation)
    assertString(fileName)
    assertString(printPath)
    assertIntegerish(width)
    assertIntegerish(height)
    if(length(A) != length(B)) {
      stop("ERROR: length(A) != length(B)")
    }
    if(class(A) == "numeric" && class(B) == "numeric") {
      if(length(intersect(names(A), names(B))) != length(A)) {
        stop("ERROR: A and B should have names, and names(sort(A)) and names(sort(B)) should be the same")
      }
      type <- "numeric"
    } else if(class(A) == "character" && class(B) == "character") {
      if(length(intersect(A, B)) != length(A)) {
        stop("ERROR: sort(A) and sort(B) should be the same")
      }
      type <- "character"
      names(A) <- A
      A <- length(A):1
      names(B) <- B
      B <- length(B):1
    } else {
      stop("ERROR: A and B should be either named numeric vectors or character vectors")
    }
    
    ### order the A and B
    if(ordering == "decreasing") {
      A <- A[order(-A)]
      B <- B[order(-B)]
    } else {
      A <- A[order(A)]
      B <- B[order(B)]
    }
    
    ### data frame of A and B
    df <- data.frame(idx=1:length(A), A=A[order(names(A))], B=B[order(names(B))],
                     stringsAsFactors = FALSE, check.names = FALSE)
    
    ### give colors to the df
    colors <- rep("lightgray", nrow(df))
    names(colors) <- rownames(df)
    colors[names(A)[1:top]] <- "skyblue"
    colors[names(B)[1:top]] <- "pink"
    colors[intersect(names(A)[1:top], names(B)[1:top])] <- "mediumpurple"
    df$colors <- factor(colors, levels = c("skyblue", "pink", "lightgray", "mediumpurple"))
    
    ### data frame for the line graph
    df2 <- data.frame(X=c(Reduce(function(x, y) {
      return(c(x, y, which(names(B) == names(A)[y])))
    }, 0:top)[-1],
    Reduce(function(x, y) {
      return(c(x, y, which(names(A) == names(B)[y])))
    }, 0:top)[-1]),
    Y=c(rep(c(1,0), top),
        rep(c(0,1), top)),
    Element=paste0("e", c(rbind(1:(top*2), 1:(top*2)))),
    Collection=c(rep("A", top*2),
                 rep("B", top*2)),
    stringsAsFactors = FALSE, check.names = FALSE)
    
    ### add top and bottom x-axis to the df2
    df2 <- rbind(df2, data.frame(X=c(0, length(A), 0, length(B)),
                                 Y=c(1, 1, 0, 0),
                                 Element=paste0("e", c(rbind((top*2+1):(top*2+2),(top*2+1):(top*2+2)))),
                                 Collection=rep("C", 4)))
    
    ### permutation test for getting a p-value (1-tail)
    set.seed(1234)
    permu_result <- sapply(1:(permutation-1), function(x) {
      random_hubs <- names(A)[sample(length(A), top)]
      return(sum(B[random_hubs]))
    })
    if(alternative == "greater") {
      pVal_A <- (length(which(permu_result > sum(B[names(A)[1:top]])))+1) / permutation
    } else {
      pVal_A <- (length(which(permu_result < sum(B[names(A)[1:top]])))+1) / permutation
    }
    permu_result <- sapply(1:(permutation-1), function(x) {
      random_hubs <- names(B)[sample(length(B), top)]
      return(sum(A[random_hubs]))
    })
    if(alternative == "greater") {
      pVal_B <- (length(which(permu_result > sum(A[names(B)[1:top]])))+1) / permutation
    } else {
      pVal_B <- (length(which(permu_result < sum(A[names(B)[1:top]])))+1) / permutation
    }
    
    if(ordering == "decreasing") {
      ### 1. a barplot of ranked metric of A
      p1 <- ggplot(data = df, aes(x = reorder(idx, -A), y = A, width = 1)) +
        ylab(paste(A_name, "(A)")) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = margin(t = 0, r = 115, b = 0, l = 5)) +
        geom_bar(stat = "identity", color = "skyblue")
      
      ### 3. a barplot of ranked metric of B
      p3 <- ggplot(data = df, aes(x = reorder(idx, -B), y = B, width = 1)) +
        ylab(paste(B_name, "(B)")) +
        theme(axis.title.x=element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = margin(t = 0, r = 115, b = 0, l = 5)) +
        geom_bar(stat = "identity", color = "pink")
    } else {
      ### 1. a barplot of ranked metric of A
      p1 <- ggplot(data = df, aes(x = reorder(idx, A), y = A, width = 1)) +
        ylab(paste(A_name, "(A)")) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = margin(t = 0, r = 115, b = 0, l = 5)) +
        geom_bar(stat = "identity", color = "skyblue")
      
      ### 3. a barplot of ranked metric of B
      p3 <- ggplot(data = df, aes(x = reorder(idx, B), y = B, width = 1)) +
        ylab(paste(B_name, "(B)")) +
        theme(axis.title.x=element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = margin(t = 0, r = 115, b = 0, l = 5)) +
        geom_bar(stat = "identity", color = "pink")
    }
    
    ### 2. a line plot that shows lines connect top same elements between ranked A and ranked B
    p2 <- ggplot(data = df2, aes(x = X, y = Y, group = Element, color = Collection)) +
      geom_line() +
      scale_color_manual(labels = c("A", "B", "C"), values = c("skyblue", "pink", "gray")) +
      theme_void() +
      theme(legend.position = "none", plot.margin = margin(t = 0, r = 5, b = 0, l = 5)) +
      annotate(geom = "text", x = length(A), y = 0.1,
               color = "skyblue", hjust = 0, vjust = 0, fontface = "bold",
               label = paste("Top", top, "(A)\nPermutation p-value\n=", pVal_A)) +
      annotate(geom = "text", x = length(B), y = 0.9,
               color = "pink", hjust = 0, vjust = 1, fontface = "bold",
               label = paste("Top", top, "(B)\nPermutation p-value\n=", pVal_B)) +
      expand_limits(x = length(A)*1.1)
    
    ### 4. a scatter plot of A and B
    p4 <- ggplot(data = df, aes(x=A, y=B)) +
      geom_point(aes(color=colors), size = 2) +
      labs(subtitle=paste0("Pearson Correlation = ", round(cor(df$A, df$B), 5),
                           " (P-value = ", round(cor.test(df$A, df$B)$p.value, 5), ")"),
           color=paste("Top", top)) +
      xlab(A_name) +
      ylab(B_name) +
      geom_smooth(method = lm, color="black", se=FALSE) +
      scale_color_manual(labels=c("A", "B", "Others", "Intersect(A, B)"), values = c("skyblue", "pink", "lightgray", "mediumpurple")) +
      theme_classic(base_size = 16)
    
    ### arrange the plots and print out
    g <- arrangeGrob(p1, p2, p3, p4, layout_matrix = rbind(c(1, 4),
                                                           c(2, 4),
                                                           c(3, 4)),
                     top = fileName)
    ggsave(file = paste0(printPath, fileName, ".png"), g, width = width, height = height)
    
  }
  
  #
  ### what are the DE genes between PS & NPS in each patient,
  ### and are they different from each other?
  #
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  all_gmp_last <- intersect(all_gmp_last,
                            which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  all_gmp_not_last <- intersect(all_gmp_not_last,
                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### only get the cells of interests
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[c(all_gmp_last, all_gmp_not_last)])
  
  ### set idents with the px info
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  
  ### run DE analysis for each patient
  de_results <- vector("list", length(unique(target_Seurat_Obj@meta.data$px)))
  names(de_results) <- unique(target_Seurat_Obj@meta.data$px)
  for(px in unique(target_Seurat_Obj@meta.data$px)) {
    ### get the patient data
    px_Seurat_Obj <- subset(target_Seurat_Obj, idents = px)
    
    ### set idents for DEA
    px_Seurat_Obj <- SetIdent(object = px_Seurat_Obj,
                              cells = rownames(px_Seurat_Obj@meta.data),
                              value = px_Seurat_Obj@meta.data$ALL_CARpos_Persister)
    
    ### run DEA
    de_results[[px]] <- FindMarkers(px_Seurat_Obj,
                                    ident.1 = "YES",
                                    ident.2 = "NO",
                                    min.pct = 0.1,
                                    logfc.threshold = 0.1,
                                    test.use = "wilcox")
    
    ### add order
    de_results[[px]]$rank <- rank(de_results[[px]]$p_val_adj)
  }
  
  ### pairwise comparison of DE results
  for(i in 1:(length(de_results)-1)) {
    for(j in (i+1):length(de_results)) {
      ### get common genes
      common_genes <- intersect(rownames(de_results[[i]])[which(de_results[[i]]$p_val_adj < 1)],
                                rownames(de_results[[j]])[which(de_results[[j]]$p_val_adj < 1)])
      
      ### get signatures - RANK
      sig1 <- de_results[[i]][common_genes,"rank"]
      names(sig1) <- common_genes
      sig2 <- de_results[[j]][common_genes,"rank"]
      names(sig2) <- common_genes
      
      ### compare the ranks
      if(length(sig1) > 2 && length(sig2) > 2) {
        compare_two_different_ranks(A = sig1,
                                    B = sig2,
                                    A_name = paste0("DE Gene Rank of", names(de_results)[i]),
                                    B_name = paste0("DE Gene Rank of", names(de_results)[j]),
                                    ordering = "decreasing",
                                    top = length(sig1),
                                    alternative = "greater",
                                    permutation = 10000,
                                    fileName = paste0("DE_Gene_Rank_Comparison_", names(de_results)[i], "_vs_", names(de_results)[j]),
                                    printPath = outputDir,
                                    width = 24,
                                    height = 12)
      } else {
        writeLines(paste(names(de_results)[i], "&", names(de_results)[j], "does not have enough genes for the analysis."))
      }
    }
  }
  
  #
  ### UMAP plot of every px in one - persister vs non-persister
  ### the purpose of it is to know how persisters and non-persisters distributed in each patient
  ### and compre those - each patient should be distinguished by color and the persistence
  ### should be also distinguished by shape
  #
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  all_gmp_last <- intersect(all_gmp_last,
                            which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  all_gmp_not_last <- intersect(all_gmp_not_last,
                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### only get the cells of interests
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[c(all_gmp_last, all_gmp_not_last)])
  
  ### UMAP with all those info
  p <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
               group.by = "px", shape.by = "GMP_CARpos_Persister",
               pt.size = 2) +
    ggtitle("UMAP of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Patient",
         shape="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.5
  ggsave(paste0(outputDir, "/", "UMAP_Plot_Persistence_per_Px.png"), plot = p, width = 20, height = 12, dpi = 300)
  
  ### get the same number of samples (persister: non-persister - same ratio) and draw a UMAP plot again
  set.seed(1234)
  target_Seurat_Obj@meta.data$New_Persistency <- NA
  multiplier_k <- 1
  for(px in unique(target_Seurat_Obj@meta.data$px)) {
    ### get specific indicies
    px_gmp_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                             which(target_Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES"))
    px_gmp_not_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                                 which(target_Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO"))
    
    ### sampling
    if((length(px_gmp_last) > 2) && (length(px_gmp_not_last) > length(px_gmp_last))) {
      px_gmp_not_last <- sample(px_gmp_not_last, size = length(px_gmp_last)*multiplier_k)
    }
    
    ### annotate new persistency info
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_last] <- "YES"
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_not_last] <- "NO"
  }
  
  ### set idents with the new info
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$New_Persistency)
  
  ### temporary seurat object for umap
  umap_seurat_obj <- subset(target_Seurat_Obj, idents = c("YES", "NO"))
  
  ### run pca & umap
  umap_seurat_obj <- RunPCA(umap_seurat_obj,
                            features = VariableFeatures(object = umap_seurat_obj),
                            npcs = 15)
  umap_seurat_obj <- RunUMAP(umap_seurat_obj, dims = 1:15)
  
  ### UMAP with all those info
  p <- list()
  p[[1]] <- DimPlot(object = umap_seurat_obj, reduction = "umap",
               group.by = "px", shape.by = "New_Persistency",
               pt.size = 3) +
    ggtitle("UMAP of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Patient",
         shape="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir, "/", "UMAP_Plot_Persistence_per_Px_Sampled.png"), plot = p[[1]], width = 20, height = 12, dpi = 300)
  
  ### umap for each px
  umap_seurat_obj <- SetIdent(object = umap_seurat_obj,
                              cells = rownames(umap_seurat_obj@meta.data),
                              value = umap_seurat_obj@meta.data$px)
  g <- ggplot_build(p[[1]])
  color_code <- data.frame(colours = unique(g$data[[1]]["colour"]), 
                           label = levels(g$plot$data[, "px"]),
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(color_code) <- color_code$label
  for(i in 1:length(unique(umap_seurat_obj@meta.data$px))) {
    temp_seurat_obj <- subset(umap_seurat_obj, idents = c(unique(umap_seurat_obj@meta.data$px)[i]))
    p[[i+1]] <- DimPlot(object = temp_seurat_obj, reduction = "umap",
                 cols = color_code[unique(umap_seurat_obj@meta.data$px)[i],"colour"], shape.by = "New_Persistency",
                 pt.size = 3) +
      ggtitle(unique(umap_seurat_obj@meta.data$px)[i]) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
      labs(shape="Is Persistent")
    p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  ### arrange the plots and save
  fName <- paste0("UMAP_Plot_Persistence_per_Px_Sampled")
  rowNum <- 3
  colNum <- 4
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir, fName, "2.png"), g, width = 25, height = 15, dpi = 300)
  
  ### PCA with all those info
  p <- list()
  p[[1]] <- DimPlot(object = umap_seurat_obj, reduction = "pca",
                    group.by = "px", shape.by = "New_Persistency",
                    pt.size = 3) +
    ggtitle("PCA of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Patient",
         shape="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir, "/", "PCA_Plot_Persistence_per_Px_Sampled.png"), plot = p[[1]], width = 20, height = 12, dpi = 300)
  
  ### PCA function from my own
  umap_seurat_obj <- SetIdent(object = umap_seurat_obj,
                              cells = rownames(umap_seurat_obj@meta.data),
                              value = umap_seurat_obj@meta.data$New_Persistency)
  de_result <- FindMarkers(umap_seurat_obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  pca_plot(normalizedMat = normalizeRNASEQwithVST(umap_seurat_obj@assays$RNA@counts[rownames(de_result)[1:100],] + 1,
                                                  filter_thresh = 0),
           grp = as.character(umap_seurat_obj@meta.data$px),
           grp2 = as.character(umap_seurat_obj@meta.data$New_Persistency),
           title = paste0("PCA_Plot_Persistence_per_Px_Sampled_DE_100"),
           outDir = outputDir)
  pca_plot(normalizedMat = normalizeRNASEQwithVST(umap_seurat_obj@assays$RNA@counts[rownames(de_result)[1:100],] + 1,
                                                  filter_thresh = 0),
           grp = as.character(umap_seurat_obj@meta.data$New_Persistency),
           title = paste0("PCA_Plot_Persistence_per_Px_Sampled_DE_100(2)"),
           outDir = outputDir)
  
  ### pca for each px
  umap_seurat_obj <- SetIdent(object = umap_seurat_obj,
                              cells = rownames(umap_seurat_obj@meta.data),
                              value = umap_seurat_obj@meta.data$px)
  g <- ggplot_build(p[[1]])
  color_code <- data.frame(colours = unique(g$data[[1]]["colour"]), 
                           label = levels(g$plot$data[, "px"]),
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(color_code) <- color_code$label
  for(i in 1:length(unique(umap_seurat_obj@meta.data$px))) {
    temp_seurat_obj <- subset(umap_seurat_obj, idents = c(unique(umap_seurat_obj@meta.data$px)[i]))
    p[[i+1]] <- DimPlot(object = temp_seurat_obj, reduction = "pca",
                        cols = color_code[unique(umap_seurat_obj@meta.data$px)[i],"colour"], shape.by = "New_Persistency",
                        pt.size = 3) +
      ggtitle(unique(umap_seurat_obj@meta.data$px)[i]) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
      labs(shape="Is Persistent")
    p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  ### arrange the plots and save
  fName <- paste0("PCA_Plot_Persistence_per_Px_Sampled")
  rowNum <- 3
  colNum <- 4
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir, fName, "2.png"), g, width = 25, height = 15, dpi = 300)
  
  
  ### get the same number of samples (30) and draw a UMAP plot again
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[c(all_gmp_last, all_gmp_not_last)])
  set.seed(1234)
  target_Seurat_Obj@meta.data$New_Persistency <- NA
  sample_num <- 30
  for(px in unique(target_Seurat_Obj@meta.data$px)) {
    ### get specific indicies
    px_gmp_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                             which(target_Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES"))
    px_gmp_not_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                                 which(target_Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO"))
    
    ### sampling
    if((length(px_gmp_last) > 2) && (length(px_gmp_not_last) > length(px_gmp_last))) {
      min_num <- min(length(px_gmp_last), sample_num)
      px_gmp_last <- sample(px_gmp_last, min_num)
      px_gmp_not_last <- sample(px_gmp_not_last, min_num)
    }
    
    ### annotate new persistency info
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_last] <- "YES"
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_not_last] <- "NO"
  }
  
  ### set idents with the new info
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$New_Persistency)
  
  ### temporary seurat object for umap
  umap_seurat_obj <- subset(target_Seurat_Obj, idents = c("YES", "NO"))
  
  ### run pca & umap
  umap_seurat_obj <- RunPCA(umap_seurat_obj,
                            features = VariableFeatures(object = umap_seurat_obj),
                            npcs = 15)
  umap_seurat_obj <- RunUMAP(umap_seurat_obj, dims = 1:15)
  
  ### UMAP with all those info
  p <- list()
  p[[1]] <- DimPlot(object = umap_seurat_obj, reduction = "umap",
                    group.by = "px", shape.by = "New_Persistency",
                    pt.size = 5) +
    ggtitle("UMAP of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Patient",
         shape="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir, "/", "UMAP_Plot_Persistence_per_Px_30_Sampled.png"), plot = p[[1]], width = 20, height = 12, dpi = 300)
  
  ### umap for each px
  umap_seurat_obj <- SetIdent(object = umap_seurat_obj,
                              cells = rownames(umap_seurat_obj@meta.data),
                              value = umap_seurat_obj@meta.data$px)
  g <- ggplot_build(p[[1]])
  color_code <- data.frame(colours = unique(g$data[[1]]["colour"]), 
                           label = levels(g$plot$data[, "px"]),
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(color_code) <- color_code$label
  for(i in 1:length(unique(umap_seurat_obj@meta.data$px))) {
    temp_seurat_obj <- subset(umap_seurat_obj, idents = c(unique(umap_seurat_obj@meta.data$px)[i]))
    p[[i+1]] <- DimPlot(object = temp_seurat_obj, reduction = "umap",
                        cols = color_code[unique(umap_seurat_obj@meta.data$px)[i],"colour"], shape.by = "New_Persistency",
                        pt.size = 5) +
      ggtitle(unique(umap_seurat_obj@meta.data$px)[i]) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
      labs(shape="Is Persistent")
    p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  ### arrange the plots and save
  fName <- paste0("UMAP_Plot_Persistence_per_Px_30_Sampled")
  rowNum <- 3
  colNum <- 4
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir, fName, "2.png"), g, width = 25, height = 15, dpi = 300)
  
  ### PCA with all those info
  p <- list()
  p[[1]] <- DimPlot(object = umap_seurat_obj, reduction = "pca",
                    group.by = "px", shape.by = "New_Persistency",
                    pt.size = 5) +
    ggtitle("PCA of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Patient",
         shape="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir, "/", "PCA_Plot_Persistence_per_Px_30_Sampled.png"), plot = p[[1]], width = 20, height = 12, dpi = 300)
  
  ### PCA function from my own
  umap_seurat_obj <- SetIdent(object = umap_seurat_obj,
                              cells = rownames(umap_seurat_obj@meta.data),
                              value = umap_seurat_obj@meta.data$New_Persistency)
  de_result <- FindMarkers(umap_seurat_obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  pca_plot(normalizedMat = normalizeRNASEQwithVST(umap_seurat_obj@assays$RNA@counts[rownames(de_result)[1:100],] + 1,
                                                  filter_thresh = 0),
           grp = as.character(umap_seurat_obj@meta.data$px),
           grp2 = as.character(umap_seurat_obj@meta.data$New_Persistency),
           title = paste0("PCA_Plot_Persistence_per_Px_30_Sampled_DE_100"),
           outDir = outputDir)
  pca_plot(normalizedMat = normalizeRNASEQwithVST(umap_seurat_obj@assays$RNA@counts[rownames(de_result)[1:100],] + 1,
                                                  filter_thresh = 0),
           grp = as.character(umap_seurat_obj@meta.data$New_Persistency),
           title = paste0("PCA_Plot_Persistence_per_Px_30_Sampled_DE_100(2)"),
           outDir = outputDir)
  
  ### pca for each px
  umap_seurat_obj <- SetIdent(object = umap_seurat_obj,
                              cells = rownames(umap_seurat_obj@meta.data),
                              value = umap_seurat_obj@meta.data$px)
  g <- ggplot_build(p[[1]])
  color_code <- data.frame(colours = unique(g$data[[1]]["colour"]), 
                           label = levels(g$plot$data[, "px"]),
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(color_code) <- color_code$label
  for(i in 1:length(unique(umap_seurat_obj@meta.data$px))) {
    temp_seurat_obj <- subset(umap_seurat_obj, idents = c(unique(umap_seurat_obj@meta.data$px)[i]))
    p[[i+1]] <- DimPlot(object = temp_seurat_obj, reduction = "pca",
                        cols = color_code[unique(umap_seurat_obj@meta.data$px)[i],"colour"], shape.by = "New_Persistency",
                        pt.size = 5) +
      ggtitle(unique(umap_seurat_obj@meta.data$px)[i]) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
      labs(shape="Is Persistent")
    p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  ### arrange the plots and save
  fName <- paste0("PCA_Plot_Persistence_per_Px_30_Sampled")
  rowNum <- 3
  colNum <- 4
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir, fName, "2.png"), g, width = 25, height = 15, dpi = 300)
  
  #
  ### Meeting with Jeremy on 02/23/2021,
  ### we decided to finalize the UMAP/PCA plots
  ### a) the PS & NPS should have the same # of cells
  ### b) there should be two versions for each: 1. with all the genes, 2. with DE genes only
  ### c) the first one is with all the pxs, then for individual plots, use color instead of shape
  ###
  ### Only for the Px06 (Best Predictor) & Px08 (Worst Predictor),
  ### run GSEA in each patient with the signature as logFC of DE genes from PS vs NPS
  ### and compare the GSEA results with rank comparison
  #
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  all_gmp_last <- intersect(all_gmp_last,
                            which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  all_gmp_not_last <- intersect(all_gmp_not_last,
                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### only get the cells of interests
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[c(all_gmp_last, all_gmp_not_last)])
  
  ### get the same number of samples (persister: non-persister - same ratio) and draw a UMAP plot again
  set.seed(1234)
  target_Seurat_Obj@meta.data$New_Persistency <- NA
  multiplier_k <- 1
  for(px in unique(target_Seurat_Obj@meta.data$px)) {
    ### get specific indicies
    px_gmp_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                             which(target_Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES"))
    px_gmp_not_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                                 which(target_Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO"))
    
    ### sampling
    if((length(px_gmp_last) > 2) && (length(px_gmp_not_last) > length(px_gmp_last))) {
      px_gmp_not_last <- sample(px_gmp_not_last, size = length(px_gmp_last)*multiplier_k)
    }
    
    ### annotate new persistency info
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_last] <- "YES"
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_not_last] <- "NO"
  }
  
  ### set idents with the new info
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$New_Persistency)
  
  #
  ### 1. with all the genes
  ### temporary seurat object for umap/pca
  target_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("YES", "NO"))
  
  ### run pca & umap
  target_Seurat_Obj <- RunPCA(target_Seurat_Obj,
                            features = VariableFeatures(object = target_Seurat_Obj),
                            npcs = 15)
  target_Seurat_Obj <- RunUMAP(target_Seurat_Obj, dims = 1:15)
  
  ### UMAP with all those info
  p <- list()
  p[[1]] <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
                    group.by = "px", shape.by = "New_Persistency",
                    pt.size = 3) +
    ggtitle("UMAP of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Patient",
         shape="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  
  ### umap for each px
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                              cells = rownames(target_Seurat_Obj@meta.data),
                              value = target_Seurat_Obj@meta.data$px)
  g <- ggplot_build(p[[1]])
  color_code <- data.frame(colours = unique(g$data[[1]]["colour"]), 
                           label = levels(g$plot$data[, "px"]),
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(color_code) <- color_code$label
  for(i in 1:length(unique(target_Seurat_Obj@meta.data$px))) {
    temp_seurat_obj <- subset(target_Seurat_Obj, idents = c(unique(target_Seurat_Obj@meta.data$px)[i]))
    p[[i+1]] <- DimPlot(object = temp_seurat_obj, reduction = "umap",
                        # cols = color_code[unique(target_Seurat_Obj@meta.data$px)[i],"colour"],
                        group.by = "New_Persistency",
                        pt.size = 3) +
      ggtitle(unique(target_Seurat_Obj@meta.data$px)[i]) +
      lims(x = g$layout$panel_scales_x[[1]]$range$range,
           y = g$layout$panel_scales_y[[1]]$range$range) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
      labs(color="Is Persistent")
    p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  ### arrange the plots and save
  fName <- paste0("UMAP_Plot_Persistence_per_Px_Sampled")
  rowNum <- 3
  colNum <- 4
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir, fName, "3.png"), g, width = 25, height = 15, dpi = 300)
  
  
  ### PCA with all those info
  p <- list()
  p[[1]] <- DimPlot(object = target_Seurat_Obj, reduction = "pca",
                    group.by = "px", shape.by = "New_Persistency",
                    pt.size = 3) +
    ggtitle("PCA of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Patient",
         shape="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  
  ### pca for each px
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  g <- ggplot_build(p[[1]])
  color_code <- data.frame(colours = unique(g$data[[1]]["colour"]), 
                           label = levels(g$plot$data[, "px"]),
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(color_code) <- color_code$label
  for(i in 1:length(unique(target_Seurat_Obj@meta.data$px))) {
    temp_seurat_obj <- subset(target_Seurat_Obj, idents = c(unique(target_Seurat_Obj@meta.data$px)[i]))
    p[[i+1]] <- DimPlot(object = temp_seurat_obj, reduction = "pca",
                        # cols = color_code[unique(target_Seurat_Obj@meta.data$px)[i],"colour"],
                        group.by = "New_Persistency",
                        pt.size = 3) +
      ggtitle(unique(target_Seurat_Obj@meta.data$px)[i]) +
      lims(x = g$layout$panel_scales_x[[1]]$range$range,
           y = g$layout$panel_scales_y[[1]]$range$range) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
      labs(color="Is Persistent")
    p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  ### arrange the plots and save
  fName <- paste0("PCA_Plot_Persistence_per_Px_Sampled")
  rowNum <- 3
  colNum <- 4
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir, fName, "3.png"), g, width = 25, height = 15, dpi = 300)
  
  
  #
  ### 2. with DE genes only
  #
  ### get DE genes
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$New_Persistency)
  de_result <- FindMarkers(target_Seurat_Obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### retain the DE genes only
  target_Seurat_Obj@assays$RNA@counts <- target_Seurat_Obj@assays$RNA@counts[rownames(de_result)[1:100],]
  target_Seurat_Obj@assays$RNA@data <- target_Seurat_Obj@assays$RNA@data[rownames(de_result)[1:100],]
  target_Seurat_Obj@assays$RNA@var.features <- rownames(de_result)[1:100]
  
  ### normalization
  target_Seurat_Obj <- NormalizeData(target_Seurat_Obj,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  
  ### find variable genes
  target_Seurat_Obj <- FindVariableFeatures(target_Seurat_Obj,
                                            selection.method = "vst", nfeatures = 2000)
  
  ### scaling
  target_Seurat_Obj <- ScaleData(target_Seurat_Obj,
                                 vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
  
  ### run pca & umap
  target_Seurat_Obj <- RunPCA(target_Seurat_Obj,
                              features = VariableFeatures(object = target_Seurat_Obj),
                              npcs = 15)
  target_Seurat_Obj <- RunUMAP(target_Seurat_Obj, dims = 1:15)
  
  ### UMAP with all those info
  p <- list()
  p[[1]] <- DimPlot(object = target_Seurat_Obj, reduction = "umap",
                    group.by = "px", shape.by = "New_Persistency",
                    pt.size = 3) +
    ggtitle("UMAP of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Patient",
         shape="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  
  ### umap for each px
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  g <- ggplot_build(p[[1]])
  color_code <- data.frame(colours = unique(g$data[[1]]["colour"]), 
                           label = levels(g$plot$data[, "px"]),
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(color_code) <- color_code$label
  for(i in 1:length(unique(target_Seurat_Obj@meta.data$px))) {
    temp_seurat_obj <- subset(target_Seurat_Obj, idents = c(unique(target_Seurat_Obj@meta.data$px)[i]))
    p[[i+1]] <- DimPlot(object = temp_seurat_obj, reduction = "umap",
                        # cols = color_code[unique(target_Seurat_Obj@meta.data$px)[i],"colour"],
                        group.by = "New_Persistency",
                        pt.size = 3) +
      ggtitle(unique(target_Seurat_Obj@meta.data$px)[i]) +
      lims(x = g$layout$panel_scales_x[[1]]$range$range,
           y = g$layout$panel_scales_y[[1]]$range$range) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
      labs(color="Is Persistent")
    p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  ### arrange the plots and save
  fName <- paste0("UMAP_Plot_Persistence_per_Px_Sampled_DE_Genes_Only")
  rowNum <- 3
  colNum <- 4
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir, fName, "3.png"), g, width = 25, height = 15, dpi = 300)
  
  ### PCA with all those info
  p <- list()
  p[[1]] <- DimPlot(object = target_Seurat_Obj, reduction = "pca",
                    group.by = "px", shape.by = "New_Persistency",
                    pt.size = 3) +
    ggtitle("PCA of SJCAR19 Data") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30)) +
    labs(color="Patient",
         shape="Is Persistent")
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  
  ### pca for each px
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  g <- ggplot_build(p[[1]])
  color_code <- data.frame(colours = unique(g$data[[1]]["colour"]), 
                           label = levels(g$plot$data[, "px"]),
                           stringsAsFactors = FALSE, check.names = FALSE)
  rownames(color_code) <- color_code$label
  for(i in 1:length(unique(target_Seurat_Obj@meta.data$px))) {
    temp_seurat_obj <- subset(target_Seurat_Obj, idents = c(unique(target_Seurat_Obj@meta.data$px)[i]))
    p[[i+1]] <- DimPlot(object = temp_seurat_obj, reduction = "pca",
                        # cols = color_code[unique(target_Seurat_Obj@meta.data$px)[i],"colour"],
                        group.by = "New_Persistency",
                        pt.size = 3) +
      ggtitle(unique(target_Seurat_Obj@meta.data$px)[i]) +
      lims(x = g$layout$panel_scales_x[[1]]$range$range,
           y = g$layout$panel_scales_y[[1]]$range$range) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 20)) +
      labs(color="Is Persistent")
    p[[i+1]]$layers[[1]]$aes_params$alpha <- 0.7
  }
  
  ### arrange the plots and save
  fName <- paste0("PCA_Plot_Persistence_per_Px_Sampled_DE_Genes_Only")
  rowNum <- 3
  colNum <- 4
  g <- arrangeGrob(grobs = p,
                   nrow = rowNum,
                   ncol = colNum,
                   top = textGrob(paste0(fName, "\n"), gp=gpar(fontsize=25)))
  ggsave(file = paste0(outputDir, fName, "3.png"), g, width = 25, height = 15, dpi = 300)
  
  
  #
  ### The next task
  #
  ### Only for the Px06 (Best Predictor) & Px08 (Worst Predictor),
  ### run GSEA in each patient with the signature as logFC of DE genes from PS vs NPS
  ### and compare the GSEA results with rank comparison
  
  ### the indicies of the persisters
  all_gmp_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES")
  all_gmp_not_last <- which(Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO")
  
  ### only use the CD8 cells
  all_gmp_last <- intersect(all_gmp_last,
                            which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  all_gmp_not_last <- intersect(all_gmp_not_last,
                                which(Seurat_Obj@meta.data$CD4_CD8_by_Consensus == "CD8"))
  
  ### only get the cells of interests
  target_Seurat_Obj <- subset(Seurat_Obj, cells = rownames(Seurat_Obj@meta.data)[c(all_gmp_last, all_gmp_not_last)])
  
  ### get the same number of samples (persister: non-persister - same ratio)
  set.seed(1234)
  target_Seurat_Obj@meta.data$New_Persistency <- NA
  multiplier_k <- 1
  for(px in unique(target_Seurat_Obj@meta.data$px)) {
    ### get specific indicies
    px_gmp_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                             which(target_Seurat_Obj@meta.data$GMP_CARpos_Persister == "YES"))
    px_gmp_not_last <- intersect(which(target_Seurat_Obj@meta.data$px == px),
                                 which(target_Seurat_Obj@meta.data$GMP_CARpos_Persister == "NO"))
    
    ### sampling
    if((length(px_gmp_last) > 2) && (length(px_gmp_not_last) > length(px_gmp_last))) {
      px_gmp_not_last <- sample(px_gmp_not_last, size = length(px_gmp_last)*multiplier_k)
    }
    
    ### annotate new persistency info
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_last] <- "YES"
    target_Seurat_Obj@meta.data$New_Persistency[px_gmp_not_last] <- "NO"
  }
  
  ### set idents with the new info
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$New_Persistency)
  
  ### remove NA and only contatin PS & NPS cells
  target_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("YES", "NO"))
  
  ### set idents with the px info
  target_Seurat_Obj <- SetIdent(object = target_Seurat_Obj,
                                cells = rownames(target_Seurat_Obj@meta.data),
                                value = target_Seurat_Obj@meta.data$px)
  
  ### Px06
  px06_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("SJCAR19-06"))
  
  ### get DE genes
  px06_Seurat_Obj <- SetIdent(object = px06_Seurat_Obj,
                              cells = rownames(px06_Seurat_Obj@meta.data),
                              value = px06_Seurat_Obj@meta.data$New_Persistency)
  de_result <- FindMarkers(px06_Seurat_Obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### GSEA directory
  output_dir <- paste0(outputDir, "GSEA/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  ### signature preparation
  signat <- de_result$avg_logFC
  names(signat) <- rownames(de_result)
  
  ### run GSEA
  GSEA_result <- run_gsea(gene_list = m_list, signature = list(signat), printPlot = FALSE)
  GSEA_result <- GSEA_result[order(GSEA_result$pval),]
  
  ### only get pathways that have pval < 0.01 & size > 30 & up-regulating results (enriched with important) only
  pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                      which(GSEA_result$size > 30)),
                                            which(GSEA_result$NES > 2))]
  if(length(pathways) < 5) {
    pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                        which(GSEA_result$size > 30)),
                                              which(GSEA_result$NES > 1.8))]
  }
  if(length(pathways) < 5) {
    pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                        which(GSEA_result$size > 30)),
                                              which(GSEA_result$NES > 1.6))]
  }
  if(length(pathways) > 30) {
    pathways <- GSEA_result$pathway[1:30]
  }
  
  ### run GSEA again with the significant result - plot printing
  dir.create(paste0(output_dir, "SJCAR19-06/"), showWarnings = FALSE, recursive = TRUE)
  GSEA_result2 <- run_gsea(gene_list = m_list[pathways], signature = list(signat),
                           heatmap_color_type = "absolute",
                           printPlot = TRUE, printPath = paste0(output_dir, "SJCAR19-06/"))
  
  ### write out the result file
  write.xlsx2(GSEA_result, file = paste0(output_dir, "/GSEA_Px06_Persisters_vs_Non-Persisters.xlsx"),
              sheetName = paste0("SJCAR19-06"), row.names = FALSE, append = TRUE)
  
  ### save the GSEA result
  px06_GSEA_result <- GSEA_result
  
  ### Px08
  px08_Seurat_Obj <- subset(target_Seurat_Obj, idents = c("SJCAR19-08"))
  
  ### get DE genes
  px08_Seurat_Obj <- SetIdent(object = px08_Seurat_Obj,
                              cells = rownames(px08_Seurat_Obj@meta.data),
                              value = px08_Seurat_Obj@meta.data$New_Persistency)
  de_result <- FindMarkers(px08_Seurat_Obj,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### GSEA directory
  output_dir <- paste0(outputDir, "GSEA/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  ### signature preparation
  signat <- de_result$avg_logFC
  names(signat) <- rownames(de_result)
  
  ### run GSEA
  GSEA_result <- run_gsea(gene_list = m_list, signature = list(signat), printPlot = FALSE)
  GSEA_result <- GSEA_result[order(GSEA_result$pval),]
  
  ### only get pathways that have pval < 0.01 & size > 30 & up-regulating results (enriched with important) only
  pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                      which(GSEA_result$size > 30)),
                                            which(GSEA_result$NES > 2))]
  if(length(pathways) < 5) {
    pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                        which(GSEA_result$size > 30)),
                                              which(GSEA_result$NES > 1.8))]
  }
  if(length(pathways) < 5) {
    pathways <- GSEA_result$pathway[intersect(intersect(which(GSEA_result$pval < 1e-02),
                                                        which(GSEA_result$size > 30)),
                                              which(GSEA_result$NES > 1.6))]
  }
  if(length(pathways) > 30) {
    pathways <- GSEA_result$pathway[1:30]
  }
  
  ### run GSEA again with the significant result - plot printing
  dir.create(paste0(output_dir, "SJCAR19-08/"), showWarnings = FALSE, recursive = TRUE)
  GSEA_result2 <- run_gsea(gene_list = m_list[pathways], signature = list(signat),
                           heatmap_color_type = "absolute",
                           printPlot = TRUE, printPath = paste0(output_dir, "SJCAR19-08/"))
  
  ### write out the result file
  write.xlsx2(GSEA_result, file = paste0(output_dir, "/GSEA_Px08_Persisters_vs_Non-Persisters.xlsx"),
              sheetName = paste0("SJCAR19-08"), row.names = FALSE, append = TRUE)
  
  ### save the GSEA result
  px08_GSEA_result <- GSEA_result
  
  #
  ### GSEA comparison between Px06 & Px08
  #
  ### give ranks
  px06_GSEA_result$rank <- rank(px06_GSEA_result$padj)
  px08_GSEA_result$rank <- rank(px08_GSEA_result$padj)
  
  ### set rownames
  rownames(px06_GSEA_result) <- px06_GSEA_result$pathway
  rownames(px08_GSEA_result) <- px08_GSEA_result$pathway
  
  ### get common pathways
  common_pathways <- intersect(px06_GSEA_result$pathway[which(px06_GSEA_result$padj < 1)],
                               px08_GSEA_result$pathway[which(px08_GSEA_result$padj < 1)])
  
  ### get signatures - RANK
  sig1 <- px06_GSEA_result[common_pathways,"rank"]
  names(sig1) <- common_pathways
  sig2 <- px08_GSEA_result[common_pathways,"rank"]
  names(sig2) <- common_pathways
  
  ### compare the ranks
  if(length(sig1) > 2 && length(sig2) > 2) {
    compare_two_different_ranks(A = sig1,
                                B = sig2,
                                A_name = paste0("GSEA Pathway Rank of Px06"),
                                B_name = paste0("GSEA Pathway Rank of Px08"),
                                ordering = "decreasing",
                                top = 300,
                                alternative = "greater",
                                permutation = 10000,
                                fileName = paste0("GSEA_Comparison_Px06_vs_Px08"),
                                printPath = output_dir,
                                width = 24,
                                height = 12)
  } else {
    writeLines("There aren't enough pathways for the analysis.")
  }
  
  
  
  
}
