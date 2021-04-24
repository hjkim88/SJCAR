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
#                  c) Cluster 0,1,4,8 vs others
#                  d) Pseudotime analysis on PCA
#                  e) PCA/UMAP plot of lineages with size=1 vs lineages with size > 1 (coloring differently)
#               6. Visualize some interesting genes on UMAP
#               7. Time series DE analysis & pathway analysis
#
#   Instruction
#               1. Source("Manuscript_Figures_And_Tables.R")
#               2. Run the function "manuscript_prep" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Manuscript_Figures_And_Tables.R/Manuscript_Figures_And_Tables.R")
#               > manuscript_prep(Seurat_RObj_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/data/SJCAR19_Oct2020_Seurat_Obj.RDS",
#                                 gmp_carpos_cd8_de_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/GMP_CARpos_CD8_Persisters_vs_NonPersisters.xlsx",
#                                 px_result_dir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/",
#                                 outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin/JCC212_SJCAR19/new_results/Manuscript/")
###

manuscript_prep <- function(Seurat_RObj_path="./data/NEW_SJCAR_SEURAT_OBJ/SJCAR19_Oct2020_Seurat_Obj_Total2.RDS",
                            gmp_carpos_cd8_de_path="./results/New3/GMP_CARpos_CD8_Persisters_vs_NonPersisters.xlsx",
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
  if(!require(slingshot, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("slingshot")
    require(slingshot, quietly = TRUE)
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
  
  #' @title Plot Slingshot output
  #' @name plot-SlingshotDataSet
  #' @aliases plot-SlingshotDataSet plot,SlingshotDataSet,ANY-method
  #'
  #' @description Tools for visualizing lineages inferred by \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{"curves"}, or \code{"both"} (by partial matching),
  #'   see Details for more.
  #' @param linInd integer, an index indicating which lineages should be plotted
  #'   (default is to plot all lineages). If \code{col} is a vector, it will be
  #'   subsetted by \code{linInd}.
  #' @param show.constraints logical, whether or not the user-specified initial
  #'   and terminal clusters should be specially denoted by green and red dots,
  #'   respectively.
  #' @param add logical, indicates whether the output should be added to an
  #'   existing plot.
  #' @param dims numeric, which dimensions to plot (default is \code{1:2}).
  #' @param asp numeric, the y/x aspect ratio, see \code{\link{plot.window}}.
  #' @param cex numeric, amount by which points should be magnified, see
  #'   \code{\link{par}}.
  #' @param lwd numeric, the line width, see \code{\link{par}}.
  #' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
  #' @param ... additional parameters to be passed to \code{\link{lines}}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' plot(sds, type = 'b')
  #'
  #' # add to existing plot
  #' plot(rd, col = 'grey50')
  #' lines(sds, lwd = 3)
  #'
  #' @import graphics
  #' @import grDevices
  #' @export
  setMethod(
    f = "plot",
    signature = signature(x = "SlingshotDataSet"),
    definition = function(x, type = NULL,
                          linInd = NULL,
                          show.constraints = FALSE,
                          constraints.col = NULL,
                          add = FALSE,
                          dims = seq_len(2),
                          asp = 1,
                          cex = 2,
                          lwd = 2,
                          col = 1,
                          ...) {
      col <- rep(col, length(slingLineages(x)))
      curves <- FALSE
      lineages <- FALSE
      if(is.null(type)){
        if(length(slingCurves(x)) > 0){
          type <- 'curves'
        }else if(length(slingLineages(x)) > 0){
          type <- 'lineages'
        }else{
          stop('No lineages or curves detected.')
        }
      }else{
        type <- c('curves','lineages','both')[pmatch(type,
                                                     c('curves','lineages','both'))]
        if(is.na(type)){
          stop('Unrecognized type argument.')
        }
      }
      
      if(type %in% c('lineages','both')){
        lineages <- TRUE
      }
      if(type %in% c('curves','both')){
        curves <- TRUE
      }
      
      if(lineages & (length(slingLineages(x))==0)){
        stop('No lineages detected.')
      }
      if(curves & (length(slingCurves(x))==0)){
        stop('No curves detected.')
      }
      
      if(is.null(linInd)){
        linInd <- seq_along(slingLineages(x))
      }else{
        linInd <- as.integer(linInd)
        if(!all(linInd %in% seq_along(slingLineages(x)))){
          if(any(linInd %in% seq_along(slingLineages(x)))){
            linInd.removed <-
              linInd[! linInd %in% seq_along(slingLineages(x))]
            linInd <-
              linInd[linInd %in% seq_along(slingLineages(x))]
            message('Unrecognized lineage indices (linInd): ',
                    paste(linInd.removed, collapse = ", "))
          }else{
            stop('None of the provided lineage indices',
                 ' (linInd) were found.')
          }
        }
      }
      
      if(lineages){
        X <- reducedDim(x)
        clusterLabels <- slingClusterLabels(x)
        connectivity <- slingAdjacency(x)
        clusters <- rownames(connectivity)
        nclus <- nrow(connectivity)
        centers <- t(vapply(clusters,function(clID){
          w <- clusterLabels[,clID]
          return(apply(X, 2, weighted.mean, w = w))
        }, rep(0,ncol(X))))
        rownames(centers) <- clusters
        X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
        clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                       drop = FALSE]
        linC <- slingParams(x)
        clus2include <- unique(unlist(slingLineages(x)[linInd]))
      }
      
      if(!add){
        xs <- NULL
        ys <- NULL
        if(lineages){
          xs <- c(xs, centers[,dims[1]])
          ys <- c(ys, centers[,dims[2]])
        }
        if(curves){
          npoints <- nrow(slingCurves(x)[[1]]$s)
          xs <- c(xs, as.numeric(vapply(slingCurves(x),
                                        function(c){ c$s[,dims[1]] }, rep(0,npoints))))
          ys <- c(ys, as.numeric(vapply(slingCurves(x),
                                        function(c){ c$s[,dims[2]] }, rep(0,npoints))))
        }
        plot(x = NULL, y = NULL, asp = asp,
             xlim = range(xs), ylim = range(ys),
             xlab = colnames(reducedDim(x))[dims[1]],
             ylab = colnames(reducedDim(x))[dims[2]])
      }
      
      if(lineages){
        for(i in seq_len(nclus-1)){
          for(j in seq(i+1,nclus)){
            if(connectivity[i,j]==1 &
               all(clusters[c(i,j)] %in% clus2include)){
              lines(centers[c(i,j), dims],
                    lwd = lwd, col = col[1], ...)
            }
          }
        }
        points(centers[clusters %in% clus2include, dims],
               cex = cex+1, pch = 16, col = col[1])
        if(show.constraints && !is.null(constraints.col)){
          for(const in names(constraints.col)) {
            points(centers[clusters %in% const, dims,
                           drop=FALSE], cex = cex-0.5,
                   col = constraints.col[const], pch = 16)
            text(x = centers[clusters %in% const, dims[1]]+0,
                 y = centers[clusters %in% const, dims[2]]+1.3,
                 labels = const,
                 font = 2,
                 cex = cex-0.5,
                 col = "black")
          }
        }
      }
      if(curves){
        for(ii in seq_along(slingCurves(x))[linInd]){
          c <- slingCurves(x)[[ii]]
          lines(c$s[c$ord, dims], lwd = lwd, col = col[ii], ...)
        }
      }
      invisible(NULL)
    }
  )
  
  #' @title Pairs plot of Slingshot output
  #' @name pairs-SlingshotDataSet
  #'
  #' @description A tool for quickly visualizing lineages inferred by
  #'   \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
  #'   Details for more.
  #' @param show.constraints logical, whether or not the user-specified initial
  #'   and terminal clusters should be specially denoted by green and red dots,
  #'   respectively.
  #' @param col character, color vector for points.
  #' @param pch integer or character specifying the plotting symbol, see
  #'   \code{\link{par}}.
  #' @param cex numeric, amount by which points should be magnified, see
  #'   \code{\link{par}}.
  #' @param lwd numeric, the line width, see \code{\link{par}}.
  #' @param ... additional parameters for \code{plot} or \code{axis}, see
  #'   \code{\link{pairs}}.
  #' @param labels character, the names of the variables, see \code{\link{pairs}}.
  #' @param horInd see \code{\link{pairs}}.
  #' @param verInd see \code{\link{pairs}}.
  #' @param lower.panel see \code{\link{pairs}}.
  #' @param upper.panel see \code{\link{pairs}}.
  #' @param diag.panel see \code{\link{pairs}}.
  #' @param text.panel see \code{\link{pairs}}.
  #' @param label.pos see \code{\link{pairs}}.
  #' @param line.main see \code{\link{pairs}}.
  #' @param cex.labels see \code{\link{pairs}}.
  #' @param font.labels see \code{\link{pairs}}.
  #' @param row1attop see \code{\link{pairs}}.
  #' @param gap see \code{\link{pairs}}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' pairs(sds, type = 'curves')
  #'
  #' @export
  pairs.SlingshotDataSet <-
    function (x, type = NULL, show.constraints = FALSE, col = NULL,
              constraints.col = NULL,
              pch = 16, cex=1, lwd=2, ...,
              labels, horInd = seq_len(nc), verInd = seq_len(nc),
              lower.panel = FALSE, upper.panel = TRUE,
              diag.panel = NULL, text.panel = textPanel,
              label.pos = 0.5 + has.diag/3, line.main = 3,
              cex.labels = NULL, font.labels = 1,
              row1attop = TRUE, gap = 1,
              xlim=NULL, ylim=NULL) {
      #####
      lp.sling <- lower.panel
      up.sling <- upper.panel
      panel <- points
      if(!up.sling){
        upper.panel <- NULL
      }else{
        upper.panel <- panel
      }
      if(!lower.panel){
        lower.panel <- NULL
      }else{
        lower.panel <- panel
      }
      log = ""
      sds <- x
      x <- reducedDim(sds)
      curves <- FALSE
      lineages <- FALSE
      if(is.null(type)){
        if(length(slingCurves(sds)) > 0){
          type <- 'curves'
        }else if(length(slingLineages(sds)) > 0){
          type <- 'lineages'
        }else{
          stop('No lineages or curves detected.')
        }
      }else{
        type <- c('curves','lineages','both')[pmatch(type,
                                                     c('curves','lineages',
                                                       'both'))]
        if(is.na(type)){
          stop('Unrecognized type argument.')
        }
      }
      if(type %in% c('lineages','both')){
        lineages <- TRUE
      }
      if(type %in% c('curves','both')){
        curves <- TRUE
      }
      if(lineages & (length(slingLineages(sds))==0)){
        stop('No lineages detected.')
      }
      if(curves & (length(slingCurves(sds))==0)){
        stop('No curves detected.')
      }
      if(lineages){
        forest <- slingAdjacency(sds)
        clusters <- rownames(forest)
        nclus <- nrow(forest)
        centers <- t(vapply(clusters,function(clID){
          w <- slingClusterLabels(sds)[,clID]
          return(apply(x, 2, weighted.mean, w = w))
        }, rep(0,ncol(reducedDim(sds)))))
        rownames(centers) <- clusters
        linC <- slingParams(sds)
      }
      range.max <- max(apply(x,2,function(xi){
        r <- range(xi, na.rm = TRUE)
        return(abs(r[2] - r[1]))
      }))
      plot.ranges <- apply(x,2,function(xi){
        mid <- (max(xi,na.rm = TRUE) + min(xi,na.rm = TRUE))/2
        return(c(mid - range.max/2, mid + range.max/2))
      })
      if(is.null(col)){
        if(requireNamespace("RColorBrewer", quietly = TRUE)) {
          cc <- c(RColorBrewer::brewer.pal(9, "Set1")[-c(1,3,6)],
                  RColorBrewer::brewer.pal(7, "Set2")[-2],
                  RColorBrewer::brewer.pal(6, "Dark2")[-5],
                  RColorBrewer::brewer.pal(8, "Set3")[-c(1,2)])
        } else {
          cc <- seq_len(100)
        }
        col <- cc[apply(slingClusterLabels(sds),1,which.max)]
      }
      #####
      if(doText <- missing(text.panel) || is.function(text.panel))
        textPanel <-
        function(x = 0.5, y = 0.5, txt, cex, font)
          text(x, y, txt, cex = cex, font = font)
      
      localAxis <- function(side, x, y, xpd, bg, col=NULL, lwd=NULL, main,
                            oma, ...) {
        ## Explicitly ignore any color argument passed in as
        ## it was most likely meant for the data points and
        ## not for the axis.
        xpd <- NA
        if(side %% 2L == 1L && xl[j]) xpd <- FALSE
        if(side %% 2L == 0L && yl[i]) xpd <- FALSE
        if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
      }
      
      localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
      localLowerPanel <- function(..., main, oma, font.main, cex.main)
        lower.panel(...)
      localUpperPanel <- function(..., main, oma, font.main, cex.main)
        upper.panel(...)
      localDiagPanel <- function(..., main, oma, font.main, cex.main)
        diag.panel(...)
      
      dots <- list(...); nmdots <- names(dots)
      if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for(i in seq_along(names(x))) {
          if(is.factor(x[[i]]) || is.logical(x[[i]]))
            x[[i]] <- as.numeric(x[[i]])
          if(!is.numeric(unclass(x[[i]])))
            stop("non-numeric argument to 'pairs'")
        }
      } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")
      panel <- match.fun(panel)
      if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel <- match.fun(lower.panel)
      if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel <- match.fun(upper.panel)
      if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))
        diag.panel <- match.fun( diag.panel)
      
      if(row1attop) {
        tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
        tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
      }
      
      nc <- ncol(x)
      if (nc < 2L) stop("only one column in the argument to 'pairs'")
      if(!all(horInd >= 1L & horInd <= nc))
        stop("invalid argument 'horInd'")
      if(!all(verInd >= 1L & verInd <= nc))
        stop("invalid argument 'verInd'")
      if(doText) {
        if (missing(labels)) {
          labels <- colnames(x)
          if (is.null(labels)) labels <- paste("var", 1L:nc)
        }
        else if(is.null(labels)) doText <- FALSE
      }
      oma <- if("oma" %in% nmdots) dots$oma
      main <- if("main" %in% nmdots) dots$main
      if (is.null(oma))
        oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)
      opar <- par(mfrow = c(length(horInd), length(verInd)),
                  mar = rep.int(gap/2, 4), oma = oma)
      on.exit(par(opar))
      dev.hold(); on.exit(dev.flush(), add = TRUE)
      
      xl <- yl <- logical(nc)
      if (is.numeric(log)) xl[log] <- yl[log] <- TRUE
      else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}
      for (i in if(row1attop) verInd else rev(verInd))
        for (j in horInd) {
          l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
          
          if(is.null(xlim) & !is.null(ylim))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = plot.ranges[,j], ylim=ylim)
          else if(!is.null(xlim) & is.null(ylim))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = xlim, ylim = plot.ranges[,i])
          else if(!is.null(xlim) & !is.null(ylim))
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = xlim, ylim=ylim)
          else
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ..., log = l,
                      xlim = plot.ranges[,j], ylim = plot.ranges[,i])
          
          if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
            box()
            if(i == 1  && (!(j %% 2L) || !has.upper || !has.lower ))
              localAxis(1L + 2L*row1attop, x[, j], x[, i], ...)
            if(i == nc && (  j %% 2L  || !has.upper || !has.lower ))
              localAxis(3L - 2L*row1attop, x[, j], x[, i], ...)
            if(j == 1  && (!(i %% 2L) || !has.upper || !has.lower ))
              localAxis(2L, x[, j], x[, i], ...)
            if(j == nc && (  i %% 2L  || !has.upper || !has.lower ))
              localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if(i == j) {
              if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
              if (doText) {
                par(usr = c(0, 1, 0, 1))
                if(is.null(cex.labels)) {
                  l.wid <- strwidth(labels, "user")
                  cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                }
                xlp <- if(xl[i]) 10^0.5 else 0.5
                ylp <- if(yl[j]) 10^label.pos else label.pos
                text.panel(xlp, ylp, labels[i],
                           cex = cex.labels, font = font.labels)
              }
            } else if(i < j){
              if(up.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        seg.col <- 1
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col, ...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex=2*cex)
                  if(show.constraints && is.null(constraints.col)){
                    if(any(linC$start.given)){
                      st.ind <- clusters %in%
                        linC$start.clus[linC$start.given]
                      points(centers[st.ind,j],
                             centers[st.ind,i], cex = cex,
                             col = 'green3',
                             pch = pch)
                    }
                    if(any(linC$end.given)){
                      en.ind <- clusters %in%
                        linC$end.clus[linC$end.given]
                      points(centers[en.ind,j],
                             centers[en.ind,i], cex = cex,
                             col = 'red2', pch = pch)
                    }
                  } else if(show.constraints && !is.null(constraints.col)){
                    for(const in names(constraints.col)) {
                      points(centers[clusters %in% const, j, drop=FALSE],
                             centers[clusters %in% const, i, drop=FALSE],
                             cex = cex, pch = 16,
                             col = constraints.col[const])
                    }
                  }
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)], lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            else{
              if(lp.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        if(clusters[ii] %in%
                           linC$start.clus |
                           clusters[jj] %in%
                           linC$start.clus){
                          seg.col <- 'green3'
                        }else if(clusters[ii] %in%
                                 linC$end.clus[
                                   linC$end.given] |
                                 clusters[jj] %in%
                                 linC$end.clus[
                                   linC$end.given]){
                          seg.col <- 'red2'
                        }else{
                          seg.col <- 1
                        }
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col,...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex = 2*cex)
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)],lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            if (any(par("mfg") != mfg))
              stop("the 'panel' function made a new plot")
          } else par(new = FALSE)
          
        }
      if (!is.null(main)) {
        font.main <- if("font.main" %in% nmdots){
          dots$font.main
        }else par("font.main")
        cex.main <- if("cex.main" %in%
                       nmdots) dots$cex.main else par("cex.main")
        mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main,
              font = font.main)
      }
      invisible(NULL)
    }
  
  #' @name plot3d-SlingshotDataSet
  #' @title Plot Slingshot output in 3D
  #'
  #' @description Tools for visualizing lineages inferred by \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
  #'   Details for more.
  #' @param linInd integer, an index indicating which lineages should be plotted
  #'   (default is to plot all lineages). If \code{col} is a vector, it will be
  #'   subsetted by \code{linInd}.
  #' @param add logical, indicates whether the output should be added to an
  #'   existing plot.
  #' @param dims numeric, which dimensions to plot (default is \code{1:3}).
  #' @param aspect either a logical indicating whether to adjust the aspect ratio
  #'   or a new ratio, see \code{\link[rgl:plot3d]{plot3d}}.
  #' @param size numeric, size of points for MST (default is \code{10}), see
  #'   \code{\link[rgl:plot3d]{plot3d}}.
  #' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
  #' @param ... additional parameters to be passed to \code{lines3d}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' \dontrun{
  #' library(rgl)
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' rd <- cbind(rd, rnorm(nrow(rd)))
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' plot3d(sds, type = 'b')
  #'
  #' # add to existing plot
  #' plot3d(rd, col = 'grey50', aspect = 'iso')
  #' plot3d(sds, lwd = 3, add = TRUE)
  #' }
  # #' @importFrom rgl plot3d
  #' @export
  plot3d.SlingshotDataSet <- function(x,
                                      type = NULL,
                                      linInd = NULL,
                                      add = FALSE,
                                      dims = seq_len(3),
                                      aspect = 'iso',
                                      size = 10,
                                      col = 1,
                                      col2 = NULL,
                                      ...){
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("Package 'rgl' is required for 3D plotting.",
           call. = FALSE)
    }
    col <- rep(col, length(slingLineages(x)))
    n <- nrow(reducedDim(x))
    curves <- FALSE
    lineages <- FALSE
    if(is.null(type)){
      if(length(slingCurves(x)) > 0){
        type <- 'curves'
      }else if(length(slingLineages(x)) > 0){
        type <- 'lineages'
      }else{
        stop('No lineages or curves detected.')
      }
    }else{
      type <- c('curves','lineages','both')[pmatch(type,c('curves','lineages',
                                                          'both'))]
      if(is.na(type)){
        stop('Unrecognized type argument.')
      }
    }
    
    if(type %in% c('lineages','both')){
      lineages <- TRUE
    }
    if(type %in% c('curves','both')){
      curves <- TRUE
    }
    
    if(lineages & (length(slingLineages(x))==0)){
      stop('No lineages detected.')
    }
    if(curves & (length(slingCurves(x))==0)){
      stop('No curves detected.')
    }
    
    if(is.null(linInd)){
      linInd <- seq_along(slingLineages(x))
    }else{
      linInd <- as.integer(linInd)
      if(!all(linInd %in% seq_along(slingLineages(x)))){
        if(any(linInd %in% seq_along(slingLineages(x)))){
          linInd.removed <-
            linInd[! linInd %in% seq_along(slingLineages(x))]
          linInd <-
            linInd[linInd %in% seq_along(slingLineages(x))]
          message('Unrecognized lineage indices (linInd): ',
                  paste(linInd.removed, collapse = ", "))
        }else{
          stop('None of the provided lineage indices',
               ' (linInd) were found.')
        }
      }
    }
    
    if(lineages){
      X <- reducedDim(x)
      clusterLabels <- slingClusterLabels(x)
      connectivity <- slingAdjacency(x)
      clusters <- rownames(connectivity)
      nclus <- nrow(connectivity)
      centers <- t(vapply(clusters,function(clID){
        w <- clusterLabels[,clID]
        return(apply(X, 2, weighted.mean, w = w))
      }, rep(0,ncol(X))))
      rownames(centers) <- clusters
      X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
      clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                     drop = FALSE]
      clus2include <- unique(unlist(slingLineages(x)[linInd]))
    }
    
    if(!add){
      xs <- NULL
      ys <- NULL
      zs <- NULL
      if(lineages){
        xs <- c(xs, centers[,dims[1]])
        ys <- c(ys, centers[,dims[2]])
        zs <- c(zs, centers[,dims[3]])
      }
      if(curves){
        npoints <- nrow(slingCurves(x)[[1]]$s)
        xs <- c(xs, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[1]] }, rep(0,npoints))))
        ys <- c(ys, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[2]] }, rep(0,npoints))))
        zs <- c(zs, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[3]] }, rep(0,npoints))))
      }
      rgl::plot3d(x = NULL, y = NULL, z = NULL, aspect = aspect,
                  xlim = range(xs), ylim = range(ys), zlim = range(zs),
                  xlab = colnames(reducedDim(x))[dims[1]],
                  ylab = colnames(reducedDim(x))[dims[2]],
                  zlab = colnames(reducedDim(x))[dims[3]])
    }
    
    if(lineages){
      for(i in seq_len(nclus-1)){
        for(j in seq(i+1,nclus)){
          if(connectivity[i,j]==1 &
             all(clusters[c(i,j)] %in% clus2include)){
            rgl::lines3d(x = centers[c(i,j),dims[1]],
                         y = centers[c(i,j),dims[2]],
                         z = centers[c(i,j),dims[3]],
                         col = col[1], ...)
          }
        }
      }
      rgl::points3d(centers[clusters %in% clus2include, dims],
                    size = size/2, col = col2[clusters[clusters %in% clus2include]])
      rgl::points3d(centers[clusters %in% clus2include, dims],
                    size = size, col = col[1])
    }
    if(curves){
      for(ii in seq_along(slingCurves(x))[linInd]){
        c <- slingCurves(x)[[ii]]
        rgl::lines3d(c$s[c$ord,dims], col = col[ii], ...)
      }
    }
    invisible(NULL)
  }
  
  ### the plot3d.SlingshotDataSet of the slingshot package is incomplete and too simple,
  ### so, i'm implementing a 3d plot function myself
  slingshot_3d_lineages <- function(slingshot_obj, color, title,
                                    print=FALSE, outputDir=NULL,
                                    width=1200, height=800) {
    
    ### load libraries
    if(!require(Seurat, quietly = TRUE)) {
      install.packages("Seurat")
      require(Seurat, quietly = TRUE)
    }
    if(!require(slingshot, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("slingshot")
      require(slingshot, quietly = TRUE)
    }
    if(!require(rgl, quietly = TRUE)) {
      install.packages("rgl")
      require(rgl, quietly = TRUE)
    }
    if(!require(rmarkdown, quietly = TRUE)) {
      install.packages("rmarkdown")
      require(rmarkdown, quietly = TRUE)
    }
    
    #
    ### 3D Slingshot
    #
    
    ### draw 3D PCA
    par3d(windowRect = c(50, 50, width+50, height+50))
    plot3d.SlingshotDataSet(slingshot_obj, dims = 1:3, col = "black", col2 = color, type = "lineages", add = TRUE)
    plot3d(slingshot_obj@reducedDim, col = apply(slingshot_obj@clusterLabels, 1, function(x) color[names(x)[which(x == 1)]]),
           size = 5, alpha = 0.5, aspect = FALSE, add = TRUE)
    axes3d(edges=c("x+-", "y+-", "z++"), lwd = 2,
           labels=TRUE, tick = FALSE, nticks = 3, box = TRUE, expand = 1.05)
    mtext3d(text = expression(bold("PC1")), edge="x+-", line = -2, at = min(slingshot_obj@reducedDim[,1]), pos = NA)
    mtext3d(text = expression(bold("PC2")), edge="y+-", line = -2, at = min(slingshot_obj@reducedDim[,2]), pos = NA)
    mtext3d(text = expression(bold("PC3")), edge="z++", line = -2, at = max(slingshot_obj@reducedDim[,3]), pos = NA)
    decorate3d(xlim = NULL, ylim = NULL, zlim = NULL, 
               xlab = "", ylab = "", zlab = "", 
               box = FALSE, axes = FALSE, main = title, sub = NULL,
               top = TRUE, aspect = FALSE, expand = 1.05)
    legend3d("topright", legend = names(color), title = "Clusters",
             col = color, pch = 19, cex=2)
    if(print) {
      writeWebGL(dir=outputDir, filename = paste0(outputDir, title, ".html"),
                 width=width, height = height)
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
  
  #
  ### 5. b) Clustering + UMAP
  #
  
  ### perform clustering
  sub_seurat_obj2 <- FindNeighbors(sub_seurat_obj2, dims = 1:10)
  sub_seurat_obj2 <- FindClusters(sub_seurat_obj2, resolution = 0.5)
  
  ### save the clustering result to meta.data
  sub_seurat_obj2@meta.data$clusters <- Idents(sub_seurat_obj2)
  
  ### UMAP with clusters by each patient
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "clusters", split.by = "px",
               pt.size = 3, ncol = 3) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 36) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 30),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 30))
  ggsave(paste0(outputDir2, "UMAP_CARpos_Clusters.png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### UMAP with clusters with all patients
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "clusters",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Clusters") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "UMAP_CARpos_Clusters_ALL.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  #
  ### 5. c) Cluster 0,1,4,8 vs others
  #
  
  ### set column for cluster 0,1,4,8 and the others
  sub_seurat_obj2$Cluster0148 <- "NO"
  sub_seurat_obj2$Cluster0148[which(as.character(sub_seurat_obj2@meta.data$clusters) %in% c("0", "1", "4", "8"))] <- "YES"
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$Cluster0148)
  
  ### DE analysis
  de_result <- FindMarkers(sub_seurat_obj2,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.2,
                           logfc.threshold = 0.2,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_Cluster0148_vs_Others.xlsx"),
              sheetName = "CARpos_Cluster0148_DE_Result", row.names = FALSE)
  
  
  #
  ### 6. Visualize some interesting genes on UMAP
  #
  
  ### create outputDir
  outputDir2 <- paste0(outputDir, "/6/")
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
  
  ### UMAP with px info
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "px",
               pt.size = 5) +
    ggtitle("") +
    labs(color="Patient") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.7
  ggsave(paste0(outputDir2, "UMAP_CARpos_Px_ALL.png"), plot = p, width = 15, height = 10, dpi = 400)
  
  ### lineages
  persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "YES")])
  non_persister_clones <- unique(Seurat_Obj@meta.data$clonotype_id_by_patient_one_alpha_beta[which(Seurat_Obj@meta.data$GMP_CARpos_CD8_Persister == "NO")])
  
  ### define persisters
  sub_seurat_obj2@meta.data$CD8_Persisters <- "Non-subsisters"
  sub_seurat_obj2@meta.data$CD8_Persisters[which(sub_seurat_obj2@meta.data$clonotype_id_by_patient_one_alpha_beta %in% persister_clones)] <- "Subsisters"
  
  ### UMAP with subsister info
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "CD8_Persisters",
               pt.size = 5, cols = c("lightgray", "red")) +
    ggtitle("") +
    labs(color="Patient") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.3
  ggsave(paste0(outputDir2, "UMAP_CARpos_Subsister_ALL.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### UMAP with subsister info split by each patient
  p <- DimPlot(object = sub_seurat_obj2, reduction = "umap",
               group.by = "CD8_Persisters", split.by = "px",
               pt.size = 5, cols = c("lightgray", "red"), ncol = 3) +
    ggtitle("") +
    labs(color="Patient") +
    theme_classic(base_size = 64) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48),
          axis.text.x = element_text(size = 48),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 48),
          legend.title = element_text(size = 36),
          legend.text = element_text(size = 36))
  p[[1]]$layers[[1]]$aes_params$alpha <- 0.3
  ggsave(paste0(outputDir2, "UMAP_CARpos_Subsister.png"), plot = p, width = 30, height = 20, dpi = 350)
  
  ### load DE genes
  gmp_carpos_cd8_de_genes <- read.xlsx2(gmp_carpos_cd8_de_path, sheetIndex = 1,
                                        stringsAsFactors = FALSE, check.names = FALSE)
  rownames(gmp_carpos_cd8_de_genes) <- gmp_carpos_cd8_de_genes$Gene
  
  ### See gene expressions on UMAP with the top 1:9 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[1:9], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(1).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 10:18 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[10:18], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(2).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 19:27 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[19:27], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(3).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 28:36 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[28:36], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(4).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 37:45 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = gmp_carpos_cd8_de_genes$Gene[37:45], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_GEXP_9_DE_Genes(5).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### set column for cluster 0,1,4,8 and the others
  sub_seurat_obj2$Cluster0148 <- "NO"
  sub_seurat_obj2$Cluster0148[which(as.character(sub_seurat_obj2@meta.data$clusters) %in% c("0", "1", "4", "8"))] <- "YES"
  sub_seurat_obj2 <- SetIdent(object = sub_seurat_obj2,
                              cells = rownames(sub_seurat_obj2@meta.data),
                              value = sub_seurat_obj2@meta.data$Cluster0148)
  
  ### DE analysis
  de_result <- FindMarkers(sub_seurat_obj2,
                           ident.1 = "YES",
                           ident.2 = "NO",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### See gene expressions on UMAP with the top 1:9 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[1:9], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0148_GEXP_9_DE_Genes(1).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 10:18 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[10:18], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0148_GEXP_9_DE_Genes(2).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 19:27 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[19:27], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0148_GEXP_9_DE_Genes(3).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 28:36 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[28:36], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0148_GEXP_9_DE_Genes(4).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### See gene expressions on UMAP with the top 37:45 DE genes
  p <- FeaturePlot(sub_seurat_obj2, features = rownames(de_result)[37:45], cols = c("lightgray", "red"))
  ggsave(paste0(outputDir2, "CARpos_Cluster0148_GEXP_9_DE_Genes(5).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  # CD3D - t cells
  # CD4 - cd4 t cells
  # CD8A - cd8 t cells
  # KLRF1 - NK cells
  # CD14 - monocytes/macrophages
  # FCGR3A - NK cells
  # FOXP3 - t regs
  # STAT4 - th1 cells
  # STAT6 - th2 cells
  # STAT3 - th17 cells
  # CD44 - memory
  # SELL(CD62L) - central memory cells
  # IL7R - central memory cells
  # CCR7 - central memory
  # CD27 - effector cells
  # KLRG1 - SLEC (short-lived effector cells)
  # CD69 - tissue resident memory
  # CD127 - MPEC (memory precursor effector cells)
  # CD45RO - effector memory t cells
  
  ### 1
  # CD4 - cd4 t cells
  # CD8A - cd8 t cells
  # KLRF1 - NK cells
  # CD14 - monocytes/macrophages
  # FCGR3A - NK cells
  # FOXP3 - t regs
  # STAT4 - th1 cells
  # STAT6 - th2 cells
  # STAT3 - th17 cells
  gene_set <- c("CD4", "CD8A", "KLRF1",
                "CD14", "FCGR3A", "FOXP3",
                "STAT4", "STAT6", "STAT3")
  names(gene_set) <- c("CD4", "CD8", "NK",
                       "Mono/Mcrphg", "NK", "Treg",
                       "Th1", "Th2", "Th17")
  p <- FeaturePlot(sub_seurat_obj2, features = gene_set,
              cols = c("lightgray", "red"))
  for(i in 1:length(gene_set)) {
    p[[i]]$labels$title <- paste(gene_set[i], "-", names(gene_set)[i])
  }
  ggsave(paste0(outputDir2, "CARpos_9_Interesting_Genes(1).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### 2
  # CD4 - cd4 t cells
  # CD8A - cd8 t cells
  # CD44 - memory
  # SELL(CD62L) - central memory cells
  # IL7R - central memory cells
  # CCR7 - central memory
  # CD27 - effector cells
  # KLRG1 - SLEC (short-lived effector cells)
  # CD69 - tissue resident memory
  gene_set <- c("CD4", "CD8A", "CD44",
                "SELL", "IL7R", "CCR7",
                "CD27", "KLRG1", "CD69")
  names(gene_set) <- c("CD4", "CD8", "Memory",
                       "(CD62L) Central Memory", "Centeral Memory", "Central Memory",
                       "Effector", "SLEC", "Tissue Resident Memory")
  p <- FeaturePlot(sub_seurat_obj2, features = gene_set,
                   cols = c("lightgray", "red"), max.cutoff = 3)
  for(i in 1:length(gene_set)) {
    p[[i]]$labels$title <- paste(gene_set[i], "-", names(gene_set)[i])
  }
  ggsave(paste0(outputDir2, "CARpos_9_Interesting_Genes(2).png"), plot = p, width = 15, height = 10, dpi = 350)
  
  
}
