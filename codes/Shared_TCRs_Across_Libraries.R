###
#   File name : Shared_TCRs_Across_Libraries.R
#   Author    : Hyunjin Kim
#   Date      : Sep 22, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : We would like to know how many TCRs are shared across libraries.
#
#   Instruction
#               1. Source("Shared_TCRs_Across_Libraries.R")
#               2. Run the function "shared_tcrs" - specify the input directory and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Shared_TCRs_Across_Libraries.R/Shared_TCRs_Across_Libraries.R")
#               > shared_tcrs(new_TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRredo/",
#                             Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
#                             barcode_dir="C:/Users/hkim8/SJ/SJCAR19/GEXbarcodes_5Oct2020/",
#                             outputDir="./results/TCR/")
###

shared_tcrs <- function(new_TCR_dir="C:/Users/hkim8/SJ/SJCAR19/TCRs_15Oct2020/",
                        barcode_dir="C:/Users/hkim8/SJ/SJCAR19/GEXbarcodes_15Oct2020/",
                        Seurat_RObj_path="./data/JCC212_21Feb2020Aggreg_regress_TCR_clonotyped_PROTO2.Robj",
                        outputDir="./results/TCR_1015/") {
  
  ### load libraries
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    library(gplots, quietly = TRUE)
  }
  if(!require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    require(RColorBrewer, quietly = TRUE)
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
          axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
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
          axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
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
  
  ### get new TCR info file paths
  new_TCR_data_dirs <- list.files(path = new_TCR_dir, pattern = "filtered_contig_annotations.csv$",
                                  full.names = TRUE, recursive = TRUE)
  
  ### load and combine the TCR data
  tcr <- vector("list", length = length(new_TCR_data_dirs))
  names(tcr) <- sapply(basename(new_TCR_data_dirs), function(x) {
    return(substr(x, 1, nchar(x)-36))
  }, USE.NAMES = FALSE)
  for(i in 1:length(new_TCR_data_dirs)) {
    
    ### load filtered contig annotation data
    tcr_data <- read.csv(file = new_TCR_data_dirs[i],
                         stringsAsFactors = FALSE, check.names = FALSE)
    
    ### remove the -* at the end of each barcode.
    tcr_data$barcode <- sapply(strsplit(tcr_data$barcode, split = "-", fixed = TRUE), function(x) x[1])
    
    ### remove the rows that do not have CDR3 sequence
    tcr_data <- tcr_data[which(tcr_data$cdr3 != "None"),]
    
    ### remove redundant rows
    tcr_data <- tcr_data[!duplicated(tcr_data[c("barcode", "cdr3")]),]
    
    ### remove non-productive rows
    tcr_data <- tcr_data[which(tcr_data$productive == "True"),]
    
    ### order by "chain" so that TRA rows come first than TRB rows
    ### and secondly, order by "CDR3 Nucleotide" sequence
    tcr_data <- tcr_data[order(as.character(tcr_data$chain), as.character(tcr_data$cdr3)),]
    
    ### annotate TRA & TRB info to the cdr3 sequences
    tcr_data$cdr3 <- paste0(tcr_data$chain, ":", tcr_data$cdr3)
    tcr_data$cdr3_nt <- paste0(tcr_data$chain, ":", tcr_data$cdr3_nt)
    
    ### now merge different TRA & TRB info to one row
    dups <- which(duplicated(tcr_data$barcode))
    if(length(dups) > 0) {
      temp <- tcr_data[dups,]
      tcr_data <- tcr_data[-dups,]
      rownames(tcr_data) <- tcr_data$barcode
      for(barcode in tcr_data$barcode) {
        idx <- which(temp$barcode == barcode)
        tcr_data[barcode,"cdr3"] <- paste(c(tcr_data[barcode,"cdr3"], temp[idx, "cdr3"]), collapse = ";")
        tcr_data[barcode,"cdr3_nt"] <- paste(c(tcr_data[barcode,"cdr3_nt"], temp[idx, "cdr3_nt"]), collapse = ";")
        tcr_data[barcode,"productive"] <- paste(c(tcr_data[barcode,"productive"], temp[idx, "productive"]), collapse = ";")
        tcr_data[barcode,"reads"] <- paste(c(tcr_data[barcode,"reads"], temp[idx, "reads"]), collapse = ";")
        tcr_data[barcode,"umis"] <- paste(c(tcr_data[barcode,"umis"], temp[idx, "umis"]), collapse = ";")
      }
    }
    
    ### only retain informative columns
    tcr_data <- tcr_data[,c("barcode", "raw_clonotype_id", "cdr3", "cdr3_nt", "reads", "umis", "productive")]
    colnames(tcr_data) <- c("barcode", "raw_clonotype_id", "cdr3_aa", "cdr3_nt", "tcr_reads", "tcr_umis", "tcr_productive")
    
    ### add time & patient info
    temp <- strsplit(basename(new_TCR_data_dirs[i]), split = "_", fixed = TRUE)[[1]]
    px <- temp[3]
    temp <- temp[-c(1,2,3)]
    if(temp[1] == "PB" || temp[1] == "BM") {
      tcr_data$time <- temp[2]
      tcr_data$type <- temp[1]
      tcr_data$px <- px
    } else if(grepl("GMP-redo", temp[1])) {
      tcr_data$time <- "GMP-redo"
      tcr_data$type <- "GMP"
      tcr_data$px <- paste0("SJCAR19-Donor", substring(px, 9))
    } else if(grepl("GMP", temp[1]) || grepl("GMP", px)) {
      tcr_data$time <- "GMP"
      tcr_data$type <- "GMP"
      tcr_data$px <- paste0("SJCAR19-Donor", substring(px, 9))
    } else if(temp[2] == "PB" || temp[2] == "BM") {
      tcr_data$time <- temp[1]
      tcr_data$type <- temp[2]
      tcr_data$px <- px
    } else {
      tcr_data$time <- temp[1]
      tcr_data$type <- "PB"
      tcr_data$px <- px
    }
    if(grepl("HealthyDonor", px)) {
      tcr_data$px <- paste0("SJCAR19-Donor", substring(px, 13))
    }
    
    ### save the TCR data
    tcr[[i]] <- tcr_data
    
  }
  
  ### there is a typo "DMP". change it to "GMP"
  # grep("DMP", names(tcr), fixed = TRUE)
  names(tcr)[9] <- "1662792_JCC212_GMPdonor33"
  
  ### get real tcr names
  real_tcr_names <- sapply(names(tcr), function(x) {
    if(startsWith(x, "X_")) {
      return(substring(x, 3))
    } else {
      return(substring(x, 9))
    }
  }, USE.NAMES = FALSE)
  
  ### match real_tcr_names to the names(real_barcodes)
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-02_Wk1")] <- "JCC212_SJCAR19-02_PB_Wk1"
  real_tcr_names[which(real_tcr_names == "JCC212_GMPdonor30")] <- "JCC212_SJCAR19_GMPdonor30"
  real_tcr_names[which(real_tcr_names == "JCC212_GMPdonor32")] <- "JCC212_SJCAR19_GMPdonor32"
  real_tcr_names[which(real_tcr_names == "JCC212_GMPdonor33")] <- "JCC212_SJCAR19_GMPdonor33"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_Wk-1")] <- "JCC212_SJCAR19-04_PB_Wk-1"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_Wk2")] <- "JCC212_SJCAR19_04_Wk2"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_GMP19003")] <- "JCC212_SJCAR19_04_GMP19003"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_PB_Wk-1")] <- "JCC212_SJCAR19_05_PB_Wk-1"
  real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor32_PreTrans")[1]] <- "JCC212_SJCAR19_Donor32_PreTrans"
  real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor33_PreTrans")[1]] <- "JCC212_SJCAR19_Donor33_PreTrans"
  real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor32_PreTrans")] <- "JCC212_SJCAR19_Donor32_PreTransB"
  real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor33_PreTrans")] <- "JCC212_SJCAR19_Donor33_PreTransB"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_Wk8")] <- "JCC212_SJCAR19_04_Wk8"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_Wk2")] <- "JCC212_SJCAR19_05_Wk2"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_Wk3")] <- "JCC212_SJCAR19_05_Wk3"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_PB_Wk4")] <- "JCC212_SJCAR19_05_PB_Wk4"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_BM_Wk4")] <- "JCC212_SJCAR19_05_BM_Wk4"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-03_PreTrans")] <- "JCC212_SJCAR19_03_PreTransB"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-04_PreTrans")] <- "JCC212_SJCAR19_04_PreTransB"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-00_PreTrans")] <- "JCC212_SJCAR19_00_PreTransB"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_PreTrans")] <- "JCC212_SJCAR19_05_PreTransB"
  real_tcr_names[which(real_tcr_names == "JCC212_HealthyDonor30_PreTrans")] <- "JCC212_SJCAR19_Donor30_PreTransB"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-06_Wk-1_PB")] <- "JCC212_SJCAR19-06_Wk-1"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR-06_GMP19028")] <- "JCC212_SJCAR19-06_GMP19028"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_3mo_PB")] <- "JCC212_SJCAR19-05_PB_3mo"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_3mo_BM")] <- "JCC212_SJCAR19-05_BM_3mo"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-06_PreTrans")] <- "JCC212_SJCAR19-6_PreTrans"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-06_3_mo_BM")] <- "JCC212_SJCAR19-06_3mo_BM"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-08_GMP")] <- "JCC212_SJCAR19-08_GMP19064"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-13_Wk-1")] <- "JCC212_SJCAR19-13_Wk-1_PB"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_Wk0")] <- "JCC212_SJCAR19_05_PB_Wk0"
  real_tcr_names[which(real_tcr_names == "JCC212_SJCAR19-05_Wk1")] <- "JCC212_SJCAR19_05_Wk1"
  
  ### shared TCRs
  shared_tcr_mat <- matrix(NA, length(tcr), length(tcr))
  rownames(shared_tcr_mat) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_tcr_mat) <- rownames(shared_tcr_mat)
  shared_tcr_pct <- matrix(NA, length(tcr), length(tcr))
  rownames(shared_tcr_pct) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_tcr_pct) <- rownames(shared_tcr_pct)
  
  for(i in 1:length(tcr)) {
    for(j in i:length(tcr)) {
      shared_tcr_mat[i,j] <- length(intersect(tcr[[i]]$cdr3_aa, tcr[[j]]$cdr3_aa))
      shared_tcr_pct[i,j] <- 100 * shared_tcr_mat[i,j] / length(unique(tcr[[j]]$cdr3_aa))
    }
  }
  
  ### make full matrices
  shared_tcr_mat[lower.tri(shared_tcr_mat)] <- t(shared_tcr_mat)[lower.tri(shared_tcr_mat)]
  shared_tcr_pct[lower.tri(shared_tcr_pct)] <- t(shared_tcr_pct)[lower.tri(shared_tcr_pct)]
  
  ### because of Jeremy's request to keep the original sample name
  rownames(shared_tcr_mat) <- names(tcr)
  colnames(shared_tcr_mat) <- names(tcr)
  rownames(shared_tcr_pct) <- names(tcr)
  rownames(shared_tcr_pct) <- names(tcr)
  
  ### write out the result
  write.xlsx2(data.frame(shared_tcr_mat, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "Shared_TCRs_Across_All_Libraries.xlsx"),
              sheetName = "Shared_TCR_Numbers", append = FALSE)
  write.xlsx2(data.frame(shared_tcr_pct, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "Shared_TCRs_Across_All_Libraries.xlsx"),
              sheetName = "Shared_TCR_Percentages", append = TRUE)
  
  ### set diagonal value as max
  # diag(shared_tcr_pct) <- max(shared_tcr_pct[lower.tri(shared_tcr_mat)])
  
  ### reset the row & col names
  rownames(shared_tcr_mat) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_tcr_mat) <- rownames(shared_tcr_mat)
  rownames(shared_tcr_pct) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_tcr_pct) <- rownames(shared_tcr_pct)
  
  ### hierarchical clustering functions
  dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
  hclust.ave <- function(x) hclust(x, method="average")
  
  ### get patient number
  px <- sapply(rownames(shared_tcr_pct), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]][2]
    if(grepl("onor", temp)) {
      return(substring(temp, nchar(temp)-1, nchar(temp)))
    } else {
      return(substr(temp, 9, 10))
    }
  })
  px["1779067_SJCAR-06_GMP19028"] <- "06"
  pxf <- factor(px, levels = unique(px[order(as.numeric(px))]))
  
  ### get GMP samples
  gmps <- sapply(rownames(shared_tcr_pct), function(x) {
    if(grepl("GMP", x) || grepl("DMP", x)) {
      if(grepl("redo", x)) {
        return("GMP-REDO")
      } else {
        return("GMP")
      }
    } else {
      return("NON-GMP")
    }
  })
  
  ### set side colors1 - px
  uniqueV <- levels(pxf)
  colors1 <- colorRampPalette(brewer.pal(9,"Reds"))(length(uniqueV))
  names(colors1) <- uniqueV
  
  ### set side colors2 - GMP & GMP-redo
  uniqueV <- c("NON-GMP", "GMP", "GMP-REDO")
  colors2 <- c("deeppink", "deepskyblue", "deepskyblue4")
  names(colors2) <- uniqueV
  
  ### heatmap
  png(paste0(outputDir, "Shared_TCRs_Across_All_Libraries.png"), width = 3200, height = 3000, res = 250)
  par(oma=c(15,0,0,15))
  heatmap.3(as.matrix(shared_tcr_pct),
            main = paste0("Shared_TCRs_Across_All_Libraries"),
            xlab = "", ylab = "", col=rich.colors(200, "blues"),
            scale="none", key=TRUE, keysize=0.8, density.info="density",
            dendrogram = "none", trace = "none",
            labRow = names(tcr), labCol = names(tcr),
            Rowv = TRUE, Colv = TRUE,
            distfun=dist.spear, hclustfun=hclust.ave,
            ColSideColors = cbind(colors1[as.character(px[colnames(shared_tcr_pct)])],
                                  colors2[as.character(gmps[colnames(shared_tcr_pct)])]),
            RowSideColors = t(cbind(colors2[as.character(gmps[rownames(shared_tcr_pct)])],
                                    colors1[as.character(px[rownames(shared_tcr_pct)])])),
            cexRow = 0.5, cexCol = 0.5, na.rm = TRUE)
  legend("left", inset = 0, xpd = TRUE, title = "Patient", legend = names(colors1), fill = colors1, cex = 1.1, box.lty = 0)
  legend("bottomleft", inset = -0.01, xpd = TRUE, title = "GMP", legend = names(colors2), fill = colors2, cex = 0.7, box.lty = 0)
  dev.off()
  
  # ### GMP-specific percentages
  # gmp_idx <- union(which(gmps == "GMP"), which(gmps == "GMP-REDO"))
  # non_gmp_idx <- setdiff(1:nrow(shared_tcr_mat), gmp_idx)
  # gmp_shared_tcr_mat <- shared_tcr_mat[gmp_idx, non_gmp_idx]
  # gmp_shared_tcr_pct <- shared_tcr_pct[gmp_idx, non_gmp_idx]
  # 
  # ### write out the result
  # write.xlsx2(data.frame(gmp_shared_tcr_mat, stringsAsFactors = FALSE, check.names = FALSE),
  #             file = paste0(outputDir, "GMP_Shared_TCRs_Across_All_Libraries.xlsx"),
  #             sheetName = "Shared_TCR_Numbers", append = FALSE)
  # write.xlsx2(data.frame(gmp_shared_tcr_pct, stringsAsFactors = FALSE, check.names = FALSE),
  #             file = paste0(outputDir, "GMP_Shared_TCRs_Across_All_Libraries.xlsx"),
  #             sheetName = "Shared_TCR_Percentages", append = TRUE)
  # 
  # ### heatmap
  # png(paste0(outputDir, "GMP_Shared_TCRs_Across_All_Libraries.png"), width = 2100, height = 1600, res = 130)
  # par(oma=c(0,0,0,25))
  # heatmap.3(as.matrix(gmp_shared_tcr_pct),
  #           main = paste0("Shared_TCRs_Across_All_Libraries"),
  #           xlab = "", ylab = "", col=rich.colors(200, "blues"),
  #           scale="none", key=TRUE, keysize=0.8, density.info="density",
  #           dendrogram = "none", trace = "none",
  #           labRow = rownames(gmp_shared_tcr_pct), labCol = FALSE,
  #           Rowv = TRUE, Colv = TRUE,
  #           distfun=dist.spear, hclustfun=hclust.ave,
  #           ColSideColors = cbind(colors1[as.character(px[colnames(gmp_shared_tcr_pct)])],
  #                                 colors2[as.character(gmps[colnames(gmp_shared_tcr_pct)])]),
  #           RowSideColors = t(cbind(colors2[as.character(gmps[rownames(gmp_shared_tcr_pct)])],
  #                                   colors1[as.character(px[rownames(gmp_shared_tcr_pct)])])),
  #           cexRow = 2, cexCol = 1, na.rm = TRUE)
  # legend("left", inset = -0.01, xpd = TRUE, title = "Patient", legend = names(colors1), fill = colors1, cex = 2, box.lty = 0)
  # legend("bottomleft", inset = -0.02, xpd = TRUE, title = "GMP", legend = names(colors2), fill = colors2, cex = 1, box.lty = 0)
  # dev.off()
  # 
  # ### Shared TCRs - GMP-REDO vs GMP
  # old_gmp_idx <- which(gmps == "GMP")
  # new_gmp_idx <- which(gmps == "GMP-REDO")
  # gmp_shared_tcr_mat2 <- shared_tcr_mat[old_gmp_idx, new_gmp_idx]
  # gmp_shared_tcr_pct2 <- shared_tcr_pct[old_gmp_idx, new_gmp_idx]
  # 
  # ### write out the result
  # write.xlsx2(data.frame(gmp_shared_tcr_mat2, stringsAsFactors = FALSE, check.names = FALSE),
  #             file = paste0(outputDir, "GMP-REDO_vs_GMP_Shared_TCRs.xlsx"),
  #             sheetName = "Shared_TCR_Numbers", append = FALSE)
  # write.xlsx2(data.frame(gmp_shared_tcr_pct2, stringsAsFactors = FALSE, check.names = FALSE),
  #             file = paste0(outputDir, "GMP-REDO_vs_GMP_Shared_TCRs.xlsx"),
  #             sheetName = "Shared_TCR_Percentages", append = TRUE)
  # 
  # ### heatmap
  # png(paste0(outputDir, "GMP-REDO_vs_GMP_Shared_TCRs.png"), width = 2100, height = 2000, res = 130)
  # par(oma=c(25,0,0,25))
  # heatmap.3(as.matrix(gmp_shared_tcr_pct2),
  #           main = paste0("GMP-REDO_vs_GMP_Shared_TCRs"),
  #           xlab = "", ylab = "", col=rich.colors(200, "blues"),
  #           scale="none", key=TRUE, keysize=0.8, density.info="density",
  #           dendrogram = "none", trace = "none",
  #           labRow = rownames(gmp_shared_tcr_pct2), labCol = colnames(gmp_shared_tcr_pct2),
  #           Rowv = FALSE, Colv = FALSE,
  #           distfun=dist.spear, hclustfun=hclust.ave,
  #           ColSideColors = cbind(colors1[as.character(px[colnames(gmp_shared_tcr_pct2)])],
  #                                 colors2[as.character(gmps[colnames(gmp_shared_tcr_pct2)])]),
  #           RowSideColors = t(cbind(colors2[as.character(gmps[rownames(gmp_shared_tcr_pct2)])],
  #                                   colors1[as.character(px[rownames(gmp_shared_tcr_pct2)])])),
  #           cexRow = 2, cexCol = 2, na.rm = TRUE)
  # legend("left", inset = -0.01, xpd = TRUE, title = "Patient", legend = names(colors1), fill = colors1, cex = 2, box.lty = 0)
  # legend("bottomleft", inset = -0.02, xpd = TRUE, title = "GMP", legend = names(colors2), fill = colors2, cex = 1, box.lty = 0)
  # dev.off()
  
  
  #
  ### Barcode matching with GEX
  #
  
  # ### load the Seurat object and save the object name
  # tmp_env <- new.env()
  # load(Seurat_RObj_path, tmp_env)
  # obj_name <- ls(tmp_env)
  # assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  # rm(tmp_env)
  # gc()
  # 
  # ### only retain the needed data and remove the rest
  # metadata <- Seurat_Obj@meta.data
  # rm(Seurat_Obj)
  # gc()
  # 
  # ### get patient, time, and PB/BM info from the TCRs
  # px <- sapply(rownames(shared_tcr_pct), function(x) {
  #   temp <- strsplit(x, split = "_", fixed = TRUE)[[1]][2]
  #   if(grepl("onor", temp)) {
  #     return(substring(temp, nchar(temp)-1, nchar(temp)))
  #   } else {
  #     return(substr(temp, 9, 10))
  #   }
  # })
  # px["1779067_SJCAR-06_GMP19028"] <- "06"
  # 
  # time <- sapply(rownames(shared_tcr_pct), function(x) {
  #   temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
  #   if(length(grep("Wk", temp)) == 1) {
  #     return(temp[grep("Wk", temp)])
  #   } else if(length(grep("GMP", temp)) == 1) {
  #     return("GMP")
  #   } else if(length(grep("mo", temp)) == 1) {
  #     return(temp[grep("mo", temp)])
  #   } else if(length(grep("PreTrans", temp)) == 1) {
  #     return(temp[grep("PreTrans", temp)])
  #   } else {
  #     return("ERROR")
  #   }
  # })
  # time["1977283_SJCAR19-05relapse-2ndInfusionWk1_PB"] <- "Wk1"
  # time["1977284_SJCAR19-05relapse-2ndInfusionWk2_PB"] <- "Wk2"
  # 
  # type <- sapply(rownames(shared_tcr_pct), function(x) {
  #   temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
  #   if(length(grep("PB", temp)) == 1) {
  #     return("PB")
  #   } else if(length(grep("BM", temp)) == 1) {
  #     return("BM")
  #   } else if(length(grep("GMP", temp)) == 1) {
  #     return("GMP")
  #   } else if(length(grep("PreTrans", temp)) == 1) {
  #     return("PB")
  #   } else {
  #     return("PB")
  #   }
  # })
  # 
  # ### combine the info to library
  # library <- sapply(1:nrow(shared_tcr_pct), function(x) {
  #   if(time[x] == "GMP") {
  #     paste0("JCC212_SJCAR19-", px[x], "_", time[x])
  #   } else if(time[x] == "PreTrans" || time[x] == "PreTransB") {
  #     paste0("JCC212_SJCAR19-", px[x], "_", time[x])
  #   } else {
  #     paste0("JCC212_SJCAR19-", px[x], "_", time[x], "_", type[x])
  #   }
  # })
  # library[which(library == "JCC212_SJCAR19-03_PreTrans")] <- "JCC212_SJCAR19-03_PreTransB"
  # library[which(library == "JCC212_SJCAR19-04_PreTrans")] <- "JCC212_SJCAR19-04_PreTransB"
  # library[which(library == "JCC212_SJCAR19-05_PreTrans")] <- "JCC212_SJCAR19-05_PreTransB"
  # 
  # ### shared barcodes
  # shared_barcode_mat <- matrix(NA, length(tcr), length(tcr))
  # rownames(shared_barcode_mat) <- sapply(names(tcr), function(x) {
  #   temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
  #   return(paste0(temp[-2], collapse = "_"))
  # }, USE.NAMES = FALSE)
  # colnames(shared_barcode_mat) <- rownames(shared_barcode_mat)
  # shared_barcode_pct <- matrix(NA, length(tcr), length(tcr))
  # rownames(shared_barcode_pct) <- sapply(names(tcr), function(x) {
  #   temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
  #   return(paste0(temp[-2], collapse = "_"))
  # }, USE.NAMES = FALSE)
  # colnames(shared_barcode_pct) <- rownames(shared_barcode_pct)
  # 
  # for(i in 1:length(tcr)) {
  #   for(j in i:length(tcr)) {
  #     shared_barcode_mat[i,j] <- length(intersect(rownames(tcr[[i]]), metadata$GexCellShort[which(metadata$Library == library[j])]))
  #     shared_barcode_pct[i,j] <- shared_barcode_mat[i,j] / length(union(rownames(tcr[[i]]), metadata$GexCellShort[which(metadata$Library == library[j])]))
  #   }
  # }
  # 
  # ### remove some all zero rows & cols
  # colsum <- apply(shared_barcode_mat, 2, function(x) {
  #   sum(x, na.rm = TRUE)
  # })
  # rIdx <- which(colsum == 0)
  # shared_barcode_mat <- shared_barcode_mat[-rIdx, -rIdx]
  # shared_barcode_pct <- shared_barcode_pct[-rIdx, -rIdx]
  # 
  # ### make full matrices
  # shared_barcode_mat[lower.tri(shared_barcode_mat)] <- t(shared_barcode_mat)[lower.tri(shared_barcode_mat)]
  # shared_barcode_pct[lower.tri(shared_barcode_pct)] <- t(shared_barcode_pct)[lower.tri(shared_barcode_pct)]
  # 
  # ### write out the result
  # write.xlsx2(data.frame(shared_barcode_mat, stringsAsFactors = FALSE, check.names = FALSE),
  #             file = paste0(outputDir, "Shared_Barcodes_Across_All_Libraries.xlsx"),
  #             sheetName = "Shared_Barcodes_Numbers", append = FALSE)
  # write.xlsx2(data.frame(shared_barcode_pct, stringsAsFactors = FALSE, check.names = FALSE),
  #             file = paste0(outputDir, "Shared_Barcodes_Across_All_Libraries.xlsx"),
  #             sheetName = "Shared_Barcodes_Percentages", append = TRUE)
  # 
  # ### diagonal <- set diagonal value as max
  # diag(shared_barcode_pct) <- max(shared_barcode_pct[lower.tri(shared_barcode_mat)])
  # 
  # 
  # ### get patient number
  # px <- sapply(rownames(shared_barcode_pct), function(x) {
  #   temp <- strsplit(x, split = "_", fixed = TRUE)[[1]][2]
  #   if(grepl("onor", temp)) {
  #     return(substring(temp, nchar(temp)-1, nchar(temp)))
  #   } else {
  #     return(substr(temp, 9, 10))
  #   }
  # })
  # px["1779067_SJCAR-06_GMP19028"] <- "06"
  # pxf <- factor(px, levels = unique(px[order(as.numeric(px))]))
  # 
  # ### get GMP samples
  # gmps <- sapply(rownames(shared_barcode_pct), function(x) {
  #   if(grepl("GMP", x) || grepl("DMP", x)) {
  #     if(grepl("redo", x)) {
  #       return("GMP-REDO")
  #     } else {
  #       return("GMP")
  #     }
  #   } else {
  #     return("NON-GMP")
  #   }
  # })
  # 
  # ### set side colors1 - px
  # uniqueV <- levels(pxf)
  # colors1 <- colorRampPalette(brewer.pal(9,"Reds"))(length(uniqueV))
  # names(colors1) <- uniqueV
  # 
  # ### set side colors2 - GMP & GMP-redo
  # uniqueV <- c("NON-GMP", "GMP", "GMP-REDO")
  # colors2 <- c("deeppink", "deepskyblue", "deepskyblue4")
  # names(colors2) <- uniqueV
  # 
  # ### heatmap
  # png(paste0(outputDir, "Shared_Barcodes_Across_All_Libraries.png"), width = 2000, height = 1600, res = 130)
  # par(oma=c(0,0,0,15))
  # heatmap.3(as.matrix(shared_barcode_pct),
  #           main = paste0("Shared_Barcodes_Across_All_Libraries"),
  #           xlab = "", ylab = "", col=rich.colors(200, "blues"),
  #           scale="none", key=TRUE, keysize=0.8, density.info="density",
  #           dendrogram = "none", trace = "none",
  #           labRow = rownames(shared_barcode_pct), labCol = FALSE,
  #           Rowv = TRUE, Colv = TRUE,
  #           distfun=dist.spear, hclustfun=hclust.ave,
  #           ColSideColors = cbind(colors1[as.character(px[colnames(shared_barcode_pct)])],
  #                                 colors2[as.character(gmps[colnames(shared_barcode_pct)])]),
  #           RowSideColors = t(cbind(colors2[as.character(gmps[rownames(shared_barcode_pct)])],
  #                                   colors1[as.character(px[rownames(shared_barcode_pct)])])),
  #           cexRow = 1, cexCol = 1, na.rm = TRUE)
  # legend("left", inset = 0, xpd = TRUE, title = "Patient", legend = names(colors1), fill = colors1, cex = 2, box.lty = 0)
  # legend("bottomleft", inset = -0.01, xpd = TRUE, title = "GMP", legend = names(colors2), fill = colors2, cex = 1, box.lty = 0)
  # dev.off()
  # 
  # ### GMP-specific percentages
  # gmp_idx <- union(which(gmps == "GMP"), which(gmps == "GMP-REDO"))
  # non_gmp_idx <- setdiff(1:nrow(shared_barcode_mat), gmp_idx)
  # gmp_shared_barcode_mat <- shared_barcode_mat[gmp_idx, non_gmp_idx]
  # gmp_shared_barcode_pct <- shared_barcode_pct[gmp_idx, non_gmp_idx]
  # 
  # ### write out the result
  # write.xlsx2(data.frame(gmp_shared_barcode_mat, stringsAsFactors = FALSE, check.names = FALSE),
  #             file = paste0(outputDir, "GMP_Shared_Barcodes_Across_All_Libraries.xlsx"),
  #             sheetName = "Shared_Barcode_Numbers", append = FALSE)
  # write.xlsx2(data.frame(gmp_shared_barcode_pct, stringsAsFactors = FALSE, check.names = FALSE),
  #             file = paste0(outputDir, "GMP_Shared_Barcodes_Across_All_Libraries.xlsx"),
  #             sheetName = "Shared_Barcode_Percentages", append = TRUE)
  # 
  # ### heatmap
  # png(paste0(outputDir, "GMP_Shared_Barcodes_Across_All_Libraries.png"), width = 2100, height = 1600, res = 130)
  # par(oma=c(0,0,0,25))
  # heatmap.3(as.matrix(gmp_shared_barcode_pct),
  #           main = paste0("Shared_Barcodes_Across_All_Libraries"),
  #           xlab = "", ylab = "", col=rich.colors(200, "blues"),
  #           scale="none", key=TRUE, keysize=0.8, density.info="density",
  #           dendrogram = "none", trace = "none",
  #           labRow = rownames(gmp_shared_barcode_pct), labCol = FALSE,
  #           Rowv = TRUE, Colv = TRUE,
  #           distfun=dist.spear, hclustfun=hclust.ave,
  #           ColSideColors = cbind(colors1[as.character(px[colnames(gmp_shared_barcode_pct)])],
  #                                 colors2[as.character(gmps[colnames(gmp_shared_barcode_pct)])]),
  #           RowSideColors = t(cbind(colors2[as.character(gmps[rownames(gmp_shared_barcode_pct)])],
  #                                   colors1[as.character(px[rownames(gmp_shared_barcode_pct)])])),
  #           cexRow = 2, cexCol = 1, na.rm = TRUE)
  # legend("left", inset = -0.01, xpd = TRUE, title = "Patient", legend = names(colors1), fill = colors1, cex = 2, box.lty = 0)
  # legend("bottomleft", inset = -0.02, xpd = TRUE, title = "GMP", legend = names(colors2), fill = colors2, cex = 1, box.lty = 0)
  # dev.off()
  
  
  #
  ### get real barcodes
  #
  
  ### get new barcode info file paths
  barcode_file_paths <- list.files(path = barcode_dir, pattern = "barcodes.tsv.gz$",
                                   full.names = TRUE, recursive = TRUE)
  
  ### load and combine the TCR data
  real_barcodes <- vector("list", length = length(barcode_file_paths))
  names(real_barcodes) <- sapply(basename(barcode_file_paths), function(x) {
    return(substr(x, 1, nchar(x)-31))
  }, USE.NAMES = FALSE)
  for(i in 1:length(barcode_file_paths)) {
    ### load barcode data
    real_barcodes[[i]] <- read.table(file = gzfile(barcode_file_paths[i]),
                                     header = FALSE,
                                     stringsAsFactors = FALSE, check.names = FALSE)
    
    ### remove the -* at the end of each barcode.
    real_barcodes[[i]]$V1 <- sapply(strsplit(real_barcodes[[i]]$V1, split = "-", fixed = TRUE), function(x) x[1])
  }
  
  ### shared barcodes
  shared_barcode_mat <- matrix(NA, length(tcr), length(tcr))
  rownames(shared_barcode_mat) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_mat) <- rownames(shared_barcode_mat)
  shared_barcode_pct <- matrix(NA, length(tcr), length(tcr))
  rownames(shared_barcode_pct) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_pct) <- rownames(shared_barcode_pct)
  
  for(i in 1:length(tcr)) {
    for(j in i:length(tcr)) {
      shared_barcode_mat[i,j] <- length(intersect(rownames(tcr[[i]]), real_barcodes[[real_tcr_names[j]]]$V1))
      shared_barcode_pct[i,j] <- 100* shared_barcode_mat[i,j] / nrow(real_barcodes[[real_tcr_names[j]]])
    }
  }
  
  ### make full matrices
  shared_barcode_mat[lower.tri(shared_barcode_mat)] <- t(shared_barcode_mat)[lower.tri(shared_barcode_mat)]
  shared_barcode_pct[lower.tri(shared_barcode_pct)] <- t(shared_barcode_pct)[lower.tri(shared_barcode_pct)]
  
  ### because of Jeremy's request to keep the original sample name
  rownames(shared_barcode_mat) <- names(tcr)
  colnames(shared_barcode_mat) <- real_tcr_names
  rownames(shared_barcode_pct) <- names(tcr)
  rownames(shared_barcode_pct) <- real_tcr_names
  
  ### write out the result
  write.xlsx2(data.frame(shared_barcode_mat, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "Shared_Barcodes_TCR-GEX_Across_All_Libraries.xlsx"),
              sheetName = "Shared_Barcodes_Numbers", append = FALSE)
  write.xlsx2(data.frame(shared_barcode_pct, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "Shared_Barcodes_TCR-GEX_Across_All_Libraries.xlsx"),
              sheetName = "Shared_Barcodes_Percentages", append = TRUE)
  
  ### diagonal <- set diagonal value as max
  # diag(shared_barcode_pct) <- max(shared_barcode_pct[lower.tri(shared_barcode_mat)])
  
  ### reset the row & col names
  rownames(shared_barcode_mat) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_mat) <- rownames(shared_barcode_mat)
  rownames(shared_barcode_pct) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_pct) <- rownames(shared_barcode_pct)
  
  ### get patient number
  px <- sapply(rownames(shared_barcode_pct), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]][2]
    if(grepl("onor", temp)) {
      return(substring(temp, nchar(temp)-1, nchar(temp)))
    } else {
      return(substr(temp, 9, 10))
    }
  })
  px["1779067_SJCAR-06_GMP19028"] <- "06"
  pxf <- factor(px, levels = unique(px[order(as.numeric(px))]))
  
  ### get GMP samples
  gmps <- sapply(rownames(shared_barcode_pct), function(x) {
    if(grepl("GMP", x) || grepl("DMP", x)) {
      if(grepl("redo", x)) {
        return("GMP-REDO")
      } else {
        return("GMP")
      }
    } else {
      return("NON-GMP")
    }
  })
  
  ### set side colors1 - px
  uniqueV <- levels(pxf)
  colors1 <- colorRampPalette(brewer.pal(9,"Reds"))(length(uniqueV))
  names(colors1) <- uniqueV
  
  ### set side colors2 - GMP & GMP-redo
  uniqueV <- c("NON-GMP", "GMP", "GMP-REDO")
  colors2 <- c("deeppink", "deepskyblue", "deepskyblue4")
  names(colors2) <- uniqueV
  
  ### heatmap
  png(paste0(outputDir, "Shared_Barcodes_TCR-GEX_Across_All_Libraries.png"), width = 3000, height = 2800, res = 220)
  par(oma=c(15,0,0,15))
  heatmap.3(as.matrix(shared_barcode_pct),
            main = paste0("Shared_Barcodes_Across_All_Libraries"),
            xlab = "", ylab = "", col=rich.colors(200, "blues"),
            scale="none", key=TRUE, keysize=0.8, density.info="density",
            dendrogram = "none", trace = "none",
            labRow = names(tcr), labCol = real_tcr_names,
            Rowv = TRUE, Colv = TRUE,
            distfun=dist.spear, hclustfun=hclust.ave,
            ColSideColors = cbind(colors1[as.character(px[colnames(shared_barcode_pct)])],
                                  colors2[as.character(gmps[colnames(shared_barcode_pct)])]),
            RowSideColors = t(cbind(colors2[as.character(gmps[rownames(shared_barcode_pct)])],
                                    colors1[as.character(px[rownames(shared_barcode_pct)])])),
            cexRow = 0.5, cexCol = 0.5, na.rm = TRUE)
  legend("left", inset = 0, xpd = TRUE, title = "Patient", legend = names(colors1), fill = colors1, cex = 1, box.lty = 0)
  legend("bottomleft", inset = -0.01, xpd = TRUE, title = "GMP", legend = names(colors2), fill = colors2, cex = 0.8, box.lty = 0)
  dev.off()
  
  #
  ### barcodes matching TCR-TCR
  #
  
  ### shared barcodes
  shared_barcode_mat <- matrix(NA, length(tcr), length(tcr))
  rownames(shared_barcode_mat) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_mat) <- rownames(shared_barcode_mat)
  shared_barcode_pct <- matrix(NA, length(tcr), length(tcr))
  rownames(shared_barcode_pct) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_pct) <- rownames(shared_barcode_pct)
  
  for(i in 1:length(tcr)) {
    for(j in i:length(tcr)) {
      shared_barcode_mat[i,j] <- length(intersect(rownames(tcr[[i]]), rownames(tcr[[j]])))
      shared_barcode_pct[i,j] <- 100* shared_barcode_mat[i,j] / nrow(tcr[[j]])
    }
  }
  
  ### make full matrices
  shared_barcode_mat[lower.tri(shared_barcode_mat)] <- t(shared_barcode_mat)[lower.tri(shared_barcode_mat)]
  shared_barcode_pct[lower.tri(shared_barcode_pct)] <- t(shared_barcode_pct)[lower.tri(shared_barcode_pct)]
  
  ### because of Jeremy's request to keep the original sample name
  rownames(shared_barcode_mat) <- names(tcr)
  colnames(shared_barcode_mat) <- names(tcr)
  rownames(shared_barcode_pct) <- names(tcr)
  rownames(shared_barcode_pct) <- names(tcr)
  
  ### write out the result
  write.xlsx2(data.frame(shared_barcode_mat, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "Shared_Barcodes_TCR-TCR_Across_All_Libraries.xlsx"),
              sheetName = "Shared_Barcodes_Numbers", append = FALSE)
  write.xlsx2(data.frame(shared_barcode_pct, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "Shared_Barcodes_TCR-TCR_Across_All_Libraries.xlsx"),
              sheetName = "Shared_Barcodes_Percentages", append = TRUE)
  
  ### reset the row & col names
  rownames(shared_barcode_mat) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_mat) <- rownames(shared_barcode_mat)
  rownames(shared_barcode_pct) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_pct) <- rownames(shared_barcode_pct)
  
  ### heatmap
  png(paste0(outputDir, "Shared_Barcodes_TCR-TCR_Across_All_Libraries.png"), width = 3000, height = 2800, res = 220)
  par(oma=c(15,0,0,15))
  heatmap.3(as.matrix(shared_barcode_pct),
            main = paste0("Shared_Barcodes_Across_All_Libraries"),
            xlab = "", ylab = "", col=rich.colors(200, "blues"),
            scale="none", key=TRUE, keysize=0.8, density.info="density",
            dendrogram = "none", trace = "none",
            labRow = names(tcr), labCol = names(tcr),
            Rowv = TRUE, Colv = TRUE,
            distfun=dist.spear, hclustfun=hclust.ave,
            ColSideColors = cbind(colors1[as.character(px[colnames(shared_barcode_pct)])],
                                  colors2[as.character(gmps[colnames(shared_barcode_pct)])]),
            RowSideColors = t(cbind(colors2[as.character(gmps[rownames(shared_barcode_pct)])],
                                    colors1[as.character(px[rownames(shared_barcode_pct)])])),
            cexRow = 0.5, cexCol = 0.5, na.rm = TRUE)
  legend("left", inset = 0, xpd = TRUE, title = "Patient", legend = names(colors1), fill = colors1, cex = 1, box.lty = 0)
  legend("bottomleft", inset = -0.01, xpd = TRUE, title = "GMP", legend = names(colors2), fill = colors2, cex = 0.8, box.lty = 0)
  dev.off()
  
  
  #
  ### barcodes matching GEX-GEX
  #
  
  ### shared barcodes
  shared_barcode_mat <- matrix(NA, length(tcr), length(tcr))
  rownames(shared_barcode_mat) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_mat) <- rownames(shared_barcode_mat)
  shared_barcode_pct <- matrix(NA, length(tcr), length(tcr))
  rownames(shared_barcode_pct) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_pct) <- rownames(shared_barcode_pct)
  
  for(i in 1:length(tcr)) {
    for(j in i:length(tcr)) {
      shared_barcode_mat[i,j] <- length(intersect(real_barcodes[[real_tcr_names[i]]]$V1, real_barcodes[[real_tcr_names[j]]]$V1))
      shared_barcode_pct[i,j] <- 100* shared_barcode_mat[i,j] / nrow(real_barcodes[[real_tcr_names[j]]])
    }
  }
  
  ### make full matrices
  shared_barcode_mat[lower.tri(shared_barcode_mat)] <- t(shared_barcode_mat)[lower.tri(shared_barcode_mat)]
  shared_barcode_pct[lower.tri(shared_barcode_pct)] <- t(shared_barcode_pct)[lower.tri(shared_barcode_pct)]
  
  ### because of Jeremy's request to keep the original sample name
  rownames(shared_barcode_mat) <- real_tcr_names
  colnames(shared_barcode_mat) <- real_tcr_names
  rownames(shared_barcode_pct) <- real_tcr_names
  rownames(shared_barcode_pct) <- real_tcr_names
  
  ### write out the result
  write.xlsx2(data.frame(shared_barcode_mat, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "Shared_Barcodes_GEX-GEX_Across_All_Libraries.xlsx"),
              sheetName = "Shared_Barcodes_Numbers", append = FALSE)
  write.xlsx2(data.frame(shared_barcode_pct, stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "Shared_Barcodes_GEX-GEX_Across_All_Libraries.xlsx"),
              sheetName = "Shared_Barcodes_Percentages", append = TRUE)
  
  ### reset the row & col names
  rownames(shared_barcode_mat) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_mat) <- rownames(shared_barcode_mat)
  rownames(shared_barcode_pct) <- sapply(names(tcr), function(x) {
    temp <- strsplit(x, split = "_", fixed = TRUE)[[1]]
    return(paste0(temp[-2], collapse = "_"))
  }, USE.NAMES = FALSE)
  colnames(shared_barcode_pct) <- rownames(shared_barcode_pct)
  
  ### heatmap
  png(paste0(outputDir, "Shared_Barcodes_GEX-GEX_Across_All_Libraries.png"), width = 3000, height = 2800, res = 220)
  par(oma=c(15,0,0,15))
  heatmap.3(as.matrix(shared_barcode_pct),
            main = paste0("Shared_Barcodes_Across_All_Libraries"),
            xlab = "", ylab = "", col=rich.colors(200, "blues"),
            scale="none", key=TRUE, keysize=0.8, density.info="density",
            dendrogram = "none", trace = "none",
            labRow = real_tcr_names, labCol = real_tcr_names,
            Rowv = TRUE, Colv = TRUE,
            distfun=dist.spear, hclustfun=hclust.ave,
            ColSideColors = cbind(colors1[as.character(px[colnames(shared_barcode_pct)])],
                                  colors2[as.character(gmps[colnames(shared_barcode_pct)])]),
            RowSideColors = t(cbind(colors2[as.character(gmps[rownames(shared_barcode_pct)])],
                                    colors1[as.character(px[rownames(shared_barcode_pct)])])),
            cexRow = 0.5, cexCol = 0.5, na.rm = TRUE)
  legend("left", inset = 0, xpd = TRUE, title = "Patient", legend = names(colors1), fill = colors1, cex = 1, box.lty = 0)
  legend("bottomleft", inset = -0.01, xpd = TRUE, title = "GMP", legend = names(colors2), fill = colors2, cex = 0.8, box.lty = 0)
  dev.off()
  
}
