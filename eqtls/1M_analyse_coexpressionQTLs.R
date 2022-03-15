############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_analyse_coexpressionQTLs.R
# Function: analyse the output of running the co-eQTL mapping
############################################################################################################################


###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(methods)
library(pbapply)
library(Matrix)
require(broom)
library(Seurat)
library(ggplot2)
library(data.table)
library(meta)
require("heatmap.plus")
library(RColorBrewer)
library(UpSetR)
library(foreach)
library(doMC)
library(cowplot)
library(lattice)
library(grid)
library(gridExtra)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################


####################
# Functions        #
####################



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
                      to_na=NULL,
                      extreme=NULL,
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
  if(!is.null(to_na)){
    x[x == to_na] <- NA
  }
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
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
    if(!is.null(extreme)){
      extreme <- extreme
      breaks <- seq(-extreme, extreme, length = breaks)
    }
    else if (!symbreaks)
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
  return(retval)
}



plot_interaction <- function(seurat_object, gene1, gene2, genotype, snp.name, version_chem, output_loc, p.value = NULL, r.value = NULL, sign.cutoff = NULL, use_SCT=T, sct_data=F, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), set_axis_by_max=T, height=5, width=5, extention='.png', remove_prepend=NULL, snp_rename=NULL, xlim=NULL, ylim=NULL){
  # go through all the conditions
  for(condition in intersect(unique(as.character(seurat_object@meta.data$timepoint)), conditions)){
    plot.matrix <- NULL
    # subset to only this condition
    seurat_object_condition <- seurat_object[,seurat_object@meta.data['timepoint'] == condition]
    for(sample in unique(seurat_object_condition@meta.data$assignment)){
      # subset to only this participant
      seurat_object_condition_participant <- seurat_object_condition[, seurat_object_condition@meta.data$assignment == sample]
      # grab the normalized counts
      sample.matrix <- NULL
      if(use_SCT){
        DefaultAssay(seurat_object_condition_participant) <- 'SCT'
        if(sct_data){
          sample.matrix <- data.frame(seurat_object_condition_participant@assays$SCT@data[gene1, ], seurat_object_condition_participant@assays$SCT@data[gene2, ])
        }
        else{
          sample.matrix <- data.frame(seurat_object_condition_participant@assays$SCT@counts[gene1, ], seurat_object_condition_participant@assays$SCT@counts[gene2, ])
        }
      }
      else{
        DefaultAssay(seurat_object_condition_participant) <- 'RNA'
        sample.matrix <- data.frame(seurat_object_condition_participant@assays$RNA@data[gene1, ], seurat_object_condition_participant@assays$RNA@data[gene2, ])
      }
      # remove the rownames, we don't need them
      rownames(sample.matrix) <- NULL
      # set the new colnames
      colnames(sample.matrix) <- c('eqtl', 'interaction')
      # add the participant name
      sample.matrix$sample.name <- sample
      # add the genotype of the participant
      sample.matrix$snp <- as.character(genotype[[sample]])
      if(is.null(plot.matrix)){
        plot.matrix <- sample.matrix
      }
      else{
        plot.matrix <- rbind(plot.matrix, sample.matrix)
      }
    }
    max_y <- max(plot.matrix$eqtl)
    max_x <- max(plot.matrix$interaction)
    # change the SNP into letters if possible
    if(!is.null(snp_rename)){
      plot.matrix$snp <- unlist(snp_rename[plot.matrix$snp])
    }
    # turn the SNP into a factor
    plot.matrix$snp <- as.factor(plot.matrix$snp)
    print(head(plot.matrix))
    # store separate plots
    plots <- list()
    # do the actual plotting
    for (i in 1:length(levels(plot.matrix$snp))) {
      genotype.to.plot <- plot.matrix[plot.matrix$snp == levels(plot.matrix$snp)[i],]

      color.high <- c("lightgreen", "orange", "lightblue")
      color.low <- c("darkgreen", "red3", "blue")

      plots[[i]] <- ggplot(genotype.to.plot, aes(y=eqtl, x=interaction, color=snp, group = sample.name)) +
        geom_point(size=0.3, alpha = 0.2) +
        geom_smooth(method="lm", se=F, size = 0.5) +
        #scale_colour_gradient(name = "sample.name", guide = F,
        #                      low = color.low[i], high = color.high[i]) +
        theme_minimal() +
        theme(panel.background = element_rect(fill = "white", colour = "grey"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "grey", fill=NA, size=2),
              plot.title = element_text(colour = c("#57a350", "#fd7600", "#383bfe")[i], face='bold')) +
        #scale_x_continuous(breaks = round(seq(min(genotype.to.plot$interaction), max(genotype.to.plot$interaction), by = 0.1),1)) +
        scale_color_manual(values=c(c("#57a350", "#fd7600", "#383bfe")[i]), guide = F) +
        ylab("") + xlab("") +
        ggtitle(levels(plot.matrix$snp)[i])
      if(set_axis_by_max){
        plots[[i]] <- plots[[i]] + xlim(c(0, max_x)) + ylim(0, max_y)
      }
      if(!is.null(xlim)){
        plots[[i]] <- plots[[i]] + xlim(xlim)
      }
      if(!is.null(ylim)){
        plots[[i]] <- plots[[i]] + ylim(ylim)
      }
    }
    title <- paste(gene1, "/", gene2, "/", snp.name, '/', condition, '/', version_chem)
    if(!is.null(remove_prepend)){
      title <- paste(gene1, "/", gene2, "/", snp.name, '/', sub(paste('^', remove_prepend, sep = ''),'', condition), '/', version_chem)
    }
    #if(!is.null(p.value) & !is.null(r) & !is.null(sign.cutoff)){
    #  title <- paste(gene1, "/", gene2, "/", snp.name, "\n", signif(p.value, 3), "/", signif(r, 3), "/ cutoff=", sign.cutoff)
    #}

    plot.all <- ggplot(plot.matrix, aes(y=eqtl, x=interaction, color = snp)) +
      geom_point(size=0.7, alpha = 0.4) +
      theme_minimal(base_size = 16) +
      theme(panel.background = element_rect(fill = "white", colour = "grey"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "grey", fill=NA, size=2)) +
      scale_color_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
      ylab(paste(gene1, 'expression')) +
      xlab(paste(gene2, 'expression')) +
      geom_smooth(method="lm") +
      ggtitle(title)

    if(set_axis_by_max){
      plot.all <- plot.all + xlim(c(0, max_x)) + ylim(0, max_y)
    }
    if(!is.null(xlim)){
      plot.all <- plot.all + xlim(xlim)
    }
    if(!is.null(ylim)){
      plot.all <- plot.all + ylim(ylim)
    }

    if(!is.null(p.value)){
      plot.all <- plot.all + annotation_custom(
        grobTree(textGrob(
          label = paste('P = ',p.value[[condition]]), x=0.5, y=0.95, gp=gpar(col="gray", fontsize=16, just=0))
        )
      )
    }
    if(!is.null(r.value)){
      plot.all <- plot.all + annotation_custom(
        grobTree(textGrob(
          label = paste('r = ',r.value[[condition]]), x=0.5, y=0.92, gp=gpar(col="gray", fontsize=16, just=0))
        )
      )
    }

    if (length(plots) == 3){
      right.col <- plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol=1)
      print(plot_grid(plot.all, right.col, rel_widths = c(2,1)))
    } else {
      right.col <- plot_grid(plots[[1]], plots[[2]], ncol=1)
      print(plot_grid(plot.all, right.col, rel_widths = c(2,1)))
    }
    ggsave(paste(output_loc, "interaction", condition, gene1, gene2, version_chem, extention, sep = ''), width=width, height=height)
  }
}

plot_boxplot_and_gt_interactions <- function(genotype_data, mapping_folder, dataset, seurat_object, gene_a, gene_b, snp, monniker='_meta_', cell_type='monocyte', condition='UT', na_to_zero=T, to_numeric=T, snp_rename=NULL, p.value=NULL, r.value=NULL, ylim=NULL, xlim=NULL, ylims=NULL, use_SCT=T, sct_data=F, cor_table=NULL, plot_points=F, violin=F, return_table=F){
  # subset to condition and cell type
  seurat_object_condition <- seurat_object[, (seurat_object@meta.data[['timepoint']] == condition & seurat_object@meta.data[['cell_type_lowerres']] == cell_type)]
  # create plot matrix
  plot.matrix <- NULL
  for(sample in unique(seurat_object_condition@meta.data$assignment)){
    # subset to only this participant
    seurat_object_condition_participant <- seurat_object_condition[, seurat_object_condition@meta.data$assignment == sample]
    # grab the normalized counts
    sample.matrix <- NULL
    if(use_SCT){
      DefaultAssay(seurat_object_condition_participant) <- 'SCT'
      if(sct_data){
        sample.matrix <- data.frame(seurat_object_condition_participant@assays$SCT@data[gene_a, ], seurat_object_condition_participant@assays$SCT@data[gene_b, ])
      }
      else{
        sample.matrix <- data.frame(seurat_object_condition_participant@assays$SCT@counts[gene_a, ], seurat_object_condition_participant@assays$SCT@counts[gene_b, ])
      }
    }
    else{
      DefaultAssay(seurat_object_condition_participant) <- 'RNA'
      sample.matrix <- data.frame(seurat_object_condition_participant@assays$RNA@data[gene_a, ], seurat_object_condition_participant@assays$RNA@data[gene_b, ])
    }
    # remove the rownames, we don't need them
    rownames(sample.matrix) <- NULL
    # set the new colnames
    colnames(sample.matrix) <- c('eqtl', 'interaction')
    # add the participant name
    sample.matrix$sample.name <- sample
    # add the genotype of the participant
    sample.matrix$snp <- as.character(genotype_data[snp,sample])
    if(is.null(plot.matrix)){
      plot.matrix <- sample.matrix
    }
    else{
      plot.matrix <- rbind(plot.matrix, sample.matrix)
    }
  }
  max_y <- max(plot.matrix$eqtl)
  max_x <- max(plot.matrix$interaction)
  # change the SNP into letters if possible
  if(!is.null(snp_rename)){
    plot.matrix$snp <- unlist(snp_rename[plot.matrix$snp])
  }
  # turn the SNP into a factor
  plot.matrix$snp <- as.factor(plot.matrix$snp)
  print(head(plot.matrix))
  # store separate plots
  plots <- list()
  # do the actual plotting
  for (i in 1:length(levels(plot.matrix$snp))) {
    genotype.to.plot <- plot.matrix[plot.matrix$snp == levels(plot.matrix$snp)[i],]

    color.high <- c("lightgreen", "orange", "lightblue")
    color.low <- c("darkgreen", "red3", "blue")

    plots[[i]] <- ggplot(genotype.to.plot, aes(y=eqtl, x=interaction, color=snp, group = sample.name)) + theme_minimal()
    if(!plot_points){
      plots[[i]] <- plots[[i]] + geom_point(alpha = 0.0)
    }
    else{
      plots[[i]] <- plots[[i]] + geom_point(size=0.3, alpha = 0.2)
    }


      plots[[i]] <- plots[[i]] + geom_smooth(method="lm", se=F, size = 0.5) +
      #scale_colour_gradient(name = "sample.name", guide = F,
      #                      low = color.low[i], high = color.high[i]) +
      theme(panel.background = element_rect(fill = "white", colour = "grey"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "grey", fill=NA, size=2),
            plot.title = element_text(colour = c("#57a350", "#fd7600", "#383bfe")[i], face='bold')) +
      #scale_x_continuous(breaks = round(seq(min(genotype.to.plot$interaction), max(genotype.to.plot$interaction), by = 0.1),1)) +
      scale_color_manual(values=c(c("#57a350", "#fd7600", "#383bfe")[i]), guide = F) +
      ylab("") + xlab("") +
      ggtitle(levels(plot.matrix$snp)[i]) +
      geom_smooth(method='lm', data=genotype.to.plot, mapping=aes(y=eqtl, x=interaction), inherit.aes=F, color = "black", se=F)
    if(!is.null(xlim)){
      plots[[i]] <- plots[[i]] + xlim(xlim)
    }
    if(!is.null(ylim)){
      plots[[i]] <- plots[[i]] + ylim(ylim)
    }
  }
  cor_data <- NULL
  if(!is.null(cor_table)){
    cor_data <- data.frame(cor_table)
  }
  else{
    # get the location of the correlated output
    cor_data_loc <- paste(mapping_folder, 'correlationMatrices/', condition, '_', cell_type, '_correlation_matrix_', gene_a, '.', dataset, '.txt', sep = '')
    # read the correlation data
    cor_data <- read.table(cor_data_loc, header = T, row.names = 1)
  }
  # get what we have data for in both
  common_data <- intersect(colnames(cor_data), colnames(genotype_data))
  # grab the data for the gene
  cor_gene <- as.vector(unlist(cor_data[gene_b, common_data]))

  # set to zero if set so
  if(na_to_zero){
    cor_gene[is.na(cor_gene)] <- 0
  }
  # grab the data for the SNP
  snps <- as.vector(unlist(genotype_data[snp, common_data]))
  if(to_numeric){
    #snps <- as.numeric(as.factor(snps))
  }
  # merge the SNP and correlation
  plot_data <- data.frame(participant=common_data, snp=snps, correlation=cor_gene)
  # rename SNPs if requested
  if(!is.null(snp_rename)){
    plot_data$snp <- unlist(snp_rename[plot_data$snp])
    plot_data$snp <- as.factor(plot_data$snp)
  }
  if(return_table){
    # we need to supply the table data, so we can also just return the plot data
    return(list(correlation=plot_data, interaction=plot.matrix))
  }
  else{
    # do plotting
    p <- ggplot(data=plot_data, aes(x=snp,y=correlation, fill = snp), pch=21, size=0.75, alpha=1)
    if(violin){
      p <- p + geom_violin(outlier.shape = NA)
    }
    else{
      p <- p + geom_boxplot(outlier.shape = NA)
    }
    p <- p + theme_bw() +
      ggtitle(paste(snp, gene_a, gene_b, condition, sep = ' ')) +
      scale_fill_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
      ylab(paste(gene_a, '-', gene_b, ' Spearman correlation', sep='')) +
      xlab('genotype') + ylim(c(-0.3, 0.7))
    if(length(unique(plot_data$correlation)) > 1){
      p <- p + geom_jitter(width = 0.1, alpha = 0.2)
    }
    if(!is.null(p.value)){
      p <- p + annotation_custom(
        grobTree(textGrob(
          label = paste('P = ',p.value[[i]]), x=0.5, y=0.95, gp=gpar(col="gray", fontsize=13, just=0))
        )
      )
    }
    if(!is.null(r.value)){
      p <- p + annotation_custom(
        grobTree(textGrob(
          label = paste('r = ',r.value[[i]]), x=0.5, y=0.92, gp=gpar(col="gray", fontsize=13, just=0))
        )
      )
    }
    if(!is.null(ylims)){
      p <- p + ylim(ylims)
    }
    if (length(plots) == 3){
      #txtboxplot(plot_data[plot_data$snp == unique(plot_data$snp)[1], 'correlation'], plot_data[plot_data$snp == unique(plot_data$snp)[2], 'correlation'], plot_data[plot_data$snp == unique(plot_data$snp)[3], 'correlation'])
      top.row <- plot_grid(plots[[1]], plots[[2]], plots[[3]], nrow=1, ncol=3)
      print(plot_grid(top.row, p, rel_heights  = c(1,2), nrow=2, ncol=1))
    } else {
      #txtboxplot(plot_data[plot_data$snp == unique(plot_data$snp)[1], 'correlation'], plot_data[plot_data$snp == unique(plot_data$snp)[2], 'correlation'])
      top.row <- plot_grid(plots[[1]], plots[[2]], nrow=1, ncol=2)
      print(plot_grid(top.row, p, rel_heights = c(1,2), nrow=2, ncol=1))
    }
  }
}


plot_coexpression_qtl <- function(genotype_data, mapping_folder, gene_name, snp, monniker='_meta_', cell_type='monocyte', condition='UT', gene_b=NULL, na_to_zero=T, to_numeric=T, snp_rename=NULL, p.value=NULL, r.value=NULL, ylims=NULL, datasets=c(1,2)){
  # use the supplied gene if possible
  gene_to_use <- gene_b
  if(is.null(gene_to_use)){
    # get the most significant one otherwise
    p_loc <- paste(mapping_folder, gene_name, monniker, cell_type, '_p.tsv', sep = '')
    # read output
    p_vals <- read.table(p_loc, sep = '\t', header = T, row.names = 1)
    # remove the significance threshold one
    p_vals <- p_vals[!rownames(p_vals) %in% c('significance_threshold'), ]
    # order by p value of the condition
    p_vals <- p_vals[order(p_vals[[condition]]), ]
    # grab the first one
    gene_to_use <- rownames(p_vals)[1]
  }
  # we have these per dataset
  plots <- list()
  # we have a meta-analysis of two sets
  for(i in datasets){
    # get the location of the correlated output
    cor_data_loc <- paste(mapping_folder, 'correlationMatrices/', condition, '_', cell_type, '_correlation_matrix_', gene_name, '.', i, '.txt', sep = '')
    # read the correlation data
    cor_data <- read.table(cor_data_loc, header = T, row.names = 1)
    # get what we have data for in both
    common_data <- intersect(colnames(cor_data), colnames(genotype_data))
    # grab the data for the gene
    cor_gene <- as.vector(unlist(cor_data[gene_to_use, common_data]))
    # set to zero if set so
    if(na_to_zero){
      cor_gene[is.na(cor_gene)] <- 0
    }
    # grab the data for the SNP
    snps <- as.vector(unlist(genotype_data[snp, common_data]))
    if(to_numeric){
      #snps <- as.numeric(as.factor(snps))
    }
    # merge the SNP and correlation
    plot_data <- data.frame(participant=common_data, snp=snps, correlation=cor_gene)
    # rename SNPs if requested
    if(!is.null(snp_rename)){
      plot_data$snp <- unlist(snp_rename[plot_data$snp])
      plot_data$snp <- as.factor(plot_data$snp)
    }
    # do plotting
    p <- ggplot(data=plot_data, aes(x=snp,y=correlation, fill = snp), pch=21, size=0.75, alpha=1) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.1, alpha = 0.2) +
      theme_bw() +
      ggtitle(paste(snp, gene_name, gene_to_use, condition, sep = ' ')) +
      scale_fill_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
      ylab(paste(gene_name, '-', gene_to_use, ' Spearman correlation', sep='')) +
      xlab('genotype') + ylim(c(-0.3, 0.7))
    if(!is.null(p.value)){
      p <- p + annotation_custom(
        grobTree(textGrob(
          label = paste('P = ',p.value[[i]]), x=0.5, y=0.95, gp=gpar(col="gray", fontsize=13, just=0))
        )
      )
    }
    if(!is.null(r.value)){
      p <- p + annotation_custom(
        grobTree(textGrob(
          label = paste('r = ',r.value[[i]]), x=0.5, y=0.92, gp=gpar(col="gray", fontsize=13, just=0))
        )
      )
    }
    if(!is.null(ylims)){
      p <- p + ylim(ylims)
    }
    # add to plots
    plots[[i]] <- p
  }
  return(plots)
}

plot_top_hit_per_condition <- function(genotype_data, mappings_folder, mapping_folder_prepend, mapping_folder_append, plot_output_loc, genes, snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), na_to_zero=T, to_numeric=T){
  # get a summary of the numbers
  coeqtl_summary <- summarize_coeqtl_tsvs(paste(mappings_folder, mapping_folder_prepend, sep = ''), paste(mapping_folder_append, '/', sep = ''), genes, cell_types=c('monocyte'), conditions=conditions)[[cell_type]]
  # I like dataframes more than matrices
  coeqtl_summary <- data.frame(coeqtl_summary)
  # can only plot what we actually have output for
  genes_to_plot <- intersect(genes, rownames(coeqtl_summary))
  # plot for each geneA
  for(gene in genes_to_plot){
    # build the path to this specific mappings folder
    #mapping_folder <- paste(mappings_folder, mapping_folder_prepend, monniker, gene, mapping_folder_append, sep = '')
    # really dislike this, but don't want to fix this inconsistency now
    mapping_folder <- paste(mappings_folder, mapping_folder_prepend, gene, substr(monniker, 1, nchar(monniker) - 1), mapping_folder_append, sep = '')
    # get the snp belonging to the gene
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
    # plot for each condition
    for(condition in conditions){
      # check if there is any coeqtl for this condition and gene
      if(!is.na(coeqtl_summary[gene, condition]) & coeqtl_summary[gene, condition] > 0){
        # try to get the plot
        top_gene_condition_plots <- plot_coexpression_qtl(genotype_data=genotype_data, mapping_folder=mapping_folder, gene_name=gene, snp=snp, monniker=monniker, cell_type=cell_type, condition=condition, gene_b=NULL, na_to_zero=na_to_zero, to_numeric=to_numeric)
        # create the plot out locations
        #plot_out_v2 <- paste(plot_output_loc, 'co-expressionQTL_', gene, '_tophit_', condition, '_V2.png', sep = '')
        #top_gene_condition_plots[[1]]
        #ggsave(plot_out_v2)
        # twice of course
        #plot_out_v3 <- paste(plot_output_loc, 'co-expressionQTL_', gene, '_tophit_', condition, '_V3.png', sep = '')
        #top_gene_condition_plots[[2]]
        #ggsave(plot_out_v3)
        # plot both
        ggsave(paste(plot_output_loc, 'co-expressionQTL_', gene, '_tophit_', condition, '.png', sep = ''), arrangeGrob(grobs = top_gene_condition_plots))
      }
    }
  }
}

plot_bottom_hit_per_condition <- function(genotype_data, mappings_folder, mapping_folder_prepend, mapping_folder_append, plot_output_loc, genes, snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), na_to_zero=T, to_numeric=T){
  # get a summary of the numbers
  coeqtl_summary <- summarize_coeqtl_tsvs(paste(mappings_folder, mapping_folder_prepend, sep = ''), paste(mapping_folder_append, '/', sep = ''), genes, cell_types=c('monocyte'), conditions=conditions)[[cell_type]]
  # I like dataframes more than matrices
  coeqtl_summary <- data.frame(coeqtl_summary)
  # can only plot what we actually have output for
  genes_to_plot <- intersect(genes, rownames(coeqtl_summary))
  # plot for each geneA
  for(gene in genes_to_plot){
    # build the path to this specific mappings folder
    #mapping_folder <- paste(mappings_folder, mapping_folder_prepend, monniker, gene, mapping_folder_append, sep = '')
    # really dislike this, but don't want to fix this inconsistency now
    mapping_folder <- paste(mappings_folder, mapping_folder_prepend, gene, substr(monniker, 1, nchar(monniker) - 1), mapping_folder_append, sep = '')
    # get the snp belonging to the gene
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
    # plot for each condition
    for(condition in conditions){
      # check if there is any coeqtl for this condition and gene
      if(!is.na(coeqtl_summary[gene, condition]) & coeqtl_summary[gene, condition] > 0){
        # get the most significant one otherwise
        p_loc <- paste(mapping_folder, gene, monniker, cell_type, '_p.tsv', sep = '')
        # read output
        p_vals <- read.table(p_loc, sep = '\t', header = T, row.names = 1)
        # grab the significance threshold
        sig_thres <- p_vals['significance_threshold', condition]
        # now remove the threshold itself
        p_vals <- p_vals[!rownames(p_vals) %in% c('significance_threshold'), ]
        # get the gene that has the significance threshold
        gene_b <- rownames(p_vals[!is.na(p_vals[[condition]]) & p_vals[[condition]] == sig_thres, ])[1]
        print(gene_b)
        # try to get the plot
        top_gene_condition_plots <- plot_coexpression_qtl(genotype_data=genotype_data, mapping_folder=mapping_folder, gene_name=gene, snp=snp, monniker=monniker, cell_type=cell_type, condition=condition, gene_b=gene_b, na_to_zero=na_to_zero, to_numeric=to_numeric)
        ggsave(paste(plot_output_loc, 'co-expressionQTL_', gene, '_tophit_', condition, '.png', sep = ''), arrangeGrob(grobs = top_gene_condition_plots))
      }
    }
  }
}



plot_top_hit_per_interaction <- function(genotype_data, mappings_folder, mapping_folder_prepend, mapping_folder_append, plot_output_loc, genes, v2_object, v3_object, snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # get a summary of the numbers
  coeqtl_summary <- summarize_coeqtl_tsvs(paste(mappings_folder, mapping_folder_prepend, sep = ''), paste(mapping_folder_append, '/', sep = ''), genes, cell_types=c('monocyte'), conditions=conditions)[[cell_type]]
  # I like dataframes more than matrices
  coeqtl_summary <- data.frame(coeqtl_summary)
  # can only plot what we actually have output for
  genes_to_plot <- intersect(genes, rownames(coeqtl_summary))
  # plot for each geneA
  for(gene in genes_to_plot){
    # build the path to this specific mappings folder
    #mapping_folder <- paste(mappings_folder, mapping_folder_prepend, monniker, gene, mapping_folder_append, sep = '')
    # really dislike this, but don't want to fix this inconsistency now
    mapping_folder <- paste(mappings_folder, mapping_folder_prepend, gene, substr(monniker, 1, nchar(monniker) - 1), mapping_folder_append, sep = '')
    # get the snp belonging to the gene
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
    # plot for each condition
    for(condition in conditions){
      # check if there is any coeqtl for this condition and gene
      if(!is.na(coeqtl_summary[gene, condition]) & coeqtl_summary[gene, condition] > 0){
        # get the most significant one otherwise
        p_loc <- paste(mapping_folder, gene, monniker, cell_type, '_p.tsv', sep = '')
        # read output
        p_vals <- read.table(p_loc, sep = '\t', header = T, row.names = 1)
        # remove the significance threshold one
        p_vals <- p_vals[!rownames(p_vals) %in% c('significance_threshold'), ]
        # order by p value of the condition
        p_vals <- p_vals[order(p_vals[[condition]]), ]
        # grab the first one
        gene_to_use <- rownames(p_vals)[1]
        # try to get the plot
        plot_interaction(v2_object, gene, gene_to_use, genotype_data[snp, ], snp, 'v2', plot_output_loc, p.value = NULL, r = NULL, sign.cutoff = NULL, use_SCT=T, conditions=c(condition))
        plot_interaction(v3_object, gene, gene_to_use, genotype_data[snp, ], snp, 'v3', plot_output_loc, p.value = NULL, r = NULL, sign.cutoff = NULL, use_SCT=T, conditions=c(condition))
      }
    }
  }
}


plot_bottom_hit_per_interaction <- function(genotype_data, mappings_folder, mapping_folder_prepend, mapping_folder_append, plot_output_loc, genes, v2_object, v3_object, snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # get a summary of the numbers
  coeqtl_summary <- summarize_coeqtl_tsvs(paste(mappings_folder, mapping_folder_prepend, sep = ''), paste(mapping_folder_append, '/', sep = ''), genes, cell_types=c('monocyte'), conditions=conditions)[[cell_type]]
  # I like dataframes more than matrices
  coeqtl_summary <- data.frame(coeqtl_summary)
  # can only plot what we actually have output for
  genes_to_plot <- intersect(genes, rownames(coeqtl_summary))
  # plot for each geneA
  for(gene in genes_to_plot){
    # really dislike this, but don't want to fix this inconsistency now
    mapping_folder <- paste(mappings_folder, mapping_folder_prepend, gene, substr(monniker, 1, nchar(monniker) - 1), mapping_folder_append, sep = '')
    # get the snp belonging to the gene
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
    # plot for each condition
    for(condition in conditions){
      # check if there is any coeqtl for this condition and gene
      if(!is.na(coeqtl_summary[gene, condition]) & coeqtl_summary[gene, condition] > 0){
        # get the most significant one otherwise
        p_loc <- paste(mapping_folder, gene, monniker, cell_type, '_p.tsv', sep = '')
        # read output
        p_vals <- read.table(p_loc, sep = '\t', header = T, row.names = 1)
        # grab the significance threshold
        sig_thres <- p_vals['significance_threshold', condition]
        # now remove the threshold itself
        p_vals <- p_vals[!rownames(p_vals) %in% c('significance_threshold'), ]
        # get the gene that has the significance threshold
        gene_to_use <- rownames(p_vals[!is.na(p_vals[[condition]]) & p_vals[[condition]] == sig_thres, ])[1]
        print(gene_to_use)
        # plot for that specific gene
        plot_interaction(v2_object, gene, gene_to_use, genotype_data[snp, ], snp, 'v2', plot_output_loc, p.value = NULL, r = NULL, sign.cutoff = NULL, use_SCT=T, conditions=c(condition))
        plot_interaction(v3_object, gene, gene_to_use, genotype_data[snp, ], snp, 'v3', plot_output_loc, p.value = NULL, r = NULL, sign.cutoff = NULL, use_SCT=T, conditions=c(condition))
      }
    }
  }
}

output_rds_to_tsv <- function(output_loc, tsv_output_prepend, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte')){
  # check each cell type
  for(cell_type in cell_types){
    combined_ct_p <- NULL
    combined_ct_r <- NULL
    # check each condition
    for(condition in conditions){
      tryCatch({
        # read the RDS
        rds <- readRDS(paste(output_loc, condition, '_', cell_type, '.rds', sep=''))
        # grab the P values
        p_df <- data.frame(rds[[2]])
        # set the condition as colname of the single column df
        colnames(p_df) <- c(condition)
        # add to existing df possible
        if(is.null(combined_ct_p)){
          combined_ct_p <- p_df
        }
        else{
          # merge if not the first condition
          combined_ct_p <- merge(combined_ct_p, p_df, by=0, all=T)
          # merge does this thing with rownames we need to correct
          rownames(combined_ct_p) <- combined_ct_p$Row.names
          combined_ct_p$Row.names <- NULL
        }
        # check if there are r values
        if(!is.null(rds[[1]])){
          r_df <- data.frame(rds[[1]])
          # set the condition as colname of the single column df
          colnames(r_df) <- c(condition)
          # add to existing df possible
          if(is.null(combined_ct_r)){
            combined_ct_r <- r_df
          }
          else{
            # merge if not the first condition
            combined_ct_r <- merge(combined_ct_r, r_df, by=0, all=T)
            # merge does this thing with rownames we need to correct
            rownames(combined_ct_r) <- combined_ct_r$Row.names
            combined_ct_r$Row.names <- NULL
          }
        }
      }, error=function(cond) {
        print(paste('issue with', condition, cell_type))
        message(cond)
      })
    }
    print(paste(tsv_output_prepend, cell_type, '_p.tsv', sep=''))
    print(head(combined_ct_p))
    write.table(combined_ct_p, paste(tsv_output_prepend, cell_type, '_p.tsv', sep=''), row.names=T, col.names=T, sep='\t')
    if(!is.null(combined_ct_r)){
      write.table(combined_ct_r, paste(tsv_output_prepend, cell_type, '_r.tsv', sep=''), row.names=T, col.names=T, sep='\t')
    }
  }
}

summarize_coeqtl_tsvs <- function(parent_output_dir_prepend, parent_output_dir_append, genes, snp_probe_mapping, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # we will have multiple summary matrices, one for each cell type
  summary_matrices <- list()
  # check for each cell type
  for(cell_type in cell_types){
    # create an empty matrix
    summary_matrix <- matrix(, ncol = length(conditions), nrow = length(genes), dimnames = list(genes, conditions))
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      v2_out <- paste(parent_output_dir_prepend, gene,'_v2', parent_output_dir_append, gene, '_v2_', cell_type, sep = '')
      v3_out <- paste(parent_output_dir_prepend, gene,'_v3', parent_output_dir_append, gene, '_v3_', cell_type, sep = '')
      meta_out <- paste(parent_output_dir_prepend, gene,'_meta', parent_output_dir_append, gene, '_meta_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      meta_p_loc <- paste(meta_out, '_p.tsv', sep = '')
      v2_r_loc <- paste(v2_out, '_r.tsv', sep = '')
      v3_r_loc <- paste(v3_out, '_r.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        print(meta_p_loc)
        meta_p <- read.table(meta_p_loc, sep = '\t', header = T, row.names = 1)
        #v2_r <- read.table(v2_r_loc, sep = '\t', header = T, row.names = 1)
        #v3_r <- read.table(v3_r_loc, sep = '\t', header = T, row.names = 1)
        # check each condition for this gene
        for(condition in conditions){
          # we can only do something if we have data for this condition
          if(condition %in% colnames(meta_p)){
            # check the significance threshold for this condition
            significance_threshold <- meta_p['significance_threshold', condition]
            # check how many are equal to or smaller than the significance threshold (substracting one for the significance threshold row itself)
            nr_significant <- nrow(meta_p[!is.na(meta_p[[condition]]) & meta_p[[condition]] <= significance_threshold, ]) - 1
            # add this to the matrix
            summary_matrix[gene, condition] <- nr_significant
          }
        }
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })

    }
    # add this 'complete' summary to the list
    summary_matrices[[cell_type]] <- summary_matrix
  }
  # return the result
  return(summary_matrices)
}


read_tsv_results <- function(parent_output_dir_prepend, parent_output_dir_append, genes, cell_types=c('monocyte')){
  # will store it all in a list
  results_per_ct_per_gene <- list()
  # check for each cell type
  for(cell_type in cell_types){
    # make a new list for this cell type
    results_per_ct_per_gene[[cell_type]] <- list()
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      meta_out <- paste(parent_output_dir_prepend, gene,'_meta', parent_output_dir_append, gene, '_meta_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      meta_p_loc <- paste(meta_out, '_p.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        # read the table
        meta_p <- read.table(meta_p_loc, sep = '\t', header = T, row.names = 1)
        # put it in the list
        results_per_ct_per_gene[[cell_type]][[gene]] <- meta_p
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })
    }
  }
  return(results_per_ct_per_gene)
}

read_tsv_results_r <- function(parent_output_dir_prepend, parent_output_dir_append, genes, version, cell_types=c('monocyte')){
  # will store it all in a list
  results_per_ct_per_gene <- list()
  # check for each cell type
  for(cell_type in cell_types){
    # make a new list for this cell type
    results_per_ct_per_gene[[cell_type]] <- list()
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      r_out <- paste(parent_output_dir_prepend, gene,'_', version, parent_output_dir_append, gene, '_', version, '_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      r_loc <- paste(r_out, '_r.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        # read the table
        r <- read.table(r_loc, sep = '\t', header = T, row.names = 1)
        # put it in the list
        results_per_ct_per_gene[[cell_type]][[gene]] <- r
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })
    }
  }
  return(results_per_ct_per_gene)
}

check_ut_overlaps <- function(parent_output_dir_prepend, parent_output_dir_append, genes, version, cell_types=c('monocyte'), stim_conditions=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # we will have multiple summary matrices, one for each cell type
  summary_matrices <- list()
  # first get the results
  results_p <- read_tsv_results(parent_output_dir_prepend, parent_output_dir_append, genes, cell_types)
  # check for each cell type
  for(cell_type in cell_types){
    # create an empty matrix
    summary_matrix <- matrix(, ncol = length(stim_conditions), nrow = length(genes), dimnames = list(genes, stim_conditions))
    # check the output of each gene
    for(gene in genes){
      # grab from the list
      p_table <- results_p[[cell_type]][[gene]]
      # check each condition against UT
      for(stim_condition in stim_conditions){
        # check if both UT and that stim condition are in there
        if('UT' %in% colnames(p_table) & stim_condition %in% colnames(p_table)){
          # get the significance thresholds
          sig_ut_threshold <- p_table['significance_threshold', 'UT']
          sig_stim_threshold <- p_table['significance_threshold', stim_condition]
          # get genes significant for each condition
          sig_ut <- rownames(p_table[!is.na(p_table[['UT']]) & p_table[['UT']] <= sig_ut_threshold, ])
          sig_stim <- rownames(p_table[!is.na(p_table[[stim_condition]]) & p_table[[stim_condition]] <= sig_stim_threshold, ])
          # check which are in both
          sig_both <- intersect(sig_ut, sig_stim)
          # check how many, doing -1, because 'significance_threshold' is always in there
          sig_both_number <- length(sig_both) - 1
          # put that as a result in the matrix
          summary_matrix[gene, stim_condition] <- sig_both_number
        }
      }
    }
    # add this 'complete' overlap summary to the list
    summary_matrices[[cell_type]] <- summary_matrix
  }
  # return the result
  return(summary_matrices)
}


write_significant_genes <- function(parent_output_dir_prepend, parent_output_dir_append, sig_gene_output_loc, genes, version, cell_types=c('monocyte'), conditions=c('UT','X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  for(cell_type in cell_types){
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      meta_out <- paste(parent_output_dir_prepend, gene,'_', version, parent_output_dir_append, gene, '_', version, '_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      meta_p_loc <- paste(meta_out, '_p.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        meta_p <- read.table(meta_p_loc, sep = '\t', header = T, row.names = 1)
        # check each condition for this gene
        for(condition in conditions){
          # we can only do something if we have data for this condition
          if(condition %in% colnames(meta_p)){
            # check the significance threshold for this condition
            significance_threshold <- meta_p['significance_threshold', condition]
            # check how many are equal to or smaller than the significance threshold (substracting one for the significance threshold row itself)
            significant_genes <- rownames(meta_p[!is.na(meta_p[[condition]]) & meta_p[[condition]] <= significance_threshold, ])
            # set output loc
            sig_output_loc <- paste(sig_gene_output_loc, gene, '_', version, '_', cell_type, '_', condition, '_sig_genes.txt', sep = '')
            # write the table
            write.table(significant_genes, sig_output_loc, quote = F, col.names=F, row.names=F)
          }
        }
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })

    }
  }
}

get_gene_lists <- function(tsv_output_loc, use_threshold=T){
  # read the table
  result <- read.table(tsv_output_loc, sep = '\t', header = T, row.names = 1)
  # create a list to store the genes
  sig_genes_per_condition <- list()
  # check each condition, which is per column
  for(condition in colnames(result)){
    if(use_threshold){
      # get for the condition in the column, if it was tested, and if it was below the significance threshold for this condition
      sig_genes <- rownames(result[!is.na(result[[condition]]) & result[[condition]] <= result['significance_threshold', condition], ])
    }
    else{
      sig_genes <- rownames(result[!is.na(result[[condition]]) & result[[condition]] < 0.05, ])
    }
    # remove significance threshold itself
    sig_genes <- sig_genes[!(sig_genes %in% c('significance_threshold'))]
    # add to the list under the condition name
    sig_genes_per_condition[[condition]] <- sig_genes
  }
  return(sig_genes_per_condition)
}


get_gene_list_geneAs <- function(path_prepend, path_append, geneAs){
  # get the significant genes per coeqtl
  sigs_per_geneA <- list()
  for(gene in geneAs){
    # combine to get the path
    full_path <- paste(path_prepend, gene, path_append, sep = '')
    # get the sig genes per condition
    sigs_geneA <- get_gene_lists(full_path)
    # add to the list
    sigs_per_geneA[[gene]] <- sigs_geneA
  }
  return(sigs_per_geneA)
}

get_correlationceofs_geneAs <- function(path_prepend, midpend, path_append, geneAs, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # list for the coefs per geneA
  coefs_per_geneA <- list()
  # check each gene
  for(gene in geneAs){
    # list for coefs per condition
    rs_per_condition <- list()
    # check each condition
    for(condition in conditions){
      # read the correlation coefficients
      rs_loc <- paste(path_prepend, gene, midpend, condition, path_append, sep='')
      try({
        rs <- read.table(rs_loc, sep = '\t', header = T, row.names = 1)
        # add to the list
        rs_per_condition[[condition]] <- rs
      })
    }
    # add to the list
    coefs_per_geneA[[gene]] <- rs_per_condition
  }
  return(coefs_per_geneA)
}


get_geneB_cor_matrix <- function(seurat_object, path_prepend, path_append, geneA, conditions = c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), condition.column = 'timepoint', nthreads = 1, method = 'spearman', verbose = T, cache_loc = './'){
  # combine to get the path
  full_path <- paste(path_prepend, geneA, path_append, sep = '')
  # get the sig genes per condition
  sigs_geneA <- get_gene_lists(full_path)
  # get unique genes in any condition
  sigs_geneA <- unique(as.vector(unlist(sigs_geneA)))
  # create the correlation matrix
  cor_matrix <- get_cor_matrix_per_cond(seurat_object = seurat_object, genes1 = sigs_geneA, genes2 = sigs_geneA, conditions = conditions, condition.column = condition.column, nthreads = nthreads, method = method, verbose = verbose, cache_loc = cache_loc)
  return(cor_matrix)
}


plot_coexpression_vs_coeqtl_direction <- function(output_dir, matrix_prepend, gene, cell_type, col, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # get the most significant one otherwise
  p_loc <- paste(output_dir, gene, '_meta_', cell_type, '_p.tsv', sep = '')
  # get the sig genes per condition
  sigs_per_cond <- get_gene_lists(p_loc)
  # init plot data
  plot_data <- NULL
  # check each condition
  for(condition in conditions){
    try({
      # get the output dir of the condition
      r_output_loc <- paste(output_dir, gene, '_', cell_type, '_', condition, '_rs.tsv', sep='')
      # read the r vals
      rs <- read.table(r_output_loc, sep = '\t', header = T, row.names = 1)
      # get significant gene
      rs <- rs[sigs_per_cond[[condition]], ]
      # get the output dir of the matrix
      matrix_dir <- paste(matrix_prepend, '', gene, '_', cell_type, '_', condition, '_cor_coeqtlgenes.tsv', sep='')
      # read the file
      coex_matrix <- read.table(matrix_dir, sep = '\t', header = T, row.names = 1)
      # of course the colnames in R replace dashes with dots, so we need to do the same
      sigs_geneA_rcolsafe <- gsub('-', '.', sigs_per_cond[[condition]])
      # subset to the genes that were coeqtls
      matrix_coeqtls <- coex_matrix[sigs_per_cond[[condition]], sigs_geneA_rcolsafe]
      # get the coeqtl genes with a positive correlation coefficient
      rs_pos_genes <- rownames(rs[rs[[col]] > 0, ])
      rs_pos_genes_rcolsafe <- gsub('-', '.', rs_pos_genes)
      # get the coeqtl genes with a negative correlation coefficient
      rs_neg_genes <- rownames(rs[rs[[col]] < 0, ])
      rs_neg_genes_rcolsafe <- gsub('-', '.', rs_neg_genes)
      # grab the correlations from the matrix checking coeqtl genes with a positive or negative r
      pos_pos <- as.vector(unlist(matrix_coeqtls[rs_pos_genes, rs_pos_genes_rcolsafe]))
      pos_neg <- as.vector(unlist(matrix_coeqtls[rs_pos_genes, rs_neg_genes_rcolsafe]))
      neg_pos <- as.vector(unlist(matrix_coeqtls[rs_neg_genes, rs_pos_genes_rcolsafe]))
      neg_neg <- as.vector(unlist(matrix_coeqtls[rs_neg_genes, rs_neg_genes_rcolsafe]))
      # remove perfect correlations, as these are with the gene against themselves
      pos_pos <- pos_pos[pos_pos != 1]
      neg_neg <- neg_neg[neg_neg != 1]
      # turn into plot data
      plot_data_condition <- data.frame(correlation=pos_pos, geno_directions=rep('pos-pos', times=length(pos_pos)))
      plot_data_condition <- rbind(plot_data_condition, data.frame(correlation=pos_neg, geno_directions=rep('pos-neg', times=length(pos_neg))))
      plot_data_condition <- rbind(plot_data_condition, data.frame(correlation=neg_pos, geno_directions=rep('neg-pos', times=length(neg_pos))))
      plot_data_condition <- rbind(plot_data_condition, data.frame(correlation=neg_neg, geno_directions=rep('neg-neg', times=length(neg_neg))))
      # add condition
      plot_data_condition$condition <- condition
      # add to the plot data
      if(is.null(plot_data)){
        plot_data <- plot_data_condition
      }
      else{
        plot_data <- rbind(plot_data, plot_data_condition)
      }
    })
  }
  p <- ggplot(data=plot_data, aes(x=geno_directions, y=correlation, fill=geno_directions)) +
    geom_boxplot() +
    geom_jitter() +
    facet_grid(. ~ condition)
  return(p)
}


significant_coeqtl_genes_to_file <- function(input_path_prepend, input_path_append, output_path_prepend, output_path_append, geneAs){
  # first get these genes
  sigs_per_geneA <- get_gene_list_geneAs(input_path_prepend, input_path_append, geneAs)
  # now start writing each geneA set
  for(geneA in names(sigs_per_geneA)){
    # and we need to write per condition
    for(condition in names(sigs_per_geneA[[geneA]])){
      # get those genes
      sig_genes <- sigs_per_geneA[[geneA]][[condition]]
      # set up the output location
      output_loc <- paste(output_path_prepend, 'coeqtls_', geneA, '_', condition, output_path_append, sep='')
      # add the extention if required
      if(!(endsWith('.txt', output_loc))){
        output_loc <- paste(output_loc, '.txt', sep = '')
      }
      # put the significant genes in a dataframe
      sig_genes_df <- data.frame(genes=sig_genes)
      # finally write the table
      write.table(sig_genes_df, output_loc, quote = F, col.names=F, row.names=F)
    }
  }
}

# fix this method, as it's ugly
significant_coeqtl_geneBs_to_file_per_dir <- function(input_path_prepend, input_path_append, output_path_loc, geneAs){
  # first get these genes
  sigs_per_geneA <- get_gene_list_geneAs(input_path_prepend, input_path_append, geneAs)
  # now start writing each geneA set
  for(geneA in names(sigs_per_geneA)){
    # and we need to write per condition
    for(condition in names(sigs_per_geneA[[geneA]])){
      # get those genes
      sig_genes <- sigs_per_geneA[[geneA]][[condition]]
      if(length(sig_genes) > 0){
        # get the output dir of the condition
        r_output_loc <- paste(input_path_prepend, geneA, '_', 'monocyte', '_', condition, '_rs.tsv', sep='')
        # read the Rs
        rs <- read.table(r_output_loc, header = T, row.names = 1, sep = '\t')
        # subset to the significant ones
        rs <- rs[sig_genes, ]
        # add the mean r to the rs
        rs$meanR <- apply(rs, 1, mean)
        # grab the positive direction genes
        pos_genes <- rownames(rs[rs$meanR > 0, ])
        # grab the negative direction genes
        neg_genes <- rownames(rs[rs$meanR < 0, ])
        # write to files
        output_loc_pos <- paste(output_path_loc, 'coeqtls_', geneA, '_', condition, '_pos_meta_monocyte.txt', sep = '')
        output_loc_neg <- paste(output_path_loc, 'coeqtls_', geneA, '_', condition, '_neg_meta_monocyte.txt', sep = '')
        write.table(data.frame(gene=pos_genes), output_loc_pos, row.names = F, col.names = F, quote=F)
        write.table(data.frame(gene=neg_genes), output_loc_neg, row.names = F, col.names = F, quote=F)
      }
    }
  }
}


get_significant_gene_overlap <- function(input_path_prepend, input_path_append, geneAs){
  # first get these genes
  sigs_per_geneA <- get_gene_list_geneAs(input_path_prepend, input_path_append, geneAs)
  # then get the plots per gene
  plots_per_geneA <- list()
  # now start writing each geneA set
  for(geneA in names(sigs_per_geneA)){
    # and we need to write per condition
    sigs_per_condition <- sigs_per_geneA[[geneA]]
    upsetplot <- upset(fromList(sigs_per_condition), nsets = length(names(sigs_per_condition)), order.by = 'freq')
    plots_per_geneA[[geneA]] <- upsetplot
  }
  return(plots_per_geneA)
}

get_r_values <- function(input_path_prepend, input_path_append, gene, snp, snps, cell_type='monocyte', to_numeric=T, sig_only=T, cell_counts=NULL, na_to_zero=T){
  # build the input directory
  output_dir <- paste(input_path_prepend, gene, '_meta', input_path_append, sep = '')
  # get the most significant one otherwise
  p_loc <- paste(output_dir, gene, '_meta_', cell_type, '_p.tsv', sep = '')
  # get the gene lists
  sigs_per_cond <- get_gene_lists(p_loc)
  # grab specifically this SNP
  specific_snp <- snps[snp, , drop=F]
  # save R matrix per condition
  matrix_per_cond <- list()
  for(i in 1:2){
    # check per condition
    for(condition in names(sigs_per_cond)){
      # read the correlation table
      cor_i_loc <- paste(output_dir, 'correlationMatrices/', condition, '_', cell_type, '_correlation_matrix_', gene, '.', i, '.txt', sep = '')
      cor_i <- read.table(cor_i_loc, header = T, row.names = 1)
      # get the snps
      snp_i <- unlist(specific_snp[,match(colnames(cor_i), colnames(specific_snp))])
      if(to_numeric){
        snp_i <- as.numeric(as.factor(snp_i)) - 1
      }
      if(na_to_zero){
        cor_i[is.na(cor_i)] <- 0
      }
      # subset to the genes we care about
      if(sig_only){
        cor_i <- cor_i[sigs_per_cond[[condition]], ]
      }
      # do the interaction analysis
      interaction.statistics <- interaction.regression(cor.matrix = cor_i, eqtl.gene = gene, snp = snp_i, cell.counts = as.vector(unlist(cell_counts[[condition]][[i]][colnames(cor_i)])))
      r.matrix <- (interaction.statistics$statistic / sqrt(length(snp_i) - 2 + interaction.statistics$statistic ** 2))
      r.matrix <- data.frame(r.matrix)
      if(sig_only){
        rownames(r.matrix) <- sigs_per_cond[[condition]]
      }
      else{
        rownames(r.matrix) <- rownames(cor_i)
      }
      colnames(r.matrix) <- i
      # add to existing r matrix if possible
      if(condition %in% names(matrix_per_cond) & sig_only){ # this is way faster, so do this if possible
        matrix_per_cond[[condition]] <- cbind(matrix_per_cond[[condition]], r.matrix)
        print('cbinding')
      }
      else if(condition %in% names(matrix_per_cond) & sig_only == F){ # if checking non-sig, we might not always have the same number of genes
        matrix_per_cond[[condition]] <- merge(matrix_per_cond[[condition]], r.matrix, by=0)
        rownames(matrix_per_cond[[condition]]) <- matrix_per_cond[[condition]]$Row.names
        matrix_per_cond[[condition]]$Row.names <- NULL
      }
      else{
        matrix_per_cond[[condition]] <- r.matrix
        print('setting')
      }
    }
  }
  return(matrix_per_cond)
}


write_r_values_per_gene_and_condition <- function(input_path_prepend, input_path_append, genes, snp_probe_mapping, snps, cell_type='monocyte', to_numeric=F, sig_only=T, cell_counts=NULL){
  # check each gene
  for(gene in genes){
    # get the matching snp
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
    # do the r catching
    r_values_per_per_cond <- get_r_values(input_path_prepend=input_path_prepend, input_path_append=input_path_append, gene=gene, snp=snp, snps=snps, cell_type=cell_type, to_numeric=to_numeric, sig_only=sig_only, cell_counts=cell_counts)
    # build the output dir
    output_dir <- paste(input_path_prepend, gene, '_meta', input_path_append, sep = '')
    # check each condition
    for(condition in names(r_values_per_per_cond)){
      # get the rs for this condition
      rs <- r_values_per_per_cond[[condition]]
      # set the output location
      r_output_loc <- paste(output_dir, gene, '_', cell_type, '_', condition, '_rs.tsv', sep='')
      if(sig_only == F){
        r_output_loc <- paste(output_dir, gene, '_', cell_type, '_', condition, '_rs_full.tsv', sep='')
      }
      # write the result
      write.table(rs, r_output_loc, row.names = T, col.names = T, sep = '\t')
    }
  }
}


write_r_plots_per_gene_and_condition <- function(input_path_prepend, input_path_append, genes, plot_output_loc, cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), sig_only=T){
  # check each gene
  for(gene in genes){
    # get base output dir for gene
    base_output_dir <- paste(input_path_prepend, gene, '_meta', input_path_append, sep = '')
    # get the most significant one otherwise
    p_loc <- paste(base_output_dir, gene, '_meta_', cell_type, '_p.tsv', sep = '')
    # get the gene lists
    sigs_per_cond <- get_gene_lists(p_loc)
    # go through the conditions
    plot_per_condition <- list()
    i <- 1
    for(condition in conditions){
      try({
        # get the output dir of the condition
        r_output_loc <- paste(base_output_dir, gene, '_', cell_type, '_', condition, '_rs.tsv', sep='')
        # read the r vals
        rs <- read.table(r_output_loc, sep = '\t', header = T, row.names = 1)
        if(!is.null(rs) & !is.null(dim(rs)) & nrow(rs) > 0 & ncol(rs) > 1 & sig_only){
          rs[sigs_per_cond[[condition]], ]
        }
        if(!is.null(rs) & !is.null(dim(rs)) & nrow(rs) > 0 & ncol(rs) > 1){
          # add column to denote whether or not the coeqtl was significant
          rs$significant <- 'no'
          if(nrow(rs[rownames(rs) %in% sigs_per_cond[[condition]], ]) > 0){
            rs[rownames(rs) %in% sigs_per_cond[[condition]], ]$significant <- 'yes'
          }
          rs$significant <- as.factor(rs$significant)
          # order by the X1 column
          rs <- rs[order(rs$X1), ]
          p <- ggplot(data=rs, aes(x=X1, y=X2))
          if(sig_only){
            p <- p + geom_point()
          }
          else{
            p <- p + geom_point(aes(colour = significant))
          }
          p <- p + ggtitle(paste(gene, condition, 'v2 vs v3')) + labs(x = 'v2', y = 'v3') + coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) + geom_rect(data=data.frame(X2=c(-1,1), X1=c(-1,1)),aes(xmin = -1, xmax = 0, ymin = -1, ymax = 0), fill = "green", alpha = 0.1) + geom_rect(data=data.frame(X2=c(-1,1), X1=c(-1,1)),aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = "green", alpha = 0.1) + geom_rect(data=data.frame(X2=c(-1,1), X1=c(-1,1)),aes(xmin = -1, xmax = 0, ymin = 0, ymax = 1), fill = "pink", alpha = 0.1) + geom_rect(data=data.frame(X2=c(-1,1), X1=c(-1,1)),aes(xmin = 0, xmax =1, ymin = -1, ymax = 0), fill = "pink", alpha = 0.1)
          plot_per_condition[[i]] <- p
          i <- i + 1
        }
      })
    }
    if(length(plot_per_condition) > 0){
      ggsave(paste(plot_output_loc, 'co-expressionQTL_', gene, '_Rs_', '.png', sep = ''), arrangeGrob(grobs = plot_per_condition), width=12, height=8)
    }
  }
}

interaction.regression <- function(cor.matrix, eqtl.gene, snp, cell.counts) {
  interaction.statistics <- do.call("rbind", apply(cor.matrix, 1, function(x) {
    model <- NULL
    if(is.null(cell.counts)){
      model <- lm(formula = x~snp)
    }
    else{
      model <- lm(formula = x~snp, weights = sqrt(cell.counts))
    }
    return(tidy(model)[2,])
  }))
  return(interaction.statistics)
}

get_most_varying_from_df <- function(dataframe, top_so_many=10, dont_get_least_varying=T){
  # now calculate the sd over this set of columns
  sds <- apply(dataframe, 1, sd, na.rm=T)
  # add the sds as a column
  dataframe$sds <- sds
  # order by the sd
  dataframe <- dataframe[order(dataframe$sds, decreasing = dont_get_least_varying), ]
  # we will return the rownames
  most_varied <- NULL
  # we need to make sure we can return as many rownames as requested
  if(nrow(dataframe) < top_so_many){
    print(paste('requested ', top_so_many, ', but dataframe only has ', nrow(most_varied), ' rows', sep = ''))
    most_varied <- rownames(dataframe)
  }
  else{
    most_varied <- rownames(dataframe)[1:top_so_many]
  }
  return(most_varied)
}

coeqt_gene_pathways_to_df <- function(output_path_prepend, output_path_append, genes, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T, toppfun=T, reactome=F){
  # put all results in a shared DF
  pathway_df <- NULL
  # check each cell type
  for(gene in genes){
    # check each stim
    for(cell_type in cell_types){
      #
      for(condition in conditions){
        try({
          pathways <- NULL
          if(toppfun){
            # paste the filepath together
            filepath <- paste(output_path_prepend, '', gene, '_', condition, '_meta_', cell_type, output_path_append, sep = '')
            # read the file
            pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
            pathways$id_name <- paste(pathways$ID, pathways$Name, sep = '_')
          }
          else if(reactome){
            # paste the filepath together
            filepath <- paste(output_path_prepend, '', gene, '_', condition, '_meta_', cell_type, output_path_append, sep = '')
            # read the file
            pathways <- read.table(filepath, header = T, sep = ',', stringsAsFactors = F, dec = '.') # , colClasses = c('character', 'character', 'character', 'double', 'double', 'double', 'double', 'double', 'double', 'character')
            pathways$id_name <- paste('1', pathways$term_name, sep = '_')
          }
          # create column name
          newcolname <- paste(gene, '_', condition, '_', cell_type, sep = '')
          # get the log2 of the significance value
          if(use_ranking){
            pathways[[newcolname]] <- as.numeric(rownames(pathways))
          }
          else{
            pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
          }

          # reduce to the only two columns we care about
          pathways <- pathways[, c('id_name', newcolname)]
          # join with other pathway files
          if(is.null(pathway_df)){
            # just set as df if the first round through
            pathway_df <- pathways
            pathway_df <- data.table(pathway_df, key = c('id_name'))
          }
          else{
            # otherwise, merge with existing pathways
            pathway_df <- merge(pathway_df, data.table(pathways, key = c('id_name')), by.x='id_name', by.y='id_name', all=T)
          }
        })
      }
    }
  }
  # turn into regular df
  pathway_df <- setDF(pathway_df)
  # set all NA to zero
  pathway_df[is.na(pathway_df)] <- 0
  # set rownames
  rownames(pathway_df) <- pathway_df$id_name
  pathway_df$id_name <- NULL
  # remove rows which amount to zero
  pathway_df <- pathway_df[apply(pathway_df[,-1], 1, function(x) !all(x==0)),]
  return(pathway_df)
}


get_most_shared_pathways <- function(pathway_table, top_so_many=10, use_sd_method=F){
  most_shared <- c()
  # there are different methods to get the most shared pathways
  if(use_sd_method){
    # the SD methods just gets the pathways with the lowest standard deviation
    most_shared <- get_most_varying_from_df(pathway_table, top_so_many = top_so_many, dont_get_least_varying = F)
  }
  else{
    # first get the sum of the rankings over the rows, lower overall ranking means more shared
    summed_rank <- apply(pathway_table, 1, sum)
    # get the top X so many from the pathway table, ordered by this sum of rankings
    most_shared <- rownames(pathway_table[order(summed_rank), ])[1:top_so_many]

    # get the sum of the 3h conditions
    timepoints_3h <- colnames(pathway_table)[grep('3h', colnames(pathway_table))]
    # only if there are two or more columns, can we do this
    if(length(timepoints_3h) > 1){
      summed_rank_3h <- apply(pathway_table[, timepoints_3h], 1, sum)
      most_shared <- c(most_shared, rownames(pathway_table[order(summed_rank_3h), ])[1:top_so_many])
    }
    else{
      print('one or less 3h timepoints')
    }
    # get the sum of the 24h conditions
    timepoints_24h <- colnames(pathway_table)[grep('24h', colnames(pathway_table))]
    # only if there are two or more columns
    if(length(timepoints_24h) > 1){
      summed_rank_24h <- apply(pathway_table[, timepoints_24h], 1, sum)
      most_shared <- c(most_shared, rownames(pathway_table[order(summed_rank_24h), ])[1:top_so_many])
    }
    else{
      print('one or less 24h timepoints')
    }
    # make unique of course
    most_shared <- unique(most_shared)
  }
  return(most_shared)
}


get_most_varied_pathways <- function(pathway_table, top_so_many=10, use_sd_method=F){
  most_varied <- c()
  # most varied in 3h condition
  timepoints_3h <- colnames(pathway_table)[grep('3h', colnames(pathway_table))]
  # most varied in 24h condition
  timepoints_24h <- colnames(pathway_table)[grep('24h', colnames(pathway_table))]
  # most varied in PA condition
  timepoints_pa <- colnames(pathway_table)[grep('hPA', colnames(pathway_table))]
  # most varied in PA condition
  timepoints_ca <- colnames(pathway_table)[grep('hCA', colnames(pathway_table))]
  # most varied in PA condition
  timepoints_mtb <- colnames(pathway_table)[grep('hMTB', colnames(pathway_table))]
  if(use_sd_method){
    # first overall most variation
    most_varied <- get_most_varying_from_df(pathway_table, top_so_many = top_so_many)
    if(length(timepoints_3h) > 1){
      most_varied <- c(most_varied, get_most_varying_from_df(pathway_table[, timepoints_3h], top_so_many = top_so_many))
    }
    if(length(timepoints_24h) > 1){
      most_varied <- c(most_varied, get_most_varying_from_df(pathway_table[, timepoints_24h], top_so_many = top_so_many))
    }
  }
  else{
    if(length(timepoints_3h) > 0 & length(timepoints_24h) > 0){
      # large difference between 3h and 24h
      mean_3h_ranks <- apply(pathway_table[, timepoints_3h], 1, mean)
      mean_24h_ranks <- apply(pathway_table[, timepoints_24h], 1, mean)
      # absolute difference
      mean_rank_diff <- abs(mean_3h_ranks - mean_24h_ranks)
      # grab by varied over abs mean difference
      most_varied <- rownames(pathway_table[order(mean_rank_diff, decreasing = T), ])[1:top_so_many]
    }
    if(length(timepoints_ca) > 0 & length(timepoints_mtb) > 0 & length(timepoints_pa) > 0){
      # get pathogen mean ranks
      mean_pa_ranks <- apply(pathway_table[, timepoints_pa], 1, mean)
      mean_ca_ranks <- apply(pathway_table[, timepoints_ca], 1, mean)
      mean_mtb_ranks <- apply(pathway_table[, timepoints_mtb], 1, mean)
      # turn into separate df
      path_df <- data.frame(pa=mean_pa_ranks, ca=mean_ca_ranks, mtb=mean_mtb_ranks)
      rownames(path_df) <- rownames(pathway_table)
      # use sd method
      most_varied <- c(most_varied, get_most_varying_from_df(path_df, top_so_many = top_so_many))
    }
  }
  print(most_varied)
  most_varied <- unique(most_varied)
  return(most_varied)
}

plot_pathway_sharing <- function(output_path_prepend, output_path_append, gene, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T, toppfun=T, reactome=F){
  # get the pathways
  pathways <- coeqt_gene_pathways_to_df(output_path_prepend=output_path_prepend, output_path_append=output_path_append, genes=c(gene), cell_types=cell_types, conditions=conditions, use_ranking = use_ranking, toppfun =toppfun, reactome = reactome)
  # set the zeroes to max+1
  pathways[pathways==0] <- max(pathways)+1
  # get most shared pathways
  most_shared_pathways <- get_most_shared_pathways(pathways)
  print(most_shared_pathways)
  # subset to those pathways
  pathways <- pathways[most_shared_pathways, ]
  if(length(cell_types) == 1){
    colnames(pathways) <- gsub(paste('_', cell_types[1], sep=''), '', colnames(pathways))
  }
  # plot as heatmap
  heatmap.3(pathways, to_na = max(pathways), dendrogram = 'none', margins=c(10,15))
}


plot_pathway_unique <- function(output_path_prepend, output_path_append, gene, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T, use_sd_method=F){
  # get the pathways
  pathways <- coeqt_gene_pathways_to_df(output_path_prepend=output_path_prepend, output_path_append=output_path_append, genes=c(gene), cell_types=cell_types, conditions=conditions, use_ranking = use_ranking)
  # set the zeroes to max+1
  pathways[pathways==0] <- max(pathways)+1
  # get most shared pathways
  most_unique_pathways <- get_most_varied_pathways(pathways, use_sd_method=use_sd_method)
  print(most_unique_pathways)
  # subset to those pathways
  pathways <- pathways[most_unique_pathways, ]
  if(length(cell_types) == 1){
    colnames(pathways) <- gsub(paste('_', cell_types[1], sep=''), '', colnames(pathways))
  }
  # plot as heatmap
  heatmap.3(pathways, to_na = max(pathways), dendrogram = 'none', margins=c(10,20))
}


write_pathway_plot_tables <- function(output_path_prepend, output_path_append, output_loc, genes, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T, top_so_many=NULL, reactome = T, toppfun = F){
  # each cell type
  for(cell_type in cell_types){
    # check each gene
    for(gene in genes){
      # grab the pathways for that co-eQTL gene
      pathway_df <- plot_pathway_all(output_path_prepend, output_path_append, gene, cell_types=c(cell_type), conditions=conditions, use_ranking=use_ranking, top_so_many=top_so_many, reactome = reactome, toppfun = toppfun, return_table=T)
      # paste together the full path
      output_loc_full <- paste(output_loc, 'pathwas_table_', cell_type, '_', gene, '.tsv', sep = '')
      # write the table
      write.table(pathway_df, output_loc_full, sep = '\t', row.names = T, col.names = T)
    }
  }
}


plot_pathway_all <- function(output_path_prepend, output_path_append, gene, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T, top_so_many=NULL, stringwrap=F, reactome = T, toppfun = F, margins = c(6,25.1), return_table=F){
  # get the pathways
  pathways <- coeqt_gene_pathways_to_df(output_path_prepend=output_path_prepend, output_path_append=output_path_append, genes=c(gene), cell_types=cell_types, conditions=conditions, use_ranking = use_ranking, reactome = reactome, toppfun = toppfun)
  # set the zeroes to max+1
  pathways[pathways==0] <- max(pathways)+1
  if(length(cell_types) == 1){
    colnames(pathways) <- gsub(paste('_', cell_types[1], sep=''), '', colnames(pathways))
  }
  # subset to top so many if requested
  if(!is.null(top_so_many)){
    # which we will use
    top_pathways <- c()
    # check for each column
    for(colname in colnames(pathways)){
      # we don't care about max/na
      nona_pathways <- pathways[!is.na(pathways[[colname]]) & pathways[[colname]] != max(pathways), ]
      # order by this column
      nona_pathways_ordered <- pathways[order(pathways[[colname]]), ]
      # get all if there are less than the top we can get
      if(nrow(nona_pathways_ordered) < top_so_many){
        # add to list
        top_pathways <- c(top_pathways, rownames(nona_pathways_ordered))
      }
      # otherwise we will get so many
      else{
        top_pathways <- c(top_pathways, rownames(nona_pathways_ordered)[1:top_so_many])
      }
    }
    top_pathways <- unique(top_pathways)
    pathways <- pathways[top_pathways, ]
  }
  colnames(pathways) <- gsub(paste(gene, '_', sep = ''), '', colnames(pathways))
  colnames(pathways) <- gsub('^X', '', colnames(pathways))
  rownames(pathways) <- gsub('^\\d+_', '', rownames(pathways))
  if(stringwrap){
    rownames(pathways) <- sapply(rownames(pathways),function(x){paste(strwrap(x,70), collapse="\n")})
  }
  if(return_table){
    # we were asked to supply the tables the figures are based on
    return(pathways)
  }
  else{
    # plot as heatmap
    heatmap.3(pathways, to_na = max(pathways), dendrogram = 'none', margins=margins, KeyValueName = 'ranking of pathway', main = paste('pathways in', gene, '\nco-expressionQTLs'), xlab = 'conditions', ylab = 'pathways')
    
  }
 }

plot_pathway_sharing_genes <- function(output_path_prepend, output_path_append, genes, plot_out_loc, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T){
  # check each gene
  for(gene in genes){
    try({
      # build the output location fully
      full_plot_output_loc <- paste(plot_out_loc, 'coeqtl_pathways_shared_', gene, '.pdf', sep = '')
      # create the plot
      pdf(full_plot_output_loc)
      plot_pathway_sharing(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, gene=gene, cell_types=cell_types, conditions=conditions)
      dev.off()
    })
  }
}

plot_pathway_unique_genes <- function(output_path_prepend, output_path_append, genes, plot_out_loc, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T, use_sd_method=F){
  # check each gene
  for(gene in genes){
    try({
      # build the output location fully
      full_plot_output_loc <- paste(plot_out_loc, 'coeqtl_pathways_unique_', gene, '.pdf', sep = '')
      # create the plot
      pdf(full_plot_output_loc)
      plot_pathway_unique(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, gene=gene, cell_types=cell_types, conditions=conditions, use_sd_method=use_sd_method)
      dev.off()
    })
  }
}

plot_pathway_all_genes <- function(output_path_prepend, output_path_append, genes, plot_out_loc, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T){
  # check each gene
  for(gene in genes){
    try({
      # build the output location fully
      full_plot_output_loc <- paste(plot_out_loc, 'coeqtl_pathways_all_', gene, '.pdf', sep = '')
      # create the plot
      pdf(full_plot_output_loc)
      plot_pathway_all(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, gene=gene, cell_types=cell_types, conditions=conditions)
      dev.off()
    })
  }
}

plot_pathway_all_genes_top <- function(output_path_prepend, output_path_append, genes, plot_out_loc, top_so_many=10, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T){
  # check each gene
  for(gene in genes){
    try({
      # build the output location fully
      full_plot_output_loc <- paste(plot_out_loc, 'coeqtl_pathways_all_', gene, '_top_', top_so_many,'.pdf', sep = '')
      # create the plot
      pdf(full_plot_output_loc)
      plot_pathway_all(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, gene=gene, cell_types=cell_types, conditions=conditions, top_so_many = top_so_many)
      dev.off()
    })
  }
}

get_cors_per_part <- function(seurat_object, gene, pathway, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), condition.column='timepoint', assignment.column='assignment', cell.type.column='cell_type_lowerres'){
  # init table
  cor_score_table <- NULL
  # check each cell type
  for(cell_type in intersect(cell_types, unique(as.character(seurat_object@meta.data[[cell.type.column]])))){
    print(cell_type)
    # subset to cell type
    seurat_object_ct <- seurat_object[, seurat_object@meta.data[[cell.type.column]] == cell_type]
    # check each condition
    for(condition in intersect(conditions, unique(as.character(seurat_object_ct@meta.data[[condition.column]])))){
      print(condition)
      # subset to condition
      seurat_object_ct_cond <- seurat_object_ct[, seurat_object_ct@meta.data[[condition.column]] == condition]
      # finally check for each participant
      for(participant in unique(seurat_object_ct_cond@meta.data[[assignment.column]])){
        print(participant)
        # subset to participant
        seurat_object_ct_cond_part <- seurat_object_ct_cond[, seurat_object_ct_cond[[assignment.column]] == participant]
        # get the correlation
        gene_exp <- as.vector(unlist(seurat_object_ct_cond_part@assays$SCT@counts[gene, ]))
        pathway_score <- as.vector(seurat_object_ct_cond_part@meta.data[, pathway])
        correlation <- 0
        try({
          correlation <- cor(gene_exp, pathway_score, method='spearman')
        })
        # make table of result
        this_cor_score_table <- data.frame(cell_type=c(cell_type), condition=c(condition), participant=c(participant), cor=c(correlation), stringsAsFactors=F)
        # add to overarching table if possible
        if(is.null(cor_score_table)){
          cor_score_table <- this_cor_score_table
        }
        else{
          cor_score_table <- rbind(cor_score_table, this_cor_score_table)
        }
      }
    }
    return(cor_score_table)
  }
}

add_gt_to_cors <- function(cors_per_part_df, snps, gene, snp_probe_mapping, assignment.column='participant', new_gt_colname=NULL){
  # get the gene that belongs to the SNP
  snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
  # get the participants, not only unique, so we can grab and paste easily later on
  participants <- cors_per_part_df[[assignment.column]]
  # grab for those participants, that SNP
  snps_participants <- as.vector(unlist(snps[snp, participants]))
  # should be the same length as correlations, so we can just cbind
  cors_per_part_df[[snp]] <- snps_participants
  # if requested, we'll name this column differently
  if(!is.null(new_gt_colname)){
    # add with new name
    cors_per_part_df[[new_gt_colname]] <- cors_per_part_df[[snp]]
    # remove old column
    cors_per_part_df[[snp]] <- NULL
  }
  return(cors_per_part_df)
}

plot_gene_to_pathway <- function(cors_per_part_df, cor_column='cor', gt_column='snp', title='pathway vs gene'){
  # create plot data, new df so we have consistent column names
  plot_data <- data.frame(snp=cors_per_part_df[[gt_column]], correlation=cors_per_part_df[[cor_column]])
  # create plot
  p <- ggplot(data=plot_data, aes(x=snp,y=correlation, fill = snp)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.2) +
    theme_bw() +
    ggtitle(title)
  # add to plots
  return(p)
}


create_gene_correlations <- function(seurat_object, genes=NULL, method='spearman'){
  # get the genes to check
  genes_to_check <- genes
  # if no genes are supplied, then check all of them
  if(is.null(genes_to_check)){
    genes_to_check <- rownames(seurat_object)
  }
    # create the matrix to fill in
    correlations <- matrix(, ncol=length(genes_to_check), nrow=length(genes_to_check), dimnames=list(genes_to_check, genes_to_check))
    # check each gene against each gene
    for(i in 1:length(genes_to_check)){
      for(i2 in i:length(genes_to_check)){
        # calculate the correlations
        try({
          # grab the gene names
          gene1 <- genes_to_check[i]
          gene2 <- genes_to_check[i2]
          # check the correlation
          correlation <- cor(as.vector(unlist(seurat_object@assays$SCT@counts[gene1, ])), as.vector(unlist(seurat_object@assays$SCT@counts[gene2, ])), method=method)
          # and add that to the correlation matrix
          correlations[gene1, gene2] <- correlation
          correlations[gene2, gene1] <- correlation
        })
      }
    }
}

create_gene_correlations_mt <- function(seurat_object, genes=NULL, nthreads=1, method='spearman', verbose=F){
  # we'll try this in parallel
  registerDoMC(nthreads)
  # get the genes to check
  genes_to_check <- genes
  # if no genes are supplied, then check all of them
  if(is.null(genes_to_check)){
    genes_to_check <- rownames(seurat_object)
  }
  # create the matrix to fill in
  correlations <- matrix(, ncol=length(genes_to_check), nrow=length(genes_to_check), dimnames=list(genes_to_check, genes_to_check))
  # check each gene against each gene
  foreach(i=1:length(genes_to_check)) %dopar% {
    for(i2 in i:length(genes_to_check)){
      if(verbose & i2 %% 1000 == 0){
        print(paste('doing', i, i2))
      }
      # calculate the correlations
      try({
        # grab the gene names
        gene1 <- genes_to_check[i]
        gene2 <- genes_to_check[i2]
        # check the correlation
        correlation <- cor(as.vector(unlist(seurat_object@assays$SCT@counts[gene1, ])), as.vector(unlist(seurat_object@assays$SCT@counts[gene2, ])), method = method)
        # and add that to the correlation matrix
        correlations[gene1, gene2] <- correlation
        correlations[gene2, gene1] <- correlation
      })
    }
  }
  return(correlations)
}

create_cor_matrix_mt <- function(seurat_object, genes1=NULL, genes2=NULL, nthreads=1, method='spearman', verbose=F, cache_loc='./'){
  DefaultAssay(seurat_object) <- 'SCT'
  # we'll try this in parallel
  registerDoMC(nthreads)
  # set the genes to check
  genes_x <- genes1
  genes_y <- genes2
  # check if they were supplied, otherwise grab use the genes available in the object
  if(is.null(genes_x)){
    genes_x <- rownames(seurat_object)
  }
  if(is.null(genes_y)){
    genes_y <- rownames(seurat_object)
  }
  # we'll cache to a file, let's create a random string, so we won't get into trouble running multiple instances
  random_string <- get_random_strings(1)
  # check the genes
  foreach(i=1:length(genes_x)) %dopar% {
    # grab the gene1
    gene1 <- genes_x[i]
    # create the matrix to fill in
    correlations <- matrix(, ncol=length(genes_y), nrow=length(1), dimnames=list(c(gene1), genes_y))
    # check against each other gene
    for(i2 in 1:length(genes_y)){
      if(verbose & i2 %% 1000 == 0){
        print(paste('doing', i, i2))
      }
      # grab the gene2
      gene2 <- genes_y[i2]
      # check the correlation
      try({
        #correlation <- cor(as.vector(unlist(seurat_object@assays$SCT@counts[gene1, ])), as.vector(unlist(seurat_object@assays$SCT@counts[gene2, ])), method = method)
        correlation <- rcorr(as.vector(unlist(seurat_object@assays$SCT@counts[gene1, ])), as.vector(unlist(seurat_object@assays$SCT@counts[gene2, ])), type = method)$r['x', 'y']
        # and add that to the correlation matrix
        correlations[gene1, gene2] <- correlation
      })
    }
    # write for this gene to a file
    cache_file_loc <- paste(cache_loc, random_string, '_', gene1, '.tsv', sep = '')
    write.table(correlations, cache_file_loc, col.names=T, row.names=T, sep = '\t')
  }
  # merge the genes together again
  full_cor_table <- NULL
  for(gene in genes_x){
    # construct the location again
    cache_file_loc <- paste(cache_loc, random_string, '_', gene, '.tsv', sep = '')
    # read the table
    correlations <- read.table(cache_file_loc, header=T, row.names=1, sep='\t')
    # merge
    if(is.null(full_cor_table)){
      full_cor_table <- correlations
    }
    else{
      full_cor_table <- rbind(full_cor_table, correlations)
    }
  }
  return(full_cor_table)
}

get_random_strings <- function(n = 1) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}


get_cor_matrix_per_cond <- function(seurat_object, genes1=NULL, genes2=NULL, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), condition.column='timepoint', nthreads=1, method='spearman', verbose=F, cache_loc='./'){
  # store correlations per condition
  cors_per_cond <- list()
  # check each condition
  for(condition in conditions){
    # subset to condition
    seurat_object_condition <- seurat_object[, seurat_object@meta.data[[condition.column]] == condition]
    # create correlation matrix
    cors_condition <- create_cor_matrix_mt(seurat_object_condition, genes1=genes1, genes2=genes2, nthreads=nthreads, method=method, verbose=verbose, cache_loc=cache_loc)
    # store in list
    cors_per_cond[[condition]] <- cors_condition
  }
  return(cors_per_cond)
}


get_avg_cor_per_gt <- function(cor_matrix_loc, genotype_data, snp, na_to_zero=T){
  # read the correlation matrix
  cor_matrix <- read.table(cor_matrix_loc, header=T, row.names=1)
  # subset to the snp
  geno_snp <- genotype_data[snp, colnames(cor_matrix)]
  geno_snp <- data.frame(t(geno_snp))
  # create dataframe to get the average cor per get
  avg_cor_per_gts <- NULL
  # get the alleles
  alleles_combos <- unique(geno_snp[[snp]])
  # check each genotype
  for(alleles in alleles_combos){
    # get the participants with this genotype
    parts_alleles <- rownames(geno_snp[as.character(geno_snp[[snp]]) == alleles, , drop = F])
    # we can only do this if there are people with this genotype
    if(length(parts_alleles) > 0){
      # subset to that allele
      cor_matrix_genotype <- cor_matrix[, parts_alleles, drop = F]
      if(na_to_zero){
        cor_matrix_genotype[is.na(cor_matrix_genotype)] <- 0
      }
      # get the mean expression of each gene in this group
      avg_cor <- apply(cor_matrix_genotype, 1, mean)
      # add to the cor df
      if(is.null(avg_cor_per_gts)){
        avg_cor_per_gts <- avg_cor
      }
      else{
        avg_cor_per_gts <- cbind(avg_cor_per_gts, avg_cor)
      }
    }
  }
  # set names
  rownames(avg_cor_per_gts) <- rownames(cor_matrix)
  colnames(avg_cor_per_gts) <- alleles_combos
  return(avg_cor_per_gts)
}

write_avg_cor_per_gt_per_cond <- function(cor_matrices_loc, gene, genotype_data, snp, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), sets=c(1, 2), na_to_zero=F){
  # check each cell type
  for(cell_type in cell_types){
    # check each condition
    for(condition in conditions){
      # check each set
      for(set in sets){
        # paste together the reading path
        input_loc <- paste(cor_matrices_loc, condition, '_', cell_type, '_correlation_matrix_', gene, '.', set, '.txt', sep = '')
        # paste together the output loc
        output_loc <- paste(cor_matrices_loc, condition, '_', cell_type, '_avg_cor_matrix_', gene, '.', set, '.tsv', sep = '')
        # get the average correlation matrix
        avg_cor_matrix <- get_avg_cor_per_gt(input_loc, genotype_data, snp, na_to_zero=na_to_zero)
        # write the result
        write.table(avg_cor_matrix, output_loc, sep = '\t', row.names = T, col.names = T, quote = F)
      }
    }
  }
}


write_avg_cor_per_coeqtl_gene <- function(output_loc_prepend, output_loc_append, coeqtl_genes, genotype_data, snp_probe_mapping, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), sets=c(1, 2), na_to_zero=F){
  # check each coeqtl gene
  for(coeqtl_gene in coeqtl_genes){
    # get the SNP that belongs to that gene
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == coeqtl_gene, ]$snp[1]
    # paste together the location of the matrices
    cor_matrices_loc <- paste(output_loc_prepend, coeqtl_gene, output_loc_append, '/correlationMatrices/', sep = '')
    # put in the work for the average correlation
    write_avg_cor_per_gt_per_cond(cor_matrices_loc = cor_matrices_loc, gene = coeqtl_gene, genotype_data = genotype_data, snp = snp, cell_types = cell_types, conditions = conditions, sets = sets, na_to_zero=na_to_zero)
  }
}

plot_baseline_correlations <- function(correlations_matrix_per_condition, geneA, geneBs, conditions=NULL, title='correlations', xlab='condition', ylab='correlation'){
  conditions.to.check <- names(correlations_matrix_per_condition)
  # subset to requested genes if asked to do so
  if(!is.null(conditions)){
    conditions.to.check <- intersect(conditions.to.check, conditions)
  }
  # init plot data
  plot_data <- NULL
  # get correlations for each condition
  for(condition in conditions.to.check){
    # get the correct correlation matrix
    cor_matrix_condition <- correlations_matrix_per_condition[[condition]]
    # get the specific correlations we wanted
    correlations <- as.vector(unlist(cor_matrix_condition[geneA, intersect(gsub('-', '.', geneBs), colnames(cor_matrix_condition))]))
    # create plot data
    plot_data_condition <- data.frame(correlation=correlations, condition=rep(condition, times=length(correlations)))
    # add to overarching plot data
    if(is.null(plot_data)){
      plot_data <- plot_data_condition
    }
    else{
      plot_data <- rbind(plot_data, plot_data_condition)
    }
  }
  plot_data$condition <- as.factor(plot_data$condition)
  levels(plot_data$condition) <- conditions
  cc <- get_color_coding_dict()
  #colScale <- scale_fill_manual(name = plot_data$condition, values = unlist(cc[conditions]))
  ggplot(plot_data, aes(x=condition, y=correlation, fill=condition)) +
    geom_boxplot() +
    #colScale +
    ggtitle(title) +
    labs(y = ylab, x=ylab) +
    geom_jitter(width = 0.1, alpha = 0.2)
}

plot_coeqtl_baseline_correlations <- function(geneA, tsv_loc, coexp_loc_prepend, coexp_append='_cor_coeqtlgenes.tsv', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), do_ggparcoord=F){
  # get the significant genes per condition
  sigs_per_cond <- get_gene_lists(tsv_loc)
  # subset to conditions we care about
  sigs_per_cond <- sigs_per_cond[intersect(names(sigs_per_cond), conditions)]
  # merge these into a single vector
  sigs <- unlist(sigs_per_cond)
  # make unique, because some genes are a coeqtl in multiple conditions
  sigs <- unique(sigs)
  # read the coexpression matrices
  correlations_matrix_per_condition <- list()
  for(condition in conditions){
    # paste together the path
    coexp_path <- paste(coexp_loc_prepend, condition, coexp_append, sep='')
    # read and add to list
    try({
      coexp_condition <- read.table(coexp_path, header=T, row.names=1, sep = '\t')
      correlations_matrix_per_condition[[condition]] <- coexp_condition
    })
  }
  # actually plot now
  if(do_ggparcoord){
    plot_baseline_correlations_ggparcoord(correlations_matrix_per_condition, geneA, sigs, conditions)
  }
  else{
    plot_baseline_correlations(correlations_matrix_per_condition, geneA, sigs, conditions)
  }

}

plot_baseline_correlations_ggparcoord <- function(correlations_matrix_per_condition, geneA, geneBs, conditions=NULL){
  full_matrix <- NULL
  conditions.to.check <- names(correlations_matrix_per_condition)
  # subset to requested genes if asked to do so
  if(!is.null(conditions)){
    conditions.to.check <- intersect(conditions.to.check, conditions)
  }
  # get correlations for each condition
  for(condition in conditions.to.check){
    # get the correct correlation matrix
    cor_matrix_condition <- correlations_matrix_per_condition[[condition]]
    # get the specific correlations we wanted
    correlations <- as.vector(unlist(cor_matrix_condition[geneA, gsub('-', '.', geneBs)]))
    # turn into df
    condition_df <- data.frame(correlation=correlations)
    colnames(condition_df) <- condition
    # add to dataframe
    if(is.null(full_matrix)){
      full_matrix <- condition_df
    }
    else{
      full_matrix <- cbind(full_matrix, condition_df)
    }
  }
  ggparcoord(full_matrix)
}

plot_coeqtl_baseline_correlations_per_gt <- function(geneA, tsv_loc, i, prepend, midpend='_monocyte_avg_cor_matrix_', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), avg_pos_only=F, avg_neg_only=F, coef_path_prepend='', coef_path_append='_rs.tsv', positive_gen_only=F, negative_gen_only=F, na_to_zero=F, ut_sig_only=F, stim_sig_only=F, stim_sig=F, ut_sig=F, stim_specific_only=F, color_by_condition=F, color_by_genotype=F, remove_prepend='^X', sig_change_loc=NULL, snp_rename=NULL, apply_paper_settings=T, xlab=NULL, ylab=NULL, remove_outliers=F, ylim=NULL){
  # get the significant genes per condition
  sigs_per_cond <- get_gene_lists(tsv_loc)
  # subset to conditions we care about
  sigs_per_cond <- sigs_per_cond[intersect(names(sigs_per_cond), conditions)]
  # merge these into a single vector
  sigs <- unlist(sigs_per_cond)
  # make unique, because some genes are a coeqtl in multiple conditions
  sigs <- unique(sigs)
  # init plot data
  plot_data <- NULL
  title <- 'average gene expression'
  runthrough1 <- T
  for(condition in conditions){
    try({
      # paste together the path
      coexp_path <- paste(prepend, condition, midpend, geneA, '.', i, '.tsv', sep='')
      # load file
      cors <- read.table(coexp_path, header=T, row.names=1, sep='\t')
      # turn NA into zero if requested
      if(na_to_zero){
        cors[is.na(cors)] <- 0
      }
      # copy for subsetting
      sigs_per_cond_nout <- sigs_per_cond
      sigs_per_cond_nout[['UT']] <- NULL
      # subset to only positive or negative correlations
      if(avg_pos_only){
        mean_row <- apply(cors, 1, mean)
        cors <- cors[mean_row >= 0, ]
        if(runthrough1){title <- paste(title, 'positive correlation')}
      }
      if(avg_neg_only){
        mean_row <- apply(cors, 1, mean)
        cors <- cors[mean_row <= 0, ]
        if(runthrough1){title <- paste(title, 'negative correlation')}
      }
      if(positive_gen_only){
        dirs <- read.table(paste(coef_path_prepend, condition, coef_path_append, sep = ''), header=T, row.names=1, sep='\t')
        cors <- cors[rownames(dirs[dirs[[i]] > 0, ]), ]
        if(runthrough1){title <- paste(title, 'positive gt effect')}
      }
      if(negative_gen_only){
        dirs <- read.table(paste(coef_path_prepend, condition, coef_path_append, sep = ''), header=T, row.names=1, sep='\t')
        cors <- cors[rownames(dirs[dirs[[i]] < 0, ]), ]
        if(runthrough1){title <- paste(title, 'negative gt effect')}
      }
      if(ut_sig_only){
        cors <- cors[rownames(cors) %in% sigs_per_cond[['UT']] & !(rownames(cors) %in% as.vector(unlist(sigs_per_cond_nout[conditions]))), ]
        if(runthrough1){title <- paste(title, 'significant in UT only')}
      }
      if(stim_sig_only){
        cors <- cors[!(rownames(cors) %in% sigs_per_cond[['UT']]) & rownames(cors) %in% as.vector(unlist(sigs_per_cond_nout[conditions])), ]
        if(runthrough1){title <- paste(title, 'significant in stim only')}
      }
      if(stim_specific_only){
        cors <- cors[!(rownames(cors) %in% sigs_per_cond[['UT']]) & rownames(cors) %in% sigs_per_cond[[condition]], ]
        if(runthrough1){title <- paste(title, 'significant specific stim only')}
      }
      if(stim_sig){
        cors <- cors[(rownames(cors) %in% sigs_per_cond[[condition]]), ]
        if(runthrough1){title <- paste(title, 'significant in stim')}
      }
      if(ut_sig){
        cors <- cors[(rownames(cors) %in% sigs_per_cond[['UT']]), ]
        if(runtrough1){title <- paste(title, 'significant in UT')}
      }
      if(!is.null(sig_change_loc)){
        if(!condition == 'UT'){
          changes <- read.table(sig_change_loc, sep = '\t', header = T, row.names = 1)
          sig_chances_cond <- rownames(changes[!is.na(changes[[paste('UT_vs_', condition, sep = '')]]) & changes[[paste('UT_vs_', condition, sep = '')]], ])
          cors <- cors[rownames(cors) %in% sig_chances_cond, ]
        }
        else{
          all_sigs <- unique(as.vector(unlist(get_gene_lists(sig_change_loc, use_threshold = F))))
          cors <- cors[rownames(cors) %in% all_sigs, ]
        }
        if(runthrough1){title <- paste(title, 'nominal effect change')}
      }
      # check each genotype
      for(geno in colnames(cors)){
        cors_geno_and_genes <- cors[sigs, geno]
        # turn into dataframe
        cors_geno_and_genes_df <- data.frame(cor=cors_geno_and_genes, stringsAsFactors = F)
        cors_geno_and_genes_df$condition <- condition
        cors_geno_and_genes_df$geno <- geno
        cors_geno_and_genes_df$gene <- sigs
        # add to df
        if(is.null(plot_data)){
          plot_data <- cors_geno_and_genes_df
        }
        else{
          plot_data <- rbind(plot_data, cors_geno_and_genes_df)
        }
      }
      runthrough1 <- F
    })
  }
  # some magic to order the names and remove prepends
  condition_names <- conditions
  if(!is.null(remove_prepend)){
    plot_data$condition <- gsub(remove_prepend, '', plot_data$condition)
    condition_names <- gsub(remove_prepend, '', condition_names)
  }
  plot_data$condition <- factor(plot_data$condition, levels=condition_names)
  # rename snps if requested
  if(!is.null(snp_rename)){
    plot_data$geno <- unlist(snp_rename[plot_data$geno])
    plot_data$geno <- as.factor(plot_data$geno)
  }
  # sort conditions
  p <- ggplot(plot_data, aes(condition, cor, fill=geno)) + geom_boxplot()
  if(remove_outliers){
    p <- ggplot(plot_data, aes(condition, cor, fill=geno)) + geom_boxplot(outlier.shape = NA)
  }
  if(color_by_condition){
    p <- ggplot(plot_data, aes(geno, cor, fill=condition))
    if(remove_outliers){
      p <- p + geom_boxplot(outlier.shape = NA) + facet_grid(. ~ condition)
    }
    else{
      p <- p + geom_boxplot() + facet_grid(. ~ condition)
    }
    cc <- get_color_coding_dict()
    colScale <- scale_fill_manual(values=unlist(cc[plot_data$condition]))
    p <- p + colScale
    p <- p + theme(legend.position = 'none',
                   #axis.ticks.x = element_blank(),
                   #axis.text.x = element_blank(),
                   #axis.title.x = element_blank(),
                   #axis.title.y = element_blank(),
                   plot.title = element_text(size = 10)) +
                  ggtitle(title)
  }
  if(color_by_genotype){
    #colScale <- scale_fill_manual(values = unlist(cc[plot_data$conditions]))
    #p <- p + colScale
  }
  if(apply_paper_settings){
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + geom_jitter(aes(x=geno, y=cor, fill="none"), pch=21, size=0.5, alpha=0.3, dodge.width=0.4) + scale_color_manual(values="lightgrey")
  }
  if(!is.null(xlab)){
    p <- p + xlab(xlab)
  }
  if(!is.null(ylab)){
    p <- p + ylab(ylab)
  }
  if(!is.null(ylim)){
    p <- p + ylim(ylim)
  }
  return(p)
}


scatter_baseline_correlations <- function(geneA, tsv_loc, condition.1, condition.2, coexp_loc_prepend, coexp_append='_cor_coeqtlgenes.tsv'){
  # get significant genes per condition
  sigs_per_cond <- get_gene_lists(tsv_loc)
  # grab for both conditions and make unique
  cond1_genes <- sigs_per_cond[[condition.1]]
  cond2_genes <- sigs_per_cond[[condition.2]]
  cond1_genes <- gsub('-', '.', cond1_genes)
  cond2_genes <- gsub('-', '.', cond2_genes)
  interested_genes <- unique(c(cond1_genes, cond2_genes))
  # read the two coexpression matrices
  coexp_path_cond1 <- paste(coexp_loc_prepend, condition.1, coexp_append, sep='')
  coexp_condition1 <- read.table(coexp_path_cond1, header=T, row.names=1, sep = '\t')
  coexp_path_cond2 <- paste(coexp_loc_prepend, condition.2, coexp_append, sep='')
  coexp_condition2 <- read.table(coexp_path_cond2, header=T, row.names=1, sep = '\t')
  # we can only plot what is in both matrices
  interested_genes <- intersect(interested_genes, intersect(colnames(coexp_condition1), colnames(coexp_condition2)))
  # grab the values from each matrix
  plot_data <- data.frame(x=as.vector(unlist(coexp_condition1[geneA, interested_genes])), y=as.vector(unlist(coexp_condition2[geneA, interested_genes])))
  rownames(plot_data) <- interested_genes
  # order by x to make the plot easier on the eyes
  plot_data <- plot_data[order(plot_data$x), ]
  # add extra column to note where this was a coeqtl
  plot_data$coeqtl <- 'neither'
  if(nrow(plot_data[rownames(plot_data) %in% cond1_genes, ]) > 0){
    plot_data[rownames(plot_data) %in% cond1_genes, ]$coeqtl <- 'cond1'
  }
  if(nrow(plot_data[rownames(plot_data) %in% cond2_genes, ]) > 0){
    plot_data[rownames(plot_data) %in% cond2_genes, ]$coeqtl <- 'cond2'
  }
  if(nrow(plot_data[(rownames(plot_data) %in% cond1_genes & rownames(plot_data) %in% cond2_genes), ]) > 0){
    plot_data[rownames(plot_data) %in% cond1_genes & rownames(plot_data) %in% cond2_genes, ]$coeqtl <- 'both'
  }
  plot_data$coeqtl <- as.factor(plot_data$coeqtl)
  # actually plot now
  ggplot(plot_data, aes(x=x, y=y, color=coeqtl)) + geom_point() + geom_smooth(method=lm) + coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) + labs(y = condition.2, x=condition.1)
}


create_scatter_per_condition_combination <- function(geneA, tsv_loc, plot_save_loc, coexp_loc_prepend, coexp_append='_cor_coeqtlgenes.tsv', conditions=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  condition_plots <- list()
  for(condition in conditions){
    try({
      condition_plot <- scatter_baseline_correlations(geneA, tsv_loc, condition.1='UT', condition.2=condition, coexp_loc_prepend=coexp_loc_prepend, coexp_append=coexp_append)
      condition_plots[[condition]] <- condition_plot
    })
  }
  ggsave(plot_save_loc, arrangeGrob(grobs = condition_plots), width=12, height=8)
}

combined_sigs_to_file <- function(sigs_per_geneA, coefs_per_geneA, genes, combinations, output_loc){
  # check each gene
  for(gene in intersect(names(sigs_per_geneA), genes)){
    # we'll store the positive and negative direction genes
    pos_dir_geneBs <- c()
    neg_dir_geneBs <- c()
    # grab per condition
    for(condition in intersect(combinations, names(sigs_per_geneA[[gene]])) ){
      # grab these conditions
      geneBs <- sigs_per_geneA[[gene]][[condition]]
      # get the coef matrix
      coefs <- coefs_per_geneA[[gene]][[condition]]
      # get the mean coef
      coefs$mean <- apply(coefs, 1, mean)
      # grab the positive and negative genes
      coefs_pos <- rownames(coefs[!is.na(coefs$mean) &coefs$mean > 0, ])
      coefs_neg <- rownames(coefs[!is.na(coefs$mean) &coefs$mean < 0, ])
      # get the significant ones with the direction
      sigs_pos <- intersect(geneBs, coefs_pos)
      sigs_neg <- intersect(geneBs, coefs_neg)
      # add to the list
      pos_dir_geneBs <- c(pos_dir_geneBs, sigs_pos)
      neg_dir_geneBs <- c(neg_dir_geneBs, sigs_neg)
    }
    # make unique
    pos_dir_geneBs <- unique(pos_dir_geneBs)
    neg_dir_geneBs <- unique(neg_dir_geneBs)
    # paste combinations togeter
    combs_string <- paste(combinations, collapse = '')
    # set output locs
    pos_dir_out <- paste(output_loc, gene, '_', combs_string, '_pos.txt', sep = '')
    neg_dir_out <- paste(output_loc, gene, '_', combs_string, '_neg.txt', sep = '')
    # write the lists of genes
    write.table(data.frame(gene=pos_dir_geneBs), pos_dir_out, row.names = F, col.names = F, quote=F)
    write.table(data.frame(gene=neg_dir_geneBs), neg_dir_out, row.names = F, col.names = F, quote=F)
  }
}

create_re_coeqtl_matrices <- function(gene, cell_type, input_path_prepend, input_path_append, output_path_prepend, output_path_append, unstim_condition='UT', stim_conditions=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), allowed_missingness=0, na_to_zero=T){
  # build the input directory
  input_dir <- paste(input_path_prepend, gene, '_meta', input_path_append, sep = '')
  output_dir <- paste(output_path_prepend, gene, '_meta', output_path_append, sep = '')
  for(i in 1:2){
    cor_ut_loc <- paste(input_dir, 'correlationMatrices/', unstim_condition, '_', cell_type, '_correlation_matrix_', gene, '.', i, '.txt', sep = '')
    cor_ut <- read.table(cor_ut_loc, header = T, row.names = 1)
    # check per condition
    for(condition in stim_conditions){
      # read the correlation table
      cor_i_loc <- paste(input_dir, 'correlationMatrices/', condition, '_', cell_type, '_correlation_matrix_', gene, '.', i, '.txt', sep = '')
      cor_i <- read.table(cor_i_loc, header = T, row.names = 1)
      # check for missingness
      if(allowed_missingness <= 0){
        cor_i <- cor_i[apply(cor_i, 1, function(x){!any(is.na(x))}),]
      }
      # or if a number was given, then subset to the genes that have correlation
      else{
        cor_i <- cor_i[apply(cor_i, 1, function(x){sum(is.na(x))/length(x) <= allowed_missingness}),]
      }
      if(na_to_zero){
        cor_i[is.na(cor_i)] <- 0
      }
      # find common genes
      common_genes <- intersect(rownames(cor_ut), rownames(cor_i))
      # find common participants
      common_participants <- intersect(colnames(cor_ut), colnames(cor_i))
      # create new matrix
      diff_matrix <- matrix(, nrow=length(common_genes), ncol=length(common_participants), dimnames = list(common_genes, common_participants))
      # fill the matrix
      for(cor_gene in common_genes){
        for(participant in common_participants){
          ut_val <- cor_ut[cor_gene, participant]
          i_val <- cor_i[cor_gene, participant]
          diff_val <- diff(c(ut_val, i_val))
          diff_matrix[cor_gene, participant] <- diff_val
        }
      }
      output_full <- paste(output_dir, unstim_condition, '_vs_', condition, '_', cell_type, '_correlation_matrix_', gene, '.', i, '.txt', sep = '')
      write.table(diff_matrix, output_full, row.names=T, col.names=T, quote=F)
    }
  }
}

do_unpermuted_coeqtl_mapping_from_dir <- function(folder_loc, gene, cell_type, genotype_info, snp, cell_counts, condition_combinations=c('UT_vs_X3hCA', 'UT_vs_X24hCA', 'UT_vs_X3hMTB', 'UT_vs_X24hMTB', 'UT_vs_X3hPA', 'UT_vs_X24hPA'), sets=c(1,2)){
  # create overall p value matrix
  pval_df_full <- NULL
  # check each condition combination
  for(condition in condition_combinations){
    print(paste('starting condition', condition))
    dataset_per_set <- list()
    common_genes <- NULL
    # check each dataset
    for(set in sets){
      # load set
      set_loc <- paste(folder_loc, 'correlationMatrices/', condition, '_', cell_type, '_correlation_matrix_', gene, '.', set, '.txt', sep = '')
      set_i <- read.table(set_loc, header=T, row.names=1)
      # add to list
      dataset_per_set[[as.character(set)]] <- set_i
      # check the common genes
      if(is.null(common_genes)){
        common_genes <- rownames(set_i)
      }
      else{
        common_genes <- intersect(common_genes, rownames(set_i))
      }
    }
    # store the P values
    pvals <- c()
    # check the common genes
    for(common_gene in common_genes){
      # create the vectors to meta-analyse on
      betas <- c()
      stdes <- c()
      # check in each set
      for(set in names(dataset_per_set)){
        # grab the dataset
        set_i <- dataset_per_set[[set]]
        # get the common genotype data
        parts_common <- intersect(colnames(set_i), colnames(genotype_info))
        # grab the SNPs
        snps <- as.vector(unlist(genotype_info[snp, parts_common]))
        snps <- as.numeric(as.factor(snps))
        # grab the correlation difference
        cors <- as.vector(unlist(set_i[common_gene, parts_common]))
        # grab the cell type numbers
        cond.1 <- strsplit(condition, '_vs_')[[1]][1]
        cond.2 <- strsplit(condition, '_vs_')[[1]][2]
        cell.counts.1 <- as.vector(unlist(cell_counts[[cond.1]][[set]][parts_common]))
        cell.counts.2 <- as.vector(unlist(cell_counts[[cond.2]][[set]][parts_common]))
        cell.counts <- cell.counts.1 + cell.counts.2
        # do the model
        model.1 <- lm(formula = cors~snps, weights = sqrt(cell.counts))
        modelSummary.1 <- summary(model.1)
        modelCoeffs.1 <- modelSummary.1$coefficients
        beta.estimate.1 <- modelCoeffs.1[2, "Estimate"]
        std.error.1 <- modelCoeffs.1[2, "Std. Error"]
        # add to vectors
        betas <- c(betas, beta.estimate.1)
        stdes <- c(stdes, std.error.1)
      }
      # do the meta analysis
      metaAnalysis <- metagen(TE=betas, seTE = stdes, studlab = names(dataset_per_set))
      # grab the p value
      pval <- metaAnalysis$pval.random
      # add to pvals
      pvals <- c(pvals, pval)
    }
    # create pval df
    pval_df <- data.frame(p.value=pvals)
    rownames(pval_df) <- common_genes
    colnames(pval_df) <- condition
    # add to overarching df
    if(is.null(pval_df_full)){
      pval_df_full <- pval_df
    }
    else{
      pval_df_full <- merge(pval_df_full, pval_df, by=0, all=T)
      rownames(pval_df_full) <- pval_df_full$Row.names
      pval_df_full$Row.names <- NULL
    }
  }
  return(pval_df_full)
}

plot_coeqtl_numbers <- function(sigs_per_geneA, paper_style = F, table_instead = F){
  # init the table
  coeqtl_numbers <- NULL
  # check each gene
  for(gene in names(sigs_per_geneA)){
    # check each coeqtl
    for(condition in names(sigs_per_geneA[[gene]])){
      group <- '-'
      if(startsWith(condition, 'X3h')){
        group <- '3h'
      }
      else if(startsWith(condition, 'X24h')){
        group <- '24h'
      }
      else if(startsWith(condition, 'UT')){
        group <- 'UT'
      }
      nr_of_coeqtls <- length(sigs_per_geneA[[gene]][[condition]])
      # convert the label
      condition_pretty <- label_dict()[[condition]]

      row <- data.frame(gene=c(gene), condition=c(condition_pretty), timepoint=c(group), number=c(nr_of_coeqtls), stringsAsFactors = F)
      # add to df
      if(is.null(coeqtl_numbers)){
        coeqtl_numbers <- row
      }
      else{
        coeqtl_numbers <- rbind(coeqtl_numbers, row)
      }
    }
  }
  # set the order
  coeqtl_numbers$timepoint <- factor(coeqtl_numbers$timepoint, levels=c('UT', '3h', '24h'))
  # make the plot
  p <- ggplot(data=coeqtl_numbers, aes(x=timepoint, y=number, fill=condition)) + geom_bar(position='stack', stat='identity') + ggtitle(paste('')) + facet_wrap(. ~ gene, ncol=ceiling(sqrt(length(unique(coeqtl_numbers$gene)))))
  # create colour palette for conditions
  cc <- get_color_coding_dict()
  colScale <- scale_fill_manual(name = 'condition', values = unlist(cc[coeqtl_numbers$condition]))
  p <- p + colScale
  # add better labels
  p <- p + labs(x =  'Timepoint', y = 'number of co-expressionQTLs')
  # add the paper style (ugh)
  if(paper_style){
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if(table_instead){
    # we needed to supply the tables the figures were based on, that is implemented here
    return(coeqtl_numbers)
  }
  else{
    return(p)
  }
}

get_genes_per_pathway_per_condition <- function(output_path_prepend, output_path_append, gene, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', '24hCA', '3hMTB', '24hMTB', 'X3hPA', 'X24hPA'), sig_col='q.value.Bonferroni', sig_cutoff=0.05){
  all_lists <- list()
  for(cell_type in cell_types){
    list_ct <- list()
    for(condition in conditions){
      try({
      # paste the filepath together
      filepath <- paste(output_path_prepend, '', gene, '_', condition, '_meta_', cell_type, output_path_append, sep = '')
      # read the table
      pathways <- read.table(filepath, header=T, sep = '\t', quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
      # subset to significant stuff
      pathways <- pathways[pathways[[sig_col]] < sig_cutoff, ]
      # get the genes
      genes <- apply(pathways, 1, function(x){as.vector(unlist(strsplit(x['Hit.in.Query.List'], split = ',')))})
      # the pathways are the keys
      names(genes) <- pathways$Name
      # add to the list of condition
      list_ct[[condition]] <- genes
      })
    }
    # add to overarching list
    all_lists[[cell_type]] <- list_ct
  }
  return(all_lists)
}

plot_ng2018_rps26 <- function(ng2018_loc, coeqtl_output_loc, gene, cell_type, set, condition='UT', allele_dif=F, paper_style=F, verbose=T){
  # paste together the coeqtl output loc for p
  coeqtl_1m_p_loc <- paste(coeqtl_output_loc, gene, '_meta_', cell_type, '_p.tsv', sep = '')
  coeqtl_1m_r_loc <- paste(coeqtl_output_loc, gene, '_', cell_type, '_', condition, '_rs.tsv', sep = '')
  # read both tables
  coeqtl_1m_p <- read.table(coeqtl_1m_p_loc, header = T, sep = '\t', row.names = 1)
  coeqtl_1m_r <- read.table(coeqtl_1m_r_loc, header = T, sep = '\t', row.names = 1)
  # grab the r values
  sig_genes_1m <- rownames(coeqtl_1m_p[!is.na(coeqtl_1m_p[[condition]]) & coeqtl_1m_p[[condition]] < coeqtl_1m_p['significance_threshold', condition], , drop = F])
  # then grab those r values
  coeqtl_1m_r <- coeqtl_1m_r[sig_genes_1m, ]
  # change direction if other allele assessed
  if(allele_dif){
    coeqtl_1m_r <- coeqtl_1m_r * -1
  }
  # read the ng2018 file
  coeqtl_ng2018 <- read.table(ng2018_loc, header = T, sep = '\t', row.names = 1)
  # get the genes that are in both
  genes_both <- intersect(rownames(coeqtl_1m_r), rownames(coeqtl_ng2018))
  # create a dataframe with both
  coeqtl_r_both <- data.frame(ng2018=coeqtl_ng2018[genes_both, 'r'], m1=coeqtl_1m_r[genes_both, set], row.names = genes_both)
  # calculate the concordance
  concordance <- nrow(coeqtl_r_both[(coeqtl_r_both$ng2018 > 0 & coeqtl_r_both$m1 > 0) | (coeqtl_r_both$ng2018 < 0 & coeqtl_r_both$m1 < 0), ]) / nrow(coeqtl_r_both)
  # order by the first column
  coeqtl_r_both <- coeqtl_r_both[order(coeqtl_r_both$ng2018), ]
  # print verbose if requested
  if(verbose){
    print(paste('ng2018 coeqtls: ', nrow(coeqtl_ng2018), ', ',
                '1M coeqtls: ', nrow(coeqtl_1m_r), ', ',
                'coeqtls in both: ', nrow(coeqtl_r_both), ', ',
                'coeqtls concordant: ', nrow(coeqtl_r_both[(coeqtl_r_both$ng2018 > 0 & coeqtl_r_both$m1 > 0) | (coeqtl_r_both$ng2018 < 0 & coeqtl_r_both$m1 < 0), ]),
                sep = ''
          )
    )
  }
  # plot
  p <- ggplot(data=coeqtl_r_both, aes(x=ng2018, y=m1)) +
    xlim(c(-1, 1)) + ylim(c(-1, 1)) +
    xlab('r Nature 2018 CD4+ T') + ylab(paste('r 1M', cell_type, condition)) +
    ggtitle(paste('concordance of ', gene, ' co-eQTLs', sep = '')) +
    annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=0, fill="green", alpha=0.1) +
    annotate("rect", xmin=0, xmax=Inf, ymin=0, ymax=Inf, fill="green", alpha=0.1) +
    annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf, fill="red", alpha=0.1) +
    annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=0, fill="red", alpha=0.1) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    geom_vline(xintercept=0, linetype="dashed", color = "black") +
    geom_point() +
    annotate("rect", xmin=0.5, xmax=1.0, ymin=-1, ymax=-0.75, fill="white", alpha=1) +
    annotate("segment", x = 0.5, xend = 1, y = -0.75, yend = -0.75, colour = "black") +
    annotate("segment", x = 0.5, xend = 1, y = -1, yend = -1, colour = "black") +
    annotate("segment", x = 0.5, xend = 0.5, y = -0.75, yend = -1, colour = "black") +
    annotate("segment", x = 1, xend = 1, y = -0.75, yend = -1, colour = "black") +
    annotate("text", x = 0.75, y = -0.875, label = paste('', round(concordance*100, digits = 1), "% concordance", sep=''), fontface='bold')
  if(paper_style){
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  return(p)
}


write_coeqtl_boxplot_interaction_tables <- function(output_loc, tables_per_condition, make_anonymous=T){
  # we will write two tables
  correlation_table <- NULL
  interaction_table <- NULL
  # get the tables
  for(table_name in names(tables_per_condition)){
    # get the tables
    if('correlation' %in% names(tables_per_condition[[table_name]])){
      correlations <- tables_per_condition[[table_name]][['correlation']]
      correlations[['condition']] <- table_name
      # add to complete table
      if(is.null(correlation_table)){
        correlation_table <- correlations
      }
      else{
        correlation_table <- rbind(correlation_table, correlations)
      }
    }
    if('interaction' %in% names(tables_per_condition[[table_name]])){
      interactions <- tables_per_condition[[table_name]][['interaction']]
      interactions[['condition']] <- table_name
      # add to complete table
      if(is.null(interaction_table)){
        interaction_table <- interactions
      }
      else{
        interaction_table <- rbind(interaction_table, interactions)
      }
    }
  }
  if(make_anonymous){
    correlation_table[['participant']] <- NULL
    interaction_table[['sample.name']] <- NULL
  }
  # write results
  output_loc_correlations <- paste(output_loc, 'correlations.tsv', sep = '')
  write.table(correlation_table, output_loc_correlations, sep = '\t', col.names = T, row.names = F)
  output_loc_interactions <- paste(output_loc, 'interactions.tsv', sep = '')
  write.table(interaction_table, output_loc_interactions, sep = '\t', col.names = T, row.names = F)
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- "lightgrey"
  color_coding[["3hCA"]] <- "darkolivegreen2"
  color_coding[["24hCA"]] <- "forestgreen"
  color_coding[["3hMTB"]] <- "lightskyblue"
  color_coding[["24hMTB"]] <- "deepskyblue3"
  color_coding[["3hPA"]] <- "sandybrown"
  color_coding[["24hPA"]] <- "darkorange1"
  color_coding[["X3hCA"]] <- "darkolivegreen2"
  color_coding[["X24hCA"]] <- "forestgreen"
  color_coding[["X3hMTB"]] <- "lightskyblue"
  color_coding[["X24hMTB"]] <- "deepskyblue3"
  color_coding[["X3hPA"]] <- "sandybrown"
  color_coding[["X24hPA"]] <- "darkorange1"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  color_coding[["CD4/CD8+ T"]] <- "#0B6799"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["hemapoietic stem"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["megakaryocyte"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["plasma B"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}

get_color_coding_dict_darker <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- "lightgrey"
  color_coding[["3hCA"]] <- "limegreen"
  color_coding[["24hCA"]] <- "darkgreen"
  color_coding[["3hMTB"]] <- "royalblue"
  color_coding[["24hMTB"]] <- "darkblue"
  color_coding[["3hPA"]] <- "chocolate"
  color_coding[["24hPA"]] <- "brown"
  color_coding[["X3hCA"]] <- "darkolivegreen2"
  color_coding[["X24hCA"]] <- "forestgreen"
  color_coding[["X3hMTB"]] <- "lightskyblue"
  color_coding[["X24hMTB"]] <- "deepskyblue3"
  color_coding[["X3hPA"]] <- "sandybrown"
  color_coding[["X24hPA"]] <- "darkorange1"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  color_coding[["CD4/CD8+ T"]] <- "#0B6799"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["hemapoietic stem"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["megakaryocyte"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["plasma B"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}

get_text_colour_dict <- function(){
  # set the condition colors
  color_coding <- list()
  # set the cell type colors
  color_coding[["Bulk"]] <- "white"
  color_coding[["CD4T"]] <- "white"
  color_coding[["CD8T"]] <- "black"
  color_coding[["monocyte"]] <- "black"
  color_coding[["NK"]] <- "black"
  color_coding[["B"]] <- "black"
  color_coding[["DC"]] <- "white"
  color_coding[["CD4+ T"]] <- "white"
  color_coding[["CD8+ T"]] <- "black"
  # other cell type colors
  color_coding[["HSPC"]] <- "black"
  color_coding[["hemapoietic stem"]] <- "black"
  color_coding[["platelet"]] <- "black"
  color_coding[["megakaryocyte"]] <- "black"
  color_coding[["plasmablast"]] <- "black"
  color_coding[["plasma B"]] <- "black"
  color_coding[["other T"]] <- "black"
  return(color_coding)
}

label_dict <- function(){
  label_dict <- list()
  label_dict[["UT"]] <- "UT"
  label_dict[["X3hCA"]] <- "3hCA"
  label_dict[["X24hCA"]] <- "24hCA"
  label_dict[["X3hMTB"]] <- "3hMTB"
  label_dict[["X24hMTB"]] <- "24hMTB"
  label_dict[["X3hPA"]] <- "3hPA"
  label_dict[["X24hPA"]] <- "24hPA"
  label_dict[["3hCA"]] <- "3hCA"
  label_dict[["24hCA"]] <- "24hCA"
  label_dict[["3hMTB"]] <- "3hMTB"
  label_dict[["24hMTB"]] <- "24hMTB"
  label_dict[["3hPA"]] <- "3hPA"
  label_dict[["24hPA"]] <- "24hPA"
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["hemapoietic stem"]] <- "hemapoietic stem"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["plasma B"]] <- "plasma B"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["megakaryocyte"]] <- "megakaryocyte"
  label_dict[["T_other"]] <- "other T"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  return(label_dict)
}


# location of the coeqtl output
coeqtl_out_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/'
# the coeqtl 'geneA' genes we looked at
coeqtl_genes <- c('HLA-DQA1', 'TMEM176B','TMEM176A', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26', 'HLA-DQA2')
# we need to paste together the whole thing
prepend <- coeqtl_out_loc
append <- '_meta_monocyte_p.tsv'
# get the sigs per geneA
sigs_per_geneA <- get_gene_list_geneAs(prepend, append, coeqtl_genes)
# set the location of where to store the significant genes
significant_coeqtl_genes_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/significant_genes_20210208_numeric_snps/'
# write these genes
significant_coeqtl_genes_to_file(prepend, append, significant_coeqtl_genes_loc, '_meta_monocyte', coeqtl_genes)
# we'll plot the numbers
plot_coeqtl_numbers(sigs_per_geneA, paper_style = T, table_instead = F)
# we want the table as well
coeqtl_numbers_numbers <- plot_coeqtl_numbers(sigs_per_geneA, paper_style = T, table_instead = T)

# a look at the pathways as well
plot_pathway_all(output_path_prepend = '/data/scRNA/eQTL_mapping/coexpressionQTLs/pathways/pathway_out_20210208_numeric_snps/coeqtls_', output_path_append = '_pathways.txt', gene = 'CLEC12A', toppfun = T, top_so_many = 5)
plot_pathway_all(output_path_prepend = '/data/scRNA/eQTL_mapping/coexpressionQTLs/pathways/pathway_out_20210208_numeric_snps/coeqtls_', output_path_append = '_pathways.txt', gene = 'CLEC12A', toppfun = T, top_so_many = 5)




# get the genotype data
vcf <- fread('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/genotypes/LL_trityper_plink_converted.vcf.gz')
genotypes_all <- as.data.frame(vcf[, 10:ncol(vcf)])
rownames(genotypes_all) <- vcf$ID

# mapping for the probes that belong to the genes
snp_probe_mapping_location <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/snp_gene_mapping_20201113.tsv'
# get the mapping of the probe to the cis SNP
snp_probe_mapping <- read.table(snp_probe_mapping_location, sep = '\t', header=T, stringsAsFactors = F)

# location of coeqtls output and correlation matrices
coeqtl_out_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/'
# prep and append we will use for now
prepend <- 'output_'
append <- '_meta_mono_missingness05replacena100permzerogenebnumeric_1/'
# location of plots
coeqtl_plot_loc <- paste(coeqtl_out_loc, 'plots/', sep = '')
# geneAs to plot
coeqtl_genes <- c('SSU72','NDUFA12','NMB','PLGRKT','SEPHS2','BTN3A2','PGD','TNFAIP6','HLA-B','ZFAND2A','HEBP1','CTSC','TMEM109','NUCB2','HIP1','AP2S1','CD52','PPID','RPS26','TMEM176B','ERAP2','HLA-DQA2','TMEM176A','CLEC12A','MAP3K7CL','BATF3','MRPL54','LILRA3','NAAA','PRKCB','SMDT1','LGALS9','KIAA1598','UBE2D1','SCO2','DNAJC15','NDUFA10','NAA38','HLA-DQA1','ROGDI','RBP7','SDCCAG8','CFD','GPX1','PRDX2','C6orf48','RBBP8','IQGAP2','PTK2B','SMAP1')

# make a plot of CLEC12A interactions
plot_interaction(seurat_object=v3,
                 gene1='CLEC12A',
                 gene2='PML', 
                 genotype=genotypes_all['rs12230244', ], 
                 snp.name='rs12230244', 
                 version_chem='v3', 
                 output_loc='/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/', p.value = NULL, r.value = NULL, sign.cutoff = NULL, use_SCT=T, sct_data=F, conditions=c('UT', 'X3hCA', 'X24hCA'), set_axis_by_max=T, height=5, width=5, extention='.png', remove_prepend=NULL, snp_rename=NULL, xlim=NULL, ylim=NULL)


# load a correlation matrix
cor_data_X24hCA <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_CLEC12A_meta_mono_missingness05replacena100permzerogenebnumeric_1/correlationMatrices/X24hCA_monocyte_correlation_matrix_CLEC12A.1.txt')
cor_data_X3hCA <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_CLEC12A_meta_mono_missingness05replacena100permzerogenebnumeric_1/correlationMatrices/X3hCA_monocyte_correlation_matrix_CLEC12A.1.txt')
cor_data_UT <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_CLEC12A_meta_mono_missingness05replacena100permzerogenebnumeric_1/correlationMatrices/UT_monocyte_correlation_matrix_CLEC12A.2.txt')
# print the boxplot
X24hCA_CLEC12A_tables <- plot_boxplot_and_gt_interactions(genotype_data=genotypes_all['rs12230244', ],
                                 mapping_folder='/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_CLEC12A_meta_mono_missingness05replacena100permzerogenebnumeric_1/', 
                                 dataset='1', 
                                 seurat_object=v2, 
                                 gene_a='CLEC12A', 
                                 gene_b='PML', 
                                 snp='rs12230244', 
                                 monniker='_meta_', 
                                 cell_type='monocyte', 
                                 condition='X24hCA', 
                                 na_to_zero=T, 
                                 to_numeric=T, 
                                 snp_rename=NULL, 
                                 p.value=NULL, 
                                 r.value=NULL, 
                                 ylim=NULL, 
                                 xlim=NULL, 
                                 ylims=NULL, 
                                 use_SCT=T, 
                                 sct_data=F, 
                                 cor_table=cor_data_X24hCA, 
                                 plot_points=F, 
                                 violin=F,
                                 return_table=T)
X3hCA_CLEC12A_tables <- plot_boxplot_and_gt_interactions(genotype_data=genotypes_all['rs12230244', ], mapping_folder='/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_CLEC12A_meta_mono_missingness05replacena100permzerogenebnumeric_1/', dataset='1', seurat_object=v2, gene_a='CLEC12A', gene_b='PML', snp='rs12230244', monniker='_meta_', cell_type='monocyte', condition='X3hCA', na_to_zero=T, to_numeric=T, snp_rename=NULL, p.value=NULL, r.value=NULL, ylim=NULL, xlim=NULL, ylims=NULL, use_SCT=T, sct_data=F, cor_table=cor_data_X3hCA, plot_points=F, violin=F, return_table=T) 
UT_CLEC12A_tables <- plot_boxplot_and_gt_interactions(genotype_data=genotypes_all['rs12230244', ], mapping_folder='/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_CLEC12A_meta_mono_missingness05replacena100permzerogenebnumeric_1/', dataset='2', seurat_object=v3, gene_a='CLEC12A', gene_b='PML', snp='rs12230244', monniker='_meta_', cell_type='monocyte', condition='UT', na_to_zero=T, to_numeric=T, snp_rename=NULL, p.value=NULL, r.value=NULL, ylim=NULL, xlim=NULL, ylims=NULL, use_SCT=T, sct_data=F, cor_table=cor_data_UT, plot_points=F, violin=F, return_table=T) 
write_coeqtl_boxplot_interaction_tables(output_loc = '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/figure_tables/figure_4b_', list(UT = UT_CLEC12A_tables, X3hCA = X3hCA_CLEC12A_tables, X24hCA = X24hCA_CLEC12A_tables))

# plot CLEC12A
plot_pathway_all(output_path_prepend = '/data/scRNA/eQTL_mapping/coexpressionQTLs/pathways/pathway_out_20210208_numeric_snps/coeqtls_', output_path_append = '_pathways.txt', gene = 'CLEC12A', toppfun = T, top_so_many = 5)
# print the tables used for the pathway heatmaps
write_pathway_plot_tables(output_path_prepend = '/data/scRNA/eQTL_mapping/coexpressionQTLs/pathways/pathway_out_20210208_numeric_snps/coeqtls_', output_path_append = '_pathways.txt', output_loc='~/Desktop/figureS6/', c('NDUFA12','CTSC','RPS26' ,'TMEM176B','CLEC12A','DNAJC15','HLA-DQA1'), cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), top_so_many=5, reactome = F, toppfun = T)
  

# do the top hit plotting
plot_top_hit_per_condition(genotype_data=genotypes_all, mappings_folder=coeqtl_out_loc, mapping_folder_prepend=prepend, mapping_folder_append='_mono_missingness05replacena100permzerogeneb_1/', plot_output_loc=coeqtl_plot_loc, genes=coeqtl_genes, snp_probe_mapping=snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), na_to_zero=T)
# top interactions as well
plot_top_hit_per_interaction(genotype_data=genotypes_all, mappings_folder=coeqtl_out_loc, mapping_folder_prepend=prepend, mapping_folder_append='_mono_missingness05replacena100permzerogenebnumeric_1/', plot_output_loc=paste(coeqtl_plot_loc, '_num_', sep=''), genes=ff_coeqtl_genes_less, v2_object=v2_mono, v3_object=v3_mono, snp_probe_mapping=snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))

# check the Rs
get_r_values(input_path_prepend='/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', input_path_append='_mono_missingness05replacena100permzerogenebnumeric_1/', gene='RPS26', snp='rs1131017', snps=genotypes_all, cell_type='monocyte', to_numeric=T, sig_only = F, cell_counts = cell_counts)
# write the Rs
write_r_values_per_gene_and_condition(input_path_prepend='/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', input_path_append='_mono_missingness05replacena100permzerogenebnumeric_1/', genes=coeqtl_genes, snp_probe_mapping=snp_probe_mapping, snps=genotypes_all, cell_type='monocyte', to_numeric=T)
# write the R plots
write_r_plots_per_gene_and_condition(input_path_prepend='/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', input_path_append='_mono_missingness05replacena100permzerogenebnumeric_1/', genes=coeqtl_genes, plot_output_loc=coeqtl_plot_loc, cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
