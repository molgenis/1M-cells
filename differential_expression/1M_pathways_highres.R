############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_MAST_highres.R
# Function: read enrichr pathway analysis output and compare high resolution vs low resolution
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
require("heatmap.plus")
library(RColorBrewer)
library(VennDiagram)


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

get_pathway_table <- function(pathway_output_loc, append='_sig_up_pathways.txt', sig_val_to_use = 'q.value.FDR', sig_val_to_store=NULL, significance_cutoff=0.05, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB'), use_ranking=F){
  # put all results in a list
  pathways_analysis <- list()
  # put all results in a shared DF
  pathway_df <- NULL
  # check each cell type
  for(cell_type in cell_types){
    # check each stim
    for(stim in stims){
      try({
        print(paste(cell_type, stim, sep = ' '))
        # paste the filepath together
        #filepath <- paste(pathway_output_loc, cell_type, 'UT',stim,'_sig_pathways.txt', sep = '')
        filepath <- paste(pathway_output_loc, cell_type, 'UT',stim,append, sep = '')
        # read the file
        pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
        # create column name
        newcolname <- paste(cell_type, 'UT', stim, sep = '')
        # get the log2 of the significance value
        pathways <- pathways[pathways[[sig_val_to_use]] < significance_cutoff, ]
        if(use_ranking){
          pathways[[newcolname]] <- as.numeric(rownames(pathways))
        }
        else{
          # get the significance value to put in the table, the default is the same as the significance value to 
          sig_val_to_for_storing <- sig_val_to_use
          if(!is.null(sig_val_to_store)){
            sig_val_to_for_storing <- sig_val_to_store
          }
          pathways[[newcolname]] <- log10(pathways[[sig_val_to_for_storing]])
          #pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
          # log transform the P values
          #pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
        }
        pathways$id_name <- paste(pathways$ID, pathways$Name, sep = '_')
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
          #pathway_df[[newcolname]] <- pathways[[newcolname]][match(pathway_df$Name, pathways$Name)]
          #pathway_df <- left_join(pathway_df, pathways)
          
        }
      })
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


get_pathway_tables <- function(pathway_output_loc, append='_sig_up_pathways.txt', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
  # put all results in a list
  pathways_analysis <- list()
  # put all results in a shared DF
  pathway_df <- NULL
  # check each cell type
  for(cell_type in cell_types){
    # check each stim
    for(stim in stims){
      try({
        print(paste(cell_type, stim, sep = ' '))
        # paste the filepath together
        filepath <- paste(pathway_output_loc, cell_type, 'UT',stim, append, sep = '')
        # read the file
        pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
        # put in in the list
        pathways_analysis[[paste(cell_type, 'UT', stim, sep='')]] <- pathways
      })
    }
  }
  return(pathways_analysis)
}

filter_pathway_on_category <- function(pathway_table, filtered_names){
  pathway_table_filtered <- pathway_table[pathway_table$Name %in% filtered_names, ]
  return(pathway_table_filtered)
}

filter_pathways_on_category <- function(pathway_list, filtered_names){
  # new list for filtered pathway
  filtered_pathway_list <- list()
  # check each pathway df in the result list
  for(key in names(pathway_list)){
    pathway_df <- pathway_list[[key]]
    filtered_pathway_df <- filter_pathway_on_category(pathway_df, filtered_names)
    # put in the filtered list
    filtered_pathway_list[[key]] <- filtered_pathway_df
  }
  return(filtered_pathway_list)
}

get_filtered_pathway_names <- function(pathway_table, relation_table, starting_id){
  # get all of the children of the starting ID
  all_children <- get_children(relation_table, starting_id)
  # get the names of the pathways that are children
  pathway_names <- as.character(pathway_table[pathway_table$V1 %in% all_children, ]$V2)
  return(pathway_names)
}

get_children <- function(relation_table, starting_id){
  # get all of the children of the starting ID
  children <- as.character(relation_table[relation_table$V1 == starting_id, 'V2'])
  # these children are all family
  family <- children
  # see if there were any children
  if(length(children) > 0){
    # if there were children, we need to get their children as well
    for(child in children){
      # get the grandchildren and add these to the family
      grand_children <- get_children(relation_table, child)
      family <- c(family, grand_children)
    }
  }
  return(family)
}

get_top_pathways <- function(pathway_table, nr_of_top_genes, is_ranked=F){
  # init pathways list
  pathways <- c()
  # go through the columns
  for(col in colnames(pathway_table)){
    # order by that column
    ordered <- pathway_table[order(pathway_table[[col]], decreasing = T), ]
    if(is_ranked){
      ordered <- pathway_table[order(pathway_table[[col]], decreasing = F), ]
    }
    # get those top ones
    top_col <- rownames(ordered)[1:nr_of_top_genes]
    pathways <- c(pathways, top_col)
  }
  # limit to those top pathways now
  pathway_table_smaller <- pathway_table[rownames(pathway_table) %in% pathways, ]
  return(pathway_table_smaller)
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
    summed_rank_3h <- apply(pathway_table[, timepoints_3h], 1, sum)
    most_shared <- c(most_shared, rownames(pathway_table[order(summed_rank_3h), ])[1:top_so_many])
    
    # get the sum of the 24h conditions
    timepoints_24h <- colnames(pathway_table)[grep('24h', colnames(pathway_table))]
    summed_rank_24h <- apply(pathway_table[, timepoints_24h], 1, sum)
    most_shared <- c(most_shared, rownames(pathway_table[order(summed_rank_24h), ])[1:top_so_many])
    
    # make unique of course
    most_shared <- unique(most_shared)
    
  }
  return(most_shared)
}

get_most_varied_pathways <- function(pathway_table, top_so_many=10, use_sd_method=F, disregard_3h_vs_24h=F){
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
    most_varied <- c()
    # because 3h has sharing and 24h has sharing, looking overall, we might see 3h vs 24h effects, we can choose to ignore this
    if(disregard_3h_vs_24h == F){
      most_varied <- get_most_varying_from_df(pathway_table, top_so_many = top_so_many)
    }
    # most varied in 3h and most varied in 24h
    most_varied <- c(most_varied, get_most_varying_from_df(pathway_table[, timepoints_3h], top_so_many = top_so_many))
    most_varied <- c(most_varied, get_most_varying_from_df(pathway_table[, timepoints_24h], top_so_many = top_so_many))
  }
  else{
    # large difference between 3h and 24h
    mean_3h_ranks <- apply(pathway_table[, timepoints_3h], 1, mean)
    mean_24h_ranks <- apply(pathway_table[, timepoints_24h], 1, mean)
    # absolute difference
    mean_rank_diff <- abs(mean_3h_ranks - mean_24h_ranks)
    # grab by varied over abs mean difference
    most_varied <- rownames(pathway_table[order(mean_rank_diff, decreasing = T), ])[1:top_so_many]
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
  most_varied <- unique(most_varied)
  return(most_varied)
}


get_most_varying_from_df <- function(dataframe, top_so_many=10, dont_get_least_varying=T, minimal_difference=NULL){
  # sometimes we want to set a cutoff of a minimal value
  if(!is.null(minimal_difference)){
    max_difference <- apply(dataframe, 1, function(x){
      # extract the minimum and maximum variable
      min_val <- min(as.vector(unlist(x)))
      max_val <- max(as.vector(unlist(x)))
      # check the difference
      diffr_max <- diff(c(min_val, max_val))
      # turn it into an absolute value (easier for the next step)
      diffr_abs <- abs(diffr_max)
      return(diffr_abs)
    })
    # can now filtering on this
    dataframe <- dataframe[max_difference >= minimal_difference, ]
  }
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


get_number_of_hits_pathways <- function(pathway_output_loc, append='_sig_up_pathways.txt', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
  # get pathways per cell type and condition combination
  pathways_per_ct_tp <- get_pathway_tables(pathway_output_loc=pathway_output_loc, append=append, cell_types=cell_types, stims=stims)
  # store the percentage somewhere
  per_cell_type_pathways <- list()
  # check each cell type
  for(key in names(pathways_per_ct_tp)){
    # grab the data
    pathways <- pathways_per_ct_tp[[key]]
    # get the fractions
    fractions <- as.list(pathways$Hit.Count.in.Query.List/pathways$Hit.Count.in.Genome)
    # the pathway names will be the keys
    keys <- pathways$Name
    # if the names are not unique, we'll have to add the ID as well
    if(length(keys) != length(unique(keys))){
      keys <- paste(pathways$Name, pathways$ID, sep = '_')
    }
    # set the keys on the list
    names(fractions) <- keys
    # now put these into the list
    per_cell_type_pathways[[key]] <- fractions
  }
  return(per_cell_type_pathways)
}


pathway_hits_to_table <- function(pathway_per_cond_and_ct){
  # there will be an overarching table
  all_pathways <- NULL
  # we will check each key, which is a combination of the cell type and condition combination
  for(key in names(pathway_per_cond_and_ct)){
    # get the list
    this_combination <- pathway_per_cond_and_ct[[key]]
    # put into data table
    this_pathways_table <- data.table(pathways=names(this_combination), values=as.vector(unlist(this_combination)))
    # set the column name to contain the key, it's what we use in most other places
    colnames(this_pathways_table) <- c('pathway', key)
    # merge together if not the first entry
    if(is.null(all_pathways)){
      all_pathways <- this_pathways_table
    }
    else{
      all_pathways <- merge(all_pathways, this_pathways_table, by='pathway', all = T)
    }
  }
  # turn into dataframe object we usually use
  all_pathways <- data.frame(all_pathways)
  # we usually have the name of the pathway as rownames
  rownames(all_pathways) <- all_pathways$pathway
  # remove the column with pathways, as it is the rowname now
  all_pathways$pathway <- NULL
  return(data.frame(all_pathways))
}


create_dotplot_per_condition_and_celltypes <- function(pathway_p_or_ranking_table, pathway_fraction_table, stims=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), normal='UT',
                                                       cell_type_sets=list('CD4T' = c('Naive CD4+ T', 'Memory CD4+ T'), 'DC' = c('pDC', 'mDC'), 'CD8T' = c('Naive CD8+ T', 'Memory CD8+ T'), 'NK' = c('NKbright', 'NKdim'), 'monocyte' = c('cMono', 'ncMono')), use_most_varied=F, top_so_many=10, minimal_difference=NULL, order_by_diff=F){
  # we'll store per stimulation
  plot_frames_per_stim <- get_plot_frames_per_condition_and_celltypes(pathway_p_or_ranking_table, pathway_fraction_table, stims=stims, normal=normal, cell_type_sets=cell_type_sets, use_most_varied=use_most_varied, top_so_many=top_so_many, minimal_difference = minimal_difference, order_by_diff = order_by_diff)
  # create list to store the lists of ggplot plots
  plots_per_stim <- list()
  # check each stim
  for(stim in stims){
    # create a list to store the ggplots
    plots_per_ct <- list()
    # check each cell type
    for(cell_type in names(cell_type_sets)){
      # fetch the plot data
      plot_data_stim_ct <- plot_frames_per_stim[[stim]][[cell_type]]
      # create the plot
      p <- ggplot(data=plot_data_stim_ct, mapping = aes(x=cell_type, y=pathway, color=p_or_rank, size=fraction)) +
        geom_point() +
        scale_color_gradient(low='red', high='yellow') + 
        theme(text = element_text(size = 15)) + 
        xlab('cell type') + 
        ylab('Pathways') + 
        labs(color = '-log10(p)')
      # put into list
      plots_per_ct[[cell_type]] <-p
    }
    # put in list
    plots_per_stim[[stim]] <- plots_per_ct
  }
  return(plots_per_stim)
}


get_plot_frames_per_condition_and_celltypes <- function(pathway_p_or_ranking_table, pathway_fraction_table, stims=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), normal='UT',
                                                        cell_type_sets=list('CD4T' = c('Naive CD4+ T', 'Memory CD4+ T'), 'DC' = c('pDC', 'mDC'), 'CD8T' = c('Naive CD8+ T', 'Memory CD8+ T'), 'NK' = c('NKbright', 'NKdim'), 'monocyte' = c('cMono', 'ncMono')), use_most_varied=F, top_so_many=10, minimal_difference=NULL,order_by_diff=F){
  # we'll store per stimulation
  plot_frames_per_stim <- list()
  # check each stim
  for(stim in stims){
    # init the list we will use
    plot_frame_ct <- list()
    # check each major cell type
    for(cell_type_major in names(cell_type_sets)){
      print(paste('plotframe:', stim, cell_type_major, sep = ''))
      # grab the subcell types
      sub_cell_type <- cell_type_sets[[cell_type_major]]
      # init the table
      p_and_frac_tbl_ct <- NULL
      # create columns of the sub cell types we will use
      append_colnames <- paste(normal, stim, sep = '')
      ct_and_tp_cols <- paste(sub_cell_type, append_colnames, sep = '')
      # filter by existing columns
      ct_and_tp_cols <- intersect(ct_and_tp_cols, colnames(pathway_p_or_ranking_table))
      # check each subceltype
      for(i in 1:length(ct_and_tp_cols)){
        # paste the column together (they are always combined with the normal condition)
        ct_and_tp_col <- ct_and_tp_cols[i]
        # grab ranking or p
        pathway_p_or_ranking_ct_tp <- pathway_p_or_ranking_table[, c(ct_and_tp_col), drop = F]
        # set a name for the column with the p or rank
        colnames(pathway_p_or_ranking_ct_tp) <- c('p_or_rank')
        # the fractions as well
        pathway_fraction_table_ct_tp <- pathway_fraction_table[, c(ct_and_tp_col), drop = F]
        # set a name for the column fraction
        colnames(pathway_fraction_table_ct_tp) <- c('fraction')
        # join them by the rownames
        pathway_p_or_ranking_and_fraction <- merge(pathway_fraction_table_ct_tp, pathway_p_or_ranking_ct_tp, by=0, all=T)
        # when merging on rownames with data.frame objects, the row names are added under Row.names. This is rediculous, so let's fix that
        pathway_p_or_ranking_and_fraction$pathway <- as.vector(unlist(pathway_p_or_ranking_and_fraction$Row.names))
        pathway_p_or_ranking_and_fraction$Row.names <- NULL
        # set the cell type
        pathway_p_or_ranking_and_fraction$cell_type <- sub_cell_type[i]
        # add to the existing table
        if(is.null(p_and_frac_tbl_ct)){
          p_and_frac_tbl_ct <- pathway_p_or_ranking_and_fraction
        }
        else{
          p_and_frac_tbl_ct <- rbind(p_and_frac_tbl_ct, pathway_p_or_ranking_and_fraction)
        }
      }
      # subset to most varied if requested
      if(use_most_varied & length(ct_and_tp_cols) > 1){
        # we will do most varied for the sub cell types, so that's those
        vars_to_check <- pathway_p_or_ranking_table[, c(ct_and_tp_cols)]
        # convert zero to max if ranked
        if(sum(vars_to_check) > 0){
          vars_to_check[vars_to_check==0] <- max(vars_to_check)
        }
        # sign flip if p value
        else{
          vars_to_check <- vars_to_check*-1
        }
        # remove completely empty rows
        vars_to_check <- vars_to_check[rowSums(vars_to_check) != 0, ]
        # get most variation
        most_varied <- get_most_varying_from_df(vars_to_check, top_so_many = top_so_many, minimal_difference = minimal_difference)
        # subset to these most varied
        p_and_frac_tbl_ct <- p_and_frac_tbl_ct[p_and_frac_tbl_ct$pathway %in% most_varied, ]
      }
      else if(use_most_varied & length(ct_and_tp_cols) == 1){
        # sort by the variable
        p_and_frac_tbl_ct <- p_and_frac_tbl_ct[order(p_and_frac_tbl_ct[[1]]), ,drop = F]
        # confine to smallest values if there are more that the top x we requestedd
        if(nrow(p_and_frac_tbl_ct) > top_so_many){
          p_and_frac_tbl_ct <- p_and_frac_tbl_ct[1:top_so_many, ,drop = F]
        }
      }
      if(order_by_diff & length(ct_and_tp_cols) > 1){
        # check which pathways we have left
        pathways_left <- p_and_frac_tbl_ct$pathway
        # now get back the original values
        vars_to_check <- pathway_p_or_ranking_table[, c(ct_and_tp_cols)]
        # remove the ones we already removed
        vars_to_check[rownames(vars_to_check) %in% pathways_left, ]
        # get the differences in the row
        max_difference <- apply(vars_to_check, 1, function(x){
          # extract the minimum and maximum variable
          min_val <- min(as.vector(unlist(x)))
          max_val <- max(as.vector(unlist(x)))
          # check the difference
          diffr_max <- diff(c(min_val, max_val))
          # turn it into an absolute value (easier for the next step)
          diffr_abs <- abs(diffr_max)
          return(diffr_abs)
        })
        # use these to order
        val_order <- order(max_difference)
        # get the pathways for those
        pathway_ordered <- rownames(vars_to_check[val_order, ])
        # now turn into factor
        p_and_frac_tbl_ct$pathway <- factor(as.character(p_and_frac_tbl_ct$pathway), levels=as.character(pathway_ordered), ordered=T)
      }
      # add to list of each cell type
      plot_frame_ct[[cell_type_major]] <- p_and_frac_tbl_ct
    }
    plot_frames_per_stim[[stim]] <- plot_frame_ct
  }
  return(plot_frames_per_stim)
}


print_listed_ggplot_objects <- function(list_in_lists, output_loc, extention='pdf', width = 10, height = 10, xlab=NULL, labs=NULL){
  # check each element
  for(key in names(list_in_lists)){
    # get the actual element
    element <- list_in_lists[[key]]
    # paste together the output loc
    output_loc_nest <- paste(output_loc, key, '.', sep = '')
    # check if we are at the last leaves, or need to nest further
    if(is.ggplot(element)){
      # finish the path
      output_loc_final <- paste(output_loc_nest, extention, sep = '')
      # get the object
      p <- element
      # check for extra parameters
      if(!is.null(xlab)){
        p <- p + xlab
      }
      if(!is.null(labs)){
        p <- p + labs
      }
      print(output_loc_final)
      # save the object
      ggsave(output_loc_final, p, width = width, height = height)
    }
    else{
      # if we are nested further, we can recursively use this function again
      print_listed_ggplot_objects(element, output_loc_nest, extention = extention, width = width, height = height, xlab = xlab, labs = labs)
    }
  }
}

write_dataframe_per_condition_and_celltypes <- function(pathway_p_or_ranking_table, pathway_fraction_table, output_loc, stims=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), normal='UT', cell_type_sets=list('CD4T' = c('Naive CD4+ T', 'Memory CD4+ T'), 'DC' = c('pDC', 'mDC'), 'CD8T' = c('Naive CD8+ T', 'Memory CD8+ T'), 'NK' = c('NKbright', 'NKdim'), 'monocyte' = c('cMono', 'ncMono')), use_most_varied=F, top_so_many=10, minimal_difference=NULL, order_by_diff=F){
  # we'll store per stimulation
  plot_frames_per_stim <- get_plot_frames_per_condition_and_celltypes(pathway_p_or_ranking_table, pathway_fraction_table, stims=stims, normal=normal, cell_type_sets=cell_type_sets, use_most_varied=use_most_varied, top_so_many=top_so_many, minimal_difference = minimal_difference, order_by_diff = order_by_diff)
  # check each stim
  for(stim in names(plot_frames_per_stim)){
    # check each cell type
    for(celltype in names(plot_frames_per_stim[[stim]])){
      # get those variables
      df_stim_ct <- plot_frames_per_stim[[stim]][[celltype]]
      # paste together an output location
      output_loc_full <- paste(output_loc, stim, '_', celltype, '.tsv', sep = '')
      # write the table
      write.table(df_stim_ct, output_loc_full, sep = '\t', row.names = F, col.names = T)
    }
  }
}

get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[['UT']] <- 'lightgrey'
  color_coding[["3hCA"]] <- "darkolivegreen2"
  color_coding[["24hCA"]] <- "forestgreen"
  color_coding[["3hMTB"]] <- "lightskyblue"
  color_coding[["24hMTB"]] <- "deepskyblue3"
  color_coding[["3hPA"]] <- "sandybrown"
  color_coding[["24hPA"]] <- "darkorange1"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  return(color_coding)
}


####################
# Main Code        #
####################

# where we write and read
write_partition <- 'tmp01'
read_partition <- 'tmp01'

# this is the reactome ID for the immune system
immune_system_reactome_id <- 'R-HSA-168256'
# load the pathways
pathways <- read.table('/data/scRNA/pathways/ReactomePathways.tsv', sep='\t')
# subset to just human to speed up the search
pathways <- pathways[pathways$V3 == 'Homo sapiens', ]
# load the pathway mapping
pathway_mappings <- read.table('/data/scRNA/pathways/ReactomePathwaysRelation.tsv', sep = '\t')
# get the filtered names
filtered_names <- get_filtered_pathway_names(pathways, pathway_mappings, 'R-HSA-168256')

# location of the lower resolution output
pathway_output_up_lowres_loc <- paste('/data/scRNA/pathways/mast/meta_paired_lores_lfc01minpct01_20201106/rna/sigs_pos/')
# read the 
pathway_up_df_low_monocyte <- get_pathway_table(pathway_output_up_lowres_loc, cell_types = c('monocyte'), append = '_sig_up_pathways.txt', use_ranking = T)
# location of the higher resolution output
pathway_output_up_highres_loc <- paste('/data/scRNA/pathways/mast/meta_paired_highres_lfc01minpct01_20210905/rna/sigs_pos/')
# read the 
pathway_up_df_high_monocyte <- get_pathway_table(pathway_output_up_highres_loc, cell_types = c('cMono', 'ncMono'), append = '_reactome.txt', use_ranking = T)
# add the rownames as a column, to allow the data tables to be merged
pathway_up_df_low_monocyte$pathway <- rownames(pathway_up_df_low_monocyte)
pathway_up_df_high_monocyte$pathway <- rownames(pathway_up_df_high_monocyte)
# merge together
pathway_up_df_monocyte <- merge(pathway_up_df_low_monocyte, pathway_up_df_high_monocyte, by='pathway', all = T)
# add back the rownames from the column we put them in
rownames(pathway_up_df_monocyte) <- pathway_up_df_monocyte$pathway
# we no longer need the rownames as a column, so let's remove it
pathway_up_df_monocyte$pathway <- NULL
# to zero transformation after merging
pathway_up_df_monocyte[is.na(pathway_up_df_monocyte)] <- 0
# remove the number prepend, it should not be necessary if using only one data source
rownames(pathway_up_df_monocyte) <- gsub('\\d+_', '', rownames(pathway_up_df_monocyte))
# plot all mono
heatmap.3(pathway_up_df_monocyte, to_na=0, dendrogram = 'none', labRow = NA)
# now for each pathogen
pdf('~/Desktop/all_pathways_monocyte.pdf', width=9, height=12)
for(timepoint in c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # subset to that timepoint
  pathway_up_df_monocyte_tp <- pathway_up_df_monocyte[, colnames(pathway_up_df_monocyte)[grep(timepoint, colnames(pathway_up_df_monocyte))]]
  # remove the labels of the cell type
  colnames(pathway_up_df_monocyte_tp) <- gsub(paste('UT', timepoint, sep = ''), '', colnames(pathway_up_df_monocyte_tp))
  # remove completely empty rows
  pathway_up_df_monocyte_tp <- pathway_up_df_monocyte_tp[rowSums(pathway_up_df_monocyte_tp) != 0, ]
  # plot
  heatmap.3(pathway_up_df_monocyte_tp, to_na=0, dendrogram = 'none', labRow = NA, main = paste('monocyte', 'UT', 'vs', timepoint), margins = c(12,4))
}
dev.off()
# do the same, but select on variation first
pdf('~/Desktop/varying_pathways_monocyte.pdf', width=9, height=12)
for(timepoint in c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # subset to that timepoint
  pathway_up_df_monocyte_tp <- pathway_up_df_monocyte[, colnames(pathway_up_df_monocyte)[grep(timepoint, colnames(pathway_up_df_monocyte))]]
  # remove the labels of the cell type
  colnames(pathway_up_df_monocyte_tp) <- gsub(paste('UT', timepoint, sep = ''), '', colnames(pathway_up_df_monocyte_tp))
  # remove completely empty rows
  pathway_up_df_monocyte_tp <- pathway_up_df_monocyte_tp[rowSums(pathway_up_df_monocyte_tp) != 0, ]
  # for variation selection, we can't leave the zeroes in
  pathway_up_df_monocyte_tp_zerosafe <- pathway_up_df_monocyte_tp
  pathway_up_df_monocyte_tp_zerosafe[pathway_up_df_monocyte_tp_zerosafe == 0] <- max(pathway_up_df_monocyte_tp_zerosafe)
  # select by variation
  most_varied_pathways <- get_most_varying_from_df(pathway_up_df_monocyte_tp_zerosafe, top_so_many = 20)
  pathway_up_df_monocyte_tp_vary <- pathway_up_df_monocyte_tp[most_varied_pathways, ]
  # plot
  heatmap.3(pathway_up_df_monocyte_tp_vary, to_na=0, dendrogram = 'none', main = paste('monocyte', 'UT', 'vs', timepoint), margins = c(12,20))
}
dev.off()

# now only look specifically at immune pathways
pathway_up_df_monocyte_immuneonly <- pathway_up_df_monocyte[rownames(pathway_up_df_monocyte) %in% filtered_names, ]
# now for each pathogen
pdf('~/Desktop/immune_pathways_monocyte.pdf', width=9, height=12)
for(timepoint in c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # subset to that timepoint
  pathway_up_df_monocyte_immuneonly_tp <- pathway_up_df_monocyte_immuneonly[, colnames(pathway_up_df_monocyte_immuneonly)[grep(timepoint, colnames(pathway_up_df_monocyte))]]
  # remove the labels of the cell type
  colnames(pathway_up_df_monocyte_immuneonly_tp) <- gsub(paste('UT', timepoint, sep = ''), '', colnames(pathway_up_df_monocyte_immuneonly_tp))
  # remove completely empty rows
  pathway_up_df_monocyte_immuneonly_tp <- pathway_up_df_monocyte_immuneonly_tp[rowSums(pathway_up_df_monocyte_immuneonly_tp) != 0, ]
  # plot
  heatmap.3(pathway_up_df_monocyte_immuneonly_tp, to_na=0, dendrogram = 'none', labRow = NA, main = paste('monocyte', 'UT', 'vs', timepoint), margins = c(12,4))
}
dev.off()
# do the same, but select on variation first
pdf('~/Desktop/immune_varying_pathways_monocyte.pdf', width=9, height=12)
for(timepoint in c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # subset to that timepoint
  pathway_up_df_monocyte_immuneonly_tp <- pathway_up_df_monocyte_immuneonly[, colnames(pathway_up_df_monocyte_immuneonly)[grep(timepoint, colnames(pathway_up_df_monocyte_immuneonly))]]
  # remove the labels of the cell type
  colnames(pathway_up_df_monocyte_immuneonly_tp) <- gsub(paste('UT', timepoint, sep = ''), '', colnames(pathway_up_df_monocyte_immuneonly_tp))
  # remove completely empty rows
  pathway_up_df_monocyte_immuneonly_tp <- pathway_up_df_monocyte_immuneonly_tp[rowSums(pathway_up_df_monocyte_immuneonly_tp) != 0, ]
  # for variation selection, we can't leave the zeroes in
  pathway_up_df_monocyte_immuneonly_tp_zerosafe <- pathway_up_df_monocyte_immuneonly_tp
  pathway_up_df_monocyte_immuneonly_tp_zerosafe[pathway_up_df_monocyte_immuneonly_tp_zerosafe == 0] <- max(pathway_up_df_monocyte_immuneonly_tp_zerosafe)
  # select by variation
  most_varied_pathways <- get_most_varying_from_df(pathway_up_df_monocyte_immuneonly_tp_zerosafe[ ,c('cMono', 'ncMono')], top_so_many = 10)
  pathway_up_df_monocyte_immuneonly_tp_vary <- pathway_up_df_monocyte_immuneonly_tp[most_varied_pathways, ]
  # plot
  heatmap.3(pathway_up_df_monocyte_immuneonly_tp_vary, to_na=0, dendrogram = 'none', main = paste('monocyte', 'UT', 'vs', timepoint), margins = c(12,20))
}
dev.off()




# read the 
pathway_up_df_high_monocyte_p <- get_pathway_table(pathway_output_up_highres_loc, cell_types = c('cMono', 'ncMono'), append = '_reactome.txt', use_ranking = F)
# to zero transformation after merging
pathway_up_df_monocyte_p <- pathway_up_df_high_monocyte_p
pathway_up_df_monocyte_p[is.na(pathway_up_df_monocyte_p)] <- 0
# remove the number prepend, it should not be necessary if using only one data source
rownames(pathway_up_df_monocyte_p) <- gsub('\\d+_', '', rownames(pathway_up_df_monocyte_p))
# get the fraction of genes in each pathway
pathway_up_df_monocyte_fractions_list <- get_number_of_hits_pathways('/data/scRNA/pathways/mast/meta_paired_highres_lfc01minpct01_20210905/rna/sigs_pos/', append = '_reactome.txt', cell_types = c('cMono', 'ncMono'))
# convert to a table
pathway_up_df_monocyte_fractions <- pathway_hits_to_table(pathway_up_df_monocyte_fractions_list)
pathway_up_df_monocyte_fractions[is.na(pathway_up_df_monocyte_fractions)] <- 0
# remove full zero rows
pathway_up_df_monocyte_p <- pathway_up_df_monocyte_p[rowSums(pathway_up_df_monocyte_p) != 0, ]
pathway_up_df_monocyte_fractions <- pathway_up_df_monocyte_fractions[rowSums(pathway_up_df_monocyte_fractions) != 0, ]
# confine to only immune
pathway_up_df_monocyte_immuneonly_p <- pathway_up_df_monocyte_p[rownames(pathway_up_df_monocyte_p) %in% filtered_names, ]
pathway_up_df_monocyte_fractions_immuneonly <- pathway_up_df_monocyte_fractions[rownames(pathway_up_df_monocyte_fractions) %in% filtered_names, ]
# now we'll plot everything
pdf('~/Desktop/immune_varying_pathways_monocyte_p.pdf', width=9, height=12)
for(timepoint in c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # subset to that timepoint
  pathway_up_df_monocyte_immuneonly_tp_p <- pathway_up_df_monocyte_immuneonly_p[, colnames(pathway_up_df_monocyte_immuneonly_p)[grep(timepoint, colnames(pathway_up_df_monocyte_immuneonly_p))]]
  pathway_up_df_monocyte_fractions_immuneonly_tp <- pathway_up_df_monocyte_fractions_immuneonly[, colnames(pathway_up_df_monocyte_fractions_immuneonly)[grep(timepoint, colnames(pathway_up_df_monocyte_fractions_immuneonly))]]
  # remove the labels of the cell type
  #colnames(pathway_up_df_monocyte_immuneonly_tp_p) <- gsub(paste('UT', timepoint, sep = ''), '', colnames(pathway_up_df_monocyte_immuneonly_tp_p))
  # remove completely empty rows
  pathway_up_df_monocyte_immuneonly_tp_p <- pathway_up_df_monocyte_immuneonly_tp_p[rowSums(pathway_up_df_monocyte_immuneonly_tp_p) != 0, ]
  pathway_up_df_monocyte_fractions_immuneonly_tp <- pathway_up_df_monocyte_fractions_immuneonly_tp[rowSums(pathway_up_df_monocyte_fractions_immuneonly_tp) != 0, ]
  # for variation selection, we can't leave the zeroes in
  pathway_up_df_monocyte_immuneonly_tp_zerosafe_p <- pathway_up_df_monocyte_immuneonly_tp_p*-1
  #pathway_up_df_monocyte_immuneonly_tp_zerosafe_p[pathway_up_df_monocyte_immuneonly_tp_zerosafe_p == 0] <- max(pathway_up_df_monocyte_immuneonly_tp_zerosafe_p)
  # select by variation
  most_varied_pathways_p <- get_most_varying_from_df(pathway_up_df_monocyte_immuneonly_tp_zerosafe_p, top_so_many = 10)
  pathway_up_df_monocyte_immuneonly_tp_vary_p <- pathway_up_df_monocyte_immuneonly_tp_p[most_varied_pathways_p, ]
  pathway_up_df_monocyte_fractions_immuneonly_tp_vary <- pathway_up_df_monocyte_fractions_immuneonly_tp[most_varied_pathways_p, ]
  # plot
  #heatmap.3(pathway_up_df_monocyte_immuneonly_tp_vary_p*-1, to_na=0, dendrogram = 'none', main = paste('monocyte', 'UT', 'vs', timepoint), margins = c(12,20), breaks = seq(from = 0, to = (max(pathway_up_df_monocyte_immuneonly_tp_vary_p*-1)), length.out = 10), col=(brewer.pal(9,"YlOrRd")))
  plots_condition_celltype <- create_dotplot_per_condition_and_celltypes(pathway_up_df_monocyte_immuneonly_tp_vary_p, pathway_up_df_monocyte_fractions_immuneonly_tp_vary, cell_type_sets = list('monocyte' = c('ncMono', 'cMono')), stims = c(timepoint))
  plots_condition_celltype[[timepoint]][['monocyte']]
}
dev.off()


pathway_up_df_monocyte_immuneonly_p <- pathway_up_df_monocyte_p[rownames(pathway_up_df_monocyte_p) %in% filtered_names, ]
pathway_up_df_monocyte_fractions_immuneonly <- pathway_up_df_monocyte_fractions[rownames(pathway_up_df_monocyte_fractions) %in% filtered_names, ]
plots_per_cond_p <- create_dotplot_per_condition_and_celltypes(pathway_up_df_monocyte_immuneonly_p, pathway_up_df_monocyte_fractions_immuneonly, cell_type_sets = list('monocyte' = c('ncMono', 'cMono')), use_most_varied = T)

pathway_up_df_high_monocyte_rank <- pathway_up_df_high_monocyte
rownames(pathway_up_df_high_monocyte_rank) <- gsub('\\d+_', '',rownames(pathway_up_df_high_monocyte_rank))
pathway_up_df_high_monocyte_immuneonly_rank <- pathway_up_df_high_monocyte_rank[rownames(pathway_up_df_high_monocyte_rank) %in% filtered_names, ]
plots_per_cond_rank <- create_dotplot_per_condition_and_celltypes(pathway_up_df_high_monocyte_immuneonly_rank, pathway_up_df_monocyte_fractions_immuneonly, cell_type_sets = list('monocyte' = c('ncMono', 'cMono')), use_most_varied = T)

# read the 
pathway_up_df_monocyte_p_nominal <- get_pathway_table(pathway_output_up_highres_loc, cell_types = c('cMono', 'ncMono'), append = '_reactome.txt', use_ranking = F, sig_val_to_store = 'p.value')
# remove the number prepend, it should not be necessary if using only one data source
rownames(pathway_up_df_monocyte_p_nominal) <- gsub('\\d+_', '', rownames(pathway_up_df_monocyte_p_nominal))
# filter by immune
pathway_up_df_monocyte_immuneonly_p_nominal <- pathway_up_df_monocyte_p_nominal[rownames(pathway_up_df_monocyte_p_nominal) %in% filtered_names, ]
# to zero transformation after merging
pathway_up_df_monocyte_immuneonly_p_nominal[is.na(pathway_up_df_monocyte_immuneonly_p_nominal)] <- 0
# plot
plots_per_cond_nominal <- create_dotplot_per_condition_and_celltypes(pathway_up_df_monocyte_immuneonly_p_nominal, pathway_up_df_monocyte_fractions_immuneonly, cell_type_sets = list('monocyte' = c('ncMono', 'cMono')), use_most_varied = T)


# read the 
pathway_up_df_all_p_fdr <- get_pathway_table(pathway_output_up_highres_loc, append = '_reactome.txt', use_ranking = F, cell_types = c('CD4 Naive', 'CD4 Memory', 'pDC', 'mDC','CD8 Naive', 'CD8 Memory', 'NKbright', 'NKdim','cMono', 'ncMono'), sig_val_to_use = 'q.value.FDR.B.Y')
rownames(pathway_up_df_all_p_fdr) <- gsub('\\d+_', '', rownames(pathway_up_df_all_p_fdr))
pathway_up_df_all_p_fdr <- pathway_up_df_all_p_fdr[rownames(pathway_up_df_all_p_fdr) %in% filtered_names, ]
# spaces to dots
colnames(pathway_up_df_all_p_fdr) <- gsub(' ', '.', colnames(pathway_up_df_all_p_fdr))
# get the fraction of genes in each pathway
pathway_up_df_all_fractions_list <- get_number_of_hits_pathways('/data/scRNA/pathways/mast/meta_paired_highres_lfc01minpct01_20210905/rna/sigs_pos/', append = '_reactome.txt', cell_types = c('CD4 Naive', 'CD4 Memory', 'pDC', 'mDC','CD8 Naive', 'CD8 Memory', 'NKbright', 'NKdim','cMono', 'ncMono'))
pathway_up_df_all_fractions_list <- pathway_hits_to_table(pathway_up_df_all_fractions_list)
pathway_up_df_all_fractions_list[is.na(pathway_up_df_all_fractions_list)] <- 0
# to zero transformation after merging
pathway_up_df_all_p_fdr[is.na(pathway_up_df_all_p_fdr)] <- 0
plots_per_cond_fdr <- create_dotplot_per_condition_and_celltypes(pathway_up_df_all_p_fdr, pathway_up_df_all_fractions_list, use_most_varied = T, cell_type_sets=list('CD4T' = c('CD4.Naive', 'CD4.Memory'), 'DC' = c('pDC', 'mDC'), 'CD8T' = c('CD8.Naive', 'CD8.Memory'), 'NK' = c('NKbright', 'NKdim'), 'monocyte' = c('cMono', 'ncMono')), minimal_difference = 1)

print_listed_ggplot_objects(plots_per_cond_fdr, '/data/scRNA/pathways/plots/de_pathways/meta_paired_highres_lfc01minpct01_20210905/rna/sigs_pos/', extention='pdf', width = 10, height = 10, xlab=xlab('cell type'), labs=labs(size='fraction', color='-log P value'))


write_dataframe_per_condition_and_celltypes(pathway_up_df_all_p_fdr, pathway_up_df_all_fractions_list, '/Users/royoelen/Desktop/figureS3/', use_most_varied = T, cell_type_sets=list('CD4T' = c('CD4.Naive', 'CD4.Memory'), 'DC' = c('pDC', 'mDC'), 'CD8T' = c('CD8.Naive', 'CD8.Memory'), 'NK' = c('NKbright', 'NKdim'), 'monocyte' = c('cMono', 'ncMono')), minimal_difference = 1)

