
######################
# libraries          #
######################

library(data.table)
require("heatmap.plus")
library(RColorBrewer)
library(VennDiagram)


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
  return(retval)
}


get_pathway_tables <- function(pathway_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
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
        filepath <- paste(pathway_output_loc, cell_type, 'UT',stim,'_sig_up_pathways.txt', sep = '')
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

get_genes_pathways_lists <- function(pathway_list, sig_col='q.value.Bonferroni', sig_cutoff=0.05){
  # create list to put everything in
  genes_lists <- list()
  # check each pathway
  for(key in names(pathway_list)){
    pathway_df <- pathway_list[[key]]
    print(key)
    # get the genes
    genes_pathways <- get_genes_pathways(pathway_df, sig_col, sig_cutoff)
    # put result in gene list
    genes_lists[[key]] <- genes_pathways
  }
  return(genes_lists)
}

get_genes_pathways <- function(pathway_table, sig_col='q.value.Bonferroni', sig_cutoff=0.05){
  # check each row
  genes <- apply(pathway_table, 1, function(x){
    # get the P we want
    sig_val <- as.numeric(as.vector(unlist(x[sig_col])))
    # get the gene list
    genes_line <- x['Hit.in.Query.List']
    # check if below the cutoff
    if(sig_val < sig_cutoff){
      # flatten value
      genes_line <- as.vector(unlist(genes_line))
      # remove the quotes
      genes_line <- gsub('\"','', genes_line)
      # split by comma
      genes_in_line <- strsplit(genes_line, ',')
      return(genes_in_line[[1]])
    }
    else{
      return(c())
    }
  })
  # turn into flat vector
  genes <- unlist(as.vector(genes))
  # get only the unique ones
  genes <- unique(genes)
  return(genes)
}

get_lfcs_for_genes <- function(mast_output_loc, genes_to_use, pval_column='metap_bonferroni', sig_pval=0.05, lfc_column='metafc', na_to_zero=F, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
  # init dataframe
  lfc_df <- NULL
  # check each cell type
  for(cell_type in cell_types){
    # check each stim
    for(stim in stims){
      try({
        print(paste(cell_type, stim, sep = ' '))
        # get the location of the MAST output file
        output_loc <- paste(mast_output_loc, cell_type, 'UT', stim, '.tsv', sep = '')
        # read the table
        mast_output <- read.table(output_loc, header = T, row.names = 1, sep = '\t')
        # filter on significant
        mast_output <- mast_output[mast_output[[pval_column]] < sig_pval, ]
        # grab the genes I care about
        mast_output_genes <- mast_output[genes_to_use, ]
        # R makes the rownames also NA for rows it can't find, so we need to redo those crownames
        rownames(mast_output_genes) <- genes_to_use
        # I only care about the LFC
        mast_output_genes_lfc <- mast_output_genes[, lfc_column, drop=F] # drop=F because R does not like dfs with one column
        # set the colname as the comparison
        colnames(mast_output_genes_lfc) <- c(paste(cell_type, 'UT', stim, sep = ''))
        # add to dataframe
        if(is.null(lfc_df)){
          lfc_df <- mast_output_genes_lfc
        }
        else{
          # we select the genes each time, the order is the same, so a cbind is safe
          lfc_df <- cbind(lfc_df, mast_output_genes_lfc)
        }
      })
    }
    # convert NA to zero if requested
    if(na_to_zero){
      lfc_df[is.na(lfc_df)] <- 0
    }
  }
  return(lfc_df)
}


get_top_vary_genes <- function(de_table, use_tp=T, use_pathogen=T, use_ct=T, sd_cutoff=0.5, use_dynamic_sd=F, top_so_many=10, must_be_positive_once=F, pathogens=c("CA", "MTB", "PA"), timepoints=c("3h", "24h"), cell_types=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")){
  top_vary_de <- c()
  cols_to_loop <- NULL
  # grab the appriate grep
  if(use_tp & use_pathogen){
    # we want a combo of pathogen and timepoints, so 3hCA for example
    cols_to_loop <- paste(rep(timepoints, each = length(pathogens)), pathogens, sep = "")
  }
  else if(use_pathogen & use_ct){
    # cell type and pathogen, so monocyte3hCA and monocyte24hCA for example
    cols_to_loop <- paste(rep(cell_types, each = length(pathogens)), pathogens, sep = ".*")
  }
  else if(use_tp & use_ct){
    # cell type at a timepoint, so monocyte3hCA and monocyte3hPA and monocyte3hMTB for example
    cols_to_loop <- paste(rep(cell_types, each = length(timepoints)), timepoints, sep = ".*")
  }
  else if(use_pathogen){
    cols_to_loop <- pathogens
  }
  else if(use_tp){
    cols_to_loop <- timepoints
  }
  else if(use_ct){
    cols_to_loop <- cell_types
  }
  # go through our group of columns
  for(col_grep in cols_to_loop){
    # grab the column names that have this in their name
    appropriate_columns <- colnames(de_table)[(grep(col_grep, colnames(de_table)))]
    print('getting most varying out of: ')
    print(appropriate_columns)
    # now subset the frame to only have these columns
    sub_de_table <- de_table[, appropriate_columns]
    # subset to only the genes that were upregulated at least once, if requested
    if(must_be_positive_once){
      sub_de_table <- sub_de_table[apply(sub_de_table,1,min) < 0,]
    }
    # we will return the rownames
    varying_genes <- NULL
    # either use a set SD or grab so many genes
    if(use_dynamic_sd){
      varying_genes <- get_most_varying_from_df(sub_de_table, top_so_many)
    }
    else{
      # now calculate the sd over this set of columns
      sds <- apply(sub_de_table, 1, sd, na.rm=T)
      # then grab the genes that are 'this' varied
      varying_genes <- rownames(sub_de_table[sds > sd_cutoff,])
    }
    
    # and add them to the list
    top_vary_de <- c(top_vary_de, varying_genes)
  }
  # constrain to the unique genes
  top_vary_de <- unique(top_vary_de)
  top_vary_de <- sort(top_vary_de)
  return(top_vary_de)
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


avg_exp_table_to_hm_table <- function(expression_table){
  # initialise table
  hm_table <- NULL
  # go through the conditions
  for(condition in unique(expression_table$condition)){
    # get the expression for that table
    expression_table_cond <- expression_table[expression_table$condition == condition, c('gene', 'average')]
    # set colnames so that we can merge these later
    colnames(expression_table_cond) <- c('gene', condition)
    # convert to data.table for efficient merging
    expression_table_cond <- data.table(expression_table_cond)
    # try to merge if necessary
    if(is.null(hm_table)){
      hm_table <- expression_table_cond
    }
    else{
      hm_table <- merge(hm_table, expression_table_cond, by='gene')
    }
  }
  # convert back to regular dataframe
  hm_table <- data.frame(hm_table)
  # set rownames
  rownames(hm_table) <- hm_table$gene
  # remove the old gene column
  hm_table$gene <- NULL
  return(hm_table)
}


get_gene_list_from_hm_branch <- function(heatmap, branch_directions, use_col=T){
  branch <- NULL
  # grab the row or column 
  if(use_col){
    branch <- heatmap$colDendrogram
  }
  else{
    branch <- heatmap$rowDendrogram
  }
  # go through the branch depths to get to the specific branch
  for(branch_direction in branch_directions){
    # the branch direction is 1 for left and 2 for right
    branch <- branch[[branch_direction]]
  }
  # now get all the children of this branch
  genes <- get_child_genes(branch)
  return(genes)
}

get_child_genes <- function(heatmap_branch){
  genes <- c()
  # if the height is zero, we are at the leaf and we can get the gene
  if(attr(heatmap_branch, 'height') == 0){
    genes <- names(attr(heatmap_branch, 'value'))
  }
  # if we are not at zero, there are more branches or leaves and we need to go further down both directions
  else{
    genes_branch1 <- get_child_genes(heatmap_branch[[1]])
    genes_branch2 <- get_child_genes(heatmap_branch[[2]])
    genes <- c(genes_branch1, genes_branch2)
  }
  return(genes)
}


pathways_to_hm_colors <- function(expression_heatmap, pathways_named_lists){
  colors_df <- NULL
  for(pathway_name in names(pathways_named_lists)){
    # get the pathway genes
    pathway.genes <- read.table(pathways_named_lists[[pathway_name]], header=F)
    pathway.genes <- as.character(pathway.genes$V1)
    pathway.genes <- pathway.genes[pathway.genes %in% rownames(expression_heatmap)]
    pathway.annotation <- rep("gray97", nrow(expression_heatmap))
    pathway.annotation[rownames(expression_heatmap) %in% pathway.genes] <- "gray55"
    # add to colors df
    if(is.null(colors_df)){
      colors_df <- data.frame(pathway.annotation)
      colnames(colors_df) <- pathway_name
    }
    else{
      colors_df[[pathway_name]] <- pathway.annotation
    }
  }
  # transform to matrix
  colors_m <- as.matrix(colors_df)
  return(colors_m)
}


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

# the location of the pathway output
pathways_up_loc <- '/data/scRNA/pathways/mast/meta_paired_lores_lfc01minpct01_20201106/rna/sigs_pos/'
# the location of the DE output
mast_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/'
# location of the expression
v2_exp_loc <- '/data/scRNA/expression/1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029_avgexp_rna.tsv'
v3_exp_loc <- '/data/scRNA/expression/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106_avgexp_rna.tsv'

# get all the pathway output
pathway_list <- get_pathway_tables(pathways_up_loc, cell_types = c('monocyte'))
# filter on immune
filtered_pathway_list <- filter_pathways_on_category(pathway_list, filtered_names)
# get the genes in the lists
pathway_genes_per_cond <- get_genes_pathways_lists(filtered_pathway_list)
# get just all the genes
de_pathway_genes <- unique(as.vector(unlist(pathway_genes_per_cond)))
# create an lfc dataframe
lfc_de_df <- get_lfcs_for_genes(mast_output_loc, de_pathway_genes, cell_types = c('monocyte'), na_to_zero = T)
# get the most varying
lfc_de_most_vary <- get_top_vary_genes(lfc_de_df, use_ct=T, use_pathogen = F, use_tp = F, use_dynamic_sd = T, top_so_many = 100, cell_types = c('monocyte'))
# 
lfc_de_varying_df <- lfc_de_df[lfc_de_most_vary, ]
# remove the monocyte and UTX monniker from the colnames
colnames(lfc_de_varying_df) <- gsub('monocyteUTX', '', colnames(lfc_de_varying_df))
# make heatmap based on LFC
heatmap.3(t(lfc_de_varying_df), col=rev(brewer.pal(10,"RdBu")), margins=c(6,8), to_na = 0)

# read the expression table
v2_expression <- read.table(v2_exp_loc, header = T, sep = '\t')
# subset to only the cell type and genes we care about
v2_expression_mono_de <- v2_expression[v2_expression$gene %in% lfc_de_most_vary & v2_expression$cell_type == 'monocyte', ]
# turn into the format for the heatmap
v2_expression_mono_de_hm <- avg_exp_table_to_hm_table(v2_expression_mono_de)
# scale to the max of the row
v2_expression_mono_de_hm <- apply(v2_expression_mono_de_hm, 1, function(x){x <- x/max(x)})

# grab some pathway genes to use for annotation
pathways_list <- list()
# 3h
pathways_list[['Interferon Signalling']] <- "/data/scRNA/pathways/REACTOME_Interferon_Signaling_genes.txt"
# 3h
pathways_list[['Interleukin-1 signalling']] <- "/data/scRNA/pathways/REACTOME_Interleukin-1_signaling.txt"
# 3h
#pathways_list[['MyD88 cascade initiated on plasma membrane']] <- "/data/scRNA/pathways/REACTOME_MyD88_cascade_initiated_on_plasma_membrane.txt"
# 3h
#pathways_list[['MyD88 dependent cascade initiated on endosome']] <- "/data/scRNA/pathways/REACTOME_MyD88_dependent_cascade_initiated_on_endosome.txt"
# 3h
pathways_list[['Toll-Like Receptors Cascades']] <- "/data/scRNA/pathways/REACTOME_Toll-Like_Receptors_Cascades.txt"
# 24h
#pathways_list[['Dectin-1 mediated noncanonical NF-kB signaling']] <- "/data/scRNA/pathways/REACTOME_Dectin-1_mediated_noncanonical_NF-kB_signaling.txt"
# 24h, bit in 3h
#pathways_list[['CLEC7A (Dectin-1) signaling']] <- "/data/scRNA/pathways/REACTOME_CLEC7A_(Dectin-1)_signaling.txt"
# 24h
#pathways_list[['antigen presenting']] <- "/data/scRNA/pathways/REACTOME_Antigen_processing-Cross_presentation.txt"
# 24h
#pathways_list[['Cross-presentation of soluble exogenous antigens']] <- "/data/scRNA/pathways/REACTOME_Cross-presentation_of_soluble_exogenous_antigens.txt"
# 24h
#pathways_list[['Regulation of RAS by GAPs']] <- "/data/scRNA/pathways/REACTOME_Regulation_of_RAS_by_GAPs.txt"
# 24h
pathways_list[['C-type lectin receptors']] <- "/data/scRNA/pathways/REACTOME_C-type_lectin_receptors.txt"
# all
#pathways_list[['Neutrophil degranulation']] <- "/data/scRNA/pathways/REACTOME_Neutrophil_degranulation.txt"
# all
#pathways_list[['Cytokine signalling']] <- "/data/scRNA/pathways/REACTOME_Cytokine_Signaling_in_Immune_system_genes.txt"


# add all the genes from the pathways togeter
pathway_genes <- c()
for(pathway in names(pathways_list)){
  genes_pathway_loc <- pathways_list[[pathway]]
  genes_pathway <- read.table(genes_pathway_loc, header=F, stringsAsFactors = F)$V1
  pathway_genes <- c(pathway_genes, genes_pathway)
}
pathway_genes <- unique(pathway_genes)

# transform the pathways to colors
colors_pathways_v2_de_hm <- pathways_to_hm_colors(t(v2_expression_mono_de_hm), pathways_list)

# plot that stuff
heatmap.3(v2_expression_mono_de_hm, col=rev(brewer.pal(10,"RdBu")), margins=c(6,8), to_na = 0, dendrogram = 'none', labCol = NA, ColSideColors = colors_pathways_v2_de_hm)




# now with all DE genes
v2_expression <- read.table(v2_exp_loc, header = T, sep = '\t')
v2_expression_mono_de_all <- v2_expression[v2_expression$gene %in% rownames(lfc_de_df) & v2_expression$cell_type == 'monocyte', ]
v2_expression_mono_de_all_hm <- avg_exp_table_to_hm_table(v2_expression_mono_de_all)
v2_expression_mono_de_all_hm <- data.frame(t(apply(v2_expression_mono_de_all_hm, 1, function(x){x <- x/max(x)})))
colnames(v2_expression_mono_de_all_hm) <- gsub('X', '', colnames(v2_expression_mono_de_all_hm))
colors_pathways_v2_de_all_hm <- pathways_to_hm_colors(v2_expression_mono_de_all_hm, pathways_list)
heatmap.3(t(v2_expression_mono_de_all_hm), col=rev(brewer.pal(10,"RdBu")), margins=c(6,8), to_na = 0, dendrogram = 'none', labCol = NA, ColSideColors = colors_pathways_v2_de_all_hm, ColSideColorsSize = 3, main = 'Differentially Expressed Genes', xlab = 'genes', ylab = 'conditions', cexRow = 1.5, side.height.fraction = 0.6, KeyValueName = 'expression')

# and V3
v3_expression <- read.table(v3_exp_loc, header = T, sep = '\t')
v3_expression_mono_de_all <- v3_expression[v3_expression$gene %in% rownames(lfc_de_df) & v3_expression$cell_type == 'monocyte', ]
v3_expression_mono_de_all_hm <- avg_exp_table_to_hm_table(v3_expression_mono_de_all)
v3_expression_mono_de_all_hm <- data.frame(t(apply(v3_expression_mono_de_all_hm, 1, function(x){x <- x/max(x)})))
colnames(v3_expression_mono_de_all_hm) <- gsub('X', '', colnames(v3_expression_mono_de_all_hm))
colors_pathways_v3_de_all_hm <- pathways_to_hm_colors(v3_expression_mono_de_all_hm, pathways_list)
heatmap.3(t(v3_expression_mono_de_all_hm), col=rev(brewer.pal(10,"RdBu")), margins=c(6,8), to_na = 0, dendrogram = 'none', labCol = NA, ColSideColors = colors_pathways_v3_de_all_hm, ColSideColorsSize = 3, main = 'Differentially Expressed Genes', xlab = 'genes', ylab = 'conditions', cexRow = 1.5, side.height.fraction = 0.6, KeyValueName = 'expression')




