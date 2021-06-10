require("heatmap.plus")
library(RColorBrewer)
library(VennDiagram)

get_mast_meta_output_overlap <- function(mast_meta_output_1, mast_meta_output_2, venn_output_loc='./', only_significant=T, pval_column = 'p_val_adj', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='avg_logFC', only_ribosomal=F, symbols.to.ensg.mapping='./', group1name='group 1', group2name='group 2'){
  # list the files in directory 1
  files <- list.files(mast_meta_output_1)
  # check the files
  for(file in files){
    # set the full path
    mast1_loc <- paste(mast_meta_output_1, '/', file, sep='')
    mast2_loc <- paste(mast_meta_output_2, '/', file, sep='')
    try({
      # try to read both tables
      mast1 <- read.table(mast1_loc, header=T, row.names = 1, sep = '\t')
      mast1[[lfc_column]] <- mast1[[lfc_column]] * -1
      mast2 <- read.table(mast2_loc, header=T, row.names = 1, sep = '\t')
      # need to add the hgnc names
      if('hgnc.names' %in% colnames(mast1)){
        mast1$hgnc.names <- as.character(mast1$hgnc.names)
      }
      else{
        mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
        mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
        mast1$hgnc.names <- rownames(mast1)
        rownames(mast1) <- mapping[match(rownames(mast1), mapping$V2),"V1"]
      }
      if('hgnc.names' %in% colnames(mast2)){
        mast2$hgnc.names <- as.character(mast2$hgnc.names)
      }
      else{
        mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
        mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
        mast2$hgnc.names <- rownames(mast2)
        rownames(mast2) <- mapping[match(rownames(mast2), mapping$V2),"V1"]
      }
      # confine to ribosomal if requested
      if(only_ribosomal){
        mast1 <- mast1[startsWith(mast1$hgnc.names, 'RPS') | startsWith(mast1$hgnc.names, 'RPL'), ]
        mast2 <- mast2[startsWith(mast2$hgnc.names, 'RPS') | startsWith(mast2$hgnc.names, 'RPL'), ]
      }
      # remove insignificant results if requested
      if(only_significant){
        mast1 <- mast1[mast1[[pval_column]] <= 0.05, ]
        mast2 <- mast2[mast2[[pval_column]] <= 0.05, ]
      }
      # only positive if requested
      if(only_positive){
        mast1 <- mast1[mast1[[lfc_column]] < 0, ]
        mast2 <- mast2[mast2[[lfc_column]] < 0, ]
      }
      # only negateive if requested
      if(only_negative){
        mast1 <- mast1[mast1[[lfc_column]] > 0, ]
        mast2 <- mast2[mast2[[lfc_column]] > 0, ]
      }
      # report something about the P-values
      quantile_mast1 <- quantile(mast1[[pval_column]])
      quantile_mast2 <- quantile(mast2[[pval_column]])
      print(paste(substr(file, 1, regexpr("\\.[^\\.]*$", file)-1), group1name, 'quantiles:'))
      print(quantile_mast1)
      print(paste(substr(file, 1, regexpr("\\.[^\\.]*$", file)-1), group2name, 'quantiles:'))
      print(quantile_mast2)
      # grab the genes
      mast1_genes <- rownames(mast1)
      mast2_genes <- rownames(mast2)
      # grab the name of the file without the extention
      myCol <- brewer.pal(3, "Pastel2")[1:2]
      output_file <- substr(file, 1, regexpr("\\.[^\\.]*$", file)-1)
      venn.diagram(x = list(mast1_genes, mast2_genes),
                   main = substr(file, 1, regexpr("\\.[^\\.]*$", file)-1),
                   category.names = c(group1name, group2name),
                   filename = paste(venn_output_loc, output_file, '.png', sep = ''),
                   imagetype="png" ,
                   height = 600 , 
                   width = 600 , 
                   resolution = 300,
                   compression = "lzw",
                   lwd = 2,
                   lty = 'blank',
                   fill = myCol,
                   cex = .6,
                   fontface = "bold",
                   fontfamily = "sans",
                   cat.cex = 0.6,
                   cat.fontface = "bold",
                   cat.default.pos = "outer",
                   cat.pos = c(-27, 27),
                   cat.dist = c(0.055, 0.055),
                   cat.fontfamily = "sans")
    })
  }
}

plot_concordances <- function(mast_meta_output_1, mast_meta_output_2, plot_output_loc='./', only_significant=T, pval_column = 'p_val_adj', sig_pval=0.05, lfc_column='avg_logFC', only_ribosomal=F, symbols.to.ensg.mapping='./', group1name='group 1', group2name='group 2'){
  # list the files in directory 1
  files <- list.files(mast_meta_output_1)
  # check the files
  for(file in files){
    # set the full path
    mast1_loc <- paste(mast_meta_output_1, '/', file, sep='')
    mast2_loc <- paste(mast_meta_output_2, '/', file, sep='')
    try({
      # try to read both tables
      mast1 <- read.table(mast1_loc, header=T, row.names = 1, sep = '\t')
      mast1[[lfc_column]] <- mast1[[lfc_column]] * -1
      mast2 <- read.table(mast2_loc, header=T, row.names = 1, sep = '\t')
      # need to add the hgnc names
      if('hgnc.names' %in% colnames(mast1)){
        mast1$hgnc.names <- as.character(mast1$hgnc.names)
      }
      else{
        mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
        mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
        mast1$hgnc.names <- rownames(mast1)
        rownames(mast1) <- mapping[match(rownames(mast1), mapping$V2),"V1"]
      }
      if('hgnc.names' %in% colnames(mast2)){
        mast2$hgnc.names <- as.character(mast2$hgnc.names)
      }
      else{
        mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
        mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
        mast2$hgnc.names <- rownames(mast2)
        rownames(mast2) <- mapping[match(rownames(mast2), mapping$V2),"V1"]
      }
      # confine to ribosomal if requested
      if(only_ribosomal){
        mast1 <- mast1[startsWith(mast1$hgnc.names, 'RPS') | startsWith(mast1$hgnc.names, 'RPL'), ]
        mast2 <- mast2[startsWith(mast2$hgnc.names, 'RPS') | startsWith(mast2$hgnc.names, 'RPL'), ]
      }
      # remove insignificant results if requested
      if(only_significant){
        mast1 <- mast1[mast1[[pval_column]] <= 0.05, ]
        mast2 <- mast2[mast2[[pval_column]] <= 0.05, ]
      }
      mast_merged <- merge(mast1, mast2, by=0)
      rownames(mast_merged) <- mast_merged$row.names
      mast_merged$row.names <- NULL
      mast_merged <- mast_merged[order(mast_merged[[paste(lfc_column, '.x', sep='')]]), ]
      plot(mast_merged[[paste(lfc_column, '.x', sep='')]], mast_merged[[paste(lfc_column, '.y', sep='')]], xlab = paste(group1name, 'logFC'), ylab=paste(group2name, 'logFC'), main=file, xlim=c(-2, 2), ylim=c(-2, 2))
      abline(lm(mast_merged[[paste(lfc_column, '.x', sep='')]]~mast_merged[[paste(lfc_column, '.y', sep='')]]))
    })
  }
}


# MAST output locations
cells_1M_MAST_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/v2_paired_lores_lfc01minpct01_20201106/rna/'
pilot4_MAST_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/pilot4_DE_output/'
symbols.to.ensg.mapping <- '/data/scRNA/differential_expression/genesymbol_to_ensid.tsv'
venn_overlap_loc <- '/data/scRNA/differential_expression/seurat_MAST/overlap/paired_lores_rna_lfc01minpct01_20201106_v2_vs_pilot4/rna/'
venn_overlap_negative_loc <- paste(venn_overlap_loc, 'negative/', sep='')
venn_overlap_positive_loc <- paste(venn_overlap_loc, 'positive/', sep='')
venn_overlap_negative_ribosomal_loc <- paste(venn_overlap_loc, 'negative_ribosomal/', sep='')
venn_overlap_positive_ribosomal_loc <- paste(venn_overlap_loc, 'positive_ribosomal/', sep='')
lfc_concordance_loc <- '/data/scRNA/differential_expression/seurat_MAST/concordance/paired_lores_rna_lfc01minpct01_20201106_v2_vs_pilot4/rna/'

# create overlap plots
get_mast_meta_output_overlap(pilot4_MAST_output_loc, cells_1M_MAST_output_loc, venn_output_loc=venn_overlap_positive_loc, only_significant=T, pval_column = 'p_val_adj', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=T, only_negative=F, lfc_column='avg_logFC', only_ribosomal=F, symbols.to.ensg.mapping=symbols.to.ensg.mapping, group1name='p4', group2name='1m')
get_mast_meta_output_overlap(pilot4_MAST_output_loc, cells_1M_MAST_output_loc, venn_output_loc=venn_overlap_negative_loc, only_significant=T, pval_column = 'p_val_adj', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=T, lfc_column='avg_logFC', only_ribosomal=F, symbols.to.ensg.mapping=symbols.to.ensg.mapping, group1name='p4', group2name='1m')
get_mast_meta_output_overlap(pilot4_MAST_output_loc, cells_1M_MAST_output_loc, venn_output_loc=venn_overlap_positive_ribosomal_loc, only_significant=T, pval_column = 'p_val_adj', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=T, only_negative=F, lfc_column='avg_logFC', only_ribosomal=T, symbols.to.ensg.mapping=symbols.to.ensg.mapping, group1name='p4', group2name='1m')
get_mast_meta_output_overlap(pilot4_MAST_output_loc, cells_1M_MAST_output_loc, venn_output_loc=venn_overlap_negative_ribosomal_loc, only_significant=T, pval_column = 'p_val_adj', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=T, lfc_column='avg_logFC', only_ribosomal=T, symbols.to.ensg.mapping=symbols.to.ensg.mapping, group1name='p4', group2name='1m')

plot_concordances(pilot4_MAST_output_loc, cells_1M_MAST_output_loc, plot_output_loc=lfc_concordance_loc, only_significant=T, pval_column = 'p_val_adj', sig_pval=0.05, lfc_column='avg_logFC', only_ribosomal=T, symbols.to.ensg.mapping=symbols.to.ensg.mapping, group1name='p4', group2name='1m')
  
