
library(ggplot2)
library(ggpubr)


get_overlapping_DE_reQTL_genes_overlap <- function(condition1, condition2, cell_type, eQTL_output_loc, MAST_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, only_positive=F, only_negative=F, reqtl=T, symbols.to.ensg.mapping.loc=NULL){
  # get the full path to the differential expression results
  MAST_output_loc_full <- paste(MAST_output_loc, cell_type, condition1, 'X', condition2, '.tsv', sep = '')
  # read eqtl or reqtl output
  eQTL_output_loc_full <- NULL
  if(reqtl){
    # get the full path to the eQTL mapping results
    eQTL_output_loc_full <- paste(eQTL_output_loc, condition1, '_vs_', condition2, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
  }
  else{
    # get the full path to the eQTL mapping results
    eQTL_output_loc_full <- paste(eQTL_output_loc, condition2, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
  }
  # init value
  overlap_percentage <- 0
  # try to read the two files
  try({
    # read the MAST output
    mast <- read.table(MAST_output_loc_full, sep = '\t', header = T, row.names = 1)
    # filter to only include the significant results
    mast <- mast[mast[[pval_column]] <= sig_pval, ]
    sig_de_genes <- rownames(mast)
    # read the eQTL output
    reQTL <- read.table(eQTL_output_loc_full, header = T, sep = '\t')
    # filter to only include the significant results
    reQTL <- reQTL[reQTL$FDR < sig_pval, ]
    # do mapping dependent on what is available
    sig_reqtl_genes <- NULL
    if(is.null(symbols.to.ensg.mapping.loc)){
      sig_reqtl_genes <- as.character(reQTL$HGNCName)
    }
    else{
      genes <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
      sig_reqtl_genes <- genes[match(as.character(reQTL$Probe), genes$ens),"gene"]
    }
    # check which are in both
    sig_both <- intersect(sig_de_genes, sig_reqtl_genes)
    # put this into numbers
    DE_genes_number <- length(sig_de_genes)
    reQTL_genes_number <- length(sig_reqtl_genes)
    sig_both_number <- length(sig_both)
    # get the percentage of reQTLs that are DE
    overlap_percentage <- sig_both_number/reQTL_genes_number
  })
  return(overlap_percentage)
}


get_overlapping_DE_reQTL_genes_overlap_per_condition <- function(eQTL_output_loc, MAST_output_loc, unstim_condition='UT', stim_conditions=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), pval_column='metap_bonferroni', sig_pval=0.05){
  # grab color codings
  colors <- get_color_coding_dict()
  # check the stim conditions
  for(stim in stim_conditions){
    # init values list
    percentages <- list()
    # check each cell type
    for(cell_type in cell_types){
      percentage <- get_overlapping_DE_reQTL_genes_overlap(unstim_condition, stim, cell_type, eQTL_output_loc, MAST_output_loc, pval_column, sig_pval)
      percentages[[cell_type]] <- percentage
    }
    barplot(unlist(percentages), main = 'reQTL DE genes fraction', sub = paste(unstim_condition, 'vs', stim), col = unlist(colors[names(percentages)]), ylim = c(0, 1))
  }
}

get_overlapping_DE_reQTL_genes_overlap_per_condition_ggplot <- function(eQTL_output_loc, MAST_output_loc, unstim_condition='UT', stim_conditions=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), pval_column='metap_bonferroni', sig_pval=0.05, reqtl=T, symbols.to.ensg.mapping.loc=NULL){
  plot_per_stim <- list()
  # grab color codings
  colors <- get_color_coding_dict()
  # check the stim conditions
  for(stim in stim_conditions){
    # init values list
    percentages <- list()
    # check each cell type
    for(cell_type in cell_types){
      percentage <- get_overlapping_DE_reQTL_genes_overlap(unstim_condition, stim, cell_type, eQTL_output_loc, MAST_output_loc, pval_column, sig_pval, reqtl, symbols.to.ensg.mapping.loc)
      percentages[[cell_type]] <- percentage
    }
    pc_df <- data.frame(cell_type=names(percentages), percentage=unlist(percentages))
    colScale <- scale_fill_manual(name = pc_df$cell_type, values = unlist(colors[cell_types]))
    p<-ggplot(pc_df, aes(x=cell_type, y=percentage, fill=cell_type)) +
      geom_bar(stat="identity")+
      colScale+
      theme_minimal() +
      ylim(0, 1) +
      theme(legend.position='none', axis.title.x=element_blank()) +
      ggtitle(stim)
    #
    plot_per_stim[[stim]] <- p
  }
  return(plot_per_stim)
}


get_overlapping_DE_reQTL_genes_overlap_per_cell_type <- function(eQTL_output_loc, MAST_output_loc, unstim_condition='UT', stim_conditions=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), pval_column='metap_bonferroni', sig_pval=0.05){
  # grab color codings
  colors <- get_color_coding_dict()
  # check each cell type
  for(cell_type in cell_types){
    # init values list
    percentages <- list()
    # check each condition
    for(stim in stim_conditions){
      percentage <- get_overlapping_DE_reQTL_genes_overlap(unstim_condition, stim, cell_type, eQTL_output_loc, MAST_output_loc, pval_column, sig_pval)
      percentages[[stim]] <- percentage
    }
    barplot(unlist(percentages), main = 'reQTL DE genes fraction', sub = paste(cell_type), col = unlist(colors[names(percentages)]), ylim = c(0, 1))
  }
}

get_overlapping_DE_reQTL_genes_overlap_per_cell_type_ggplot <- function(eQTL_output_loc, MAST_output_loc, unstim_condition='UT', stim_conditions=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), pval_column='metap_bonferroni', sig_pval=0.05, reqtl=T, symbols.to.ensg.mapping.loc=NULL){
  plot_per_stim <- list()
  # grab color codings
  colors <- get_color_coding_dict()
  # check each cell type
  for(cell_type in cell_types){
    # init values list
    percentages <- list()
    # check each condition
    for(stim in stim_conditions){
      percentage <- get_overlapping_DE_reQTL_genes_overlap(unstim_condition, stim, cell_type, eQTL_output_loc, MAST_output_loc, pval_column, sig_pval, reqtl, symbols.to.ensg.mapping.loc)
      percentages[[stim]] <- percentage
    }
    pc_df <- data.frame(condition=names(percentages), percentage=unlist(percentages))
    colScale <- scale_fill_manual(name = pc_df$condition, values = unlist(colors[conditions]))
    p<-ggplot(pc_df, aes(x=condition, y=percentage, fill=condition)) +
      geom_bar(stat="identity")+
      colScale+
      theme_minimal() +
      ylim(0, 1) +
      theme(legend.position='none', axis.title.x=element_blank()) +
      ggtitle(cell_type)
    #
    plot_per_stim[[cell_type]] <- p
  }
  return(plot_per_stim)
}

# 
get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["3hCA"]] <- "khaki2"
  color_coding[["24hCA"]] <- "khaki4"
  color_coding[["3hMTB"]] <- "paleturquoise1"
  color_coding[["24hMTB"]] <- "paleturquoise3"
  color_coding[["3hPA"]] <- "rosybrown1"
  color_coding[["24hPA"]] <- "rosybrown3"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  return(color_coding)
}



# location of the eQTL output
eQTL_output_loc <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/'
# location of the meta analysed MAST output
mast_meta_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/'

# the conditions to look at
conditions <- c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')
#  the cell types to look at
cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')

# gene to ensemble id mapping
symbols.to.ensg.mapping.loc <- '/data/scRNA/gene_to_ensemble.tsv'

# put in the work and call the function
de_reqtl_overlap_plot <- get_overlapping_DE_reQTL_genes_overlap_per_condition_ggplot(eQTL_output_loc, mast_meta_output_loc)
annotate_figure(ggarrange(de_reqtl_overlap_plot[['3hCA']], de_reqtl_overlap_plot[['3hMTB']],de_reqtl_overlap_plot[['3hPA']],
          de_reqtl_overlap_plot[['24hCA']], de_reqtl_overlap_plot[['24hMTB']], de_reqtl_overlap_plot[['24hPA']],ncol=3, nrow=2), top = 'reQTLs that have DE genes')
de_reqtl_overlap_plot_ct <- get_overlapping_DE_reQTL_genes_overlap_per_cell_type_ggplot(eQTL_output_loc, mast_meta_output_loc)
annotate_figure(ggarrange(de_reqtl_overlap_plot_ct[['bulk']], de_reqtl_overlap_plot_ct[['B']],de_reqtl_overlap_plot_ct[['CD4T']], de_reqtl_overlap_plot_ct[['CD8T']],
                          de_reqtl_overlap_plot_ct[['DC']], de_reqtl_overlap_plot_ct[['monocyte']],ncol=4, nrow=2), top = 'reQTLs that have DE genes')
# create the plots per eQTLs
de_eqtl_overlap_plot <- get_overlapping_DE_reQTL_genes_overlap_per_condition_ggplot(eQTL_output_loc, mast_meta_output_loc, reqtl=F, symbols.to.ensg.mapping.loc=NULL)
# create the arranged plot
annotate_figure(ggarrange(de_eqtl_overlap_plot[['3hCA']], de_eqtl_overlap_plot[['3hMTB']],de_eqtl_overlap_plot[['3hPA']],
                          de_eqtl_overlap_plot[['24hCA']], de_eqtl_overlap_plot[['24hMTB']], de_eqtl_overlap_plot[['24hPA']],ncol=3, nrow=2), top = 'eQTLs that have DE genes')

