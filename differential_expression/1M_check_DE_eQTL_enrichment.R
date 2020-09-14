


get_overlapping_DE_reQTL_genes_overlap <- function(condition1, condition2, cell_type, eQTL_output_loc, MAST_output_loc, pval_column='metap_bonferroni', sig_pval=0.05){
  # get the full path to the differential expression results
  MAST_output_loc_full <- paste(MAST_output_loc, cell_type, condition1, 'X', condition2, '.tsv', sep = '')
  # get the full path to the eQTL mapping results
  eQTL_output_loc_full <- paste(eQTL_output_loc, condition1, '_vs_', condition2, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
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
    sig_reqtl_genes <- as.character(reQTL$HGNCName)
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
  color_coding[["monocyte"]] <- "#E64B50"
  color_coding[["NK"]] <- "#EDBA1B"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  return(color_coding)
}



# location of the eQTL output
eQTL_output_loc <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/'
# location of the meta analysed MAST output
mast_meta_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/meta_paired_lores_lfc01minpct01_20200713/rna/'

# the conditions to look at
conditions <- c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')
#  the cell types to look at
cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')

# put in the work and call the function
get_overlapping_DE_reQTL_genes_overlap_per_condition(eQTL_output_loc, mast_meta_output_loc)
get_overlapping_DE_reQTL_genes_overlap_per_cell_type(eQTL_output_loc, mast_meta_output_loc)

