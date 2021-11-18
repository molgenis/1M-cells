############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_higres_vs_lowres_eqtls.R
# Function: compare eQTL results between higher and lower resolution eQTL runs
############################################################################################################################


####################
# libraries        #
####################

library(UpSetR)
library(VennDiagram)

####################
# Functions        #
####################

compare_overlap_low_to_high_res_egenes <- function(upset_plot, low_res_output_loc, high_res_output_loc, cell_type_matches, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), sig_column='FDR', sig_cutoff=0.05, output_loc=NULL){
  result <- NULL
  if(upset_plot){
    result <- compare_overlap_low_to_high_res_egenes_upset(low_res_output_loc=low_res_output_loc, high_res_output_loc=high_res_output_loc, cell_type_matches=cell_type_matches, conditions=conditions, sig_column=sig_column, sig_cutoff=sig_cutoff)
  }
  else{
    result <- compare_overlap_low_to_high_res_egenes_venn(low_res_output_loc=low_res_output_loc, high_res_output_loc=high_res_output_loc, cell_type_matches=cell_type_matches, conditions=conditions, sig_column=sig_column, sig_cutoff=sig_cutoff, output_loc=output_loc)
  }
  return(result)
}



compare_overlap_low_to_high_res_egenes_upset <- function(low_res_output_loc, high_res_output_loc, cell_type_matches, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), sig_column='FDR', sig_cutoff=0.05){
  # we will have plot for each condition
  plots_per_condition <- list()
  # check each condition
  for(condition in conditions){
    # we will also have a plot per major cell type
    plot_per_major <- list()
    # each name in the cell type matches is the major cell type
    for(major in names(cell_type_matches)){
      # save the gene lists in a list we can use for upset
      gene_lists <- list()
      # paste together the path to the major MAST output
      major_output_loc <- paste(low_res_output_loc, condition, '/', major, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
      # read the output
      major_output <- read.table(major_output_loc, sep = '\t', header = T, stringsAsFactors = F)
      # get the genes that are significant
      sig_genes_major <- major_output[major_output[[sig_column]] < sig_cutoff, ][['ProbeName']]
      # add to the list
      gene_lists[[major]] <- sig_genes_major
      # get the high resolution cell types this major type is linked to
      minor_matches <- cell_type_matches[[major]]
      # check each of these minor cell types
      for(minor in minor_matches){
        try({
          # paste together the path to the major MAST output
          minor_output_loc <- paste(high_res_output_loc, condition, '/', minor, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
          # read the output
          minor_output <- read.table(minor_output_loc, sep = '\t', header = T, stringsAsFactors = F)
          # get the genes that are significant
          sig_genes_minor <- minor_output[minor_output[[sig_column]] < sig_cutoff, ][['ProbeName']]
          # add to list
          gene_lists[[minor]] <- sig_genes_minor
        })
      }
      # create the upset plot
      p <- upset(data = fromList(gene_lists), nsets = length(names(gene_lists)), order.by = 'freq')
      # add to the list
      plot_per_major[[major]] <- p
    }
    # add list of plots to condition
    plots_per_condition[[condition]] <- plot_per_major
  }
  return(plots_per_condition)
}


compare_overlap_low_to_high_res_egenes_venn <- function(low_res_output_loc, high_res_output_loc, cell_type_matches, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), sig_column='FDR', sig_cutoff=0.05, output_loc=NULL){
  # we will have the genes per cell types
  genes_per_major_cell_type <- list()
  # and we will have plot for each condition as well
  plot_per_major_cell_type <- list()
  # check each condition
  for(major in names(cell_type_matches)){
    # the genes of the major cell type
    genes_per_major_cell_type[[major]] <- list()
    # save the major one separately
    genes_per_major_cell_type[[major]][[major]]
    # check each condition
    for(condition in conditions){
      # paste together the path to the major MAST output
      major_output_loc <- paste(low_res_output_loc, condition, '/', major, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
      # read the output
      major_output <- read.table(major_output_loc, sep = '\t', header = T, stringsAsFactors = F)
      # get the genes that are significant
      sig_genes_major <- major_output[major_output[[sig_column]] < sig_cutoff, ][['ProbeName']]
      # add to the list
      genes_per_major_cell_type[[major]][[major]] <- c(genes_per_major_cell_type[[major]][[major]], sig_genes_major)
    }
    # make unique (same genes in mapped in some conditions)
    genes_per_major_cell_type[[major]][[major]] <- unique(genes_per_major_cell_type[[major]][[major]])
    # check for the minors
    minor_matches <- cell_type_matches[[major]]
    for(minor in minor_matches){
      # we'll save the results 
      plots_per_minor_cell_type <- genes_per_major_cell_type[[major]][[minor]] <- c()
      # get for each condition
      for(condition in conditions){
        # check each condition
        minor_output_loc <- paste(high_res_output_loc, condition, '/', minor, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
        try({
          # read the output
          minor_output <- read.table(minor_output_loc, sep = '\t', header = T, stringsAsFactors = F)
          # get the genes that are significant
          sig_genes_minor <- minor_output[minor_output[[sig_column]] < sig_cutoff, ][['ProbeName']]
          # stuff them in the list with the rest
          genes_per_major_cell_type[[major]][[minor]] <- c(genes_per_major_cell_type[[major]][[minor]], sig_genes_minor)
        })
      }
      genes_per_major_cell_type[[major]][[minor]] <- unique(genes_per_major_cell_type[[major]][[minor]])
    }
    # create the venn data
    plot_per_major_cell_type[[major]] <- list(x = genes_per_major_cell_type[[major]], filename=NULL, category.names=names(genes_per_major_cell_type[[major]]), 
                                              cat.col = c("#440154ff", '#21908dff', '#fde725ff'), col=c("#440154ff", '#21908dff', '#fde725ff'),
                                              fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)))
    # plot straight away if requested
    if(!is.null(output_loc)){
      result_loc <- paste(output_loc, '/', major, '.tiff', sep = '')
      print(paste('writing', result_loc))
      venn.diagram(x = genes_per_major_cell_type[[major]], filename=result_loc, category.names=names(genes_per_major_cell_type[[major]]), 
                   cat.col = c("#440154ff", '#21908dff', '#fde725ff'), col=c("#440154ff", '#21908dff', '#fde725ff'),
                   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)))
    }
  }
  return(plot_per_major_cell_type)
}




plot_venns <- function(venn_data_per_condition){
  
}


####################
# Main Code        #
####################

# set locations of output
lowres_eQTL_output <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/'
highres_eQTL_output <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_highres_20210905_confine_lead_snp_gene/results/'
highres_T_eQTL_output <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_hires_20211008_reclassified_T_eqtlgenlead/results/'

plots_per_condition <- compare_overlap_low_to_high_res_egenes(T, lowres_eQTL_output, highres_eQTL_output, list('DC' = c('mDC', 'pDC'), 'monocyte' = c('ncMono', 'cMono'), 'NK' = c('NKdim', 'NKbright')))

pdf('eQTL_overlap_monocyte.pdf', width = 5, height = 5)
plots_per_condition[['UT']][['monocyte']]
plots_per_condition[['3hCA']][['monocyte']]
plots_per_condition[['24hCA']][['monocyte']]
plots_per_condition[['3hMTB']][['monocyte']]
plots_per_condition[['24hMTB']][['monocyte']]
plots_per_condition[['3hPA']][['monocyte']]
plots_per_condition[['24hPA']][['monocyte']]
dev.off()

pdf('eQTL_overlap_NK.pdf', width = 5, height = 5)
plots_per_condition[['UT']][['NK']]
plots_per_condition[['3hCA']][['NK']]
plots_per_condition[['24hCA']][['NK']]
plots_per_condition[['3hMTB']][['NK']]
plots_per_condition[['24hMTB']][['NK']]
plots_per_condition[['3hPA']][['NK']]
plots_per_condition[['24hPA']][['NK']]
dev.off()

pdf('eQTL_overlap_DC.pdf', width = 5, height = 5)
plots_per_condition[['UT']][['DC']]
plots_per_condition[['3hCA']][['DC']]
plots_per_condition[['24hCA']][['DC']]
plots_per_condition[['3hMTB']][['DC']]
plots_per_condition[['24hMTB']][['DC']]
plots_per_condition[['3hPA']][['DC']]
plots_per_condition[['24hPA']][['DC']]
dev.off()

venn_data_per_condition <- compare_overlap_low_to_high_res_egenes(F, lowres_eQTL_output, highres_eQTL_output, list('DC' = c('mDC', 'pDC'), 'monocyte' = c('ncMono', 'cMono'), 'NK' = c('NKdim', 'NKbright')), output_loc='/data/scRNA/eQTL_mapping/major_vs_minor_overlap/')
venn_data_per_condition <- compare_overlap_low_to_high_res_egenes(F, lowres_eQTL_output, highres_T_eQTL_output, list('CD4T' = c('CD4_Naive', 'CD4_Memory'), 'CD8T' = c('CD8_Naive', 'CD8_Memory')), output_loc='/data/scRNA/eQTL_mapping/major_vs_minor_overlap/')


