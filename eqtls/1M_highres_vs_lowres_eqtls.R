############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_higres_vs_lowres_eqtls.R
# Function: compare eQTL results between higher and lower resolution eQTL runs
############################################################################################################################


####################
# libraries        #
####################

library(UpSetR)

####################
# Functions        #
####################

compare_overlap_low_to_high_res_egenes <- function(low_res_output_loc, high_res_output_loc, cell_type_matches, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), sig_column='FDR', sig_cutoff=0.05){
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


####################
# Main Code        #
####################

# set locations of output
lowres_eQTL_output <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/'
highres_eQTL_output <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_highres_20210905_confine_lead_snp_gene/results/'

plots_per_condition <- compare_overlap_low_to_high_res_egenes(lowres_eQTL_output, highres_eQTL_output, list('DC' = c('mDC', 'pDC'), 'monocyte' = c('ncMono', 'cMono'), 'NK' = c('NKdim', 'NKbright')))

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


