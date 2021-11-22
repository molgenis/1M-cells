



compare_overlap_low_to_high_de <- function(low_res_output_loc, high_res_output_loc, cell_type_matches, conditions=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), sig_column='metap_bonferroni', sig_cutoff=0.05){
  # we will have the genes per cell types
  all_genes_major <- c()
  # and we will have plot for each condition as well
  all_genes_minor <- c()
  # check each condition
  for(major in names(cell_type_matches)){
    for(condition in conditions){
      try({
        # paste path together
        output_loc_low <- paste(low_res_output_loc, major, 'UTX', condition, '.tsv', sep = '')
        # read the output
        output_low <- read.table(output_loc_low, header = T, row.names = 1, sep = '\t', stringsAsFactors = F)
        # subset to what is significant
        output_low <- output_low[output_low[[sig_column]] < sig_cutoff, ]
        # add to the genes vector
        all_genes_major <- c(all_genes_major, rownames(output_low))
      })
    }
    # next check tne minors
    for(minor in cell_type_matches[[major]]){
      for(condition in conditions){
        try({
          # paste path together
          output_loc_high <- paste(high_res_output_loc, minor, 'UTX', condition, '.tsv', sep = '')
          # read the output
          output_high <- read.table(output_loc_high, header = T, row.names = 1, sep = '\t', stringsAsFactors = F)
          # subset to what is significant
          output_high <- output_high[output_high[[sig_column]] < sig_cutoff, ]
          # add to the genes vector
          all_genes_minor <- c(all_genes_minor, rownames(output_high))
        })
      }
    }
  }
  
  # make the minor and major ones unique
  all_genes_major <- unique(all_genes_major)
  all_genes_minor <- unique(all_genes_minor)
  
  # setdiff of what you want to get
  only_major <- setdiff(all_genes_major, all_genes_minor)
  only_minor <- setdiff(all_genes_minor, all_genes_major)
  both <- intersect(all_genes_minor, all_genes_major)
  print(paste(
    'total major: ', (length(all_genes_major)),
    ', total minor: ', (length(all_genes_minor)),
    ', only major: ', (length(only_major)), 
    ', only_minor: ', (length(only_minor)),
    ', both: ', (length(both)),
    sep = ''
  ))
}

de_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/'
lowres_de_output_loc <- paste(de_output_loc, '/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/', sep = '')
highres_de_output_loc <- paste(de_output_loc, '/paired_highres_and_hires_t_lfc01minpct01_20210905_merged/meta_paired_highres_and_hires_t_lfc01minpct01_20210905_merged/rna/', sep = '')

compare_overlap_low_to_high_de(lowres_de_output_loc, highres_de_output_loc, list('CD4T' = c('CD4 Naive', 'CD4 Memory'), 'CD8T' = c('CD8 Naive', 'CD8 Memory'), 'DC' = c('mDC', 'pDC'), 'monocyte' = c('ncMono', 'cMono'), 'NK' = c('NKdim', 'NKbright'), 'B' = c()),
                               )
