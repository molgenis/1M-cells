
write_significant_reqtl_genes <- function(reqtl_output_loc, gene_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('3hCA', '24hCA', '3hPA', '24hPA', '3hMTB', '24hMTB'), sig_column='FDR', sig_threshold=0.05, gene_column='ProbeName'){
  # get the output
  reqtls_per_ct_and_cond <- get_reqtl_per_cell_type_and_conditions(reqtl_output_loc, cell_types = cell_types, stims = stims, sig_column = sig_column, sig_threshold = sig_threshold)
  # loop the cell types
  for(cell_type in names(reqtls_per_ct_and_cond)){
    # loop through the conditions
    for(condition in names(reqtls_per_ct_and_cond[[cell_type]])){
      # grab the table
      output <- reqtls_per_ct_and_cond[[cell_type]][[condition]]
      # see if there is anything significant left
      if(nrow(output) > 0){
        # grab the genes
        sig_genes <- output[[gene_column]]
        # turn into dataframe
        sig_genes_df <- data.frame(gene=sig_genes)
        # paste together location
        full_out_loc <- paste(gene_output_loc, cell_type, '_', condition, '.txt', sep = '')
        # write the table
        write.table(sig_genes, full_out_loc, row.names = F, col.names = F, quote = F)
      }
    }
  }
}

get_reqtl_per_cell_type_and_conditions <- function(reqtl_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('3hCA', '24hCA', '3hPA', '24hPA', '3hMTB', '24hMTB'), only_significant=T, sig_column='FDR', sig_threshold=0.05){
  # create list per cell type
  reqtl_per_ct <- list()
  # read each cell type
  for(cell_type in cell_types){
    # create a list per condition
    reqtl_per_condition <- list()
    # read each condition
    for(condition in stims){
      # paste together the location
      reqtl_loc <- paste(reqtl_output_loc, 'UT_vs_', condition, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      # read the table
      reqtl_output <- read.table(reqtl_loc, sep = '\t', header = T)
      # filter for significant if requested
      if(only_significant){
        reqtl_output <- reqtl_output[reqtl_output[[sig_column]] < sig_threshold, ]
      }
      # put it in the list
      reqtl_per_condition[[condition]] <- reqtl_output
    }
    # put it in the list
    reqtl_per_ct[[cell_type]] <- reqtl_per_condition
  }
  return(reqtl_per_ct)
}
