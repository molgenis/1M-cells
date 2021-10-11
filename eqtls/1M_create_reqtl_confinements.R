############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_create_reqtl_confinements.R
# Function: read eQTL results and use those make confinement files for reQTL mapping
############################################################################################################################


####################
# libraries        #
####################

# none required, only base R

####################
# Functions        #
####################

write_confinements <- function(eqtl_output_loc, output_prepend, per_ct=F, file_name='eQTLsFDR0.05-ProbeLevel.txt.gz', output_append='_confinement.txt', conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_append='_expression', significance_column='FDR', significance_cutoff=0.05, snp_column='SNPName', gene_column='HGNCName', other_columns=c('OverallZScore'), top_esnp_only=F, sep='\t', header=T, row.names=NULL, col_names_write=NULL, sep_write='\t'){
  # get the full output
  eqtl_lists <- get_eqtls_celltypes_and_conditions(eqtl_output_loc, file_name=file_name, conditions=conditions, cell_types=cell_types, cell_append=cell_append, significance_column=significance_column, significance_cutoff=significance_cutoff, snp_column=snp_column, gene_column=gene_column, other_columns=other_columns, top_esnp_only=top_esnp_only, sep=sep, header=header, row.names=row.names)
  # write per cell type if required
  if(per_ct){
    # check each cell type then
    for(cell_type in names(eqtl_lists)){
      # grab that output
      eqtls_ct <- get_snp_gene_combinations(eqtl_lists[[cell_type]], gene_column = gene_column)
      # get a file name without spaces, because that is just easier to work with
      cell_type_safe <- gsub(' ', '_', cell_type)
      # paste together the output location
      confinement_file_loc <- paste(output_prepend, cell_type_safe, output_append, sep = '')
      # write the result
      write_confinement(eqtls_ct, confinement_file_loc, col_names_write = col_names_write, sep_write = sep_write)
    }
  }
  # if not per cell type, basically do the same, but with all the output at once
  else{
    # grab everything at once, regardless of cell type
    eqtls <- get_snp_gene_combinations(eqtl_lists, gene_column = gene_column)
    # paste together the output location
    confinement_file_loc <- paste(output_prepend, output_append, sep = '')
    # write the result
    write_confinement(eqtls, confinement_file_loc, col_names_write = col_names_write, sep_write = sep_write)
  }
}


write_confinement <- function(confinement, confinement_file_loc, col_names_write=NULL, sep_write='\t'){
  # write with colnames if requested
  if(!is.null(col_names_write)){
    # set those column names
    colnames(confinement) <- col_names_write
    write.table(confinement, confinement_file_loc, sep = sep_write, col.names = T, row.names = F, quote = F)
  }
  else{
    write.table(confinement, confinement_file_loc, sep = sep_write, col.names = F, row.names = F, quote = F)
  }
}


get_snp_gene_combinations <- function(eqtl_output, snp_column='SNPName', gene_column='HGNCName'){
  # initialize variable
  snp_gene_combinations <- NULL
  # it is a table, we can just grab the data
  if(is.data.frame(eqtl_output) | is.matrix(eqtl_output)){
    # if there is no list, it must be a table (coercing to dataframe)
    snp_gene_combinations <- data.frame(eqtl_output[, c(snp_column, gene_column)])
  }
  # if the output is a list, we need to recursively go deeper to get to the tables
  else if(is.list(eqtl_output)){
    # check each list element
    for(element_name in names(eqtl_output)){
      # recursively call this function again
      this_eqtl_output <- get_snp_gene_combinations(eqtl_output[[element_name]], snp_column = snp_column, gene_column = gene_column)
      # add to existing output
      if(is.null(snp_gene_combinations)){
        snp_gene_combinations <- this_eqtl_output
      }
      else{
        snp_gene_combinations <- rbind(snp_gene_combinations, this_eqtl_output)
      }
    }
  }
  # anything else is a mystery to me
  else{
    warning(paste('skipping because not list/data.frame/matrix, but', typeof(eqtl_output)))
  }
  # make unique, confinements only need the snp-probe combination once
  snp_gene_combinations <- unique(snp_gene_combinations)
  return(snp_gene_combinations)
}


get_eqtls_celltypes_and_conditions <- function(eqtl_output_loc, file_name='eQTLsFDR0.05-ProbeLevel.txt.gz', conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_append='_expression', significance_column='FDR', significance_cutoff=0.05, snp_column='SNPName', gene_column='HGNCName', other_columns=c('OverallZScore'), top_esnp_only=F, sep='\t', header=T, row.names=NULL){
  # store output per cell type
  sig_eqtls_per_cell_type <- list()
  # check each cell type
  for(cell_type in cell_types){
    # get list of tables for this cell type
    sig_eqtl_table_cell_type <- get_eqtls_conditions(eqtl_output_loc, file_name=file_name, conditions=conditions, cell_type, cell_append=cell_append, significance_column=significance_column, significance_cutoff=significance_cutoff, snp_column=snp_column, gene_column=gene_column, other_columns=other_columns, top_esnp_only=top_esnp_only, sep=sep, header=header, row.names=row.names)
    # add to the list with the cell types as keys
    sig_eqtls_per_cell_type[[cell_type]] <- sig_eqtl_table_cell_type
  }
  return(sig_eqtls_per_cell_type)
}


get_eqtls_conditions <- function(eqtl_output_loc, file_name='eQTLsFDR0.05-ProbeLevel.txt.gz', conditions=c(), cell_type, cell_append='_expression', significance_column='FDR', significance_cutoff=0.05, snp_column='SNPName', gene_column='HGNCName', other_columns=c('OverallZScore'), top_esnp_only=F, sep='\t', header=T, row.names=NULL){
  # store output per condition
  sig_eqtl_table_per_condition <- list()
  # check each condition
  for(condition in conditions){
    # paste together the location of the file
    eqtl_output_loc_full <- paste(eqtl_output_loc, '/', condition, '/', cell_type, cell_append, '/', file_name, sep = '')
    # get the significant eQTLs
    sig_eqtls_condition <- get_significant_snp_gene_combos(eqtl_output_loc_full, significance_column=significance_column, significance_cutoff=significance_cutoff, snp_column=snp_column, gene_column=gene_column, other_columns=other_columns, top_esnp_only=top_esnp_only, sep=sep, header=header, row.names=row.names)
    # add these to the list
    sig_eqtl_table_per_condition[[condition]] <- sig_eqtls_condition
  }
  return(sig_eqtl_table_per_condition)
}


get_significant_snp_gene_combos <- function(eqtl_output_loc, significance_column='FDR', significance_cutoff=0.05, snp_column='SNPName', gene_column='HGNCName', other_columns=c('OverallZScore'), top_esnp_only=F, sep='\t', header=T, row.names=NULL){
  eqtl_table <- tryCatch({
    # read the table
    eqtl_table <- read.table(eqtl_output_loc, sep = sep, header = header, row.names = row.names)
    # subset to what is significant and what we care about
    eqtl_table <- eqtl_table[eqtl_table[[significance_column]] < significance_cutoff, c(snp_column, gene_column, other_columns, significance_column)]
    # filter for only the top effect if requested
    if(top_esnp_only){
      # get the unique snp-gene pairs
      eqtl_table$snp_gene <- paste(eqtl_table[[snp_column]], eqtl_table[[gene_column]], sep = '_')
      snp_probe_combos <- unique(eqtl_table$snp_gene)
      # now sort by p value
      eqtl_table <- eqtl_table[order(eqtl_table[[significance_column]]), ]
      # now do a match. match stops at the first hit, which due to sorting, should also be the most significat
      eqtl_table <- eqtl_table[match(snp_probe_combos, eqtl_table$snp_gene), ]
      # remove the extra column we added
      eqtl_table$snp_gene <- NULL
    }
    return(eqtl_table)
  }, warning = function(w){
    print(paste('warning for file ', eqtl_output_loc, sep = ''))
    print(w)
    return(NULL)
  }, error = function(e){
    print(paste('skipped file ', eqtl_output_loc, sep = ''))
    print(e)
    return(NULL)
  })
  return(eqtl_table)
}


####################
# Main code        #
####################

# set the location of the of the eqtl output
eqtl_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_hires_20211008_reclassified_T_eqtlgenlead/results/'
# set the output of the confinement file
full_confinement_prepend <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/sig_hires_20211008_reclassified_T_eqtlgenlead'
# cell type specific prepend file
ct_specific <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/sig_hires_20211008_reclassified_T_eqtlgenlead_'
# these are the labels that are possible
labels_t_azimuth <- c('CD8 Naive', 'CD4 CTL', 'CD4 Naive', 'CD4 TCM', 'CD8 TEM', 'MAIT', 'CD8 TEM', 'dnT', 'CD4 TEM', 'CD8 TCM', 'CD8 Proliferating')
# we corrected for spaces
labels_t_azimuth_nospace <- gsub(' ', '_', labels_t_azimuth)
# write the confinement files
write_confinements(eqtl_output_loc, full_confinement_prepend, cell_types = labels_t_azimuth_nospace, gene_column = 'ProbeName')
write_confinements(eqtl_output_loc, ct_specific, per_ct = T, cell_types = labels_t_azimuth_nospace, gene_column = 'ProbeName')

# the same for unconfined now
eqtl_output_unconfined_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_hires_20211008_reclassified_T_unconfined/results/'
# set the output of the confinement file
full_unconfined_confinement_prepend <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/sig_hires_20211008_reclassified_T_unconfined'
# cell type specific prepend file
ct_specific_unconfined <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/sig_hires_20211008_reclassified_T_unconfined_'
# write the confinement files
write_confinements(eqtl_output_unconfined_loc, full_unconfined_confinement_prepend, cell_types = labels_t_azimuth_nospace, top_esnp_only = T, gene_column = 'ProbeName')
write_confinements(eqtl_output_unconfined_loc, ct_specific_unconfined, per_ct = T, cell_types = labels_t_azimuth_nospace, top_esnp_only = T, gene_column = 'ProbeName')



