####################
# libraries        #
####################


####################
# Functions        #
####################


harmonize_PRS_to_gt_format <- function(prs_table){
  # the genotype data has four digits, so those we at least need
  corrected_id_column <- apply(prs_table, 1, function(x){
    # get correct column
    id <- x['IID']
    # subset to the last four digits
    last_digits <- substr(id, nchar(id)-3, nchar(id))
    # get what is in front
    first_string <- substr(id, 1, nchar(id)-3)
    # reduce extra zeroes we don't need
    first_string_ll <- gsub('LL0+', '1_LLDeep_', first_string)
    # put it together again
    final_id <- paste(first_string_ll, last_digits, sep = '')
    # give that back
    return(final_id)
  })
  # convert to vector, as that is easier
  id_column_vector <- as.vector(unlist(corrected_id_column))
  # replace the IID column
  prs_table[['IID']] <- id_column_vector
  return(prs_table)
}


perform_eqts <- function(prs_table, expression_table, method='spearman', score_col='SCORESUM', keep_fails=T){
  # confine to data we have in both
  participants_both <- intersect(prs_table[['IID']], colnames(expression_table))
  prs_table <- prs_table[match(participants_both, prs_table[['IID']]), ]
  expression_table <- expression_table[, match(participants_both, colnames(expression_table))]
  # remove zero variation expression as that is useless to test
  expression_table <- expression_table[apply(expression_table, 1, var) != 0, ]
  # check the expression of each gene
  res <- apply(expression_table, 1, function(row){
    result_row <- NULL
    # get the number of participants
    n <- length(row)
    # depending on the data, sometimes doing the analysis, we need to catch these instances
    tryCatch({
      # calculate the correlation
      if(method == 'spearman'){
        result_gene <- cor.test(as.vector(unlist(row)), prs_table[[score_col]], method = method)
        result_row <- list( 'method' = method,
                          'alternative' = result_gene$alternative, 
                          'n' = n,
                          'S' = as.vector(unlist(result_gene$statistic['S'])), 
                          'rho' = as.vector(unlist(result_gene$estimate['rho'])), 
                          'p.value' = result_gene$p.value)
      }
      else{
        result_gene <- cor.test(as.vector(unlist(row)), prs_table[[score_col]], method = method)
        result_row <- list( 'method' = method,
                          'alternative' = result_gene$alternative, 
                          'n' = n,
                          'S' = as.vector(unlist(result_gene$statistic['S'])), 
                          'rho' = as.vector(unlist(result_gene$estimate['rho'])), 
                          'p.value' = result_gene$p.value)
      }
    },
    # so in case of an error, we will return an mostly empty row
    error = function(e){
      message('analysis impossible for a gene')
      print(e)
      # create empty entry
      result_row <- list( 'method' = method,
                          'alternative' = NA, 
                          'n' = n,
                          'S' = NA, 
                          'rho' = NA, 
                          'p.value' = NA)
    })
    return(result_row)
  })
  # turn into a dataframe
  res_df <- do.call(rbind.data.frame, res)
  # order by the p value
  res_df <- res_df[order(res_df[['p.value']]), ]
  # maybe we want to remove what we could not test, this is important to note when talking about the MTC!
  if(!keep_fails){
    res_df <- res_df[!is.na(res_df[['p.value']]), ]
  }
  # add MTC
  res_df$BH <- p.adjust(res_df$p.value, method = 'BH')
  res_df$bonferroni <- p.adjust(res_df$p.value, method = 'bonferroni')
  return(res_df)
}


perform_eqts_cts <- function(expression_loc, prs_table, result_loc, cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), method='spearman', score_col='SCORESUM', exp_file_append='_expression.tsv', report_sig_col=NULL, report_sig_cutoff=NULL){
  # check each cell type
  for(cell_type in cell_types){
    # paste together the location
    output_loc_cell_type <- paste(expression_loc, cell_type, exp_file_append, sep = '')
    # read the file
    expression_cell_type <- read.table(output_loc_cell_type, header = T, row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
    # perform the analysis
    result_cell_type <- perform_eqts(prs_table = prs_table, expression_table = expression_cell_type, method = method, score_col = score_col)
    # setup an output location
    output_loc_result <- paste(result_loc, '/', cell_type, '.tsv', sep = '')
    # write the result
    write.table(result_cell_type, output_loc_result, col.names = T, sep = '\t', quote = F)
    # report the findings if requested
    if(!is.null(report_sig_col)){
      number_of_sigs <- nrow(result_cell_type[result_cell_type[[report_sig_col]] < report_sig_cutoff, ])
      print(paste('found', as.character(number_of_sigs), 'in', output_loc_result))
    }
  }
}


perform_eqts_conditions <- function(expression_loc, prs_table, result_loc, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), method='spearman', score_col='SCORESUM', exp_file_append='_expression.tsv', report_sig_col=NULL, report_sig_cutoff=NULL){
  # check each timepoint
  for(condition in conditions){
    # paste the path together
    output_loc_condition <- paste(result_loc, '/', condition, '/', sep = '')
    # set the expression path
    expression_loc_condition <- paste(expression_loc, '/', condition, '/', sep = '')
    # do this analysis for this condition
    perform_eqts_cts(expression_loc_condition, prs_table, output_loc_condition, cell_types = cell_types, method=method, score_col=score_col, exp_file_append=exp_file_append, report_sig_col=report_sig_col, report_sig_cutoff=report_sig_cutoff)
  }
}


plot_eqts <- function(expression_table, prs_table, gene_name, score_col='SCORESUM'){
  # subset to what is available in both
  participants_both <- intersect(prs_table[['IID']], colnames(expression_table))
  # get the entries, in order
  prs_table <- prs_table[match(participants_both, prs_table[['IID']]), ]
  expression_table <- expression_table[, match(participants_both, colnames(expression_table))]
  # extract the gene
  expression_row  <- expression_table[gene_name, ]
  expression  <- as.vector(unlist(expression_row))
  # turn into a plotting table
  plot_table <- data.frame(part = participants_both, prs = prs_table[[score_col]], expression = expression, stringsAsFactors = F)
  # create the plot
  p <- ggplot(data = plot_table, mapping = aes(x = prs, y = expression)) +
    geom_point() +
    #geom_smooth(method = 'lm', formula = expression~prs) +
    xlab(paste('PRS')) + 
    ylab(paste('expression of', gene_name)) +
    ggtitle(paste('eQTS of PRS and', gene_name))
  
  return(p)
}


met_analyse_eqts <- function(named_eqts_outputs, output_loc){
  
}


get_exp_now <- function(condition, cell_type, path=NULL, chem=NULL){
  full_output_loc <- NULL
  # check if the user gave us a path
  if(!is.null(path)){
    full_output_loc <- paste(path, '/', condition, '/', cell_type, '_expression.tsv', sep = '')
  }
  else if(is.null(path) & !is.null(chem)){
    if(chem == 'v2'){
      full_output_loc <- paste(v2_exp_l, '/', condition, '/', cell_type, '_expression.tsv', sep = '')
    }
    else if(chem == 'v3'){
      full_output_loc <- paste(v3_exp_l, '/', condition, '/', cell_type, '_expression.tsv', sep = '')
    }
    else{
      stop('expression chemistry must be v2 or v3')
    }
  }
  else{
    stop('supply either the base expression path, or the chemistry (v2 or v3)')
  }
  # read the expression matrix
  expression <- read.table(full_output_loc, sep = '\t', header = T, row.names=1, check.names = F)
  return(expression)
}


get_res_now <- function(condition, cell_type, path=NULL, chem=NULL){
  full_output_loc <- NULL
  # check if the user gave us a path
  if(!is.null(path)){
    full_output_loc <- paste(path, '/', condition, '/', cell_type, '.tsv', sep = '')
  }
  else if(is.null(path) & !is.null(chem)){
    if(chem == 'v2'){
      full_output_loc <- paste(v2_res, '/', condition, '/', cell_type, '.tsv', sep = '')
    }
    else if(chem == 'v3'){
      full_output_loc <- paste(v3_res, '/', condition, '/', cell_type, '.tsv', sep = '')
    }
    else{
      stop('expression chemistry must be v2 or v3')
    }
  }
  else{
    stop('supply either the base expression path, or the chemistry (v2 or v3)')
  }
  # read the expression matrix
  result <- read.table(full_output_loc, sep = '\t', header = T, row.names=1, check.names = F)
  return(result)
}

# utility function to save me from passing large variable names around every time
set_global_feature_loc <- function(v2_exp_l='/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v2_sct_mqc_demux_lores_20201029/', v3_exp_l='/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v3_sct_mqc_demux_lores_20201106/'){
  # set the global paths of a couple of variables we keep using
  v2_exp_l <<- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v2_sct_mqc_demux_lores_20201029/'
  v3_exp_l <<- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v3_sct_mqc_demux_lores_20201106/'
  v2_res <<-'/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/PRS/v2_sct_mqc_demux_lores_20201029/'
  v3_res <<-'/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/PRS/v3_sct_mqc_demux_lores_20201106/'
}

#########################
# main code             #
#########################

# location of the files
prs_ra_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/PRS/RA_GWASmeta_European_v2/ra_Okada2013_full_lld_pgs_20211104.txt'
prs_sle_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/PRS/ea_imputed_mergedChroms_removedNa_removedDuplicates_selectedColumns/sle_Langefeld2017_full_lld_pgs_20211104.txt'
expression_data_root_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v2_sct_mqc_demux_lores_20201029/'
expression_data_root_loc2 <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v3_sct_mqc_demux_lores_20201106/'

# location of the eQTS output
eqts_sle_output_loc_v2 <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/PRS/SLE/v2_sct_mqc_demux_lores_20201029/'
eqts_sle_output_loc_v3 <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/PRS/SLE/v3_sct_mqc_demux_lores_20201106/'
eqts_ra_output_loc_v2 <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/PRS/RA/v2_sct_mqc_demux_lores_20201029/'
eqts_ra_output_loc_v3 <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/PRS/RA/v3_sct_mqc_demux_lores_20201106/'


# set up which cell types and conditions we want to use
cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
conditions <- c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')

# read the prs table
prs_ra_table <- read.table(prs_ra_loc, sep = '\t', header = T, stringsAsFactors = F)
# turn into compatible format
#prs_ra_table <- harmonize_PRS_to_gt_format(prs_ra_table)
# do the eqts stuff
perform_eqts_conditions(expression_data_root_loc, prs_ra_table, eqts_ra_output_loc_v2, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), method='spearman', score_col='ra_pgs_phiauto', exp_file_append='_expression.tsv', report_sig_col=NULL, report_sig_cutoff=NULL)
perform_eqts_conditions(expression_data_root_loc2, prs_ra_table, eqts_ra_output_loc_v3, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), method='spearman', score_col='ra_pgs_phiauto', exp_file_append='_expression.tsv', report_sig_col=NULL, report_sig_cutoff=NULL)


# read the prs table
prs_sle_table <- read.table(prs_sle_loc, sep = '\t', header = T, stringsAsFactors = F)
# turn into compatible format
#prs_sle_table <- harmonize_PRS_to_gt_format(prs_sle_table)
# do the eqts stuff
perform_eqts_conditions(expression_data_root_loc, prs_sle_table, eqts_sle_output_loc_v2, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), method='spearman', score_col='SCORESUM', exp_file_append='_expression.tsv', report_sig_col=NULL, report_sig_cutoff=NULL)
perform_eqts_conditions(expression_data_root_loc2, prs_sle_table, eqts_sle_output_loc_v3, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), method='spearman', score_col='SCORESUM', exp_file_append='_expression.tsv', report_sig_col=NULL, report_sig_cutoff=NULL)

