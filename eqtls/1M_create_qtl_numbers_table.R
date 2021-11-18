############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_create_qtl_numbers_table.R
# Function: get tables describing the number of QTLs in each condition and cell type
############################################################################################################################


####################
# libraries        #
####################



####################
# Functions        #
####################

get_number_of_egenes <- function(qtl_output_loc, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_append='_expression', output_file='eQTLsFDR0.05-ProbeLevel.txt.gz', sig_column='FDR', sig_cutoff=0.05, probe_column='ProbeName'){
  # initialize the dataframe to store the results in
  numbers <- data.frame(condition = rep(conditions, each = length(cell_types)), cell_type = rep(cell_types, times = length(conditions)), egenes = rep(NA, times = length(conditions) * length(cell_types)), stringsAsFactors = F)
  # go through the conditions
  for(condition in conditions){
    # paste the more specific path together
    qtl_output_condition_loc <- paste(qtl_output_loc, condition, '/', sep = '')
    # check each cell type
    for(cell_type in cell_types){
      # paste together the output for this specific cell type
      qtl_output_condition_celltype_loc <- paste(qtl_output_condition_loc, '/', cell_type, cell_type_append, '/', output_file, sep = '')
      # it's possible that we were unable to run something, so we might not be able to find an output file
      tryCatch({
        # read the output file
        qtl_output_file <- read.table(qtl_output_condition_celltype_loc, sep = '\t', header = T, stringsAsFactors = F)
        # set initial variable for number of hits
        nr_unique_probes <- 0
        # check if there is anything significant
        if(nrow(qtl_output_file <- qtl_output_file[qtl_output_file[[sig_column]] < sig_cutoff, ]) > 0){
          # filter by what is considered significant
          qtl_output_file <- qtl_output_file[qtl_output_file[[sig_column]] < sig_cutoff, ]
          # count the number of unique genes
          nr_unique_probes <- length(unique(qtl_output_file[[probe_column]]))
        }
        # add this entry
        numbers[numbers$condition == condition & numbers$cell_type == cell_type, 'egenes'] <- nr_unique_probes
      }, 
      error = function(cond){
        # let the user know there was no output file
        message(paste('error on file', qtl_output_condition_celltype_loc, '\n'))
        message(cond)
      })
    }
  }
  return(numbers)
}


# transform tall table
#1        UT         B     50
#2        UT      CD4T    610
#3        UT      CD8T    420
#4        UT        DC    144
#5        UT  monocyte    560
#6        UT        NK    289
#
# to wide table
#
tall_to_wide_table <- function(table, column_as_rownames, column_as_colnames, data_variable){
  # get the unique variables for each columns to use as row and column names
  row_names <- unique(table[[column_as_rownames]])
  col_names <- unique(table[[column_as_colnames]])
  # initialize as matrix
  wide_table <- matrix(, ncol=length(col_names), nrow=length(row_names), dimnames = list(row_names, col_names))
  # now check each entry in the first column
  for(row_entry in row_names){
    # check against each entry in the other column
    for(col_entry in col_names){
      # see if there is a variable
      if(nrow(table[table[[column_as_rownames]] == row_entry & table[[column_as_colnames]] == col_entry, ]) > 0){
        # get what is hopefully one entry
        variable <- table[table[[column_as_rownames]] == row_entry & table[[column_as_colnames]] == col_entry, data_variable]
        # warn if there is more than one matched variable
        if(length(variable) > 1){
          base::warning(message = 'there is more than one match for this combination of column and row names, using the first one')
        }
        # use the first entry
        variable <- variable[1]
        # set that as value
        wide_table[row_entry, col_entry] <- variable
      }
      else{
        base::warning(message = paste('there is no entry for', row_entry ,'and', col_entry, '\nstays NA variable'))
      }
    }
  }
  # convert to data.table because I prefer those
  wide_table <- data.frame(wide_table, stringsAsFactors = F)
  return(wide_table)
}


####################
# Main Code        #
####################

# location of the output
eqtl_meta_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/'
eqtl_meta_lowres_loc <- paste(eqtl_meta_loc, 'sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/', sep = '')
eqtl_meta_highres_loc <- paste(eqtl_meta_loc, 'sct_mqc_demux_highres_20210905_20211008_reclassified_T_combined_confine_lead_snp_gene/results/', sep = '')


# get the results for the lower resolution table
eqtl_meta_lowres_table <- get_number_of_egenes(eqtl_meta_lowres_loc, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA', 'UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA'))
# convert to other format if that one is preferable
eqtl_meta_lowres_table_wide <- tall_to_wide_table(eqtl_meta_lowres_table, 'condition', 'cell_type', 'egenes')
# higher resolution as well
eqtl_meta_highres_table <- get_number_of_egenes(eqtl_meta_highres_loc, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4_Memory', 'CD4_Naive', 'CD8_Memory', 'CD8_Naive', 'NKbright',  'NKdim', 'cMono', 'ncMono', 'mDC', 'pDC', 'plasma_B'))
# and other format
eqtl_meta_highres_table_wide <- tall_to_wide_table(eqtl_meta_highres_table, 'condition', 'cell_type', 'egenes')




