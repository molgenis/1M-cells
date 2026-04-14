#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_merge_qtl_outputs.R
# Function: merge all QTL tables
# 
############################################################################################################################

###################
# libraries       #
###################

library(data.table)


###################
# Functions       #
###################

merge_qtl_tables <- function(qtl_output_loc, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=NULL, qtl_file='eQTLsFDR-ProbeLevel.txt.gz', dirsnip=NULL, subdirsnip=NULL) {
  # we'll store each table
  tables_each <- list()
  # list the directories present
  dirs_present <- list.dirs(qtl_output_loc, recursive=FALSE, full.names = F)
  # overlap with conditions if given
  if (!is.null(conditions)) {
    dirs_present <- dirs_present[dirs_present %in% conditions]
  }
  # check each of the directories
  for (directory in dirs_present) {
    # paste the directory together
    qtl_dir_loc <- paste(qtl_output_loc, directory, sep = '/')
    # list the directories here
    subdirs_present <- list.dirs(qtl_dir_loc, recursive=FALSE, full.names = F)
    # overlap with cell types if given
    if (!is.null(cell_types)) {
      subdirs_present <- subdirs_present[subdirs_present %in% cell_types]
    }
    # check each of these subdirectories
    for (subdirector in subdirs_present) {
      # paste filepath together to QTL file
      qtl_file_loc <- paste(qtl_dir_loc, subdirector, qtl_file, sep = '/')
      # check if that one exists
      if (file.exists(qtl_file_loc)) {
        # read the table
        qtl_table <- fread(qtl_file_loc)
        # get condition
        condition <- directory
        if (!is.null(dirsnip)) {
          condition <- sub(dirsnip, '', condition)
        }
        # and cell type
        cell_type <- subdirector
        if (!is.null(subdirsnip)) {
          cell_type <- sub(subdirsnip, '', cell_type)
        }
        # add columns for condition and cell type
        qtl_table <- cbind(data.frame('condition' = rep(condition, times = nrow(qtl_table)), 'celltype' = rep(cell_type, times = nrow(qtl_table))), qtl_table)
        # store in list
        tables_each[[paste(condition, cell_type)]] <- qtl_table
      } else {
        print(paste('File not found:', qtl_file_loc))
      }
    }
  }
  # merge all
  merged_table <- rbindlist(tables_each)
  return(merged_table)
}


###################
# Settings        #
###################


###################
# Main code       #
###################

# location of the eQTL output of the confined mapping per condition
eqtl_confined_loc <- '/groups/umcg-franke-scrna/tmp04/releases/wijst-2020-hg19/v1/QTL_mapping/eQTL/output/sct_mqc_demux_lores_20201106_eqtlgenlead_metabeta/'
# get table
eqtl_confined <- merge_qtl_tables(eqtl_confined_loc, subdirsnip = '_expression')
# location of reqtls
reqtl_confined_loc <- '/groups/umcg-franke-scrna/tmp04/releases/wijst-2020-hg19/v1/QTL_mapping/eQTL/output/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_metabeta/'
# get table
reqtl_confined <- merge_qtl_tables(reqtl_confined_loc, subdirsnip = '_expression', conditions = c('UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA'))
# save the result
eqtl_confined_merged_loc <- '/groups/umcg-franke-scrna/tmp04/releases/wijst-2020-hg19/v1/QTL_mapping/eQTL/output/sct_mqc_demux_lores_20201106_eqtlgenlead_metabeta_merged.txt.xz'
# save this
write.table(rbind(eqtl_confined, reqtl_confined), xzfile(eqtl_confined_merged_loc, compression = 9), quote = F, row.names = F, col.names = T, sep = '\t')
# get unconfined mapping location
eqtl_unconfined_loc <- '/groups/umcg-franke-scrna/tmp04/releases/wijst-2020-hg19/v1/QTL_mapping/eQTL/output/sct_mqc_demux_lores_20201106_no_confine_metabeta/'
# get table
eqtl_unconfined <- merge_qtl_tables(eqtl_unconfined_loc, subdirsnip = '_expression')
# save the result
eqtl_unconfined_merged_loc <- '/groups/umcg-franke-scrna/tmp04/releases/wijst-2020-hg19/v1/QTL_mapping/eQTL/output/sct_mqc_demux_lores_20201106_no_confine_metabeta_merged.txt.xz'
# save this
write.table(eqtl_unconfined, xzfile(eqtl_unconfined_merged_loc, compression = 9), quote = F, row.names = F, col.names = T, sep = '\t')
