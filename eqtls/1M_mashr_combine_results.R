######################
# libraries          #
######################

# mashr libs
library(stringr)
library(mashr)
# image creation libs
library(VennDiagram)
library(UpSetR)

####################
# Functions        #
####################

get_conditions_mash <- function(eqtl_output_loc, conditions = c(), cell_type){
  betas = NULL
  ses = NULL
  # go through the conditions
  for(condition in conditions){
    # there might be some tables not present, so let's set to a default
    table <- NULL
    # read the table
    table_loc <- paste(eqtl_output_loc, condition, '/', cell_type, '_expression/', 'eQTLsFDR-ProbeLevel.txt.gz',  sep = '')
    # show progress
    print(paste('reading', table_loc))
    try({
      table <- read.table(table_loc, sep = '\t', header = T)
    })
    if(is.null(table)){
      print(paste('table not present for:', table_loc))
    } else{
      # grab the variables
      table_b1 <- as.numeric(str_match(table$Beta..SE., "(.+) \\x28.+\\x29;")[,2])
      table_se1 <- as.numeric(str_match(table$Beta..SE., ".+ \\x28(.+)\\x29;")[,2])
      table_b2 <- as.numeric(str_match(table$Beta..SE., ";(.+) \\x28.+\\x29")[,2])
      table_se2 <- as.numeric(str_match(table$Beta..SE., ";.+ \\x28(.+)\\x29")[,2])
      # add the SNP>probe as rownames
      table_b1 <- data.frame(table_b1, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
      table_se1 <- data.frame(table_se1, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
      table_b2 <- data.frame(table_b2, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
      table_se2 <- data.frame(table_se2, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
      # add the colnames
      colnames(table_b1) <- c(paste(condition, 'v2', sep = '.'))
      colnames(table_se1) <- c(paste(condition, 'v2', sep = '.'))
      colnames(table_b2) <- c(paste(condition, 'v3', sep = '.'))
      colnames(table_se2) <- c(paste(condition, 'v3', sep = '.'))
      # banana
      if(is.null(ses)){
        ses <- merge(table_se1, table_se2, by="row.names",all.x=TRUE)
        rownames(ses) <- ses$Row.names
        ses$Row.names <- NULL
      } else{
        ses <- merge(ses, table_se1, by="row.names",all.x=TRUE)
        rownames(ses) <- ses$Row.names
        ses$Row.names <- NULL
        ses <- merge(ses, table_se2, by="row.names",all.x=TRUE)
        rownames(ses) <- ses$Row.names
        ses$Row.names <- NULL
      }
      if(is.null(betas)){
        betas <- merge(table_b1, table_b2, by="row.names",all.x=TRUE)
        rownames(betas) <- betas$Row.names
        betas$Row.names <- NULL
      } else{
        betas <- merge(betas, table_b1, by="row.names",all.x=TRUE)
        rownames(betas) <- betas$Row.names
        betas$Row.names <- NULL
        betas <- merge(betas, table_b2, by="row.names",all.x=TRUE)
        rownames(betas) <- betas$Row.names
        betas$Row.names <- NULL
      }
    }
  }
  # in case reading failed, or we have only one column we won't continue
  data  <- NULL
  if(!is.null(ses) & !is.null(betas) & ncol(ses) > 1 & ncol(betas) > 1){
    # put it in a list
    ses <- ses[complete.cases(ses), ]
    betas <- betas[complete.cases(betas), ]
    data = mash_set_data(as.matrix(betas), as.matrix(ses))
  }
  # return the mash data set
  return(data)
}

get_conditions_mash_meta_z <- function(eqtl_output_loc, conditions = c(), cell_type){
  betas = NULL
  ses = NULL
  # go through the conditions
  for(condition in conditions){
    # there might be some tables not present, so let's set to a default
    table <- NULL
    # read the table
    table_loc <- paste(eqtl_output_loc, condition, '/', cell_type, '_expression/', 'eQTLsFDR-ProbeLevel.txt.gz',  sep = '')
    # show progress
    print(paste('reading', table_loc))
    try({
      table <- read.table(table_loc, sep = '\t', header = T)
    })
    if(is.null(table)){
      print(paste('table not present for:', table_loc))
    } else{
      # grab the variables
      table_b1 <- table$OverallZScore
      # add the SNP>probe as rownames
      table_b1 <- data.frame(table_b1, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
      # add the colnames
      colnames(table_b1) <- c(condition)
      # create the SEs
      table_se1 <- data.frame(c(rep(1, nrow(table_b1))), row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
      colnames(table_se1) <- c(condition)
      # banana
      if(is.null(ses)){
        ses <- table_se1
      } else{
        ses <- merge(ses, table_se1, by="row.names",all.x=TRUE)
        rownames(ses) <- ses$Row.names
        ses$Row.names <- NULL
      }
      if(is.null(betas)){
        betas <- table_b1
      } else{
        betas <- merge(betas, table_b1, by="row.names",all.x=TRUE)
        rownames(betas) <- betas$Row.names
        betas$Row.names <- NULL
      }
    }
  }
  # in case reading failed, or we have only one column we won't continue
  data  <- NULL
  if(!is.null(ses) & !is.null(betas) & ncol(ses) > 1 & ncol(betas) > 1){
    # put it in a list
    ses <- ses[complete.cases(ses), ]
    betas <- betas[complete.cases(betas), ]
    data = mash_set_data(as.matrix(betas), as.matrix(ses))
  }
  # return the mash data set
  return(data)
}


get_conditions_cell_types_mash <- function(eqtl_output_loc, result_output_loc, cell_types, use_z=F){
  for(cell_type in cell_types){
    # do 3h
    data.x3h <- NULL
    if(use_z){
      data.x3h <- get_conditions_mash_meta_z(eqtl_output_loc, conditions = c('3hCA', '3hPA', '3hMTB'), cell_type)
    } else{
      data.x3h <- get_conditions_mash(eqtl_output_loc, conditions = c('3hCA', '3hPA', '3hMTB'), cell_type)
    }
    # result might be NULL if there is not enough data to compare
    if(!is.null(data.x3h)){
      # run MASH
      m.1by1.x3h = mash_1by1(data.x3h)
      U.pca.x3h = cov_pca(data.x3h,ncol(data.x3h$Bhat)) # 6PCS due to number of conditions
      U.ed.x3h = cov_ed(data.x3h, U.pca.x3h)
      U.c.x3h = cov_canonical(data.x3h)
      m.x3h = mash(data.x3h, c(U.c.x3h,U.ed.x3h))
      lfsr.x3h <- m.x3h$result$lfsr
      # save result
      result_output_loc_3h <- paste(result_output_loc, cell_type, '_3h.tsv', sep = '')
      write.table(lfsr.x3h, result_output_loc_3h, sep = '\t', row.names = T, col.names = T)
    }

    # do 24h
    data.x24h <- NULL
    if(use_z){
      data.x24h <- get_conditions_mash_meta_z(eqtl_output_loc, conditions = c('24hCA', '24hPA', '24hMTB'), cell_type)
    } else{
      data.x24h <- get_conditions_mash(eqtl_output_loc, conditions = c('24hCA', '24hPA', '24hMTB'), cell_type)
    }
    # result might be NULL if there is not enough data to compare
    if(!is.null(data.x24h)){
      # run MASH
      m.1by1.x24h = mash_1by1(data.x24h)
      U.pca.x24h = cov_pca(data.x24h,ncol(data.x24h$Bhat)) # 6PCS due to number of conditions
      U.ed.x24h = cov_ed(data.x24h, U.pca.x24h)
      U.c.x24h = cov_canonical(data.x24h)
      m.x24h = mash(data.x24h, c(U.c.x24h,U.ed.x24h))
      lfsr.x24h <- m.x24h$result$lfsr
      # save result
      result_output_loc_24h <- paste(result_output_loc, cell_type, '_24h.tsv', sep = '')
      write.table(lfsr.x24h, result_output_loc_24h, sep = '\t', row.names = T, col.names = T)
    }
  }
}

get_conditions_cell_types_mash_combined <- function(eqtl_output_loc, result_output_loc, cell_types, use_z=F){
  for(cell_type in cell_types){
    # do 3h and 24h
    data <- NULL
    if(use_z){
      data <- get_conditions_mash_meta_z(eqtl_output_loc, conditions = c('3hCA', '3hPA', '3hMTB', '24hCA', '24hPA', '24hMTB'), cell_type)
    } else{
      data <- get_conditions_mash(eqtl_output_loc, conditions = c('3hCA', '3hPA', '3hMTB', '24hCA', '24hPA', '24hMTB'), cell_type)
    }
    # result might be NULL if there is not enough data to compare
    if(!is.null(data)){
      # run MASH
      m.1by1 = mash_1by1(data)
      U.pca = cov_pca(data,ncol(data$Bhat)) # 6PCS due to number of conditions
      U.ed = cov_ed(data, U.pca)
      U.c = cov_canonical(data)
      m = mash(data, c(U.c,U.ed))
      lfsr <- m$result$lfsr
      # save result
      result_output_loc_3and24h <- paste(result_output_loc, cell_type, '_3and24h.tsv', sep = '')
      write.table(lfsr, result_output_loc_3and24h, sep = '\t', row.names = T, col.names = T)
    }
  }
}


show_significant_overlaps <- function(table, sig_cutoff=0.05){
  # get the conditions from the column names
  conditions <- colnames(table)
  # initialize a list to store the values per condition
  lists <- list()
  for(condition in conditions){
    # grab the significant results
    sigs <- rownames(table[table[[condition]] <= sig_cutoff,])
    # put it in the list
    lists[[condition]] <- sigs
  }
  # put it in the actual diagram
  upset(fromList(lists), order.by = 'freq', nsets = ncol(table))
}


show_significant_overlaps_cell_types <- function(cell_types, result_loc, plot_loc, sig_cutoff=0.05, timepoints = c('3h', '24h')){
  # make plot for the cell types
  for(cell_type in cell_types){
    # make plot for that cell type per timepoint
    for(timepoint in timepoints){
      # get the location of the result file
      result_loc_ct_tp <- paste(result_loc, cell_type, '_', timepoint, '.tsv', sep = '')
      # read the table into memory
      table <- read.table(result_loc_ct_tp, sep = '\t', header = T, row.names = 1)
      # setup the place to save the plot
      plot_loc_ct_tp <- paste(plot_loc, cell_type, '_', timepoint, '.png', sep = '')
      # setup saving
      png(plot_loc_ct_tp, width = 2000, height = 2000, pointsize = 20)
      # call method that makes image
      show_significant_overlaps(table, sig_cutoff = sig_cutoff)
      # close connection
      #dev.off()
      # print progress
      print(paste('created plot at', plot_loc_ct_tp))
    }
  }
}

####################
# main code        #
####################

cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'mDC', 'megakaryocyte', 'monocyte', 'NK', 'pDC', 'plasma_B', 'hemapoietic_stem')
eqtl_result_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_lores_confine_1m_ut_all_cell_types/results/'
result_output_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/mashr/sct_mqc_lores_confine_1m_ut_all_cell_types/data_and_canonical/'
get_conditions_cell_types_mash(eqtl_result_loc, result_output_loc, cell_types)
result_output_loc_meta_z <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/mashr/sct_mqc_lores_confine_1m_ut_all_cell_types/data_and_canonical_meta_z/'
get_conditions_cell_types_mash(eqtl_result_loc, result_output_loc_meta_z, cell_types, use_z = T)
get_conditions_cell_types_mash_combined(eqtl_result_loc, result_output_loc_meta_z, cell_types, use_z = T)


# make plots
plot_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/mashr/figures/sct_mqc_lores_confine_1m_ut_all_cell_types/data_and_canonical/'
show_significant_overlaps_cell_types(cell_types, result_loc = result_output_loc, plot_loc = plot_loc)
plot_loc_meta_z <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/mashr/figures/sct_mqc_lores_confine_1m_ut_all_cell_types/data_and_canonical_meta_z/'
show_significant_overlaps_cell_types(cell_types, result_loc = result_output_loc_meta_z, plot_loc = plot_loc_meta_z)
