############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_MAST_meta.R
# Function: meta analyse output of MAST v2 and v3 chemistries, and determine significant DE genes
############################################################################################################################


######################
# libraries          #
######################

library(metap)
library(MetaVolcanoR)
library(stringr)
library(data.table)
require("heatmap.plus")
library(RColorBrewer)
library(VennDiagram)
library(Matrix)

####################
# Functions        #
####################

write_meta_mast <- function(condition_info, mast_output_loc_prepend, mast_output_loc_append, mast_meta_output_loc_prepend){
  # go through the conditions
  for(condition in c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
    # check for each cell type
    for(cell_type in c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'NK', 'hemapoietic stem', 'megakaryocyte', 'monocyte', 'plasma B')){
      # get the number of cells
      #cond1_v2_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition1 == 'UT' & condition_info$chem == 'V2', ]$nr_of_cells_condition1[1]
      #cond1_v3_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition1 == 'UT' & condition_info$chem == 'V3', ]$nr_of_cells_condition1[1]
      #cond2_v2_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition2 == condition & condition_info$chem == 'V2', ]$nr_of_cells_condition2[1]
      #cond2_v3_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition2 == condition & condition_info$chem == 'V3', ]$nr_of_cells_condition2[1]
      # get the mast output
      mast_loc_v2 <- paste(mast_output_loc_prepend, '2', mast_output_loc_append, cell_type, 'UT', condition, '.tsv', sep = '')
      mast_loc_v3 <- paste(mast_output_loc_prepend, '3', mast_output_loc_append, cell_type, 'UT', condition, '.tsv', sep = '')
      try({
        # read the mast output
        mast_v2 <- read.table(mast_loc_v2, header=T)
        mast_v3 <- read.table(mast_loc_v3, header=T)
        # get the genes that are in both
        genes_both <- intersect(rownames(mast_v2), rownames(mast_v3))
        # select only those genes
        mast_v2 <- mast_v2[rownames(mast_v2) %in% genes_both,]
        mast_v3 <- mast_v3[rownames(mast_v3) %in% genes_both,]
        # morph P val to minimum
        if(nrow(mast_v2[mast_v2$p_val == 0, ]) > 0){
          mast_v2[mast_v2$p_val == 0, ]$p_val <- .Machine$double.xmin
        }
        if(nrow(mast_v3[mast_v3$p_val == 0, ]) > 0){
          mast_v3[mast_v3$p_val == 0, ]$p_val <- .Machine$double.xmin
        }
        # add the gene name also in a column
        mast_v2$gene <- rownames(mast_v2)
        mast_v3$gene <- rownames(mast_v3)
        # add the mast results
        masts <- list()
        masts$v2 <- mast_v2
        masts$v3 <- mast_v3
        # perform the metavolcanor approach
        meta_degs_comb <- combining_mv(diffexp=masts, pcriteria='p_val', foldchangecol='avg.logFC', genenamecol = 'gene', collaps = T)
        # grab the result we care about
        volcanometa <- meta_degs_comb@metaresult
        volcanometa$metap_bonferroni <- volcanometa$metap*length(genes_both)
        # add the genes as rownames
        rownames(volcanometa) <- volcanometa$gene
        # add a colname append
        colnames(mast_v2) <- paste(colnames(mast_v2), 'v2', sep = '_')
        colnames(mast_v3) <- paste(colnames(mast_v3), 'v3', sep = '_')
        # merge the frames
        mast <- merge(mast_v2, mast_v3, by=0, all=TRUE)
        rownames(mast) <- mast$Row.names
        mast$Row.names <- NULL
        # get the meta p values using stouffers method
        #stouffers <- rep(NA, times = nrow(mast))
        #for(i in 1:nrow(mast)){
        #  # get the p-values
        #  p_vals <- c(mast[i, 'p_val_v2'], mast[i, 'p_val_v3'])
        #  # the weights are based on the number of cells
        #  weights <- c(sqrt(cond1_v2_cells + cond2_v2_cells), sqrt(cond1_v3_cells + cond2_v3_cells))
        #  # get the result from the Stouffer's method
        #  stouffers_res <- sumz(p = p_vals, weights = weights)
        #  if(!is.na(stouffers_res)){
        #    stouffers[i] <- stouffers_res$p[1,1]*length(genes_both) #bonferroni correct by multiplying by number of tests
        #  }
        #}
        # add the value
        #mast$stouffers_p <- stouffers
        #mast[mast$stouffers_p > 1 & !is.na(mast$stouffers_p), ]$stouffers_p <- 1
        # also add the volcanometa stuff
        mast <- merge(mast, volcanometa, by=0, all=TRUE)
        rownames(mast) <- mast$Row.names
        mast$Row.names <- NULL
        if(nrow(mast[mast$metap_bonferroni > 1, ]) > 0){
          mast[mast$metap_bonferroni > 1, ]$metap_bonferroni <- 1
        }
        # write the result
        output_loc <- paste(mast_meta_output_loc_prepend, cell_type, 'UT', condition, '.tsv', sep = '')
        write.table(mast, output_loc, sep = '\t')
      })
    }
  }
}

get_unique_condition_DE_genes <- function(mast_meta_output_loc, unique_output_loc, stims=c('CA', 'MTB', 'PA'), timepoints=c('3h', '24h'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), pval_column='metap_bonferroni', sig_pval=0.05, lfc_column='metafc', only_positive=F, only_negative=F){
  # check each cell type
  for(cell_type in cell_types){
    # check for each timepoint
    for(timepoint in timepoints){
      stim_de_genes <- list()
      # combine with each stim
      for(stim in stims){
        # create the output loc string
        output_loc <- paste(mast_meta_output_loc, cell_type, 'UTX', timepoint, stim, '.tsv', sep = '')
        # read the output
        mast <- read.table(output_loc, header=T)
        # filter to only include the significant results
        mast <- mast[mast[[pval_column]] <= sig_pval, ]
        # filter for only the positive lfc if required
        if(only_positive){
          mast <- mast[mast[[lfc_column]] < 0, ]
        }
        # filter for only the positive lfc if required
        if(only_negative){
          mast <- mast[mast[[lfc_column]] > 0, ]
        }
        # the genes are the rows
        de_genes <- rownames(mast)
        # add to the list
        stim_de_genes[[paste('UTX', timepoint, stim, sep = '')]] <- de_genes
      }
      # now check the combinations
      for(condition in names(stim_de_genes)){
        # make a copy of the list with all conditions
        stim_de_genes_copy <- stim_de_genes
        # grab the genes for this condition
        de_genes_condition <- stim_de_genes_copy[[condition]]
        # remove these genes from the copied list
        stim_de_genes_copy[[condition]] <- NULL
        # now grab what is left in the list, so anything that was not this condition
        de_genes_other_conditions <- as.vector(unlist(stim_de_genes_copy))
        # do a setdiff to get what is left
        de_unique <- setdiff(de_genes_condition, de_genes_other_conditions)
        # write these out
        output_loc_full <- paste(unique_output_loc, cell_type, condition, '_unique_to_pathogen_across_timepoint.txt', sep = '')
        write.table(data.frame(de_unique), output_loc_full, row.names = F, col.names = F, quote = F)
      }
    }
  }
}


get_shared_condition_DE_genes <- function(mast_meta_output_loc, shared_output_loc, stims=c('CA', 'MTB', 'PA'), timepoints=c('3h', '24h'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), pval_column='metap_bonferroni', sig_pval=0.05, lfc_column='metafc', only_positive=F, only_negative=F){
  # check each cell type
  for(cell_type in cell_types){
    # check for each timepoint
    for(timepoint in timepoints){
      stim_de_genes <- NULL
      # combine with each stim
      for(stim in stims){
        # create the output loc string
        output_loc <- paste(mast_meta_output_loc, cell_type, 'UTX', timepoint, stim, '.tsv', sep = '')
        # read the output
        mast <- read.table(output_loc, header=T)
        # filter to only include the significant results
        mast <- mast[mast[[pval_column]] <= sig_pval, ]
        # filter for only the positive lfc if required
        if(only_positive){
          mast <- mast[mast[[lfc_column]] < 0, ]
        }
        # filter for only the positive lfc if required
        if(only_negative){
          mast <- mast[mast[[lfc_column]] > 0, ]
        }
        # the genes are the rows
        de_genes <- rownames(mast)
        # add to the list
        if(is.null(stim_de_genes)){
          stim_de_genes <- de_genes
        }
        else{
          stim_de_genes <- intersect(stim_de_genes, de_genes)
        }
      }
      # write these out
      output_loc_full <- paste(shared_output_loc, cell_type, timepoint, '_shared_at_timepoint_across_pathogens.txt', sep = '')
      write.table(data.frame(stim_de_genes), output_loc_full, row.names = F, col.names = F, quote = F)
    }
  }
}

get_unique_timepoint_DE_genes <- function(mast_meta_output_loc, unique_output_loc, stims=c('CA', 'MTB', 'PA'), timepoints=c('3h', '24h'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), pval_column='metap_bonferroni', sig_pval=0.05, lfc_column='metafc', only_positive=F, only_negative=F){
  # check each cell type
  for(cell_type in cell_types){
    # check for each timepoint
    for(stim in stims){
      timepoint_de_genes <- list()
      # combine with each stim
      for(timepoint in timepoints){
        # create the output loc string
        output_loc <- paste(mast_meta_output_loc, cell_type, 'UTX', timepoint, stim, '.tsv', sep = '')
        # read the output
        mast <- read.table(output_loc, header=T)
        # filter to only include the significant results
        mast <- mast[mast[[pval_column]] <= sig_pval, ]
        # filter for only the positive lfc if required
        if(only_positive){
          mast <- mast[mast[[lfc_column]] < 0, ]
        }
        # filter for only the positive lfc if required
        if(only_negative){
          mast <- mast[mast[[lfc_column]] > 0, ]
        }
        # the genes are the rows
        de_genes <- rownames(mast)
        # add to the list
        timepoint_de_genes[[paste('UTX', timepoint, stim, sep = '')]] <- de_genes
      }
      # now check the combinations
      for(timepoint in names(timepoint_de_genes)){
        # make a copy of the list with all conditions
        timepoint_de_genes_copy <- timepoint_de_genes
        # grab the genes for this condition
        de_genes_timepoint <- timepoint_de_genes_copy[[timepoint]]
        # remove these genes from the copied list
        timepoint_de_genes_copy[[timepoint]] <- NULL
        # now grab what is left in the list, so anything that was not this condition
        de_genes_other_timepoints <- as.vector(unlist(timepoint_de_genes_copy))
        # do a setdiff to get what is left
        de_unique <- setdiff(de_genes_timepoint, de_genes_other_timepoints)
        # write these out
        output_loc_full <- paste(unique_output_loc, cell_type, timepoint, '_unique_to_timepoint_same_pathogen.txt', sep = '')
        write.table(data.frame(de_unique), output_loc_full, row.names = F, col.names = F, quote = F)
      }
    }
  }
}

get_shared_timepoint_DE_genes <- function(mast_meta_output_loc, shared_output_loc, stims=c('CA', 'MTB', 'PA'), timepoints=c('3h', '24h'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), pval_column='metap_bonferroni', sig_pval=0.05, lfc_column='metafc', only_positive=F, only_negative=F){
  # check each cell type
  for(cell_type in cell_types){
    # check for each timepoint
    for(stim in stims){
      timepoint_de_genes <- NULL
      # combine with each stim
      for(timepoint in timepoints){
        # create the output loc string
        output_loc <- paste(mast_meta_output_loc, cell_type, 'UTX', timepoint, stim, '.tsv', sep = '')
        # read the output
        mast <- read.table(output_loc, header=T)
        # filter to only include the significant results
        mast <- mast[mast[[pval_column]] <= sig_pval, ]
        # filter for only the positive lfc if required
        if(only_positive){
          mast <- mast[mast[[lfc_column]] < 0, ]
        }
        # filter for only the positive lfc if required
        if(only_negative){
          mast <- mast[mast[[lfc_column]] > 0, ]
        }
        # the genes are the rows
        de_genes <- rownames(mast)
        # add to the list
        if(is.null(timepoint_de_genes)){
          timepoint_de_genes <- de_genes
        }
        else{
          timepoint_de_genes <- intersect(timepoint_de_genes, de_genes)
        }
      }
      # write these out
      output_loc_full <- paste(shared_output_loc, cell_type, stim, '_shared_across_timepoints_same_pathogen.txt', sep = '')
      write.table(data.frame(timepoint_de_genes), output_loc_full, row.names = F, col.names = F, quote = F)
    }
  }
}



write_meta_mast_3hvs24h <- function(condition_info, mast_output_loc_prepend, mast_output_loc_append, mast_meta_output_loc_prepend){
  # go through the conditions
  for(condition in c('CA', 'PA', 'MTB')){
    # check for each cell type
    for(cell_type in c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'hemapoietic stem', 'megakaryocyte', 'monocyte', 'plasma B')){
      # get the mast output
      mast_loc_v2 <- paste(mast_output_loc_prepend, '2', mast_output_loc_append, cell_type, 'X3h', condition, 'X24h', condition, '.tsv', sep = '')
      mast_loc_v3 <- paste(mast_output_loc_prepend, '3', mast_output_loc_append, cell_type, 'X3h', condition, 'X24h', condition, '.tsv', sep = '')
      try({
        # read the mast output
        mast_v2 <- read.table(mast_loc_v2, header=T)
        mast_v3 <- read.table(mast_loc_v3, header=T)
        # get the genes that are in both
        genes_both <- intersect(rownames(mast_v2), rownames(mast_v3))
        # select only those genes
        mast_v2 <- mast_v2[rownames(mast_v2) %in% genes_both,]
        mast_v3 <- mast_v3[rownames(mast_v3) %in% genes_both,]
        # morph P val to minimum
        if(nrow(mast_v2[mast_v2$p_val == 0, ]) > 0){
          mast_v2[mast_v2$p_val == 0, ]$p_val <- .Machine$double.xmin
        }
        if(nrow(mast_v3[mast_v3$p_val == 0, ]) > 0){
          mast_v3[mast_v3$p_val == 0, ]$p_val <- .Machine$double.xmin
        }
        # add the gene name also in a column
        mast_v2$gene <- rownames(mast_v2)
        mast_v3$gene <- rownames(mast_v3)
        # add the mast results
        masts <- list()
        masts$v2 <- mast_v2
        masts$v3 <- mast_v3
        # perform the metavolcanor approach
        meta_degs_comb <- combining_mv(diffexp=masts, pcriteria='p_val', foldchangecol='avg.logFC', genenamecol = 'gene', collaps = T)
        # grab the result we care about
        volcanometa <- meta_degs_comb@metaresult
        volcanometa$metap_bonferroni <- volcanometa$metap*length(genes_both)
        # add the genes as rownames
        rownames(volcanometa) <- volcanometa$gene
        # add a colname append
        colnames(mast_v2) <- paste(colnames(mast_v2), 'v2', sep = '_')
        colnames(mast_v3) <- paste(colnames(mast_v3), 'v3', sep = '_')
        # merge the frames
        mast <- merge(mast_v2, mast_v3, by=0, all=TRUE)
        rownames(mast) <- mast$Row.names
        mast$Row.names <- NULL
        # also add the volcanometa stuff
        mast <- merge(mast, volcanometa, by=0, all=TRUE)
        rownames(mast) <- mast$Row.names
        mast$Row.names <- NULL
        if(nrow(mast[mast$metap_bonferroni > 1, ]) > 0){
          mast[mast$metap_bonferroni > 1, ]$metap_bonferroni <- 1
        }
        # write the result
        output_loc <- paste(mast_meta_output_loc_prepend, cell_type, 'X3h', condition, 'X24h', condition,'.tsv', sep = '')
        write.table(mast, output_loc, sep = '\t')
      })
    }
  }
}

write_background_meta <- function(background_loc_prepend, background_loc_append, meta_output_loc, to_ens = F, symbols.to.ensg.mapping = 'genes.tsv'){
  # go through the conditions
  for(condition in c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
    # check for each cell type
    for(cell_type in c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'hemapoietic stem', 'megakaryocyte', 'monocyte', 'plasma B')){
      # set the paths
      background_v2_loc <- paste(background_loc_prepend, '2/', cell_type, '_UT_', condition, background_loc_append, sep = '')
      background_v3_loc <- paste(background_loc_prepend, '3/', cell_type, '_UT_', condition, background_loc_append, sep = '')
      # read the genes
      v2_genes <- read.table(background_v2_loc, header = F)
      v3_genes <- read.table(background_v3_loc, header = F)
      # intersect to get the genes that are in both
      meta_genes <- intersect(v2_genes$V1, v3_genes$V1)
      # convert the symbols to ensemble IDs
      if (to_ens) {
        mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
        mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
        meta_genes <- mapping[match(meta_genes, mapping$V2),"V1"]
      }
      # otherwise change the Seurat replacement back
      else{
        #meta_genes <- gsub("-", "_", meta_genes)
      }
      # write the result
      background_meta_loc <- paste(meta_output_loc, '/', cell_type, '_UT_', condition, background_loc_append, sep = '')
      write.table(meta_genes, background_meta_loc, col.names = F, row.names = F, quote = F)
    }
  }
}

get_significant_genes <- function(mast_output_loc, sig_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv'){
  # get the files
  files <- list.files(mast_output_loc)
  # try to read each file
  for(file in files){
    try({
      # read the mast output
      mast <- read.table(paste(mast_output_loc, file, sep = ''), header=T)
      # filter to only include the significant results
      mast <- mast[mast[[pval_column]] <= sig_pval, ]
      # filter for only the positive lfc if required
      if(only_positive){
        mast <- mast[mast[[lfc_column]] < 0, ]
      }
      # filter for only the positive lfc if required
      if(only_negative){
        mast <- mast[mast[[lfc_column]] > 0, ]
      }
      # confine in some way if reporting a max number of genes
      if(!is.null(max)){
        # by p if required
        if(max_by_pval){
          mast <- mast[order(mast[[p_val_column]]), ]
        }
        # by lfc otherwise
        else{
          mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
        }
        # subset to the number we requested if max was set
        mast <- mast[1:max,]
      }
      # grab the genes from the column names
      genes <- rownames(mast)
      # convert the symbols to ensemble IDs
      if (to_ens) {
        mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
        mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
        genes <- mapping[match(genes, mapping$V2),"V1"]
      }
      # otherwise change the Seurat replacement back
      else{
        #genes <- gsub("-", "_", genes)
      }
      # create a regex to get the last index of .
      last_dot_pos <- "\\.[^\\.]*$"
      # this allows us to remove the filename extention
      file_no_ext <- substr(file, 1, regexpr(last_dot_pos,file)-1)
      # create output location
      sig_output <- paste(sig_output_loc, file_no_ext, '.txt', sep = '')
      # write the genes
      write.table(genes, sig_output, sep = '\t', quote = F, row.names = F, col.names = F)
    })
  }
}

get_pathway_table <- function(pathway_output_loc, sig_val_to_use = 'q.value.Bonferroni', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB'), use_ranking=F){
  # put all results in a list
  pathways_analysis <- list()
  # put all results in a shared DF
  pathway_df <- NULL
  # check each cell type
  for(cell_type in cell_types){
    # check each stim
    for(stim in stims){
      try({
        print(paste(cell_type, stim, sep = ' '))
        # paste the filepath together
        #filepath <- paste(pathway_output_loc, cell_type, 'UT',stim,'_sig_pathways.txt', sep = '')
        filepath <- paste(pathway_output_loc, cell_type, 'UT',stim,'_sig_up_pathways.txt', sep = '')
        # read the file
        pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
        # create column name
        newcolname <- paste(cell_type, 'UT', stim, sep = '')
        # get the log2 of the significance value
        #pathways[[newcolname]] <- log2(pathways[[sig_val_to_use]])
        if(use_ranking){
          pathways[[newcolname]] <- as.numeric(rownames(pathways))
        }
        else{
          pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
        }
        pathways$id_name <- paste(pathways$ID, pathways$Name, sep = '_')
        # reduce to the only two columns we care about
        pathways <- pathways[, c('id_name', newcolname)]
        # join with other pathway files
        if(is.null(pathway_df)){
          # just set as df if the first round through
          pathway_df <- pathways
          pathway_df <- data.table(pathway_df, key = c('id_name'))
        }
        else{
          # otherwise, merge with existing pathways
          pathway_df <- merge(pathway_df, data.table(pathways, key = c('id_name')), by.x='id_name', by.y='id_name', all=T)
          #pathway_df[[newcolname]] <- pathways[[newcolname]][match(pathway_df$Name, pathways$Name)]
          #pathway_df <- left_join(pathway_df, pathways)

        }
      })
    }
  }
  # turn into regular df
  pathway_df <- setDF(pathway_df)
  # set all NA to zero
  pathway_df[is.na(pathway_df)] <- 0
  # set rownames
  rownames(pathway_df) <- pathway_df$id_name
  pathway_df$id_name <- NULL
  # remove rows which amount to zero
  pathway_df <- pathway_df[apply(pathway_df[,-1], 1, function(x) !all(x==0)),]
  return(pathway_df)
}


get_top_pathways <- function(pathway_table, nr_of_top_genes, is_ranked=F){
  # init pathways list
  pathways <- c()
  # go through the columns
  for(col in colnames(pathway_table)){
    # order by that column
    ordered <- pathway_table[order(pathway_table[[col]], decreasing = T), ]
    if(is_ranked){
      ordered <- pathway_table[order(pathway_table[[col]], decreasing = F), ]
    }
    # get those top ones
    top_col <- rownames(ordered)[1:nr_of_top_genes]
    pathways <- c(pathways, top_col)
  }
  # limit to those top pathways now
  pathway_table_smaller <- pathway_table[rownames(pathway_table) %in% pathways, ]
  return(pathway_table_smaller)
}


get_most_shared_pathways <- function(pathway_table, top_so_many=10, use_sd_method=F){
  most_shared <- c()
  # there are different methods to get the most shared pathways
  if(use_sd_method){
    # the SD methods just gets the pathways with the lowest standard deviation
    most_shared <- get_most_varying_from_df(pathway_table, top_so_many = top_so_many, dont_get_least_varying = F)
  }
  else{
    # first get the sum of the rankings over the rows, lower overall ranking means more shared
    summed_rank <- apply(pathway_table, 1, sum)
    # get the top X so many from the pathway table, ordered by this sum of rankings
    most_shared <- rownames(pathway_table[order(summed_rank), ])[1:top_so_many]

    # get the sum of the 3h conditions
    timepoints_3h <- colnames(pathway_table)[grep('3h', colnames(pathway_table))]
    summed_rank_3h <- apply(pathway_table[, timepoints_3h], 1, sum)
    most_shared <- c(most_shared, rownames(pathway_table[order(summed_rank_3h), ])[1:top_so_many])

    # get the sum of the 24h conditions
    timepoints_24h <- colnames(pathway_table)[grep('24h', colnames(pathway_table))]
    summed_rank_24h <- apply(pathway_table[, timepoints_24h], 1, sum)
    most_shared <- c(most_shared, rownames(pathway_table[order(summed_rank_24h), ])[1:top_so_many])

    # make unique of course
    most_shared <- unique(most_shared)

  }
  return(most_shared)
}

get_most_varied_pathways <- function(pathway_table, top_so_many=10, use_sd_method=F, disregard_3h_vs_24h=F){
  most_varied <- c()
  # most varied in 3h condition
  timepoints_3h <- colnames(pathway_table)[grep('3h', colnames(pathway_table))]
  # most varied in 24h condition
  timepoints_24h <- colnames(pathway_table)[grep('24h', colnames(pathway_table))]
  # most varied in PA condition
  timepoints_pa <- colnames(pathway_table)[grep('hPA', colnames(pathway_table))]
  # most varied in PA condition
  timepoints_ca <- colnames(pathway_table)[grep('hCA', colnames(pathway_table))]
  # most varied in PA condition
  timepoints_mtb <- colnames(pathway_table)[grep('hMTB', colnames(pathway_table))]
  if(use_sd_method){
    # first overall most variation
    most_varied <- c()
    # because 3h has sharing and 24h has sharing, looking overall, we might see 3h vs 24h effects, we can choose to ignore this
    if(disregard_3h_vs_24h == F){
      most_varied <- get_most_varying_from_df(pathway_table, top_so_many = top_so_many)
    }
    # most varied in 3h and most varied in 24h
    most_varied <- c(most_varied, get_most_varying_from_df(pathway_table[, timepoints_3h], top_so_many = top_so_many))
    most_varied <- c(most_varied, get_most_varying_from_df(pathway_table[, timepoints_24h], top_so_many = top_so_many))
  }
  else{
    # large difference between 3h and 24h
    mean_3h_ranks <- apply(pathway_table[, timepoints_3h], 1, mean)
    mean_24h_ranks <- apply(pathway_table[, timepoints_24h], 1, mean)
    # absolute difference
    mean_rank_diff <- abs(mean_3h_ranks - mean_24h_ranks)
    # grab by varied over abs mean difference
    most_varied <- rownames(pathway_table[order(mean_rank_diff, decreasing = T), ])[1:top_so_many]
    # get pathogen mean ranks
    mean_pa_ranks <- apply(pathway_table[, timepoints_pa], 1, mean)
    mean_ca_ranks <- apply(pathway_table[, timepoints_ca], 1, mean)
    mean_mtb_ranks <- apply(pathway_table[, timepoints_mtb], 1, mean)
    # turn into separate df
    path_df <- data.frame(pa=mean_pa_ranks, ca=mean_ca_ranks, mtb=mean_mtb_ranks)
    rownames(path_df) <- rownames(pathway_table)
    # use sd method
    most_varied <- c(most_varied, get_most_varying_from_df(path_df, top_so_many = top_so_many))
  }
  most_varied <- unique(most_varied)
  return(most_varied)
}

get_mast_meta_output_overlap <- function(mast_meta_output_1, mast_meta_output_2, venn_output_loc='./', only_significant=T, pval_column = 'metap_bonferroni', group1name='group 1', group2name='group 2'){
  # list the files in directory 1
  files <- list.files(mast_meta_output_1)
  # check the files
  for(file in files){
    # set the full path
    mast1_loc <- paste(mast_meta_output_1, '/', file, sep='')
    mast2_loc <- paste(mast_meta_output_2, '/', file, sep='')
    try({
      # try to read both tables
      mast1 <- read.table(mast1_loc, header=T, row.names = 1, sep = '\t')
      mast2 <- read.table(mast2_loc, header=T, row.names = 1, sep = '\t')
      # remove insignificant results if requested
      if(only_significant){
        mast1 <- mast1[mast1[[pval_column]] <= 0.05, ]
        mast2 <- mast2[mast2[[pval_column]] <= 0.05, ]
      }
      # report something about the P-values
      quantile_mast1 <- quantile(mast1[[pval_column]])
      quantile_mast2 <- quantile(mast2[[pval_column]])
      print(paste(substr(file, 1, regexpr("\\.[^\\.]*$", file)-1), group1name, 'quantiles:'))
      print(quantile_mast1)
      print(paste(substr(file, 1, regexpr("\\.[^\\.]*$", file)-1), group2name, 'quantiles:'))
      print(quantile_mast2)
      # grab the genes
      mast1_genes <- rownames(mast1)
      mast2_genes <- rownames(mast2)
      # grab the name of the file without the extention
      myCol <- brewer.pal(3, "Pastel2")[1:2]
        output_file <- substr(file, 1, regexpr("\\.[^\\.]*$", file)-1)
        venn.diagram(x = list(mast1_genes, mast2_genes),
                     main = substr(file, 1, regexpr("\\.[^\\.]*$", file)-1),
                     category.names = c(group1name, group2name),
                     filename = paste(venn_output_loc, output_file, '.png', sep = ''),
                     imagetype="png" ,
                     height = 600 ,
                     width = 600 ,
                     resolution = 300,
                     compression = "lzw",
                     lwd = 2,
                     lty = 'blank',
                     fill = myCol,
                     cex = .6,
                     fontface = "bold",
                     fontfamily = "sans",
                     cat.cex = 0.6,
                     cat.fontface = "bold",
                     cat.default.pos = "outer",
                     cat.pos = c(-27, 27),
                     cat.dist = c(0.055, 0.055),
                     cat.fontfamily = "sans")
    })
  }
}


get_combined_meta_de_table <- function(meta_output_loc, must_be_positive_once=F, convert_insignificant_p_to_lfc0=F, pval_column='metap_bonferroni', lfc_column='metafc', pval_significance_threshold=0.05, pathogens=c("CA", "MTB", "PA"),timepoints=c("3h", "24h"), cell_types_to_use=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")){
  deg_meta_combined <- NULL
  for(pathogen in pathogens) {
    for (timepoint in timepoints) {
      for (cell_type in cell_types_to_use) {
        deg_table <- read.table(paste0(meta_output_loc, cell_type, "UTX", timepoint, pathogen, ".tsv"), stringsAsFactors = F, sep = "\t")
        # subset to only the significant ones if so selected
        if(convert_insignificant_p_to_lfc0){
          deg_table <- deg_table[deg_table[[pval_column]] < pval_significance_threshold, ]
        }
        deg_table <- deg_table[lfc_column]
        colnames(deg_table) <- c(paste(cell_type, "UTX", timepoint, pathogen, sep=''))
        deg_table$genes <- rownames(deg_table)
        deg_table <- data.table(deg_table)
        print(head(deg_table))
        if(is.null(deg_meta_combined)){
          deg_meta_combined <- deg_table
        }
        else{
          deg_meta_combined <- merge(deg_meta_combined, deg_table, by.x='genes', by.y='genes', all=TRUE)
          #deg_meta_combined[,paste(cell_type, timepoint, pathogen, sep = "_")] <- deg_table$metafc
        }

      }
    }
  }

  deg_meta_combined <- data.frame(deg_meta_combined)
  rownames(deg_meta_combined) <- deg_meta_combined$genes
  deg_meta_combined$genes <- NULL
  deg_meta_combined[is.na(deg_meta_combined)] <- 0
  # limit to those that were upregulated at least once if requested
  if(must_be_positive_once){
    deg_meta_combined <- deg_meta_combined[apply(deg_meta_combined,1,min) < 0,]
  }
  return(deg_meta_combined)
}

plot_loose_overlap_cell_types <- function(mast_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB'), pval_column='metap_bonferroni', lfc_column='metafc', pval_significance_threshold=0.05, pval_column_other='metap_bonferroni', pval_significance_threshold_other=0.1){
  # get the overlaps
  loose_overlaps <- get_loose_overlap_cell_types(mast_output_loc, cell_types=cell_types, stims=stims, pval_column=pval_column, lfc_column=lfc_column, pval_significance_threshold=pval_significance_threshold, pval_column_other=pval_column_other, pval_significance_threshold_other=pval_significance_threshold_other)
  # create the plot df
  plot_df <- NULL
  # go through the non-unique ones
  for(number in unique(loose_overlaps$count[loose_overlaps$count > 1])){
    # get the number of genes with that number of overlap
    number_of_genes <- nrow(loose_overlaps[loose_overlaps$count == number, ])
    # create the row
    plot_row <- data.frame(number=c(number), number_of_genes=c(number_of_genes), cell_type=c('mixed'), stringsAsFactors = F)
    # add to the dataframe
    if(is.null(plot_df)){
      plot_df <- plot_row
    }
    else{
      plot_df <- rbind(plot_df, plot_row)
    }
  }
  # add the unique ones as well
  for(ct in unique(loose_overlaps[loose_overlaps$count == 1, 'specificity'])){
    # get number of genes specific to that cell type
    number_of_genes <- nrow(loose_overlaps[loose_overlaps$specificity == ct, ])
    # create the row
    plot_row <- data.frame(number=c(1), number_of_genes=c(number_of_genes), cell_type=c(ct), stringsAsFactors = F)
    # add to the dataframe
    if(is.null(plot_df)){
      plot_df <- plot_row
    }
    else{
      plot_df <- rbind(plot_df, plot_row)
    }
  }
  # grab the colours
  cc <- get_color_coding_dict()
  # add the 'mixed' condition
  cc[['mixed']] <- 'gray'
  fillScale <- scale_fill_manual(name = "cell type",values = unlist(cc[plot_df$cell_type]))
  # make the plot finally
  p <- ggplot(plot_df, aes(fill=cell_type, y=number_of_genes, x=number)) +
    geom_bar(position='stack', stat='identity') +
    labs(x='number', y='number of DE genes') +
    ggtitle(paste('DE genes and cell type specificity')) +
    labs(fill = "Found in") +
    fillScale
  return(p)
}

get_loose_overlap_cell_types <- function(mast_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB'), pval_column='metap_bonferroni', lfc_column='metafc', pval_significance_threshold=0.05, pval_column_other='metap_bonferroni', pval_significance_threshold_other=0.1){
  # get the output
  all_output <- get_mast_per_cell_type_and_condition(mast_output_loc = mast_output_loc, cell_types = cell_types, stims = stims)
  # get the truly significant and loosely significant entries
  sig_per_ct_truly <- get_significant_genes_per_cell_type_and_condition(all_output, pval_column=pval_column, pval_significance_threshold=pval_significance_threshold, cell_types=cell_types, stims=stims)
  sig_per_ct_loose <- get_significant_genes_per_cell_type_and_condition(all_output, pval_column=pval_column_other, pval_significance_threshold=pval_significance_threshold_other, cell_types=cell_types, stims=stims)
  # reduce the conditions
  sig_per_ct_truly_red <- list()
  sig_per_ct_loose_red <- list()
  for(cell_type in intersect(cell_types, names(sig_per_ct_truly))){
    # create a vector to put the genes in, regardless of condition
    truly_vector_all_conds <- c()
    loose_vector_all_conds <- c()
    # add for this condition
    for(condition in intersect(stims, names(sig_per_ct_loose[[cell_type]]))){
      truly_vector_all_conds <- c(truly_vector_all_conds, sig_per_ct_truly[[cell_type]][[condition]])
      loose_vector_all_conds <- c(loose_vector_all_conds, sig_per_ct_loose[[cell_type]][[condition]])
    }
    # make unique
    truly_vector_all_conds <- unique(truly_vector_all_conds)
    loose_vector_all_conds <- unique(loose_vector_all_conds)
    # add to the lists
    sig_per_ct_truly_red[[cell_type]] <- truly_vector_all_conds
    sig_per_ct_loose_red[[cell_type]] <- loose_vector_all_conds
  }
  # initialize the dataframe
  count_df <- NULL
  # now check each gene, and see if it is in any other sets
  for(cell_type in names(sig_per_ct_truly_red)){
    # we only need to check genes we havent't already done
    for(gene in setdiff(sig_per_ct_truly_red[[cell_type]], count_df$gene)){
      # we start counting at 1
      count <- 1
      # check the other cell types
      for(cell_type_2 in setdiff(names(sig_per_ct_loose_red), c(cell_type))){
        # it it's also loosely found, increase the count
        if(gene %in% sig_per_ct_loose_red[[cell_type_2]]){
          count <- count + 1
        }
      }
      # create the row
      count_row <- data.frame(gene=c(gene), count=c(count), stringsAsFactors = F)
      # depending on if the count is one, set the type
      if(count == 1){
        count_row$specificity <- cell_type
      }
      else if(count > 1){
        count_row$specificity <- 'mixed'
      }
      # add the count for this gene
      if(is.null(count_df)){
        count_df <- count_row
      }
      else{
        count_df <- rbind(count_df, count_row)
      }
    }
  }
  return(count_df)
}


get_significant_genes_per_cell_type_and_condition <- function(output_per_ct_and_condition, pval_column='metap_bonferroni', pval_significance_threshold=0.05, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
  # get the truly significant and loosely significant entries
  sig_per_ct <- list()
  # check specificity per cell type
  for(cell_type in intersect(names(output_per_ct_and_condition), cell_types)){
    # list per condition
    sig_per_cond <- list()
    # for each stimulation
    for(stim in intersect(names(output_per_ct_and_condition[[cell_type]]), stims)){
      # grab the entry
      output <- output_per_ct_and_condition[[cell_type]][[stim]]
      # grab the significant ones
      sigs <- rownames(output[output[[pval_column]] < pval_significance_threshold, ])
      # add this to list
      sig_per_cond[[stim]] <- sigs
    }
    # add to list of lists
    sig_per_ct[[cell_type]] <- sig_per_cond
  }
  return(sig_per_ct)
}

get_mast_per_cell_type_and_condition <- function(mast_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
  # create list per cell type
  mast_per_ct <- list()
  # read each cell_type
  for(cell_type in cell_types){
    # create list per condition
    mast_per_condition <- list()
    # read each condition
    for(condition in stims){
      try({
        # paste together the location
        mast_loc <- paste(mast_output_loc, cell_type, 'UT', condition, '.tsv', sep = '')
        # read the table
        mast_output <- read.table(mast_loc, header = T, row.names = 1)
        # put it in the list
        mast_per_condition[[condition]] <- mast_output
      })
    }
    # put in the list
    mast_per_ct[[cell_type]] <- mast_per_condition
  }
  return(mast_per_ct)
}

create_overlap_pathway_df_cell_type <- function(cell_type, locations, use_ranking=T){
  # put all results in a list
  pathways_analysis <- list()
  # put all results in a shared DF
  pathway_df <- NULL
  # check each cell type
  for(location in locations){
    # get the output of the specific cell type
    files_to_read <- list.files(location, pattern = cell_type)
    # check each stim
    for(file in files_to_read){
      try({
        print(file)
        # paste the filepath together
        #filepath <- paste(pathway_output_loc, cell_type, 'UT',stim,'_sig_pathways.txt', sep = '')
        filepath <- paste(location, file, sep = '')
        # read the file
        pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
        # create column name
        newcolname <- substr(file, 1, nchar(file)-4)
        # get the log2 of the significance value
        #pathways[[newcolname]] <- log2(pathways[[sig_val_to_use]])
        if(use_ranking){
          pathways[[newcolname]] <- as.numeric(rownames(pathways))
        }
        else{
          pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
        }
        pathways$id_name <- paste(pathways$ID, pathways$Name, sep = '_')
        # reduce to the only two columns we care about
        pathways <- pathways[, c('id_name', newcolname)]
        # join with other pathway files
        if(is.null(pathway_df)){
          # just set as df if the first round through
          pathway_df <- pathways
          pathway_df <- data.table(pathway_df, key = c('id_name'))
        }
        else{
          # otherwise, merge with existing pathways
          pathway_df <- merge(pathway_df, data.table(pathways, key = c('id_name')), by.x='id_name', by.y='id_name', all=T)
        }
      })
    }
  }
  # turn into regular df
  pathway_df <- setDF(pathway_df)
  # set all NA to zero
  pathway_df[is.na(pathway_df)] <- 0
  # set rownames
  rownames(pathway_df) <- pathway_df$id_name
  pathway_df$id_name <- NULL
  # remove rows which amount to zero
  pathway_df <- pathway_df[apply(pathway_df[,-1], 1, function(x) !all(x==0)),]
  return(pathway_df)
}


get_top_vary_genes <- function(de_table, use_tp=T, use_pathogen=T, use_ct=T, sd_cutoff=0.5, use_dynamic_sd=F, top_so_many=10, must_be_positive_once=F, pathogens=c("CA", "MTB", "PA"), timepoints=c("3h", "24h"), cell_types=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")){
  top_vary_de <- c()
  cols_to_loop <- NULL
  # grab the appriate grep
  if(use_tp & use_pathogen){
    # we want a combo of pathogen and timepoints, so 3hCA for example
    cols_to_loop <- paste(rep(timepoints, each = length(pathogens)), pathogens, sep = "")
  }
  else if(use_pathogen & use_ct){
    # cell type and pathogen, so monocyte3hCA and monocyte24hCA for example
    cols_to_loop <- paste(rep(cell_types, each = length(pathogens)), pathogens, sep = ".*")
  }
  else if(use_tp & use_ct){
    # cell type at a timepoint, so monocyte3hCA and monocyte3hPA and monocyte3hMTB for example
    cols_to_loop <- paste(rep(cell_types, each = length(timepoints)), timepoints, sep = ".*")
  }
  else if(use_pathogen){
    cols_to_loop <- pathogens
  }
  else if(use_tp){
    cols_to_loop <- timepoints
  }
  else if(use_ct){
    cols_to_loop <- cell_types
  }
  # go through our group of columns
  for(col_grep in cols_to_loop){
    # grab the column names that have this in their name
    appropriate_columns <- colnames(de_table)[(grep(col_grep, colnames(de_table)))]
    print('getting most varying out of: ')
    print(appropriate_columns)
    # now subset the frame to only have these columns
    sub_de_table <- de_table[, appropriate_columns]
    # subset to only the genes that were upregulated at least once, if requested
    if(must_be_positive_once){
      sub_de_table <- sub_de_table[apply(sub_de_table,1,min) < 0,]
    }
    # we will return the rownames
    varying_genes <- NULL
    # either use a set SD or grab so many genes
    if(use_dynamic_sd){
      varying_genes <- get_most_varying_from_df(sub_de_table, top_so_many)
    }
    else{
      # now calculate the sd over this set of columns
      sds <- apply(sub_de_table, 1, sd, na.rm=T)
      # then grab the genes that are 'this' varied
      varying_genes <- rownames(sub_de_table[sds > sd_cutoff,])
    }

    # and add them to the list
    top_vary_de <- c(top_vary_de, varying_genes)
  }
  # constrain to the unique genes
  top_vary_de <- unique(top_vary_de)
  top_vary_de <- sort(top_vary_de)
  return(top_vary_de)
}

get_most_varying_from_df <- function(dataframe, top_so_many=10, dont_get_least_varying=T){
  # now calculate the sd over this set of columns
  sds <- apply(dataframe, 1, sd, na.rm=T)
  # add the sds as a column
  dataframe$sds <- sds
  # order by the sd
  dataframe <- dataframe[order(dataframe$sds, decreasing = dont_get_least_varying), ]
  # we will return the rownames
  most_varied <- NULL
  # we need to make sure we can return as many rownames as requested
  if(nrow(dataframe) < top_so_many){
    print(paste('requested ', top_so_many, ', but dataframe only has ', nrow(most_varied), ' rows', sep = ''))
    most_varied <- rownames(dataframe)
  }
  else{
    most_varied <- rownames(dataframe)[1:top_so_many]
  }
  return(most_varied)
}

filter_pathway_df_on_starting_id <- function(pathway_df, filtered_pathway_names, remove_id_from_pathway_name=T){
  # get the ones now in the pathway df
  pathway_names_with_id <- rownames(pathway_df)
  last_dash_pos <- "\\_[^\\_]*$"
  # get the pathway names by skipping from the underscore
  pathway_names <- substr(pathway_names_with_id, regexpr(last_dash_pos, pathway_names_with_id)+1, nchar(pathway_names_with_id))
  print(head(pathway_names))
  # filter the pathway df
  pathway_df_filtered <- pathway_df[pathway_names %in% filtered_pathway_names, ]
  # remove the ID from the pathway name if asked
  if(remove_id_from_pathway_name){
    rownames(pathway_df_filtered) <- substr(rownames(pathway_df_filtered), regexpr(last_dash_pos, rownames(pathway_df_filtered))+1, nchar(rownames(pathway_df_filtered)))
  }
  return(pathway_df_filtered)
}

get_filtered_pathway_names <- function(pathway_table, relation_table, starting_id){
  # get all of the children of the starting ID
  all_children <- get_children(relation_table, starting_id)
  # get the names of the pathways that are children
  pathway_names <- as.character(pathway_table[pathway_table$V1 %in% all_children, ]$V2)
  return(pathway_names)
}

get_children <- function(relation_table, starting_id){
  # get all of the children of the starting ID
  children <- as.character(relation_table[relation_table$V1 == starting_id, 'V2'])
  # these children are all family
  family <- children
  # see if there were any children
  if(length(children) > 0){
    # if there were children, we need to get their children as well
    for(child in children){
      # get the grandchildren and add these to the family
      grand_children <- get_children(relation_table, child)
      family <- c(family, grand_children)
    }
  }
  return(family)
}

plot_de_vs_ct_numbers <- function(ct_number_loc, mast_output_loc, cell_types_to_use=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC"), conditions=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), pval_column='metap_bonferroni'){
  ct_numbers <- read.table(ct_number_loc, header = T, sep = '\t')
  # create table to store results
  numbers_table <- NULL
  # check each cell type
  for(cell_type in cell_types_to_use){
    # check each condition
    for(condition in conditions){
      # paste the output loc together
      file <- paste(cell_type, 'UT', 'X', condition, '.tsv', sep='')
      # read the mast output
      mast <- read.table(paste(mast_output_loc, file, sep = ''), header=T, sep='\t', row.names=1)
      # filter to only include the significant results
      mast <- mast[mast[[pval_column]] <= 0.05, ]
      # get number of DE genes
      de_gene_count <- nrow(mast)
      # get the number of cells in chem2
      nr_of_cells_chem2 <- ct_numbers[ct_numbers$chem == 'V2' & ct_numbers$cell_type == cell_type & ct_numbers$condition1 == 'UT' & ct_numbers$condition2 == paste('X', condition, sep = ''), 'nr_of_cells_condition1'] +
        ct_numbers[ct_numbers$chem == 'V2' & ct_numbers$cell_type == cell_type & ct_numbers$condition1 == 'UT' & ct_numbers$condition2 == paste('X', condition, sep = ''), 'nr_of_cells_condition2']
      # and in chem 3
      nr_of_cells_chem3 <- ct_numbers[ct_numbers$chem == 'V3' & ct_numbers$cell_type == cell_type & ct_numbers$condition1 == 'UT' & ct_numbers$condition2 == paste('X', condition, sep = ''), 'nr_of_cells_condition1'] +
        ct_numbers[ct_numbers$chem == 'V3' & ct_numbers$cell_type == cell_type & ct_numbers$condition1 == 'UT' & ct_numbers$condition2 == paste('X', condition, sep = ''), 'nr_of_cells_condition2']
      # sum these
      nr_of_cells_both <- nr_of_cells_chem2 + nr_of_cells_chem3
      # put into dataframe
      print(paste(condition, cell_type, nr_of_cells_both, de_gene_count))
      numbers_combination <- data.frame(c(condition), c(cell_type), c(nr_of_cells_both), c(de_gene_count), stringsAsFactors = F)
      colnames(numbers_combination) <- c('condition', 'cell_type', 'nr_of_cells', 'nr_of_DE_genes')
      # add to larger df
      if(is.null(numbers_table)){
        numbers_table <- numbers_combination
      }
      else{
        numbers_table <- rbind(numbers_table, numbers_combination)
      }
    }
  }
  # the label is different depending on whether we show the proportion or not
  xlab <- 'nr of cells'
  cc <- get_color_coding_dict()
  colScale <- scale_colour_manual(name = "cell_type",values = unlist(cc[cell_types_to_use]))
  ggplot(numbers_table, aes(x=as.numeric(nr_of_cells), y=as.numeric(nr_of_DE_genes), shape=condition, color=cell_type)) +
    geom_point(size=3) +
    labs(x = xlab, y = 'nr of DE genes', title = 'cells vs nr of DE genes') +
    colScale
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
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  return(color_coding)
}


get_gene_list_from_hm_branch <- function(heatmap, branch_directions, use_col=T){
  branch <- NULL
  # grab the row or column
  if(use_col){
    branch <- heatmap$colDendrogram
  }
  else{
    branch <- heatmap$rowDendrogram
  }
  # go through the branch depths to get to the specific branch
  for(branch_direction in branch_directions){
    # the branch direction is 1 for left and 2 for right
    branch <- branch[[branch_direction]]
  }
  # now get all the children of this branch
  genes <- get_child_genes(branch)
  return(genes)
}

get_child_genes <- function(heatmap_branch){
  genes <- c()
  # if the height is zero, we are at the leaf and we can get the gene
  if(attr(heatmap_branch, 'height') == 0){
    genes <- names(attr(heatmap_branch, 'value'))
  }
  # if we are not at zero, there are more branches or leaves and we need to go further down both directions
  else{
    genes_branch1 <- get_child_genes(heatmap_branch[[1]])
    genes_branch2 <- get_child_genes(heatmap_branch[[2]])
    genes <- c(genes_branch1, genes_branch2)
  }
  return(genes)
}

get_average_gene_expression_per_ct_and_tp <- function(seurat_object, condition.column = 'timepoint', cell.type.column = 'cell_type_lowerres', cell_types_to_use=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC"), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), assay='RNA'){
  exp_df <- NULL
  # calculate for each condition
  for(condition in conditions){
    # subset to just the cells of this condition
    seurat_object_condition <- seurat_object[,seurat_object@meta.data[condition.column] == condition]
    # calculate for each cell_type
    for(cell_type in cell_types_to_use){
      # subset to just the cells of the cell type
      seurat_object_cell_type <- seurat_object_condition[,seurat_object_condition@meta.data[cell.type.column] == cell_type]
      # calculate the relevant matrix from the relevant assay
      exp_df_ct_cond <- NULL
      if(assay == 'RNA'){
        DefaultAssay(seurat_object_cell_type) <- 'RNA'
        averages <- apply(seurat_object_cell_type$RNA@data, 1, mean)
        exp_df_ct_cond <- data.frame(condition=rep(condition, times = length(averages)), cell_type=rep(cell_type, times = length(averages)), gene=rownames(seurat_object_cell_type$RNA@data), average=averages)
      }
      else if(assay == 'SCT'){
        DefaultAssay(seurat_object_cell_type) <- 'SCT'
        averages <- apply(seurat_object_cell_type$SCT@counts, 1, mean)
        exp_df_ct_cond <- data.frame(condition=rep(condition, times = length(averages)), cell_type=rep(cell_type, times = length(averages)), gene=rownames(seurat_object_cell_type$SCT@counts), average=averages)
      }
      # paste to the overall df
      if(is.null(exp_df)){
        exp_df <- exp_df_ct_cond
      }
      else{
        exp_df <- rbind(exp_df, exp_df_ct_cond)
      }
    }
  }
  return(exp_df)
}

avg_exp_to_table <- function(avg_exp_table, cell_type, genes){
  subsetted_avg_exp <- avg_exp_table[avg_exp_table$cell_type %in% c(cell_type) & avg_exp_table$gene %in% genes, ]
  avg_exp_matrix <- matrix(, nrow=length(unique(subsetted_avg_exp$gene)), ncol=length(unique(subsetted_avg_exp$condition)), dimnames=list(unique(subsetted_avg_exp$gene), unique(subsetted_avg_exp$condition)))
  for(gene in unique(subsetted_avg_exp$gene)){
    for(condition in unique(subsetted_avg_exp$condition)){
      avg_exp_matrix[gene, condition] <- subsetted_avg_exp[subsetted_avg_exp$condition ==  condition & subsetted_avg_exp$gene == gene, 'average']
    }
  }
  return(avg_exp_matrix)
}

get_normed_to_max_expression <- function(expression_table, cell_type){
  # subset to table of cell_type
  expression_table_ct <- expression_table[expression_table$cell_type == cell_type, ]
  # check each gene
  for(gene in unique(expression_table_ct$gene)){
    # get the highest value
    max_expression <- max(expression_table_ct[expression_table_ct$gene == gene, ]$average)
    # divide the expression in all conditions by that max value
    expression_table_ct[expression_table_ct$gene == gene, ]$average <- expression_table_ct[expression_table_ct$gene == gene, ]$average / max_expression
  }
  return(expression_table_ct)
}

subset_expression_table_by_de <- function(expression_table, mast_output_loc, cell_types_to_use=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC"), conditions=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), pval_column='metap_bonferroni', sig_pval=0.05, only_positive=F, only_negative=F, lfc_column='metafc'){
  genes_de <- c()
  # get the files
  for(cell_type in cell_types_to_use){
    # try to read each file
    for(condition in conditions){
      # paste together the file
      file <- paste(cell_type, 'UT', condition, '.tsv', sep = '')
      try({
        # read the mast output
        mast <- read.table(paste(mast_output_loc, file, sep = ''), header=T)
        # filter to only include the significant results
        mast <- mast[mast[[pval_column]] <= 0.05, ]
        # filter for only the positive lfc if required
        if(only_positive){
          mast <- mast[mast[[lfc_column]] < 0, ]
        }
        # filter for only the positive lfc if required
        if(only_negative){
          mast <- mast[mast[[lfc_column]] > 0, ]
        }
        # grab the genes from the column names
        genes <- rownames(mast)
        # add to our overarching list
        genes_de <- c(genes_de, genes)
      })
    }
  }
  # remove duplicates
  genes_de <- unique(genes_de)
  # subset to genes that are in the DE genes list
  expression_table_filtered <- expression_table[as.character(expression_table$gene) %in% genes_de, ]
  return(expression_table_filtered)
}

avg_exp_table_to_hm_table <- function(expression_table){
  # initialise table
  hm_table <- NULL
  # go through the conditions
  for(condition in unique(expression_table$condition)){
    # get the expression for that table
    expression_table_cond <- expression_table[expression_table$condition == condition, c('gene', 'average')]
    # set colnames so that we can merge these later
    colnames(expression_table_cond) <- c('gene', condition)
    # convert to data.table for efficient merging
    expression_table_cond <- data.table(expression_table_cond)
    # try to merge if necessary
    if(is.null(hm_table)){
      hm_table <- expression_table_cond
    }
    else{
      hm_table <- merge(hm_table, expression_table_cond, by='gene')
    }
  }
  # convert back to regular dataframe
  hm_table <- data.frame(hm_table)
  # set rownames
  rownames(hm_table) <- hm_table$gene
  # remove the old gene column
  hm_table$gene <- NULL
  return(hm_table)
}


pathways_to_hm_colors <- function(expression_heatmap, pathways_named_lists){
  colors_df <- NULL
  for(pathway_name in names(pathways_named_lists)){
    # get the pathway genes
    pathway.genes <- read.table(pathways_named_lists[[pathway_name]], header=F)
    pathway.genes <- as.character(pathway.genes$V1)
    pathway.genes <- pathway.genes[pathway.genes %in% rownames(expression_heatmap)]
    pathway.annotation <- rep("gray97", nrow(expression_heatmap))
    pathway.annotation[rownames(expression_heatmap) %in% pathway.genes] <- "gray55"
    # add to colors df
    if(is.null(colors_df)){
      colors_df <- data.frame(pathway.annotation)
      colnames(colors_df) <- pathway_name
    }
    else{
      colors_df[[pathway_name]] <- pathway.annotation
    }
  }
  # transform to matrix
  colors_m <- as.matrix(colors_df)
  return(colors_m)
}

transform.to.zscore <- function(x){
  if (as.numeric(x[3]) < 5e-324*2){
    x[3] <- 5e-324*2
  }
  z.score <- qnorm(as.numeric(x[3])/2)
  if (as.numeric(x[2]) > 0){
    return(z.score * -1)
  }
  else {
    return(z.score)
  }
}

transform.to.zscore.log10 <- function(x){
	z.score <- -log10(as.numeric(x[3]))
	if (is.infinite(z.score)){
		z.score <- -log10(5e-324)
	}
	if (as.numeric(x[2]) > 0){
		return(z.score)
	}
	else {
		return(z.score * -1)
	}
}

# get the locations of the DE output
mast_output_prepend <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/v'
mast_output_append <- '_paired_lores_lfc01minpct01_20201106/rna/'
# write the location of the combined output
mast_meta_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/'

# write meta output
write_meta_mast(NULL, mast_output_prepend, mast_output_append, mast_meta_output_loc)

# do the same for the 3h vs 24h stuff
mast_meta_output_3h24h_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106_3h24h/rna/'
write_meta_mast_3hvs24h(NULL, mast_output_prepend, mast_output_append, mast_meta_output_3h24h_loc)


# we can go from gene symbols to ensemble IDs with this file
gene_to_ens_mapping <- "/data/scRNA/differential_expression/genesymbol_to_ensid.tsv"
# set the location to write the significant genes
sig_output_loc <- '/data/scRNA/differential_expression/sigs/meta_paired_lores_lfc01minpct01_20201106_ensid/rna/'
# write the significant genes
get_significant_genes(mast_meta_output_loc, sig_output_loc, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_up_output_loc <- '/data/scRNA/differential_expression/sigs_pos/meta_paired_lores_lfc01minpct01_20201106_ensid/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc, sig_up_output_loc, only_positive = T, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_down_output_loc <- '/data/scRNA/differential_expression/sigs_neg/meta_paired_lores_lfc01minpct01_20201106_ensid/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc, sig_down_output_loc, only_negative = T, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)

# set the location to write the significant genes
sig_output_loc <- '/data/scRNA/differential_expression/sigs/meta_paired_lores_lfc01minpct01_20201106/rna/'
# write the significant genes
get_significant_genes(mast_meta_output_loc, sig_output_loc, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_up_output_loc <- '/data/scRNA/differential_expression/sigs_pos/meta_paired_lores_lfc01minpct01_20201106/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc, sig_up_output_loc, only_positive = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_down_output_loc <- '/data/scRNA/differential_expression/sigs_neg/meta_paired_lores_lfc01minpct01_20201106/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc, sig_down_output_loc, only_negative = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)

# do all this stuff for the 3h vs 24h as well
sig_output_3h24h_loc <- '/data/scRNA/differential_expression/sigs/meta_paired_lores_lfc01minpct01_20201106_3h24h/rna/'
get_significant_genes(mast_meta_output_3h24h_loc, sig_output_3h24h_loc, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
sig_up_output_3h24h_loc <- '/data/scRNA/differential_expression/sigs_pos/meta_paired_lores_lfc01minpct01_20201106_3h24h/rna/'
get_significant_genes(mast_meta_output_3h24h_loc, sig_up_output_3h24h_loc, only_positive = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
sig_down_output_3h24h_loc <- '/data/scRNA/differential_expression/sigs_neg/meta_paired_lores_lfc01minpct01_20201106_3h24h/rna/'
get_significant_genes(mast_meta_output_3h24h_loc, sig_down_output_3h24h_loc, only_negative = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
