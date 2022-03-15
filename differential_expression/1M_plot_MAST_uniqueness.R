




####################
# libraries        #
####################

library(ggplot2)
library(ggbeeswarm)

####################
# Functions        #
####################


plot_de_gene_uniqueness_condition <- function(base_mast_output_path, marked_hours=T, condition_combinations=c('UTX3hCA', 'UTX24hCA', 'UTX3hMTB', 'UTX24hMTB', 'UTX3hPA', 'UTX24hPA'), cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), hours=list('3h' = c('UTX3hCA', 'UTX3hMTB', 'UTX3hPA'), '24h' = c('UTX24hCA', 'UTX24hMTB', 'UTX24hPA')),  pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc',use_color_coding_dict=F, to_ens=F, symbols.to.ensg.mapping='genes.tsv', return_table=F){
  # we want to get all the DE and where they were significant
  gene_conditions_df <- NULL
  # so, check each condition
  for(condition in condition_combinations){
    # we'll store the genes
    genes_condition <- c()
    # then check each cell type
    for(cell_type in cell_types_to_use){
      # paste together the filename
      full_mast_path <- paste(base_mast_output_path, cell_type, condition, '.tsv', sep = '')
      # get the genes
      genes_cell_type <- get_de_genes_from_mast_file(full_mast_path, pval_column=pval_column, sig_pval=sig_pval, max=max, max_by_pval=max_by_pval, only_positive=only_positive, only_negative=only_negative, lfc_column=lfc_column, to_ens=to_ens, symbols.to.ensg.mapping=symbols.to.ensg.mapping)
      # add to the list of genes for this condition
      genes_condition <- c(genes_condition, genes_cell_type)
    }
    # make the genes unique, as the cell types have overlap
    genes_condition <- unique(genes_condition)
    # set as a dataframe
    condition_df <- data.frame(genes_condition, rep(T, times=length(genes_condition)), stringsAsFactors = F)
    # set the condition as the column name
    colnames(condition_df) <- c('gene', condition)
    # merge with existing gene and conditions
    if(is.null(gene_conditions_df)){
      gene_conditions_df <- condition_df
    }
    else{
      gene_conditions_df <- merge(gene_conditions_df, condition_df, by = 'gene', all = T)
    }
  }
  # genes that were not found before merging, will be NA, those are not found and thus of course F
  gene_conditions_df[is.na(gene_conditions_df)] <- F
  # disregard the gene column (makes the applies easier)
  gene_conditions_df <- gene_conditions_df[, condition_combinations]
  # get for each genes how many conditions it was significant in
  number_of_times_sig <- apply(gene_conditions_df, 1, sum)
  # turn into plot df
  plot_df <- NULL
  # we need a bar for unique number that the DE genes were unique
  for(number_sig in unique(number_of_times_sig)){
    # get the number of times this was the case
    number_times_this_sig <- sum(number_of_times_sig == number_sig)
    # grab for this the data in the gene conditions
    gene_conditions_df_this_number <- gene_conditions_df[number_of_times_sig == number_sig, ]
    # we will calculate how many are specific
    specific <- 0
    # check how many of these are specific to a timepoint
    for(timepoint_group in names(hours)){
      # get the group column names
      columns_included <- hours[[timepoint_group]]
      # get all the other groups
      columns_groups_other <- setdiff(names(hours), timepoint_group)
      # and grab all those columns
      columns_other <- as.vector(unlist(hours[columns_groups_other]))
      # check which are in this timepoint group, but not in of the others
      specificity_yes <- apply(gene_conditions_df_this_number, 1, function(x){
        specific_timepoint <- (sum(x[columns_included]) > 0 & sum(x[columns_other]) == 0)
        return(specific_timepoint)
      })
      # check how many
      specificity_yes_number  <- sum(specificity_yes)
      # increase the number of specific ones
      specific <- specific + specificity_yes_number
      # add to plot frame
      row <- data.frame(number_sig=c(number_sig), number_times_this_sig=c(specificity_yes_number), condition=c(timepoint_group), stringsAsFactors = F)
      # add or set
      if(is.null(plot_df)){
        plot_df <- row
      }
      else{
        plot_df <- rbind(plot_df, row)
      }
    }
    # calculate what is left
    non_specific <- number_times_this_sig - specific
    # create appropriate row
    row <- data.frame(number_sig=number_sig, number_times_this_sig=non_specific, condition=c('mixed'), stringsAsFactors = F)
    # add or set
    if(is.null(plot_df)){
      plot_df <- row
    }
    else{
      plot_df <- rbind(plot_df, row)
    }
  }
  # set the condition of the numbers
  #plot_df$condition <- 'mixed'
  # we have to supply the tables used to make the figures, which we'll do in this manner
  if(return_table){
    return(plot_df)
  }
  else{
    # set the order I like for the legend, but setting the factor order
    plot_df$condition <- factor(plot_df$condition, levels=c('mixed', setdiff(unique(plot_df$condition), 'mixed')))
    
    # make the plot finally
    p <- ggplot(plot_df, aes(fill=condition, y=number_times_this_sig, x=number_sig)) +
      geom_bar(position='stack', stat='identity') +
      labs(x='number of conditions a gene is differentially expressed in', y='Number of significant DE genes') +
      ggtitle('overlap of DE genes in condition combinations') +
      labs(fill = "Found in")
    # use standard colors if requested
    if(use_color_coding_dict){
      # grab the colours
      cc <- get_color_coding_dict()
      # add the 'mixed' condition
      cc[['mixed']] <- 'gray'
      fillScale <- scale_fill_manual(name = "condition",values = unlist(cc[c(names(hours), 'mixed')]))
      p <- p + fillScale
    }
    
    return(p)
  }
}

plot_de_gene_uniqueness_celltype <- function(base_mast_output_path, marked_singles=T, condition_combinations=c('UTX3hCA', 'UTX24hCA', 'UTX3hMTB', 'UTX24hMTB', 'UTX3hPA', 'UTX24hPA'), cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv', paper_style=F, return_table=F){
  # we want to get all the DE and where they were significant
  gene_cell_type_df <- NULL
  # so, check each condition
  for(cell_type in cell_types_to_use){
    # we'll store the genes
    genes_cell_type <- c()
    # then check each cell type
    for(condition in condition_combinations){
      # paste together the filename
      full_mast_path <- paste(base_mast_output_path, cell_type, condition, '.tsv', sep = '')
      # get the genes
      genes_condition <- get_de_genes_from_mast_file(full_mast_path, pval_column=pval_column, sig_pval=sig_pval, max=max, max_by_pval=max_by_pval, only_positive=only_positive, only_negative=only_negative, lfc_column=lfc_column, to_ens=to_ens, symbols.to.ensg.mapping=symbols.to.ensg.mapping)
      # add to the list of genes for this condition
      genes_cell_type <- c(genes_cell_type, genes_condition)
    }
    # make the genes unique, as the cell types have overlap
    genes_cell_type <- unique(genes_cell_type)
    # set as a dataframe
    cell_type_df <- data.frame(genes_cell_type, rep(T, times=length(genes_cell_type)), stringsAsFactors = F)
    # set the condition as the column name
    colnames(cell_type_df) <- c('gene', cell_type)
    # merge with existing gene and conditions
    if(is.null(gene_cell_type_df)){
      gene_cell_type_df <- cell_type_df
    }
    else{
      gene_cell_type_df <- merge(gene_cell_type_df, cell_type_df, by = 'gene', all = T)
    }
  }
  # genes that were not found before merging, will be NA, those are not found and thus of course F
  gene_cell_type_df[is.na(gene_cell_type_df)] <- F
  # disregard the gene column (makes the applies easier)
  gene_cell_type_df <- gene_cell_type_df[, cell_types_to_use]
  # get for each genes how many conditions it was significant in
  number_of_times_sig <- apply(gene_cell_type_df, 1, sum)
  # turn into plot df
  plot_df <- NULL
  # we need a bar for unique number that the DE genes were unique
  for(number_sig in unique(number_of_times_sig)){
    # get the number of times this was the case
    number_times_this_sig <- sum(number_of_times_sig == number_sig)
    # create appropriate row
    row <- data.frame(number_sig=number_sig, number_times_this_sig=number_times_this_sig, stringsAsFactors = F)
    # add or set
    if(is.null(plot_df)){
      plot_df <- row
    }
    else{
      plot_df <- rbind(plot_df, row)
    }
  }
  # set the condition of the numbers
  plot_df$cell_type <- 'mixed'
  # for the singles, so in one condition only, we might want to see the proportions
  if(marked_singles){
    # remove the singles, as we're overwriting those
    plot_df <- plot_df[plot_df$number_sig != 1, ]
    # subset to the singles
    gene_cell_type_df_singles <- gene_cell_type_df[number_of_times_sig == 1, ]
    # now that we have only the singles, we can sum over the columns, to get the number of genes unique to the condition
    number_unique_per_cell_type <- apply(gene_cell_type_df_singles, 2, sum)
    # make that into a df
    plot_df_uniques <- data.frame(number_sig=rep(1, times=length(cell_types_to_use)), number_times_this_sig=number_unique_per_cell_type, cell_type=cell_types_to_use)
    # add to the current plot
    plot_df <- rbind(plot_df, plot_df_uniques)
  }
  # set the order I like for the legend, but setting the factor order
  plot_df$cell_type <- factor(plot_df$cell_type, levels=c('mixed', cell_types_to_use))
  # we are supposed to supply the tables to make the figures, we can do that like this
  if(return_table){
    return(plot_df)
  }
  else{
    # grab the colours
    cc <- get_color_coding_dict()
    # add the 'mixed' condition
    cc[['mixed']] <- 'gray'
    fillScale <- scale_fill_manual(name = "cell type",values = unlist(cc[c(cell_types_to_use, 'mixed')]))
    # make the plot finally
    p <- ggplot(plot_df, aes(fill=cell_type, y=number_times_this_sig, x=number_sig)) +
      geom_bar(position='stack', stat='identity') +
      labs(x='number of cell types a gene is differentially expressed in', y='Number of significant DE genes') +
      ggtitle('overlap of DE genes in cell types') +
      labs(fill = "Found in") +
      fillScale
    if(paper_style){
      p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
    }
    return(p)
  }
}


get_de_genes_from_mast_file <- function(mast_full_file_path, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv'){
  significant_genes <- c()
  try({
    # read the mast output
    mast <- read.table(mast_full_file_path, header=T, row.names = 1)
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
    # confine in some way if reporting a max number of genes
    if(!is.null(max)){
      # by p if required
      if(max_by_pval){
        mast <- mast[order(mast[[pval_column]]), ]
      }
      # by lfc otherwise
      else{
        mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
        # if we only have the upregulated genes, order the other way around
        if(only_positive){
          mast <- mast[order(mast[[lfc_column]], decreasing = F), ]
        }
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
    significant_genes <- genes
  })
  return(significant_genes)
}


plot_de_gene_numbers <- function(base_mast_output_path, condition_combinations=c('UTX3hCA', 'UTX24hCA', 'UTX3hMTB', 'UTX24hMTB', 'UTX3hPA', 'UTX24hPA'), cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), condition_to_hour=list('UTX3hCA'='3h', 'UTX3hMTB' = '3h', 'UTX3hPA' = '3h', 'UTX24hCA' = '24h', 'UTX24hMTB' = '24h', 'UTX24hPA' = '24h'), pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', paper_style=F, return_table=F){
  # we want to get all the DE and where they were significant
  gene_number_table <- NULL
  # so, check each condition
  for(cell_type in cell_types_to_use){
    # then check each cell type
    for(condition in condition_combinations){
      # paste together the filename
      full_mast_path <- paste(base_mast_output_path, cell_type, condition, '.tsv', sep = '')
      # get the genes
      genes_condition <- get_de_genes_from_mast_file(full_mast_path, pval_column=pval_column, sig_pval=sig_pval, max=max, max_by_pval=max_by_pval, only_positive=only_positive, only_negative=only_negative, lfc_column=lfc_column, to_ens=F, symbols.to.ensg.mapping=NULL)
      # count the number of genes
      nr_genes <- length(genes_condition)
      # get the time
      time <- condition_to_hour[[condition]]
      # add to a row
      row <- data.frame(cell_type = c(cell_type), condition = c(condition), time=c(time), nr_genes=c(nr_genes))
      # add to table
      if(is.null(gene_number_table)){
        gene_number_table <- row
      }
      else{
        gene_number_table <- rbind(gene_number_table, row)
      }
    }
  }
  # we are supposed to supply the tables of the figures, we'll do that like this
  if(return_table){
    return(gene_number_table)
  }
  else{
    p  <- ggplot(gene_number_table) + 
      geom_quasirandom(aes(x=time, y=nr_genes), size=2.6) +
      geom_quasirandom(aes(x=time, y=nr_genes, color=condition), size=2.4) +
      scale_y_continuous(minor_breaks = seq(500 , 3000, 500), breaks = seq(500, 3000, 500))
    if(paper_style){
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
    }
    #scale_color_manual(values=condition.colors) +
    p <- p + ylim(0, max(gene_number_table$nr_genes)) +
    guides(size=FALSE, fill=FALSE, shape=FALSE) +
    labs(x="Time point", y="Number of significant DE genes", guide="Stimulation") +
    ggtitle("Number of DE genes per cell type and condition")
    
    p <- p + facet_grid(cols=vars(cell_type))
    return(p)
  }
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  # set the condition colors
  color_coding <- list()
  color_coding[['UT']] <- 'lightgrey'
  color_coding[["3hCA"]] <- "darkolivegreen2"
  color_coding[["24hCA"]] <- "forestgreen"
  color_coding[["3hMTB"]] <- "lightskyblue"
  color_coding[["24hMTB"]] <- "deepskyblue3"
  color_coding[["3hPA"]] <- "sandybrown"
  color_coding[["24hPA"]] <- "darkorange1"
  color_coding[['3h']] <- 'green'
  color_coding[['24h']] <- 'blue'
  color_coding[['3h specific']] <- 'green'
  color_coding[['24h specific']] <- 'blue'
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}

label_dict <- function(){
  label_dict <- list()
  label_dict[["UT"]] <- "UT"
  label_dict[["X3hCA"]] <- "3hCA"
  label_dict[["X24hCA"]] <- "24hCA"
  label_dict[["X3hMTB"]] <- "3hMTB"
  label_dict[["X24hMTB"]] <- "24hMTB"
  label_dict[["X3hPA"]] <- "3hPA"
  label_dict[["X24hPA"]] <- "24hPA"
  label_dict[["3hCA"]] <- "3hCA"
  label_dict[["24hCA"]] <- "24hCA"
  label_dict[["3hMTB"]] <- "3hMTB"
  label_dict[["24hMTB"]] <- "24hMTB"
  label_dict[["3hPA"]] <- "3hPA"
  label_dict[["24hPA"]] <- "24hPA"
  label_dict[['3h']] <- '3h specific'
  label_dict[['24h']] <- '24h specific'
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["hemapoietic stem"]] <- "hemapoietic stem"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["plasma B"]] <- "plasma B"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["megakaryocyte"]] <- "megakaryocyte"
  label_dict[["T_other"]] <- "other T"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  return(label_dict)
}


####################
# Main Code        #
####################


# DE output locations
mast_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/'
# make plot of celltype uniqueness
plot_de_gene_uniqueness_celltype(mast_output_loc)
# make table instead
de_overlap_celltype_table <- plot_de_gene_uniqueness_celltype(mast_output_loc, return_table = T)
# make plot of the condition uniqueness
plot_de_gene_uniqueness_condition(mast_output_loc)
# make table instead
de_overlap_condition_table <- plot_de_gene_uniqueness_condition(mast_output_loc, return_table = T)
# plot the number of DE genes
plot_de_gene_numbers(mast_output_loc)
# make table instead
de_number_table <- plot_de_gene_numbers(mast_output_loc, return_table = T)

