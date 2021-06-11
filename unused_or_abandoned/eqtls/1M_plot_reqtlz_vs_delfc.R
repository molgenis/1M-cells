
plot_reqtl_z_vs_de_lfc <- function(reqtl_output_location, mast_output_location, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), only_significant_de=T, only_significant_reqtl=T, sig_column_de='metap_bonferroni', sig_column_reqtl='FDR', sig_threshold_de=0.05, sig_threshold_reqtl=0.05, lfc_column='metafc', gene_to_ens_mapping=NULL){
  # read all the differentially expressed data
  stims_with_prepend <- paste('X', stims, sep = '')
  all_de <- get_mast_per_cell_type_and_condition(mast_output_location, cell_types=cell_types, stims=stims_with_prepend, only_significant = only_significant_de, sig_column = sig_column_de, sig_threshold = sig_threshold_de)
  # read all the re-QTL output
  all_reqtl <- get_reqtl_per_cell_type_and_conditions(reqtl_output_location, cell_types=cell_types, stims=stims, only_significant = only_significant_reqtl, sig_column = sig_column_reqtl, sig_threshold = sig_threshold_reqtl)
  # make plot list
  plots_per_ct <- list()
  # make a plot per cell type, which we can only do for cell types with both DE and reQTL results
  for(cell_type in intersect(names(all_reqtl), names(all_de))){
    # make plot list
    plots_per_condition <- list()
    # check each condition
    for(stim in stims){
      # only if they are in both
      if(stim %in% names(all_reqtl[[cell_type]]) & paste('X', stim, sep = '') %in% names(all_de[[cell_type]])){
        # grab the DE table
        de_output <- all_de[[cell_type]][[paste('X', stim, sep = '')]]
        # grab the reQTL table
        reqtl_output <- all_reqtl[[cell_type]][[stim]]
        # make DE results all positive
        de_output[de_output[[lfc_column]] < 0, lfc_column] <- de_output[de_output[[lfc_column]] < 0, lfc_column] * -1
        # make Z-scores all positive
        reqtl_output[reqtl_output[['OverallZScore']] < 0, 'OverallZScore'] <- reqtl_output[reqtl_output[['OverallZScore']] < 0, 'OverallZScore'] * -1
        # plot
        de_vs_reqtl_plot_df <- NULL
        # 
        if(!is.null(gene_to_ens_mapping)){
          # change the rownames of the DE to the ensemble IDs
          gene_to_ens_mapping$gene <- gsub("_", "-", (make.unique(gene_to_ens_mapping$gene)))
          print(head(gene_to_ens_mapping[match(rownames(de_output), gene_to_ens_mapping$gene), 'ens']))
          rownames(de_output) <- gene_to_ens_mapping[match(rownames(de_output), gene_to_ens_mapping$gene), 'ens']
          # get the common genes
          common_genes <- intersect(as.character(reqtl_output$ProbeName), rownames(de_output))
          # make a dataframe of the overlap
          de_vs_reqtl_plot_df <- data.frame(gene=common_genes, lfc=de_output[common_genes, lfc_column], Z=reqtl_output[match(common_genes, reqtl_output$ProbeName), 'OverallZScore'])
        }
        else{
          # get the common genes
          common_genes <- intersect(as.character(reqtl_output$HGNCName), rownames(de_output))
          # make a dataframe of the overlap
          de_vs_reqtl_plot_df <- data.frame(gene=common_genes, lfc=de_output[common_genes, lfc_column], Z=reqtl_output[match(common_genes, reqtl_output$HGNCName), 'OverallZScore'])
        }
        # sort by the DE genes
        de_vs_reqtl_plot_df <- de_vs_reqtl_plot_df[order(de_vs_reqtl_plot_df$lfc), ]
        # make the actual plot
        p <- ggplot(data=de_vs_reqtl_plot_df, aes(x=lfc, y=Z)) + geom_point() + xlim(c(0, 3)) + ylim(c(0, 10)) + ggtitle(paste('DE vs LFC of', cell_type, 'in UT vs', stim)) + geom_smooth(method='lm')
        # add to plot list
        plots_per_condition[[stim]] <- p
      }
      # add to plot list
      plots_per_ct[[cell_type]] <- plots_per_condition
    }
  }
  return(plots_per_ct)
}

plot_reqtl_z_vs_de_lfc_combined <- function(reqtl_output_location, mast_output_location, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), only_significant_de=T, only_significant_reqtl=T, sig_column_de='metap_bonferroni', sig_column_reqtl='FDR', sig_threshold_de=0.05, sig_threshold_reqtl=0.05, lfc_column='metafc', gene_to_ens_mapping=NULL){
  # read all the differentially expressed data
  stims_with_prepend <- paste('X', stims, sep = '')
  all_de <- get_mast_per_cell_type_and_condition(mast_output_location, cell_types=cell_types, stims=stims_with_prepend, only_significant = only_significant_de, sig_column = sig_column_de, sig_threshold = sig_threshold_de)
  # read all the re-QTL output
  all_reqtl <- get_reqtl_per_cell_type_and_conditions(reqtl_output_location, cell_types=cell_types, stims=stims, only_significant = only_significant_reqtl, sig_column = sig_column_reqtl, sig_threshold = sig_threshold_reqtl)
  # make plot list
  plots_per_ct <- list()
  # make a plot per cell type, which we can only do for cell types with both DE and reQTL results
  for(cell_type in intersect(names(all_reqtl), names(all_de))){
    # create the plot DF
    de_vs_reqtl_plot_df <- NULL
    # check each condition
    for(stim in stims){
      # only if they are in both
      if(stim %in% names(all_reqtl[[cell_type]]) & paste('X', stim, sep = '') %in% names(all_de[[cell_type]])){
        # grab the DE table
        de_output <- all_de[[cell_type]][[paste('X', stim, sep = '')]]
        # grab the reQTL table
        reqtl_output <- all_reqtl[[cell_type]][[stim]]
        # make DE results all positive
        de_output[de_output[[lfc_column]] < 0, lfc_column] <- de_output[de_output[[lfc_column]] < 0, lfc_column] * -1
        # make Z-scores all positive
        reqtl_output[reqtl_output[['OverallZScore']] < 0, 'OverallZScore'] <- reqtl_output[reqtl_output[['OverallZScore']] < 0, 'OverallZScore'] * -1
        #
        de_vs_reqtl_plot_df_stim <- NULL
        # 
        if(!is.null(gene_to_ens_mapping)){
          # change the rownames of the DE to the ensemble IDs
          gene_to_ens_mapping$gene <- gsub("_", "-", (make.unique(gene_to_ens_mapping$gene)))
          rownames(de_output) <- gene_to_ens_mapping[match(rownames(de_output), gene_to_ens_mapping$gene), 'ens']
          # get the common genes
          common_genes <- intersect(as.character(reqtl_output$ProbeName), rownames(de_output))
          # make a dataframe of the overlap
          de_vs_reqtl_plot_df_stim <- data.frame(gene=common_genes, lfc=de_output[common_genes, lfc_column], Z=reqtl_output[match(common_genes, reqtl_output$ProbeName), 'OverallZScore'])
        }
        else{
          # get the common genes
          common_genes <- intersect(as.character(reqtl_output$HGNCName), rownames(de_output))
          # make a dataframe of the overlap
          de_vs_reqtl_plot_df_stim <- data.frame(gene=common_genes, lfc=de_output[common_genes, lfc_column], Z=reqtl_output[match(common_genes, reqtl_output$HGNCName), 'OverallZScore'])
        }
        if(is.null(de_vs_reqtl_plot_df)){
          de_vs_reqtl_plot_df <- de_vs_reqtl_plot_df_stim
        }
        else{
          de_vs_reqtl_plot_df <- rbind(de_vs_reqtl_plot_df, de_vs_reqtl_plot_df_stim)
        }
      }
    }
    # sort by lfc
    de_vs_reqtl_plot_df <- de_vs_reqtl_plot_df[order(de_vs_reqtl_plot_df$lfc), ]
    # make the actual plot
    p <- ggplot(data=de_vs_reqtl_plot_df, aes(x=lfc, y=Z)) + geom_point() + xlim(c(0, 3)) + ylim(c(0, 10)) + ggtitle(paste('DE vs LFC of', cell_type)) + geom_smooth(method='lm')
    # add to plot list
    plots_per_ct[[cell_type]] <- p
  }
  return(plots_per_ct)
}


get_mast_per_cell_type_and_condition <- function(mast_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB'), only_significant=T, sig_column='metap_bonferroni', sig_threshold=0.05){
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
        # filter for significant if requested
        if(only_significant){
          mast_output <- mast_output[mast_output[[sig_column]] < sig_threshold, ]
        }
        # put it in the list
        mast_per_condition[[condition]] <- mast_output
      })
    }
    # put in the list
    mast_per_ct[[cell_type]] <- mast_per_condition
  }
  return(mast_per_ct)
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

# the mapping used for ensemble ID to gene symbol in QTL analysis
gene_to_ens_mapping_loc <- "/data/scRNA/gene_to_ensemble.tsv"
gene_to_ens_mapping <- read.table(gene_to_ens_mapping_loc, header = T, sep = '\t', stringsAsFactors = F)
# location of the reQTL mapping
reqtl_loc <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/'
# location of the DE analysis
de_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/'
# get the plots
lfc_z_plots <- plot_reqtl_z_vs_de_lfc(reqtl_loc, de_loc, gene_to_ens_mapping = gene_to_ens_mapping)

