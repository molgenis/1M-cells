####################
# libraries        #
####################

library(Seurat)
library(ggplot2)


####################
# Functions        #
####################

# plot the module score in a boxplot
plot_average_expression <- function(seurat_object, module_score_column_name, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_type_column='cell_type_lowerres', condition_column='timepoint', title='pathway', color_by_ct=T, to_per_part=F, participant_column='participant'){
  # get the metadata
  metadata <- seurat_object@meta.data
  # limit the cell types
  metadata <- metadata[metadata[[cell_type_column]] %in% cell_types, ]
  # limit the conditions
  metadata <- metadata[metadata[[condition_column]] %in% conditions, ]
  # set some hardcoded column names
  metadata$ct_column <- metadata[[cell_type_column]]
  metadata$tp_column <- metadata[[condition_column]]
  metadata$msc_column <- metadata[[module_score_column_name]]
  print(head(metadata))
  metadata$cttp_column <- paste(metadata$ct_column, metadata$tp_column, sep = '\n')
  # convert to per_participant if requested
  if(to_per_part){
    metadata <- metadata_to_average_module_per_part(metadata, cttp_column='cttp_column', msc_column='msc_column', ct_column='ct_column', condition_column='tp_column', participant_column=participant_column)
  }
  # grab our dict for colors
  cc <- get_color_coding_dict()
  # if the colouring is by cell type
  if(color_by_ct){
    colScale <- scale_fill_manual(name = metadata$ct_column, values = unlist(cc[cell_types]))
    ggplot(metadata, aes(x=cttp_column, y=msc_column, fill=ct_column)) +
      geom_boxplot() +
      colScale +
      ggtitle(title) +
      labs(y = 'module scores', x='cell type and condition')
  }
  # if the colouring is by the condition
  else{
    ct_tp_order <- c()
    for(cell_type in cell_types){
      ct_tp_order <- c(ct_tp_order, paste(rep(cell_type, each = length(conditions)), conditions, sep = "\n"))
    }
    colScale <- scale_fill_manual(name = metadata$tp_column, values = unlist(cc[conditions]))
    ggplot(metadata, aes(x=cttp_column, y=msc_column, fill=tp_column)) +
      geom_boxplot() +
      colScale +
      ggtitle(title) +
      labs(y = 'module scores', x='cell type and condition') +
      scale_x_discrete(limits = ct_tp_order)
  }
}

metadata_to_average_module_per_part <- function(metadata, cttp_column='cttp_column', msc_column='msc_column', ct_column='ct_column', condition_column='tp_column', participant_column='participant'){
  metadata_averaged <- NULL
  for(cttp in unique(metadata[[cttp_column]])){
    # subset to cell type and condition combo
    metadata_cttp <- metadata[metadata[[cttp_column]] == cttp ,]
    # check for each participant here
    for(participant in unique(metadata_cttp[[participant_column]])){
      # calculate the mean
      mean_participant <- mean(metadata_cttp[metadata_cttp[[participant_column]], msc_column])
      # get other data
      ct <- metadata_cttp[metadata_cttp[[participant_column]], ct_column][1]
      tp <- metadata_cttp[metadata_cttp[[participant_column]], condition_column][1]
      # set as a dataframe
      metadata_averaged_row <- data.frame(cttp_column=c(cttp), msc_column=c(mean_participant), ct_column=c(ct), condition_column=c(tp), participant_column=c(participant))
      # set as initial or append
      if(is.null(metadata_averaged)){
        metadata_averaged <- metadata_averaged_row
      }
      else{
        metadata_averaged <- rbind(metadata_averaged, metadata_averaged_row)
      }
    }
  }
  # set colnames as before
  colnames(metadata_averaged) <- c(cttp_column, msc_column, ct_column, condition_column, participant_column)
  return(metadata_averaged)
}

# add the module score from a table of genes (now using as pathway plot)
add_module_score_from_table <- function(pathway_gene_table_loc, pathway_name, seurat_object, and_plot=T, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'Baseline', 't24h', 't8w'), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = T, to_per_part=F, participant_column='participant.final'){
  # get the cytokine genes
  pathway_df <- read.table(pathway_gene_table_loc)
  pathway_genes <- pathway_df$V1
  pathway_genes_in_seurat_object <- intersect(rownames(seurat_object), pathway_genes)
  # add to list
  pathway_list_seurat_object <- list()
  pathway_list_seurat_object[[pathway_name]] <- pathway_genes_in_seurat_object
  # add the module score
  seurat_object <- AddModuleScore(seurat_object, features = pathway_list_seurat_object, name = pathway_name)
  if(and_plot){
    plot_average_expression(seurat_object, module_score_column_name=paste(pathway_name, '1', sep = ''), cell_types=cell_types, conditions=conditions, cell_type_column=cell_type_column, condition_column=condition_column, title=pathway_name, color_by_ct = color_by_ct, to_per_part=to_per_part, participant_column=participant_column)
  }
  return(seurat_object)
}

# mapping of conditions and cell types to colours
get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- 'grey'
  color_coding[["3hCA"]] <- "khaki2"
  color_coding[["24hCA"]] <- "khaki4"
  color_coding[["3hMTB"]] <- "paleturquoise1"
  color_coding[["24hMTB"]] <- "paleturquoise3"
  color_coding[["3hPA"]] <- "rosybrown1"
  color_coding[["24hPA"]] <- "rosybrown3"
  color_coding[["X3hCA"]] <- "khaki2"
  color_coding[["X24hCA"]] <- "khaki4"
  color_coding[["X3hMTB"]] <- "paleturquoise1"
  color_coding[["X24hMTB"]] <- "paleturquoise3"
  color_coding[["X3hPA"]] <- "rosybrown1"
  color_coding[["X24hPA"]] <- "rosybrown3"
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

####################
# Main code        #
####################

# locations of objects
v2_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds'
v3_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds'

# location of the pathway genes
pathway_gene_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/pathways/'
cytokine_signalling_loc <- paste(pathway_gene_loc, 'REACTOME_Cytokine_Signaling_in_Immune_system_genes.txt', sep = '')
interferon_signalling_loc <- paste(pathway_gene_loc, 'REACTOME_Interferon_Signaling_genes.txt', sep = '')
antigen_processing_cross_presentation_loc <- paste(pathway_gene_loc, 'REACTOME_Antigen_processing-Cross_presentation.txt', sep = '')
class_I_MHC_mediated_antigen_processing_and_presentation_loc <- paste(pathway_gene_loc, 'REACTOME_Class_I_MHC_mediated_antigen_processing_and_presentation.txt', sep = '')
CLEC7A_Dectin_1_signaling_loc <- paste(pathway_gene_loc, 'REACTOME_CLEC7A_(Dectin-1)_signaling.txt', sep = '')
DAP12_signaling_loc <- paste(pathway_gene_loc, 'REACTOME_DAP12_signaling.txt', sep = '')
interferon_alpha_or_beta_signaling <- paste(pathway_gene_loc, 'REACTOME_Interferon_alpha_or_beta_signaling.txt', sep = '')
interleukin_2_signaling_loc <- paste(pathway_gene_loc, 'REACTOME_Interleukin-2_signaling.txt', sep = '')
interleukin_10_signaling_loc <- paste(pathway_gene_loc, 'REACTOME_Interleukin-10_signaling.txt', sep = '')
ye_interferon_response_loc <- paste(pathway_gene_loc, 'YE_INF_shared.txt', sep = '')
ye_interferon1_response_loc <- paste(pathway_gene_loc, 'YE_INF_type1.txt', sep = '')
ye_interferon2_response_loc <- paste(pathway_gene_loc, 'YE_INF_type2.txt', sep = '')

# read object
v2 <- readRDS(v2_loc)

# get the cytokine genes
cytokine_df <- read.table(cytokine_signalling_loc)
cytokine_genes <- cytokine_df$V1
cytokine_genes_in_v2 <- intersect(rownames(v2), cytokine_genes)
# add to list
cytokines_list_v2 <- list()
cytokines_list_v2[['cytokine_signalling']] <- cytokine_genes_in_v2
# add the module score
v2 <- AddModuleScore(v2, features = cytokines_list_v2, name = 'cytokine_signalling')
plot_average_expression(v2, module_score_column_name='cytokine_signalling1', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'X3hCA'), cell_type_column='cell_type_lowerres', condition_column='timepoint', title='cytokine_signalling')

# get the cytokine genes
interferon_df <- read.table(interferon_signalling_loc)
interferon_genes <- interferon_df$V1
interferon_genes_in_v2 <- intersect(rownames(v2), interferon_genes)
# add to list
interferon_list_v2 <- list()
interferon_list_v2[['interferon_signalling']] <- interferon_genes_in_v2
# add the module score
v2 <- AddModuleScore(v2, features = interferon_list_v2, name = 'interferon_signalling')
plot_average_expression(v2, module_score_column_name='interferon_signalling1', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'X3hCA'), cell_type_column='cell_type_lowerres', condition_column='timepoint')

# get the cytokine genes
antigen_processing_cross_presentation_df <- read.table(antigen_processing_cross_presentation_loc)
antigen_processing_cross_presentation_genes <- antigen_processing_cross_presentation_df$V1
antigen_processing_cross_presentation_genes_in_v2 <- intersect(rownames(v2), antigen_processing_cross_presentation_genes)
# add to list
antigen_processing_cross_presentation_list_v2 <- list()
antigen_processing_cross_presentation_list_v2[['antigen_processing_cross_presentation']] <- antigen_processing_cross_presentation_genes_in_v2
# add the module score (manually, doing automatically later)
v2 <- AddModuleScore(v2, features = antigen_processing_cross_presentation_list_v2, name = 'antigen_processing_cross_presentation')
plot_average_expression(v2, module_score_column_name='antigen_processing_cross_presentation1', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'X24hCA'), cell_type_column='cell_type_lowerres', condition_column='timepoint', title='antigen_processing_cross_presentation')
plot_average_expression(v2, module_score_column_name='antigen_processing_cross_presentation1', cell_types=c('monocyte'), cell_type_column='cell_type_lowerres', condition_column='timepoint', title='antigen_processing_cross_presentation', color_by_ct = F)

# doing it in one go, basically all I did before
add_module_score_from_table(CLEC7A_Dectin_1_signaling_loc, 'CLEC7A_Dectin_1_signaling', v2, cell_types=c('monocyte'), cell_type_column='cell_type_lowerres', condition_column='timepoint', color_by_ct = F)
add_module_score_from_table(DAP12_signaling_loc, 'DAP12_signaling', v2, and_plot=T, cell_types=c('monocyte'), cell_type_column='cell_type_lowerres', condition_column='timepoint', color_by_ct = F)


