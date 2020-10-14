library(Seurat)
library(ggplot2)


plot_average_expression <- function(seurat_object, module_score_column_name, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_type_column='cell_type_lowerres', condition_column='timepoint', title='pathway', color_by_ct=T){
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
  
  cc <- get_color_coding_dict()
  if(color_by_ct){
    colScale <- scale_fill_manual(name = metadata$ct_column, values = unlist(cc[cell_types]))
    ggplot(metadata, aes(x=cttp_column, y=msc_column, fill=ct_column)) +
      geom_boxplot() +
      colScale +
      ggtitle(title) +
      labs(y = 'module scores', x='cell type and condition')
  }
  else{
    ct_tp_order <- paste(rep(cell_type, each = length(conditions)), conditions, sep = "\n")
    colScale <- scale_fill_manual(name = metadata$tp_column, values = unlist(cc[conditions]))
    ggplot(metadata, aes(x=cttp_column, y=msc_column, fill=tp_column)) +
      geom_boxplot() +
      colScale +
      ggtitle(title) +
      labs(y = 'module scores', x='cell type and condition') +
      scale_x_discrete(limits = ct_tp_order)
  }
}

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




# locations of objects
v2_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds'
v3_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds'

# location of the pathway genes
pathway_gene_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/pathways/'
cytokine_signalling_loc <- paste(pathway_gene_loc, 'REACTOME_Cytokine_Signaling_in_Immune_system_genes.txt', sep = '')
interferon_signalling_loc <- paste(pathway_gene_loc, 'REACTOME_Interferon_Signaling_genes.txt', sep = '')

# read object
v2 <- readRDS(v2_loc)
v2 <- v2[,!is.na(v2@meta.data$timepoint)]
v2 <- v2[,!is.na(v2@meta.data$assignment)]

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

