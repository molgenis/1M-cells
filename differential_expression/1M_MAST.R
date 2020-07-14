####################
# libraries        #
####################

library(MAST)
library(Seurat)

####################
# Functions        #
####################

perform_mast <- function(seurat_object, output_loc, condition.1, condition.2, split.column = 'timepoint', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25, features = NULL, latent.vars=NULL){
  # set the assay we wish to analyse
  DefaultAssay(seurat_object) <- assay
  # a hard column name is required for the model.matrix step, so let's add that
  seurat_object <- AddMetaData(seurat_object, seurat_object@meta.data[split.column], 'condition')
  Idents(seurat_object) <- seurat_object@meta.data$condition
  # create an output directory
  output_loc_final <- paste(output_loc, condition.1, condition.2, '.tsv', sep = '')
  # perform MAST if possible
  result <- NULL
  # we need to 'try' here, as MAST ungraciously exits if it cannot perform
  try({
    result <- FindMarkers(object = seurat_object, ident.1 = condition.1, ident.2 = condition.2, test.use = 'MAST', min.pct = min.pct, logfc.threshold = logfc.threshold, assay = assay, features = features, latent.vars = latent.vars)
  })
  # write the actual table if possible
  if(is.null(result)){
    print(paste('nothing to compare due to threshold for', output_loc_final))
  } else{
    write.table(result, output_loc_final, sep = '\t', row.names = T)
  }
}

perform_mast_per_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL){
  # do subselection based on features
  features = NULL
  # grab the top expressed genes if that is what the user wanted
  if(!is.null(use_top_expressed)){
    features <- get_top_expressed_features(seurat_object, use_top_expressed)
  }
  # go through the cell types
  for(cell_type in unique(seurat_object@meta.data[[cell.type.column]])){
    # make a more specific path
    output_loc_cell_type <- paste(output_loc, cell_type, sep = '')
    # grab the subset
    seurat_object_cell_type <- seurat_object[,seurat_object@meta.data[cell.type.column] == cell_type]
    # do the MAST
    stims_to_do <- c('CA', 'PA', 'MTB')
    # if user want to do other stims, that's fine
    if(!is.null(stims)){
      stims_to_do <- stims
    }
    # for the three stims (looping here, so don't have to subset the celltype multiple times)
    for(stim in stims_to_do){
      # paste together the conditions
      tp3h <- paste('X3h', stim, sep = '')
      tp24h <- paste('X24h', stim, sep = '')
      try({perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'UT', condition.2 = tp3h, split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, features = features, latent.vars = latent.vars)})
      try({perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'UT', condition.2 = tp24h ,split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, features = features, latent.vars = latent.vars)})
      try({perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = tp3h, condition.2 = tp24h ,split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, features = features, latent.vars = latent.vars)})
    }
  }
  # finally do a bulk analysis as well
  output_loc_bulk <- paste(output_loc, 'bulk', sep = '')
  for(stim in c('CA', 'PA', 'MTB')){
    # paste together the conditions
    tp3h <- paste('X3h', stim, sep = '')
    tp24h <- paste('X24h', stim, sep = '')
    try({perform_mast(seurat_object, output_loc_bulk, condition.1 = 'UT', condition.2 = tp3h, split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, latent.vars = latent.vars)})
    try({perform_mast(seurat_object, output_loc_bulk, condition.1 = 'UT', condition.2 = tp24h, split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, latent.vars = latent.vars)})
    try({perform_mast(seurat_object, output_loc_bulk, condition.1 = tp3h, condition.2 = tp24h, split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, latent.vars = latent.vars)})
  }
}

get_top_expressed_features <- function(seurat_object, top_genes_number){
  # go to RNA assay
  DefaultAssay(seurat_object) <- 'RNA'
  # check most variable features again
  seurat_object <- FindVariableFeatures(seurat_object)
  # get the feature means
  seurat_object_means <- seurat_object[['RNA']]@meta.features
  # sort by mean expression
  seurat_object_means <- seurat_object_means[with(seurat_object_means, order(-vst.mean)),]
  # get the top 1000 expressed
  seurat_object_top_expressed <- rownames(seurat_object_means)[1:top_genes_number]
  # return those genes
  return(seurat_object_top_expressed)
}

####################
# Main Code        #
####################

# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
#object_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
#object_loc <- '/data/p287578/1M_cells_scRNAseq/scanpy_preprocess_samples/objects/'
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')
#object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617_anon.rds', sep = '')
#object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617_anon.rds', sep = '')

# DE output locations
#mast_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/pairwise_DE_comparison/'
#mast_output_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/pairwise_DE_comparison/'
#mast_output_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/pairwise_DE_comparison_20200617_top1000/'
#mast_output_loc <- '/data/p287578/1M_cells_scRNAseq/differential_expression/seurat_MAST/output/pairwise_DE_comparison_20200617/'
#mast_output_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/pairwise_DE_comparison_20200617_top1000/'
mast_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/'

# for a MAST comparison, also do only paired comparisons
mast_output_paired_loc_v2 <- paste(mast_output_loc, 'v2_paired/', sep = '')
mast_output_paired_loc_v3 <- paste(mast_output_loc, 'v3_paired/', sep = '')
mast_output_paired_lores_loc_v2 <- paste(mast_output_loc, 'v2_paired_lores_lfc025minpct01ncountrna_20200713/', sep = '')
mast_output_paired_lores_loc_v3 <- paste(mast_output_loc, 'v3_paired_lores_lfc025minpct01ncountrna_20200713/', sep = '')
mast_output_paired_loc_v2_rna <- paste(mast_output_paired_loc_v2, 'rna/', sep = '')
mast_output_paired_loc_v3_rna <- paste(mast_output_paired_loc_v3, 'rna/', sep = '')
mast_output_paired_loc_v2_sct <- paste(mast_output_paired_loc_v2, 'sct/', sep = '')
mast_output_paired_loc_v3_sct <- paste(mast_output_paired_loc_v3, 'sct/', sep = '')
mast_output_paired_lores_loc_v2_rna <- paste(mast_output_paired_lores_loc_v2, 'rna/', sep = '')
mast_output_paired_lores_loc_v3_rna <- paste(mast_output_paired_lores_loc_v3, 'rna/', sep = '')
mast_output_paired_lores_loc_v2_sct <- paste(mast_output_paired_lores_loc_v2, 'sct/', sep = '')
mast_output_paired_lores_loc_v3_sct <- paste(mast_output_paired_lores_loc_v3, 'sct/', sep = '')


# put in the work for v2
v2 <- readRDS(object_loc_v2)
#perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_lores_loc_v2_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
#perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_lores_loc_v2_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA', min.pct = 0, logfc.threshold = 0, stims = c('CA'))
#perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_lores_loc_v2_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25)
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_lores_loc_v2_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25, latent.vars=c('nCount_RNA'))
#perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_loc_v2_sct, cell.type.column = 'cell_type', assay = 'SCT')
#perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_loc_v2_rna, cell.type.column = 'cell_type', assay = 'RNA')
rm(v2)

# put in the work for v3
v3 <- readRDS(object_loc_v3)
#perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_lores_loc_v3_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
#perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_lores_loc_v3_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA', min.pct = 0, logfc.threshold = 0, stims = c('CA'))
#perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_lores_loc_v3_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25)
perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_lores_loc_v3_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25, latent.vars=c('nCount_RNA'))
#perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_loc_v3_sct, cell.type.column = 'cell_type', assay = 'SCT')
#perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_loc_v3_rna, cell.type.column = 'cell_type', assay = 'RNA')
rm(v3)
