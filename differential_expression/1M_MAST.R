####################
# libraries        #
####################

library(MAST)
library(Seurat)

####################
# Functions        #
####################

perform_mast <- function(seurat_object, output_loc, condition.1, condition.2, split.column = 'timepoint', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25){
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
    result <- FindMarkers(object = seurat_object, ident.1 = condition.1, ident.2 = condition.2, test.use = 'MAST', min.pct = min.pct, logfc.threshold = logfc.threshold, assay = assay)
  })
  # write the actual table if possible
  if(is.null(result)){
    print(paste('nothing to compare due to threshold for', output_loc_final))
  } else{
    write.table(result, output_loc_final, sep = '\t', row.names = T)
  }
}

perform_mast_per_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', min.pct = 0.1){
  # go through the cell types
  for(cell_type in unique(seurat_object@meta.data[[cell.type.column]])){
    # make a more specific path
    output_loc_cell_type <- paste(output_loc, cell_type, sep = '')
    # grab the subset
    seurat_object_cell_type <- seurat_object[,seurat_object@meta.data[cell.type.column] == cell_type]
    # do variable feature detection, why Hilde?
    # seurat_object <- FindVariableFeatures(object = seurat_object, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
    # do the MAST
    # for the three stims (looping here, so don't have to subset the celltype multiple times)
    for(stim in c('CA', 'PA', 'MTB')){
      # paste together the conditions
      tp3h <- paste('X3h', stim, sep = '')
      tp24h <- paste('X24h', stim, sep = '')
      try({perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'UT', condition.2 = tp3h, split.column = split.column, assay = assay, min.pct = min.pct)})
      try({perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'UT', condition.2 = tp24h ,split.column = split.column, assay = assay, min.pct = min.pct)})
      try({perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = tp3h, condition.2 = tp24h ,split.column = split.column, assay = assay, min.pct = min.pct)})
    }
  }
  # finally do a bulk analysis as well
  output_loc_bulk <- paste(output_loc, 'bulk', sep = '')
  for(stim in c('CA', 'PA', 'MTB')){
    # paste together the conditions
    tp3h <- paste('X3h', stim, sep = '')
    tp24h <- paste('X24h', stim, sep = '')
    try({perform_mast(seurat_object, output_loc_bulk, condition.1 = 'UT', condition.2 = tp3h, split.column = split.column, assay = assay, min.pct = min.pct)})
    try({perform_mast(seurat_object, output_loc_bulk, condition.1 = 'UT', condition.2 = tp24h, split.column = split.column, assay = assay, min.pct = min.pct)})
    try({perform_mast(seurat_object, output_loc_bulk, condition.1 = tp3h, condition.2 = tp24h, split.column = split.column, assay = assay, min.pct = min.pct)})
  }
}

####################
# Main Code        #
####################

# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20200520.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200520.rds', sep = '')

# DE output locations
#mast_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/pairwise_DE_comparison/'
mast_output_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/pairwise_DE_comparison/'
# for a MAST comparison, also do only paired comparisons
mast_output_paired_loc_v2 <- paste(mast_output_loc, 'v2_paired/', sep = '')
mast_output_paired_loc_v3 <- paste(mast_output_loc, 'v3_paired/', sep = '')
mast_output_paired_lores_loc_v2 <- paste(mast_output_loc, 'v2_paired_lores/', sep = '')
mast_output_paired_lores_loc_v3 <- paste(mast_output_loc, 'v3_paired_lores/', sep = '')
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
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_lores_loc_v2_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_lores_loc_v2_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA')
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_loc_v2_sct, cell.type.column = 'cell_type', assay = 'SCT')
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_loc_v2_rna, cell.type.column = 'cell_type', assay = 'RNA')
rm(v2)

# put in the work for v3
v3 <- readRDS(object_loc_v3)
perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_lores_loc_v3_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_lores_loc_v3_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA')
perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_loc_v3_sct, cell.type.column = 'cell_type', assay = 'SCT')
perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_loc_v3_rna, cell.type.column = 'cell_type', assay = 'RNA')
rm(v3)
