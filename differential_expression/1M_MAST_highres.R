############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_MAST_highres.R
# Function: perform MAST DE analysis on higher resolution cell types
############################################################################################################################


####################
# libraries        #
####################

library(MAST)
library(Seurat)
library(Matrix)

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

perform_mast_per_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, cell_types_to_use=NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL){
  # do subselection based on features
  features = NULL
  # grab the top expressed genes if that is what the user wanted
  if(!is.null(use_top_expressed)){
    features <- get_top_expressed_features(seurat_object, use_top_expressed)
  }
  cell_types <- unique(as.character(seurat_object@meta.data[[cell.type.column]]))
  # confine if requested
  if(!is.null(cell_types_to_use)){
    cell_types <- intersect(cell_types, cell_types_to_use)
  }
  # go through the cell types
  for(cell_type in cell_types){
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
}


####################
# Main Code        #
####################

# object locations
object_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029_azimuth.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106_azimuth.rds', sep = '')
object_loc_v2_new <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20210905.rds', sep = '')
object_loc_v3_new <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20210905.rds', sep = '')
# DE output locations
mast_output_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/'

# for a MAST comparison, also do only paired comparisons
mast_output_paired_highres_loc_v2 <- paste(mast_output_loc, 'v2_paired_highres_lfc01minpct01_20210905/', sep = '')
mast_output_paired_highres_loc_v3 <- paste(mast_output_loc, 'v3_paired_highres_lfc01minpct01_20210905/', sep = '')
# we'll use the RNA assay
mast_output_paired_highres_rna_loc_v2 <- paste(mast_output_paired_highres_loc_v2, 'rna/', sep = '')
mast_output_paired_highres_rna_loc_v3 <- paste(mast_output_paired_highres_loc_v3, 'rna/', sep = '')


# read the object
v2 <- readRDS(object_loc_v2)
DefaultAssay(v2) <- 'RNA'
# we've done some refinements at the marker gene level, let's make those changes permanent
v2@meta.data[v2@meta.data$cell_type == 'NK', 'cell_type'] <- 'NKdim'
# clean up
v2 <- v2[, !is.na(v2@meta.data$cell_type) & !is.na(v2@meta.data$assignment) & !is.na(v2@meta.data$timepoint)]
# write the new object
saveRDS(v2, object_loc_v2_new)
# do the mapping
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_highres_rna_loc_v2, cell_types_to_use = c('NKdim', 'NKbright', 'mDC', 'pDC', 'mono 1', 'mono 2', 'mono 3', 'mono 4'), logfc.threshold = 0.1)


# read the object
v3 <- readRDS(object_loc_v3)
DefaultAssay(v3) <- 'RNA'
# we've done some refinements at the marker gene level, let's make those changes permanent
v3@meta.data[v2@meta.data$cell_type == 'NK', 'cell_type'] <- 'NKdim'
# clean up
v3 <- v3[, !is.na(v3@meta.data$cell_type) & !is.na(v3@meta.data$assignment) & !is.na(v3@meta.data$timepoint)]
# write the new object
saveRDS(v3, object_loc_v3_new)
# do the mapping
perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_highres_rna_loc_v3, cell_types_to_use = c('NKdim', 'NKbright', 'mDC', 'pDC', 'mono 1', 'mono 2', 'mono 3', 'mono 4'), logfc.threshold = 0.1)


