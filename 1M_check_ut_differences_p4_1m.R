

library(Seurat)
library(MAST)


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

perform_mast_per_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL, condition.1 = '1M_cells', condition.2 = 'pilot4_unstimulated'){
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
    try({perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = condition.1, condition.2 = condition.2 ,split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, features = features, latent.vars = latent.vars)})
  }
  # finally do a bulk analysis as well
  output_loc_bulk <- paste(output_loc, 'bulk', sep = '')
  try({perform_mast(seurat_object, output_loc_bulk, condition.1 = condition.1, condition.2 = condition.2, split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, latent.vars = latent.vars)})
}

perform_mast_per_celltype_subsampled <- function(seurat_object, seurat_object2, output_loc, subsample_size, subsample_times=10, subsample_column='assignment', split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL){
  for(i in 1:subsample_times){
    # set the location
    output_loc_ss <- paste(output_loc, i, '_', sep = '')
    # grab a subsample
    subsampled_assignments <- sample(unique(seurat_object@meta.data[[subsample_column]]), subsample_size)
    # grab the subsampled object
    seurat_object_subsampled <- seurat_object[, seurat_object@meta.data[[subsample_column]] %in% subsampled_assignments]
    # merge the two objects
    seurat_object_merged <- merge(seurat_object_subsampled, seurat_object2)
    # normalize data
    seurat_object_merged <- NormalizeData(seurat_object_merged)
    # actually perform
    perform_mast_per_celltype(seurat_object_merged, output_loc_ss, split.column = split.column, cell.type.column = cell.type.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, use_top_expressed = use_top_expressed, latent.vars=latent.vars)
  }
}

perform_mast_per_celltype_subsampled_same <- function(seurat_object, output_loc, subsample_size, subsample_times=10, subsample_column='assignment', split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL){
  for(i in 1:subsample_times){
    # set the location
    output_loc_ss <- paste(output_loc, i, '_', sep = '')
    # grab a subsample
    subsampled_assignments <- as.character(sample(unique(seurat_object@meta.data[[subsample_column]]), subsample_size))
    subsampled_assignments_split <- split(subsampled_assignments, sample(2, length(subsampled_assignments), repl = TRUE) )
    subsampled_assignments_1 <- subsampled_assignments_split[[1]]
    subsampled_assignments_2 <- subsampled_assignments_split[[2]]
    # grab the subsampled object
    seurat_object_subsampled <- seurat_object[, seurat_object@meta.data[[subsample_column]] %in% subsampled_assignments]
    seurat_object_subsampled@meta.data[seurat_object_subsampled@meta.data[[subsample_column]] %in% subsampled_assignments_1, ][[split.column]] <- '1'
    seurat_object_subsampled@meta.data[seurat_object_subsampled@meta.data[[subsample_column]] %in% subsampled_assignments_2, ][[split.column]] <- '2'
    # normalize data
    seurat_object_subsampled <- NormalizeData(seurat_object_subsampled)
    # actually perform
    perform_mast_per_celltype(seurat_object_subsampled, output_loc_ss, split.column = split.column, cell.type.column = cell.type.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, use_top_expressed = use_top_expressed, latent.vars=latent.vars, condition.1 = '1', condition.2 = '2')
  }
}



# get the old object
pilot4_old_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/ut_compare/objects/pilot4.Rds'
pilot4_old <- readRDS(pilot4_old_loc)
# add cell_types
pilot4_old@meta.data$cell_type <- pilot4_old@ident
# convert to gene symbols
gene_to_ens_mapping <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/resources/features_v3.tsv"
genes <- read.table(gene_to_ens_mapping, header = F, stringsAsFactors = F)
genes$V2 <- gsub("_", "-", make.unique(genes$V2))
gene_symbols <- genes[match(rownames(pilot4_old@data), genes$V1),"V2"]
rownames(pilot4_old@data) <- gene_symbols
rownames(pilot4_old@raw.data) <- gene_symbols
rownames(pilot4_old@scale.data) <- gene_symbols
# recreate seurat object
pilot4 <- CreateSeuratObject(counts = pilot4_old@raw.data)
pilot4@meta.data <- pilot4_old@meta.data
# update to new version
#pilot4 <- UpdateSeuratObject(pilot4_old)
# grab just UT
pilot4_ut <- pilot4[, !is.na(pilot4@meta.data$stimulation) & pilot4@meta.data$stimulation == 'unstimulated']
# clear some memory
rm(pilot4)
# get the new v2 object
v2_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds'
v2 <- readRDS(v2_loc)
# grab just UT
v2_ut <- v2[, !is.na(v2@meta.data$timepoint) & v2@meta.data$timepoint == 'UT']
# clear some memory
rm(v2)
# harmonise the cell types
v2_ut@meta.data$cell_type <- v2_ut@meta.data$cell_type_lowerres
#pilot4_ut@meta.data$cell_type <- pilot4_ut@active.ident
levels(pilot4_ut@meta.data$cell_type) <- c(levels(pilot4_ut@meta.data$cell_type), 'CD4T', 'CD8T', 'DC', 'monocyte')
pilot4_ut@meta.data[pilot4_ut@meta.data$cell_type == 'CD4+_T', ]$cell_type <- 'CD4T'
pilot4_ut@meta.data[pilot4_ut@meta.data$cell_type == 'CD8+_T', ]$cell_type <- 'CD8T'
pilot4_ut@meta.data[pilot4_ut@meta.data$cell_type == 'pDC', ]$cell_type <- 'DC'
pilot4_ut@meta.data[pilot4_ut@meta.data$cell_type == 'Monocyte', ]$cell_type <- 'monocyte'
pilot4_ut@meta.data$cell_type <- droplevels(pilot4_ut@meta.data$cell_type)

# merge the objects
ut_merged <- merge(pilot4_ut, v2_ut)
#normalize
ut_merged <- NormalizeData(ut_merged)

# set the location of the output
mast_output_loc_full <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/ut_compare/mast_output/full_'
mast_output_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/ut_compare/mast_output/'
# do all vs all
perform_mast_per_celltype(seurat_object = ut_merged, output_loc = mast_output_loc_full, split.column = 'orig.ident', cell.type.column = 'cell_type', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25)
#
perform_mast_per_celltype_subsampled(v2_ut, pilot4_ut, output_loc = mast_output_loc, subsample_size=10, subsample_times=20, subsample_column='assignment', split.column = 'orig.ident', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL)
perform_mast_per_celltype_subsampled_same(v2_ut, output_loc = mast_output_loc, subsample_size=20, subsample_times=20, subsample_column='assignment', split.column = 'orig.ident', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL)

