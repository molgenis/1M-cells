####################
# Libraries        #
####################

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)


####################
# Functions        #
####################

# create the feature plots
plot_celltype_markers <- function(seurat_object, grouping_column='seurat_clusters', celltype_marker_genes=c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34"), assay = "RNA", slot="data", plot_dir = "./"){
  # set the identity
  Idents(seurat_object) <- grouping_column
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # set a standard value
  seurat_object$clusters <- seurat_object[[grouping_column]]
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$clusters
      LabelClusters(plot = p, id = 'clusters')
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}


# create violin plots
plot_celltype_violins <- function(seurat_object, grouping_column='seurat_clusters', celltype_marker_genes=c("CCR7","S100A4","CD3E","CD4","CD8A", "CD8B", "FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34"), assay = "RNA", slot="data", plot_dir = "./"){
  # set the identity
  Idents(seurat_object) <- grouping_column
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), assay=assay, slot=slot)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}


plot_dimplots <- function(seurat_object, plot_prepend, grouping_columns=c('seurat_clusters'), reduction='umap'){
  # check each of the grouping columns that we have
  for(column in grouping_columns){
    # set up where to save
    output_loc <- paste(plot_prepend, column, '.pdf', sep = '')
    # make the plot
    p <- DimPlot(seurat_object, group.by = column, reduction = reduction)
    # save this plot
    ggsave(output_loc, plot = p)
  }
}

subset_and_recluster <- function(seurat_object, subset_column, subset_values, use_sct=T, dims=1:30, resolution=1.2){
  # subset to the specific variables in the variable column
  seurat_subset <- seurat_object[, seurat_object@meta.data[[subset_column]] %in% subset_values]
  # clear any normalization that was already done
  try({seurat_subset@assays$SCT <- NULL})
  try({seurat_subset@assays$RNA@data <- NULL})
  try({seurat_subset@assays$RNA@scale.data <- NULL})
  try({seurat_subset@graphs <- NULL})
  try({seurat_subset@neighbors <- NULL})
  try({seurat_subset@reductions <- NULL})
  # we will set the default data to operate on
  assay <- 'SCT'
  # normalize the data
  if(use_sct){
    seurat_subset <- SCTransform(seurat_object)
  }
  else{
    seurat_subset <- NormalizeData(seurat_object)
    assay <- 'RNA'
  }
  # first calculate the PCs
  seurat_subset <- RunPCA(seurat_subset, assay = assay)
  # do dimensional reduction to 2d space, based on the PCs
  seurat_subset <- RunUMAP(seurat_subset, dims = dims)
  # use k-nearest neighbour to determine how different cells are
  seurat_subset <- FindNeighbors(seurat_subset)
  # see how they form discrete clusters
  seurat_subset <- FindClusters(seurat_subset, resolution = resolution)
  # return our new object
  return(seurat_subset)
}


split_integrate_and_cluster <- function(seurat_object, split_column, dims=1:30, resolution=1.2){
  # split by the column provided
  seurat_objects_split <- SplitObject(seurat_object, split.by = split_column)
  # normalize each set separately
  seurat_objects_split <- lapply(X = seurat_objects_split, FUN = SCTransform)
  # select the genes to use for integration
  features <- SelectIntegrationFeatures(object.list = seurat_objects_split, nfeatures = 3000)
  # prepare before the integration
  seurat_objects_split <- PrepSCTIntegration(object.list = seurat_objects_split, anchor.features = features)
  # get the anchors
  anchors <- FindIntegrationAnchors(object.list = seurat_objects_split, normalization.method = 'SCT', anchor.features = features)
  # perform the actual integration
  seurat_objects_integrated <- IntegrateData(anchorset = anchors, normalization.method = 'SCT')
  # clear up some memory
  rm(seurat_objects_split)
  rm(anchors)
  # do the normal clustering workflow
  seurat_objects_integrated <- RunPCA(seurat_objects_integrated)
  seurat_objects_integrated <- RunUMAP(seurat_objects_integrated, dims = dims)
  seurat_objects_integrated <- FindNeighbors(seurat_objects_integrated)
  seurat_objects_integrated <- FindClusters(seurat_objects_integrated, resolution = resolution)
  # because we need it later on, let's also do the RNA Normalization
  DefaultAssay(seurat_objects_integrated) <- 'RNA'
  seurat_objects_integrated <- NormalizeData(seurat_objects_integrated)
  return(seurat_objects_integrated)
}


split_integrate_and_cluster_old <- function(seurat_object, split_column, sets=list('CA' = c('UT', 'X3hCA', 'X24hCA'), 'MTB' = c('UT', 'X3hMTB', 'X24hMTB'), 'PA' = c('UT', 'X3hPA', 'X24hPA')), dims=1:30, resolution=1.2){
  # clear any normalization that was already done
  try({seurat_object@assays$SCT <- NULL})
  try({seurat_object@assays$RNA@data <- NULL})
  try({seurat_object@assays$RNA@scale.data <- NULL})
  try({seurat_object@graphs <- NULL})
  try({seurat_object@neighbors <- NULL})
  try({seurat_object@reductions <- NULL})
  # init list to save results
  integrated_objects <- list()
  # save an integrated object per set
  for(set_name in names(sets)){
    # grab the set to use
    set <- sets[[set_name]]
    # subset to that set
    seurat_object_set <- seurat_object[, seurat_object@meta.data[[split_column]] %in% set]
    # split the object by the column to separate objects
    seurat_object_set.list <- SplitObject(seurat_object_set, split.by = split_column)
    # normalize and identify variable features for each dataset independently
    seurat_object_set.list <- lapply(X = seurat_object_set.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = seurat_object_set.list)
    # get the anchors
    anchors <- FindIntegrationAnchors(object.list = seurat_object_set.list, anchor.features = features)
    # perform the actual integration
    seurat_objects_integrated <- IntegrateData(anchorset = anchors)
    # clear up memory
    rm(seurat_object_set.list)
    # use integrated assay
    DefaultAssay(seurat_objects_integrated) <- "integrated"
    # Run the standard workflow for visualization and clustering
    seurat_objects_integrated <- ScaleData(seurat_objects_integrated, verbose = FALSE)
    seurat_objects_integrated <- RunPCA(seurat_objects_integrated, npcs = length(dims), verbose = FALSE)
    seurat_objects_integrated <- RunUMAP(seurat_objects_integrated, reduction = "pca", dims = dims)
    seurat_objects_integrated <- FindNeighbors(seurat_objects_integrated, reduction = "pca", dims = dims)
    seurat_objects_integrated <- FindClusters(seurat_objects_integrated, resolution = resolution)
    # switch to RNA slot
    DefaultAssay(seurat_objects_integrated) <- 'RNA'
    seurat_objects_integrated <- NormalizeData(seurat_objects_integrated)
    # add this to the list we were keeping track of
    integrated_objects[[set_name]] <- seurat_objects_integrated
  }
  return(integrated_objects)
}


# add metadata that is based on existing metadata in the seurat object
add_imputed_meta_data <- function(seurat_object, column_to_transform, column_to_reference, column_to_create){
  # add the column
  seurat_object@meta.data[column_to_create] <- NA
  # go through the grouping we have for the entire object
  for(group in unique(seurat_object@meta.data[[column_to_transform]])){
    # subset to get only this group
    seurat_group <- seurat_object[,seurat_object@meta.data[[column_to_transform]] == group]
    best_group <- 'unknown'
    best_number <- 0
    # check against the reference column
    for(reference in unique(seurat_group@meta.data[[column_to_reference]])){
      # we don't care for the NA reference, if we had all data, we wouldn't need to do this anyway
      if(is.na(reference) == F){
        # grab the number of cells in this group, with this reference
        number_of_reference_in_group <- nrow(seurat_group@meta.data[seurat_group@meta.data[[column_to_reference]] == reference & is.na(seurat_group@meta.data[[column_to_reference]]) == F,])
        correctpercent <- number_of_reference_in_group/ncol(seurat_group)
        print(paste(group,"matches",reference,correctpercent,sep=" "))
        # update numbers if better match
        if(number_of_reference_in_group > best_number){
          best_number <- number_of_reference_in_group
          best_group <- reference
        }
      }
    }
    print(paste("setting ident:",best_group,"for group", group, sep=" "))
    # set this best identity
    seurat_object@meta.data[seurat_object@meta.data[[column_to_transform]] == group,][column_to_create] <- best_group
    # force cleanup
    rm(seurat_group)
  }
  return(seurat_object)
}



####################
# Main Code        #
####################

# object locations
object_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds', sep = '')
object_loc_v3_t <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106_T.rds', sep = '')
# plot loc
plot_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/sub_classification/plots/'
plot_loc_violin <- paste(plot_loc, 'violin_plots/', sep = '')
plot_loc_feature <- paste(plot_loc, 'feature_plots/', sep = '')
plot_loc_dim <- paste(plot_loc, 'dimplots/', sep = '')
# classifications
classifications_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/seurat_multimodal/classifications/'
classifications_loc_v2 <- paste(classifications_loc, 'cluster_azi_t_classifications_v2.tsv', sep = '')
classifications_loc_v3 <- paste(classifications_loc, 'cluster_azi_t_classifications_v3.tsv', sep = '')

# read the object
v3 <- readRDS(object_loc_v3)

# subset and recluster
v3_t <- subset_and_recluster(v3, 'cell_type_lowerres', c('CD4T', 'CD8T'))
# save the object
saveRDS(v3_t, object_loc_v3_t)

v3_t  <- NormalizeData(v3_t)
plot_dimplots(seurat_object = v3_t, paste(plot_loc_dim, 'v3_t_20201106_rna_data_', sep = ''), grouping_columns = c('seurat_clusters', 'timepoint'))
plot_dimplots(seurat_object = v3_t, paste(plot_loc_dim, 'v3_t_20201106_rna_data_old_', sep = ''), grouping_columns = c('cell_type', 'cell_type_lowerres'))
plot_celltype_violins(v3_t, plot_dir = paste(plot_loc_violin, 'v3_t_20201106_rna_data_', sep = ''))
plot_celltype_markers(v3_t, plot_dir = paste(plot_loc_feature, 'v3_t_20201106_rna_data_', sep = ''))



# let's try again, but this time with integration
object_loc <- '/data/p287578/1M_cells_scRNAseq/seurat_preprocess_samples/objects/'
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds', sep = '')
object_loc_v3_t_int <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20210914_T_int.rds', sep = '')
# read the object
v3 <- readRDS(object_loc_v3)
# subset to what we are interested in
v3_t <- v3[, as.character(v3@meta.data$cell_type_lowerres) %in% c('CD4T', 'CD8T')]
# clear up memory
rm(v3)
# integrated and recluster
v3_t <- split_integrate_and_cluster(v3_t, 'timepoint')
# save the result
saveRDS(v3_t, object_loc_v3_t_int)
# also plot here
plot_dimplots(seurat_object = v3_t, paste(plot_loc_dim, 'v3_t_20210914_int_rna_data_', sep = ''), grouping_columns = c('seurat_clusters', 'timepoint'))
plot_dimplots(seurat_object = v3_t, paste(plot_loc_dim, 'v3_t_20210914_int_rna_data_old_', sep = ''), grouping_columns = c('cell_type', 'cell_type_lowerres'))
plot_celltype_violins(v3_t, plot_dir = paste(plot_loc_violin, 'v3_t_20210914_int_rna_data_', sep = ''))
plot_celltype_markers(v3_t, plot_dir = paste(plot_loc_feature, 'v3_t_20210914_int_rna_data_', sep = ''))



# and of course, using the old style of integration
object_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds', sep = '')
object_loc_v3_t_int_perpat <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106_T_int_perpat.rds', sep = '')
# read the object
v3 <- readRDS(object_loc_v3)
# subset to what we are interested in
v3_t <- v3[, as.character(v3@meta.data$cell_type_lowerres) %in% c('CD4T', 'CD8T')]
# integrate and recluster
v3_t <- split_integrate_and_cluster_old(v3_t, 'timepoint')
# save the result
saveRDS(v3_t, object_loc_v3_t_int_perpat)
# create the plots that show the measure of integration
plot_dimplots(seurat_object = v3_t[['CA']], paste(plot_loc_dim, 'v3_t_CA_20201106_int_rna_data_', sep = ''), grouping_columns = c('seurat_clusters', 'timepoint'))
plot_dimplots(seurat_object = v3_t[['MTB']], paste(plot_loc_dim, 'v3_t_MTB_20201106_int_rna_data_', sep = ''), grouping_columns = c('seurat_clusters', 'timepoint'))
plot_dimplots(seurat_object = v3_t[['PA']], paste(plot_loc_dim, 'v3_t_PA_20201106_int_rna_data_', sep = ''), grouping_columns = c('seurat_clusters', 'timepoint'))
# plot the cell type violins as well
plot_celltype_violins(v3_t[['CA']], plot_dir = paste(plot_loc_violin, 'v3_t_CA_20201106_int_rna_data_', sep = ''))
plot_celltype_violins(v3_t[['MTB']], plot_dir = paste(plot_loc_violin, 'v3_t_CA_20201106_int_rna_data_', sep = ''))
plot_celltype_violins(v3_t[['PA']], plot_dir = paste(plot_loc_violin, 'v3_t_CA_20201106_int_rna_data_', sep = ''))
# and finally the feature plots
plot_celltype_markers(v3_t[['CA']], plot_dir = paste(plot_loc_feature, 'v3_t_CA_20201106_int_rna_data_', sep = ''))
plot_celltype_markers(v3_t[['MTB']], plot_dir = paste(plot_loc_feature, 'v3_t_MTB_20201106_int_rna_data_', sep = ''))
plot_celltype_markers(v3_t[['PA']], plot_dir = paste(plot_loc_feature, 'v3_t_PA_20201106_int_rna_data_', sep = ''))
# check the old assignments as well
plot_dimplots(seurat_object = v3_t[['CA']], paste(plot_loc_dim, 'v3_t_CA_20201106_int_rna_data_old_', sep = ''), grouping_columns = c('cell_type', 'cell_type_lowerres'))
plot_dimplots(seurat_object = v3_t[['MTB']], paste(plot_loc_dim, 'v3_t_MTB_20201106_int_rna_data_old_', sep = ''), grouping_columns = c('cell_type', 'cell_type_lowerres'))
plot_dimplots(seurat_object = v3_t[['PA']], paste(plot_loc_dim, 'v3_t_PA_20201106_int_rna_data_old_', sep = ''), grouping_columns = c('cell_type', 'cell_type_lowerres'))

# read the Azimuth prediction we did before
object_loc_v3_azt <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106_azimuth.rds', sep = '')
v3_predicted <- readRDS(object_loc_v3_azt)

# add the Azimuth labels
v3_t[['CA']] <- AddMetaData(v3_t[['CA']], v3_predicted@meta.data['predicted.celltype.l2'])
v3_t[['MTB']] <- AddMetaData(v3_t[['MTB']], v3_predicted@meta.data['predicted.celltype.l2'])
v3_t[['PA']] <- AddMetaData(v3_t[['PA']], v3_predicted@meta.data['predicted.celltype.l2'])

# these are the labels that are possible
labels_t_azimuth <- c('CD8 Naive', 'CD4 CTL', 'CD4 Naive', 'CD4 TCM', 'CD8 TEM', 'MAIT', 'CD8 TEM', 'dnT', 'CD4 TEM', 'CD8 TCM', 'CD8 Proliferating')
labels_t_highres <- c('th2 CD4T', 'naive CD4T', 'th1 CD4T', 'memory CD8T', 'naive CD8T', 'naive CD4T transitioning to stim', 'reg CD4T', 'memory CD8T left and naive CD8T right', 'cyto CD4T', 'double negative T', 'memory CD8T', 'memory CD8T left', 'T helper')

v3_t[['CA']]@meta.data$cell_type_t <- v3_t[['CA']]@meta.data$cell_type
levels(v3_t[['CA']]@meta.data$cell_type_t) <- c(levels(v3_t[['CA']]@meta.data$cell_type_t), 'other')
v3_t[['CA']]@meta.data[!(v3_t[['CA']]@meta.data$cell_type_t %in% labels_t_highres), 'cell_type_t'] <- 'other'
v3_t[['CA']]@meta.data$predicted.celltype.l2.t <- v3_t[['CA']]@meta.data$predicted.celltype.l2
levels(v3_t[['CA']]@meta.data$predicted.celltype.l2.t) <- c(levels(v3_t[['CA']]@meta.data$predicted.celltype.l2.t), 'other')
v3_t[['CA']]@meta.data[!(v3_t[['CA']]@meta.data$predicted.celltype.l2.t %in% labels_t_azimuth), 'predicted.celltype.l2.t'] <- 'other'

v3_t[['MTB']]@meta.data$cell_type_t <- v3_t[['MTB']]@meta.data$cell_type
levels(v3_t[['MTB']]@meta.data$cell_type_t) <- c(levels(v3_t[['MTB']]@meta.data$cell_type_t), 'other')
v3_t[['MTB']]@meta.data[!(v3_t[['MTB']]@meta.data$cell_type_t %in% labels_t_highres), 'cell_type_t'] <- 'other'
v3_t[['MTB']]@meta.data$predicted.celltype.l2.t <- v3_t[['MTB']]@meta.data$predicted.celltype.l2
levels(v3_t[['MTB']]@meta.data$predicted.celltype.l2.t) <- c(levels(v3_t[['MTB']]@meta.data$predicted.celltype.l2.t), 'other')
v3_t[['MTB']]@meta.data[!(v3_t[['MTB']]@meta.data$predicted.celltype.l2.t %in% labels_t_azimuth), 'predicted.celltype.l2.t'] <- 'other'

v3_t[['PA']]@meta.data$cell_type_t <- v3_t[['PA']]@meta.data$cell_type
levels(v3_t[['PA']]@meta.data$cell_type_t) <- c(levels(v3_t[['PA']]@meta.data$cell_type_t), 'other')
v3_t[['PA']]@meta.data[!(v3_t[['PA']]@meta.data$cell_type_t %in% labels_t_highres), 'cell_type_t'] <- 'other'
v3_t[['PA']]@meta.data$predicted.celltype.l2.t <- v3_t[['PA']]@meta.data$predicted.celltype.l2
levels(v3_t[['PA']]@meta.data$predicted.celltype.l2.t) <- c(levels(v3_t[['PA']]@meta.data$predicted.celltype.l2.t), 'other')
v3_t[['PA']]@meta.data[!(v3_t[['PA']]@meta.data$predicted.celltype.l2.t %in% labels_t_azimuth), 'predicted.celltype.l2.t'] <- 'other'


v3_t_ca_p1 <- DimPlot(v3_t[['CA']], group.by = 'cell_type_t')
v3_t_ca_p2 <- DimPlot(v3_t[['CA']], group.by = 'predicted.celltype.l2.t')
plot_grid(v3_t_ca_p1, v3_t_ca_p2, ncol=1, nrow=2)
ggsave(paste(plot_loc_dim, 'v3_t_CA_20201106_int_rna_data_cttvspcl2t.pdf', sep = ''), width = 20, height = 20)
v3_t_mtb_p1 <- DimPlot(v3_t[['MTB']], group.by = 'cell_type_t')
v3_t_mtb_p2 <- DimPlot(v3_t[['MTB']], group.by = 'predicted.celltype.l2.t')
plot_grid(v3_t_mtb_p1, v3_t_mtb_p2, ncol=1, nrow=2)
ggsave(paste(plot_loc_dim, 'v3_t_MTB_20201106_int_rna_data_cttvspcl2t.pdf', sep = ''), width = 20, height = 20)

# classify the T cells, by assigning the clusters based on the cell type with the largest proportion in that cluster
v3_t[['CA']] <- add_imputed_meta_data(v3_t[['CA']], 'seurat_clusters', 'predicted.celltype.l2.t', 'clustered.celltype.l2.t')
v3_t[['MTB']] <- add_imputed_meta_data(v3_t[['MTB']], 'seurat_clusters', 'predicted.celltype.l2.t', 'clustered.celltype.l2.t')
v3_t[['PA']] <- add_imputed_meta_data(v3_t[['PA']], 'seurat_clusters', 'predicted.celltype.l2.t', 'clustered.celltype.l2.t')

# create a table with just the predictions, making sure we don't have any double entries, and not overwriting the T classifications
v3_t_table <- v3_t[['CA']]@meta.data['clustered.celltype.l2.t']
v3_t_table <- rbind(v3_t_table, v3_t[['MTB']]@meta.data[v3_t[['MTB']]@meta.data[['clustered.celltype.l2.t']] != 'other' & !(rownames(v3_t[['MTB']]@meta.data) %in% rownames(v3_t_table)), ]['clustered.celltype.l2.t'])
v3_t_table <- rbind(v3_t_table, v3_t[['PA']]@meta.data[v3_t[['PA']]@meta.data[['clustered.celltype.l2.t']] != 'other' & !(rownames(v3_t[['PA']]@meta.data) %in% rownames(v3_t_table)), ]['clustered.celltype.l2.t'])
# add the cell barcodes as a column (I just like it more that way)
v3_t_table$barcode <- rownames(v3_t_table)
# and reorder the way we like
v3_t_table <- v3_t_table[, c('barcode', 'clustered.celltype.l2.t')]
# write the result somewhere
t_classification_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/seurat_multimodal/classifications/'
t_classification_loc_v3 <- paste(t_classification_loc, 'v3_azi_to_cluster_l2_cell_types.tsv', sep = '')
write.table(v3_t_table, t_classification_loc_v3, sep = '\t', row.names = F, col.names = T, quote = F)


# now everything for v2 as well
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds', sep = '')
object_loc_v2_t_int_perpat <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201029_T_int_perpat.rds', sep = '')
v2 <- readRDS(object_loc_v2)
# subset to what we are interested in
v2_t <- v2[, as.character(v2@meta.data$cell_type_lowerres) %in% c('CD4T', 'CD8T')]
# clear memory
rm(v2)
# integrate and recluster
v2_t <- split_integrate_and_cluster_old(v2_t, 'timepoint')
# save the result
saveRDS(v2_t, object_loc_v2_t_int_perpat)

# read the Azimuth prediction we did before
object_loc_v2_azt <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029_azimuth.rds', sep = '')
v2_predicted <- readRDS(object_loc_v2_azt)

# add the Azimuth labels
v2_t[['CA']] <- AddMetaData(v2_t[['CA']], v2_predicted@meta.data['predicted.celltype.l2'])
v2_t[['MTB']] <- AddMetaData(v2_t[['MTB']], v2_predicted@meta.data['predicted.celltype.l2'])
v2_t[['PA']] <- AddMetaData(v2_t[['PA']], v2_predicted@meta.data['predicted.celltype.l2'])

# clear memory
rm(v2_predicted)

v2_t[['CA']]@meta.data$predicted.celltype.l2.t <- v2_t[['CA']]@meta.data$predicted.celltype.l2
levels(v2_t[['CA']]@meta.data$predicted.celltype.l2.t) <- c(levels(v2_t[['CA']]@meta.data$predicted.celltype.l2.t), 'other')
v2_t[['CA']]@meta.data[!(v2_t[['CA']]@meta.data$predicted.celltype.l2.t %in% labels_t_azimuth), 'predicted.celltype.l2.t'] <- 'other'
v2_t[['MTB']]@meta.data$predicted.celltype.l2.t <- v2_t[['MTB']]@meta.data$predicted.celltype.l2
levels(v2_t[['MTB']]@meta.data$predicted.celltype.l2.t) <- c(levels(v2_t[['MTB']]@meta.data$predicted.celltype.l2.t), 'other')
v2_t[['MTB']]@meta.data[!(v2_t[['MTB']]@meta.data$predicted.celltype.l2.t %in% labels_t_azimuth), 'predicted.celltype.l2.t'] <- 'other'
v2_t[['PA']]@meta.data$predicted.celltype.l2.t <- v2_t[['PA']]@meta.data$predicted.celltype.l2
levels(v2_t[['PA']]@meta.data$predicted.celltype.l2.t) <- c(levels(v2_t[['PA']]@meta.data$predicted.celltype.l2.t), 'other')
v2_t[['PA']]@meta.data[!(v2_t[['PA']]@meta.data$predicted.celltype.l2.t %in% labels_t_azimuth), 'predicted.celltype.l2.t'] <- 'other'


# these are the labels that are possible
labels_t_azimuth <- c('CD8 Naive', 'CD4 CTL', 'CD4 Naive', 'CD4 TCM', 'CD8 TEM', 'MAIT', 'CD8 TEM', 'dnT', 'CD4 TEM', 'CD8 TCM', 'CD8 Proliferating')

# classify the T cells, by assigning the clusters based on the cell type with the largest proportion in that cluster
v2_t[['CA']] <- add_imputed_meta_data(v2_t[['CA']], 'seurat_clusters', 'predicted.celltype.l2.t', 'clustered.celltype.l2.t')
v2_t[['MTB']] <- add_imputed_meta_data(v2_t[['MTB']], 'seurat_clusters', 'predicted.celltype.l2.t', 'clustered.celltype.l2.t')
v2_t[['PA']] <- add_imputed_meta_data(v2_t[['PA']], 'seurat_clusters', 'predicted.celltype.l2.t', 'clustered.celltype.l2.t')

# create a table with just the predictions, making sure we don't have any double entries, and not overwriting the T classifications
v2_t_table <- v2_t[['CA']]@meta.data['clustered.celltype.l2.t']
v2_t_table <- rbind(v2_t_table, v2_t[['MTB']]@meta.data[v2_t[['MTB']]@meta.data[['clustered.celltype.l2.t']] != 'other' & !(rownames(v2_t[['MTB']]@meta.data) %in% rownames(v2_t_table)), ]['clustered.celltype.l2.t'])
v2_t_table <- rbind(v2_t_table, v2_t[['PA']]@meta.data[v2_t[['PA']]@meta.data[['clustered.celltype.l2.t']] != 'other' & !(rownames(v2_t[['PA']]@meta.data) %in% rownames(v2_t_table)), ]['clustered.celltype.l2.t'])
# add the cell barcodes as a column (I just like it more that way)
v2_t_table$barcode <- rownames(v2_t_table)
# and reorder the way we like
v2_t_table <- v2_t_table[, c('barcode', 'clustered.celltype.l2.t')]
# write the result somewhere
t_classification_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/seurat_multimodal/classifications/'
t_classification_loc_v2 <- paste(t_classification_loc, 'v2_azi_to_cluster_l2_cell_types.tsv', sep = '')
write.table(v2_t_table, t_classification_loc_v2, sep = '\t', row.names = F, col.names = T, quote = F)




