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
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object[[grouping_column]]
      LabelClusters(plot = p, id = grouping_column)
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
    rm(seurat_object)
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

# read the object
v3 <- readRDS(object_loc_v3)

# subset and recluster
v3_t <- subset_and_recluster(v3, 'cell_type_lowerres', c('CD4T', 'CD8T'))
# save the object
saveRDS(v3_t, object_loc_v3_t)

v3_t  <- NormalizeData(v3_t)
plot_celltype_violins(v3_t, plot_dir = paste(plot_loc_violin, 'v3_t_20210909_rna_data_', sep = ''))



# let's try again, but this time with integration
object_loc <- '/data/p287578/1M_cells_scRNAseq/seurat_preprocess_samples/objects/'
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds', sep = '')
object_loc_v3_t_int <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106_T_int.rds', sep = '')
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
plot_dimplots(seurat_object = v3_t[['CA']], paste(plot_loc_dim, 'v3_t_CA_20210914_int_rna_data_', sep = ''), grouping_columns = c('seurat_clusters', 'timepoint'))
plot_dimplots(seurat_object = v3_t[['MTB']], paste(plot_loc_dim, 'v3_t_MTB_20210914_int_rna_data_', sep = ''), grouping_columns = c('seurat_clusters', 'timepoint'))
plot_dimplots(seurat_object = v3_t[['PA']], paste(plot_loc_dim, 'v3_t_PA_20210914_int_rna_data_', sep = ''), grouping_columns = c('seurat_clusters', 'timepoint'))
# plot the cell type violins as well
plot_celltype_violins(v3_t[['CA']], plot_dir = paste(plot_loc_violin, 'v3_t_CA_20210909_int_rna_data_', sep = ''))
plot_celltype_violins(v3_t[['MTB']], plot_dir = paste(plot_loc_violin, 'v3_t_MTB_20210909_int_rna_data_', sep = ''))
plot_celltype_violins(v3_t[['PA']], plot_dir = paste(plot_loc_violin, 'v3_t_PA_20210909_int_rna_data_', sep = ''))
# and finally the feature plots
plot_celltype_markers(v3_t[['CA']], plot_dir = paste(plot_loc_feature, 'v3_t_CA_20210909_int_rna_data_', sep = ''))
plot_celltype_markers(v3_t[['MTB']], plot_dir = paste(plot_loc_feature, 'v3_t_MTB_20210909_int_rna_data_', sep = ''))
plot_celltype_markers(v3_t[['PA']], plot_dir = paste(plot_loc_feature, 'v3_t_PA_20210909_int_rna_data_', sep = ''))



