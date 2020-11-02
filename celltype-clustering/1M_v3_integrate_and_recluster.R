library(Seurat)
library(ggplot2)

# add metadata that is based on existing incomplete metadata in the seurat object
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


# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')
# plot locations
dimplot_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/plots/dimplots/'
v3_CA_dimplot_ct_loc <- paste(dimplot_loc, 'v3_CA_reclus_origct_20201029.png', sep = '')
v3_MTB_dimplot_ct_loc <- paste(dimplot_loc, 'v3_MTB_reclus_origct_20201029.png', sep = '')
v3_PA_dimplot_ct_loc <- paste(dimplot_loc, 'v3_PA_reclus_origct_20201029.png', sep = '')

# read the v3 object
v3 <- readRDS(object_loc_v3)
# grab the UT timepoint we will use for all integrations
v3_UT <- subset(v3, subset = timepoint == "UT")
# grab the stimmed unstimmed samples
v3_UT <- v3_UT[, !(v3_UT@meta.data$assignment %in% c('LLDeep_1058','LLDeep_1229', 'LLDeep_1179', 'LLDeep_1247', 'LLDeep_1016', 'LLDeep_1067', 'LLDeep_0747', 'LLDeep_0906'))]
# try to do merge the CA timepoints with the untreated samples
v3_X3hCA <- subset(v3, subset = timepoint == "X3hCA")
v3_X24hCA <- subset(v3, subset = timepoint == "X24hCA")
# get a list of the untreated and the two CA conditions
v3_CA_list <- list(v3_UT,v3_X3hCA,v3_X24hCA)
# grab the variable features in each condition
v3_CA_list <- lapply(X = v3_CA_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# calculate the anchors
v3_CA_anchors <- FindIntegrationAnchors(object.list = v3_CA_list, dims = 1:20)
# merge the untreated with the CA conditions based on the anchors
v3_CA <- IntegrateData(anchorset = v3_CA_anchors, dims = 1:20)
# set the default assay to the integrated one for clustering etc
DefaultAssay(v3_CA) <- "integrated"
# do scaling, required for integrated sets
v3_CA <- ScaleData(v3_CA, verbose = FALSE)
# do PCA
v3_CA <- RunPCA(v3_CA, npcs = 30, verbose = FALSE)
# UMAP and Clustering
v3_CA <- RunUMAP(v3_CA, reduction = "pca", dims = 1:20)
v3_CA <- FindNeighbors(v3_CA, reduction = "pca", dims = 1:20)
v3_CA <- FindClusters(v3_CA, resolution = 0.5)
# clear up some memory
rm(v3_X3hCA)
rm(v24_X3hCA)
# save combined file
saveRDS(v3_CA, "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_CA_cca_integrated_classicnormscale_20201029.rds")
# make the ct plot
DimPlot(v3_CA, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(v3_CA_dimplot_ct_loc, width = 10, height = 10)
# clear memory after saving
rm(v3_CA)

# try to do merge the PA timepoints with the untreated samples 
v3_X3hPA <- subset(v3, subset = timepoint == "X3hPA")
v3_X24hPA <- subset(v3, subset = timepoint == "X24hPA")
# get a list of the untreated and the two CA conditions
v3_PA_list <- list(v3_UT,v3_X3hPA,v3_X24hPA)
# grab the variable features in each condition
v3_PA_list <- lapply(X = v3_PA_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# calculate the anchors
v3_PA_anchors <- FindIntegrationAnchors(object.list = v3_PA_list, dims = 1:20)
# merge the untreated with the CA conditions based on the anchors
v3_PA <- IntegrateData(anchorset = v3_PA_anchors, dims = 1:20)
# set the default assay to the integrated one for clustering etc
DefaultAssay(v3_PA) <- "integrated"
# do scaling, required for integrated sets
v3_PA <- ScaleData(v3_PA, verbose = FALSE)
# do PCA
v3_PA <- RunPCA(v3_PA, npcs = 30, verbose = FALSE)
# UMAP and Clustering
v3_PA <- RunUMAP(v3_PA, reduction = "pca", dims = 1:20)
v3_PA <- FindNeighbors(v3_PA, reduction = "pca", dims = 1:20)
v3_PA <- FindClusters(v3_PA, resolution = 0.5)
# clear up some memory
rm(v3_X3hPA)
rm(v24_X3hPA)
# save combined file
saveRDS(v3_PA, "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_PA_cca_integrated_classicnormscale_20201029.rds")
# make the ct plot
DimPlot(v3_PA, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(v3_PA_dimplot_ct_loc, width = 10, height = 10)
# clear memory after saving
rm(v3_PA)

# try to do merge the PA timepoints with the untreated samples 
v3_X3hMTB <- subset(v3, subset = timepoint == "X3hMTB")
v3_X24hMTB <- subset(v3, subset = timepoint == "X24hMTB")
# get a list of the untreated and the two CA conditions
v3_MTB_list <- list(v3_UT,v3_X3hMTB,v3_X24hMTB)
# grab the variable features in each condition
v3_MTB_list <- lapply(X = v3_MTB_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# calculate the anchors
v3_MTB_anchors <- FindIntegrationAnchors(object.list = v3_MTB_list, dims = 1:20)
# merge the untreated with the CA conditions based on the anchors
v3_MTB <- IntegrateData(anchorset = v3_MTB_anchors, dims = 1:20)
# set the default assay to the integrated one for clustering etc
DefaultAssay(v3_MTB) <- "integrated"
# do scaling, required for integrated sets
v3_MTB <- ScaleData(v3_MTB, verbose = FALSE)
# do PCA
v3_MTB <- RunPCA(v3_MTB, npcs = 30, verbose = FALSE)
# UMAP and Clustering
v3_MTB <- RunUMAP(v3_MTB, reduction = "pca", dims = 1:20)
v3_MTB <- FindNeighbors(v3_MTB, reduction = "pca", dims = 1:20)
v3_MTB <- FindClusters(v3_MTB, resolution = 0.5)
# clear up some memory
rm(v3_X3hMTB)
rm(v24_X3hMTB)
# save combined file
saveRDS(v3_MTB, "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_MTB_cca_integrated_classicnormscale_20201029.rds")
# make the ct plot
DimPlot(v3_MTB, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(v3_MTB_dimplot_ct_loc, width = 10, height = 10)
# clear memory after saving
rm(v3_MTB)

# read the CA integrated object
v3_CA <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_CA_cca_integrated_classicnormscale_20201029.rds")
# create the clusters plot
DimPlot(v3_CA, reduction = 'umap', group.by = 'seurat_clusters')
ggsave(paste(dimplot_loc, 'v3_CA_reclus_clus_20201029.png', sep = ''), width = 10, height = 10)
# backup current cell type
v3_CA@meta.data$cell_type_lowerres_old <- v3_CA@meta.data$cell_type_lowerres
v3_CA@meta.data$cell_type_lowerres <- NULL
# add imputed cell types
v3_CA <- add_imputed_meta_data(v3_CA, 'seurat_clusters', 'cell_type_lowerres_old', 'cell_type_lowerres')
# make new plot
DimPlot(v3_CA, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(paste(dimplot_loc, 'v3_CA_reclus_ctnew_20201029.png', sep = ''), width = 10, height = 10)
# save new object
saveRDS(v3_CA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_CA_cca_integrated_classicnormscale_newct_20201029.rds")

# read the MTB integrated object
v3_MTB <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_MTB_cca_integrated_classicnormscale_20201029.rds")
# create the clusters plot
DimPlot(v3_MTB, reduction = 'umap', group.by = 'seurat_clusters')
ggsave(paste(dimplot_loc, 'v3_MTB_reclus_clus_20201029.png', sep = ''), width = 10, height = 10)
# backup current cell type
v3_MTB@meta.data$cell_type_lowerres_old <- v3_MTB@meta.data$cell_type_lowerres
v3_MTB@meta.data$cell_type_lowerres <- NULL
# add imputed cell types
v3_MTB <- add_imputed_meta_data(v3_MTB, 'seurat_clusters', 'cell_type_lowerres_old', 'cell_type_lowerres')
# make new plot
DimPlot(v3_MTB, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(paste(dimplot_loc, 'v3_MTB_reclus_ctnew_20201029.png', sep = ''), width = 10, height = 10)
# save new object
saveRDS(v3_MTB, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_MTB_cca_integrated_classicnormscale_newct_20201029.rds")


# read the PA integrated object
v3_PA <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_PA_cca_integrated_classicnormscale_20201029.rds")
# create the clusters plot
DimPlot(v3_PA, reduction = 'umap', group.by = 'seurat_clusters')
ggsave(paste(dimplot_loc, 'v3_PA_reclus_clus_20201029.png', sep = ''), width = 10, height = 10)
# backup current cell type
v3_PA@meta.data$cell_type_lowerres_old <- v3_PA@meta.data$cell_type_lowerres
v3_PA@meta.data$cell_type_lowerres <- NULL
# add imputed cell types
v3_PA <- add_imputed_meta_data(v3_PA, 'seurat_clusters', 'cell_type_lowerres_old', 'cell_type_lowerres')
# make new plot
DimPlot(v3_PA, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(paste(dimplot_loc, 'v3_PA_reclus_ctnew_20201029.png', sep = ''), width = 10, height = 10)
# save new object
saveRDS(v3_PA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_PA_cca_integrated_classicnormscale_newct_20201029.rds")

# create the cell type tsv
v3_CA_cell_types <- v3_CA@meta.data['cell_type_lowerres']
v3_MTB_cell_types <- v3_MTB@meta.data['cell_type_lowerres']
v3_PA_cell_types <- v3_PA@meta.data['cell_type_lowerres']
v3_cell_types <- rbind(v3_CA_cell_types, v3_MTB_cell_types)
v3_cell_types <- rbind(v3_cell_types, v3_PA_cell_types)
# read the v3 object to put the identitities in
v3 <- readRDS(object_loc_v3)
v3@meta.data$cell_type_lowerres_old <- v3@meta.data$cell_type_lowerres
v3@meta.data$cell_type_lowerres <- NULL
# add the cell types
v3 <- AddMetaData(v3, v3_cell_types)
# remove the UT samples that were incorrectly stimulated
v3 <- v3[, !(v3@meta.data$assignment %in% c('LLDeep_1058','LLDeep_1229', 'LLDeep_1179', 'LLDeep_1247', 'LLDeep_1016', 'LLDeep_1067', 'LLDeep_0747', 'LLDeep_0906') & v3@meta.data$timepoint == 'UT')]
# redo normalization for RNA assay, since those are gone
DefaultAssay(v3) <- 'RNA'                                                  
v3 <- NormalizeData(v3)
# write the result
saveRDS(v3, '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201029.rds')

