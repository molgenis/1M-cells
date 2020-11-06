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

# create the feature plots
plot_celltype_markers <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$seurat_clusters 
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}

# create the feature plots
plot_celltype_markers_param <- function(seurat_object, markers, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes <- markers
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$seurat_clusters
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}

# create violin plots
plot_celltype_violins <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes <- c("CCR7","S100A4","CD3E","CD4","CD8A", "CD8B", "FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = seurat_object$integrated_snn_res.1.2)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10, assay=assay, slot=slot)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}




# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')
# plot locations
dimplot_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/plots/dimplots/'
v3_CA_dimplot_ct_loc <- paste(dimplot_loc, 'v3_CA_reclus_origct_20201105.png', sep = '')
v3_MTB_dimplot_ct_loc <- paste(dimplot_loc, 'v3_MTB_reclus_origct_20201105.png', sep = '')
v3_PA_dimplot_ct_loc <- paste(dimplot_loc, 'v3_PA_reclus_origct_20201105.png', sep = '')

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
v3_CA <- RunUMAP(v3_CA, reduction = "pca", dims = 1:30)
v3_CA <- FindNeighbors(v3_CA, reduction = "pca", dims = 1:30)
v3_CA <- FindClusters(v3_CA, resolution = 1.2)
# clear up some memory
rm(v3_X3hCA)
rm(v24_X3hCA)
# save combined file
saveRDS(v3_CA, "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_CA_cca_integrated_classicnormscale_20201105.rds")
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
v3_PA <- RunUMAP(v3_PA, reduction = "pca", dims = 1:30)
v3_PA <- FindNeighbors(v3_PA, reduction = "pca", dims = 1:30)
v3_PA <- FindClusters(v3_PA, resolution = 1.2)
# clear up some memory
rm(v3_X3hPA)
rm(v24_X3hPA)
# save combined file
saveRDS(v3_PA, "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_PA_cca_integrated_classicnormscale_20201105.rds")
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
v3_MTB <- RunUMAP(v3_MTB, reduction = "pca", dims = 1:30)
v3_MTB <- FindNeighbors(v3_MTB, reduction = "pca", dims = 1:30)
v3_MTB <- FindClusters(v3_MTB, resolution = 1.2)
# clear up some memory
rm(v3_X3hMTB)
rm(v24_X3hMTB)
# save combined file
saveRDS(v3_MTB, "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_MTB_cca_integrated_classicnormscale_20201105.rds")
# clear memory after saving
rm(v3_MTB)


# do T reclustering based on 20201105 1.2 resolution output
#v3_3hMTB_T <- subset(v3_MTB, subset = (cell_type_lowerres == 'CD4T' | cell_type_lowerres == 'CD8T') & timepoint == 'X3hMTB')
#v3_24hMTB_T <- subset(v3_MTB, subset = (cell_type_lowerres == 'CD4T' | cell_type_lowerres == 'CD8T') & timepoint == 'X24hMTB')
#v3_MTB_UT_T <- subset(v3_MTB, subset = (cell_type_lowerres == 'CD4T' | cell_type_lowerres == 'CD8T') & timepoint == 'UT')
#v3_3hCA_T <- subset(v3_CA, subset = (cell_type_lowerres == 'CD4T' | cell_type_lowerres == 'CD8T') & timepoint == 'X3hCA')
#v3_24hCA_T <- subset(v3_CA, subset = (cell_type_lowerres == 'CD4T' | cell_type_lowerres == 'CD8T') & timepoint == 'X24hCA')
#v3_CA_UT_T <- subset(v3_CA, subset = (cell_type_lowerres == 'CD4T' | cell_type_lowerres == 'CD8T') & timepoint == 'UT')
#v3_3hPA_T <- subset(v3_PA, subset = (cell_type_lowerres == 'CD4T' | cell_type_lowerres == 'CD8T') & timepoint == 'X3hPA')
#v3_24hPA_T <- subset(v3_PA, subset = (cell_type_lowerres == 'CD4T' | cell_type_lowerres == 'CD8T') & timepoint == 'X24hPA')
#v3_PA_UT_T <- subset(v3_PA, subset = (cell_type_lowerres == 'CD4T' | cell_type_lowerres == 'CD8T') & timepoint == 'UT')
# grab T cells from any integration with UT
#v3_UT_T <- v3_MTB_UT_T
#v3_CA_UT_T_notinmtb <- v3_CA_UT_T[, !(rownames(v3_CA_UT_T@meta.data) %in% rownames(v3_UT_T@meta.data))]
# the SCT assay seems to give issues, so just NULLing that one
#v3_UT_T@assays$SCT <- NULL
#v3_CA_UT_T@assays$SCT <- NULL
#v3_UT_T <- merge(v3_CA_UT_T_notinmtb, v3_UT_T)
#v3_PA_UT_T_notinmtborca <- v3_PA_UT_T[, !(rownames(v3_PA_UT_T@meta.data) %in% rownames(v3_UT_T@meta.data))]
#v3_PA_UT_T_notinmtborca@assays$SCT <- NULL
#v3_UT_T <- merge(v3_PA_UT_T_notinmtborca, v3_UT_T)
# add everything into a list
#v3_T_list <- list(v3_UT_T, v3_3hMTB_T, v3_24hMTB_T, v3_3hCA_T, v3_24hCA_T, v3_3hPA_T, v3_24hPA_T)
#v3_T_list <- lapply(X = v3_T_list, FUN = function(x) {
#  DefaultAssay(x) <- "RNA"
#  x <- NormalizeData(x)
#  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#})
#saveRDS(v3_T_list, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_T_list_20201105.rds")
# calculate the anchors
#v3_T_anchors <- FindIntegrationAnchors(object.list = v3_T_list, dims = 1:20)
# merge the untreated with the CA conditions based on the anchors
#v3_T <- IntegrateData(anchorset = v3_T_anchors, dims = 1:20)
# set the default assay to the integrated one for clustering etc
#DefaultAssay(v3_T) <- "integrated"
# do scaling, required for integrated sets
#v3_T <- ScaleData(v3_T, verbose = FALSE)
# do PCA
#v3_T <- RunPCA(v3_T, npcs = 30, verbose = FALSE)
# UMAP and Clustering
#v3_T <- RunUMAP(v3_T, reduction = "pca", dims = 1:30)
#v3_T <- FindNeighbors(v3_T, reduction = "pca", dims = 1:30)
#v3_T <- FindClusters(v3_T, resolution = 1.2)
# create the plot
#DimPlot(v3_T, reduction = 'umap', group.by = 'seurat_clusters')
#ggsave(paste(dimplot_loc, 'v3_T_reclus_clus_20201105.png', sep = ''), width = 10, height = 10)
# switch to RNA assay to prepare for featureplots
#DefaultAssay(v3_T) <- 'RNA'
#v3_T <- NormalizeData(v3_CA)
#saveRDS(v3_T, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_T_cca_integrated_classicnorm_20201105.rds")
# plot
#plot_celltype_markers(seurat_object = v3_T, assay = 'RNA', slot = 'data', plot_dir = paste(features_base_plot, 'v3_T/RNA_data_20201105/', sep = ''))

v3_CA <- readRDS('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_CA_cca_integrated_classicnormscale_20201105.rds')
DefaultAssay(v3_CA) <- 'integrated'
v3_CA <- FindClusters(v3_CA, resolution = 1.5)
DimPlot(v3_CA)
dimplot_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/plots/dimplots/'
ggsave(paste(dimplot_loc, 'v3_CA_reclus_clusts_20201106.png', sep=''), width = 10, height = 10)
DefaultAssay(v3_CA) <- 'RNA'
v3_CA <- NormalizeData(v3_CA)
plot_celltype_markers_param(seurat_object = v3_CA, markers=c('SIRPA', 'ITGAM'),  assay = 'RNA', slot = 'data', plot_dir = paste(features_base_plot, 'v3_CA/RNA_data_20201106/', sep = ''))
# backup current cell type
v3_CA@meta.data$cell_type_lowerres_old <- v3_CA@meta.data$cell_type_lowerres
v3_CA@meta.data$cell_type_lowerres <- NULL
# add imputed cell types
v3_CA <- add_imputed_meta_data(v3_CA, 'seurat_clusters', 'cell_type_lowerres_old', 'cell_type_lowerres')
# make new plot
DimPlot(v3_CA, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(paste(dimplot_loc, 'v3_CA_reclus_ctnew_20201106.png', sep = ''), width = 10, height = 10)
plot_celltype_markers(seurat_object = v3_CA,  assay = 'RNA', slot = 'data', plot_dir = paste(features_base_plot, 'v3_CA/RNA_data_20201106/', sep = ''))


v3_PA <- readRDS('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_PA_cca_integrated_classicnormscale_20201105.rds')
DefaultAssay(v3_PA) <- 'integrated'
v3_PA <- FindClusters(v3_PA, resolution = 1.5)
DimPlot(v3_PA)
dimplot_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/plots/dimplots/'
ggsave(paste(dimplot_loc, 'v3_PA_reclus_clusts_20201106.png', sep=''), width = 10, height = 10)
DefaultAssay(v3_PA) <- 'RNA'
v3_PA <- NormalizeData(v3_PA)
plot_celltype_markers_param(seurat_object = v3_PA, markers=c('SIRPA', 'ITGAM'),  assay = 'RNA', slot = 'data', plot_dir = paste(features_base_plot, 'v3_PA/RNA_data_20201106/', sep = ''))
# backup current cell type
v3_PA@meta.data$cell_type_lowerres_old <- v3_PA@meta.data$cell_type_lowerres
v3_PA@meta.data$cell_type_lowerres <- NULL
# add imputed cell types
v3_PA <- add_imputed_meta_data(v3_PA, 'seurat_clusters', 'cell_type_lowerres_old', 'cell_type_lowerres')
# make new plot
DimPlot(v3_PA, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(paste(dimplot_loc, 'v3_PA_reclus_ctnew_20201106.png', sep = ''), width = 10, height = 10)
plot_celltype_markers(seurat_object = v3_PA,  assay = 'RNA', slot = 'data', plot_dir = paste(features_base_plot, 'v3_PA/RNA_data_20201106/', sep = ''))


v3_MTB <- readRDS('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_MTB_cca_integrated_classicnormscale_20201105.rds')
DefaultAssay(v3_MTB) <- 'integrated'
v3_MTB <- FindClusters(v3_MTB, resolution = 1.5)
DimPlot(v3_MTB)
dimplot_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/plots/dimplots/'
ggsave(paste(dimplot_loc, 'v3_MTB_reclus_clusts_20201106.png', sep=''), width = 10, height = 10)
DefaultAssay(v3_MTB) <- 'RNA'
v3_MTB <- NormalizeData(v3_MTB)
plot_celltype_markers_param(seurat_object = v3_MTB, markers=c('SIRPA', 'ITGAM'),  assay = 'RNA', slot = 'data', plot_dir = paste(features_base_plot, 'v3_MTB/RNA_data_20201106/', sep = ''))
# backup current cell type
v3_MTB@meta.data$cell_type_lowerres_old <- v3_MTB@meta.data$cell_type_lowerres
v3_MTB@meta.data$cell_type_lowerres <- NULL
# add imputed cell types
v3_MTB <- add_imputed_meta_data(v3_MTB, 'seurat_clusters', 'cell_type_lowerres_old', 'cell_type_lowerres')
# make new plot
DimPlot(v3_MTB, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(paste(dimplot_loc, 'v3_MTB_reclus_ctnew_20201106.png', sep = ''), width = 10, height = 10)
plot_celltype_markers(seurat_object = v3_MTB,  assay = 'RNA', slot = 'data', plot_dir = paste(features_base_plot, 'v3_MTB/RNA_data_20201106/', sep = ''))

# save objects
saveRDS(v3_CA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_CA_cca_integrated_classicnormscale_newct_20201106.rds")
saveRDS(v3_PA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_PA_cca_integrated_classicnormscale_newct_20201106.rds")
saveRDS(v3_MTB, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_MTB_cca_integrated_classicnormscale_newct_20201106.rds")

# backup cell types
v3_CA@meta.data$cell_type_lowerres_new <- v3_CA@meta.data$cell_type_lowerres
v3_MTB@meta.data$cell_type_lowerres_new <- v3_MTB@meta.data$cell_type_lowerres
v3_PA@meta.data$cell_type_lowerres_new <- v3_PA@meta.data$cell_type_lowerres

# just change cluster 38 to a megakaryocyte
v3_MTB@meta.data[v3_MTB@meta.data$seurat_clusters == 38, ]$cell_type_lowerres <- 'megakaryocyte'
# and 34 to CD4T
v3_MTB@meta.data[v3_MTB@meta.data$seurat_clusters == 33, ]$cell_type_lowerres <- 'CD4T'
DimPlot(v3_MTB, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(paste(dimplot_loc, 'v3_MTB_reclus_ctnewmega_20201106.png', sep = ''), width = 10, height = 10)

# just change cluster 36 to megakaryocyte
v3_CA@meta.data[v3_CA@meta.data$seurat_clusters == 36, ]$cell_type_lowerres <- 'megakaryocyte'
# determine PF4 expression 
featureplot_PF4_CA <- FeaturePlot(v3_CA, features = c("PF4"))
v3_CA <- AddMetaData(v3_CA, featureplot_PF4_CA$data["PF4"], "PF4_relative_expression")
v3_CA@meta.data$above_PF4_threshold <- F
v3_CA@meta.data[v3_CA@meta.data$PF4_relative_expression > 2, ]$above_PF4_threshold <- T
DimPlot(v3_CA, group.by='above_PF4_threshold')
ggsave(paste(dimplot_loc, 'v3_CA_reclus_PF4mega_20201106.png', sep = ''), width = 10, height = 10)
# change cell type based on the PF4 threshold
v3_CA@meta.data[v3_CA@meta.data$above_PF4_threshold == T, ]$cell_type_lowerres <- 'megakaryocyte'
DimPlot(v3_CA, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(paste(dimplot_loc, 'v3_CA_reclus_ctnewmega_20201106.png', sep = ''), width = 10, height = 10)

# determine PF4 expression 
featureplot_PF4_PA <- FeaturePlot(v3_PA, features = c("PF4"))
v3_PA <- AddMetaData(v3_PA, featureplot_PF4_PA$data["PF4"], "PF4_relative_expression")
# only 37 to mega
v3_PA@meta.data[v3_PA@meta.data$seurat_clusters == 37, ]$cell_type_lowerres <- 'megakaryocyte'
DimPlot(v3_PA, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave(paste(dimplot_loc, 'v3_PA_reclus_ctnewmega_20201106.png', sep = ''), width = 10, height = 10)

# save objects once again
saveRDS(v3_CA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_CA_cca_integrated_classicnormscale_newctmega_20201106.rds")
saveRDS(v3_PA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_PA_cca_integrated_classicnormscale_newctmega_20201106.rds")
saveRDS(v3_MTB, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_MTB_cca_integrated_classicnormscale_newctmega_20201106.rds")

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
saveRDS(v3, '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds')


