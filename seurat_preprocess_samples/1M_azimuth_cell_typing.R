####################
# Libraries        #
####################

library(SeuratDisk)
library(ggplot2)
library(patchwork)


####################
# Functions        #
####################

add_cell_type_predictions <- function(reference, query, normalization.method='SCT', reference.reduction="spca", dims=1:50, recompute.residuals=FALSE, reduction.model="wnn.umap",
                                      refdata = list(
                                        celltype.l1 = "celltype.l1",
                                        celltype.l2 = "celltype.l2",
                                        predicted_ADT = "ADT"
                                      )){
  # find transfer anchors between the reference and the query, the query is your dataset
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = normalization.method,
    reference.reduction = reference.reduction,
    dims = dims,
    recompute.residuals = recompute.residuals
  )
  
  # map the query to the reference, this will project your dataset onto the reference, and copy the labels with a certain amount of certainty (the score)
  query <- MapQuery(
    anchorset = anchors,
    query = query,
    reference = reference,
    refdata = refdata,
    reference.reduction = reference.reduction,
    reduction.model = reduction.model
  )
  
  return(query)
}

####################
# Main Code        #
####################

# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds', sep = '')
object_loc_v2_azt <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029_azimuth.rds', sep = '')
object_loc_v3_azt <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106_azimuth.rds', sep = '')
reference_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/seurat_multimodal/references/pbmc_multimodal.h5seurat'


# load the reference Seurat object
reference <- LoadH5Seurat(reference_loc)
# load the v2 object
v2 <- readRDS(object_loc_v2)
# match assay to reference
DefaultAssay(v2) <- 'SCT'
# do the cell typing
v2 <-add_cell_type_predictions(reference = reference, query = v2)
# save the result
saveRDS(v2, object_loc_v2_azt)

# and the whole shebang for v2 as well
v3 <- readRDS(object_loc_v3)
DefaultAssay(v3) <- 'SCT'
v3 <-add_cell_type_predictions(reference = reference, query = v3)
saveRDS(v3, object_loc_v3_azt)

