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

convert_cell_labels <- function(seurat_object, high_res_pred_col='predicted.celltype.l2', high_res_col='cell_type_pred', low_res_col='cell_type_low_pred'){
  # convert to character, to make easier to use
  seurat_object@meta.data[[high_res_col]] <- as.character(seurat_object@meta.data[[high_res_pred_col]])
  # harmonize the cell types again
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] == 'NK', ][[high_res_col]] <- 'NKdim'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] == 'NK_CD56bright', ][[high_res_col]] <- 'NKbright'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] == 'CD14 Mono', ][[high_res_col]] <- 'cMono'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] == 'CD16 Mono', ][[high_res_col]] <- 'ncMono'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] == 'Plasmablast', ][[high_res_col]] <- 'plasmablast'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] == 'Platelet', ][[high_res_col]] <- 'platelet'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] == 'Eryth', ][[high_res_col]] <- 'eryth'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] == 'Doublet', ][[high_res_col]] <- 'doublet'
  # spaces in variables is inconvenient in a lot of places, so we'll replace these with underscores
  seurat_object@meta.data[[high_res_col]] <- gsub(' ', '_', seurat_object@meta.data[[high_res_col]])
  # determine the higher resolution cell types
  # we will define a lower resolution cell type as well, we need to create some groupings for this
  cd4t <- c('Treg', 'CD4_Naive', 'CD4_TCM', 'CD4_TEM', 'CD4_CTL', 'CD4_Proliferating')
  cd8t <- c('MAIT', 'CD8_Naive', 'CD8_TCM', 'CD8_TEM', 'CD8_Proliferating')
  t_other <- c('dnT', 'gdT', 'ILC')
  nk <- c('NKdim', 'NKbright', 'NK_Proliferating')
  monocyte <- c('cMono', 'ncMono')
  dc <- c('cDC1', 'cDC2', 'pDC', 'ASDC')
  b <- c('B_naive', 'B_intermediate', 'B_memory')
  # add the new column by copying the higher res first, 
  seurat_object@meta.data[[low_res_col]] <- seurat_object@meta.data[[high_res_col]]
  # in this new column, overwrite them to have the lower resolution
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] %in% cd4t, ][[low_res_col]] <- 'CD4T'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] %in% cd8t, ][[low_res_col]] <- 'CD8T'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] %in% t_other, ][[low_res_col]] <- 'T_other'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] %in% nk, ][[low_res_col]] <- 'NK'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] %in% monocyte, ][[low_res_col]] <- 'monocyte'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] %in% dc, ][[low_res_col]] <- 'DC'
  seurat_object@meta.data[seurat_object@meta.data[[high_res_col]] %in% b, ][[low_res_col]] <- 'B'
  # change to factors because that is cleaner
  seurat_object@meta.data[[low_res_col]] <- as.factor(seurat_object@meta.data[[low_res_col]])
  seurat_object@meta.data[[high_res_col]] <- as.factor(seurat_object@meta.data[[high_res_col]])
  return(seurat_object)
}


create_confusion_matrix <- function(assignment_table, truth_column, prediction_column, truth_column_label=NULL, prediction_column_label=NULL){
  # convert both of the columns to characters, to make everything easier
  assignment_table[[truth_column]] <- as.character(assignment_table[[truth_column]])
  assignment_table[[prediction_column]] <- as.character(assignment_table[[prediction_column]])
  # init the table
  confusion_table <- NULL
  # check each truth
  for(truth in unique(assignment_table[[truth_column]])){
    # get these truths
    truth_rows <- assignment_table[assignment_table[[truth_column]] == truth, ]
    # check now many have this truth
    this_truth_number <- nrow(truth_rows)
    # check what was predicted for these truths
    #for(prediction in unique(truth_rows[[prediction_column]])){
    for(prediction in unique(assignment_table[[prediction_column]])){
      # check the number of this prediction
      this_prediction_number <- nrow(truth_rows[truth_rows[[prediction_column]] == prediction, ])
      # init variable
      fraction <- NULL
      # we can only calculate a fraction if the result is not zero
      if(this_prediction_number > 0){
        # calculate the fraction
        fraction <- this_prediction_number / this_truth_number
      }
      # otherwise we just set it to zero
      else{
        fraction <- 0
      }
      # turn into row
      this_row <- data.frame(truth=c(truth), prediction=c(prediction), freq=c(fraction), stringsAsFactors = F)
      # add this entry to the dataframe
      if(is.null(confusion_table)){
        confusion_table <- this_row
      }
      else{
        confusion_table <- rbind(confusion_table, this_row)
      }
    }
  }
  # round the frequency off to a sensible cutoff
  confusion_table$freq <- round(confusion_table$freq, digits=2)
  # turn into plot
  p <- ggplot(data=confusion_table, aes(x=truth, y=prediction, fill=freq)) + geom_tile() + scale_fill_gradient(low='red', high='blue') + geom_text(aes(label=freq))
  # some options
  if(!is.null(truth_column_label)){
    p <- p + xlab(truth_column_label)
  }
  if(!is.null(prediction_column_label)){
    p <- p + ylab(prediction_column_label)
  }
  return(p)
}


####################
# Main Code        #
####################

# object locations
object_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds', sep = '')
object_loc_v2_azt <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029_azimuth.rds', sep = '')
object_loc_v3_azt <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106_azimuth.rds', sep = '')
reference_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/seurat_multimodal/references/pbmc_multimodal.h5seurat'
# plot loc
plot_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/seurat_multimodal/plots/'


# load the reference Seurat object
reference <- LoadH5Seurat(reference_loc)
# load the v2 object
v2 <- readRDS(object_loc_v2)
# match assay to reference
DefaultAssay(v2) <- 'SCT'
# do the cell typing
v2 <-add_cell_type_predictions(reference = reference, query = v2)
# clean up the column names and entries
v2 <- convert_cell_labels(v2)
# save the result
saveRDS(v2, object_loc_v2_azt)

# and the whole shebang for v2 as well
v3 <- readRDS(object_loc_v3)
DefaultAssay(v3) <- 'SCT'
v3 <-add_cell_type_predictions(reference = reference, query = v3)
v3 <- convert_cell_labels(v3)
saveRDS(v3, object_loc_v3_azt)

# create a plot of marker based assignment, and azimuth based assignment
v2_confusion_matrix <- create_confusion_matrix(v2@meta.data, 'cell_type_lowerres', 'cell_type_low_pred', 'markers', 'azimuth')
ggsave(paste(plot_loc, 'v2_confusion_matrix.pdf', sep = ''), plot = v2_confusion_matrix, width = 10 ,height = 10)

# create a plot of marker based assignment, and azimuth based assignment
v3_confusion_matrix <- create_confusion_matrix(v3@meta.data, 'cell_type_lowerres', 'cell_type_low_pred', 'markers', 'azimuth')
ggsave(paste(plot_loc, 'v3_confusion_matrix.pdf', sep = ''), plot = v3_confusion_matrix, width = 10 ,height = 10)

