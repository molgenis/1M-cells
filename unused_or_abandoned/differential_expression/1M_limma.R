####################
# libraries        #
####################

library(Seurat)
library(limma)

####################
# Functions        #
####################

perform_limma_trend <- function(seurat_object, output_loc, split.column = 'timepoint', assay = 'RNA', lfc = 0.25){
  # set output loc
  output_loc_final <- paste(output_loc, '.tsv', sep = '')
  print(paste('trying analysis for', output_loc_final))
  # set the assay we wish to analyse
  DefaultAssay(seurat_object) <- assay
  # a hard column name is required for the model.matrix step, so let's add that
  seurat_object <- AddMetaData(seurat_object, seurat_object@meta.data[split.column], 'condition')
  # get expression from object
  expr <- as.matrix(GetAssayData(seurat_object))
  # get zero-expression genes
  bad <- which(rowSums(expr) == 0)
  # remove zero-expression genes
  expr <- expr[-bad,]
  # create model
  mm <- model.matrix(~0 + condition, data = seurat_object@meta.data)
  # create fit
  fit <- lmFit(expr, mm)
  # get DE with limma trend (lfc is significant logfold change)
  #fit <- treat(fit, lfc=log2(lfc), trend = T)
  fit <- treat(fit, lfc=lfc, trend = T)
  # get all the results
  result <- topTreat(fit, n = nrow(expr), sort.by = 'P')
  # save the result
  write.table(result, output_loc_final, sep = '\t')
  # remove largest object
  rm(expr)
}

perform_limma_trend_per_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', lfc = 0.25){
  # go through the cell types
  for(cell_type in unique(seurat_object@meta.data[[cell.type.column]])){
    # make a more specific path
    output_loc_cell_type <- paste(output_loc, cell_type, sep = '')
    # grab the subset
    seurat_object_cell_type <- seurat_object[,seurat_object@meta.data[cell.type.column] == cell_type]
    # do the limma
    try({
      perform_limma_trend(seurat_object_cell_type, output_loc_cell_type, split.column = split.column, assay = assay, lfc = lfc)
    })
      # clear memory
    rm(seurat_object_cell_type)
  }
  # finally do a bulk analysis as well
  output_loc_bulk <- paste(output_loc, 'bulk', sep = '')
  try({
    perform_limma_trend(seurat_object, output_loc_bulk, split.column = split.column, assay = assay, lfc = lfc)
  })
}

perform_limma_trend_per_set_and_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', lfc = 0.25){
  # grab CA combo and do limma
  seurat_object_ca_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'UT' | seurat_object@meta.data[split.column] == 'X3hCA' | seurat_object@meta.data[split.column] == 'X24hCA']
  output_loc_ca <- paste(output_loc, 'CA', sep = '')
  perform_limma_trend_per_celltype(seurat_object_ca_compare, output_loc_ca, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_ca_compare)
  # grab PA combo and do limma
  seurat_object_pa_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'UT' | seurat_object@meta.data[split.column] == 'X3hPA' | seurat_object@meta.data[split.column] == 'X24hPA']
  output_loc_pa <- paste(output_loc, 'PA', sep = '')
  perform_limma_trend_per_celltype(seurat_object_pa_compare, output_loc_pa, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_pa_compare)
  # grab MBT combo and do limma
  seurat_object_mtb_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'UT' | seurat_object@meta.data[split.column] == 'X3hMTB' | seurat_object@meta.data[split.column] == 'X24hMTB']
  output_loc_mtb <- paste(output_loc, 'MTB', sep = '')
  perform_limma_trend_per_celltype(seurat_object_mtb_compare, output_loc_mtb, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_mtb_compare)
}
  
perform_limma_trend_per_set_and_celltype_paired <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', lfc = 1.2){
  # grab CA combo and do limma
  seurat_object_ca_ut_3h_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'UT' | seurat_object@meta.data[split.column] == 'X3hCA']
  output_loc_ca <- paste(output_loc, 'CA_UT_3hCA', sep = '')
  perform_limma_trend_per_celltype(seurat_object_ca_ut_3h_compare, output_loc_ca, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_ca_ut_3h_compare)
  # grab MBT combo and do limma
  seurat_object_ca_ut_24h_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'UT' | seurat_object@meta.data[split.column] == 'X24hCA']
  output_loc_ca <- paste(output_loc, 'CA_UT_24hCA', sep = '')
  perform_limma_trend_per_celltype(seurat_object_ca_ut_24h_compare, output_loc_ca, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_ca_ut_24h_compare)
  # grab MBT combo and do limma
  seurat_object_ca_3h_24h_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'X3hCA' | seurat_object@meta.data[split.column] == 'X24hCA']
  output_loc_ca <- paste(output_loc, 'CA_3hCA_24hCA', sep = '')
  perform_limma_trend_per_celltype(seurat_object_ca_3h_24h_compare, output_loc_ca, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_ca_3h_24h_compare)
  
  # grab PA combo and do limma
  seurat_object_pa_ut_3h_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'UT' | seurat_object@meta.data[split.column] == 'X3hPA']
  output_loc_pa <- paste(output_loc, 'PA_UT_3hPA', sep = '')
  perform_limma_trend_per_celltype(seurat_object_pa_ut_3h_compare, output_loc_pa, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_pa_ut_3h_compare)
  # grab PA combo and do limma
  seurat_object_pa_ut_24h_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'UT' | seurat_object@meta.data[split.column] == 'X24hPA']
  output_loc_pa <- paste(output_loc, 'PA_UT_24hPA', sep = '')
  perform_limma_trend_per_celltype(seurat_object_pa_ut_24h_compare, output_loc_pa, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_pa_ut_24h_compare)
  # grab PA combo and do limma
  seurat_object_pa_3h_24h_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'X3hPA' | seurat_object@meta.data[split.column] == 'X24hPA']
  output_loc_pa <- paste(output_loc, 'PA_3hPA_24hPA', sep = '')
  perform_limma_trend_per_celltype(seurat_object_pa_3h_24h_compare, output_loc_pa, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_pa_3h_24h_compare)
  
  # grab MTB combo and do limma
  seurat_object_mtb_ut_3h_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'UT' | seurat_object@meta.data[split.column] == 'X3hMTB']
  output_loc_mtb <- paste(output_loc, 'MTB_UT_3hMTB', sep = '')
  perform_limma_trend_per_celltype(seurat_object_mtb_ut_3h_compare, output_loc_mtb, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_mtb_ut_3h_compare)
  # grab MTB combo and do limma
  seurat_object_mtb_ut_24h_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'UT' | seurat_object@meta.data[split.column] == 'X24hMTB']
  output_loc_mtb <- paste(output_loc, 'MTB_UT_24hMTB', sep = '')
  perform_limma_trend_per_celltype(seurat_object_mtb_ut_24h_compare, output_loc_mtb, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_mtb_ut_24h_compare)
  # grab MTB combo and do limma
  seurat_object_mtb_3h_24h_compare <- seurat_object[,seurat_object@meta.data[split.column] == 'X3hMTB' | seurat_object@meta.data[split.column] == 'X24hMTB']
  output_loc_mtb <- paste(output_loc, 'MTB_3hMTB_24hMTB', sep = '')
  perform_limma_trend_per_celltype(seurat_object_mtb_3h_24h_compare, output_loc_mtb, split.column = split.column, cell.type.column = cell.type.column, assay = assay, lfc = lfc)
  rm(seurat_object_mtb_3h_24h_compare)
}

####################
# Main Code        #
####################

# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20200427.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200427.rds', sep = '')
#object_loc <- '/data/p287578/1M_cells_scRNAseq/scanpy_preprocess_samples/objects/'
#object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20200427_anon.rds', sep = '')
#object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200427_anon.rds', sep = '')


# DE output locations
limma_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/limma/results/'
#limma_output_loc <- '/data/p287578/1M_cells_scRNAseq/differential_expression/limma/results/'
limma_output_loc_v2 <- paste(limma_output_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20200427/', sep = '')
limma_output_loc_v3 <- paste(limma_output_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200427/', sep = '')
limma_output_loc_v2_rna <- paste(limma_output_loc_v2, 'rna/', sep = '')
limma_output_loc_v3_rna <- paste(limma_output_loc_v3, 'rna/', sep = '')
limma_output_loc_v2_sct <- paste(limma_output_loc_v2, 'sct/', sep = '')
limma_output_loc_v3_sct <- paste(limma_output_loc_v3, 'sct/', sep = '')

# put in the work for v2
v2 <- readRDS(object_loc_v2)
perform_limma_trend_per_set_and_celltype(seurat_object = v2, output_loc = limma_output_loc_v2_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_limma_trend_per_set_and_celltype(seurat_object = v2, output_loc = limma_output_loc_v2_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA')

# paired as well
perform_limma_trend_per_set_and_celltype_paired(seurat_object = v2, output_loc = limma_output_loc_v2_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_limma_trend_per_set_and_celltype_paired(seurat_object = v2, output_loc = limma_output_loc_v2_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA')


# clear up memory
rm(v2)

# put in the work for v3
v3 <- readRDS(object_loc_v3)
perform_limma_trend_per_set_and_celltype(seurat_object = v3, output_loc = limma_output_loc_v3_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_limma_trend_per_set_and_celltype(seurat_object = v3, output_loc = limma_output_loc_v3_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA')

# paired as well
perform_limma_trend_per_set_and_celltype_paired(seurat_object = v3, output_loc = limma_output_loc_v3_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_limma_trend_per_set_and_celltype_paired(seurat_object = v3, output_loc = limma_output_loc_v3_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA')

# clear up memory
rm(v3)
