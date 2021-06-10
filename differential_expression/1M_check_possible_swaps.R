library(Seurat)
require("heatmap.plus")
library(RColorBrewer)
library(ggpubr)

do_MAST_per_participant <- function(seurat_object, condition.1, condition.2, output_loc, cell_types_to_use=NULL, split.column = 'timepoint', cell.type.column = 'cell_type', assignment.column='assignment', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.1){
  cell_types <- cell_types_to_use
  # if the cell types were not supplied, 
  if(is.null(cell_types_to_use)){
    cell_types <- unique(seurat_object@meta.data[[cell.type.column]])
  }
  # do the MAST per cell type
  for(cell_type in cell_types){
    # subset to the cell type
    seurat_object_cell_type <- seurat_object[,seurat_object@meta.data[cell.type.column] == cell_type]
    # check for each participant
    for(participant in unique(seurat_object_cell_type@meta.data[[assignment.column]])){
      try({
        Idents(seurat_object_cell_type) <- seurat_object_cell_type@meta.data$timepoint
        # subset to the participant
        seurat_object_cell_type_part <- seurat_object_cell_type[,seurat_object_cell_type@meta.data[assignment.column] == participant]
        # do the actual DE analysis
        result <- FindMarkers(object = seurat_object_cell_type_part, ident.1 = condition.1, ident.2 = condition.2, test.use = 'MAST', min.pct = min.pct, logfc.threshold = logfc.threshold, assay = assay)
        # set up the output location
        output_loc_full <- paste(output_loc, cell_type, '_', condition.1, condition.2, '_', participant, '.tsv', sep = '')
        # write the result
        write.table(result, output_loc_full, sep = '\t', row.names = T, col.names = T)
      })
    }
  }
}

do_MAST <- function(seurat_object, condition.1, condition.2, output_loc, cell_types_to_use=NULL, split.column = 'timepoint', cell.type.column = 'cell_type', assignment.column='assignment', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.1){
  cell_types <- cell_types_to_use
  # if the cell types were not supplied, 
  if(is.null(cell_types_to_use)){
    cell_types <- unique(seurat_object@meta.data[[cell.type.column]])
  }
  # do the MAST per cell type
  for(cell_type in cell_types){
    Idents(seurat_object_cell_type) <- seurat_object_cell_type@meta.data$timepoint
    # subset to the cell type
    seurat_object_cell_type <- seurat_object[,seurat_object@meta.data[cell.type.column] == cell_type]
    Idents(seurat_object_cell_type) <- seurat_object_cell_type@meta.data$timepoint
    # do the actual DE analysis
    result <- FindMarkers(object = seurat_object_cell_type, ident.1 = condition.1, ident.2 = condition.2, test.use = 'MAST', min.pct = min.pct, logfc.threshold = logfc.threshold, assay = assay)
    # set up the output location
    output_loc_full <- paste(output_loc, cell_type, '_', condition.1, condition.2, '.tsv', sep = '')
    # write the result
    write.table(result, output_loc_full, sep = '\t', row.names = T, col.names = T)
  }
}

combine_MAST_outputs <- function(mast_output_loc, merged_output_loc, condition.1 = 'UT', condition.2 = 'X3hPA', cell_types_to_check=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
  for(cell_type in cell_types_to_check){
    # create regex to list the files
    list_dir_regex <- paste(cell_type, '_', condition.1, condition.2, '_', '.*\\.tsv', sep = '')
    # list the files
    files <- list.files(mast_output_loc, list_dir_regex)
    # init DE result
    de_table_ct <- NULL
    # check each file
    for(i in 1:length(files)){
      # get the file
      file_loc <- paste(mast_output_loc, files[[i]], sep = '')
      # read the table
      de_results <- read.table(file_loc, sep = '\t', header = T, row.names = 1)
      # extract the participant id, which is after the other parameters such as ct and cond, but before '.tsv'
      participant <- substring(files[[i]], nchar(paste(cell_type, '_', condition.1, condition.2, sep = ''))+2, nchar(files[[i]])-4)
      # add prepend
      colnames(de_results) <- paste(participant, colnames(de_results), sep = '_')
      # add gene as column name
      de_results$gene <- rownames(de_results)
      # convert to datatable
      de_results <- data.table(de_results)
      # add to table or create new one if non-existant
      if(is.null(de_table_ct)){
        de_table_ct <- de_results
      }
      else{
        de_table_ct <- merge(de_table_ct, de_results, by='gene', all=T)
      }
    }
    merged_output_loc_ct <- paste(merged_output_loc, cell_type, condition.1, condition.2, '.tsv', sep = '')
    write.table(de_table_ct, merged_output_loc_ct, sep = '\t', col.names = T, row.names = F)
  }
}

plot_lfc_per_participant_hm <- function(combined_de_table){
  lfc_table <- NULL
  # grab the logfold column names
  logfc_cols <- colnames(mono_combined_table)[grep('_avg_logFC', colnames(mono_combined_table))]
  # replace the '_avg_logFC' portion with nothing so that we end up with just the participants
  participants <- gsub('_avg_logFC', '', logfc_cols)
  # now do each participant
  for(participant in participants){
    # grab the logfc column
    logfcs <- combined_de_table[[paste(participant, 'avg_logFC', sep = '_')]]
    # grab the pvals
    pvals <- combined_de_table[[paste(participant, 'p_val_adj', sep = '_')]]
    # make all the logfcs that are insignificant or NA zero
    logfcs[is.na(pvals) | pvals > 0.05 ] <- 0
    # create table
    part_table <- data.frame(logfcs)
    colnames(part_table) <- participant
    rownames(part_table) <- rownames(combined_de_table)
    # add to dataframe
    if(is.null(lfc_table)){
      lfc_table <- part_table
    }
    else{
      lfc_table <- cbind(lfc_table, part_table)
    }
  }
  # remove any gene that is never DE
  lfc_table <- lfc_table[rowSums(lfc_table) != 0, ]
  heatmap.3(t(lfc_table), labCol = NA, col=(brewer.pal(10,"RdBu")))
}

get_most_varying_from_df <- function(dataframe, top_so_many=10){
  # now calculate the sd over this set of columns
  sds <- apply(dataframe, 1, sd, na.rm=T)
  # add the sds as a column
  dataframe$sds <- sds
  # order by the sd
  dataframe <- dataframe[order(dataframe$sds, decreasing = T), ]
  # we will return the rownames
  most_varied <- NULL
  # we need to make sure we can return as many rownames as requested
  if(nrow(dataframe) < top_so_many){
    print(paste('requested ', top_so_many, ', but dataframe only has ', nrow(most_varied), ' rows', sep = ''))
    most_varied <- rownames(dataframe)
  }
  else{
    most_varied <- rownames(dataframe)[1:top_so_many]
  }
  return(most_varied)
}

plot_per_participant_concordance <- function(){
  
}


# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')
# read the V3 object
v3 <- readRDS(object_loc_v3)
# we want to look specifically at UT vs 3hPA
v3_ut_3hpa <- subset(v3, subset = timepoint == 'UT' | timepoint == 'X3hPA')
# setup the location of the per-participant DE for UT vs 3hPA
per_participant_de_v3_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/differential_expression/output/MAST/per_participant_de_paired_lores_lfc01minpct01_20200617/v3_per_participant_de_paired_lores_lfc01minpct01_20200617/'
# do the DE
do_MAST_per_participant(v3_ut_3hpa, condition.1='UT', condition.2='X3hPA', output_loc=per_participant_de_v3_output_loc, cell_types_to_use=c('monocyte'), split.column = 'timepoint', cell.type.column = 'cell_type_lowerres', assignment.column='assignment', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.1)
# two lane DE of the lanes that seem to have issues
v3_lane_specific <- subset(v3, subset = lane == '190123_lane1' | lane == '190123_lane2')
# this lane specific output loc
lane_specific_de_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/differential_expression/output/MAST/190123_lane1_190123_lane2_specific_de_paired_lores_lfc01minpct01_20200617/'
# do the DE
do_MAST(seurat_object=v3_lane_specific, condition.1='UT', condition.2='X3hPA', output_loc=lane_specific_de_output_loc, cell_types_to_use=c('monocyte'), split.column = 'timepoint', cell.type.column = 'cell_type_lowerres', assignment.column='assignment', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.1)

# redo the clustering and dimension reduction to check
DefaultAssay(v3_ut_3hpa) <- 'SCT'
v3_ut_3hpa <- RunPCA(v3_ut_3hpa, verbose = FALSE)
v3_ut_3hpa <- RunUMAP(v3_ut_3hpa, dims = 1:30, verbose = FALSE)
v3_ut_3hpa <- FindNeighbors(v3_ut_3hpa, dims = 1:30, verbose = FALSE)
v3_ut_3hpa <- FindClusters(v3_ut_3hpa, verbose = FALSE)
# save the object
v3_ut_3hpa_output_loc <- paste(object_loc, '1M_v3_UTX3hPA_mediumQC_ctd_rnanormed_demuxids_20201027.rds', sep = '')
saveRDS(v3_ut_3hpa, v3_ut_3hpa_output_loc)
# create some plots
DimPlot(v3_ut_3hpa, reduction = 'umap', group.by = 'seurat_clusters')
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_clus.png', width = 10, height = 10)
DimPlot(v3_ut_3hpa, reduction = 'umap', group.by = 'cell_type_lowerres')
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_ct.png', width = 10, height = 10)
DimPlot(v3_ut_3hpa, reduction = 'umap', group.by = 'timepoint')
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_tp.png', width = 10, height = 10)
DimPlot(v3_ut_3hpa, reduction = 'umap', group.by = 'lane')
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_lane.png', width = 10, height = 10)
DimPlot(v3_ut_3hpa, reduction = 'umap', group.by = 'participant')
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_part.png', width = 10, height = 10)
# add a new column to indicate the specific samples that are possibly swapped
v3_ut_3hpa@meta.data$participant_timepoint <- paste(v3_ut_3hpa@meta.data$assignment, v3_ut_3hpa@meta.data$timepoint, sep = '_')
# plot only the suspected participants and their timepoint
v3_ut_3hpa@meta.data$participant_timepoint_specific <- v3_ut_3hpa@meta.data$participant_timepoint
v3_ut_3hpa@meta.data[!(v3_ut_3hpa@meta.data$assignment %in% c('LLDeep_1058','LLDeep_1229', 'LLDeep_1179', 'LLDeep_1247', 'LLDeep_1016', 'LLDeep_1067', 'LLDeep_0747', 'LLDeep_0906')), ]$participant_timepoint_specific <- NA
DimPlot(v3_ut_3hpa[, !is.na(v3_ut_3hpa@meta.data$participant_timepoint_specific)], reduction='umap', group.by='participant_timepoint_specific')
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_parttp.png', width = 10, height = 10)
# plot only the non-suspected participants I think are not swapped and their timepoints
v3_ut_3hpa@meta.data$participant_timepoint_specific <- v3_ut_3hpa@meta.data$participant_timepoint
v3_ut_3hpa@meta.data[!(v3_ut_3hpa@meta.data$assignment %in% c('LLDeep_1133','LLDeep_1318','LLDeep_0853','LLDeep_1191', 'LLDeep_0768','LLDeep_0717','LLDeep_1045','LLDeep_0526')), ]$participant_timepoint_specific <- NA
DimPlot(v3_ut_3hpa[, !is.na(v3_ut_3hpa@meta.data$participant_timepoint_specific)], reduction='umap', group.by='participant_timepoint_specific')
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_parttp_unsuspect.png', width = 10, height = 10)
# plot only the suspected participants, but only the UT timepoint
v3_ut_3hpa@meta.data$participant_timepoint_specific <- v3_ut_3hpa@meta.data$participant_timepoint
v3_ut_3hpa@meta.data[!(v3_ut_3hpa@meta.data$assignment %in% c('LLDeep_1058','LLDeep_1229', 'LLDeep_1179', 'LLDeep_1247', 'LLDeep_1016', 'LLDeep_1067', 'LLDeep_0747', 'LLDeep_0906') & v3_ut_3hpa@meta.data$timepoint == 'UT'), ]$participant_timepoint_specific <- NA
DimPlot(v3_ut_3hpa[, !is.na(v3_ut_3hpa@meta.data$participant_timepoint_specific)], reduction='umap', group.by='participant_timepoint_specific', order=c('LLDeep_1058_UT','LLDeep_1229_UT', 'LLDeep_1179_UT', 'LLDeep_1247_UT', 'LLDeep_1016_UT', 'LLDeep_1067_UT', 'LLDeep_0747_UT', 'LLDeep_0906_UT'))
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_parttp_ut.png', width = 10, height = 10)
# plot only the suspected participants, but only the X3hPA timepoint
v3_ut_3hpa@meta.data$participant_timepoint_specific <- v3_ut_3hpa@meta.data$participant_timepoint
v3_ut_3hpa@meta.data[!(v3_ut_3hpa@meta.data$assignment %in% c('LLDeep_1058','LLDeep_1229', 'LLDeep_1179', 'LLDeep_1247', 'LLDeep_1016', 'LLDeep_1067', 'LLDeep_0747', 'LLDeep_0906') & v3_ut_3hpa@meta.data$timepoint == 'X3hPA'), ]$participant_timepoint_specific <- NA
DimPlot(v3_ut_3hpa[, !is.na(v3_ut_3hpa@meta.data$participant_timepoint_specific)], reduction='umap', group.by='participant_timepoint_specific', order=c('LLDeep_1058_X3hPA','LLDeep_1229_X3hPA', 'LLDeep_1179_X3hPA', 'LLDeep_1247_X3hPA', 'LLDeep_1016_X3hPA', 'LLDeep_1067_X3hPA', 'LLDeep_0747_X3hPA', 'LLDeep_0906_X3hPA'))
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_parttp_x3hpa.png', width = 10, height = 10)
# plot only the non-suspected participants I think are not swapped, but only the UT timepoint
v3_ut_3hpa@meta.data$participant_timepoint_specific <- v3_ut_3hpa@meta.data$participant_timepoint
v3_ut_3hpa@meta.data[!(v3_ut_3hpa@meta.data$assignment %in% c('LLDeep_1133','LLDeep_1318','LLDeep_0853','LLDeep_1191', 'LLDeep_0768','LLDeep_0717','LLDeep_1045','LLDeep_0526')  & v3_ut_3hpa@meta.data$timepoint == 'UT'), ]$participant_timepoint_specific <- NA
DimPlot(v3_ut_3hpa[, !is.na(v3_ut_3hpa@meta.data$participant_timepoint_specific)], reduction='umap', group.by='participant_timepoint_specific')
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_parttp_unsuspect_ut.png', width = 10, height = 10)
# plot only the non-suspected participants I think are not swapped, but only the 3hPA timepoint
v3_ut_3hpa@meta.data$participant_timepoint_specific <- v3_ut_3hpa@meta.data$participant_timepoint
v3_ut_3hpa@meta.data[!(v3_ut_3hpa@meta.data$assignment %in% c('LLDeep_1133','LLDeep_1318','LLDeep_0853','LLDeep_1191', 'LLDeep_0768','LLDeep_0717','LLDeep_1045','LLDeep_0526')  & v3_ut_3hpa@meta.data$timepoint == 'X3hPA'), ]$participant_timepoint_specific <- NA
DimPlot(v3_ut_3hpa[, !is.na(v3_ut_3hpa@meta.data$participant_timepoint_specific)], reduction='umap', group.by='participant_timepoint_specific')
ggsave('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v3_ut_3hpa_parttp_unsuspect_x3hpa.png', width = 10, height = 10)

# grab the v2 3hCA unintegrated object I made
v2_ut_3hca <- readRDS(paste(object_loc, '1M_v2_UTX3hCA_mediumQC_ctd_rnanormed_demuxids_20201027.rds', sep = ''))
# grab only the UT
v2_ut_3hca_utonly <- subset(v2_ut_3hca, subset = timepoint == 'UT')
# do this per lane
for(i in 1:9){
  from <- (i*8)-7
  to <- i*8
  v2_ut_3hca_utonly_i <- DimPlot(v2_ut_3hca_utonly[, v2_ut_3hca_utonly@meta.data$assignment %in% unique(v2_ut_3hca_utonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_3hca_utonly_',i,'.png', sep=''), width=10, height=10)
}
# grab the v2 3hCA unintegrated object I made
v2_ut_3hpa <- readRDS(paste(object_loc, '1M_v2_UTX3hPA_mediumQC_ctd_rnanormed_demuxids_20201027.rds', sep = ''))
# grab only the UT
v2_ut_3hpa_utonly <- subset(v2_ut_3hpa, subset = timepoint == 'UT')
# do this per lane
for(i in 1:9){
  from <- (i*8)-7
  to <- i*8
  v2_ut_3hpa_utonly_i <- DimPlot(v2_ut_3hpa_utonly[, v2_ut_3hpa_utonly@meta.data$assignment %in% unique(v2_ut_3hpa_utonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_3hpa_utonly_',i,'.png', sep=''), width=10, height=10)
}
# grab the v2 3hCA unintegrated object I made
v2_ut_3hmtb <- readRDS(paste(object_loc, '1M_v2_UTX3hMTB_mediumQC_ctd_rnanormed_demuxids_20201027.rds', sep = ''))
# grab only the UT
v2_ut_3hmtb_utonly <- subset(v2_ut_3hmtb, subset = timepoint == 'UT')
# do this per lane
for(i in 1:9){
  from <- (i*8)-7
  to <- i*8
  v2_ut_3hmtb_utonly_i <- DimPlot(v2_ut_3hpa_utonly[, v2_ut_3hmtb_utonly@meta.data$assignment %in% unique(v2_ut_3hmtb_utonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_3hmtb_utonly_',i,'.png', sep=''), width=10, height=10)
}
# grab the v2 24hCA unintegrated object I made
v2_ut_24hca <- readRDS(paste(object_loc, '1M_v2_UTX24hCA_mediumQC_ctd_rnanormed_demuxids_20201027.rds', sep = ''))
# grab only the UT
v2_ut_24hca_utonly <- subset(v2_ut_24hca, subset = timepoint == 'UT')
# do this per lane
for(i in 1:9){
  from <- (i*8)-7
  to <- i*8
  v2_ut_24hca_utonly_i <- DimPlot(v2_ut_24hca_utonly[, v2_ut_24hca_utonly@meta.data$assignment %in% unique(v2_ut_24hca_utonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_24hca_utonly_',i,'.png', sep=''), width=10, height=10)
}
# grab the v2 24hCA unintegrated object I made
v2_ut_24hpa <- readRDS(paste(object_loc, '1M_v2_UTX24hPA_mediumQC_ctd_rnanormed_demuxids_20201027.rds', sep = ''))
# grab only the UT
v2_ut_24hpa_utonly <- subset(v2_ut_24hpa, subset = timepoint == 'UT')
# do this per lane
for(i in 1:9){
  from <- (i*8)-7
  to <- i*8
  v2_ut_24hpa_utonly_i <- DimPlot(v2_ut_24hpa_utonly[, v2_ut_24hpa_utonly@meta.data$assignment %in% unique(v2_ut_24hpa_utonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_24hpa_utonly_',i,'.png', sep=''), width=10, height=10)
}
# grab the v2 24hCA unintegrated object I made
v2_ut_24hmtb <- readRDS(paste(object_loc, '1M_v2_UTX24hMTB_mediumQC_ctd_rnanormed_demuxids_20201027.rds', sep = ''))
# grab only the UT
v2_ut_24hmtb_utonly <- subset(v2_ut_24hmtb, subset = timepoint == 'UT')
# do this per lane
for(i in 1:9){
  from <- (i*8)-7
  to <- i*8
  v2_ut_24hmtb_utonly_i <- DimPlot(v2_ut_24hpa_utonly[, v2_ut_24hmtb_utonly@meta.data$assignment %in% unique(v2_ut_24hmtb_utonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_24hmtb_utonly_',i,'.png', sep=''), width=10, height=10)
}

# grab only the UT
v2_ut_3hca_3hcaonly <- subset(v2_ut_3hca, subset = timepoint == 'X3hCA')
# do this per lane
for(i in 1:10){
  from <- (i*8)-7
  to <- i*8
  v2_ut_3hca_3hca_i <- DimPlot(v2_ut_3hca_3hcaonly[, v2_ut_3hca_3hcaonly@meta.data$assignment %in% unique(v2_ut_3hca_3hcaonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_3hca_3hca_',i,'.png', sep=''), width=10, height=10)
}
# grab only the UT
v2_ut_3hpa_3hpaonly <- subset(v2_ut_3hpa, subset = timepoint == 'X3hPA')
# do this per lane
for(i in 1:11){
  from <- (i*8)-7
  to <- i*8
  v2_ut_3hpa_3hpa_i <- DimPlot(v2_ut_3hpa_3hpaonly[, v2_ut_3hpa_3hpaonly@meta.data$assignment %in% unique(v2_ut_3hpa_3hpaonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_3hpa_3hpaonly_',i,'.png', sep=''), width=10, height=10)
}
# grab only the UT
v2_ut_3hmtb_3hmtbonly <- subset(v2_ut_3hmtb, subset = timepoint == 'X3hMTB')
# do this per lane
for(i in 1:9){
  from <- (i*8)-7
  to <- i*8
  v2_ut_3hmtb_3hmtb_i <- DimPlot(v2_ut_3hmtb_3hmtbonly[, v2_ut_3hmtb_3hmtbonly@meta.data$assignment %in% unique(v2_ut_3hmtb_3hmtbonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_3hmtb_3hmtbonly_',i,'.png', sep=''), width=10, height=10)
}

# grab only the UT
v2_ut_24hca_24hcaonly <- subset(v2_ut_24hca, subset = timepoint == 'X24hCA')
# do this per lane
for(i in 1:10){
  from <- (i*8)-7
  to <- i*8
  v2_ut_24hca_24hca_i <- DimPlot(v2_ut_24hca_24hcaonly[, v2_ut_24hca_24hcaonly@meta.data$assignment %in% unique(v2_ut_24hca_24hcaonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_24hca_24hca_',i,'.png', sep=''), width=10, height=10)
}
# grab only the UT
v2_ut_24hpa_24hpaonly <- subset(v2_ut_24hpa, subset = timepoint == 'X24hPA')
# do this per lane
for(i in 1:10){
  from <- (i*8)-7
  to <- i*8
  v2_ut_24hpa_24hpa_i <- DimPlot(v2_ut_24hpa_24hpaonly[, v2_ut_24hpa_24hpaonly@meta.data$assignment %in% unique(v2_ut_24hpa_24hpaonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_24hpa_24hpaonly_',i,'.png', sep=''), width=10, height=10)
}
# grab only the UT
v2_ut_24hmtb_24hmtbonly <- subset(v2_ut_24hmtb, subset = timepoint == 'X24hMTB')
# do this per lane
for(i in 1:10){
  from <- (i*8)-7
  to <- i*8
  v2_ut_24hmtb_24hmtb_i <- DimPlot(v2_ut_24hmtb_24hmtbonly[, v2_ut_24hmtb_24hmtbonly@meta.data$assignment %in% unique(v2_ut_24hmtb_24hmtbonly@meta.data$assignment)[from:to]], reduction = 'umap', group.by = 'assignment')
  ggsave(paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/swap_check/clustering/plots/v2_ut_24hmtb_24hmtbonly_',i,'.png', sep=''), width=10, height=10)
}


# get the local output folder of the per participant output
lane_specific_de_output_loc <- '/data/scRNA/swap_check/differential_expression/output/MAST/per_participant_de_paired_lores_lfc01minpct01_20200617/v3_per_participant_de_paired_lores_lfc01minpct01_20200617/'
# and where to output the merged result
lane_merged_de_output_loc <- '/data/scRNA/swap_check/differential_expression/output/MAST/per_participant_de_paired_lores_lfc01minpct01_20200617_merged/v3_per_participant_de_paired_lores_lfc01minpct01_20200617_merged/'
# write the merged table
combine_MAST_outputs(lane_specific_de_output_loc, lane_merged_de_output_loc, cell_types_to_check = c('monocyte'))
# get the mono combined table
mono_combined_table_loc <- '/data/scRNA/swap_check/differential_expression/output/MAST/per_participant_de_paired_lores_lfc01minpct01_20200617_merged/v3_per_participant_de_paired_lores_lfc01minpct01_20200617_merged/monocyteUTX3hPA.tsv'
mono_combined_table <- read.table(mono_combined_table_loc, sep = '\t', header = T, row.names = 1)
plot_lfc_per_participant_hm(mono_combined_table)
# try with a subset of the genes, the image stays the same
mono_combined_table_most_varied_genes <- get_most_varying_from_df(mono_combined_table, top_so_many = 250)
mono_combined_table_most_varied <- mono_combined_table[mono_combined_table_most_varied_genes, ]
plot_lfc_per_participant_hm(mono_combined_table_most_varied)
# try with specifically the interferon pathways
pathway_gene_loc <- '/data/scRNA/pathways/'
interferon_signalling_loc <- paste(pathway_gene_loc, 'REACTOME_Interferon_Signaling_genes.txt', sep = '')
interferon_df <- read.table(interferon_signalling_loc)
interferon_genes <- interferon_df$V1
mono_combined_table_interferon <- mono_combined_table[rownames(mono_combined_table) %in% interferon_genes, ]
plot_lfc_per_participant_hm(mono_combined_table_interferon)
# get the DE for the entirety of UTX3hPA
mono_ut_3hpa_full <- read.table('/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v3_paired_lores_lfc01minpct01_20200713/rna/monocyteUTX3hPA.tsv', sep = '\t', header = T, row.names = 1)
# set insignificant logfcs to zero
mono_ut_3hpa_full[mono_ut_3hpa_full$p_val_adj > 0.05, ]$avg_logFC <- 0
# get the genes in both the full and the per-participant logfc tables
genes_full_and_perpart <- intersect(rownames(mono_ut_3hpa_full), rownames(mono_combined_table))
# 


