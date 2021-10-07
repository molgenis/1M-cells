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
library(MetaVolcanoR)
library(UpSetR)

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


write_meta_mast <- function(mast_output_loc_prepend, mast_output_loc_append, mast_meta_output_loc_prepend, cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'NK', 'hemapoietic stem', 'megakaryocyte', 'monocyte', 'plasma B'), conditions=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB'), mtc_method='bonferroni'){
  # go through the conditions
  for(condition in conditions){
    # check for each cell type
    for(cell_type in cell_types){
      # get the mast output
      mast_loc_v2 <- paste(mast_output_loc_prepend, '2', mast_output_loc_append, cell_type, 'UT', condition, '.tsv', sep = '')
      mast_loc_v3 <- paste(mast_output_loc_prepend, '3', mast_output_loc_append, cell_type, 'UT', condition, '.tsv', sep = '')
      try({
        # read the mast output
        mast_v2 <- read.table(mast_loc_v2, header=T)
        mast_v3 <- read.table(mast_loc_v3, header=T)
        # get the genes that are in both
        genes_both <- intersect(rownames(mast_v2), rownames(mast_v3))
        # there need to be genes left for this
        if(length(genes_both) == 0){
          print(paste('skipped', condition, cell_type, 'due to there not being gene overlapping in meta-analysis'))
        }
        else{
          # select only those genes
          mast_v2 <- mast_v2[rownames(mast_v2) %in% genes_both,]
          mast_v3 <- mast_v3[rownames(mast_v3) %in% genes_both,]
          # morph P val to minimum
          if(nrow(mast_v2[mast_v2$p_val == 0, ]) > 0){
            mast_v2[mast_v2$p_val == 0, ]$p_val <- .Machine$double.xmin
          }
          if(nrow(mast_v3[mast_v3$p_val == 0, ]) > 0){
            mast_v3[mast_v3$p_val == 0, ]$p_val <- .Machine$double.xmin
          }
          # add the gene name also in a column
          mast_v2$gene <- rownames(mast_v2)
          mast_v3$gene <- rownames(mast_v3)
          # add the mast results
          masts <- list()
          masts$v2 <- mast_v2
          masts$v3 <- mast_v3
          # perform the metavolcanor approach
          meta_degs_comb <- combining_mv(diffexp=masts, pcriteria='p_val', foldchangecol='avg_log2FC', genenamecol = 'gene', collaps = T)
          # grab the result we care about
          volcanometa <- meta_degs_comb@metaresult
          #volcanometa$metap_bonferroni <- volcanometa$metap*length(genes_both)
          volcanometa[[paste('metap', mtc_method, sep = '_')]] <- p.adjust(volcanometa$metap, method = mtc_method)
          # add the genes as rownames
          rownames(volcanometa) <- volcanometa$gene
          # add a colname append
          colnames(mast_v2) <- paste(colnames(mast_v2), 'v2', sep = '_')
          colnames(mast_v3) <- paste(colnames(mast_v3), 'v3', sep = '_')
          # merge the frames
          mast <- merge(mast_v2, mast_v3, by=0, all=TRUE)
          rownames(mast) <- mast$Row.names
          mast$Row.names <- NULL
          # also add the volcanometa stuff
          mast <- merge(mast, volcanometa, by=0, all=TRUE)
          rownames(mast) <- mast$Row.names
          mast$Row.names <- NULL
          #if(nrow(mast[mast$metap_bonferroni > 1, ]) > 0){
          #  mast[mast$metap_bonferroni > 1, ]$metap_bonferroni <- 1
          #}
          # write the result
          output_loc <- paste(mast_meta_output_loc_prepend, cell_type, 'UT', condition, '.tsv', sep = '')
          write.table(mast, output_loc, sep = '\t')
        }
      })
    }
  }
}


compare_overlap_low_to_high_res_de_genes <- function(low_res_output_loc, high_res_output_loc, cell_type_matches, conditions=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), sig_column='metap', sig_cutoff=0.05){
  # we will have plot for each condition
  plots_per_condition <- list()
  # check each condition
  for(condition in conditions){
    # we will also have a plot per major cell type
    plot_per_major <- list()
    # each name in the cell type matches is the major cell type
    for(major in names(cell_type_matches)){
      # save the gene lists in a list we can use for upset
      gene_lists <- list()
      # paste together the path to the major MAST output
      major_output_loc <- paste(low_res_output_loc, major, 'UTX', condition, '.tsv', sep = '')
      # read the output
      major_output <- read.table(major_output_loc, sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
      # get the genes that are significant
      sig_genes_major <- rownames(major_output[major_output[[sig_column]] < sig_cutoff, ])
      # add to the list
      gene_lists[[major]] <- sig_genes_major
      # get the high resolution cell types this major type is linked to
      minor_matches <- cell_type_matches[[major]]
      # check each of these minor cell types
      for(minor in minor_matches){
        try({
          # paste together the path to the minor MAST output
          minor_output_loc <- paste(high_res_output_loc, minor, 'UTX', condition, '.tsv', sep = '')
          # read the output
          minor_output <- read.table(minor_output_loc, sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
          # get the genes that are significant
          sig_genes_minor <- rownames(minor_output[minor_output[[sig_column]] < sig_cutoff, ])
          # add to list
          gene_lists[[minor]] <- sig_genes_minor
        })
      }
      # create the upset plot
      p <- upset(data = fromList(gene_lists), nsets = length(names(gene_lists)), order.by = 'freq')
      # add to the list
      plot_per_major[[major]] <- p
    }
    # add list of plots to condition
    plots_per_condition[[condition]] <- plot_per_major
  }
  return(plots_per_condition)
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
mast_output_paired_highres_loc_v2 <- paste(mast_output_loc, 'paired_highres_lfc01minpct01_20210905/v2_paired_highres_lfc01minpct01_20210905/', sep = '')
mast_output_paired_highres_loc_v3 <- paste(mast_output_loc, 'paired_highres_lfc01minpct01_20210905/v3_paired_highres_lfc01minpct01_20210905/', sep = '')
# we'll use the RNA assay
mast_output_paired_highres_rna_loc_v2 <- paste(mast_output_paired_highres_loc_v2, 'rna/', sep = '')
mast_output_paired_highres_rna_loc_v3 <- paste(mast_output_paired_highres_loc_v3, 'rna/', sep = '')
# additionally, we will do a separate T cell DE analysis
mast_output_paired_t_loc_v2 <- paste(mast_output_loc, 'paired_highres_lfc01minpct01_20210905/v2_paired_t_lfc01minpct01_20210905/', sep = '')
mast_output_paired_t_loc_v3 <- paste(mast_output_loc, 'paired_highres_lfc01minpct01_20210905/v3_paired_t_lfc01minpct01_20210905/', sep = '')
# and the assignments for these T cells needs to come from somewhere as well
t_classifications_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/seurat_multimodal/classifications/'
t_classifications_loc_v2 <- paste(t_classifications, 'v2_azi_to_cluster_l2_cell_types.tsv', sep = '')
t_classifications_loc_v3 <- paste(t_classifications, 'v3_azi_to_cluster_l2_cell_types.tsv', sep = '')

# we'll need some plots as well
mast_overlap_plot_loc <- '/data/scRNA/differential_expression/seurat_MAST/overlap/meta_paired_highres_lfc01minpct01_20210905_meta_paired_lores_lfc01minpct01_20201106/'

# read the object
v2 <- readRDS(object_loc_v2)
DefaultAssay(v2) <- 'RNA'
# we've done some refinements at the marker gene level, let's make those changes permanent
v2 <- v2[, !is.na(v2@meta.data$timepoint)]
v2 <- v2[, !is.na(v2@meta.data$assignment)]
v2 <- v2[, !is.na(v2@meta.data$cell_type)]
# we've done some refinements at the marker gene level, let's make those changes permanent
v2@meta.data[v2@meta.data$cell_type == 'NK', 'cell_type'] <- 'NKdim'
levels(v2@meta.data$cell_type) <- c(levels(v2@meta.data$cell_type), 'cMono', 'ncMono')
v2@meta.data[v2@meta.data$cell_type %in% c('mono 1', 'mono 4'), 'cell_type'] <- 'cMono'
v2@meta.data[v2@meta.data$cell_type %in% c('mono 2'), 'cell_type'] <- 'ncMono'
v2@meta.data$cell_type <- droplevels(v2@meta.data$cell_type)
# clean up
v2 <- v2[, !is.na(v2@meta.data$cell_type) & !is.na(v2@meta.data$assignment) & !is.na(v2@meta.data$timepoint)]
# write the new object
saveRDS(v2, object_loc_v2_new)
# do the mapping
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_highres_rna_loc_v2, cell_types_to_use = c('NKdim', 'NKbright', 'mDC', 'pDC', 'cMono', 'ncMono'), logfc.threshold = 0.1)
# now for our T cell stuff
v2 <- readRDS(object_loc_v2)
DefaultAssay(v2) <- 'RNA'
v2 <- v2[, !is.na(v2@meta.data$timepoint)]
v2 <- v2[, !is.na(v2@meta.data$assignment)]
# add the T classifications
t_classifications_v2 <- read.table(t_classifications_loc_v2, sep = '\t', row.names = 1, stringsAsFactors = F)
v2 <- AddMetaData(v2, t_classifications_loc_v2['clustered.celltype.l2'])
# remove what we don't care about
v2 <- v2[, !is.na(v2@meta.data$clustered.celltype.l2)]
v2 <- v2[, v2@meta.data$clustered.celltype.l2 != 'unknown']
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_t_rna_loc_v2, cell_types_to_use = NULL, logfc.threshold = 0.1)


# read the object
v3 <- readRDS(object_loc_v3)
DefaultAssay(v3) <- 'RNA'
# we've done some refinements at the marker gene level, let's make those changes permanent
v3@meta.data[v3@meta.data$cell_type == 'NK', 'cell_type'] <- 'NKdim'
levels(v3@meta.data$cell_type) <- c(levels(v3@meta.data$cell_type), 'cMono', 'ncMono')
v3@meta.data[v3@meta.data$cell_type %in% c('mono 1', 'mono 4'), 'cell_type'] <- 'cMono'
v3@meta.data[v3@meta.data$cell_type %in% c('mono 2'), 'cell_type'] <- 'ncMono'
v3@meta.data$cell_type <- droplevels(v3@meta.data$cell_type)
# clean up
v3 <- v3[, !is.na(v3@meta.data$cell_type) & !is.na(v3@meta.data$assignment) & !is.na(v3@meta.data$timepoint)]
# write the new object
saveRDS(v3, object_loc_v3_new)
# do the mapping
perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_highres_rna_loc_v3, cell_types_to_use = c('NKdim', 'NKbright', 'mDC', 'pDC', 'cMono', 'ncMono'), logfc.threshold = 0.1)
# now for our T cell stuff
v3 <- readRDS(object_loc_v3)
DefaultAssay(v3) <- 'RNA'
v3 <- v3[, !is.na(v3@meta.data$timepoint)]
v3 <- v3[, !is.na(v3@meta.data$assignment)]
# add the T classifications
t_classifications_v3 <- read.table(t_classifications_loc_v3, sep = '\t', row.names = 1, stringsAsFactors = F)
v3 <- AddMetaData(v3, t_classifications_loc_v3['clustered.celltype.l2'])
# remove what we don't care about
v3 <- v3[, !is.na(v3@meta.data$clustered.celltype.l2)]
v3 <- v3[, v3@meta.data$clustered.celltype.l2 != 'unknown']
perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_t_rna_loc_v3, cell_types_to_use = NULL, logfc.threshold = 0.1)


# also perform the meta analysis
mast_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/'
mast_output_prepend <- paste(mast_output_loc, 'paired_highres_lfc01minpct01_20210905/v', sep = '')
mast_output_append <- '_paired_highres_lfc01minpct01_20210905/rna/'
# write the location of the combined output
mast_meta_output_loc <- paste(mast_output_loc, 'paired_highres_lfc01minpct01_20210905//meta_paired_highres_lfc01minpct01_20210905/rna/', sep = '')

# write meta output
write_meta_mast(mast_output_prepend, mast_output_append, mast_meta_output_loc, cell_types = c('NKdim', 'NKbright', 'cMono', 'ncMono', 'mDC', 'pDC'))

# create plots
overlap_plots <- compare_overlap_low_to_high_res_de_genes('/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/', '/data/scRNA/differential_expression/seurat_MAST/output/paired_highres_lfc01minpct01_20210905/meta_paired_highres_lfc01minpct01_20210905/rna/', list('DC' = c('mDC', 'pDC'), 'monocyte' = c('cMono', 'ncMono'), 'NK' = c('NKdim', 'NKbright')))
# save the result
saveRDS(overlap_plots, paste(mast_overlap_plot_loc,  'overlap_lores_highres.rds', sep = ''))

pdf(paste(mast_overlap_plot_loc,  'overlap_lores_highres_monocyte.pdf', sep = ''), width = 5, height = 5)
overlap_plots[['UT']][['monocyte']]
overlap_plots[['3hCA']][['monocyte']]
overlap_plots[['24hCA']][['monocyte']]
overlap_plots[['3hMTB']][['monocyte']]
overlap_plots[['24hMTB']][['monocyte']]
overlap_plots[['3hPA']][['monocyte']]
overlap_plots[['24hPA']][['monocyte']]
dev.off()

pdf(paste(mast_overlap_plot_loc,  'overlap_lores_highres_NK.pdf', sep = ''), width = 5, height = 5)
overlap_plots[['UT']][['NK']]
overlap_plots[['3hCA']][['NK']]
overlap_plots[['24hCA']][['NK']]
overlap_plots[['3hMTB']][['NK']]
overlap_plots[['24hMTB']][['NK']]
overlap_plots[['3hPA']][['NK']]
overlap_plots[['24hPA']][['NK']]
dev.off()

pdf(paste(mast_overlap_plot_loc,  'overlap_lores_highres_DC.pdf', sep = ''), width = 5, height = 5)
overlap_plots[['UT']][['DC']]
overlap_plots[['3hCA']][['DC']]
overlap_plots[['24hCA']][['DC']]
overlap_plots[['3hMTB']][['DC']]
overlap_plots[['24hMTB']][['DC']]
overlap_plots[['3hPA']][['DC']]
overlap_plots[['24hPA']][['DC']]
dev.off()
