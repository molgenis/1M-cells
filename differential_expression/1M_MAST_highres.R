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


get_significant_genes <- function(mast_output_loc, sig_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv'){
  # get the files
  files <- list.files(mast_output_loc)
  # try to read each file
  for(file in files){
    try({
      # read the mast output
      mast <- read.table(paste(mast_output_loc, file, sep = ''), header=T)
      # filter to only include the significant results
      mast <- mast[mast[[pval_column]] <= sig_pval, ]
      # filter for only the positive lfc if required
      if(only_positive){
        mast <- mast[mast[[lfc_column]] < 0, ]
      }
      # filter for only the positive lfc if required
      if(only_negative){
        mast <- mast[mast[[lfc_column]] > 0, ]
      }
      # confine in some way if reporting a max number of genes
      if(!is.null(max)){
        # by p if required
        if(max_by_pval){
          mast <- mast[order(mast[[p_val_column]]), ]
        }
        # by lfc otherwise
        else{
          mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
        }
        # subset to the number we requested if max was set
        mast <- mast[1:max,]
      }
      # grab the genes from the column names
      genes <- rownames(mast)
      # convert the symbols to ensemble IDs
      if (to_ens) {
        mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
        mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
        genes <- mapping[match(genes, mapping$V2),"V1"]
      }
      # otherwise change the Seurat replacement back
      else{
        #genes <- gsub("-", "_", genes)
      }
      # create a regex to get the last index of .
      last_dot_pos <- "\\.[^\\.]*$"
      # this allows us to remove the filename extention
      file_no_ext <- substr(file, 1, regexpr(last_dot_pos,file)-1)
      # create output location
      sig_output <- paste(sig_output_loc, file_no_ext, '.txt', sep = '')
      # write the genes
      write.table(genes, sig_output, sep = '\t', quote = F, row.names = F, col.names = F)
    })
  }
}

plot_DE_sharing_per_celltype <- function(condition_combination, mast_output_loc, cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), pval_column='metap_bonferroni', sig_pval=0.05, only_positive=F, only_negative=F, lfc_column='metafc', use_label_dict=T, use_color_dict=T){
  DE_genes_per_ct <- list()
  # get the DE genes for each cell type
  for(cell_type in cell_types_to_use){
    # build the full path
    full_mast_path <- paste(mast_output_loc, cell_type, condition_combination, '.tsv', sep = '')
    # grab the significant genes
    try({
      # read the mast output
      mast <- read.table(full_mast_path, header=T, row.names = 1, sep = '\t')
      # filter to only include the significant results
      mast <- mast[mast[[pval_column]] <= 0.05, ]
      # filter for only the positive lfc if required
      if(only_positive){
        mast <- mast[mast[[lfc_column]] < 0, ]
      }
      # filter for only the positive lfc if required
      if(only_negative){
        mast <- mast[mast[[lfc_column]] > 0, ]
      }
      # we just care about the gene names
      sig_genes <- rownames(mast)
      # store these for the cell type
      DE_genes_per_ct[[cell_type]] <- sig_genes
    })
  }
  if(use_label_dict){
    names(DE_genes_per_ct) <- label_dict()[names(DE_genes_per_ct)]
  }
  queries <- NULL
  sets.bar.color <- 'black'
  if(use_color_dict){
    # we can fill a query list now
    queries <- list()
    # create df to store the number of each set, so we know how to order
    nrs_df <- NULL
    # add the colors for the cell types
    for(i in 1:length(names(DE_genes_per_ct))){
      cell_type <- names(DE_genes_per_ct)[i]
      # add for the singles in the intersection sizes
      ct_list <- list(
        query = intersects,
        params = list(cell_type),
        color = get_color_coding_dict()[[cell_type]],
        active = T)
      queries[[i]] <- ct_list
      # add for the DF to order the set sizes
      numbers_row <- data.frame(ct=c(cell_type), nr=c(length(DE_genes_per_ct[[cell_type]])), stringsAsFactors = F)
      if(is.null(nrs_df)){
        nrs_df <- numbers_row
      }
      else{
        nrs_df <- rbind(nrs_df, numbers_row)
      }
    }
    # get the order of the sets
    ordered_cts <- nrs_df[order(nrs_df$nr, decreasing = T), 'ct']
    # add the colors for the sets
    sets.bar.color <- unlist(get_color_coding_dict()[ordered_cts])
  }
  upset(fromList(DE_genes_per_ct), order.by = 'freq', nsets = length(DE_genes_per_ct), queries = queries, sets.bar.color=sets.bar.color	)
  #return(DE_genes_per_ct)
}

plot_DE_sharing_per_celltype_meh <- function(condition_combination, mast_output_loc, cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), pval_column='metap_bonferroni', sig_pval=0.05, only_positive=F, only_negative=F, lfc_column='metafc', use_label_dict=T, use_color_dict=T, flip_shared_unique=F){
  DE_genes_per_ct <- list()
  # get the DE genes for each cell type
  for(cell_type in cell_types_to_use){
    # build the full path
    full_mast_path <- paste(mast_output_loc, cell_type, condition_combination, '.tsv', sep = '')
    # grab the significant genes
    try({
      # read the mast output
      mast <- read.table(full_mast_path, header=T, row.names = 1, sep = '\t')
      # filter to only include the significant results
      mast <- mast[mast[[pval_column]] <= 0.05, ]
      # filter for only the positive lfc if required
      if(only_positive){
        mast <- mast[mast[[lfc_column]] < 0, ]
      }
      # filter for only the positive lfc if required
      if(only_negative){
        mast <- mast[mast[[lfc_column]] > 0, ]
      }
      # we just care about the gene names
      sig_genes <- rownames(mast)
      # store these for the cell type
      DE_genes_per_ct[[cell_type]] <- sig_genes
    })
  }
  if(use_label_dict){
    names(DE_genes_per_ct) <- label_dict()[names(DE_genes_per_ct)]
  }
  # init the dataframe
  numbers <- NULL
  for(cell_type in names(DE_genes_per_ct)){
    # make a copy of everything
    DE_genes_per_ct_not_this_ct <- DE_genes_per_ct
    # remove this specific cell type from the copy
    DE_genes_per_ct_not_this_ct[[cell_type]] <- NULL
    # combine all the genes (which dont containt this cell type specific genes anymore)
    DE_genes_not_this_ct <- unique(as.vector(unlist(DE_genes_per_ct_not_this_ct)))
    # get the genes specific for this cell type
    DE_genes_ct <- DE_genes_per_ct[[cell_type]]
    # check which are unique to the cell type
    DE_genes_ct_unique <- setdiff(DE_genes_ct, DE_genes_not_this_ct)
    # put it into numbers
    DE_genes_ct_unique_number <- length(DE_genes_ct_unique)
    DE_genes_ct_shared_number <- length(DE_genes_ct) - DE_genes_ct_unique_number
    # turn into rows
    rows <- data.frame(cell_type=c(cell_type, cell_type), unique=c(cell_type, 'shared'), number=c(DE_genes_ct_unique_number, DE_genes_ct_shared_number), stringsAsFactors = F)
    # add to dataframe
    if(is.null(numbers)){
      numbers <- rows
    }
    else{
      numbers <- rbind(numbers, rows)
    }
  }
  # set the order I like for the legend, but setting the factor order
  numbers$unique <- factor(numbers$unique, levels=c('shared', names(DE_genes_per_ct)))
  if(flip_shared_unique){
    numbers$unique <- factor(numbers$unique, levels=c(names(DE_genes_per_ct), 'shared'))
  }
  # set default label
  condition_combination_label <- condition_combination
  # make a prettier one if we can
  if(use_label_dict){
    # get a nice label for the condition combination
    condition_combination_label <- label_dict()[[condition_combination]]
  }
  # make the plot finally
  p <- ggplot(numbers, aes(fill=unique, y=number, x=cell_type)) +
    geom_bar(position='stack', stat='identity') +
    labs(x='cell type', y='number of DE genes') +
    ggtitle(paste('DE genes and cell type specificity in', condition_combination_label)) +
    labs(fill = "Found in")
  if(use_color_dict){
    # grab the colours
    cc <- get_color_coding_dict()
    # add the 'mixed' condition
    cc[['shared']] <- 'gray'
    fillScale <- scale_fill_manual(name = "cell type",values = unlist(cc[c('shared', names(DE_genes_per_ct))]))
    p <- p + fillScale
  }
  return(p)
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
mast_output_paired_t_loc_v2 <- paste(mast_output_loc, 'paired_t_lfc01minpct01_20210905/v2_paired_t_lfc01minpct01_20210905/', sep = '')
mast_output_paired_t_loc_v3 <- paste(mast_output_loc, 'paired_t_lfc01minpct01_20210905/v3_paired_t_lfc01minpct01_20210905/', sep = '')
mast_output_paired_t_rna_loc_v2 <- paste(mast_output_paired_highres_loc_v2, 'rna/', sep = '')
mast_output_paired_t_rna_loc_v3 <- paste(mast_output_paired_highres_loc_v3, 'rna/', sep = '')
# and the assignments for these T cells needs to come from somewhere as well
t_classifications_loc <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/seurat_multimodal/classifications/'
t_classifications_loc_v2 <- paste(t_classifications_loc, 'v2_azi_to_cluster_l2_cell_types.tsv', sep = '')
t_classifications_loc_v3 <- paste(t_classifications_loc, 'v3_azi_to_cluster_l2_cell_types.tsv', sep = '')

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
t_classifications_v2 <- read.table(t_classifications_loc_v2, sep = '\t', row.names = 1, stringsAsFactors = F, header = T)
v2 <- AddMetaData(v2, t_classifications_v2['clustered.celltype.l2.t'])
# remove what we don't care about
v2 <- v2[, !is.na(v2@meta.data$clustered.celltype.l2.t)]
v2 <- v2[, v2@meta.data$clustered.celltype.l2.t != 'unknown']
# reclassify the memory cells together
v2@meta.data$clustered.celltype.l2.t.merged <- v2@meta.data$clustered.celltype.l2.t
v2@meta.data[v2@meta.data$clustered.celltype.l2.t %in% c('CD4 TCM', 'CD4 TEM'), 'clustered.celltype.l2.t.merged'] <- 'CD4 Memory'
v2@meta.data[v2@meta.data$clustered.celltype.l2.t %in% c('CD8 TCM', 'CD8 TEM'), 'clustered.celltype.l2.t.merged'] <- 'CD8 Memory'
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_t_rna_loc_v2, cell.type.column = 'clustered.celltype.l2.t.merged', cell_types_to_use = NULL, logfc.threshold = 0.1)
rm(v2)

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
t_classifications_v3 <- read.table(t_classifications_loc_v3, sep = '\t', row.names = 1, stringsAsFactors = F, header = T)
v3 <- AddMetaData(v3, t_classifications_v3['clustered.celltype.l2.t'])
# remove what we don't care about
v3 <- v3[, !is.na(v3@meta.data$clustered.celltype.l2.t)]
v3 <- v3[, v3@meta.data$clustered.celltype.l2.t != 'unknown']
# reclassify the memory cells together
v3@meta.data$clustered.celltype.l2.t.merged <- v3@meta.data$clustered.celltype.l2.t
v3@meta.data[v3@meta.data$clustered.celltype.l2.t %in% c('CD4 TCM', 'CD4 TEM'), 'clustered.celltype.l2.t.merged'] <- 'CD4 Memory'
v3@meta.data[v3@meta.data$clustered.celltype.l2.t %in% c('CD8 TCM', 'CD8 TEM'), 'clustered.celltype.l2.t.merged'] <- 'CD8 Memory'
perform_mast_per_celltype(seurat_object = v3, output_loc = mast_output_paired_t_rna_loc_v3, cell.type.column = 'clustered.celltype.l2.t.merged', cell_types_to_use = NULL, logfc.threshold = 0.1)


# also perform the meta analysis
mast_output_local_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/'
mast_output_prepend <- paste(mast_output_loc, 'paired_highres_lfc01minpct01_20210905/v', sep = '')
mast_output_append <- '_paired_highres_lfc01minpct01_20210905/rna/'
# write the location of the combined output
mast_meta_output_loc <- paste(mast_output_loc, 'paired_highres_lfc01minpct01_20210905//meta_paired_highres_lfc01minpct01_20210905/rna/', sep = '')
mast_meta_output_local_loc <- paste(mast_output_local_loc, 'paired_highres_lfc01minpct01_20210905//meta_paired_highres_lfc01minpct01_20210905/rna/', sep = '')

# write meta output
write_meta_mast(mast_output_prepend, mast_output_append, mast_meta_output_loc, cell_types = c('NKdim', 'NKbright', 'cMono', 'ncMono', 'mDC', 'pDC'))

# now for T as well
t_mast_output_prepend <- paste(mast_meta_output_local_loc, 'paired_t_lfc01minpct01_20210905/v', sep = '')
t_mast_output_append <- '_paired_t_lfc01minpct01_20210905/rna/'
t_mast_meta_output_loc <- paste(mast_output_local_loc, 'paired_t_lfc01minpct01_20210905/meta_paired_t_lfc01minpct01_20210905/rna/', sep = '')
# these are the labels we want to look at
labels_t_azimuth <- c('CD8 Naive', 'CD4 CTL', 'CD4 Naive', 'CD4 TCM', 'CD8 TEM', 'MAIT', 'CD8 TEM', 'dnT', 'CD4 TEM', 'CD8 TCM', 'CD8 Proliferating')
labels_memory <- c('CD4 Memory', 'CD8 Memory')
# write the meta results
write_meta_mast(t_mast_output_prepend, t_mast_output_append, t_mast_meta_output_loc, cell_types = labels_t_azimuth)


sig_output_loc_hires <- '/data/scRNA/differential_expression/sigs/meta_paired_highres_lfc01minpct01_20210905/rna/'
sig_down_output_loc_hires <- '/data/scRNA/differential_expression/sigs_pos/meta_paired_highres_lfc01minpct01_20210905/rna/'
sig_up_output_loc_hires <- '/data/scRNA/differential_expression/sigs_neg/meta_paired_highres_lfc01minpct01_20210905/rna/'
# write the sinificant genes to files
get_significant_genes(mast_meta_output_local_loc, sig_output_loc_hires)
get_significant_genes(mast_meta_output_local_loc, sig_down_output_loc_hires, only_negative = T)
get_significant_genes(mast_meta_output_local_loc, sig_up_output_loc_hires, only_positive = T)
# we will write the significant T genes to files
sig_output_loc_t <- '/data/scRNA/differential_expression/sigs/meta_paired_t_lfc01minpct01_20210905/rna/'
sig_down_output_loc_t <- '/data/scRNA/differential_expression/sigs_pos/meta_paired_t_lfc01minpct01_20210905/rna/'
sig_up_output_loc_t <- '/data/scRNA/differential_expression/sigs_neg/meta_paired_t_lfc01minpct01_20210905/rna/'
# write the sinificant genes to files
get_significant_genes(t_mast_meta_output_loc, sig_output_loc_t)
get_significant_genes(t_mast_meta_output_loc, sig_down_output_loc_t, only_negative = T)
get_significant_genes(t_mast_meta_output_loc, sig_up_output_loc_t, only_positive = T)


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

# and for T
overlap_plots_t <- compare_overlap_low_to_high_res_de_genes('/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/', '/data/scRNA/differential_expression/seurat_MAST/output/paired_t_lfc01minpct01_20210905/meta_paired_t_lfc01minpct01_20210905/rna/', list('CD4T' = c('CD4 CTL', 'CD4 Naive', 'CD4 TCM'), 'CD8T' = c('CD8 Naive', 'CD8 TEM', 'CD8 TEM', 'CD8 TCM', 'CD8 Proliferating')))

pdf(paste(mast_overlap_plot_loc,  'overlap_lores_highres_CD4T.pdf', sep = ''), width = 5, height = 5)
overlap_plots_t[['UT']][['CD4T']]
overlap_plots_t[['3hCA']][['CD4T']]
overlap_plots_t[['24hCA']][['CD4T']]
overlap_plots_t[['3hMTB']][['CD4T']]
overlap_plots_t[['24hMTB']][['CD4T']]
overlap_plots_t[['3hPA']][['CD4T']]
overlap_plots_t[['24hPA']][['CD4T']]
dev.off()

pdf(paste(mast_overlap_plot_loc,  'overlap_lores_highres_CD8T.pdf', sep = ''), width = 5, height = 5)
overlap_plots_t[['UT']][['CD8T']]
overlap_plots_t[['3hCA']][['CD8T']]
overlap_plots_t[['24hCA']][['CD8T']]
overlap_plots_t[['3hMTB']][['CD8T']]
overlap_plots_t[['24hMTB']][['CD8T']]
overlap_plots_t[['3hPA']][['CD8T']]
overlap_plots_t[['24hPA']][['CD8T']]
dev.off()

# plot the amount of sharing
pdf(paste(mast_overlap_plot_loc, 'hires_t_sharing.pdf', sep = ''), width = 5, height = 5)
plot_DE_sharing_per_celltype(mast_output_loc = t_mast_meta_output_loc, condition_combination = 'UTX3hCA', cell_types_to_use = labels_t_azimuth, use_label_dict = F, use_color_dict = F)
plot_DE_sharing_per_celltype(mast_output_loc = t_mast_meta_output_loc, condition_combination = 'UTX24hCA', cell_types_to_use = labels_t_azimuth, use_label_dict = F, use_color_dict = F)
plot_DE_sharing_per_celltype(mast_output_loc = t_mast_meta_output_loc, condition_combination = 'UTX3hMTB', cell_types_to_use = labels_t_azimuth, use_label_dict = F, use_color_dict = F)
plot_DE_sharing_per_celltype(mast_output_loc = t_mast_meta_output_loc, condition_combination = 'UTX24hMTB', cell_types_to_use = labels_t_azimuth, use_label_dict = F, use_color_dict = F)
plot_DE_sharing_per_celltype(mast_output_loc = t_mast_meta_output_loc, condition_combination = 'UTX3hPA', cell_types_to_use = labels_t_azimuth, use_label_dict = F, use_color_dict = F)
plot_DE_sharing_per_celltype(mast_output_loc = t_mast_meta_output_loc, condition_combination = 'UTX24hPA', cell_types_to_use = labels_t_azimuth, use_label_dict = F, use_color_dict = F)
dev.off()



