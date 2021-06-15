

library(Seurat)
library(MAST)
library(UpSetR)
library(data.table)


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

perform_mast_per_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL, condition.1 = '1', condition.2 = '2'){
  # do subselection based on features
  features = NULL
  # grab the top expressed genes if that is what the user wanted
  if(!is.null(use_top_expressed)){
    features <- get_top_expressed_features(seurat_object, use_top_expressed)
  }
  # go through the cell types
  for(cell_type in unique(seurat_object@meta.data[[cell.type.column]])){
    # make a more specific path
    output_loc_cell_type <- paste(output_loc, cell_type, sep = '')
    # grab the subset
    seurat_object_cell_type <- seurat_object[,seurat_object@meta.data[cell.type.column] == cell_type]
    # do the MAST
    try({perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = condition.1, condition.2 = condition.2 ,split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, features = features, latent.vars = latent.vars)})
  }
  # finally do a bulk analysis as well
  output_loc_bulk <- paste(output_loc, 'bulk', sep = '')
  try({perform_mast(seurat_object, output_loc_bulk, condition.1 = condition.1, condition.2 = condition.2, split.column = split.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, latent.vars = latent.vars)})
}

perform_mast_per_celltype_subsampled <- function(seurat_object, seurat_object2, output_loc, subsample_size, subsample_times=10, subsample_column='assignment', split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL, condition.1 = '1', condition.2 = '2'){
  for(i in 1:subsample_times){
    # set the location
    output_loc_ss <- paste(output_loc, i, '_', sep = '')
    # grab a subsample
    subsampled_assignments <- sample(unique(seurat_object@meta.data[[subsample_column]]), subsample_size)
    # grab the subsampled object
    seurat_object_subsampled <- seurat_object[, seurat_object@meta.data[[subsample_column]] %in% subsampled_assignments]
    # merge the two objects
    seurat_object_merged <- merge(seurat_object_subsampled, seurat_object2)
    # normalize data
    seurat_object_merged <- NormalizeData(seurat_object_merged)
    # actually perform
    perform_mast_per_celltype(seurat_object_merged, output_loc_ss, split.column = split.column, cell.type.column = cell.type.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, use_top_expressed = use_top_expressed, latent.vars=latent.vars, condition.1 = condition.1, condition.2 = condition.2)
  }
}

perform_mast_per_celltype_subsampled_same <- function(seurat_object, output_loc, subsample_size, subsample_times=10, subsample_column='assignment', split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL, condition.1='1', condition.2='2'){
  for(i in 1:subsample_times){
    # set the location
    output_loc_ss <- paste(output_loc, i, '_', sep = '')
    # grab a subsample
    subsampled_assignments <- as.character(sample(unique(seurat_object@meta.data[[subsample_column]]), subsample_size))
    subsampled_assignments_split <- split(subsampled_assignments, sample(2, length(subsampled_assignments), repl = TRUE) )
    subsampled_assignments_1 <- subsampled_assignments_split[[1]]
    subsampled_assignments_2 <- subsampled_assignments_split[[2]]
    # grab the subsampled object
    seurat_object_subsampled <- seurat_object[, seurat_object@meta.data[[subsample_column]] %in% subsampled_assignments]
    seurat_object_subsampled@meta.data[seurat_object_subsampled@meta.data[[subsample_column]] %in% subsampled_assignments_1, ][[split.column]] <- condition.1
    seurat_object_subsampled@meta.data[seurat_object_subsampled@meta.data[[subsample_column]] %in% subsampled_assignments_2, ][[split.column]] <- condition.2
    # normalize data
    seurat_object_subsampled <- NormalizeData(seurat_object_subsampled)
    # actually perform
    perform_mast_per_celltype(seurat_object_subsampled, output_loc_ss, split.column = split.column, cell.type.column = cell.type.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, use_top_expressed = use_top_expressed, latent.vars=latent.vars, condition.1 = condition.1, condition.2 = condition.2)
  }
}

perform_mast_per_celltype_subsampled_tp2 <- function(seurat_object, output_loc, subsample_size, subsample_times=10, subsample_column='assignment', split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL, condition.1='1', condition.2='2'){
  for(i in 1:subsample_times){
    # set the location
    output_loc_ss <- paste(output_loc, i, '_', sep = '')
    # grab a subsample
    subsampled_assignments <- as.character(sample(unique(seurat_object@meta.data[[subsample_column]]), subsample_size))
    subsampled_assignments_split <- split(subsampled_assignments, sample(2, length(subsampled_assignments), repl = TRUE) )
    subsampled_assignments_1 <- subsampled_assignments_split[[1]]
    subsampled_assignments_2 <- subsampled_assignments_split[[2]]
    # grab the subsampled object
    seurat_object_subsampled <- seurat_object[, seurat_object@meta.data[[subsample_column]] %in% subsampled_assignments]
    # subsample to the two conditions
    seurat_object_subsampled <- seurat_object_subsampled[, ((seurat_object_subsampled@meta.data[[subsample_column]] %in% subsampled_assignments_1 & seurat_object_subsampled@meta.data[[split.column]] == condition.1) | (seurat_object_subsampled@meta.data[[subsample_column]] %in% subsampled_assignments_2 & seurat_object_subsampled@meta.data[[split.column]] == condition.2))]
    # normalize data
    seurat_object_subsampled <- NormalizeData(seurat_object_subsampled)
    # actually perform
    perform_mast_per_celltype(seurat_object_subsampled, output_loc_ss, split.column = split.column, cell.type.column = cell.type.column, assay = assay, min.pct = min.pct, logfc.threshold = logfc.threshold, use_top_expressed = use_top_expressed, latent.vars=latent.vars, condition.1 = condition.1, condition.2 = condition.2)
  }
}


combine_mast_results_ss <- function(mast_output_loc, merged_output_loc, condition.1 = '1M_cells', condition.2 = 'pilot4_unstimulated', cell_types_to_check=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
  for(cell_type in cell_types_to_check){
    # create regex to list the files
    list_dir_regex <- paste('*_', cell_type, condition.1, condition.2, '.tsv', sep = '')
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
      # add prepend
      colnames(de_results) <- paste(i, colnames(de_results), sep = '_')
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

get_nr_of_times_DE <- function(output_loc, cell_type){
  # create regex to list the files
  list_dir_regex <- cell_type
  # list the files
  files <- list.files(output_loc)
  # filter
  files <- files[grepl(list_dir_regex, files)]
  # there should be only one result
  file <- files[[1]]
  # append file loc
  file_loc <- paste(output_loc, file, sep = '')
  # read the file
  de_results <- read.table(file_loc, sep = '\t', header = T, row.names = 1)
  # get the distribution of DE genes
  nr_of_times_de <- apply(de_results, 1, function(x){
    # grab the p_val_adj columns
    x <- x[names(x)[grep('p_val_adj', names(x))]]
    # check how often they are significant
    x_sig <- sum(!is.na(x) & x < 0.05)
    return(x_sig)
  })
  return(nr_of_times_de)
}

plot_DE_distributions <- function(merged_mast_output_locs, comparison_names, cell_types_to_check=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
  # make a plot per cell type
  for(cell_type in cell_types_to_check){
    # display all the plots for the cell type
    par(mfrow=c(1,length(merged_mast_output_locs)))
    # read each output file
    for(i in 1:length(merged_mast_output_locs)){
      output_loc <- merged_mast_output_locs[i]
      # grab the numer of times each gene was DE in all the subsampling
      nr_of_times_de <- get_nr_of_times_DE(output_loc, cell_type)
      # add a plot for this specific comparison
      hist(nr_of_times_de, main = paste(cell_type, comparison_names[i]))
    }
  }
}



plot_DE_overlaps <- function(merged_mast_output_locs, comparison_names, cell_types_to_check=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), true_de_cutoff=15){
  # make a plot per cell type
  for(cell_type in cell_types_to_check){
    # create a list to put the DE genes in
    de_per_comparison <- list()
    for(i in 1:length(merged_mast_output_locs)){
      # grab the output location
      output_loc <- merged_mast_output_locs[i]
      # grab the numer of times each gene was DE in all the subsampling
      nr_of_times_de <- get_nr_of_times_DE(output_loc, cell_type)
      # grab the gene names that were DE more times than given in the threshold
      true_de_comparison <- names(nr_of_times_de)[nr_of_times_de >= true_de_cutoff]
      # grab the name of the comparison
      comparison_name <- comparison_names[i]
      # add these genes to the list
      de_per_comparison[[comparison_name]] <- true_de_comparison
    }
    upset(fromList(de_per_comparison), order.by = 'freq')
  }
}

get_DE_overlap_plot_ct<- function(merged_mast_output_locs, comparison_names, cell_type, true_de_cutoff=15){
  # create a list to put the DE genes in
  de_per_comparison <- list()
  for(i in 1:length(merged_mast_output_locs)){
    # grab the output location
    output_loc <- merged_mast_output_locs[i]
    # grab the numer of times each gene was DE in all the subsampling
    nr_of_times_de <- get_nr_of_times_DE(output_loc, cell_type)
    # grab the gene names that were DE more times than given in the threshold
    true_de_comparison <- names(nr_of_times_de)[nr_of_times_de >= true_de_cutoff]
    # grab the name of the comparison
    comparison_name <- comparison_names[i]
    # add these genes to the list
    de_per_comparison[[comparison_name]] <- true_de_comparison
  }
  upset(fromList(de_per_comparison), order.by = 'freq')
}


output_shared_and_exclusive_de_genes <- function(merged_mast_output_loc1, merged_mast_output_loc2, comparison_name1, comparison_name2, output_loc, cell_types_to_check=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), true_de_cutoff=15){
  # make a file per cell type
  for(cell_type in cell_types_to_check){
    # grab the numer of times each gene was DE in all the subsampling
    nr_of_times_de1 <- get_nr_of_times_DE(merged_mast_output_loc1, cell_type)
    # grab the gene names that were DE more times than given in the threshold
    true_de_comparison1 <- names(nr_of_times_de1)[nr_of_times_de1 >= true_de_cutoff]
    # grab the numer of times each gene was DE in all the subsampling
    nr_of_times_de2 <- get_nr_of_times_DE(merged_mast_output_loc2, cell_type)
    # grab the gene names that were DE more times than given in the threshold
    true_de_comparison2 <- names(nr_of_times_de2)[nr_of_times_de2 >= true_de_cutoff]
    # these were the common DE genes
    common_de <- intersect(true_de_comparison1, true_de_comparison2)
    # these were the set1 exclusive genes
    set1_exclusive_de <- setdiff(true_de_comparison1, true_de_comparison2)
    # these were the set2 exclusive genes
    set2_exclusive_de <- setdiff(true_de_comparison2, true_de_comparison1)
    # create the output paths
    common_de_output_path <- paste(output_loc, cell_type, '_', comparison_name1, '_vs_', comparison_name2, '_common.txt', sep = '')
    set1_de_output_path <- paste(output_loc, cell_type, '_', comparison_name1, '_vs_', comparison_name2, '_', comparison_name1, '_exclusive.txt', sep = '')
    set2_de_output_path <- paste(output_loc, cell_type, '_', comparison_name1, '_vs_', comparison_name2, '_', comparison_name2, '_exclusive.txt', sep = '')
    # write the files
    write.table(common_de, common_de_output_path, col.names = F, row.names = F, quote = F)
    write.table(set1_exclusive_de, set1_de_output_path, col.names = F, row.names = F, quote = F)
    write.table(set2_exclusive_de, set2_de_output_path, col.names = F, row.names = F, quote = F)
  }
}




# get the old object
pilot4_old_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/ut_compare/objects/pilot4.Rds'
pilot4_old <- readRDS(pilot4_old_loc)
# add cell_types
pilot4_old@meta.data$cell_type <- pilot4_old@ident
# convert to gene symbols
gene_to_ens_mapping <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/resources/features_v3.tsv"
genes <- read.table(gene_to_ens_mapping, header = F, stringsAsFactors = F)
genes$V2 <- gsub("_", "-", make.unique(genes$V2))
gene_symbols <- genes[match(rownames(pilot4_old@data), genes$V1),"V2"]
rownames(pilot4_old@data) <- gene_symbols
rownames(pilot4_old@raw.data) <- gene_symbols
rownames(pilot4_old@scale.data) <- gene_symbols
# recreate seurat object
pilot4 <- CreateSeuratObject(counts = pilot4_old@raw.data)
pilot4@meta.data <- pilot4_old@meta.data
# update to new version
#pilot4 <- UpdateSeuratObject(pilot4_old)
# grab just UT
pilot4_ut <- pilot4[, !is.na(pilot4@meta.data$stimulation) & pilot4@meta.data$stimulation == 'unstimulated']
# clear some memory
rm(pilot4)
# get the new v2 object
v2_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds'
v2 <- readRDS(v2_loc)
# grab just UT
v2_ut <- v2[, !is.na(v2@meta.data$timepoint) & v2@meta.data$timepoint == 'UT']
# grab UT and 24hMTB
v2_ut_24hmtb <- v2[, !is.na(v2@meta.data$timepoint) & (v2@meta.data$timepoint == 'UT' | v2@meta.data$timepoint == 'X24hMTB')]
# clear some memory
rm(v2)
# harmonise the cell types
v2_ut@meta.data$cell_type <- v2_ut@meta.data$cell_type_lowerres
#pilot4_ut@meta.data$cell_type <- pilot4_ut@active.ident
levels(pilot4_ut@meta.data$cell_type) <- c(levels(pilot4_ut@meta.data$cell_type), 'CD4T', 'CD8T', 'DC', 'monocyte')
pilot4_ut@meta.data[pilot4_ut@meta.data$cell_type == 'CD4+_T', ]$cell_type <- 'CD4T'
pilot4_ut@meta.data[pilot4_ut@meta.data$cell_type == 'CD8+_T', ]$cell_type <- 'CD8T'
pilot4_ut@meta.data[pilot4_ut@meta.data$cell_type == 'pDC', ]$cell_type <- 'DC'
pilot4_ut@meta.data[pilot4_ut@meta.data$cell_type == 'Monocyte', ]$cell_type <- 'monocyte'
pilot4_ut@meta.data$cell_type <- droplevels(pilot4_ut@meta.data$cell_type)

# merge the objects
ut_merged <- merge(pilot4_ut, v2_ut)
#normalize
ut_merged <- NormalizeData(ut_merged)

# set the location of the output
mast_output_loc_full <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/ut_compare/mast_output/full_'
mast_output_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/ut_compare/mast_output/'
# do all vs all
perform_mast_per_celltype(seurat_object = ut_merged, output_loc = mast_output_loc_full, split.column = 'orig.ident', cell.type.column = 'cell_type', assay = 'RNA', min.pct = 0.1, logfc.threshold = 0.25)
# perform with subsampling, taking ten participants each time
perform_mast_per_celltype_subsampled(v2_ut, pilot4_ut, output_loc = mast_output_loc, subsample_size=10, subsample_times=20, subsample_column='assignment', split.column = 'orig.ident', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL, condition.1 = 'pilot4_unstimulated', condition.2 = '1M_cells')
# to compare to differences in the dataset itself, do the same by subsampling from the 1M object and comparing then
perform_mast_per_celltype_subsampled_same(v2_ut, output_loc = mast_output_loc, subsample_size=20, subsample_times=20, subsample_column='assignment', split.column = 'orig.ident', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL)
# now compare two actual timepoints in the same dataset
perform_mast_per_celltype_subsampled_tp2(v2_ut_24hmtb, output_loc = mast_output_loc, subsample_size=20, subsample_times=20, subsample_column='assignment', split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', stims = NULL, min.pct = 0.1, logfc.threshold = 0.25, use_top_expressed = NULL, latent.vars=NULL, condition.1 = 'UT', condition.2 = 'X24hMTB')


# set location to store the combined MAST output
mast_output_loc_merged <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/ut_compare/mast_output_merged/'
# created the merged columns
combine_mast_results_ss(mast_output_loc, mast_output_loc_merged, condition.1 = 'pilot4_unstimulated', condition.2 = '1M_cells',  cell_types_to_check=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))

# set location to store the combined MAST output for the subset of the same 1M
mast_output_loc_merged_ss <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/ut_compare/mast_output_merged_ss_self/'
combine_mast_results_ss(mast_output_loc, mast_output_loc_merged_ss, condition.1 = '1', condition.2 = '2', cell_types_to_check=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))

# set location to store the combined MAST output for the subset of the same 1M
mast_output_loc_merged_2tp <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/ut_compare/mast_output_merged_ss_2tp/'
combine_mast_results_ss(mast_output_loc, mast_output_loc_merged_2tp, condition.1 = 'UT', condition.2 = 'X24hMTB', cell_types_to_check=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))

# set vectors of locations and the names
merged_mast_output_locs <- c(mast_output_loc_merged, mast_output_loc_merged_ss, mast_output_loc_merged_2tp)
merged_mast_output_locs <- c('/data/scRNA/ut_compare/mast_output_merged/', '/data/scRNA/ut_compare/mast_output_merged_ss_self/', '/data/scRNA/ut_compare/mast_output_merged_ss_2tp/') # I copied these locally
comparison_names <- c('Plos UT vs 1M UT subsampled 20*10', '1M UT vs 1M UT subsampled 20*10', '1M UT vs 24hMTB subsampled 20*10')


# plot the number of times a DE gene was found back to be DE in subsampling
plot_DE_distributions(merged_mast_output_locs, comparison_names)

# upsetR has issues with plotting in loops, so unfortunately we will have to manually create the plots per cell type
merged_mast_output_locs <- c(mast_output_loc_merged, mast_output_loc_merged_ss, mast_output_loc_merged_2tp)
merged_mast_output_locs <- c('/data/scRNA/ut_compare/mast_output_merged/', '/data/scRNA/ut_compare/mast_output_merged_ss_self/', '/data/scRNA/ut_compare/mast_output_merged_ss_2tp/') # I copied these locally
comparison_names <- c('vsp4', 'vsself', 'vsreal')
# so per cell type
get_DE_overlap_plot_ct(merged_mast_output_locs, comparison_names, cell_type='B')
get_DE_overlap_plot_ct(merged_mast_output_locs, comparison_names, cell_type='CD4T')
get_DE_overlap_plot_ct(merged_mast_output_locs, comparison_names, cell_type='CD8T')
get_DE_overlap_plot_ct(merged_mast_output_locs, comparison_names, cell_type='DC')
get_DE_overlap_plot_ct(merged_mast_output_locs, comparison_names, cell_type='monocyte')
get_DE_overlap_plot_ct(merged_mast_output_locs, comparison_names, cell_type='NK')

# next we should check to see which genes are specifically shared or exclusive
gene_list_output_loc <- '/data/scRNA/ut_compare/de_gene_sharing/'

# 
output_shared_and_exclusive_de_genes('/data/scRNA/ut_compare/mast_output_merged/', '/data/scRNA/ut_compare/mast_output_merged_ss_2tp/', 'P4UT1MUT', '1MUT1M24hMTB', gene_list_output_loc)


