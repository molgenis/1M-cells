######################
# libraries          #
######################

library(Seurat)
library(Matrix)

####################
# Functions        #
####################

create_features_files <- function(seurat_object,  output_loc, cell_types_to_output = NULL, symbols.to.ensg=F, symbols.to.ensg.mapping = "genes.tsv", sample.id.column.name="assignment", cell_type_column = "cell_type", assay = "SCT", prepend_1 = T){
  # set the assay
  DefaultAssay(seurat_object) <- assay
  # get the individuals in the object
  individuals <- unique(seurat_object@meta.data[,sample.id.column.name])
  # remove any na individuals (sometimes caused if using demux assignments)
  individuals <- individuals[!is.na(individuals)]
  # if the user did not specify the cell types to output, just do all of them
  if(is.null(cell_types_to_output)){
    cell_types_to_output = unique(seurat_object@meta.data[[cell_type_column]])
  }
  # we will want to report on the number of cells per cell type, so that needs to be stored in a matrix
  cell_counts <- matrix(nrow=length(cell_types_to_output), ncol = length(individuals), 
                        dimnames = list(cell_types_to_output, individuals))
  # go through the cell types we want to output
  for (celltype in cell_types_to_output) {
    # grab the cells specific to that cell type
    cells_cell_type <- seurat_object[,seurat_object@meta.data[cell_type_column] == celltype]
    # create a matrix where we will store the mean expressions
    mean_expression_matrix <- matrix(nrow=nrow(cells_cell_type), ncol = length(individuals), 
                                     dimnames = list(rownames(cells_cell_type), individuals))
    # go through the individuals
    #for (individual in individuals) {
    cells_individuals <- unique(cells_cell_type@meta.data[[sample.id.column.name]])
    cells_individuals <- cells_individuals[!is.na(cells_individuals)]
    for (individual in cells_individuals) {
      # take into account that there might be zero expression
      if (sum(is.na(cells_cell_type@meta.data[,sample.id.column.name]) | cells_cell_type@meta.data[,sample.id.column.name] == individual) == 0) {
        mean_expression_matrix[,individual] <- 0
        cell_counts[celltype,individual] <- 0
      }
      else {
        # grab the cells of this cell type and this individual
        cells_cell_type_individual <- cells_cell_type[,cells_cell_type@meta.data[,sample.id.column.name] == individual]
        # calculate the mean expression for each gene
        mean_expression_matrix[,individual] <- rowMeans(cells_cell_type_individual)
        # get the number of cells of this cell type for this individual
        cell_counts[celltype,individual] <- ncol(cells_cell_type_individual)
      }
    }
    # remap symbols to ensg numbers if necessary
    if (symbols.to.ensg) {
      genes <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
      genes$V2 <- gsub("_", "-", make.unique(genes$V2))
      rownames(mean_expression_matrix) <- genes[match(rownames(mean_expression_matrix), genes$V2),"V1"]
    }
    if(prepend_1){
      # prepend the '1_' to the id 
      colnames(mean_expression_matrix) <- paste0("1_", colnames(mean_expression_matrix))
    }
    # write our table of means for this cell type
    write.table(mean_expression_matrix, 
                file = paste0(output_loc, celltype, "_expression", ".tsv"),
                quote = F, sep = "\t", col.names = NA)
  }
  # write out table of 
  write.table(cell_counts, 
              file = paste0(output_loc, "cell_counts.txt"),
              quote = F, sep = "\t", col.names = NA)
}

create_feature_files_per_condition <- function(seurat_object, output_loc, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint", sample.id.column.name="assignment", cell_type_column = "cell_type", assay = "SCT", symbols.to.ensg = F, symbols.to.ensg.mapping="genes.tsv", prepend_1 = T){
  # grab the conditions
  conditions <- unique(seurat_object@meta.data[[condition.column.name]])
  # remove any na conditions (sometimes caused if using demux assignments)
  conditions <- conditions[!is.na(conditions)]
  # unless we wnat only specific conditions, we'll just to those
  if(!(is.null(conditions_to_output))){
    conditions <- conditions_to_output
  }
  # go through the conditions
  for(condition in conditions){
    # there should be a subdirectory per condition
    output_loc_condition <- paste0(output_loc, condition,"/")
    # there might be a prepended 'X' that we might want to remove
    if(startsWith(condition, "X")){
      output_loc_condition <- paste0(output_loc, substr(condition, 2, nchar(condition)),"/")
    }
    # grab only the cells with this condition
    print(paste("grabbing condition", condition))
    cells_condition <- seurat_object[,seurat_object@meta.data[condition.column.name] == condition]
    print(paste("finished grabbing condition", condition))
    # do the feature file creation for this condition
    create_features_files(cells_condition, cell_types_to_output = NULL, output_loc = output_loc_condition, symbols.to.ensg=symbols.to.ensg, symbols.to.ensg.mapping=symbols.to.ensg.mapping, sample.id.column.name=sample.id.column.name, cell_type_column=cell_type_column, assay=assay, prepend_1 = prepend_1)
  }
}


create_features_files_bulk <- function(seurat_object,  output_loc, symbols.to.ensg=F, symbols.to.ensg.mapping = "genes.tsv", sample.id.column.name="assignment", assay = "SCT", prepend_1 = T){
  # set the assay
  DefaultAssay(seurat_object) <- assay
  # get the individuals in the object
  individuals <- unique(seurat_object@meta.data[,sample.id.column.name])
  # remove any na individuals (sometimes caused if using demux assignments)
  individuals <- individuals[!is.na(individuals)]
  
  # create a matrix where we will store the mean expressions
  mean_expression_matrix <- matrix(nrow=nrow(seurat_object), ncol = length(individuals),
                                   dimnames = list(rownames(seurat_object), individuals))
  # go through the individuals
  for (individual in individuals) {
    # take into account that there might be zero expression
    if (sum(is.na(seurat_object@meta.data[,sample.id.column.name]) | seurat_object@meta.data[,sample.id.column.name] == individual) == 0) {
      mean_expression_matrix[,individual] <- 0
    }
    else {
      # grab the cells of this cell type and this individual
      cells_individual <- seurat_object[,!(is.na(seurat_object@meta.data[,sample.id.column.name])) & seurat_object@meta.data[,sample.id.column.name] == individual]
      # calculate the mean expression for each gene
      mean_expression_matrix[,individual] <- rowMeans(cells_individual)
    }
  }
  # remap symbols to ensg numbers if necessary
  if (symbols.to.ensg) {
    genes <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
    genes$V2 <- gsub("_", "-", make.unique(genes$V2))
    rownames(mean_expression_matrix) <- genes[match(rownames(mean_expression_matrix), genes$V2),"V1"]
  }
  # prepend the '1_' to the id
  if(prepend_1){
    colnames(mean_expression_matrix) <- paste0("1_", colnames(mean_expression_matrix))
  }
  # write our table of means for this cell type
  write.table(mean_expression_matrix,
              file = paste0(output_loc, "bulk", "_expression", ".tsv"),
              quote = F, sep = "\t", col.names = NA)
}


create_feature_files_per_condition_bulk <- function(seurat_object, output_loc, condition.column.name = "timepoint", symbols.to.ensg=F, symbols.to.ensg.mapping = "genes.tsv", sample.id.column.name="assignment", assay = "SCT", cell_type_column = "cell_type"){
  # grab the conditions
  conditions <- unique(seurat_object@meta.data[[condition.column.name]])
  # remove any na conditions (sometimes caused if using demux assignments)
  conditions <- conditions[!is.na(conditions)]
  # go through the conditions
  for(condition in conditions){
    # there should be a subdirectory per condition
    output_loc_condition <- paste0(output_loc, condition,"/")
    # there might be a prepended 'X' that we might want to remove
    if(startsWith(condition, "X")){
      output_loc_condition <- paste0(output_loc, substr(condition, 2, nchar(condition)),"/")
    }
    # grab only the cells with this condition
    cells_condition <- seurat_object[,seurat_object@meta.data[condition.column.name] == condition]
    # do the feature file creation for this condition
    create_features_files_bulk(cells_condition, output_loc_condition, symbols.to.ensg=symbols.to.ensg, symbols.to.ensg.mapping = symbols.to.ensg.mapping, sample.id.column.name=sample.id.column.name, assay=assay)
  }
}
# where to write, tmp04 for Calculon, scr01 for Bender
write_partition <- 'tmp04'

# locations of objects
object_dir <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
v2_object_loc <- paste(object_dir, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds', sep = '')
#v2_object_loc <- paste(object_dir, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')
v3_object_loc <- paste(object_dir, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201029.rds', sep = '')

# eqtl output dirs
features_dir <- paste('/groups/umcg-bios/', write_partition, '/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/', sep = '')
v2_features_demux_dir <- paste(features_dir, 'v2_sct_mqc_demux_lores_20201029/', sep = '')
v3_features_demux_dir <- paste(features_dir, 'v3_sct_mqc_demux_lores_20201106/', sep = '')
v2_features_soup_dir <- paste(features_dir, 'v2_sct_mqc_soup_lores_20201029/', sep = '')
v3_features_soup_dir <- paste(features_dir, 'v3_sct_mqc_soup_lores_20201106/', sep = '')

# for the inhouse eQTL-mapping pipeline, we currently have ENSG numbers instead of gene numbers
gene_to_ens_mapping <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/resources/features_v3.tsv"

# read object
v2 <- readRDS(v2_object_loc)
v2 <- v2[, !is.na(v2@meta.data$timepoint)]
v2 <- v2[, !is.na(v2@meta.data$assignment)]
v2 <- v2[, !is.na(v2@meta.data$cell_type_lowerres)]
# write the bulk data with demux identities
create_feature_files_per_condition_bulk(seurat_object=v2, output_loc=v2_features_demux_dir, condition.column.name = "timepoint", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment", assay = "SCT", cell_type_column = "cell_type_lowerres")
# create the per-ct data with demux identities
create_feature_files_per_condition(seurat_object=v2, output_loc=v2_features_demux_dir, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint", sample.id.column.name="assignment", cell_type_column = "cell_type_lowerres", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = T)
# write the bulk data with the soup identities
create_feature_files_per_condition_bulk(seurat_object=v2, output_loc=v2_features_soup_dir, condition.column.name = "timepoint.soup", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment.soup", assay = "SCT", cell_type_column = "cell_type_lowerres")
# create the per-ct data with soup identities
create_feature_files_per_condition(seurat_object=v2, output_loc=v2_features_demux_dir, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.soup", sample.id.column.name="assignment.soup", cell_type_column = "cell_type_lowerres", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = T)
# clear up memory
rm(v2)

# read object
v3 <- readRDS(v3_object_loc)
v3 <- v3[, !is.na(v3@meta.data$timepoint)]
v3 <- v3[, !is.na(v3@meta.data$assignment)]
v3 <- v3[, !is.na(v3@meta.data$cell_type_lowerres)]
# write the bulk data with demux identities
create_feature_files_per_condition_bulk(seurat_object=v3, output_loc=v3_features_demux_dir, condition.column.name = "timepoint", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment", assay = "SCT", cell_type_column = "cell_type_lowerres")
# create the per-ct data with demux identities
create_feature_files_per_condition(seurat_object=v3, output_loc=v3_features_demux_dir, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint", sample.id.column.name="assignment", cell_type_column = "cell_type_lowerres", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = T)
# write the bulk data with the soup identities
create_feature_files_per_condition_bulk(seurat_object=v3, output_loc=v3_features_soup_dir, condition.column.name = "timepoint.soup", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment.soup", assay = "SCT", cell_type_column = "cell_type_lowerres")
# create the per-ct data with soup identities
create_feature_files_per_condition(seurat_object=v3, output_loc=v3_features_demux_dir, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.soup", sample.id.column.name="assignment.soup", cell_type_column = "cell_type_lowerres", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = T)
# clear up memory
rm(v3)
