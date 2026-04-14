#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: 1m_differential_expression_limma_parameterised.R
# Function: perform differential gene expression analysis using Limma Dream
############################################################################################################################


####################
# libraries        #
####################

# read the object
library(Seurat)
# plotting
library(ggplot2)
# limma DE dependencies
library(variancePartition)
library(edgeR)
library(BiocParallel)
# convert count matrices
library(Matrix)
library(Matrix.utils) # NOT IN CONTAINER! als grr or remotes::install_github("cvarrichio/Matrix.utils")
library(optparse) # NOT IN CONTAINER


####################
# Functions        #
####################

#' perform pseudobulk limma
#' 
#' @param seurat_object The Seurat object to add the sample assignment to
#' @param output_loc the location to write the output tables to
#' @param condition_combinations named list of vectors which to pairwise compare e.g. list('condition_final' = c('24hCA', 'UT'))
#' @param aggregates the columns to aggregate on, so pseudobulk per donor and inflammation is c('donor', 'inflammation')
#' @param fixed_effects vector of fixed effects to include in the model
#' @param random_effects vector or random effects to include in the model
#' @param minimal_cells the minimal number of cells that needs to be present for a sample to not be excluded
#' @param min_numi the minimal number of UMIs a cell must have to be used for the pseudobulk
#' @param minimal_complexity the minimal number of cells transcripts a pseudobulk needs to be based on for a sample to not be excluded
#' @param verbose whether to print progress messages
#' @returns 0 if successful
#' dream_pairwise(pbmc, './bulk_test/', list('inflammation'=c('AI','NI')))
dream_pairwise <- function(seurat_object, output_loc, condition_combinations, aggregates=c('assignment', 'timepoint'), fixed_effects=c('timepoint'), random_effects=c('assignment'), minimal_cells=0, min_numi=200, minimal_complexity=5000, verbose=T){
  # grab the countmatrix
  countMatrix <- NULL
  # subset the object based on the minimal number of umis if requested
  if (min_numi > 0) {
    seurat_object <- seurat_object[, seurat_object@meta.data[['nFeature_SCT']] >= min_numi]
  }
  # depending on the version
  if (grepl('^3|4.0', seurat_object@version)) {
    countMatrix <- seurat_object@assays$SCT@counts
  }
  else if (grepl('^4.9|5', seurat_object@version)) {
    countMatrix <- seurat_object@assays$SCT@counts
  }
  # get the metadata
  metadata <- seurat_object@meta.data
  # create the groups to aggregate on, here it's on sample and inflammation status usually
  groups <- metadata[, unique(c(aggregates, fixed_effects, random_effects))]
  # create an aggregated counts matrix
  aggregate_countMatrix <- t(aggregate.Matrix(t(countMatrix), groupings = groups, fun = 'sum'))
  # create aggregated metadata
  aggregate_metadata <- unique(metadata[, unique(c(aggregates, fixed_effects, random_effects))])
  # set the rownames of the aggregate metadata
  rownames_to_set_agg_metadata <- aggregate_metadata[[unique(c(aggregates, fixed_effects, random_effects))[1]]]
  for(i in 2:length(unique(c(aggregates, fixed_effects, random_effects)))){
    rownames_to_set_agg_metadata <- paste(rownames_to_set_agg_metadata, aggregate_metadata[[unique(c(aggregates, fixed_effects, random_effects))[i]]], sep='_')
  }
  rownames(aggregate_metadata) <- rownames_to_set_agg_metadata
  # set in the same order as the count matrix
  aggregate_metadata <- aggregate_metadata[colnames(aggregate_countMatrix), ]
  
  # next get the cell numbers for each observation, these will follow the order of the original aggregated metadata
  cell_numbers <- get_nr_cells_aggregate_combination(aggregate_metadata[, unique(c(aggregates, fixed_effects, random_effects))], metadata)
  
  # filter by the number of cells if requested
  if (minimal_cells > 0) {
    # get the indices of where the cell numbers are above this
    indices_above_threshold <- which(cell_numbers[['nr']] >= minimal_cells)
    # report how many
    if (verbose) {
      message(paste('of', as.character(nrow(cell_numbers)), 'entries, ', length(indices_above_threshold), 'contained more cells than the', as.character(minimal_cells), 'threshold'))
    }
    # do the actual filtering
    cell_numbers <- cell_numbers[indices_above_threshold, ]
    aggregate_countMatrix <- aggregate_countMatrix[, indices_above_threshold]
    aggregate_metadata <- aggregate_metadata[indices_above_threshold, ]
  }
  # filter by the complexity if requested
  if (minimal_complexity > 0) {
    # calculate the complexity first
    complexity <- data.frame((colSums(aggregate_countMatrix)))
    # set better column names
    colnames(complexity) <- c('complexity')
    # get the indices of where the complexity is above the threshold
    indices_above_complexity <- which(complexity[['complexity']] >= minimal_complexity)
    # report how many
    if (verbose) {
      message(paste('of', as.character(nrow(complexity)), 'entries, ', length(indices_above_complexity), 'contained more UMIs than the', as.character(minimal_complexity), 'threshold'))
    }
    # and filter
    cell_numbers <- cell_numbers[indices_above_complexity, ]
    aggregate_countMatrix <- aggregate_countMatrix[, indices_above_complexity]
    aggregate_metadata <- aggregate_metadata[indices_above_complexity, ]
  }
  
  # filter genes by number of counts
  isexpr = rowSums(cpm(aggregate_countMatrix)>0.1) >= 5
  
  # Standard usage of limma/voom
  geneExpr = DGEList( aggregate_countMatrix[isexpr,] )
  geneExpr = calcNormFactors( geneExpr )
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  register(param)
  
  if(verbose){
    print('aggregated metadata head:')
    print(head(aggregate_metadata))
    print('aggregated counts head:')
    print(head(geneExpr[['counts']]))
    print('cell numbers head')
    print(head(cell_numbers[order(cell_numbers[['nr']]), ]))
  }
  
  # paste together the model
  model_formula <- '~ 0 '
  for(fixed_effect in fixed_effects){
    model_formula <- paste(model_formula, fixed_effect, sep = ' + ')
  }
  for(random_effect in random_effects){
    model_formula <- paste(model_formula, ' + (1|', random_effect, ')', sep = '')
  }
  
  if(verbose){
    print(paste('formula:', model_formula))
  }
  
  # and turn into a formula
  form <- as.formula(model_formula)
  tryCatch(
    {
      # estimate weights using linear mixed model of dream
      vobjDream = voomWithDreamWeights( counts = geneExpr, formula = form, data = aggregate_metadata, weights = cell_numbers[['nr']] ) # the cell numbers are in the same order as the metadata, and as such can be passed like this
      
      # do each combination
      for(combination_name in names(condition_combinations)){
        # grab the combination
        combination <- condition_combinations[[combination_name]]
        
        if(verbose){
          print(paste('doing combination:', combination_name, sep = ' ', collapse = ' '))
        }
        
        tryCatch(
          {
            # define and then cbind contrasts
            L = getContrast( vobjDream, form, aggregate_metadata, c(paste(combination_name, combination[1], sep=''), paste(combination_name, combination[2], sep='')))
            
            # fit contrast
            fit = dream( vobjDream, form, aggregate_metadata, L)
            
            # grab the exact fit
            limma_result <- topTable(fit, coef='L1', number=length(fit$F.p.value))
            
            
            # add bonferroni adjustment
            limma_result[['p.bonferroni']] <- p.adjust(limma_result[['P.Value']])
            
            # add some statistics
            result_stats_list <- list()
            # check each condition
            result_stats_list[['combination']] <- data.frame(combination=rep(paste(combination_name, paste(combination, collapse='-'), sep = '.'), times = nrow(limma_result)))
            for (condition in combination) {
              # get the cell numbers for the condition
              cell_numbers_condition <- cell_numbers[cell_numbers[[combination_name]] == condition, ]
              result_stats_list[[paste('nsample', condition, sep = '_')]] <- data.frame(nsample=rep(nrow(cell_numbers_condition), times = nrow(limma_result)))
              result_stats_list[[paste('ncell', condition, sep = '_')]] <- data.frame(ncells=rep(
                paste(as.character(min(cell_numbers_condition[['nr']])),
                      as.character(quantile(cell_numbers_condition[['nr']])[['25%']]),
                      as.character(quantile(cell_numbers_condition[['nr']])[['50%']]),
                      as.character(quantile(cell_numbers_condition[['nr']])[['75%']]),
                      as.character(max(cell_numbers_condition[['nr']])),
                      sep = ';'
                ), times = nrow(limma_result)))
            }
            # merge all
            result_stats <- do.call('cbind', result_stats_list)
            colnames(result_stats) <- names(result_stats_list)
            
            # add combination as first column
            limma_result <- cbind(result_stats, limma_result)
            
            # finally also add the feature as an explicit column
            limma_result <- cbind(data.frame(feature = rownames(limma_result)), limma_result)
            
            # set an output location
            limma_output_loc <- gzfile(paste(output_loc, combination_name, '.tsv.gz', sep = ''))
            
            # also write the model we used
            limma_formula_loc <- paste(output_loc, combination_name, '.formula', sep = '')
            
            if(verbose){
              print(paste('writing result', paste(output_loc, combination_name, '.tsv.gz', sep = '')))
            }
            
            # write the result
            write.table(limma_result, limma_output_loc, sep = '\t', row.names = F)
            # and the formula
            write.table(model_formula, limma_formula_loc, row.names = F, col.names = F)
          }, error=function(cond) {
            print(paste('analysis failed in', combination))
            message(cond)
          }
        )
      }
    }, error=function(cond) {
      print(paste('model build failed'))
      message(cond)
    }
  )
  return(0)
}


#' perform pseudobulk limma per celltype
#' 
#' @param seurat_object The Seurat object to add the sample assignment to
#' @param output_loc the location to write the output tables to
#' @param condition_combinations named list of vectors which to pairwise compare e.g. list('condition_final' = c('24hCA', 'UT'))
#' @param celltype_column the column in the metadata that contains the cell type to perform limma on
#' @param cell_type_to_use the cell types to perform limma on. If NULL, all cell types will be used
#' @param aggregates the columns to aggregate on, so pseudobulk per donor and inflammation is c('donor', 'inflammation')
#' @param fixed_effects vector of fixed effects to include in the model
#' @param random_effects vector or random effects to include in the model
#' @param minimal_cells the minimal number of cells that needs to be present for a sample to not be excluded
#' @param min_numi the minimal number of UMIs a cell must have to be used for the pseudobulk
#' @param minimal_complexity the minimal number of cells transcripts a pseudobulk needs to be based on for a sample to not be excluded
#' @param verbose whether to print progress messages
#' @returns 0 if successful
#' do_limma_dream_pairwise_per_celltype(pbmc, './ct_test/')
do_limma_dream_pairwise_per_celltype <- function(seurat_object, output_loc, condition_combinations=list('timepoint' =  c('X24hCA', 'UT')), celltype_column='cell_type_final', cell_types_to_use=NULL, aggregates=c('assignment', 'timepoint'), fixed_effects=c('timepoint'), random_effects=c('assignment'), minimal_cells=0, min_numi=200, minimal_complexity=5000, verbose=T){
  # use the cell types supplied, or all if none are supplied
  cell_types <- cell_types_to_use
  if(is.null(cell_types_to_use)){
    cell_types <- unique(seurat_object@meta.data[[celltype_column]])
  }
  # remove NA ones
  cell_types <- cell_types[!is.na(cell_types)]
  # check each cell type
  for(cell_type in cell_types){
    if (verbose) {
      print(paste('doing', cell_type))
    }
    # get a POSIX-safe celltype name
    cell_type_safe <- gsub(' |/', '_', cell_type)
    cell_type_safe <- gsub('-', '_negative', cell_type_safe)
    cell_type_safe <- gsub('\\+', '_positive', cell_type_safe)
    cell_type_safe <- gsub('\\)', '', cell_type_safe)
    cell_type_safe <- gsub('\\(', '', cell_type_safe)
    # make that into the output prepend
    output_loc_celltype <- paste(output_loc, cell_type_safe, '_', sep = '')
    # subset to the cell type
    seurat_object_celltype <- seurat_object[, seurat_object@meta.data[[celltype_column]] == cell_type]
    # perform the analysis
    dream_pairwise(seurat_object = seurat_object_celltype, output_loc = output_loc_celltype, condition_combinations = condition_combinations, aggregates = aggregates, fixed_effects = fixed_effects, random_effects = random_effects, minimal_cells = minimal_cells, min_numi = min_numi, minimal_complexity = minimal_complexity, verbose = verbose)
  }
  return(0)
}

#' remove illegal characters from the cell type annotation
#' 
#' @param cell_types vector of cell types
#' @returns vector of celltypes with safe name
#' make_celltypes_safe(c('CD4+T', 'NK-dim'))
make_celltypes_safe <- function(cell_types){
  # get a safe file name
  cell_type_safes <- gsub(' |/', '_', cell_types)
  cell_type_safes <- gsub('-', '_negative', cell_type_safes)
  cell_type_safes <- gsub('\\+', '_positive', cell_type_safes)
  cell_type_safes <- gsub('\\)', '', cell_type_safes)
  cell_type_safes <- gsub('\\(', '', cell_type_safes)
  return(cell_type_safes)
}

#' remove illegal characters from the cell type annotation
#' 
#' @param seurat_object the seurat object to add a lower classification of cell type to
#' @param reclassification_mapping a table that has the original cell type, and what it should be changed to
#' @param mapping_original_column the column in the mapping table with the original cell types
#' @param mapping_reclass_column the column in the mapping table with what to change the cell type to
#' @param metadata_original_column the column in the metadata to get the original cell types from
#' @param metadata_reclassification_column the metadata column to add with the new mapping
#' @returns vector of celltypes with safe name
#' pmbc <- add_lower_classification(pmbc, mapping_table, 'high_res_ct', 'low_res_ct', 'cell.type', 'cell.type.low')
add_lower_classification <- function(seurat_object, reclassification_mapping, mapping_original_column, mapping_reclass_column, metadata_original_column, metadata_reclassification_column){
  # add the new column
  seurat_object@meta.data[[metadata_reclassification_column]] <- NA
  # get each cell type in the data
  metadata_original_cts <- unique(seurat_object@meta.data[[metadata_original_column]])
  # and the originals in the mapping
  reclassification_original_cts <- unique(reclassification_mapping[[mapping_original_column]])
  # we can only map what is present in both
  originals_both <- intersect(metadata_original_cts, reclassification_original_cts)
  # check what is missing
  only_metadata <- setdiff(metadata_original_cts, reclassification_original_cts)
  only_mapping <- setdiff(reclassification_original_cts, metadata_original_cts)
  # warn what is missing
  if(length(only_metadata) > 0){
    print('some celltypes only in metadata')
    print(only_metadata)
  }
  if(length(only_mapping) > 0){
    print('some celltypes only in remapping ')
    print(only_mapping)
  }
  # check each cell type
  for(celltype_original in originals_both){
    # get the appropriate remapping
    celltype_remapped <- reclassification_mapping[reclassification_mapping[[mapping_original_column]] == celltype_original, mapping_reclass_column]
    # now remap in the metadata
    seurat_object@meta.data[seurat_object@meta.data[[metadata_original_column]] == celltype_original, metadata_reclassification_column] <- celltype_remapped
  }
  return(seurat_object)
}

#' create a column which combines the values of multiple other columns
#' 
#' @param dataframe dataframe that has the original column, to which to add the aggregate column
#' @param columns the columns to combine into a new column
#' @returns the original dataframe, with a new column 'all_aggregates', which is the combination of the supplied columns
#' new_df <- combine_columns(df, c('V1', 'V3'))
combine_columns <- function(dataframe, columns) {
  # check each column
  for (column in columns) {
    # check if we were already aggregating
    if ('all_aggregates' %in% colnames(dataframe)) {
      # add the new column
      dataframe[['all_aggregates']] <- paste(dataframe[['all_aggregates']], dataframe[[column]], sep = '_')
    }
    # otherwise we need to start our aggregation
    else{
      dataframe[['all_aggregates']] <- dataframe[[column]]
    }
  }
  return(dataframe)
}

#' get the number of cells describing a combination of values from the seurat metadata
#' 
#' @param aggregate_df the dataframe with the combinations that were aggregated over
#' @param seurat_metadata the metadata of a seurat object from which to get cell numbers for the combinations
#' @returns the original dataframe, with a 'nr' column, denoting the number of cells for that combination
#' cell_numbers <- get_nr_cells_aggregate_combination(aggretates, pbmc_meta.data)
get_nr_cells_aggregate_combination <- function(aggregate_df, seurat_metadata) {
  # grab the columns we aggregated on
  aggregate_columns <- colnames(aggregate_df)
  # now use table to get the number of entries of these aggregates in the metadata
  aggregate_numbers <- data.frame(table(seurat_metadata[, aggregate_columns]))
  # add a new column for the combination of the aggregates
  aggregate_df <- combine_columns(aggregate_df, aggregate_columns)
  # in our numbers as well
  aggregate_numbers <- combine_columns(aggregate_numbers, aggregate_columns)
  # now join the frequencies from the aggragate numbers onto the aggregate columns
  aggregate_df[['nr']] <- aggregate_numbers[match(aggregate_df[['all_aggregates']], aggregate_numbers[['all_aggregates']]), 'Freq']
  # now we are sure that the cell numbers are in the same order as the aggregated metadata
  # remove the column we created
  aggregate_df[['all_aggregates']] <- NULL
  # and return the result
  return(aggregate_df)
}


do_debug <- function() {
  # fill the opt
  opt <- list()
  opt[['out']] <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation_1m/'
  opt[['file']] <- paste('/groups/umcg-franke-scrna/tmp03/releases/wijst-2020-hg19/v1/seurat/1M_both_', 'monocyte', '_UT_24hCA_seuratv5.rds', sep = '')
  opt[['cell_type_column']] <- 'cell_type_lowerres'
  opt[['min_cells']] <- 10
  opt[['min_umi']] <- 200
  opt[['min_complexity']] <- 2000
  
  # location of the DE
  limma_output_loc <- opt[['out']]
  # locations of objects
  seurat_object_object_loc <- opt[['file']]
  # celltype column
  celltype_column <- opt[['cell_type_column']]
  # minimal number of cells
  min_cells <- opt[['min_cells']]
  # minimal UMIs
  min_cell_umis <- opt[['min_umi']]
  # minimal number of UMIs of pseudobulk
  min_pseudo_umis <- opt[['min_complexity']]
  
  # read the Seurat object
  seurat_object <- readRDS(seurat_object_object_loc)
  
  # add the day to the metadata
  seurat_object@meta.data[['day']] <- gsub('_lane\\d+', '', seurat_object@meta.data[['lane']])
  
  # filter where we don't have the sex
  seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[['sex']])]
  # or the age
  seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[['age']])]
  # or the inflammation status
  seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[['timepoint']]) & seurat_object@meta.data[['timepoint']] != 'unknown']
  
  # do the bulk analysis
  do_limma_dream_pairwise_per_celltype(seurat_object, 
                                       output_loc = limma_output_loc, 
                                       celltype_column = celltype_column, 
                                       aggregates = c('timepoint', 'lane', 'assignment'), 
                                       fixed_effects = c('timepoint', 'age', 'sex'), 
                                       random_effects = c('assignment', 'lane'), 
                                       condition_combinations=list('timepoint' =  c('X24hCA', 'UT')), 
                                       minimal_cells = min_cells,
                                       min_numi = min_cell_umis, 
                                       minimal_complexity = min_pseudo_umis)
}


####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="seurat object input", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-c", "--cell_type_column"), type="character", default='cell_type_safe',
              help="column describing the celltype in the Seurat metadata [default= %default]", metavar="character"),
  make_option(c("-m", "--min_cells"), type="numeric", default=0,
              help="minimal number of cells required for a pseudobulk to consider it in the analysis [default= %default]", metavar="numeric"),
  make_option(c("-u", "--min_umi"), type="numeric", default=200,
              help="minimal UMIs for a cell to keep it for pseudobulking [default= %default]", metavar="numeric"),
  make_option(c("-l", "--min_complexity"), type="numeric", default=0,
              help="minimal number of UMIs for a pseudobulk to consider it in the analysis [default= %default]", metavar="numeric")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# we need a Seurat object
if (is.null(opt[['file']])){
  print_help(opt_parser)
  stop('no seurat object specified', call.=FALSE)
}

# location of the DE
limma_output_loc <- opt[['out']]
# locations of objects
seurat_object_object_loc <- opt[['file']]
# celltype column
celltype_column <- opt[['cell_type_column']]
# minimal number of cells
min_cells <- opt[['min_cells']]
# minimal UMIs
min_cell_umis <- opt[['min_umi']]
# minimal number of UMIs of pseudobulk
min_pseudo_umis <- opt[['min_complexity']]

# read the Seurat object
seurat_object <- readRDS(seurat_object_object_loc)

# add the day to the metadata
seurat_object@meta.data[['day']] <- gsub('_lane\\d+', '', seurat_object@meta.data[['lane']])

# filter where we don't have the sex
seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[['sex']])]
# or the age
seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[['age']])]
# or the inflammation status
seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[['timepoint']]) & seurat_object@meta.data[['timepoint']] != 'unknown']

# do the bulk analysis
do_limma_dream_pairwise_per_celltype(seurat_object, 
                                     output_loc = limma_output_loc, 
                                     celltype_column = celltype_column, 
                                     aggregates = c('timepoint', 'chem', 'lane', 'assignment'), 
                                     fixed_effects = c('timepoint', 'chem', 'age', 'sex'), 
                                     random_effects = c('assignment', 'lane'),
                                     condition_combinations=list('timepoint' =  c('X24hCA', 'UT')), 
                                     minimal_cells = min_cells,
                                     min_numi = min_cell_umis, 
                                     minimal_complexity = min_pseudo_umis)
