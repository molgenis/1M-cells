library(Seurat)

get_relevant_covariates <- function(metadata, timepoint_column='timepoint', cell_type_column='cell_type', assignment_column='assignment', timepoints=NULL, cell_types=NULL, assignments=NULL, doRNA=T, doSCT=T, verbose=F){
  # check which conditions we have
  conditions_to_use <- unique(metadata[[timepoint_column]])
  # subset if requested
  if(!is.null(timepoints)){
    conditions_to_use <- intersect(conditions_to_use, timepoints)
  }
  # get rid of NAs
  conditions_to_use <- setdiff(conditions_to_use, c(NA))
  # save the tables somewhere
  per_condition <- list()
  # check each condition
  for(condition in conditions_to_use){
    if(verbose){
      print(paste('from', as.character(nrow(metadata)), 'subsetting'))
      print(paste('doing condition', condition))
    }
    # subset to the specific condition
    cells_condition <- metadata[!is.na(metadata[[timepoint_column]]) & metadata[[timepoint_column]] == condition, ]
    
    # next, check the cell types
    cell_types_to_use <- unique(cells_condition[[cell_type_column]])
    # again, subset if requested
    if(!is.null(cell_types)){
      cell_types_to_use <- intersect(cell_types_to_use, cell_types)
    }
    # get rid of NAs
    cell_types_to_use <- setdiff(cell_types_to_use, c(NA))
    per_condition_per_celltype <- list()
    # check each cell type
    for(cell_type in cell_types_to_use){
      if(verbose){
        print(paste('from', as.character(nrow(cells_condition)), 'subsetting'))
        print(paste('doing cell type', cell_type))
      }
      # subset to specific cell type
      cells_condition_celltype <- cells_condition[!is.na(cells_condition[[cell_type_column]]) & cells_condition[[cell_type_column]] == cell_type, ]
      # get the participants to use
      participants_to_use <- unique(cells_condition_celltype[[assignment_column]])
      # subset if requested
      if(!is.null(assignments)){
        participants_to_use <- intersect(participants_to_use, assignments)
      }
      # get rid of NAs
      participants_to_use <- setdiff(participants_to_use, c(NA))
      # we will do these covariates
      covar_columns <- c()
      if(doRNA){
        covar_columns <- c(
        'mean_nCount_RNA',
        'total_nCount_RNA',
        'mean_nFeature_RNA',
        'total_nFeature_RNA',
        'mean_feature_per_cell_RNA'
        )
      }
      if(doSCT){
        covar_columns <- c(covar_columns,
        'mean_nCount_SCT',
        'total_nCount_SCT',
        'mean_nFeature_SCT',
        'total_nFeature_SCT',
        'mean_feature_per_cell_SCT'
        )
      }
      # the number of cells is always present
      covar_columns <- c('nCell', covar_columns)
        
      # will create an empty matrix to fill
      covariate_data <- matrix(, nrow=length(participants_to_use), ncol=length(covar_columns), dimnames = list(participants_to_use, covar_columns))
      
      # check each participant
      for(participant in participants_to_use){
        if(verbose){
          print(paste('doing participant', participant))
        }
        # subset to that participant
        cells_condition_celltype_participant <- cells_condition_celltype[!is.na(cells_condition_celltype[[assignment_column]]) & cells_condition_celltype[[assignment_column]] == participant, ]
        
        # grab the number of cells
        nCell <- nrow(cells_condition_celltype_participant)
        # that is a metric we always want
        covariate_data[participant, 'nCell'] <- nCell
        if(doRNA){
          # grab the relevant data
          total_nCount_RNA <- sum(cells_condition_celltype_participant$nCount_RNA)
          mean_Count_RNA <- mean(cells_condition_celltype_participant$nCount_RNA)
          total_nFeature_RNA <- sum(cells_condition_celltype_participant$nFeature_RNA)
          mean_nFeature_RNA <- mean(cells_condition_celltype_participant$nFeature_RNA)
          mean_feature_per_cell_RNA <- mean(cells_condition_celltype_participant$nFeature_RNA / nCell)
          # put into table
          covariate_data[participant, c('mean_nCount_RNA',
                                        'total_nCount_RNA',
                                        'mean_nFeature_RNA',
                                        'total_nFeature_RNA',
                                        'mean_feature_per_cell_RNA')] <- c(
                                          total_nCount_RNA,
                                          mean_Count_RNA,
                                          total_nFeature_RNA,
                                          mean_nFeature_RNA,
                                          mean_feature_per_cell_RNA
                                        )
          
        }
        if(doSCT){
          # grab the relevant data
          total_nCount_SCT <- sum(cells_condition_celltype_participant$nCount_SCT)
          mean_Count_SCT <- mean(cells_condition_celltype_participant$nCount_SCT)
          total_nFeature_SCT <- sum(cells_condition_celltype_participant$nFeature_SCT)
          mean_nFeature_SCT <- mean(cells_condition_celltype_participant$nFeature_SCT)
          mean_feature_per_cell_SCT <- mean(cells_condition_celltype_participant$nFeature_SCT / nCell)
          # put into table
          covariate_data[participant, c('mean_nCount_SCT',
                                        'total_nCount_SCT',
                                        'mean_nFeature_SCT',
                                        'total_nFeature_SCT',
                                        'mean_feature_per_cell_SCT')] <- c(
                                          total_nCount_SCT,
                                          mean_Count_SCT,
                                          total_nFeature_SCT,
                                          mean_nFeature_SCT,
                                          mean_feature_per_cell_SCT
                                        )
        }
      }
      # add to list
      per_condition_per_celltype[[cell_type]] <- covariate_data
    }
    # and to list of lists
    per_condition[[condition]] <- per_condition_per_celltype
  }
  return(per_condition)
}

# read the ng2018 file
ng2018_object_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/clustering/pilot3_seurat3_200420_sct_azimuth.rds'
ng2018 <- readRDS(ng2018_object_loc)
# add dummy UT column
ng2018@meta.data$timepoint <- 'UT'
# grab the covariates into lists
ng2018_covars <- get_relevant_covariates(ng2018@meta.data, 'timepoint', 'cell_type', 'sample', doRNA=F)

# create column for lower res ct annotation
ng2018@meta.data$cell_type_lowerres <- NA
# get a character representation of the columns because I don't like Factors
ng2018@meta.data$cell_type_char <- as.character(ng2018@meta.data$cell_type)

# create groups
plasma <- 'Plasma'
b <- c('B')
cd4t <- c('CD4+ T')
cd8t <- c('CD8+ T')
dc <- c('pDC', 'mDC')
monocyte <- c('cMonocyte', 'ncMonocyte')
nk <- c('CD56(dim) NK', 'CD56(bright) NK')
megakaryocyte <- c('Megakaryocyte')
# use groups to get lower resolution cell type
ng2018@meta.data[ng2018@meta.data$cell_type_char %in% plasma, 'cell_type_lowerres'] <- 'plasma'
ng2018@meta.data[ng2018@meta.data$cell_type_char %in% b, 'cell_type_lowerres'] <- 'B'
ng2018@meta.data[ng2018@meta.data$cell_type_char %in% cd4t, 'cell_type_lowerres'] <- 'CD4T'
ng2018@meta.data[ng2018@meta.data$cell_type_char %in% cd8t, 'cell_type_lowerres'] <- 'CD8T'
ng2018@meta.data[ng2018@meta.data$cell_type_char %in% dc, 'cell_type_lowerres'] <- 'DC'
ng2018@meta.data[ng2018@meta.data$cell_type_char %in% monocyte, 'cell_type_lowerres'] <- 'monocyte'
ng2018@meta.data[ng2018@meta.data$cell_type_char %in% nk, 'cell_type_lowerres'] <- 'NK'
ng2018@meta.data[ng2018@meta.data$cell_type_char %in% megakaryocyte, 'cell_type_lowerres'] <- 'megakaryocyte'
# remove our helper column
ng2018@meta.data$cell_type_char <- NULL


ng2018_lr_covars <- get_relevant_covariates(ng2018@meta.data, 'timepoint', 'cell_type_lowerres', 'sample', doRNA=F, verbose = T) # there is no condition, so orig.ident is used as a dummy variable


