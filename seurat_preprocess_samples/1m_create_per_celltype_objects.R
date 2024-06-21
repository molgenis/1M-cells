#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: 1m_create_per_celltype_objects.R
# Function: 
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)


####################
# Functions        #
####################


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 120 * 2000 * 2024^2)

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# location of the objects
seurat_object_loc_v2 <- '/groups/umcg-franke-scrna/tmp03/releases/wijst-2020-hg19/v1/seurat/1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds'
seurat_object_loc_v3 <- '/groups/umcg-franke-scrna/tmp03/releases/wijst-2020-hg19/v1/seurat/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds'

# load object
m1_v2 <- readRDS(seurat_object_loc_v2)
m1_v3 <- readRDS(seurat_object_loc_v3)
# merge
m1 <- merge(m1_v2, m1_v3)
# clear memory
rm(m1_v2)
rm(m1_v3)

# now go through each cell type without taking UT only
for (cell_type in unique(m1@meta.data$cell_type_lowerres)) {
  # check for NA
  if (!is.na(cell_type)) {
    # subset to this celltype and condition
    m1_ct <- m1[, !is.na(m1@meta.data$cell_type_lowerres) &
                        !is.na(m1@meta.data$timepoint_final) &
                        m1@meta.data$timepoint_final %in% c('UT', 'X24hCA')&
                        m1@meta.data$cell_type_lowerres == cell_type]
    # write this file
    saveRDS(m1_ct,
            paste('/groups/umcg-franke-scrna/tmp03/releases/wijst-2020-hg19/v1/seurat/', cell_type, 'UT_24hCA.rds', sep = ''))
    
  }
}
