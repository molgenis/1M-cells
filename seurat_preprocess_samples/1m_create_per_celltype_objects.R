#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: 1m_create_per_celltype_objects.R
# Function: create Seurat v5 object per celltype
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
# update
m1_v2 <- UpdateSeuratObject(m1_v2)
m1_v3 <- UpdateSeuratObject(m1_v3)# merge
m1 <- merge(m1_v2, m1_v3)
# clear memory
rm(m1_v2)
rm(m1_v3)

# read the age/sex file
age_sex_loc <- '/groups/umcg-franke-scrna/tmp03/releases/wijst-2020-hg19/v1/metadata/1M_exp_age_gender.tsv'
age_sex <- read.table(age_sex_loc, header = T, sep = '\t')
# add age/sex
m1@meta.data[['age']] <- age_sex[match(m1@meta.data$exp.id, age_sex$ExpNr), 'Age']
m1@meta.data[['sex']] <- age_sex[match(m1@meta.data$exp.id, age_sex$ExpNr), 'Gender']
# remove unneccessary column
m1@meta.data$cell_type_lowerres_old <- NULL
# save result
saveRDS(m1, '/groups/umcg-franke-scrna/tmp03/releases/wijst-2020-hg19/v1/seurat/1M_both_mediumQC_ctd_rnanormed_demuxids_20240621.rds')

# now go through each cell type without taking UT only
for (cell_type in unique(m1@meta.data$cell_type_lowerres)) {
  # check for NA
  if (!is.na(cell_type)) {
    # subset to this celltype and condition
    m1_ct <- m1[, !is.na(m1@meta.data$cell_type_lowerres) &
                  !is.na(m1@meta.data$timepoint) &
                  m1@meta.data$timepoint %in% c('UT', 'X24hCA')&
                  m1@meta.data$cell_type_lowerres == cell_type]
    # write this file
    saveRDS(m1_ct,
            paste('/groups/umcg-franke-scrna/tmp03/releases/wijst-2020-hg19/v1/seurat/1M_both_', cell_type, '_UT_24hCA_seuratv5.rds', sep = ''))
    
  }
}
