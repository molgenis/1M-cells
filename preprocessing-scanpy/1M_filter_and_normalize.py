###########################################################################################################################
# Authors: Dylan de Vries, Roy Oelen
# Name: 1M_filter_and_normalize.py
# Function: grab the h5ad file that had the read from cellranger and the assignments to IDs and timepoints, and do preprocessing
#
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################

import sys
import numpy as np
import os
import re
import pandas as pd
import scanpy as sc

###########################################################################################################################
#
# Functions
#
###########################################################################################################################


def preprocess_data(scanpy_object, perc_mito=8, min_cells=3, min_genes=200, max_genes=2500, use_old_doublet_detection=False):
    """
    preprocess the lanes that were loaded
    :param scanpy_object: the object to filer
    :param perc_mito: the max percent of mitochondrial RNA that may be present
    :param min_cells: the minimum number of cells there must be with expression for a gene to be used
    :param min_genes: the minimum number of genes that must be present to use a cell
    :param max_genes: the max number of genes that must be present in a cell with low singlet certainty
    :param use_old_doublet_detection: use the old style of taking the likelyhoods from demuxlet to check for doublets
    :return: the scanpy object with the filterig applied
    """
    scanpy_object.var_names_make_unique()
    # throw away genes and cells with low counts
    sc.pp.filter_cells(scanpy_object, min_genes=min_genes)
    sc.pp.filter_genes(scanpy_object, min_cells=min_cells)
    mito_genes = scanpy_object.var_names.str.startswith('MT-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    scanpy_object.obs['percent_mito'] = np.sum(
        scanpy_object[:, mito_genes].X, axis=1).A1 / np.sum(scanpy_object.X, axis=1).A1
    # add the total counts per cell as observations-annotation to scanpy_object
    scanpy_object.obs['n_counts'] = scanpy_object.X.sum(axis=1).A1
    # calculate QC metrics
    sc.pp.calculate_qc_metrics(scanpy_object)
    scanpy_object = scanpy_object[scanpy_object.obs['percent_mito'] <= perc_mito/100, :]
    # scanpy_object = scanpy_object[scanpy_object.obs['n_genes'] <= max_genes, :]
    # get possible doublets
    possible_doublets = []
    if use_old_doublet_detection:
        # the old method relies on the change of doublet vs singlet
        possible_doublets = get_possible_doublets_old_method(scanpy_object, max_genes)
    else:
        # the new method just takes the assignment from Souporcell
        possible_doublets = scanpy_object.obs['status'] == 'doublet'
    scanpy_object.obs['possible_doublets'] = possible_doublets
    # include based on singlet likelihood
    scanpy_object = scanpy_object[scanpy_object.obs['possible_doublets'] == False]
    # throw away genes and cells with low counts
    sc.pp.filter_cells(scanpy_object, min_genes=min_genes)
    sc.pp.filter_genes(scanpy_object, min_cells=min_cells)
    return scanpy_object


def get_possible_doublets_old_method(scanpy_object, max_genes):
    """
    get the possible doublets information added to a scanpy object
    :param scanpy_object: the object to get the doublet data to
    :param max_genes: the max number of genes for a cell with low singlet certainty
    :return: the doublet detection added in the order of the scanpy object
    """
    # init the array to store the possible doublets we may have missed
    doublet = []
    # go through the indices of the observations
    for i in range(len(scanpy_object.obs.index)):
        llk = float(scanpy_object.obs['LLK12_demux'][i]) - float(scanpy_object.obs['SNG_LLK1_demux'][i])
        # Seurat> pbmc <- subset(pbmc, subset = LLK12 - SNG.LLK1 < 25)
        # Seurat> pbmc <- subset(pbmc, subset = LLK12 - SNG.LLK1  < 0 | nFeature_RNA < 2000)
        if (llk < 25) and (llk < 0 | scanpy_object.obs['n_genes'][i] < max_genes):
            doublet.append(False)
        else:
            doublet.append(True)
    return doublet

###########################################################################################################################
#
# Main code
#
###########################################################################################################################

sc.settings.verbosity = 3           # verbosity: errors (0), warnings (1), info (2), hints (3)

# the output h5ad locations of each step
object_loc = "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/scanpy_preprocess_samples/objects/"
lanes_v2_assigned_loc = ''.join([object_loc, '1M_cells_scanpy_object_v2_assigned.h5ad'])
lanes_v3_assigned_loc = ''.join([object_loc, '1M_cells_scanpy_object_v3_assigned.h5ad'])
lanes_v2_preprocessed_loc = ''.join([object_loc, '1M_cells_scanpy_object_v2_pp.h5ad'])
lanes_v3_preprocessed_loc = ''.join([object_loc, '1M_cells_scanpy_object_v3_pp.h5ad'])
lanes_v2_scaled_loc = ''.join([object_loc, '1M_cells_scanpy_object_v2_scaled.h5ad'])
lanes_v3_scaled_loc = ''.join([object_loc, '1M_cells_scanpy_object_v3_scaled.h5ad'])

# load from the previous step (1M_read_and_assign)
v2_assigned = sc.read(lanes_v2_assigned_loc)
v3_assigned = sc.read(lanes_v3_assigned_loc)

# do preprocessing step
v2_processed = preprocess_data(v2_assigned, max_genes=2000)
v3_processed = preprocess_data(v3_assigned, max_genes=2000, perc_mito=15)

# write to files as backup
v2_processed.write(lanes_v2_preprocessed_loc)
v3_processed.write(lanes_v3_preprocessed_loc)

# remove zero-count genes, as that makes the regression not work
sc.pp.filter_genes(v2_processed, min_counts=1)
sc.pp.filter_genes(v3_processed, min_counts=1)

# do the normalization
#sc.pp.normalize_per_cell(v2_processed, counts_per_cell_after=1e4)
sc.pp.normalize_total(v2_processed, target_sum=1e4)
sc.pp.log1p(v2_processed)
#sc.pp.normalize_per_cell(v3_processed, counts_per_cell_after=1e4)
sc.pp.normalize_total(v3_processed, target_sum=1e4)
sc.pp.log1p(v3_processed)

# get the highly variable genes, needed later
sc.pp.highly_variable_genes(v2_processed, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(v3_processed, min_mean=0.0125, max_mean=3, min_disp=0.5)

# regress out the mitochondrial fraction, not sure whether or not to regress out the counts
sc.pp.regress_out(v2_processed, ['n_counts', 'percent_mito'])
#sc.pp.regress_out(v2_processed, ['percent_mito'])
sc.pp.regress_out(v3_processed, ['n_counts', 'percent_mito'])
#sc.pp.regress_out(v3_processed, ['percent_mito'])

# Find out what the max_value 10 actually means "Clip values exceeding standard deviation 10", not sure if needed
sc.pp.scale(v2_processed, max_value=10)
#sc.pp.scale(v2_processed)
sc.pp.scale(v3_processed, max_value=10)
#sc.pp.scale(v3_processed)

# write to files as backup
v2_processed.write(lanes_v2_scaled_loc)
v3_processed.write(lanes_v2_scaled_loc)
