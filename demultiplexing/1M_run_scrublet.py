# imports
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# base dir where each lane is located
cellranger_dir = '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/alignment/cellranger_output/'
# specific lane location
lane_append = '/outs/filtered_feature_bc_matrix/'
# grab the different lanes
lanes = os.listdir(cellranger_dir)
# we're expecting around 8% doublets, given the number of cells
expected_doublet_rate=0.08

# state where to save to
scrub_all_save_loc = '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/demultiplexing/scrublet/scrublet_assignment.tsv'

# init frame to save results
scrub_all =  None

#loop through the lanes we found in the directory
for lane in lanes:
  # get the full paths
  matrix_loc = "".join([cellranger_dir, lane, lane_append, 'matrix.mtx.gz'])
  features_loc = "".join([cellranger_dir, lane, lane_append, 'features.tsv'])
  barcodes_loc = "".join([cellranger_dir, lane, lane_append, 'barcodes.tsv'])
  # load the files
  print("".join(['reading lane ', lane]))
  counts_matrix = scipy.io.mmread(matrix_loc).T.tocsc()
  genes = np.array(scr.load_genes(features_loc, delimiter='\t', column=1))
  barcodes = pd.read_csv(barcodes_loc, sep= '\t', header=None)
  # perform scrublet
  scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate)
  # grabbing scores and assignments
  doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
  # add everything together in one frame
  print("".join(['finished scrubbing ', lane]))
  assignment = barcodes
  assignment['doublet'] = predicted_doublets
  assignment['doublet_score'] = doublet_scores
  assignment['barcode'] = assignment[0].str.slice(0, 16)
  assignment['lane_barcode'] = assignment['barcode']+'_'+lane
  # remove original column
  assignment = assignment.drop(0, axis=1)
  # add to existing dataframe
  print("".join(['merging lane ', lane]))
  if scrub_all is None:
    scrub_all = assignment
  else:
    scrub_all = pd.concat([scrub_all, assignment])

# save to file
scrub_all.to_csv(scrub_all_save_loc)
