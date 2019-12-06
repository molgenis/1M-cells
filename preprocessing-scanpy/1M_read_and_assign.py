###########################################################################################################################
# Authors: Dylan de Vries, Roy Oelen
# Name: 1M_read_and_assign.py
# Function: Read all lanes and combine them into a single object and write clusters
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


def get_tp_per_lane_and_ll(stim_mapping_loc, exp_to_ll_loc):
    """
    get the stimulation information from files
    :param stim_mapping_loc: the location of the stimulations file
    :param exp_to_ll_loc: the location of the exp-to-ll file
    :return: a dictionary with the lanes as keys, with then then another dictionary with the condition per sample in that lane
    """
    # init the mapping of the experiment to the LL id
    exp_to_ll = {}
    # init the double dictionary that contains the lanes and LLs per TP
    lanes_to_stims_mapping = {}
    try:
        # read the cell assignments file
        exp_to_ll_file = open(exp_to_ll_loc, "r")
        # check the lines
        for line in exp_to_ll_file.readlines():
            # 'someone' used multiple space instead of a tab in the file
            cleaned_line = re.sub('\s+', '\t', line).strip()
            # split by the tab
            split_line = cleaned_line.split("\t")
            # there might be an empty line at the end, so only use the lines that can evaluate to two entries
            if len(split_line) == 2:
                exp_to_ll[split_line[0]] = split_line[1]
        # map column names to indices
        column_names_to_indices = {}
        # start reading the file
        stim_mapping_file = open(stim_mapping_loc, "r")
        # state that we are reading from the first line
        first_line = True
        for line in stim_mapping_file.readlines():
            # make sure there are no spaces
            cleaned_line = re.sub('\s+', '\t', line).strip()
            # split by the tab
            split_line = cleaned_line.split("\t")
            # create the mapping if this is the first line
            if first_line:
                # map indices to column name
                for i in range(len(split_line)):
                    column_names_to_indices[split_line[i]] = i
                # we're done, so no longer on the first line
                first_line = False
            else:
                lane = ""
                tp_per_ll = {}
                for colname,index in column_names_to_indices.items():
                    # lane is special
                    if colname == "Lane":
                        lane = split_line[index]
                    # if the entry under this TP is not empty and not NA
                    elif split_line[index] != "" and split_line[index] != "NA":
                        # the entries are split by a comma
                        expNrs = split_line[index].split(",")
                        # now use the expNrs
                        for expNr in expNrs:
                            # get the LL ID instead
                            ll_id = exp_to_ll[expNr]
                            # state the TP(that was the column name) for this ll_id
                            tp_per_ll[ll_id] = colname
                # for the LLs in this lane, save their states
                lanes_to_stims_mapping[lane] = tp_per_ll
    except Exception as e:
        print(str(e))
        print("oopsy")
    return lanes_to_stims_mapping


def add_stim_tags(scanpy_object, tp_per_lane_and_ll, assignment_key="assignment_ll", batch_key="batch"):
    """
    add the stimulation data to a scanpy object
    :param batch_key:
    :param assignment_key:
    :param scanpy_object: the object to add stimulation information to
    :param tp_per_lane_and_ll: a dictionary as obtained from 'get_tp_per_lane_and_ll'
    :return: the scanpy object with the stimulations info added under .obs['tp']
    """
    # init the array to store the timepoints
    timepoints = []
    # go through the indices of the observations
    for i in range(len(scanpy_object.obs.index)):
        # grab the values we know for this observation
        lane = scanpy_object.obs[batch_key][i]
        assignment = scanpy_object.obs[assignment_key][i]
        tp = "unknown"
        # get the dictionary for this lane and then the timepoint for the sample in that lane
        if lane in tp_per_lane_and_ll.keys() and assignment in tp_per_lane_and_ll[lane].keys():
            tp = tp_per_lane_and_ll[lane][assignment]
        timepoints.append(tp)
    # set the fetched timepoints as an observation
    scanpy_object.obs['tp'] = timepoints
    return scanpy_object

def add_demux_tags(scanpy_object, demuxlet_output_loc, filename_regex='^\d{6}_lane\d{1}.best', lane_extract_regex='(^\d{6}_lane\d{1})', batch_key="batch"):
    """
    add demux assignment tags to a scanpy object
    :param scanpy_object: the object to add the demux tags to
    :param demuxlet_output_loc: the location of the .best files
    :param filename_regex: the filtering to only get the appropriate .best files
    :param lane_extract_regex: the filtering to get the lane name from the filename
    :param batch_key: the key to get the lane/batch from the observations of the scanpy object
    :return: the scanpy object with the demuxlet tags added
    """
    # read the assignments
    demux_assignments = get_demuxlet_assignments(demuxlet_output_loc, filename_regex, lane_extract_regex)
    # create the lists where we are going to store the values we get from the assignments file
    bests = []
    sng_1sts = []
    sng_llk1s = []
    llk12s = []
    # go through the indices of the observations
    for i in range(len(scanpy_object.obs.index)):
        # grab the values we know for this observation
        lane = scanpy_object.obs[batch_key][i]
        barcode_unsanitized = scanpy_object.obs.index[i]
        # we need to remove the trailing '-number'
        upto_index_to_grab = barcode_unsanitized.rindex('-')
        barcode = barcode_unsanitized[:upto_index_to_grab]
        # set the default values
        best = 'unknown'
        sng_1st = 'unknown'
        sng_llk1 = 0
        llk12 = 0
        # grab the data for this cell/barcode if possible
        if lane in demux_assignments.keys() and barcode in demux_assignments[lane].keys():
            best = demux_assignments[lane][barcode]['BEST']
            sng_1st = demux_assignments[lane][barcode]['SNG.1ST']
            sng_llk1 = demux_assignments[lane][barcode]['SNG.LLK1']
            llk12 = demux_assignments[lane][barcode]['LLK12']
        else:
            print("no assignment for {} in {}".format(barcode,lane))
        # add the values for for this cell/barcode
        bests.append(best)
        sng_1sts.append(sng_1st)
        sng_llk1s.append(sng_llk1)
        llk12s.append(llk12)
    # now add these observations
    scanpy_object.obs['BEST'] = bests
    scanpy_object.obs['SNG.1ST'] = sng_1sts
    scanpy_object.obs['SNG.LLK1'] = sng_llk1s
    scanpy_object.obs['LLK12'] = llk12s
    # return the object with the observations added
    return scanpy_object


def add_soup_tags(scanpy_object, souporcell_output_loc, filename_regex='^\d{6}_lane\d{1}.tsv', lane_extract_regex='(^\d{6}_lane\d{1})', batch_key="batch"):
    """
    add souporcell assignment tags to a scanpy object
    :param scanpy_object: the object to add the souporcell tags to
    :param demuxlet_output_loc: the location of the .tsv files
    :param filename_regex: the filtering to only get the appropriate .tsv files
    :param lane_extract_regex: the filtering to get the lane name from the filename
    :param batch_key: the key to get the lane/batch from the observations of the scanpy object
    :return: the scanpy object with the souporcell tags added
    """
    # read the assignments
    soupor_assignments = get_souporcell_assignments(souporcell_output_loc, filename_regex, lane_extract_regex)
    # create the lists where we are going to store the values we get from the assignments file
    assignments = []
    statuss = []
    # go through the indices of the observations
    for i in range(len(scanpy_object.obs.index)):
        # grab the values we know for this observation
        lane = scanpy_object.obs[batch_key][i]
        barcode_unsanitized = scanpy_object.obs.index[i]
        # we need to remove the trailing '-number'
        upto_index_to_grab = barcode_unsanitized.rindex('-')
        barcode = barcode_unsanitized[:upto_index_to_grab]
        # set the default values
        assignment = 'unknown'
        status = 'unknown'
        # grab the data for this cell/barcode if possible
        if lane in soupor_assignments.keys() and barcode in soupor_assignments[lane].keys():
            assignment = soupor_assignments[lane][barcode]['assignment_ll']
            status = soupor_assignments[lane][barcode]['status']
        else:
            print("no assignment for {} in {}".format(barcode,lane))
        # add the values for for this cell/barcode
        assignments.append(assignment)
        statuss.append(status)
    # now add these observations
    scanpy_object.obs['assignment_ll'] = assignments
    scanpy_object.obs['status'] = statuss
    # return the object with the observations added
    return scanpy_object


def get_demuxlet_assignments(demuxlet_output_loc, filename_regex='^\d{6}_lane\d{1}.best', lane_extract_regex='(^\d{6}_lane\d{1})'):
    """
    get the demuxlet assignments from a location
    :param demuxlet_output_loc: the location of the .best demuxlet files
    :param filename_regex: the regex used to just select the .best files you want
    :param lane_extract_regex: the regex used to extract the lane name from the filename
    :return: a dictionary per lane, with a dictionary per barcode, with a dictionary with the assignment data
    """
    # list the .best files
    demux_files_unfiltered = os.listdir(demuxlet_output_loc)
    demux_regex = re.compile(filename_regex) # should match 180920_lane1.best etc.
    # get only the directories matching our regex
    demux_files = filter(demux_regex.match, demux_files_unfiltered)
    # store the results per lane
    demux_per_lane = {}
    for lane in demux_files:
        # store the assignment per barcode
        assignment_per_barcode = {}
        # we need to know if we're at the first line
        at_first_line = True
        # read each demuxlet output file
        try:
            # open the demux output
            demuxlet_output = open('/'.join([demuxlet_output_loc, lane]), "r")
            # we can try to extract the lane from the filename
            lane_sanitized = lane
            lane_search = re.search(lane_extract_regex, lane)
            if lane_search:
                lane_sanitized = lane_search.group(1)
            # we'll store a mapping of the index of a column, to its column name in the header
            colname_to_index = {}
            # go through the lines of the file
            for line in demuxlet_output.readlines():
                # split by tab and remove the trailing newline
                split_line = line.strip().split("\t")
                # map the colnames to the indices if this is the first line
                if at_first_line:
                    for i in range(len(split_line)):
                        colname_to_index[split_line[i]] = i
                    # we're done, so no longer on the first line
                    at_first_line = False
                # if we're not at the first line, we can use the colname to index mapping to get the interested data
                else:
                    # fetch interested values
                    barcode = split_line[colname_to_index['BARCODE']]
                    assignment = split_line[colname_to_index["BEST"]]
                    assignment_sgn1 = split_line[colname_to_index["SNG.1ST"]]
                    assignment_sgn1_llk = float(split_line[colname_to_index["SNG.LLK1"]])
                    assignment_llk12 = float(split_line[colname_to_index["LLK12"]])
                    # put values into dictionary
                    values_this_barcode = {'BARCODE': barcode, 'BEST': assignment, 'SNG.1ST': assignment_sgn1, 'SNG.LLK1': assignment_sgn1_llk, 'LLK12': assignment_llk12}
                    # set as values for this barcode
                    assignment_per_barcode[barcode] = values_this_barcode
            # set the barcodes information we fetched for this lane
            demux_per_lane[lane_sanitized] = assignment_per_barcode
        except Exception as e:
            print(' '.join["error reading lane", lane, str(e)])
            continue
    return demux_per_lane


def get_souporcell_assignments(souporcell_output_loc, filename_regex='^\d{6}_lane\d{1}.tsv', lane_extract_regex='(^\d{6}_lane\d{1})'):
    """
    get the souporcell assignments from a location
    :param souporcell_output_loc: the location of the souporcell files with the LL assignment in the last column
    :param filename_regex: the regex used to just select the .best files you want
    :param lane_extract_regex: the regex used to extract the lane name from the filename
    :return: a dictionary per lane, with a dictionary per barcode, with a dictionary with the assignment data
    """
    # list the .best files
    soupor_files_unfiltered = os.listdir(souporcell_output_loc)
    soupor_regex = re.compile(filename_regex) # should match 180920_lane1.best etc.
    # get only the directories matching our regex
    soupor_files = filter(soupor_regex.match, soupor_files_unfiltered)
    # store the results per lane
    soupor_per_lane = {}
    for lane in soupor_files:
        # store the assignment per barcode
        assignment_per_barcode = {}
        # we need to know if we're at the first line
        at_first_line = True
        # read each demuxlet output file
        try:
            # open the demux output
            soupor_output = open('/'.join([souporcell_output_loc, lane]), "r")
            # we can try to extract the lane from the filename
            lane_sanitized = lane
            lane_search = re.search(lane_extract_regex, lane)
            if lane_search:
                lane_sanitized = lane_search.group(1)
            # we'll store a mapping of the index of a column, to its column name in the header
            colname_to_index = {}
            # go through the lines of the file
            for line in soupor_output.readlines():
                # split by tab and remove the trailing newline
                split_line = line.strip().split("\t")
                # map the colnames to the indices if this is the first line
                if at_first_line:
                    for i in range(len(split_line)):
                        colname_to_index[split_line[i]] = i
                    # we're done, so no longer on the first line
                    at_first_line = False
                # if we're not at the first line, we can use the colname to index mapping to get the interested data
                else:
                    # fetch interested values
                    barcode = split_line[colname_to_index['barcode']]
                    assignment = split_line[colname_to_index["assignment_ll"]]
                    status = split_line[colname_to_index["status"]]
                    # put values into dictionary
                    values_this_barcode = {'barcode': barcode, 'assignment_ll': assignment, 'status': status}
                    # set as values for this barcode
                    assignment_per_barcode[barcode] = values_this_barcode
            # set the barcodes information we fetched for this lane
            soupor_per_lane[lane_sanitized] = assignment_per_barcode
        except Exception as e:
            print(' '.join["error reading lane", lane, str(e)])
            continue
    return soupor_per_lane


###########################################################################################################################
#
# Main code
#
###########################################################################################################################

sc.settings.verbosity = 3           # verbosity: errors (0), warnings (1), info (2), hints (3)

# the output of the cellranger output per lane
cellranger_output_loc = "/groups/umcg-lld/tmp04/projects/1MCellRNAseq/processed/cellranger_output/"

# the output h5ad locations of each step
object_loc = "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/scanpy_preprocess_samples/objects/"
lanes_v2_raw_loc = ''.join([object_loc, '1M_cells_scanpy_object_v2_raw.h5ad'])
lanes_v3_raw_loc = ''.join([object_loc, '1M_cells_scanpy_object_v3_raw.h5ad'])
lanes_v2_assigned_loc = ''.join([object_loc, '1M_cells_scanpy_object_v2_assigned.h5ad'])
lanes_v3_assigned_loc = ''.join([object_loc, '1M_cells_scanpy_object_v3_assigned.h5ad'])

# list the output directories
lane_files_unfiltered = os.listdir(cellranger_output_loc)
lane_regex = re.compile(r'^\d{6}_lane\d{1}$') # should match 180920_lane1 etc.
# get only the directories matching our regex
lane_files = filter(lane_regex.match, lane_files_unfiltered)

# set the lists to store that data for the lanes for both chemistries
lane_data_v2 = []
lane_data_v3 = []
v3_lanes = {}
v2_lanes = {}


for i, lane in enumerate(lane_files):
    print(lane)
    # read the 10X count matrix
    lane_data = sc.read_10x_mtx(''.join([cellranger_output_loc, "/{}/outs/filtered_feature_bc_matrix".format(lane)]), var_names='gene_symbols', cache=True)
    # if the lanes is from 2018, it is a V2 chem lane
    if lane.startswith("18"):
        # so add to the V2 lanes
        v2_lanes[len(v2_lanes.keys())] = lane
        lane_data_v2.append(lane_data)
    # if the lane is from 2019, it is a V3 chem lane
    elif lane.startswith("19"):
        # so add to the V3 lanes
        v3_lanes[len(v3_lanes.keys())] = lane
        lane_data_v3.append(lane_data)

v2_combined = lane_data_v2[0].concatenate(lane_data_v2[1:])
v2_lane_assignment = []
for batch in v2_combined.obs['batch']:
    v2_lane_assignment.append(v2_lanes[int(batch)])

# set the batches
v2_combined.obs['batch'] = v2_lane_assignment


v3_combined = lane_data_v3[0].concatenate(lane_data_v3[1:])
v3_lane_assignment = []

for batch in v3_combined.obs['batch']:
    v3_lane_assignment.append(v3_lanes[int(batch)])

# set the batches
v3_combined.obs['batch'] = v3_lane_assignment

# write to files as backup
v2_combined.write(lanes_v2_raw_loc)
v3_combined.write(lanes_v3_raw_loc)

# add demuxlet assignments
demux_loc = '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/demux_relevant_samples/demux_output_cytosnp/'
v2_assigned = add_demux_tags(v2_combined, demux_loc, filename_regex='^\d{6}_lane\d{1}_sorted_hfixed.best', lane_extract_regex='(^\d{6}_lane\d{1})')
v3_assigned = add_demux_tags(v3_combined, demux_loc, filename_regex='^\d{6}_lane\d{1}_sorted_hfixed.best', lane_extract_regex='(^\d{6}_lane\d{1})')

# add souporcell assignments UPDATE THIS DIRECTORY WHEN NEW OUTPUT IS DONE
soup_loc = '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/soupor_matchgts/doublet_out_gts/'
v2_assigned = add_soup_tags(v2_assigned, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/soupor_matchgts/doublet_out_gts/', filename_regex='^\d{6}_lane\d{1}_doub_gt.tsv', lane_extract_regex='(^\d{6}_lane\d{1})')
v3_assigned = add_soup_tags(v3_assigned, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/soupor_matchgts/doublet_out_gts/', filename_regex='^\d{6}_lane\d{1}_doub_gt.tsv', lane_extract_regex='(^\d{6}_lane\d{1})')

# fetch the timepoint assignments
exp_id_file = "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/scanpy_preprocess_samples/descriptions/exp_to_ll.txt"
lane_to_tp_file = "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/scanpy_preprocess_samples/descriptions/lane_to_tp.txt"
tp_per_lane_and_ll = get_tp_per_lane_and_ll(lane_to_tp_file, exp_id_file)
# add the timepoint assignments
v2_assigned = add_stim_tags(v2_assigned, tp_per_lane_and_ll)
v3_assigned = add_stim_tags(v3_assigned, tp_per_lane_and_ll)

# remove the cells without timepoints
v2_assigned = v2_assigned[v2_assigned.obs['tp'] != 'unknown']
v3_assigned = v3_assigned[v3_assigned.obs['tp'] != 'unknown']

# write to files as backup
v2_assigned.write(lanes_v2_assigned_loc)
v3_assigned.write(lanes_v3_assigned_loc)


