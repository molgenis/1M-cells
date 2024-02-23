#!/bin/bash

###################################################################
#Script Name	  : 1M_download_qtl_results.sh
#Description	  : add the meta beta to EMP output files
#Args           : 
#Author       	: Roy Oelen
#example        : ./1M_add_metabeta_to_emp.sh /groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/QTL_mapping/eQTL/sct_mqc_demux_lores_20201106_no_confine/results/
#
###################################################################

# get as a parameter
BASE_DIR=$1
SYNC_DIR=$2

# get the folders in the base dir
conditions=$(rsync --list-only 'airlock+gearshift:'${BASE_DIR}/* | awk '{ print $5;}')
# check each condition
for cond in ${conditions[*]}
    do
        # paste the condition folder
        cond_folder=${BASE_DIR}/${cond}/
        # also create the folder locally
        local_dir=${SYNC_DIR}/${cond}/
        mkdir -p ${local_dir}
        # get every entry in this folder
        celltypes=$(rsync --list-only 'airlock+gearshift:'${cond_folder}/* | awk '{ print $5;}')
        # check each cell type
        for ct in ${celltypes[*]}
            do
                # past the cell type folder
                ct_folder=${cond_folder}/${ct}/
                # this is the file we need
                OUT_EQTL_FILE=${ct_folder}/'eQTLsFDR-ProbeLevel-metabeta.txt.gz'
                # make the folder where we will place the output
                mkdir -p ${local_dir}/${ct}/
                # now download the file
                rsync -aP 'airlock+gearshift:'${OUT_EQTL_FILE} ${local_dir}/${ct}/'eQTLsFDR-ProbeLevel-metabeta.txt.gz'
            done
    done