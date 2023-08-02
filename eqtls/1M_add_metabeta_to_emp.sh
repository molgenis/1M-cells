#!/bin/bash

###################################################################
#Script Name	  : 1M_add_metabeta_to_emp.sh
#Description	  : add the meta beta to EMP output files
#Args           : 
#Author       	: Roy Oelen
#example        : ./1M_add_metabeta_to_emp.sh /groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/QTL_mapping/eQTL/sct_mqc_demux_lores_20201106_no_confine/results/
#
###################################################################
# location of EMP
EMP_LOC='/groups/umcg-biogen/tmp01/tools/eqtl-mapping-pipeline-1.4.9b-SNAPSHOT/eqtl-mapping-pipeline.jar'
# the command
EMP_COMMAND='java -Xmx1g -jar '${EMP_LOC}

# get as a parameter
BASE_DIR=$1

# check each condition
for cond in ${BASE_DIR}/*
    do
        # paste the condition folder
        cond_folder=${cond}/
        # check each cell type
        for ct in ${cond_folder}/*
            do
                # past the cell type folder
                ct_folder=${ct}/
                # get the input and output files
                SNPQC_LOC=${ct_folder}/'SNPQCLog.txt.gz'
                MAFFILE_LOC=${ct_folder}/'MAFFILE.txt.gz'
                IN_EQTL_FILE=${ct_folder}/'eQTLsFDR-ProbeLevel.txt.gz'
                OUT_EQTL_FILE=${ct_folder}/'eQTLsFDR-ProbeLevel-metabeta.txt.gz'
                # do the full commands
                ${EMP_COMMAND} \
                --mode util --getmaffromqclog \
                --in ${SNPQC_LOC} \
                --out ${MAFFILE_LOC}
                ${EMP_COMMAND} \
                --mode util --calculatebeta \
                --in ${IN_EQTL_FILE} \
                --out ${OUT_EQTL_FILE} \
                --in2 ${MAFFILE_LOC}
            done
    done