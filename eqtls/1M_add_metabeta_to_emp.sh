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
EMP_LOC='/groups/umcg-biogen/tmp04/tools//eqtl-mapping-pipeline-1.4.9a-SNAPSHOT/eqtl-mapping-pipeline.jar'
# the command
EMP_COMMAND='java -Xmx1g -jar '${EMP_LOC}

# get as a parameter
# BASE_DIR='/groups/umcg-franke-scrna/tmp04/releases/wijst-2020-hg19/v1/QTL_mapping/eQTL/output/sct_mqc_demux_lores_20201106_no_confine/results/'
# OUT_DIR='/groups/umcg-franke-scrna/tmp04/releases/wijst-2020-hg19/v1/QTL_mapping/eQTL/output/sct_mqc_demux_lores_20201106_no_confine_metabeta/'
BASE_DIR=$1
OUT_DIR=$2

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
                # extract base ct name
                base_ct=$(basename ${ct})
                # extract base condition name
                base_cond=$(basename ${cond})
                # get the input and output files
                SNPQC_LOC=${ct_folder}/'SNPQCLog.txt.gz'
                MAFFILE_LOC=${ct_folder}/'MAFFILE.txt.gz'
                IN_EQTL_FILE=${ct_folder}/'eQTLsFDR-ProbeLevel.txt.gz'
                # make output directory
                OUT_EQTL_DIR=${OUT_DIR}/${base_cond}/${base_ct}/
                mkdir -p ${OUT_EQTL_DIR}
                # paste output filename
                OUT_EQTL_FILE=${OUT_EQTL_DIR}/'eQTLsFDR-ProbeLevel.txt.gz'
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