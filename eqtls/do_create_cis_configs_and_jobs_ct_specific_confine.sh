#!/usr/bin/env bash
#FEATURE_FOLDER="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v2_sct_mqc_demux_hires_20201029_reclassified_T/"
#FEATURE_FOLDER_2="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v3_sct_mqc_demux_hires_20201106_reclassified_T/"
FEATURE_FOLDER="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v2_sct_mqc_demux_highres_20210905/"
FEATURE_FOLDER_2="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v3_sct_mqc_demux_highres_20210905/"
OUTPUT_DIR="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_highres_20210905_eqtlgenlead_perctsig/"
CONFINEMENTS_PREPEND="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/sig_hires_20211008_reclassified_T_eqtlgenlead_"
#CONFINEMENTS_PREPEND="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/sct_mqc_demux_highres_20210905_eqtlgenlead_"
CONFINEMENTS_APPEND="_confinement.txt"

stim_dirs=${FEATURE_FOLDER}"UT_vs_*/"
for dir in ${stim_dirs}
do
    path=${dir%*/}
    stimulation=${path##*/}

    current_feature_folder=${FEATURE_FOLDER}"/"${stimulation}
    current_feature_folder_2=${FEATURE_FOLDER_2}"/"${stimulation}
    config_dir=${OUTPUT_DIR}"/config/"${stimulation}
    result_dir=${OUTPUT_DIR}"/results/"${stimulation}
    job_dir=${OUTPUT_DIR}"/jobs/"${stimulation}

    ./create_cis_configs_and_jobs_ct_specific_confine.sh \
        ${current_feature_folder} \
        ${current_feature_folder_2} \
        ${config_dir} \
        ${result_dir} \
        ${job_dir} \
        ${CONFINEMENTS_PREPEND} \
        ${CONFINEMENTS_APPEND}
done
