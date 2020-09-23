#!/usr/bin/env bash

# paths and filenames
EQTL_OUTPUT_DIR="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/"
CONDITIONS=("UT" "3hCA" "24hCA" "3hMTB" "24hMTB" "3hPA" "24hPA")
RE_CONDITIONS=("UT_vs_3hCA" "UT_vs_24hCA" "UT_vs_3hMTB" "UT_vs_24hMTB" "UT_vs_3hPA" "UT_vs_24hPA")
CELL_TYPES=("bulk" "B" "CD4T" "CD8T" "DC" "monocyte" "NK")
CELL_TYPE_APPEND="_expression"
EQTL_RESULT_FILE="eQTLsFDR0.05-ProbeLevel.txt.gz"
EQTL_TABLE_LOC="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/20200729_eqtls_numbers.txt"
REQTL_TABLE_LOC="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/20200729_reqtls_numbers.txt"

printf "\t" > ${EQTL_TABLE_LOC}
printf "%s\t" "${CELL_TYPES[@]}" >> ${EQTL_TABLE_LOC}
printf "\n" >> ${EQTL_TABLE_LOC}

# check condition
for condition in ${CONDITIONS[*]} ; do
  printf ${condition} >> ${EQTL_TABLE_LOC}
  printf "\t" >> ${EQTL_TABLE_LOC}
  # create array
  cell_type_numbers=()
  # check each cell type
  for cell_type in ${CELL_TYPES[*]} ; do
    # paste the full path together
    eqtl_result_path=${EQTL_OUTPUT_DIR}/${condition}/${cell_type}${CELL_TYPE_APPEND}/${EQTL_RESULT_FILE}
    # get the number of lines in the file
    eqtl_result_nroflines=$(zcat ${eqtl_result_path} | wc -l)
    # substract the header
    eqtl_result_nr=$((${eqtl_result_nroflines}-1))
    # add to the array
    cell_type_numbers+=(${eqtl_result_nr})
  done
  printf "%s\t" "${cell_type_numbers[@]}" >> ${EQTL_TABLE_LOC}
  printf "\n" >> ${EQTL_TABLE_LOC}
done

printf "\t" > ${REQTL_TABLE_LOC}
printf "%s\t" "${CELL_TYPES[@]}" >> ${REQTL_TABLE_LOC}
printf "\n" >> ${REQTL_TABLE_LOC}

# check condition
for re_condition in ${RE_CONDITIONS[*]} ; do
  printf ${re_condition} >> ${REQTL_TABLE_LOC}
  printf "\t" >> ${REQTL_TABLE_LOC}
  # create array
  cell_type_numbers=()
  # check each cell type
  for cell_type in ${CELL_TYPES[*]} ; do
    # paste the full path together
    reqtl_result_path=${EQTL_OUTPUT_DIR}/${re_condition}/${cell_type}${CELL_TYPE_APPEND}/${EQTL_RESULT_FILE}
    # get the number of lines in the file
    reqtl_result_nroflines=$(zcat ${reqtl_result_path} | wc -l)
    # substract the header
    reqtl_result_nr=$((${reqtl_result_nroflines}-1))
    # add to the array
    cell_type_numbers+=(${reqtl_result_nr})
  done
  printf "%s\t" "${cell_type_numbers[@]}" >> ${REQTL_TABLE_LOC}
  printf "\n" >> ${REQTL_TABLE_LOC}
done
