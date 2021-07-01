#!/usr/bin/env bash

ml R
SCRIPT_LOC="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/scripts/1M_create_substracted_features.R"
SPLIT_FEATURES_LOC=$1

UT="UT"
h3CA="3hCA"
h24CA="24hCA"
h3PA="3hPA"
h24PA="24hPA"
h3MTB="3hMTB"
h24MTB="24hMTB"

# create an array of stim conditions
stims=(${h3CA} ${h24CA} ${h3PA} ${h24PA} ${h3MTB} ${h24MTB})

# loop through files
files=${SPLIT_FEATURES_LOC}/${UT}"/*"
for f in ${files}
do
  # only create files for features
  if [[ ${f} =~ .*\.(tsv) ]]
  # ${f} is the features file full path, we also need just the name
  fbasename=$(basename "${f}" .tsv)
  then
    # create the combinations
    for stim in ${stims[*]} ; do
      ut_path=$f
      stim_path=${SPLIT_FEATURES_LOC}"/"${stim}"/"${fbasename}".tsv"
      output_folder=${SPLIT_FEATURES_LOC}"/"${UT}"_vs_"${stim}
      mkdir -p ${output_folder}
      Rscript ${SCRIPT_LOC} ${ut_path} ${stim_path} ${output_folder}"/"${fbasename}".tsv"
    done
  fi
done
