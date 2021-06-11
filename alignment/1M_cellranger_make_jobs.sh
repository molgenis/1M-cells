#!/usr/bin/env bash

FASTQ_DIR="/groups/umcg-lld/prm02/rawdata/1MCellRNAseq/fastqs"
JOB_DIR="/groups/umcg-wijmenga/scr01/umcg-hbrugge/cellranger/jobs"
TMP_OUTPUT_DIR="/dev/shm/umcg-hbrugge/tmp/cellranger"
COPY_TO_DIR="calculon:/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/"

CELLRANGER_DIR="/groups/umcg-wijmenga/tmp04/umcg-ddevries/cellranger-3.0.2/cellranger-3.0.2/"
TRANSCRIPTOME_DIR="/groups/umcg-wijmenga/tmp04/umcg-ddevries/cellranger-3.0.2/refdata-cellranger-hg19-3.0.0"

mkdir ${JOB_DIR}
mkdir -p ${TMP_OUTPUT_DIR}

for dir in "$FASTQ_DIR"/*lane*/ ; do

    lane_id=$(basename $dir)
    output_job=${JOB_DIR}/run_${lane_id}.sh

    echo -e "#!/usr/bin/env bash

	set -e

        cd ${TMP_OUTPUT_DIR}

        ${CELLRANGER_DIR}cellranger count \\
                --id=${lane_id} \\
                --transcriptome=${TRANSCRIPTOME_DIR} \\
                --fastqs=${dir} \\
                --localmem=100

        rm -rf ${TMP_OUTPUT_DIR}/${lane_id}/SC_RNA_COUNTER_CS

        rsync -rv ${TMP_OUTPUT_DIR}/${lane_id} ${COPY_TO_DIR}

        rm -rf ${TMP_OUTPUT_DIR}/${lane_id}
        rm ${output_job}

        " >> ${output_job}
done

