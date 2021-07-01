#!/usr/bin/env bash

#directory and file listings
LANES_DIR="/groups/umcg-lld/tmp04/projects/1MCellRNAseq/processed/cellranger_output/"
LANE_READGROUP_APPEND="outs/possorted_genome_bam.bam"
LANE_BARCODE_APPEND="outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
INDIVIDUAL_GENOTYPES="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/demux_relevant_samples/genotype_per_lane_cytoSNP/"
OUTPUT_DIR="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/demux_relevant_samples/demux_output_cytosnp/"
JOB_DIR="/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/demux_relevant_samples/jobs_cytosnp"
DEMUXLET_DIR="/groups/umcg-wijmenga/tmp04/umcg-hbrugge/apps/demuxlet/"

#parameters used
TAG_GROUP="CB"
TAG_UMI="UB"
FIELD="GT"

for dir in "$LANES_DIR"/*lane*/ ; do

    lane_id=$(basename $dir)
    output_job=${JOB_DIR}/demux_${lane_id}_sorted_hfixed.sh

    echo -e "#!/usr/bin/env bash

#SBATCH --job-name=demux_${lane_id}_Roelen
#SBATCH --output=demux_${lane_id}_Roelen.out
#SBATCH --error=demux_${lane_id}_Roelen.err
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

	set -e

        ${DEMUXLET_DIR}demuxlet \\
                --sam ${LANES_DIR}${lane_id}/${LANE_READGROUP_APPEND} \\
                --tag-group ${TAG_GROUP} \\
                --tag-UMI ${TAG_UMI} \\
                --field ${FIELD} \\
                --vcf ${INDIVIDUAL_GENOTYPES}${lane_id}_sorted_hfixed.vcf \\
                --out ${OUTPUT_DIR}${lane_id}_sorted_hfixed \\
                --group-list ${LANES_DIR}${lane_id}/${LANE_BARCODE_APPEND}

        " >> ${output_job}
done
