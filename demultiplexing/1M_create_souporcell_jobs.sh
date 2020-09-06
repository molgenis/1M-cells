#!/usr/bin/env bash

#directory and file listings
LANES_DIR="/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/alignment/cellranger_output/"
LANE_READGROUP_APPEND="outs/possorted_genome_bam.bam"
LANE_BARCODE_APPEND="outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
INDIVIDUAL_GENOTYPES="/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/genotype/1M_genotypes_cytoSNP_MAF_0.05_exons/"
GENOTYPE_APPEND=".vcf"
OUTPUT_DIR="/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/demultiplexing/souporcell/souporcell_output/"
JOB_DIR="/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/demultiplexing/souporcell/jobs/"
SOUPOR_IMAGE="/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/demultiplexing/souporcell/image/souporcell.sif"
GENOME_LOC="/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/external_software/refdata-cellranger-hg19-3.0.0/fasta/genome.fa"

LANES=("180920_lane1" "180920_lane2" "180921_lane1" "180921_lane2" "180925_lane1" "180925_lane2" "180926_lane1" "180926_lane2" "180927_lane1" "180927_lane2" "180928_lane1" "180928_lane2" "181003_lane1" "181003_lane2" "181003_lane3" "181004_lane1" "181004_lane2" "181010_lane1" "181010_lane2" "181011_lane1" "181011_lane2" "181017_lane1" "181017_lane2" "181018_lane1" "181018_lane2" "181022_lane1" "181022_lane2" "181023_lane1" "181023_lane2" "181024_lane1" "181024_lane2" "181024_lane3" "181025_lane1" "181025_lane2" "181029_lane1" "181029_lane2" "181030_lane1" "181030_lane2" "181105_lane1" "181105_lane2" "181106_lane1" "181106_lane2" "181107_lane1" "181107_lane2" "181108_lane1" "181108_lane2" "181114_lane1" "181114_lane2" "181115_lane1" "181115_lane2" "181121_lane1" "181121_lane2" "181122_lane1" "181122_lane2" "181128_lane1" "181128_lane2" "181129_lane1" "181129_lane2" "181213_lane1" "181213_lane2" "181213_lane3" "181214_lane1" "181214_lane2" "181218_lane1" "181218_lane2" "181218_lane3" "181218_lane4" "181219_lane1" "181219_lane2" "181219_lane3" "181231_lane1" "181231_lane2" "190101_lane1" "190101_lane2" "190109_lane1" "190109_lane2" "190110_lane1" "190110_lane2" "190114_lane1" "190114_lane2" "190115_lane1" "190115_lane2" "190116_lane1" "190116_lane2" "190117_lane1" "190117_lane2" "190121_lane1" "190121_lane2" "190122_lane1" "190122_lane2" "190123_lane1" "190123_lane2" "190124_lane1" "190124_lane2" "190129_lane1" "190129_lane2" "190130_lane1" "190130_lane2" "190131_lane1" "190131_lane2" "190201_lane1" "190201_lane2" "190204_lane1" "190204_lane2" "190205_lane1")

for dir in ${LANES[*]} ; do

    lane_id=${dir}
    output_folder=${OUTPUT_DIR}/${lane_id}/
    mkdir -p ${output_folder}
    output_job=${JOB_DIR}/soup_${lane_id}_1m_mmaf005_SBATCH.sh

    echo -e "#!/usr/bin/env bash
#SBATCH --job-name=soup_${lane_id}_mmaf005
#SBATCH --output=soup_${lane_id}_mmaf005.out
#SBATCH --error=soup_${lane_id}_mmaf005.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
        set -e
        cd ${output_folder}
        export SINGULARITY_BINDPATH=\"/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/,/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/external_software/\"
        singularity exec ${SOUPOR_IMAGE} souporcell_pipeline.py \
-i ${LANES_DIR}/${lane_id}/outs/possorted_genome_bam.bam \
-b ${LANES_DIR}/${lane_id}/outs/filtered_feature_bc_matrix/barcodes.tsv \
-f ${GENOME_LOC} \
-t 8 \
-o ${output_folder}/ \
-k 8 \
--known_genotypes ${INDIVIDUAL_GENOTYPES}/${lane_id}${GENOTYPE_APPEND} \
--skip_remap True
        " >> ${output_job}
done
