JOB_DIR='/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/LD_DB/jobs/'
POPULATION_FILE='/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/LD_DB/population_lists/1000g_eurs.txt'
FULL_VCF_DIR='/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/LD_DB/genotypes/'
SUB_VCF_DIR='/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/LD_DB/genotypes_eur/'
VCF_IN_PREPEND='ALL.chr'
VCF_IN_APPEND='.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
VCF_OUT_PREPEND='EUR.chr'
VCF_OUT_APPEND='.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

for chrom in ${CHROMS[*]} ; do


    output_job=${JOB_DIR}/subset_${VCF_OUT_PREPEND}${chrom}_SBATCH.sh

    echo -e "#!/usr/bin/env bash

#SBATCH --job-name=subset_${VCF_OUT_PREPEND}${chrom}
#SBATCH --output=${JOB_DIR}subset_${VCF_OUT_PREPEND}${chrom}.out
#SBATCH --error=${JOB_DIR}subset_${VCF_OUT_PREPEND}${chrom}.err
#SBATCH --time=5:59:59
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

set -e

ml HTSlib
ml BCFtools
bcftools view -S ${POPULATION_FILE} \
-O z \
-o ${SUB_VCF_DIR}${VCF_OUT_PREPEND}${chrom}${VCF_OUT_APPEND} \
--force-samples \
${FULL_VCF_DIR}${VCF_IN_PREPEND}${chrom}${VCF_IN_APPEND}

" >> ${output_job}
done
