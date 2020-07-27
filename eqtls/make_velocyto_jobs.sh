#!/usr/bin/env bash

lanes=( 180928_lane1 180928_lane2 181003_lane1 181003_lane2 181003_lane3 181011_lane1 181011_lane2 181022_lane1 181022_lane2 181023_lane1 181023_lane2 181024_lane1 181024_lane2 181024_lane3 181030_lane1 181030_lane2 181106_lane1 181106_lane2 181107_lane1 181107_lane2 181108_lane1 181108_lane2 181122_lane1 181122_lane2 181213_lane1 181213_lane2 181213_lane3 181214_lane1 181214_lane2 181218_lane1 181218_lane2 181219_lane1 181219_lane2 181219_lane3 190109_lane1 190109_lane2 190110_lane1 190110_lane2 190115_lane1 190115_lane2 190122_lane1 190122_lane2 190123_lane1 190123_lane2 190124_lane1 190124_lane2 190130_lane1 190130_lane2 190204_lane1 190204_lane2 )

for lane in "${lanes[@]}"
do

echo -e "#!/usr/bin/env bash
#SBATCH --job-name=velocyto${lane}
#SBATCH --output=velocyto${lane}.out
#SBATCH --error=velocyto${lane}.err
#SBATCH --chdir=/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/${lane}/
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml SAMtools
ml HDF5
ml Python
velocyto run10x -m /groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/genome_annotation/HG19_mm10_rmsk.gtf /groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/${lane} /groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/genome_annotation/cellranger-hg19-300.gtf
" >> /groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/jobs/velocyto_${lane}_SBATCH.sh
done
