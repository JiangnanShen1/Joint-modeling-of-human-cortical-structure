#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiangnan.shen@yale.edu
#SBATCH --job-name=remove_dup
#SBATCH --partition=pi_zhao,scavenge
#SBATCH --cpus-per-task=5
#SBATCH --mem=10g
#SBATCH --time=4:00:00
#SBATCH --array=1-22

chr=${SLURM_ARRAY_TASK_ID}
/gpfs/gibbs/pi/zhao/yc769/softwares/plink \
--bfile /gpfs/gibbs/pi/zhao/yz738/1000G/eur_chr${chr}_SNPmaf5 \
--exclude /gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/reference_panel/rm_rsID.txt \
--make-bed \
--out /gpfs/gibbs/pi/zhao/js4377/UKB/lava/Data/reference_panel/eur_chr${chr}_SNPmaf5_rmdu