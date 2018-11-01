#!/bin/bash
#SBATCH --job-name=impute
#SBATCH --error=imputation.err
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=55000
#SBATCH --exclude="midway2-0401"


module load bcftools
module load htslib
module load vcftools
module load R/3.4.3

Rscript tbtn_imputation_wgs.R $1 $2 $3
