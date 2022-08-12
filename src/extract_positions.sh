#!/bin/bash
#SBATCH -p hpcbio
#SBATCH -J snppos
#SBATCH --mail-user lvclark@illinois.edu
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=100G

module load R/3.6.0-IGB-gcc-8.2.0

Rscript snp_positions_from_vcf.R --args yam336.m2M2vsnps_missing0.9.recode.vcf yam336.m2M2vsnps_missing0.9.recode.rds


