#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -J all_gene_trinuc
#SBATCH --partition=general
#SBATCH -C haswell
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=30000
#SBATCH --mail-user=vincent.cannataro@yale.edu
#SBATCH --mail-type=FAIL

module load R
Rscript trinuc_gene_comp_builder.R
