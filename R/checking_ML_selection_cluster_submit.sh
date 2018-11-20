#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -J all_gene_selection
#SBATCH --partition=general
#SBATCH -C haswell
#SBATCH -t 120:00:00
#SBATCH --mem-per-cpu=30000
#SBATCH --mail-user=vincent.cannataro@yale.edu
#SBATCH --mail-type=FAIL

module load R
Rscript checking_ML_selection_results_on_cluster.R
