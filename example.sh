#!/bin/bash
#SBATCH --partition=common
#SBATCH --account=statdept
#SBATCH -n8

module load R/3.6.0
Rscript scripts/analysis_draft_v1.R
