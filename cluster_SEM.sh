#!/bin/bash
#SBATCH --partition=herringlab-low,statdept-low
#SBATCH --account=statdept
#SBATCH -c8

singularity run -B /work/sta790/ff31:/work/sta790/ff31 felpo_ex_latest.sif factor_interactions/codes
/cluster_job.R