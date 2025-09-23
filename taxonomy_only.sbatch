#!/bin/bash -l
#SBATCH -A naiss2025-22-1225
#SBATCH -p shared
#SBATCH -J tax_only
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH --output=/cfs/klemming/home/a/almcama/projects/asv_pipeline/logs/%x-%j.out

set -Eeuo pipefail
export LC_ALL=C.UTF-8 LANG=C.UTF-8 LC_CTYPE=C.UTF-8
export OMP_NUM_THREADS=1 
export MKL_NUM_THREADS=1 
export OPENBLAS_NUM_THREADS=1

ml purge >/dev/null 2>&1 || true
ml PDC >/dev/null 2>&1 || true
ml load PrgEnv-gnu/8.6.0
ml load cpeGNU/24.11
ml load R/4.4.2-cpeGNU-24.11
export R_LIBS_USER="/cfs/klemming/home/a/almcama/.R/24.11/4.4.2/library"

Rscript /cfs/klemming/home/a/almcama/projects/asv_pipeline/scripts/taxonomy_only.R \
  /cfs/klemming/home/a/almcama/projects/asv_pipeline/results/asv_out_12428927 \
  /cfs/klemming/home/a/almcama/projects/asv_pipeline/data/db/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz \
  $SLURM_CPUS_PER_TASK \
  /cfs/klemming/home/a/almcama/projects/asv_pipeline/data/db/silva/silva_v138.2_assignSpecies.fa.gz
