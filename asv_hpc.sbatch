#!/bin/bash -l
#SBATCH -A naiss2025-22-1225
#SBATCH -p shared
#SBATCH -J ASV_dada2
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --output=/cfs/klemming/home/a/almcama/projects/asv_pipeline/logs/%x-%j.out

set -Eeuo pipefail

#SLURM CPU allocation 
export LC_ALL=C.UTF-8
export LANG=C.UTF-8
export LC_CTYPE=C.UTF-8
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# ---- Modules ----
ml purge >/dev/null 2>&1 || true
ml PDC >/dev/null 2>&1 || true
ml load PrgEnv-gnu/8.6.0
ml load cpeGNU/24.11
ml load R/4.4.2-cpeGNU-24.11

# sanity checks
command -v R       >/dev/null || { echo "R not in PATH"; ml list; exit 1; }
command -v Rscript >/dev/null || { echo "Rscript not in PATH"; ml list; exit 1; }

# ---- User library path (per R version/platform) ----
export R_LIBS_USER="/cfs/klemming/home/a/almcama/.R/24.11/4.4.2/library"

# sanity check
command -v Rscript >/dev/null || { echo "Rscript not in PATH"; ml list; exit 1; }
Rscript -e 'cat("Using:", R.version.string, "\n"); print(.libPaths())'

#Paths
IN="/cfs/klemming/home/a/almcama/projects/asv_pipeline/data/raw_fastq"
OUT="/cfs/klemming/home/a/almcama/projects/asv_pipeline/results/asv_out_${SLURM_JOB_ID}"
TAX_TRAIN="/cfs/klemming/home/a/almcama/projects/asv_pipeline/data/db/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz"
TAX_SPEC="/cfs/klemming/home/a/almcama/projects/asv_pipeline/data/db/silva/silva_v138.2_assignSpecies.fa.gz"

# Prep dirs & record env
mkdir -p /cfs/klemming/home/a/almcama/projects/asv_pipeline/logs "$OUT"
ml list |& tee "$OUT/modules_loaded.txt"

# Run R script (note: file only needs read perms; not executable)
Rscript /cfs/klemming/home/a/almcama/projects/asv_pipeline/scripts/ASV_dada2.r \
  --in "$IN" \
  --out "$OUT" \
  --pattern_f "_R1" \
  --pattern_r "_R2" \
  --trunc_len_f 0 \
  --trunc_len_r 0 \
  --trunc_q 2 \
  --max_ee_f 2 \
  --max_ee_r 2 \
  --subset_err 1000000 \
  --pool pseudo \
  --n_threads "$SLURM_CPUS_PER_TASK" \
  --taxonomy_train "$TAX_TRAIN" \
  --taxonomy_species "$TAX_SPEC" \
  --make-quality-plots \
  --quality-plot-samples 16 \
  --make-error-plots \
  |& tee "$OUT/R_console_${SLURM_JOB_ID}.log"
