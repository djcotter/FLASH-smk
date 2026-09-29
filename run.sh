#!/bin/bash
#
#SBATCH -p horence
#SBATCH --time=2-00:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

set -euo pipefail

# ============================================================
# Configuration — EDIT THESE PATHS for your machine
# ============================================================
MINIFORGE_DIR="/oak/stanford/groups/horence/chester/dabs_ref/miniforge3"
FLASH_ENV="flash"

# ============================================================
# Usage:
#   ./run_flash.sh <Snakefile> [rule1 rule2 ...]
#   MODE=embeddings|genomes|ohe|umap ./run_flash.sh <Snakefile> [rule1 rule2 ...]
#
# Provide at least the Snakefile to run. For example, to reproduce the bacterial run:
#    ./run_flash.sh Snakefile
#
# MODE selects the target meta-rule (default: embeddings).
# Any rules listed after the Snakefile (space-separated) are forced with -R.
# ============================================================

MODE="${MODE:-embeddings}"
SNAKEMAKE_FILE="${1:?Usage: $0 <Snakefile> [rules_to_force...]}"
shift || true
FORCE_RULES=("$@")

# ---- activate conda environment ----
source "${MINIFORGE_DIR}/etc/profile.d/conda.sh"
conda activate "$FLASH_ENV"

# ---- pick target meta-rule for the chosen mode ----
case "$MODE" in
    embeddings) TARGET="all_embeddings" ;;
    genomes)    TARGET="all_genomes" ;;
    ohe)        TARGET="all_ohe" ;;
    umap)       TARGET="all_umap" ;;
    *) echo "ERROR: MODE must be one of: embeddings, genomes, ohe, umap (got '$MODE')" >&2; exit 1 ;;
esac

# ---- unlock and run ----
SNAKEMAKE_CMD=(snakemake --sdm conda --use-conda --conda-base-path "$MINIFORGE_DIR" --profile slurm_profile/ -s "$SNAKEMAKE_FILE")
snakemake --unlock -s "$SNAKEMAKE_FILE" || true
if [ "${#FORCE_RULES[@]}" -gt 0 ]; then
    "${SNAKEMAKE_CMD[@]}" "$TARGET" -R "${FORCE_RULES[@]}"
else
    "${SNAKEMAKE_CMD[@]}" "$TARGET"
fi