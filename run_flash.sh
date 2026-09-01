#!/bin/bash
set -euo pipefail

# ============================================================
# Configuration — EDIT THESE PATHS for your machine
# ============================================================
MINIFORGE_DIR="/path/to/miniforge3"       # your miniforge/conda install
FLASH_ENV="name_of_flash_env"             # name of the FLASH conda environment

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