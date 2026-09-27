#!/bin/bash

# This script is called by Snakemake

PYTHON_SCRIPT=$(realpath $1)
SINGULARITY_IMG=$(realpath $2)
INPUT_FILE=$(realpath $3)
OUTPUT_FILE=$(realpath $4)

# Internal paths (now baked into the container)
export MODEL_CFG="/opt/models/bacterial_128dim_config.yml"
export MODEL_CKPT="/opt/models/weights.ckpt"

echo "GPU used for Hyena embedding:"
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,uuid,driver_version --format=csv,noheader
else
    echo "WARNING: nvidia-smi is not available." >&2
fi

singularity exec --nv ${SINGULARITY_IMG} \
    python ${PYTHON_SCRIPT} \
        --model_cfg ${MODEL_CFG} \
        --ckpt_path ${MODEL_CKPT} \
        --seq_file ${INPUT_FILE} \
        --output_file ${OUTPUT_FILE} \
        --max_seqlen 128000 \
        --nlayers 4 \
        --batch_size 100
