#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

# Define the constants
HYENA_WDR="/oak/stanford/groups/horence/dcotter1/utility_files/hyena_models/hyena_wdr/"
PYTHON_SCRIPT="/oak/stanford/groups/horence/dcotter1/utility_files/hyena_models/dev_embedder.py"
MODEL_CFG="/oak/stanford/groups/horence/dcotter1/utility_files/hyena_models/aws_topES_multi_target/hg38_256dim_config.yml"
MODEL_CKPT="/oak/stanford/groups/horence/dcotter1/utility_files/hyena_models/aws_topES_multi_target/weights.ckpt"
SINGULARITY_IMG="/home/groups/horence/hyena-dna-nt6_latest.sif"

# Assign command line arguments to variables
input_file=$1
output_file=$2

# Call the Python script with the provided arguments
cd ${HYENA_WDR}
singularity run --nv --writable-tmpfs -B ${HYENA_WDR} ${SINGULARITY_IMG} python ${PYTHON_SCRIPT} --model_cfg ${MODEL_CFG} --ckpt_path ${MODEL_CKPT} --seq_file ${input_file} --output_file ${output_file} --max_seqlen 160000 --nlayers 8 --batch_size 500
