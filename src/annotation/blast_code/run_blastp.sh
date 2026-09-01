#!/bin/bash

INPUT_FASTA=$1
SPLIT_TEMP_FOLDER=$2
BLAST_OUTPUT_FOLDER=$3
OUTPUT_FILE=$4
THREADS=$5
TAXID=$6
TRANSLATION_TABLE=$7
LOCAL_BLAST_DB=$8
PROTEIN_DB=${9:-}
MINIMAL_BLAST=${10:-false}
COEFFICIENTS_FILE=${11:-}
NUM_PLOT_HITS=${12:-10}

MINIMAL_BLAST_FLAG=""
if [[ "$MINIMAL_BLAST" == "true" || "$MINIMAL_BLAST" == "True" || "$MINIMAL_BLAST" == "1" ]] ; then
  MINIMAL_BLAST_FLAG="--minimal_blast"
fi

if [[ -z $TAXID ]] ; then
  TAXID=0
fi

# if there is no local blast db provided, it will be passed as 0
if [[ $LOCAL_BLAST_DB == "0" ]] ; then
  LOCAL_BLAST_DB=""
else
  LOCAL_BLAST_DB="--local_blast_db $LOCAL_BLAST_DB"
fi

# if PROTEIN_DB is not empty, add the "--protein_db $PROTEIN_DB" flag, otherwise set it to an empty string
if [[ -z $PROTEIN_DB ]] ; then
  PROTEIN_DB_FLAG=""
else
  PROTEIN_DB_FLAG="--protein_db $PROTEIN_DB"
fi

mkdir -p $SPLIT_TEMP_FOLDER
mkdir -p $BLAST_OUTPUT_FOLDER

python src/annotation/blast_code/run_blastp.py \
  --input $INPUT_FASTA \
  --split_folder $SPLIT_TEMP_FOLDER \
  --blast_folder $BLAST_OUTPUT_FOLDER \
  --max_workers $THREADS \
  --taxid "$TAXID" \
  $MINIMAL_BLAST_FLAG \
  --coefficients "$COEFFICIENTS_FILE" \
  --num_plot_hits "$NUM_PLOT_HITS" \
  --translation_table $TRANSLATION_TABLE ${LOCAL_BLAST_DB} \
  $PROTEIN_DB_FLAG # flag will be provided as --protein_db "/path/to/db" or will be empty # flag will be provided as --local_blast_db "/path/to/db" or will be empty
 

Rscript --vanilla src/annotation/blast_code/merge_blastp_with_GO.R \
  --blast_folder $BLAST_OUTPUT_FOLDER \
  --output_file $OUTPUT_FILE \
  --max_workers $THREADS

rm -r $SPLIT_TEMP_FOLDER
rm -r $BLAST_OUTPUT_FOLDER
