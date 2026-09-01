from Bio import Entrez, SeqIO
import pandas as pd
import time
import argparse
import sys
import os
from os.path import join, basename
from concurrent.futures import ThreadPoolExecutor
from math import floor
import random
import pickle

MAX_RETRIES = 100
REQUEST_DELAY = 4  # Delay in seconds between requests
INITIAL_DELAY_RANGE = (0.1, 8)  # Range for initial random delay

BLAST_FEATURE_COLUMNS = [
    "query",
    "identity",
    "evalue",
    "qcovs",
    "features",
    "features_10000_window",
]


def blast_feature_columns(window):
    columns = BLAST_FEATURE_COLUMNS.copy()
    window_col = f"features_{window}_window"
    if window_col != "features_10000_window":
        columns[-1] = window_col
    return columns


def load_cache(cache_file):
    """Load cached records from a file."""
    if os.path.exists(cache_file):
        with open(cache_file, "rb") as f:
            return pickle.load(f)
    return {}

def save_cache(cache_file, records):
    """Save records to a cache file."""
    with open(cache_file, "wb") as f:
        pickle.dump(records, f)

def extract_unique_accessions(blast_folder):
    """Extract unique accession numbers from all BLAST output files."""
    unique_accessions = set()
    blast_outs = [join(blast_folder, f) for f in os.listdir(blast_folder) if f.endswith(".blastout.tsv")]
    
    for blast_out in blast_outs:
        try:
            df = pd.read_csv(blast_out, sep="\t", header=None)
            df.columns = ["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",
                      "q_start", "q_end", "s_start", "s_end", "sstrand", "evalue", "qcovs", "sgi",
                      "sacc", "slen", "staxids", "stitle"]
            unique_accessions.update(df["sacc"].unique())
        except pd.errors.EmptyDataError:
            print(f"File {blast_out} is empty. Skipping...")
    return unique_accessions

def fetch_sequence(seq_id, request_delay, initial_delay_range, entrez_email=None):
    Entrez.email = entrez_email
    # Introduce a small random delay at the start
    initial_delay = random.uniform(*initial_delay_range)
    time.sleep(initial_delay)

    print(f"Fetching sequence {seq_id} after initial delay of {initial_delay:.2f} seconds")
    for i in range(MAX_RETRIES):
        try:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            time.sleep(request_delay/2)  # Delay between successful requests
            return record
        except Exception as e:
            print(f"Error fetching {seq_id}: {e}. Retrying...")
            time.sleep(request_delay)  # Exponential backoff for retries
    return None

def fetch_all_sequences(unique_accessions, cache_file, max_workers=4, entrez_email=None):
    """Fetch sequences for all unique accession numbers using parallel processing."""
    # Load existing cache
    try:
        sacc_records = load_cache(cache_file)
    except Exception as e:
        print("Error loading sacc records: {e}")
        sacc_records = {}

    # Calculate delays based on the desired request rate per worker
    request_delay = (1/3) * max_workers
    initial_delay_range = (request_delay / 2, request_delay)

    def fetch_and_store(seq_id):
        if seq_id not in sacc_records:
            record = fetch_sequence(seq_id, request_delay, initial_delay_range, entrez_email)
            if record:
                sacc_records[seq_id] = record

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        executor.map(fetch_and_store, unique_accessions)

    # Save updated cache
    save_cache(cache_file, sacc_records)
    return sacc_records

def find_overlapping_features(record, window_start, window_end, strand):
    overlapping_features = []
    for feature in record.features:
        if feature.type == "source":
            continue
        feature_start = feature.location.start
        feature_end = feature.location.end

        overlap = (feature_start <= window_end) and (feature_end >= window_start)

        if overlap:
            overlapping_features.append({
                "type": feature.type,
                "start": str(feature_start),
                "end": str(feature_end),
                "gene": feature.qualifiers.get("gene"),
                "product": feature.qualifiers.get("product"),
                "protein_seq": feature.qualifiers.get("translation"),
                "protein_id": feature.qualifiers.get("protein_id"),
                "note": feature.qualifiers.get('note')
            })
    return overlapping_features

def featurize_blast_out(blast_out, window, sacc_records):
    df = pd.read_csv(blast_out, sep="\t", header=None)
    df.columns = ["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",
                  "q_start", "q_end", "s_start", "s_end", "sstrand", "evalue", "qcovs", "sgi",
                  "sacc", "slen", "staxids", "stitle"]
    df["features"] = None
    df[f"features_{window}_window"] = None

    for index, row in df.iterrows():
        record = sacc_records[row["sacc"]]
        features = find_overlapping_features(record, row["s_start"], row["s_end"], row["sstrand"])
        df.at[index, "features"] = features
        window_start = max(row["s_start"] - window, 0)
        window_end = row["s_end"] + window
        strand = row["sstrand"]
        if strand == "plus":
            strand = 1
        elif strand == "minus":
            strand = -1
        features = find_overlapping_features(record, window_start, window_end, strand)
        df.at[index, f"features_{window}_window"] = features
    
    return df[["query", "identity", "evalue", "qcovs", "features", f"features_{window}_window"]]

def process_blast_file(blast_out, blast_feat_out, blast_window, sacc_records):
    if os.path.exists(blast_feat_out) and os.path.getsize(blast_feat_out) > 0:
        print(f"Output file {blast_feat_out} exists and has data. Skipping.")
        return
    if os.path.getsize(blast_out) > 0:
        df_features = featurize_blast_out(blast_out, blast_window, sacc_records)
        df_features.to_csv(blast_feat_out, index=None, sep="\t")
        print(f"Featurize blast output complete for {blast_out}. Output file: {blast_feat_out}")
    else:
        pd.DataFrame(columns=blast_feature_columns(blast_window)).to_csv(
            blast_feat_out, index=None, sep="\t"
        )
        print(f"BLAST output {blast_out} was empty. Wrote header-only feature file: {blast_feat_out}")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Featurize BLAST output")
    parser.add_argument("--blast_folder", required=True, help="Path to the folder containing BLAST output files")
    parser.add_argument("--blast_window", type=int, default=10000, help="Window size for feature extraction")
    parser.add_argument("--output_file", required=True, help="Path to the output file for concatenated results")
    parser.add_argument("--max_workers", type=int, default=4, help="Maximum number of workers for parallel processing")
    parser.add_argument("--entrez_email", required=True, help="Email for NCBI Entrez requests")
    parser.add_argument("--temp_dir", default="tmp/", help="Path to the temporary directory for storing cached sequence records")
    args = parser.parse_args()

    CACHE_FILE = os.path.join(args.temp_dir, "sacc_records.pkl")
    blast_folder = args.blast_folder
    blast_window = args.blast_window
    output_file = args.output_file
    max_workers = args.max_workers

    if not os.path.exists(blast_folder):
        print(f"Blast folder {blast_folder} does not exist. Exiting.")
        sys.exit(0)

    # Extract unique accession numbers
    unique_accessions = extract_unique_accessions(blast_folder)
    print(f"Total unique accessions: {len(unique_accessions)}")

    # Fetch sequences for all unique accessions using parallel processing
    sacc_records = fetch_all_sequences(unique_accessions, CACHE_FILE, max_workers=max_workers, entrez_email=args.entrez_email)

    blast_outs = [join(blast_folder, f) for f in os.listdir(blast_folder) if f.endswith(".blastout.tsv")]
    blast_feat_outs = [join(blast_folder, basename(f).split(".")[0] + ".blastfeatout.tsv") for f in blast_outs]
    print(f"Total blast output files: {len(blast_outs)}")

    with ThreadPoolExecutor(max_workers=max(2, floor(max_workers / 2))) as executor:
        futures = [executor.submit(process_blast_file, blast_out, blast_feat_out, blast_window, sacc_records) for blast_out, blast_feat_out in zip(blast_outs, blast_feat_outs)]
        for future in futures:
            future.result()

    # Concatenate all .blastfeatout.tsv files into the output file
    all_feat_outs = [pd.read_csv(join(blast_folder, f), sep="\t") for f in os.listdir(blast_folder) if f.endswith(".blastfeatout.tsv")]
    if all_feat_outs:
        concatenated_df = pd.concat(all_feat_outs, ignore_index=True)
    else:
        concatenated_df = pd.DataFrame(columns=blast_feature_columns(blast_window))
    concatenated_df.to_csv(output_file, index=None, sep="\t")
    print(f"All .blastfeatout.tsv files have been concatenated into {output_file}")
