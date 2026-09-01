import os
from os.path import join, basename
import subprocess
import shutil
import sys
import argparse
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
import pandas as pd

SPLIT_THRESH = 20  # 100
SPLIT_EACH = 10  # 50
NUMERIC_PATTERN = re.compile(r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?")
CLUSTER_PATTERN = re.compile(r"(^.*cluster_\d+|\w+_kmer_\d+)(?:_|$)")


def get_first_coef_abs(coef_string):
    coefs = NUMERIC_PATTERN.findall(str(coef_string))
    if not coefs:
        return 0.0
    try:
        return abs(float(coefs[0]))
    except ValueError:
        return 0.0


def get_feature_cluster(feature):
    feature = str(feature)
    cluster_match = CLUSTER_PATTERN.search(feature)
    if cluster_match:
        return cluster_match.group(1)
    return feature.rsplit("_", 1)[0]


def get_record_cluster(record):
    return get_feature_cluster(record.id)


def get_plot_selected_clusters(coefficients_file, num_hits):
    coefficients = pd.read_csv(coefficients_file, sep="\t")
    required_columns = {"metadata_category", "feature", "coefficients"}
    missing_columns = required_columns - set(coefficients.columns)
    if missing_columns:
        raise ValueError(
            "Coefficient file is missing required columns for minimal BLAST: "
            + ", ".join(sorted(missing_columns))
        )

    coefficients = coefficients.copy()
    coefficients["cluster"] = coefficients["feature"].apply(get_feature_cluster)
    coefficients["max_coefficient"] = coefficients["coefficients"].apply(get_first_coef_abs)

    selected_clusters = set()
    for _, category_dt in coefficients.groupby("metadata_category", sort=False):
        top_dt = (
            category_dt[["cluster", "max_coefficient"]]
            .drop_duplicates()
            .sort_values("max_coefficient", ascending=False)
            .head(num_hits)
        )
        selected_clusters.update(top_dt["cluster"].astype(str))
    return selected_clusters


def filter_plot_selected_clusters(fasta_file, output_file, coefficients_file, num_hits):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    selected_clusters = get_plot_selected_clusters(coefficients_file, num_hits)
    selected = [
        record
        for record in records
        if get_record_cluster(record) in selected_clusters
    ]
    SeqIO.write(selected, output_file, "fasta")
    print(
        f"Sending {len(selected)}/{len(records)} sequences to BLAST query "
        f"(minimal BLAST: top {num_hits} plotted clusters per metadata category)."
    )
    return output_file


def format_taxids(taxid):
    taxid = str(taxid).strip().strip("\"'").strip()
    if taxid in {"", "0", "NA", "NaN", "nan", "None", "none"}:
        return ""

    taxid = taxid.strip("{}[]()")
    taxids = [
        item.strip().strip("\"'").strip()
        for item in taxid.replace(";", ",").replace("+", ",").split(",")
        if item.strip()
    ]
    if not taxids:
        return ""

    invalid_taxids = [item for item in taxids if not item.isdigit()]
    if invalid_taxids:
        raise ValueError(f"Invalid taxid value(s): {', '.join(invalid_taxids)}")

    return f"-taxids {','.join(taxids)}"


def read_fasta(fasta_file, output_type="dict"):
    """
    Read a fasta file and return a dictionary with the sequence id as key and the sequence as value.
    """
    if output_type == "list":
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record.seq)
        return sequences
    elif output_type == "dict":
        sequences = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = record.seq
        return sequences
    elif output_type == "pandas":
        sequences = []
        description = []
        ids = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record.seq)
            description.append(record.description)
            ids.append(record.id)
        return pd.DataFrame(
            {"ID": ids, "Description": description, "Sequence": sequences}
        )


def split_fasta(fasta_file, output_dir, num_seq=1):
    """
    Split a fasta file into multiple files.
    """
    os.makedirs(output_dir, exist_ok=True)
    if num_seq == 1:
        for record in SeqIO.parse(fasta_file, "fasta"):
            output_file = os.path.join(output_dir, record.id + ".fasta")
            with open(output_file, "w") as f:
                f.write(">" + record.description + "\n")
                f.write(str(record.seq) + "\n")
    else:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        for i in range(0, len(records), num_seq):
            output_file = os.path.join(output_dir, f"split_{i}.fasta")
            with open(output_file, "w") as f:
                for record in records[i : i + num_seq]:
                    f.write(">" + record.description + "\n")
                    f.write(str(record.seq) + "\n")


def run_blast(splitted_fasta, blast_folder, max_workers, taxid, local_blast_db=""):
    fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue qcovs sgi sacc slen staxids stitle"
    taxid = format_taxids(taxid)
    if local_blast_db:
        remote_flag = ""
        local_db_export = f"export BLASTDB={local_blast_db}; "
    else:
        taxid = ""  # cannot use -remote with taxids, so we skip taxid
        remote_flag = "-remote"
        local_db_export = ""

    def run_single_blast(f):
        if taxid:
            print(f"Using taxonomy flag: {taxid}")

        blast_out = join(blast_folder, basename(f).split(".")[0] + ".blastout.tsv")

        # skip if tsv file already exists and is not empty
        if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
            print(f"Skipping {f} as blast output already exists")
            return

        # create the blast command
        params = [
            local_db_export,
            "blastn",
            f"-outfmt '{fmt}'",
            f"-query {f}",
            remote_flag,
            "-db core_nt",
            f"-out {blast_out}",
            "-evalue 0.1",
            "-task blastn",
            "-dust no",
            "-word_size 24",
            "-reward 1",
            "-penalty -3",
            taxid,
            "-max_target_seqs 5",
        ]
        # join only non-empty parts to avoid extra spaces
        cmd = " ".join(p for p in params if p)

        subprocess.run(cmd, shell=True, check=True)
        print(f"Blast complete for {f}")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_single_blast, f) for f in splitted_fasta]
        for future in as_completed(futures):
            future.result()  # Raise any exceptions that occurred


def parse_args():
    parser = argparse.ArgumentParser(description="Run BLAST on input fasta file")
    parser.add_argument(
        "--input_file", required=True, help="Path to the input fasta file"
    )
    parser.add_argument(
        "--split_folder",
        required=True,
        help="Path to the folder to store split fasta files",
    )
    parser.add_argument(
        "--blast_folder", required=True, help="Path to the folder to store BLAST output"
    )
    parser.add_argument(
        "--max_workers", type=int, default=4, help="Number of concurrent BLAST commands"
    )
    parser.add_argument(
        "--taxid",
        type=str,
        default="0",
        help="What tax id to restrict to when searching BLAST",
    )
    parser.add_argument(
        "--local_blast_db",
        type=str,
        default="",
        help="Path to local BLAST database folder (if using local databases)",
    )
    parser.add_argument(
        "--minimal_blast",
        action="store_true",
        help="Only BLAST sequences from clusters that will appear in BLAST plots.",
    )
    parser.add_argument(
        "--coefficients",
        type=str,
        default="",
        help="Nonzero coefficient TSV. Required when --minimal_blast is set.",
    )
    parser.add_argument(
        "--num_plot_hits",
        type=int,
        default=10,
        help="Number of top coefficient clusters per metadata category used by BLAST plots.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.exists(args.blast_folder):
        print("Not running blast as the blast folder does not exist")
        sys.exit(0)
    query_fasta = args.input_file
    if args.minimal_blast:
        if not args.coefficients:
            raise ValueError("--coefficients is required when --minimal_blast is set")
        os.makedirs(args.split_folder, exist_ok=True)
        query_fasta = os.path.join(args.split_folder, "minimal_blast_query.fasta")
        filter_plot_selected_clusters(
            args.input_file, query_fasta, args.coefficients, args.num_plot_hits
        )

    if len(read_fasta(query_fasta)) > SPLIT_THRESH:
        split_fasta(query_fasta, args.split_folder, SPLIT_EACH)
        if query_fasta != args.input_file and os.path.exists(query_fasta):
            os.remove(query_fasta)
    else:
        if os.path.dirname(os.path.abspath(query_fasta)) != os.path.abspath(args.split_folder):
            shutil.copy(query_fasta, args.split_folder)
    splitted_fasta = [
        join(args.split_folder, f)
        for f in os.listdir(args.split_folder)
        if f.endswith(".fasta")
    ]
    run_blast(
        splitted_fasta,
        args.blast_folder,
        args.max_workers,
        args.taxid,
        args.local_blast_db,
    )
