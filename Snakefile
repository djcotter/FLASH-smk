"""
FLASH Pipeline

Daniel Cotter
Update June 17 2025

This pipeline is designed to take SPLASH results and use them to predict on phenotypic metadata.
The pipeline is written in snakemake and is best run using resources from a high-performance computing cluster.
"""

## Importing necessary modules
from pathlib import Path
import pandas as pd
import csv

## Define the config files for the pipeline
# you can define different config files for different runs of the pipeline
configfile: "config.yaml" 
dataset_table_path = Path(config["data_table"])
TEMP_DIR = Path(config["temp_dir"])

## read the dataset table with the short names as indices
## datasets here define the target of the pipeline
## NOTE: If a dataset is the result of a scSPLASH run on 10X data, it should contain "SC10X" in its short name
dataset_table = pd.read_csv(dataset_table_path, index_col = "dataset_short_name")

## Define the wildcards on which the pipeline will be run
DATASETS = list(dataset_table.index)
SELECT_TYPES = config["options"]["filters"]
CLUSTER_TYPES = config["options"]["cluster_types"]
NUM_CLUSTERS = [config["options"]["num_clusters"]]
ANCHOR_LENGTH = config["options"]["anchor_length"] # this is the length of the anchor in nucleotides
TARGET_LENGTH = config["options"]["target_length"] # this is the length of the target in nucleotides
KMER_WIDTH = [ANCHOR_LENGTH + TARGET_LENGTH] # this is Anchor Length + Target Length
KMER_STEP = [ANCHOR_LENGTH + TARGET_LENGTH] # this can be used to let the steps be variable, but for now we will use a single value
MODELS = config["options"]["models"]
NORMALIZE = config["options"]["normalize_embeddings"]
TRAIN_PROPORTION = config["options"]["train_proportion"] # this is the proportion of the data to use for training, the rest will be used for testing
raw_target_rank = config["options"]["target_rank"] # this is the rank of the target to use for prediction (1 = top target, 2 = second target, etc.)
TARGET_RANK = raw_target_rank if isinstance(raw_target_rank, (list, tuple)) else [raw_target_rank]

# whether to generate plots or to stop at the output of the prediction task
GENERATE_PLOTS = config["options"]["generate_plots"]
if GENERATE_PLOTS:
    FILE_SUFFIXES = ["nonzero_coefficients_annotated.tsv", "confusion_matrices.pdf",
                    "nonzero_coefficients_blast_annotated_plots.pdf", "nonzero_coefficients_heatmaps.pdf"]
else:
    FILE_SUFFIXES = ["nonzero_coefficients_annotated.tsv", "confusion_matrices.pdf"]


## constrain the wildcards of the pipeline
# specifically we do not want wildcards to contain underscores or spaces as they are used to 
# separate the wildcards in the output file names
wildcard_constraints:
    dataset=r"[A-Za-z\d-]+",
    model=r"[A-Za-z\d-]+",
    select_type=r"[A-Za-z\d-]+",
    cluster_type=r"[A-Za-z\d-]+",
    num_clusters=r"\d+",
    kmer_width=r"\d+",
    kmer_step=r"\d+",
    normalize=r"[A-Za-z]+"

## TARGET RULES --------------------------------
## Define Rule all for the target of the pipeline
rule all_embeddings:
    """
    Generate all summary files for the datasets defined in the dataset table
    """
    input:
        expand(Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", 
                    "{dataset}_{model}_adelie_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_{FILE}"),
               dataset=DATASETS,
               select_type=SELECT_TYPES,
               cluster_type=CLUSTER_TYPES,
               model=MODELS,
               num_clusters=NUM_CLUSTERS,
               target_rank=TARGET_RANK,
               kmer_width=KMER_WIDTH,
               kmer_step=KMER_STEP,
               normalize=NORMALIZE,
               train_proportion=TRAIN_PROPORTION,
               FILE = FILE_SUFFIXES)


rule all_genomes:
    """
    Generate prediction for genome sequences provided for a given dataset.
    """
    input:
        # all genome coefficients files
        expand(Path("results",
                    "{dataset}", 
                    "{select_type}", 
                    "{cluster_type}", 
                    "{model}", 
                    "genomes", 
                    "{normalize}", 
                    "{dataset}_{model}_adelie_genomes_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_{FILE}"),
               dataset=DATASETS,
               select_type=SELECT_TYPES,
               cluster_type=CLUSTER_TYPES,
               model=MODELS,
               num_clusters=NUM_CLUSTERS,
               kmer_width=KMER_WIDTH,
               kmer_step=KMER_STEP,
               train_proportion=TRAIN_PROPORTION,
               normalize=NORMALIZE,
               FILE = ["nonzero_coefficients_annotated.tsv", "confusion_matrices.pdf"]) 
               # no need to plot genomes blast results as they
               # will be the same as the non-genome results which generated the model for the prediction.


rule all_ohe:
    """
    Generate all one-hot encoded files for the datasets defined in the dataset table and perform
    prediction on the metadata.
    """
    input:
        expand(Path("results", "{dataset}", "{select_type}", "{cluster_type}", "ohe", 
                    "{dataset}_ohe_adelie_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_{FILE}"),
               dataset=DATASETS,
               select_type=SELECT_TYPES,
               cluster_type=CLUSTER_TYPES,
               num_clusters=NUM_CLUSTERS,
               target_rank=TARGET_RANK,
               kmer_width=KMER_WIDTH,
               kmer_step=KMER_STEP,
               train_proportion=TRAIN_PROPORTION,
               FILE = FILE_SUFFIXES)


rule all_umap:
    """
    Generate UMAP plots for the datasets defined in the dataset table. Assumes no metadata is available and colors plots by HDBSCAN cluster ID.
    The script can be modified slightly to color by metadata if available.
    """
    input: 
        expand(Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", 
                    "{dataset}_{model}_umap_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.pdf"),
               dataset=DATASETS,
               select_type=SELECT_TYPES,
               cluster_type=CLUSTER_TYPES,
               model=MODELS,
               num_clusters=NUM_CLUSTERS,
               target_rank=TARGET_RANK,
               kmer_width=KMER_WIDTH,
               kmer_step=KMER_STEP,
               normalize=NORMALIZE)


## STANDARD PIPLELINE RULES -------------------------------
# Begin Standard FLASH Pipeline Rules
rule choose_anchors:
    """
    This rule selects the top anchors from the SPLASH results based on various criteria (select_type)
    New scripts that select different anchors can be added to the config file under the "anchor_select_script" key
    """
    input:
        lambda wildcards: Path(dataset_table.loc[wildcards.dataset, "SPLASH_results"],
                               "result.after_correction.scores.tsv") # this is the path to the default SPLASH results file
    params:
        script = lambda wildcards: Path(config["scripts"]["anchor_select_script"][wildcards.select_type]),
        lookup_table = config["lookup_table_for_artifact_filtering"],
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type),
        splash_bin = config["splash_bin"],
        num_anchors = config["extended_options"]["num_anchors_to_select"],
        effect_size = config["extended_options"]["effect_size_cutoff"] # only select anchors with an effect size greater than this value
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_selected_anchors_{select_type}.txt")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --input {input} --output {output} \
        --lookup_table {params.lookup_table} --temp_dir {params.tmp_dir} --num_anchors {params.num_anchors} \
        --effect_size {params.effect_size} --splash_bin {params.splash_bin}
    """


rule cluster_anchors:
    """
    This rule clusters the selected anchors based on the selected clustering method (cluster_type)
    New scripts that cluster the anchors can be added to the config file under the "cluster_script" key
    This rule has a mix of R and Python code, so it is necessary to load all appropriate modules
    The call to the python or Rscript is stored in the config file under the "cluster_script" key
    """
    input:
        Path(TEMP_DIR, "{dataset}", "{dataset}_selected_anchors_{select_type}.txt")
    params:
        script = lambda wildcards: Path(config["scripts"]["cluster_script"][wildcards.cluster_type]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type),
        translation_table = lambda wildcards: f"--translation_table {dataset_table.loc[wildcards.dataset, "translation_table"]}"  if wildcards.cluster_type == "masked-aa-clustered" else "",
        fuzzy_params = lambda wildcards: (
            config["scripts"]["clustering_params"]["fuzzy_nt_clustering_params"]
            if wildcards.cluster_type == "masked-nucleotide-clustered"
            else (
            config["scripts"]["clustering_params"]["fuzzy_aa_clustering_params"]
            if wildcards.cluster_type == "masked-aa-clustered"
            else ""
            )
        ),
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_clustered_anchors_{select_type}_{cluster_type}.txt")
    conda:
        config["envs"]["biopython_env"]
    shell:"""
        {params.script} --input {input} --output {output} --temp_dir {params.tmp_dir} {params.translation_table} {params.fuzzy_params}
    """


rule reorder_clusters:
    """
    This rule reorders the clusters based on the SPLASH results and the selected anchors and clustering method
    One reordering example is to sort the clusters and only grab 1 anchor per cluster.
    New scripts that reorder the clusters can be added to the config file under the "reorder_script" key
    """
    input:
        clusters = Path(TEMP_DIR, "{dataset}", "{dataset}_clustered_anchors_{select_type}_{cluster_type}.txt"),
        splash_results = lambda wildcards: Path(dataset_table.loc[wildcards.dataset, "SPLASH_results"],
                                                "result.after_correction.scores.tsv")
    params:
        script = lambda wildcards: Path(config["scripts"]["reorder_script"][wildcards.cluster_type]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type),
        effect_size_threshold = config["extended_options"]["effect_size_cutoff"],
        distance_threshold = config["extended_options"]["distance_threshold"],
        max_clusters_to_process = int(config["options"]["num_clusters"] * 1.5) # process 1.5x the number of clusters to select from
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_reordered_clusters_{select_type}_{cluster_type}.txt")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --input_anchor_clusters {input.clusters} \
        --splash_stats {input.splash_results} --output {output} --temp_dir {params.tmp_dir} \
        --effect_size_cutoff {params.effect_size_threshold} --distance_threshold {params.distance_threshold} \
        --max_clusters_to_process {params.max_clusters_to_process} --num_cores {threads}
    """


rule select_N_clusters:
    """
    After reordering the clusters, this rule selects the top N clusters based on their id (since they are already reordered)
    """
    input:
        Path(TEMP_DIR, "{dataset}", "{dataset}_reordered_clusters_{select_type}_{cluster_type}.txt")
    output:
        clusters = Path(TEMP_DIR, "{dataset}", "{dataset}_selected_clusters_{select_type}_{cluster_type}_top{num_clusters}-clusters.txt"),
        anchors = Path(TEMP_DIR, "{dataset}", "{dataset}_selected_anchors_{select_type}_{cluster_type}_top{num_clusters}-clusters.txt")
    shell:"""
        awk '$1 <= {wildcards.num_clusters}' {input} > {output.clusters}
        cut -f2 {output.clusters} > {output.anchors}
    """


rule prepare_sequences:
    """
    This rule prepares the sequences for the selected clusters and anchors for the dataset.
    It formats the sequences into concatmers for each sample where all anchor-target pairs are concatenated
    together. The output is a fasta file and a tsv file with the sequences.
    """
    input:
        cluster_file = Path(TEMP_DIR, "{dataset}", "{dataset}_selected_clusters_{select_type}_{cluster_type}_top{num_clusters}-clusters.txt"),
        anchor_file = Path(TEMP_DIR, "{dataset}", "{dataset}_selected_anchors_{select_type}_{cluster_type}_top{num_clusters}-clusters.txt"),
        id_mapping = lambda wildcards: Path(dataset_table.loc[wildcards.dataset, "SPLASH_results"], "sample_name_to_id.mapping.txt"),
        sample_sheet = lambda wildcards: Path(dataset_table.loc[wildcards.dataset, "SPLASH_results"], "sample_sheet.txt")
    params:
        script = lambda wildcards: (
            Path(config["scripts"]["prepare_sequences_single_cell"])
            if ("SC10X" in wildcards.dataset or config["options"]["seqs_from_raw_data"] == False)
            else Path(config["scripts"]["prepare_sequences"])
        ),
        input_samples = lambda wildcards: (
            f"--satc_files {Path(dataset_table.loc[wildcards.dataset, 'SPLASH_results'], 'result_satc')}"
            if ("SC10X" in wildcards.dataset or config["options"]["seqs_from_raw_data"] == False)
            else f"--sample_sheet {Path(dataset_table.loc[wildcards.dataset, 'SPLASH_results'], 'sample_sheet.txt')}"
        ),
        output_prefix = lambda wildcards: Path(
            TEMP_DIR,
            f"{wildcards.dataset}",
            f"{wildcards.dataset}_prepared_sequences_{wildcards.select_type}_{wildcards.cluster_type}_top{wildcards.num_clusters}_target{wildcards.target_rank}_k{wildcards.kmer_width}_s{wildcards.kmer_step}"
        ),
        tmp_dir = lambda wildcards: Path(
            TEMP_DIR,
            wildcards.dataset,
            wildcards.select_type,
            wildcards.cluster_type,
            "target" + wildcards.target_rank,
            wildcards.num_clusters + "-clusters"
        ),
        single_cell = lambda wildcards: (
            "--single_cell" if "SC10X" in wildcards.dataset else ""
        ),
        target_rank = lambda wildcards: wildcards.target_rank,
        cluster_filter = lambda wildcards: (
            f"--apply_cluster_filter {config['options']['cluster_filter']['type']}:{config['options']['cluster_filter']['threshold']}"
            if config["options"]["cluster_filter"]["apply"] else ""
        ),
        satc_util_bin = Path(config["satc_util_bin"]) # satc utils are typically located in the same directory as the main splash binary.
    output:
        fasta = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_sample_sequences.fasta"),
        tsv = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_sample_sequences.tsv")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --anchor_file {input.anchor_file} \
        --cluster_file {input.cluster_file} --id_mapping {input.id_mapping} \
        {params.input_samples} --output_prefix {params.output_prefix} \
        --temp_dir {params.tmp_dir} --num_cores {threads} {params.single_cell} \
        --target_rank {params.target_rank} {params.cluster_filter} \
        --satc_util_bin {params.satc_util_bin}
    """


rule decompose_kmers:
    """
    Process the sample sequences to decompose them into kmers of width kmer_width and step kmer_step
    The outputs are 1) a fasta file with the unique kmers and 2) a tsv file with the ordering of the kmers
    for each sample
    """
    input:
        Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_sample_sequences.fasta")
    params:
        script = Path(config["scripts"]["decompose_kmers"]),
        output_prefix = lambda wildcards: Path(TEMP_DIR, f"{wildcards.dataset}", f"{wildcards.dataset}_decomposed_kmers_{wildcards.select_type}_{wildcards.cluster_type}_top{wildcards.num_clusters}_target{wildcards.target_rank}"),
        kmer_width = lambda wildcards: wildcards.kmer_width,
        kmer_step = lambda wildcards: wildcards.kmer_step
    output:
        unique_kmers = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta"),
        order = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv")
    conda:
        config["envs"]["biopython_env"]
    shell:"""
        python {params.script} -k {wildcards.kmer_width} -s {wildcards.kmer_step} \
        {input} {params.output_prefix}
    """


rule match_kmers_to_clusters:
    """
    Process the ordering file to produce a tsv file mapping clusters to
    their component kmers
    """
    input:
        order = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv"),
        unique_kmers = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta")
    params:
        script = Path(config["scripts"]["match_kmers_to_clusters"])
    output:
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv")
    conda:
        config["envs"]["default_r"]
    shell:"""
    Rscript --vanilla {params.script} --ordering {input.order} --kmers {input.unique_kmers} --output {output}
    """


rule embed_kmers_hyena:
    """
    Embeds the unique kmers using the Hyena model we have pretrained.
    TODO: Make this rule replaceable with other embedding models in the future
    """
    input:
        kmers = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta"),
        singularity_image = config["containers"]["hyena_embedder"]
    params:
        wrapper_script = config["scripts"]["embed_kmers_hyena"],
        python_embedder = config["scripts"]["hyena_kmer_embedder"]
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_hyena-embeddings_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv")
    shell:"""
        bash {params.wrapper_script} {params.python_embedder} {input.singularity_image} {input.kmers} {output}
    """


rule prepare_data_for_prediction_top_variance:
    """
    This rule processes the embeddings to fit into a glmnet model by grabbing the top variance embeddings per cluster
    and then saves the resulting data frame as a feather object to be used in the glmnet model.
    """
    input:
        embeddings = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_{model}-embeddings_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv"),
        ordering = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv")
    params:
        script = Path(config["scripts"]["format_embeddings_variance"]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.num_clusters + "-clusters", "target" + wildcards.target_rank,
                                         "k" + wildcards.kmer_width + "_s" + wildcards.kmer_step, wildcards.model + "_embeddings", wildcards.normalize),
        normalized_flag = lambda wildcards: "--normalized_embeddings" if wildcards.normalize =="normalized" else "",
        num_to_keep = config["extended_options"]["num_embedding_features_to_keep"]["by_variance"],
        recode_missing_flag = "--recode_missing" if config["options"]["feature_processing"]["recode_missing"] else ""
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_{normalize}.feather")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --embeddings {input.embeddings} --ordering {input.ordering} \
        --output {output} --temp_dir {params.tmp_dir} --num_threads {threads} --num_to_keep {params.num_to_keep} \
        {params.normalized_flag} {params.recode_missing_flag}
    """


rule prepare_data_for_prediction_pca:
    """
    This rule processes the embeddings to fit into a glmnet model by performing PCA on the embeddings by 
    samples matrix per cluster and then saves the resulting data frame as a feather object to be used in the glmnet model.
    The number of PCA components to keep is defined in the config file under "num_embedding_features_to_keep" -> "by_pca"
    """
    input: 
        embeddings = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_{model}-embeddings_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv"),
        ordering = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv")
    params:
        script = Path(config["scripts"]["format_embeddings_pca"]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.num_clusters + "-clusters", "target" + wildcards.target_rank,
                                         "k" + wildcards.kmer_width + "_s" + wildcards.kmer_step, wildcards.model + "_pca_embeddings", wildcards.normalize),
        normalized_flag = lambda wildcards: "--normalized" if wildcards.normalize =="normalized" else "",
        num_pcs = config["extended_options"]["num_embedding_features_to_keep"]["by_pca"]
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_pca_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_{normalize}.feather")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --embeddings {input.embeddings} \
        --ordering {input.ordering} \
        --output {output} \
        --temp_dir {params.tmp_dir} \
        --num_threads {threads} \
        --num_pcs {params.num_pcs} \
        {params.normalized_flag}
    """


rule prepare_data_for_umap_top_variance:
    """
    This rule processes the embeddings to fit into a glmnet model by grabbing the top variance embeddings per cluster
    and then saves the resulting data frame as a feather object to be used in the glmnet model.
    """
    input:
        embeddings = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_{model}-embeddings_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv"),
        ordering = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv")
    params:
        script = Path(config["scripts"]["format_embeddings_variance"]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.num_clusters + "-clusters", "target" + wildcards.target_rank,
                                         "k" + wildcards.kmer_width + "_s" + wildcards.kmer_step, wildcards.model + "_umap_embeddings", wildcards.normalize),
        normalized_flag = lambda wildcards: "--normalized_embeddings" if wildcards.normalize =="normalized" else "",
        num_to_keep = config["extended_options"]["num_embedding_features_to_keep"]["umap"]
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_umap_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_{normalize}.feather")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --embeddings {input.embeddings} --ordering {input.ordering} \
        --output {output} --temp_dir {params.tmp_dir} --num_threads {threads} --num_to_keep {params.num_to_keep} {params.normalized_flag}
    """


rule run_adelie:
    """
    This rule uses preprocessed embeddings for running the glmnet model to predict on the metadata.
    """
    input:
        embeddings = lambda wildcards: (
            Path(
                TEMP_DIR,
                "{dataset}",
                "{dataset}_{model}_pca_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_{normalize}.feather",
            )
            if config["options"]["feature_processing"]["method"] == "pca"
            else Path(
                TEMP_DIR,
                "{dataset}",
                "{dataset}_{model}_top_variance_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_{normalize}.feather",
            )
            if config["options"]["feature_processing"]["method"] == "top_variance"
            else "",
        ),
        metadata = lambda wildcards: dataset_table.loc[wildcards.dataset, "metadata_file"]
    params:
        script = Path(config["scripts"]["adelie"]),
        output_prefix = lambda wildcards: Path("results", f"{wildcards.dataset}", f"{wildcards.select_type}", f"{wildcards.cluster_type}", f"{wildcards.model}", f"{wildcards.normalize}", f"{wildcards.dataset}_{wildcards.model}_adelie_results_top{wildcards.num_clusters}_target{wildcards.target_rank}_k{wildcards.kmer_width}_s{wildcards.kmer_step}_trainProp{wildcards.train_proportion}"),
        min_samples = config["extended_options"]["min_samples_adelie"],
        grouped_flag = "--grouped" if config["options"]["grouped_model"] else "",
        alpha = config["extended_options"]["adelie_alpha"]
    output:
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_adelie_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv"),
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_adelie_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_confusion_matrices.pdf"),
    conda:
        config["envs"]["adelie_env"]
    shell:"""
        python {params.script} --data {input.embeddings} \
        --metadata {input.metadata} --output_prefix {params.output_prefix} \
        --min_samples {params.min_samples} \
        --n_threads {threads} --train_prop {wildcards.train_proportion} \
        --alpha {params.alpha} \
        {params.grouped_flag}
    """


rule prepare_data_for_prediction_ohe:
    input:
        sample_sequences = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_sample_sequences.fasta")
    params:
        script = Path(config["scripts"]["format_sequences_ohe"]),
        kmer_width = lambda wildcards: wildcards.kmer_width,
        kmer_step = lambda wildcards: wildcards.kmer_step
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_ohe_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.feather")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --input {input} --output {output} --kmer_width {params.kmer_width}
    """


rule run_adelie_ohe:
    input:
        features = Path(TEMP_DIR, "{dataset}", "{dataset}_ohe_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.feather"),
        metadata = lambda wildcards: dataset_table.loc[wildcards.dataset, "metadata_file"]
    params:
        script = Path(config["scripts"]["adelie"]),
        output_prefix = lambda wildcards: Path("results", f"{wildcards.dataset}", f"{wildcards.select_type}", f"{wildcards.cluster_type}", f"ohe", f"{wildcards.dataset}_ohe_adelie_results_top{wildcards.num_clusters}_target{wildcards.target_rank}_k{wildcards.kmer_width}_s{wildcards.kmer_step}_trainProp{wildcards.train_proportion}"),
        min_samples = config["extended_options"]["min_samples_adelie"],
        grouped_flag = "--grouped" if config["options"]["grouped_model"] else "",
        alpha = config["extended_options"]["adelie_alpha"]
    output:
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_adelie_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv"),
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_adelie_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_confusion_matrices.pdf"),
    conda:
        config["envs"]["adelie_env"]
    shell:"""
        python {params.script} --data {input.features} \
        --metadata {input.metadata} --output_prefix {params.output_prefix} \
        --min_samples {params.min_samples} \
        --n_threads {threads} --train_prop {wildcards.train_proportion} \
        --alpha {params.alpha} \
        {params.grouped_flag}
    """


## GENOME RULES --------------------------------
## the following rules are for  genome sequences
rule process_genome_to_sample_sequences:
    """
    This rule processes the genome sequences to produce a test set of sample sequences
    Takes in a clustered set of anchors and produces a fasta file with the sample sequences
    """
    input:
        cluster_file = Path(TEMP_DIR, "{dataset}", "{dataset}_selected_clusters_{select_type}_{cluster_type}_top{num_clusters}-clusters.txt"),
        genome_list = lambda wildcards: dataset_table.loc[wildcards.dataset, "genome_list"], # path to the file containing the list of genomes to process
        genome_files = lambda wildcards: dataset_table.loc[wildcards.dataset, "genomes_folder"] # path to the folder containing the genome files
    params:
        script = Path(config["scripts"]["genome_to_sample_sequences"]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.num_clusters + "-clusters", "genomes"),
        output_prefix = lambda wildcards: Path(TEMP_DIR, f"{wildcards.dataset}", f"{wildcards.dataset}_prepared_sequences_{wildcards.select_type}_{wildcards.cluster_type}_top{wildcards.num_clusters}_genomes"),
        satc_util_bin = Path(config["satc_util_bin"]) # satc utils are typically located in the same directory as the main splash binary.
    output:
        fasta = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_genomes_sample_sequences.fasta"),
        tsv = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_genomes_sample_sequences.tsv")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --cluster_file {input.cluster_file} \
        --genome_list {input.genome_list} --genome_files {input.genome_files} \
        --output_prefix {params.output_prefix} --temp_dir {params.tmp_dir} \
        --num_cores {threads} --satc_util_bin {params.satc_util_bin}
    """


rule decompose_kmers_genomes:
    """
    Process the genome sample sequences to decompose them into kmers of width kmer_width and step kmer_step.
    The outputs are 1) a fasta file with the unique kmers and 2) a tsv file with the ordering of the kmers
    for each sample.
    """
    input:
        Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_genomes_sample_sequences.fasta")
    params:
        script = Path(config["scripts"]["decompose_kmers"]),
        output_prefix = lambda wildcards: Path(TEMP_DIR, f"{wildcards.dataset}", f"{wildcards.dataset}_decomposed_kmers_genomes_{wildcards.select_type}_{wildcards.cluster_type}_top{wildcards.num_clusters}"),
        kmer_width = lambda wildcards: wildcards.kmer_width,
        kmer_step = lambda wildcards: wildcards.kmer_step
    output:
        unique_kmers = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_genomes_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta"),
        order = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_genomes_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv")
    conda:
        config["envs"]["biopython_env"]
    shell:"""
        python {params.script} -k {wildcards.kmer_width} -s {wildcards.kmer_step} \
        {input} {params.output_prefix}
    """


rule embed_kmers_hyena_genomes:
    """
    Embeds the unique kmers using the Hyena model we have pretrained.
    TODO: Make this rule replaceable with other embedding models in the future
    """
    input:
        kmers = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_genomes_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta"),
        singularity_image = config["containers"]["hyena_embedder"]
    params:
        wrapper_script = config["scripts"]["embed_kmers_hyena"],
        python_embedder = config["scripts"]["hyena_kmer_embedder"]
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_hyena-embeddings_genomes_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}.tsv")
    shell:"""
        bash {params.wrapper_script} {params.python_embedder} {input.singularity_image} {input.kmers} {output}
    """


rule prepare_data_for_prediction_genomes:
    """
    Use a specific script to take in the feather object for the original embeddings and process the 
    genome embeddings to produce a matching feather object for the glmnet model
    """
    input:
        embeddings = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_{model}-embeddings_genomes_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}.tsv"),
        ordering = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_genomes_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv"),
        original_embeddings_feather = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target1_k{kmer_width}_s{kmer_step}_{normalize}.feather"),
        nonzero_coefficients_file = Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_adelie_results_top{num_clusters}_target1_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv")
    params:
        script = Path(config["scripts"]["format_embeddings_genomes"]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.num_clusters + "-clusters",
                                         "k" + wildcards.kmer_width + "_s" + wildcards.kmer_step, wildcards.model + "_embeddings", wildcards.normalize),
        normalized_flag = lambda wildcards: "--normalized_embeddings" if wildcards.normalize =="normalized" else "",
        output_prefix = lambda wildcards: Path(TEMP_DIR, f"{wildcards.dataset}", f"{wildcards.dataset}_{wildcards.model}_top_variance_features_for_glmnet_genomes_{wildcards.select_type}_{wildcards.cluster_type}_top{wildcards.num_clusters}_target1_k{wildcards.kmer_width}_s{wildcards.kmer_step}_trainProp{wildcards.train_proportion}_{wildcards.normalize}")
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_genomes_{select_type}_{cluster_type}_top{num_clusters}_target1_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_{normalize}_ORIGINAL.feather"),
        Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_genomes_{select_type}_{cluster_type}_top{num_clusters}_target1_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_{normalize}_GENOMES.feather")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --embeddings {input.embeddings} --ordering {input.ordering} --original_feather {input.original_embeddings_feather} \
        --original_adelie_output {input.nonzero_coefficients_file} \
        --output_prefix {params.output_prefix} --temp_dir {params.tmp_dir} --num_threads {threads} {params.normalized_flag} 
    """


rule run_adelie_genomes:
    """
    Take in main embeddings as train and genome embeddings as test plus both of their metadata files
    and run the glmnet script
    """
    input:
        train_features = Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_genomes_{select_type}_{cluster_type}_top{num_clusters}_target1_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_{normalize}_ORIGINAL.feather"),
        train_metadata = lambda wildcards: dataset_table.loc[wildcards.dataset, "metadata_file"],
        test_features = Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_genomes_{select_type}_{cluster_type}_top{num_clusters}_target1_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_{normalize}_GENOMES.feather"),
        test_metadata = lambda wildcards: dataset_table.loc[wildcards.dataset, "genome_metadata_file"]
    params:
        script = Path(config["scripts"]["glmnet_genomes_script"]),
        output_prefix = lambda wildcards: Path("results", f"{wildcards.dataset}", f"{wildcards.select_type}", f"{wildcards.cluster_type}", f"{wildcards.model}", "genomes", f"{wildcards.normalize}", f"{wildcards.dataset}_{wildcards.model}_adelie_genomes_results_top{wildcards.num_clusters}_k{wildcards.kmer_width}_s{wildcards.kmer_step}_trainProp{wildcards.train_proportion}"),
        min_samples = config["extended_options"]["min_samples_adelie"],
        alpha = config["extended_options"]["adelie_alpha"],
        grouped_flag = "--grouped" if config["options"]["grouped_model"] else ""
    output:
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "genomes", "{normalize}", "{dataset}_{model}_adelie_genomes_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv"),
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "genomes", "{normalize}", "{dataset}_{model}_adelie_genomes_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_confusion_matrices.pdf"),
    conda:
        config["envs"]["adelie_env"]
    shell:"""
        python {params.script} --train_features {input.train_features} --train_metadata {input.train_metadata} \
        --test_features {input.test_features} --test_metadata {input.test_metadata} --output_prefix {params.output_prefix} \
        --n_threads {threads} \
        --min_samples {params.min_samples} \
        --alpha {params.alpha} \
        {params.grouped_flag}
    """


## ANNOTATION RULES --------------------------------
# These rules are for annotating the clusters with the lookup table and merging the annotations
# They also include rules for annotating the clusters with blast searches

rule annotate_clusters:
    """
    Annotate the clusters with the provided lookup table
    Currently, this lookup table includes annotations for antibiotic resistance genes
    """
    input:
        cluster_seqs = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv"),
        lookup_table = config["lookup_table_for_annotation"]
    params:
        script = config["scripts"]["annotate_clusters"],
        temp_dir = lambda wildcards: Path(TEMP_DIR, f"{wildcards.dataset}", f"{wildcards.select_type}", f"{wildcards.cluster_type}", f"{wildcards.num_clusters}-clusters", f"target{wildcards.target_rank}_k{wildcards.kmer_width}_s{wildcards.kmer_step}"),
        splash_bin = config["splash_bin"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target{target_rank}_k{kmer_width}_s{kmer_step}_annotated.tsv")
    conda:
        config["envs"]["biopython_env"]
    shell:"""
        python {params.script} --cluster_seqs {input.cluster_seqs} --lookup_table {input.lookup_table} \
        --output {output} --temp_dir {params.temp_dir} --splash_bin {params.splash_bin} 
    """


rule merge_annotations:
    """
    Merge the annotations for each cluster onto the non-zero coefficients file
    """
    input:
        annotations = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target{target_rank}_k{kmer_width}_s{kmer_step}_annotated.tsv"),
        coefficients = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv")
    params:
        script = config["scripts"]["merge_annotations"]
    output:
        coefs = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_annotated.tsv"),
        fasta = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences.fasta")
    conda:
        config["envs"]["default_r"]
    shell:"""
    Rscript --vanilla {params.script} --annotations {input.annotations} --coefficients {input.coefficients} --output {output.coefs}
    """


rule run_blast_nonzero_features:
    """
    Run a blast search on the significant sequences to find the closest matches in the NCBI database
    """
    input:
        fasta = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences.fasta")
    params:
        script = lambda wildcards: config["scripts"]["run_blast"],
        blast_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.model, wildcards.num_clusters, "target" + wildcards.target_rank, wildcards.normalize, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "blast"),
        split_fasta_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.model, wildcards.num_clusters, "target" + wildcards.target_rank, wildcards.normalize, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "split_fasta"),
        taxid = lambda wildcards: int(dataset_table.loc[wildcards.dataset, "taxid"]),
        entrez_email = config["entrez_email"] if config["entrez_email"] else 0,
        temp_dir = config["temp_dir"],
        blast_db_path = config["blast_db_path"] if config["blast_db_path"] else 0
    output:
        Path(TEMP_DIR, 
             "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{num_clusters}", "target" + "{target_rank}", "{normalize}", "{predictionTask}", "trainProp" + "{train_proportion}",
             "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_blast.tsv")
    conda:
        config["envs"]["biopython_env_r"]
    shell:"""
        bash {params.script} {input.fasta} {params.split_fasta_temp_dir} {params.blast_temp_dir} {output} {threads} {params.taxid} {params.entrez_email} {params.temp_dir} {params.blast_db_path}
    """


rule run_blastp_nonzero_features:
    """
    Run a blast search on the significant sequences to find the closest matches in the NCBI database
    """
    input:
        fasta = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences.fasta")
    params:
        script = lambda wildcards: config["scripts"]["run_blastp"],
        blast_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.model, wildcards.num_clusters, "target" + wildcards.target_rank, wildcards.normalize, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "blastp"),
        split_fasta_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.model, wildcards.num_clusters, "target" + wildcards.target_rank, wildcards.normalize, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "split_fasta_blastp"),
        taxid = lambda wildcards: int(dataset_table.loc[wildcards.dataset, "taxid"]),
        translation_table = lambda wildcards: dataset_table.loc[wildcards.dataset, "translation_table"],
        blast_db_path = config["blast_db_path"] if config["blast_db_path"] else 0
    output:
        Path(TEMP_DIR, 
             "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{num_clusters}", "target" + "{target_rank}", "{normalize}", "{predictionTask}", "trainProp" + "{train_proportion}",
             "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_blastp.tsv")
    conda:
        config["envs"]["biopython_env_r"]
    shell:"""
        bash {params.script} {input.fasta} {params.split_fasta_temp_dir} {params.blast_temp_dir} {output} {threads} {params.taxid} {params.translation_table} {params.blast_db_path}
    """


rule run_blastp_swissprot_nonzero_features:
    """
    Run a blast search on the significant sequences to find the closest matches in the NCBI database
    """
    input:
        fasta = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences.fasta")
    params:
        script = lambda wildcards: config["scripts"]["run_blastp"],
        blast_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.model, wildcards.num_clusters, "target" + wildcards.target_rank, wildcards.normalize, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "blastp_swissprot"),
        split_fasta_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.model, wildcards.num_clusters, "target" + wildcards.target_rank, wildcards.normalize, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "split_fasta_blastp_swissprot"),
        taxid = lambda wildcards: int(dataset_table.loc[wildcards.dataset, "taxid"]),
        translation_table = lambda wildcards: dataset_table.loc[wildcards.dataset, "translation_table"],
        blast_db_path = config["blast_db_path"] if config["blast_db_path"] else 0,
        protein_db = "swissprot"
    output:
        Path(TEMP_DIR, 
             "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{num_clusters}", "target" + "{target_rank}", "{normalize}", "{predictionTask}", "trainProp" + "{train_proportion}",
             "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_swissprot.tsv")
    conda:
        config["envs"]["biopython_env_r"]
    shell:"""
        bash {params.script} {input.fasta} {params.split_fasta_temp_dir} {params.blast_temp_dir} {output} {threads} {params.taxid} {params.translation_table} {params.blast_db_path} {params.protein_db}
    """


rule merge_blast_results:
    """
    Merge the blast results with the annotated sequences
    """
    input:
        blast_annotations = Path(TEMP_DIR, "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{num_clusters}", "target{target_rank}", "{normalize}", "{predictionTask}", "trainProp" + "{train_proportion}",
                                 "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_blast.tsv"),
        coefficients = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv")
    params:
        script = config["scripts"]["merge_blast_results"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blast_annotated.tsv")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --blast_annotations {input.blast_annotations} --coefficients {input.coefficients} --output {output}
    """


rule merge_blastp_results:
    """
    Merge the blast results with the annotated sequences
    """
    input:
        blast_annotations = Path(TEMP_DIR, "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{num_clusters}", "target{target_rank}", "{normalize}", "{predictionTask}", "trainProp" + "{train_proportion}",
                                 "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_blastp.tsv"),
        coefficients = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv")
    params:
        script = config["scripts"]["merge_blast_results"],
        translation_table = lambda wildcards: dataset_table.loc[wildcards.dataset, "translation_table"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blastp_annotated.tsv")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --blast_annotations {input.blast_annotations} --coefficients {input.coefficients} --output {output} --translation_table {params.translation_table}
    """


rule merge_blastp_swissprot_results:
    """
    Merge the blast results with the annotated sequences
    """
    input:
        blast_annotations = Path(TEMP_DIR, "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{num_clusters}", "target{target_rank}", "{normalize}", "{predictionTask}", "trainProp" + "{train_proportion}",
                                 "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_swissprot.tsv"),
        coefficients = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv")
    params:
        script = config["scripts"]["merge_blast_results"],
        translation_table = lambda wildcards: dataset_table.loc[wildcards.dataset, "translation_table"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_swissprot_annotated.tsv")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --blast_annotations {input.blast_annotations} --coefficients {input.coefficients} --output {output} --translation_table {params.translation_table}
    """


rule merge_annotations_OHE:
    """
    Merge the annotations for each cluster onto the non-zero coefficients file for OHE features
    """
    input:
        annotations = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target{target_rank}_k{kmer_width}_s{kmer_step}_annotated.tsv"),
        coefficients = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv")
    params:
        script = config["scripts"]["merge_annotations"]
    output:
        coefs = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_annotated.tsv"),
        fasta = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences.fasta")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --annotations {input.annotations} --coefficients {input.coefficients} --output {output.coefs}
    """


rule run_blast_nonzero_features_OHE:
    """
    Run a blast search on the significant sequences to find the closest matches in the NCBI database for OHE features
    """
    input:
        fasta = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences.fasta")
    params:
        script = lambda wildcards: config["scripts"]["run_blast"],
        blast_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, "ohe", wildcards.num_clusters, "target"+ wildcards.target_rank, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "blast"),
        split_fasta_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, "ohe", wildcards.num_clusters, "target"+ wildcards.target_rank, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "split_fasta"),
        taxid = lambda wildcards: int(dataset_table.loc[wildcards.dataset, "taxid"]),
        entrez_email = config["entrez_email"] if config["entrez_email"] else 0,
        temp_dir = config["temp_dir"],
        blast_db_path = config["blast_db_path"] if config["blast_db_path"] else 0
    output:
        Path(TEMP_DIR, 
             "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{num_clusters}", "target" + "{target_rank}", "{predictionTask}", "trainProp" + "{train_proportion}",
             "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_blast.tsv")
    conda:
        config["envs"]["biopython_env_r"]
    shell:"""
        bash {params.script} {input.fasta} {params.split_fasta_temp_dir} {params.blast_temp_dir} {output} {threads} {params.taxid} {params.entrez_email} {params.temp_dir} {params.blast_db_path}
    """


rule run_blastp_nonzero_features_OHE:
    """
    Run a blast search on the significant sequences to find the closest matches in the NCBI database
    """
    input:
        fasta = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences.fasta")
    params:
        script = lambda wildcards: config["scripts"]["run_blastp"],
        blast_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, "ohe", wildcards.num_clusters, "target"+ wildcards.target_rank, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "blastp"),
        split_fasta_temp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, "ohe", wildcards.num_clusters, "target"+ wildcards.target_rank, wildcards.predictionTask, "trainProp" + wildcards.train_proportion, "split_fasta_blastp"),
        taxid = lambda wildcards: int(dataset_table.loc[wildcards.dataset, "taxid"]),
        translation_table = lambda wildcards: dataset_table.loc[wildcards.dataset, "translation_table"],
        blast_db_path = config["blast_db_path"] if config["blast_db_path"] else 0
    output:
        Path(TEMP_DIR, "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{num_clusters}", "target" + "{target_rank}", "{predictionTask}", "trainProp" + "{train_proportion}",
             "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_blastp.tsv")
    conda:
        config["envs"]["biopython_env_r"]
    shell:"""
        bash {params.script} {input.fasta} {params.split_fasta_temp_dir} {params.blast_temp_dir} {output} {threads} {params.taxid} {params.translation_table} {params.blast_db_path}
    """


rule merge_blast_results_OHE:
    """
    Merge the blast results with the annotated sequences for OHE features
    """
    input:
        blast_annotations = Path(TEMP_DIR, "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{num_clusters}", "target" + "{target_rank}", "{predictionTask}", "trainProp" + "{train_proportion}",
                                 "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_blast.tsv"),
        coefficients = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv")
    params:
        script = config["scripts"]["merge_blast_results"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blast_annotated.tsv")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --blast_annotations {input.blast_annotations} --coefficients {input.coefficients} --output {output}
    """


rule merge_blastp_results_OHE:
    """
    Merge the blast results with the annotated sequences for OHE features
    """
    input:
        blast_annotations = Path(TEMP_DIR, "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{num_clusters}", "target" + "{target_rank}", "{predictionTask}", "trainProp" + "{train_proportion}",
                                 "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_significant_sequences_blastp.tsv"),
        coefficients = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv")
    params:
        script = config["scripts"]["merge_blast_results"],
        translation_table = lambda wildcards: dataset_table.loc[wildcards.dataset, "translation_table"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blastp_annotated.tsv")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --blast_annotations {input.blast_annotations} --coefficients {input.coefficients} --output {output} --translation_table {params.translation_table}
    """


rule merge_annotations_genomes:
    """
    Merge the annotations for each cluster onto the non-zero coefficients file
    """
    input:
        annotations = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target1_k{kmer_width}_s{kmer_step}_annotated.tsv"),
        coefficients = Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "genomes", "{normalize}", "{dataset}_{model}_adelie_genomes_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients.tsv")
    params:
        script = config["scripts"]["merge_annotations"]
    output:
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "genomes", "{normalize}", "{dataset}_{model}_{predictionTask}_genomes_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_annotated.tsv")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --annotations {input.annotations} --coefficients {input.coefficients} --output {output}
    """


## PLOTTING RULES --------------------------------
# These rules are for plotting the results of the blast searches and the annotated sequences

rule plot_blast_features:
    """
    Merge the blast results with the annotated sequences
    """
    input:
        nonzero_features = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blastp_annotated.tsv"),
        nonzero_features2 = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blast_annotated.tsv"),
        sequences = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_sample_sequences.tsv"),
        clusters = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv"),
        metadata = lambda wildcards: dataset_table.loc[wildcards.dataset, "metadata_file"],
        feather = Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_{normalize}.feather")
    params:
        script = config["scripts"]["plot_blast_results"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blast_annotated_plots.pdf")
    conda:
        config["envs"]["default_r"]
    resources:
        msa=1 # one downside of the msa package is that it creates a temp file that cannot be redirected, so we limit to one at a time
    shell:"""
        Rscript --vanilla {params.script} \
        --nonzero_annotations {input.nonzero_features}\
        --clusters {input.clusters} \
        --feather_file {input.feather} \
        --sample_seqs {input.sequences} \
        --metadata {input.metadata} \
        --output {output}
    """


rule plot_blast_heatmaps:
    """
    Merge the blast results with the annotated sequences
    """
    input:
        nonzero_features = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blastp_annotated.tsv"), 
        nonzero_features2 = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blast_annotated.tsv"), 
        sequences = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_sample_sequences.fasta"),
        clusters = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv"),
        metadata = lambda wildcards: dataset_table.loc[wildcards.dataset, "metadata_file"],
        feather = Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_{normalize}.feather")
    params:
        script = config["scripts"]["plot_blast_heatmaps"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_heatmaps.pdf")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} \
        --nonzero_annotations {input.nonzero_features}\
        --clusters {input.clusters} \
        --feather_file {input.feather} \
        --sample_seqs {input.sequences} \
        --metadata {input.metadata} \
        --output {output}
    """


rule plot_blast_features_OHE:
    """
    Merge the blast results with the annotated sequences
    """
    input:
        nonzero_features = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blastp_annotated.tsv"), 
        nonzero_features2 = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blast_annotated.tsv"),
        sequences = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_sample_sequences.tsv"),
        clusters = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv"),
        metadata = lambda wildcards: dataset_table.loc[wildcards.dataset, "metadata_file"],
        feather = Path(TEMP_DIR, "{dataset}", "{dataset}_ohe_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.feather")
    params:
        script = config["scripts"]["plot_blast_results"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blast_annotated_plots.pdf")
    conda:
        config["envs"]["default_r"]
    resources:
        msa=1 # one downside of the msa package is that it creates a temp file that cannot be redirected, so we limit to one at a time
    shell:"""
        Rscript --vanilla {params.script} \
        --nonzero_annotations {input.nonzero_features}\
        --clusters {input.clusters} \
        --feather_file {input.feather} \
        --sample_seqs {input.sequences} \
        --metadata {input.metadata} \
        --output {output}
    """


rule plot_blast_heatmaps_OHE:
    """
    Merge the blast results with the annotated sequences
    """
    input:
        nonzero_features = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blastp_annotated.tsv"), 
        nonzero_features2 = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_blast_annotated.tsv"), 
        sequences = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_sample_sequences.fasta"),
        clusters = Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_target{target_rank}_k{kmer_width}_s{kmer_step}.tsv"),
        metadata = lambda wildcards: dataset_table.loc[wildcards.dataset, "metadata_file"],
        feather = Path(TEMP_DIR, "{dataset}", "{dataset}_ohe_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.feather")
    params:
        script = config["scripts"]["plot_blast_heatmaps"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_{predictionTask}_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_trainProp{train_proportion}_nonzero_coefficients_heatmaps.pdf")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} \
        --nonzero_annotations {input.nonzero_features}\
        --clusters {input.clusters} \
        --feather_file {input.feather} \
        --sample_seqs {input.sequences} \
        --metadata {input.metadata} \
        --output {output}
    """


rule plot_embeddings_umap:
    """
    Plot UMAP of the top variance embedding per cluster for the given dataset.
    This rule is used to visualize the embeddings of the clusters.
    """
    input:
        embeddings = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_umap_{select_type}_{cluster_type}_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}_{normalize}.feather"),
        metadata = lambda wildcards: dataset_table.loc[wildcards.dataset, "metadata_file"]
    params:
        script = config["scripts"]["plot_embeddings_umap"],
        num_PCs = config["extended_options"]["num_PCs_umap"]
    output:
        Path('results', "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_umap_results_top{num_clusters}_target{target_rank}_k{kmer_width}_s{kmer_step}.pdf")
    conda:
        config["envs"]["default_r"]
    shell:"""
        Rscript --vanilla {params.script} --embeddings {input.embeddings} --metadata {input.metadata} \
        --output {output} --num_PCs {params.num_PCs}
    """
