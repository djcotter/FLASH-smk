# *FLASH*: FunctionaL Assigning Sequence Homing

## Table of Contents

* [Description](#description)

### Running the workflow

* [1. Install mamba, snakemake, and SPLASH](#1-install-mamba-snakemake-and-splash)
* [2. Ensure dataset_table.csv is filled out correctly](#2-ensure-dataset_tablecsv-is-filled-out-correctly)
* [3. Configure parameters in config.yaml](#3-ensure-that-the-parameters-you-want-to-use-are-specified-in-the-configyaml-file-snakemake-will-use-these-to-fill-out-the-wildcards-in-the-snakefile)
* [4. Update all paths in config.yaml](#4-update-all-paths-in-the-config-file-configyaml)
* [5. Run snakemake](#5-run-snakemake)

### Execution modes

* [One Hot Encoding](#one-hot-encoding)
* [Embedding mode or genome predictions](#embedding-mode-or-genome-predictions)

### Inputs and outputs

* [Input Files](#input-files)
* [Target rules and output files](#target-rules-and-output-files)

### Guides and examples

* [Guide: Running SPLASH and preparing inputs for FLASH](#guide-running-splash-and-preparing-inputs-for-flash)
* [Required Outputs](#required-outputs)
* [Running FLASH](#running-flash)
* [Example Run of FLASH on H5N1 Sample Data](#example-run-of-flash-on-h5n1-sample-data-step-by-step)

### Data considerations

* [Types of data that can be analyzed](#types-of-data-that-can-be-analyzed)
* [Very small datasets](#very-small-datasets)
* [Data where SPLASH has few significant anchors](#data-where-splash-has-very-few-significant-anchors-or-anchors-that-pass-the-threshold)

### Resource management

* [Resource Management, Threads, and Profiles](#resource-management-threads-and-profiles)
* [Global core limits (-j / --cores)](#global-core-limits--j----cores)
* [Using the provided Slurm profile](#using-the-provided-slurm-profile)
* [Adapting for local execution](#adapting-for-local-execution)
* [Example for a local machine](#example-for-a-local-machine)


## Description

The following pipeline uses `Snakemake` to run metadata-driven prediction and attribution on SPLASH results. Once SPLASH has been run on data and Snakemake has been installed, it should take *5-10 minutes* to configure the metadata table and check on input files, to prepare FLASH to run on a new dataset.

### IMPORTANT NOTE

This pipeline has been modified to run locally but it requires extensive resources. It can be submitted to a scheduler by providing a `--profile` that contains a job submission command. See `slurm_profile/` for an example.

## Running the workflow

### 1. Install `mamba`, `snakemake`, and `SPLASH`

You must have conda installed with the `mamba` package manager. See [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge) for info on how to install mamba. After that, create a snakemake environment using the following command:

```{bash}
mamba create -n snakemake python=3.12 snakemake -c conda-forge -c bioconda
```

To use the executor profile for cluster job submission (with sbatch) you also must run the following while the new conda environment is active (to activate: `conda activate snakemake`):

```{bash}
mamba install snakemake-executor-plugin-cluster-generic -c bioconda
```

#### Installing SPLASH

To install `SPLASH`, download the relase here: [github.com/refresh-bio/SPLASH/releases/tag/v2.11.6](https://github.com/refresh-bio/SPLASH/releases/tag/v2.11.6). Unzip it into a folder in the project root and change the `splash_bin` parameter in `config.yaml` to match. You can also run the script `get_splash.sh` if you are using x64 linux.

**The FLASH pipeline uses several SPLASH utilities throughout but the pipeline requires you to have run SPLASH on your data *before* running the full set of snakemake commands.**

### 2. Ensure `dataset_table.csv` is filled out correctly

This file contains all of the datasets on which we want to run the FLASH Snakemake pipeline. Replace the current datasets with ones you would like to run. To run on a new dataset, add a new row with a **short name**, a path to the **SPLASH run folder** (where you ***must have already run*** SPLASH with the option to dump SATC files), a path to the **metadata file**, a path to a **lookup table for artifact filtering** (you don't need to change this from the other rows), and a **translation table** (this is an integer corresponding to the correct genetic code to use for translation).

#### Inside the SPLASH run folder, the FLASH pipeline is looking for the following files specifically:

1. `sample_sheet.txt` generated in order to originally run SPLASH.
2. `result.after_correction.scores.tsv` 
3. `sample_name_to_id.mapping.txt` 
4. `result_satc/` containing all of the satc files for the run.

There is an example script for running SPLASH in the `resources/utility_scripts` folder (run_splash.sh). See below for further detail.

#### The metadata file needs to be formatted as follows

1. Column 1 contains the same sample ids as were used to run `SPLASH` and is named `sample_name` (the script will attempt to assign the first column as `sample_name` otherwise)
2. The other columns contain a simple column name describing the metadata as well as the observations of the metadata.

### 3. Ensure that the parameters you want to use are specified in the `config.yaml` file. Snakemake will use these to fill out the wildcards in the `Snakefile`

The current wildcards and parameters are defined as follows:

```{yaml}
options:
  seqs_from_raw_data: true # set to true if you want FLASH to assemble it's sequences from raw FASTQ files, false if you want to default from the pre-processed SPLASH results
  feature_processing_method: "pca" # set to either 'top_variance' or 'pca' for how to process embeddings prior to glmnet. default: 'top_variance'
  generate_plots: false # set to true to generate plots after model training
  num_clusters: 20000 # number of clusters to use for clustering anchors
  filters: ["filter1"] # define which anchor selection filters to use
  cluster_types: ["shiftDist-levFilter", "masked-nucleotide-clustered"] # define which clustering methods to use
  models: ["hyena"] # currently only 'hyena' is supported
  normalize_embeddings: ['normalized', 'unnormalized'] # define whether to use normalized or unnormalized embeddings, this will do both
  train_proportion: 0.8 # proportion of samples to use for training GLMnet models
  anchor_length: 27 # this is the length of the anchor in nucleotides, change as needed for different SPLASH runs
  target_length: 27 # this is the length of the target in nucleotides, change as needed for different SPLASH runs
  target_rank: 1 # this is the rank of the target to use when assembling the anchor-target pairs, change as needed

# you should not need to change anything below this line but you can alter these if needed
extended_options:
  num_anchors_to_select: 3000000 # number of anchors to select from each dataset for clustering
  effect_size_cutoff: 0.6 # this is the cutoff for selecting an anchor from the SPLASH results
  distance_threshold: 5 # distance threshold for filtering anchors after clustering
  num_embedding_features_to_keep:
    glmnet: 100 # number of top variance embedding features to keep for GLMnet
    umap: 1 # number of top variance embedding features to keep for UMAP
  num_PCs_umap: 10 # number of principal components to use for UMAP plotting
  min_samples_adelie: 28 # minimum number of samples in the smallest class for running adelie (requires at least N samples per class)

```

You can change these by specifying new values for the parameters.
In order to use your own filter clustering type, you must also provide new scripts (in the `scripts:` section) containing modified code for anchor selection and clustering. See the existing scripts for examples of how to do this.

#### Generating annotations and plots

To generate annotations and plots you can set the flag `generate_plots:` in the config file. This will run BLAST on the output sequences for each analysis and generate a series of visualizations. This step can be very slow as the default behavior is to run blast in `remote` mode. By pre-downloading the blast databases for `core_nt` and `refseq_protein` you can speed this step up.

To point at the folder containing the local blast databases, modify the config file, `config.yaml`, and change the field `blast_db_path` to point to the folder containing the local blast databases.

These can be downloaded with `blast+`; they are installed using `update_blastdb.pl --decompress core_nt refseq_protein`. For further instructions on downloading local copies of these databases, see <https://www.ncbi.nlm.nih.gov/books/NBK569850/>.

To enable taxonomic filtering of the blast features, using the `taxid` field in the file `dataset_table.csv`, you must A) be using local blast databases as remote blast does not support taxonomic filtering and B) You also must additionally download and decompress the file `taxdb.tar.gz` from <https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz> in the local blast database folder. Follow the instructions in the *Taxonomic filtering for BLAST databases* section of <https://www.ncbi.nlm.nih.gov/books/NBK569839/> for further details.

### 4. Update all paths in the config file, `config.yaml`

The `config.yaml` file is also used by the Snakemake workflow to specify the correct locations of all paths and some parameters.

It is necessary to change the fields below to match your setup and computing environment:

1. Change the `entrez_email` entry if you intend to run the pipeline with blast plots enabled.
2. Change the `temp_dir` entry to point to a location that has a large amount of storage and fast I/O speeds. I set this to the temporary SCRATCH storage on our file systems. Depending on the size of the input data, these intermediate files can be quite large.
3. Change `splash_bin` to point to the installed location of SPLASH. See note above on how to install SPLASH.

### 5. Run `snakemake`

#### One Hot Encoding

You can run the workflow for a given target rule (for example using One Hot Encoding) by typing the following in an interactive node:

```{bash}
export NUM_CORES=1
snakemake --sdm conda -j $NUM_CORES all_ohe
```

This will execute `snakemake` to run jobs using up to `$NUM_CORES` maximum in parallel and will produce the files specified by the rule `all_ohe`.

The flag `--sdm conda` is important because it will tell snakemake to build and use the required conda environments specified in the `envs/` directory. *Note: This can be quite slow whenever the environment is installed for the first time or is changed.*

#### Embedding mode or genome predictions

**IMPORTANT:** **Before you run in embeddings mode you will need to download the singularity image containing the model state and necessary code into the folder `containers/`. This can be done by ensuring `singularity` is installed and running `bash containers/setup.sh` to pull the container image from GitHub Container Registry.**

If you want to run FLASH using the embedding mode or genome prediction mode, you will need a GPU in the environment where you run `snakemake` or you will need to use a profile that can allocate GPUs and other resources (like the one provided for the slurm scheduler in the directory `slurm_profile/`).

**The keyword `genomes` is arbitrary and this mode can be used to perform predictions across any type of external data not used to run SPLASH or train the FLASH model. We have used this mode to run predictions on other short-read datasets, on genomes, and on long-read datasets.**


The command for running the code in an interactive environment (where a GPU is present if running embeddings) is similar to above:

```{bash}
export NUM_CORES=1
snakemake --sdm conda -j $NUM_CORES all_embeddings # for embedding mode

snakemake --sdm conda -j $NUM_CORES all_genomes # for genomes mode, requires additional files/setup
```

The code can also be run using an automatic scheduler. The included example submission script (`resources/utility_scripts/run_snakemake.sbatch`) and profile (`slurm_profile/config.v8+.yaml`) can submit the pipeline and request the required resources including GPU resources. For running on a cluster using `slurm` you can use this config using the --profile slurm_profile/config.v8+.yaml. On some clusters, Snakemake needs to know where the shared Conda/Miniforge installation lives. In that case, explicitly pass the Conda base path.

```{bash}
snakemake --sdm conda --use-conda --conda-base-path /path/to/miniforge3 --profile slurm_profile/ all_embeddings

snakemake --sdm conda --use-conda --conda-base-path /path/to/miniforge3 --profile slurm_profile/ all_embeddings
```

**NOTE: You must modify the profile to match your own partitions resources, and constraints necessary for your cluster. The provided config is an example for our local HPC resources. This can be adapted to different schedulers by modifying the profile but is currently only set to work for `slurm`.**

## Input Files

Input files and paths are detailed in `dataset_table.csv` and the columns are descriptive. Metadata paths as well as SPLASH results paths are stored there.

- `dataset_short_name`: A short name for the dataset (Don't include underscores).

- `SPLASH_results`: Path to the SPLASH run folder containing results. Should contain `result.after_correction.scores.tsv`, `sample_name_to_id.mapping.txt`, and the folder `result_satc`. (An example script for running SPLASH to produce these outputs is in the `resources/helper_scripts` folder).
  - If using the flag `seqs_from_raw_data: true` in the config file, this step will additionaly allow the script to find the `sample_sheet.txt` file. This file is the default one used to run SPLASH and contains two columns (`sample_name` and `file_path`) that map sample names to raw FASTQ files. This will allow the script to extract sequences directly from the raw data instead of from the SPLASH outputs which circumvents the limitations of SPLASH in terms of outputting sequences that appear less frequently.

- `metadata_file`: Path to the metadata file associated with the dataset.
  - The metadata file must contain a column named `sample_name` that matches the sample names used in the SPLASH run.

- `lookup_table`: Path to the lookup table for annotation.

- `translation_table`: Integer corresponding to the correct genetic code for translation.

- `taxid`: Taxonomic ID associated with the dataset. Used to restrict the BLAST-based annotation pipeline to a specific taxon.

### Additional inputs for genome predictions

For running genome-based predictions, additional fields must be filled out in the `dataset_table.csv` file. If you are not running genome-based predictions, you can ignore set these fields to NA.

- `genomes_folder`: Path to the folder containing the genome fasta files.

- `genome_list`: List of "testing data" genomes for extracting. Names must map to the filenames in `genomes folder` (e.g. name:`SAMPLE1` -> file: `SAMPLE1.fasta`)

- `genome_metadata_file`: Path to the genome metadata file. Takes the same form as the metadata file above. Column names and labels must match the original metadata file.


## Target rules and output files

Output files are available in the `results/` folder with the dataset name as a subfolder and subsequent subfolders defining different branching points along the way. The flag `GENERATE_PLOTS` can be changed to `True` to force the script to perform the annotation and plotting steps.

- The target rule **`rule_all_embeddings`** genarates all embedding-based predictions and summary files using the nucleotide language model.
- The target rule **`rule_all_genomes`** generates models based on embeddings data and then predicts on these data using test data supplied in the genome directory as specified in the `dataset_table.csv` file.
- The target rule **`rule_all_ohe`** formats the data using One Hot Encoding and predicts on the data this way. This does not require setting up the gpu-based embedding pipeline.
- The target rule **`rule_all_umap`** generates the unsupervised umap clusters using the top variance embedding per cluster in the formatted data.

## Guide: Running SPLASH and preparing inputs for FLASH

### SPLASH

You can use the following to prepare the sample data for the example. You can use these scripts to run SPLASH on your own data to prepare it for the FLASH pipeline as well. 

To run the pipeline, you must first run SPLASH on your raw sequencing data to generate the necessary input files. An example script for running SPLASH is provided in the `resources/utility_scripts` folder (run_splash.sh). This script requires a file called `sample_sheet.txt` that contains two columns: `sample_name` and `file_path`, where `sample_name` is a unique identifier for each sample and `file_path` is the path to the raw sequencing data file (in FASTQ format). Typically, when using paired end data, we run SPLASH only on the forward reads. The options in the example script are set to generate the necessary SATC files for FLASH as they are not generated by default. You can modify the anchor and target lengths as needed, but these must match the values used in the FLASH pipeline. By default, we typically use 27 for both anchor and target lengths, which corresponds to ~9 amino acids each when translated.

### Required Outputs

After running SPLASH, ensure that the output folder contains the following files:

1. `result.after_correction.scores.tsv`
2. `sample_name_to_id.mapping.txt`
3. `result_satc/` containing all of the satc files for the run.
4. `sample_sheet.txt` (necessary when using the `seqs_from_raw_data: true` option in the config file)

These files are the necessary inputs for the FLASH pipeline. You can then proceed to fill out the `dataset_table.csv` file with the path to the SPLASH output folder containing these files and to the metadata file associated with your samples. The metadata file should have a column named `sample_name` that matches the sample names used in the SPLASH run.

### Running FLASH

You should not need to modify any of the parameters in the `config.yaml` unless you change the anchor or target lengths when running SPLASH as mentioned above. If you do change these lengths, ensure that the `anchor_length` and `target_length` variables in the `config.yaml` file are updated accordingly. After confirming that all paths and parameters are correctly set, you can run the FLASH pipeline using Snakemake as described in the previous sections.

Note that if you use a very short anchor length, you should also modify the `CLUSTER_TYPES` variable in the `Snakefile` to avoid using clustering methods that rely on longer anchors, such as `shiftDist-levFilter`. This should instead be set to use `noCluster`:

```{python}
# If using default anchor length of 27, use:
ANCHOR_LENGTH = 27
TARGET_LENGTH = 27
CLUSTER_TYPES = ["shiftDist-levFilter"]

# If using a shorter anchor length, e.g., 9, use:
ANCHOR_LENGTH = 9
TARGET_LENGTH = 27
CLUSTER_TYPES = ["noCluster"]
```

## Example Run of FLASH on H5N1 Sample Data (Step-by-Step)

This section walks through a complete, reproducible example using the provided H5N1 dataset. The goal is to go from raw data → SPLASH → FLASH predictions.

### Step 0: Create and activate environment

```bash
mamba create -n flash_env python=3.12 snakemake -c conda-forge -c bioconda
conda activate flash_env
mamba install snakemake-executor-plugin-cluster-generic
```


### Step 1: Download example data

```bash
bash generate_example_data.sh
```

This will:

* Download H5N1 sequencing data into `example_data/`
* Generate `sample_sheet.txt` for SPLASH

> Requires SRA Toolkit (included in `envs/data_env.yml`)


### Step 2: Run SPLASH on the example data

```bash
bash run_splash_example.sh
```

After completion, confirm the SPLASH output folder contains:

* `result.after_correction.scores.tsv`
* `sample_name_to_id.mapping.txt`
* `result_satc/`
* `sample_sheet.txt`


### Step 3: Set up dataset table

Edit `dataset_table.csv` and add/update a row like:

```
dataset_short_name,SPLASH_results,metadata_file,lookup_table,translation_table,taxid
H5N1_example,example_data/splash_output,resources/metadata/H5N1_example_metadata.csv,resources/lookup_tables/default_lookup.tsv,1,11320
```


### Step 4: Configure `config.yaml`

At minimum, verify:

```yaml
splash_bin: /path/to/SPLASH/bin
temp_dir: /path/to/fast/storage
```

Optional but recommended:

```yaml
generate_plots: false
```


### Step 5: Download model container (REQUIRED for embeddings)

```bash
bash containers/setup.sh
```

This pulls the required Singularity container into `containers/`.


### Step 6: Create a local execution profile

Create:

```
profiles/local/config.yaml
```

With the following contents:

```yaml
cores: 8
jobs: 8
use-conda: true

rerun-incomplete: true
printshellcmds: true
latency-wait: 60

default-resources:
  - mem_mb=4000
  - disk_mb=10000

set-threads:
  all_ohe: 4
  all_embeddings: 4

set-resources:
  all_embeddings:
    mem_mb: 16000
```


### Step 7: Run FLASH

#### Option A: CPU-only (One Hot Encoding)

```bash
snakemake --profile profiles/local --sdm conda all_ohe
```

#### Option B: Embedding mode (requires GPU + container)

```bash
snakemake --profile profiles/local --sdm conda all_embeddings
```


### Step 8: Inspect outputs

Results will appear in:

```
results/H5N1_example/
```

Key outputs include:

* Model predictions
* Feature importance summaries
* Clustered anchor outputs

### Notes & Common Pitfalls

* **Container step is mandatory** for embedding mode — missing this will cause runtime failures.
* **Paths must be absolute or correct relative paths** in both `dataset_table.csv` and `config.yaml`.
* If jobs fail due to memory:

  * Reduce `cores` / `jobs`
  * Increase `mem_mb` in profile
* If SPLASH anchors are sparse:

  * Lower `effect_size_cutoff` in `config.yaml`

--- 

## Types of data that can be analyzed
In principal FLASH has been written to run on any data with sturctured phenotype/metadata labels provided. However, there are some caveats to this. 

#### Very small datasets
The FLASH paramaters are built to only include metadata categories from a metadata file with >= 28 samples in the smallest class of a metadata category. This means that only datasets with exactly 56 samples that are evenly split into 28 and 28 will pass the adelie step of the pipeline. 
This can be changed by adjusting the paramater `min_samples_adelie` in the `extended_options` section of the `config.yml` file. The default number was chosen after extensive testing in many datasets but feel free to experiment with smaller datasets. It is possilbe the the unsupervised clustering mode will still work for smaller amounts of data. 

#### Data where SPLASH has very few significant anchors or anchors that pass the threshold 
If SPLASH does not provide a good insight into the structure encoded in your set of samples, it will not provide many anchors with hgih effect size. Because we use a default of 0.6 to filter for high effect size anchors, this could mean too few anchors enter into the analysis in the first place. This can be dropped in the `config.yaml` by adjusting the `effect_size_cutoff` praramter under `extended_options`. 

## Resource Management, Threads, and Profiles

Snakemake determines how many jobs run in parallel and how many threads each rule can use through a combination of **rule definitions**, **command-line arguments**, and **execution profiles**.

### How threads are defined in the workflow

Each rule in the `Snakefile` may specify a `threads:` directive. This value represents the *maximum* number of CPU threads that a single job of that rule is allowed to use. For example:

```python
rule example:
    threads: 8
```

This does **not** guarantee that 8 threads will always be used. Instead, it defines an upper bound that Snakemake and the execution backend will respect.

### Global core limits (`-j` / `--cores`)

When running locally, the `-j` (or `--cores`) argument controls the total number of CPU cores available across all jobs:

```bash
snakemake -j 16
```

Snakemake will schedule jobs so that the sum of allocated threads does not exceed this limit. For example, two jobs requiring 8 threads each could run in parallel under `-j 16`.

### Profiles override and extend resource behavior

When using a profile (e.g., the provided Slurm profile in `slurm_profile/`), resource handling is largely delegated to the scheduler. Profiles can:

* Override `threads`
* Define additional resources (e.g., `mem_mb`, `time`, `gpus`)
* Map Snakemake resources to scheduler submission parameters (e.g., `sbatch` flags)
* Control job parallelism independently of the local `-j` flag

For example, a Slurm profile may include entries such as:

```yaml
set-threads:
  some_rule: 16

resources:
  some_rule:
    mem_mb: 64000
    runtime: 120
```

In this case, the profile—not the `Snakefile`—determines the actual resources requested at runtime.

### Using the provided Slurm profile

The example profile in `slurm_profile/` demonstrates how to:

* Submit jobs via `sbatch`
* Request CPUs, memory, and GPUs per rule
* Scale execution across a cluster

**Important:** This profile is cluster-specific. You must modify it to match:

* Available partitions/queues
* Memory and runtime limits
* GPU availability and constraints

### Adapting for local execution

If you are running FLASH on a local machine (no scheduler), you should:

1. Use the `-j` flag to reflect your available CPU cores.
2. Review the `threads:` directives in the `Snakefile`.
3. Optionally reduce thread counts for memory-intensive rules.
4. Avoid overcommitting resources (e.g., running too many high-thread jobs in parallel).

If needed, you can also define a simple local profile by copying relevant sections (such as `set-threads` and `resources`) from the Slurm profile and removing scheduler-specific fields.

### Key takeaway

* **`threads:` in rules** → defines per-job maximums
* **`-j / --cores`** → defines total parallelism (local execution)
* **Profiles** → dynamically override and map resources to your execution environment

For most users:

* On a cluster → use and adapt a profile (e.g., Slurm)
* On a workstation → tune `-j` and adjust rule threads as needed


### Example for a local machine
Here’s a minimal local profile you can drop into something like profiles/local/. This avoids any cluster or executor plugin logic and just helps you centrally manage threads and resource scaling for a single machine.

`profiles/local/config.yaml`:

```
# Maximum number of cores available to Snakemake
cores: 8

# Keep Snakemake from overloading your system
jobs: 8

# Use conda environments defined in envs/
use-conda: true

# Optional: rerun incomplete jobs
rerun-incomplete: true

# Print shell commands (useful for debugging)
printshellcmds: true

# Latency wait for filesystem (useful on network filesystems)
latency-wait: 60

# Default resources applied to rules that don't specify them
default-resources:
  - mem_mb=4000
  - disk_mb=10000

# Override threads for specific heavy rules (example)
set-threads:
  all_embeddings: 4
  all_genomes: 4

# You can also cap resources per rule if needed
set-resources:
  all_embeddings:
    mem_mb: 16000
  all_genomes:
    mem_mb: 16000
How to use it
snakemake --profile profiles/local all_ohe
```

This replaces the need to manually pass `-j`, `--cores`, and other flags every time.

#### How to adapt it to your machine

`cores / jobs`: Set this to the number of CPU cores you actually want Snakemake to use (not necessarily the maximum your system has—leave some headroom).

`set-threads`: Use this to downscale rules that are too aggressive by default. For example, if a rule asks for 16 threads but your machine only has 8 cores, cap it here.

`default-resources`: Prevents too many memory-heavy jobs from running simultaneously. Increase mem_mb if you see jobs getting killed.

#### Philosophy vs. the Slurm profile

The Slurm profile translates Snakemake jobs into scheduler submissions. This local profile simply constrains and balances execution within a single machine. It contains no submission commands, no plugins and no external dependencies.
