# process_fastqs_for_embedding_and_ohe.R
# Daniel Cotter
# 2024-09-13

# This script takes in a path to the original SPLASH sample_sheet.txt file
# in order to map back to the original data as well as as well as a list of anchors
# and a list of anchor clusters. It then dumps the anchors from the original fasta files
# and proceeds to transform these new files into a set of downstream files that will
# be formatted for input into the embedding and OHE steps.

## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(furrr))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-a", "--anchor_file"), "File with list of anchors",
    type =
      "character"
  ),
  make_option(
    c("-c", "--cluster_file"),
    "File with list of anchor clusters",
    type = "character"
  ),
  make_option(c("i", "--id_mapping"), "File with sample id mapping",
    type =
      "character"
  ),
  make_option(
    c("--sample_sheet"),
    "Path to original SPLASH sample_sheet.txt file",
    type = "character"
  ),
  make_option(c("-o", "--output_prefix"), "Output prefix.", type = "character"),
  make_option(
    c("--temp_dir"),
    "Temporary directory to store intermediate files",
    type = "character"
  ),
  make_option(
    c("--satc_util_bin"),
    "Path to SATC Util binary folder",
    type = "character",
    default = "/oak/stanford/groups/horence/dcotter1/satc_utils/"
  ),
  make_option(
    c("--num_cores"),
    "Number of cores to use",
    type = "integer",
    default = 8
  ),
  make_option(
    c("--single_cell"),
    help = "Whether to merge the first two columns of the SATC file to create sample barcode pairs",
    type = "logical",
    default = FALSE,
    action = "store_true"
  ),
  make_option(
    c("--target_rank"),
    "Rank of the target to use when assembling anchor-target pairs",
    type = "integer",
    default = 1
  ),
  make_option(
    c("--apply_cluster_filter"),
    paste(
      "If supplied, apply either a filter eliminating clusters with fewer",
      "than a certain fraction of anchors",
      "or a filter eliminating clusters with binary targets",
      "(i.e clusters with only 2 unique sequences across all samples).",
      "Format should be 'fractionMissing:N' or 'binaryTarget' where N",
      "is the fraction (i.e. 0.05) of anchors that are allowed to be missing",
      "for it to be retained."
    ),
    type = "character",
    default = NULL
  )
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set up parallel processing
plan(multicore, workers = opt$num_cores)
setDTthreads(max(1, opt$num_cores - 4))

# check that user specified all files
if (!file.exists(opt$anchor_file) |
  !file.exists(opt$cluster_file) |
  is.null(opt$sample_sheet) | is.null(opt$output_prefix)) {
  stop("Must provide anchor file, cluster file, satc files, id mapping file, and output prefix")
}

# check that the cluster filter argument is in the correct format if supplied
if (!is.null(opt$apply_cluster_filter)) {
  if (!grepl("fractionMissing:\\d+", opt$apply_cluster_filter) &&
    opt$apply_cluster_filter != "binaryTarget") {
    stop("Invalid format for --apply_cluster_filter. Must be 'fractionMissing:N' or 'binaryTarget'.")
  }
}

# create a temporary directory to store intermediate files
if (!is.null(opt$temp_dir)) {
  temp_dir <- ifelse(grepl("/$", opt$temp_dir),
    opt$temp_dir,
    paste0(opt$temp_dir, "/")
  )
  system(paste("mkdir -p", temp_dir))
} else {
  temp_dir <- file.path(dirname(opt$output_prefix), "tmp/")
  system(paste("mkdir -p", temp_dir))
}

system(paste0("rm ", temp_dir, "/*/*"))
system(paste0("rm ", temp_dir, "/*"))

# read in the anchor cluster file
anchor_clusters <- fread(opt$cluster_file,
  header = F,
  col.names = c("cluster_id", "anchor")
) %>%
  group_by(cluster_id) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# create the data frame that is used to process the satc files with fafq filter
satc_files <- fread(
  opt$sample_sheet,
  header = F,
  col.names = c("sample_id", "fastq_file")
)
satc_files <- satc_files %>% mutate(satc_file = file.path(temp_dir, "satc_files", paste0(sample_id, ".satc")))
if (!dir.exists(file.path(temp_dir, "satc_files"))) {
  system(paste("mkdir -p", file.path(temp_dir, "satc_files")))
}
# add a new column for the dumped satc files
satc_files <- satc_files %>%
  mutate(satc_dump = gsub(
    ".satc",
    ".satc.dump",
    file.path(temp_dir, "dumped", basename(satc_file))
  ))
system(paste("mkdir -p", file.path(temp_dir, "dumped")))

# declare a satc file for the output of all the dump files
all_satc_file <- file.path(temp_dir, "all_satc_merged.txt")

# get anchor and target lengths from the name of the output file
anchor_len <- nchar(readLines(opt$anchor_file, n = 1))
kmer_len <- str_extract(basename(opt$output_prefix), "(?<=_k)[0-9]+(?=_s\\d+)")
kmer_len <- as.integer(kmer_len)
# calculate target length
target_len <- kmer_len - anchor_len


# create the individual satc files for each sample using fafq filter and
# the provided anchor list
cat("Creating SATC files using fafq_filter...\n")
future_walk2(
  satc_files$fastq_file,
  satc_files$satc_file,
  \(x, y) system(
    paste0(
      file.path(opt$satc_util_bin, "fafq_filter"),
      " -d ",
      opt$anchor_file,
      " -i ",
      x,
      " -o ",
      y,
      " -n ",
      opt$target_rank + 1, # get target_rank +1 targets per anchor
      " --anchor_len ",
      anchor_len,
      " --target_len ",
      target_len
    )
  )
)

# dump the satc files
cat("Dumping SATC files...\n")
future_walk2(
  satc_files$satc_file,
  satc_files$satc_dump,
  \(x, y) system(
    paste0(
      file.path(opt$satc_util_bin, "satc_dump"),
      " --anchor_list ",
      opt$anchor_file,
      " ",
      x,
      " ",
      y
    )
  )
)

# merge all of the dumped satc files into one file
system(paste("rm", all_satc_file))
# when merging them into one file keep columns 2-3 and use sample_id as the sample name
cat("Merging dumped SATC files...\n")
walk2(
  satc_files$satc_dump,
  satc_files$sample_id,
  \(x, y) system(
    paste0(
      "awk 'BEGIN {OFS=\"\\t\"} {print \"",
      y,
      "\", $2, $3, $4}' ",
      x,
      " >> ",
      all_satc_file
    )
  )
)

# remove any lines that start or end in [ACTG]
cat("Removing invalid anchor-target pairs...\n")
system(paste(
  "grep -v '^[ACTG]{3}' ",
  all_satc_file,
  " | grep -v '[ACTG]$' > ",
  file.path(temp_dir, "all_satc_merged_no_anchor.txt")
))
system(paste(
  "mv",
  file.path(temp_dir, "all_satc_merged_no_anchor.txt"),
  all_satc_file
))


if (opt$single_cell) {
  # exit since this is not supported for single cell data
  stop(
    "Single Cell data is not supported for this script. Snakemake will default to the original process_splash_data_for_embedding.R script."
  )
}

# undump the satc files
all_satc_undumped <- file.path(temp_dir, "all_satc_merged.undumped.txt")
all_satc_temp_mapping <- file.path(temp_dir, "all_satc_merged.temp_mapping.txt")
cat("Undumping merged SATC file...\n")
system(
  paste(
    file.path(opt$satc_util_bin, "satc_undump"),
    "-i",
    all_satc_file,
    "-o",
    all_satc_undumped,
    "-m",
    all_satc_temp_mapping
  )
)


# filter the satc files
all_satc_filtered <- file.path(temp_dir, "all_satc.filtered.txt")
cat("Filtering merged SATC file...\n")
cat(paste("Using target rank:", opt$target_rank, "\n"))
system(
  paste(
    file.path(opt$satc_util_bin, "satc_filter"),
    "-i",
    all_satc_undumped,
    "-o",
    all_satc_filtered,
    "-d",
    opt$anchor_file,
    "-n",
    opt$target_rank + 1 # get target_rank +1 targets per anchor
  )
)


# redump the satc file
all_satc_filtered_dump <- file.path(temp_dir, "all_satc.filtered.dump")
cat("Redumping filtered SATC file...\n")
system(
  paste(
    file.path(opt$satc_util_bin, "satc_dump"),
    "--sample_names",
    all_satc_temp_mapping,
    all_satc_filtered,
    all_satc_filtered_dump
  )
)

cat("Reading in dumped SATC_file...\n")
# read in the dumped satc file
satc_dt <- fread(
  all_satc_filtered_dump,
  header = FALSE,
  col.names = c("sample", "anchor", "target", "count")
)

cat("Ordering SATC file by sample, anchor, desc(count)...\n")
# now filter for ONLY the specified target rank
setorder(satc_dt, sample, anchor, -count)

# add a column for target rank
cat("Adding a target rank column to the SATC...\n")
# instead we use data.table to set the target rank (after the satc_dt has been ordered)
satc_dt[, target_rank := seq_len(.N), by = .(sample, anchor)]

# satc_dt <- satc_dt %>%
#   group_by(sample, anchor) %>%
#   mutate(target_rank = row_number()) %>%
#   ungroup()

# filter for only the specified target rank
cat(paste0("Filtering to only include targets of rank ", opt$target_rank, "...\n"))
satc_dt <- satc_dt[target_rank == opt$target_rank, ]
satc_dt <- satc_dt[, .(sample, anchor, target, count)]

# satc_dt <- satc_dt %>%
#   filter(target_rank == opt$target_rank) %>%
#   ungroup() %>%
#   select(sample, anchor, target, count)

# get the target length (for filling in Ns later)
target_length <- unique(nchar(head(satc_dt$target)))

# grab the top anchor per cluster as a representative anchor
representative_anchors <- anchor_clusters %>%
  group_by(cluster_id) %>%
  distinct(cluster_id, .keep_all = T) %>%
  ungroup() %>%
  pull(anchor)

# read in the satc and pivot it wider
cat("Creating Wide SATC by merging on clusters and pivoting wider...\n")
wide_satc <- merge(satc_dt, anchor_clusters, by = "anchor", all.x = TRUE)

# drop any anchors where cluster_id is NA in case there were any extra anchors
wide_satc <- as.data.table(wide_satc)
wide_satc <- wide_satc[!is.na(cluster_id), ]

# identify clusters that are not in the satc file
missing_clusters <- setdiff(anchor_clusters$cluster_id, wide_satc$cluster_id)

# add the missing cluster for one sample with the representative anchor and Ns
if (!is_empty(missing_clusters)) {
  missing_df <- data.frame(
    sample = unique(head(wide_satc$sample))[1],
    anchor = representative_anchors[missing_clusters],
    target = strrep("N", target_length),
    count = 0,
    cluster_id = missing_clusters,
    rank = 1
  )

  wide_satc <- rbind(wide_satc, missing_df)
}

wide_satc <- as.data.table(wide_satc)
wide_satc <- wide_satc[order(cluster_id, rank)]
wide_satc <- unique(wide_satc, by = c("sample", "cluster_id"))

wide_satc <- wide_satc %>%
  mutate(seq = str_c(anchor, target, sep = "")) %>%
  select(sample, cluster_id, seq)

wide_satc <- dcast(wide_satc, sample ~ cluster_id, value.var = "seq")
wide_satc <- as.data.frame(wide_satc)

# add the representative anchors to the wide satc
wide_satc <- cbind(wide_satc[1], map2_dfc(
  wide_satc[, 2:ncol(wide_satc)],
  1:length(representative_anchors),
  \(x, y) ifelse(is.na(x), str_c(
    representative_anchors[y], strrep("N", target_length),
    sep = ""
  ), x)
))
wide_satc <- wide_satc %>% ungroup()

# if applying a cluster filter, apply it now
if (!is.null(opt$apply_cluster_filter)) {
  if (grepl("fractionMissing", opt$apply_cluster_filter)) {
    fraction_missing <- as.numeric(str_remove(opt$apply_cluster_filter, "fractionMissing:"))
    cat(paste("Applying percent missing filter with threshold:", fraction_missing, "\n"))
    cluster_filter <- function(cluster_col) {
      total_samples <- nrow(wide_satc)
      missing_samples <- sum(grepl("NNN", cluster_col))
      missing_fraction <- (missing_samples / total_samples)
      return(missing_fraction <= fraction_missing)
    }
    clusters_to_keep <- sapply(wide_satc[, -1], cluster_filter)
    wide_satc <- wide_satc[, c(TRUE, clusters_to_keep)]
  } else if (grepl("binaryTarget", opt$apply_cluster_filter)) {
    cat("Applying binary target filter...\n")
    cluster_filter <- function(cluster_col) {
      unique_seqs <- unique(cluster_col)
      unique_seqs <- unique_seqs[!grepl("NNN", unique_seqs)]
      return(length(unique_seqs) > 2)
    }
    clusters_to_keep <- sapply(wide_satc[, -1], cluster_filter)
    wide_satc <- wide_satc[, c(TRUE, clusters_to_keep)]
  } else {
    stop("Invalid format for --apply_cluster_filter. Must be 'fractionMissing:N' or 'binaryTarget'.")
  }
}

# join the columns together into one sequence and write to a tsv
seqs <- wide_satc %>% unite(seq, -sample, sep = "")
seqs %>% write_tsv(paste0(opt$output_prefix, "_sample_sequences.tsv"))

# write the data to a fasta file
dna <- Biostrings::DNAStringSet(seqs$seq)
names(dna) <- seqs$sample
writeXStringSet(dna, paste0(opt$output_prefix, "_sample_sequences.fasta"))
