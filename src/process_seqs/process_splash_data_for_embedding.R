# process_splash_data_for_embedding.R
# Daniel Cotter
# 2024-09-13

# This script takes in a path to SPLASH SATC files as well as a list of anchors
# and a list of anchor clusters. It then dumps the anchors from the SATC files
# and reformats a sequence for each sample that can be used in downstream
# analyses. The output is a fasta file and a tsv file with the sample id and
# the sequence.

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
    c("-s", "--satc_files"),
    "Path to SPLASH SATC file directory",
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
    help = "merge the first 2 cols of the SATC to create sample barcode pairs",
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

# check that user specified all files
if (!file.exists(opt$anchor_file) |
  !file.exists(opt$cluster_file) |
  is.null(opt$satc_files) |
  is.null(opt$output_prefix) | !file.exists(opt$id_mapping)) {
  stop(paste(
    "Must provide anchor file, cluster file,",
    "satc files, id mapping file, and output prefix"
  ))
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

system(paste0("rm ", temp_dir, "/*"))
system(paste0("rm ", temp_dir, "/filtered/*"))
system(paste0("rm ", temp_dir, "/dumped/*"))

# read in the anchor cluster file
anchor_clusters <- fread(opt$cluster_file,
  header = FALSE,
  col.names = c("cluster_id", "anchor")
) %>%
  group_by(cluster_id) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# list all of the .satc files in the result_satc folder
satc_files <- list.files(opt$satc_files,
  pattern = "bin\\d+.satc", full.names = TRUE
)

satc_files <- data.frame(satc_file = satc_files) %>%
  mutate(satc_filter = gsub(
    ".satc",
    "filtered.satc",
    file.path(temp_dir, "filtered", basename(satc_file))
  )) %>%
  mutate(satc_dump = gsub(
    ".satc",
    ".satc.dump",
    file.path(temp_dir, "dumped", basename(satc_file))
  ))

# create temp dir for the filtered satc files
system(paste("mkdir -p", file.path(temp_dir, "filtered")))

# create temp dir for dumped satc files
system(paste("mkdir -p", file.path(temp_dir, "dumped")))

# declare a satc file for the output of all the dump files
all_satc_file <- file.path(temp_dir, "all_satc_merged.txt")

# first filter the satc files to only include anchors in the anchor file and
# up to 1 plus the specified target rank, then dump them using the satc_util_bin
cat("Filtering SATC files...\n")
future_walk2(satc_files$satc_file, satc_files$satc_filter, function(x, y) {
  system(
    paste(
      file.path(opt$satc_util_bin, "satc_filter"),
      "-i",
      x,
      "-o",
      y,
      "-d",
      opt$anchor_file,
      "-n",
      opt$target_rank + 1 # grab all targets up to the specified rank + 1
    )
  )
})

# dump the satc files
cat("Dumping SATC files...\n")
future_walk2(satc_files$satc_filter, satc_files$satc_dump, function(x, y) {
  system(
    paste0(
      file.path(opt$satc_util_bin, "satc_dump"),
      " --anchor_list ",
      opt$anchor_file,
      " --sample_names ",
      opt$id_mapping,
      " ",
      x,
      " ",
      y
    )
  )
})

# merge them into one file
cat("Merging dumped SATC files...\n")
system(paste("rm", all_satc_file))
walk(satc_files$satc_dump, \(x) system(paste("cat", x, ">>", all_satc_file)))

# remove any lines that start or end in [ACTG]
cat("Removing invalid anchor-target pairs...\n")
system(paste(
  "grep -v '^[ACTG]{3}' ",
  all_satc_file,
  " | grep -v '[ACTG]$' > ",
  file.path(opt$temp_dir, "all_satc_merged_no_anchor.txt")
))

system(paste(
  "mv",
  file.path(opt$temp_dir, "all_satc_merged_no_anchor.txt"),
  all_satc_file
))

# if single cell, merge the first two columns to create sample barcode pairs
if (opt$single_cell) {
  cat("Single cell mode: merging first two cols...\n")
  temp_merged_file <- file.path(opt$temp_dir, "all_satc_SC_merged_columns.txt")
  system(
    paste(
      "awk 'BEGIN {OFS=\"\\t\"} {print $1\"_\"$2, $3, $4, $5}'",
      all_satc_file,
      ">",
      temp_merged_file
    )
  )
  system(paste("mv", temp_merged_file, all_satc_file))
}

# undump the satc files
cat("Undumping merged SATC file...\n")
all_satc_undumped <- file.path(opt$temp_dir, "all_satc_merged.undumped.txt")
all_satc_temp_mapping <- file.path(
  opt$temp_dir,
  "all_satc_merged.temp_mapping.txt"
)

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
cat("Filtering merged SATC file...\n")
cat(paste("Using target rank:", opt$target_rank, "\n"))
all_satc_filtered <- file.path(opt$temp_dir, "all_satc.filtered.txt")
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
    opt$target_rank + 1 # grab all targets up to the specified rank + 1
  )
)


# redump the satc file
cat("Redumping filtered SATC file for downstream processing...\n")
all_satc_filtered_dump <- file.path(opt$temp_dir, "all_satc.filtered.dump")
system(
  paste(
    file.path(opt$satc_util_bin, "satc_dump"),
    "--sample_names",
    all_satc_temp_mapping,
    all_satc_filtered,
    all_satc_filtered_dump
  )
)

# read in the dumped satc file
satc_dt <- fread(
  all_satc_filtered_dump,
  header = FALSE,
  col.names = c("sample", "anchor", "target", "count")
)

# now filter for ONLY the specified target rank
satc_dt <- satc_dt[order(sample, anchor, -count)]

# add a column for target rank
satc_dt <- satc_dt %>%
  group_by(sample, anchor) %>%
  mutate(target_rank = row_number()) %>%
  ungroup()

# filter for only the specified target rank
satc_dt <- satc_dt %>%
  filter(target_rank == opt$target_rank) %>%
  ungroup() %>%
  select(sample, anchor, target, count)

# get the target length (for filling in Ns later)
target_length <- unique(nchar(satc_dt$target))

# grab the top anchor per cluster as a representative anchor
representative_anchors <- anchor_clusters %>%
  group_by(cluster_id) %>%
  distinct(cluster_id, .keep_all = TRUE) %>%
  ungroup() %>%
  pull(anchor)

# read in the satc and pivot it wider
wide_satc <- merge(satc_dt, anchor_clusters, by = "anchor", all.x = TRUE)

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
  seq_along(representative_anchors),
  function(x, y) {
    ifelse(is.na(x), str_c(
      representative_anchors[y], strrep("N", target_length),
      sep = ""
    ), x)
  }
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
