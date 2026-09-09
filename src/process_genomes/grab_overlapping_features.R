# grab_top_embeddings_by_variance.R
# Daniel Cotter
# 2024-09-18

# This script takes in one embeddings file and the ordering file and outputs a new
# embeddings tsv with samples as rows and the top 10 embeddings by variance per cluster as columns.


## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(resample))
suppressPackageStartupMessages(library(furrr))


## parse arguments --------
# define command line arguments
# define command line arguments
option_list <- list(
  make_option(c("-e", "--embeddings"), help = "Embeddings file", type = "character"),
  make_option(c("-o", "--ordering"), help = "Ordering file", type = "character"),
  make_option(c("--original_adelie_output"), help = "Path to nonzero coefficients.", type = "character"),
  make_option(c("--original_feather"), help = "Original feather file", type = "character"),
  make_option(c("-p", "--output_prefix"), help = "Output prefix.", type = "character"),
  make_option(c("--temp_dir"),
    help = "Temporary directory to store intermediate files",
    type = "character"
  ),
  make_option(c("--num_threads"),
    help = "Number of threads to use for parallel operations",
    type = "integer", default = 1
  ),
  make_option(c("--normalized_embeddings"),
    help = "Whether to normalize the embeddings before calculating variance",
    action = "store_true", default = FALSE
  )
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that user specified all files
if (!file.exists(opt$embeddings) | !file.exists(opt$ordering) | is.null(opt$output_prefix) | !file.exists(opt$original_feather)) {
  stop("Must provide embeddings, ordering, original feather file, and output prefix")
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

# define future plans
setDTthreads(opt$num_threads)
plan(multisession, workers = opt$num_threads)
options(future.globals.maxSize = 16000 * 1024^2)

## print a summary of the arguments
cat("\n####################\n")
cat("Running grab_overlapping_features.R with the following arguments:\n")
cat("Ordering file: ", opt$ordering, "\n")
cat("Embeddings file: ", opt$embeddings, "\n")
cat("Original feather file: ", opt$original_feather, "\n")
cat("Output file: ", opt$output_prefix, "\n")
cat("Temporary directory: ", temp_dir, "\n")
cat("####################\n\n")

## load data --------
# load embeddings
cat("\nLoading embeddings...\n")
# copy the embeddings file to the temp directory to speed up I/O
embeddings_temp <- file.path(temp_dir, "genomes_raw_embeddings_temp.tsv")
system(paste("cp", opt$embeddings, embeddings_temp))
embeddings <- fread(embeddings_temp, header = F)
colnames(embeddings) <- c("kmer", paste0("embedding_", 1:(ncol(embeddings) - 1)))

# load the ordering file
cat("Loading the ordering file...\n")
# copy the ordering file to the temp directory to speed up I/O
ordering_temp <- file.path(temp_dir, "genomes_ordering_temp.tsv")
system(paste("cp", opt$ordering, ordering_temp))
ordering <- fread(ordering_temp,
  header = F, sep = "\t",
  col.names = c("sample_name", "seq", "kmer", "start", "end"),
  # first column should ALWAYS be a character
  colClasses = c("character", "character", "character", "integer", "integer")
)

ordering <- ordering %>% select(sample_name, kmer)

ordering <- ordering %>%
  group_by(sample_name) %>%
  mutate(cluster = row_number()) %>%
  mutate(cluster = cluster - 1) %>%
  ungroup()

clusters <- ordering %>%
  group_by(cluster) %>%
  group_split()

join_and_write_clusters <- function(cluster_df, filename, all_embeddings) {
  cluster_df %>%
    select(-cluster) %>%
    left_join(all_embeddings, by = "kmer") %>%
    select(-kmer) %>%
    fwrite(file = filename, nThread = 1, col.names = T)
}

temp_embeddings_dir <- file.path(temp_dir, "genome_embeddings_per_cluster/")
system(paste("mkdir -p", temp_embeddings_dir))
cluster_files <- paste0(temp_embeddings_dir, "embeddings_cluster_", 0:(length(clusters) - 1), ".csv")

cat("Writing all clusters and their embeddings out to file in: ", temp_embeddings_dir)

future_walk2(clusters, cluster_files, \(x, y) join_and_write_clusters(x, y, all_embeddings = embeddings))


# cluster_to_kmer_mapping <- ordering %>% select(cluster, kmer) %>% distinct()
# cluster_to_kmer_mapping <- cluster_to_kmer_mapping %>%
#   arrange(cluster) %>% group_by(cluster) %>%
#   summarise(kmers = str_c(kmer, collapse=",")) %>% ungroup()


cat("Formatting the embeddings for downstream use...\n")
# load the original embeddings file and grab the column names
original_dt <- feather::read_feather(opt$original_feather)
original_cols <- colnames(original_dt)

# load the nonzero coefficients column
nonzero_coef_cols <- fread(opt$original_adelie_output) %>%
  select(feature) %>%
  mutate(cluster = str_remove(feature, "_embedding_\\d+")) %>%
  distinct(cluster) %>%
  pull(cluster)

# filter the original cols to only contain the "NEW" columns
new_col_ordering <- original_cols[grepl("embedding", original_cols) & str_remove(original_cols, "_embedding_\\d+") %in% nonzero_coef_cols]
new_col_ordering <- c("sample_name", new_col_ordering)

# filter the original dt to only contain the new col ordering
original_dt <- original_dt %>% select(all_of(new_col_ordering))

# write out the original dt with the minimized column selection
new_og_dt_path <- paste0(opt$output_prefix, "_ORIGINAL.feather")
arrow::write_feather(original_dt, new_og_dt_path)

# loop through the cluster files and grab columns that match the original column names
grab_matching_columns <- function(in_file, normalized, original_cols) {
  cluster_num <- str_extract(in_file, "cluster_(\\d+).csv", group = 1) %>% as.integer()
  temp_dt <- fread(in_file, header = T, nThread = 1, colClasses = c(sample_name = "character")) %>%
    select(sample_name, starts_with("embedding")) # first filter for only one cluster
  colnames(temp_dt) <- ifelse(grepl("embedding", colnames(temp_dt)),
    yes = paste0("cluster_", cluster_num, "_", colnames(temp_dt)),
    no = colnames(temp_dt)
  )
  if (normalized) {
    temp_dt <- temp_dt %>% mutate(across(starts_with("cluster"), scale))
  }
  cols_to_keep <- intersect(colnames(temp_dt), original_cols)
  temp_dt <- temp_dt %>%
    select(sample_name, all_of(cols_to_keep)) %>%
    arrange(sample_name)
  return(temp_dt)
}

cat("Grabbing matching columns per cluster\n")
top_col_dt <- future_map(
  cluster_files,
  \(x) grab_matching_columns(x, normalized = opt$normalized_embeddings, new_col_ordering)
) %>%
  purrr::reduce(full_join, by = "sample_name") %>%
  arrange(sample_name)

# ensure the column names are in the same order as the original embeddings
top_col_dt <- top_col_dt %>% select(all_of(new_col_ordering))

# write out the embeddings matrix to a temp feather file
embeddings_feather <- paste0(opt$output_prefix, "_GENOMES.feather")
cat("Writing matching column embeddings to ", embeddings_feather, "\n")
arrow::write_feather(top_col_dt, embeddings_feather)


# cat("Grabbing matching columns per cluster\n")
# top_col_dt <- future_map_dfc(
#   cluster_files,
#   \(x) grab_matching_columns(x, normalized = opt$normalized_embeddings, original_cols)
# )

# # ensure the column names are in the same order as the original embeddings
# top_col_dt <- top_col_dt %>% select(all_of(original_cols))

# # write out the embeddings matrix to a temp feather file
# embeddings_feather <- opt$output
# cat("Writing matching column embeddings to ", embeddings_feather, "\n")
# arrow::write_feather(top_col_dt, embeddings_feather)
