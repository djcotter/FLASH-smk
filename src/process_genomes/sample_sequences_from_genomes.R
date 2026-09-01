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
suppressPackageStartupMessages(library(stringdist))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-c", "--cluster_file"), "File with list of anchor clusters", type = "character"),
  make_option(c("-g", "--genome_list"), "File with list of genomes", type = "character"),
  make_option(c("-s", "--genome_files"), "Path to directory with all genome files", type = "character"),
  make_option(c("-o", "--output_prefix"), "Output prefix.", type = "character"),
  make_option(c("--temp_dir"), "Temporary directory to store intermediate files",
    type = "character"
  ),
  make_option(c("--satc_util_bin"), "Path to SATC Util binary folder",
    type = "character", default = ""
  ),
  make_option(c("--anchor_len"), "Anchor length to pass to fafq_filter",
    type = "integer", default = 27
  ),
  make_option(c("--target_len"), "Target length to pass to fafq_filter",
    type = "integer", default = 27
  ),
  make_option(c("--num_cores"), "Number of cores to use", type = "integer", default = 8)
)


# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set up parallel processing
plan(multicore, workers = opt$num_cores)

# check that user specified all files
if (!file.exists(opt$cluster_file) | !file.exists(opt$genome_list) | is.null(opt$genome_files) | is.null(opt$output_prefix)) {
  stop("Must provide cluster file, genome list, genome files, and output prefix")
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

# read in the anchor cluster file
anchor_clusters <- fread(opt$cluster_file,
  header = F, col.names = c("cluster_id", "anchor")
) %>%
  group_by(cluster_id) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# anchor filter
anchor_filter <- file.path(opt$temp_dir, "anchor_filter.txt")
if (!file.exists(anchor_filter)) {
  anchor_clusters %>%
    pull(anchor) %>%
    unique() %>%
    write_lines(anchor_filter)
}

# list all of the genome files in the genome directory, or use explicit sample/file pairs
genome_list <- fread(opt$genome_list, header = F)
if (ncol(genome_list) >= 2) {
  genome_files <- genome_list %>%
    select(genome_name = 1, genome = 2) %>%
    mutate(genome = ifelse(str_detect(genome, "^/"), genome, file.path(opt$genome_files, genome)))
} else {
  colnames(genome_list) <- "genome"
  genome_files <- list.files(opt$genome_files, pattern = ".fna|.fasta|.fa|.fastq|.fq", full.names = T)
  genome_files <- genome_files[str_remove(basename(tools::file_path_sans_ext(str_remove(genome_files, ".gz$"))), "_1") %in% genome_list$genome]

  if (any(str_detect(genome_files, "_1"))) {
    genome_files <- genome_files[str_detect(genome_files, "_1")]
  }

  genome_files <- data.frame(genome = genome_files) %>%
    mutate(genome_name = str_remove(basename(tools::file_path_sans_ext(str_remove(genome, ".gz$"))), "_1"))
}

genome_files <- genome_files %>%
  mutate(
    genome_capitalized = file.path(temp_dir, "capitalized", genome_name),
    genome_oneline = file.path(temp_dir, "oneline", genome_name),
    genome_satc = file.path(temp_dir, "satc", paste0(genome_name, ".satc")),
    genome_dump = file.path(temp_dir, "dumped", paste0(genome_name, ".satc.dump"))
  )

cat("Genome files selected for processing: ", nrow(genome_files), "\n")
if (nrow(genome_files) == 0) {
  stop("No genome files matched --genome_list. Check sample IDs, file basenames, and genome_files path.")
}

system(paste("mkdir -p", file.path(temp_dir, "capitalized")))
system(paste("mkdir -p", file.path(temp_dir, "oneline")))
system(paste("mkdir -p", file.path(temp_dir, "satc")))
system(paste("mkdir -p", file.path(temp_dir, "dumped")))

# declare a satc file for the output of all the dump files
all_genome_file <- file.path(opt$temp_dir, "all_satc_merged.txt")

system(paste("rm", all_genome_file))
if (!file.exists(all_genome_file)) {
  if (any(str_detect(genome_files$genome, "fastq|fq"))) {
    FORMAT <- "fastq"
  } else {
    FORMAT <- "fasta"
  }

  # first process all of the genomes through Biostrings and rewrite them out to file
  future_walk2(genome_files$genome, genome_files$genome_capitalized, \(x, y) {
    Biostrings::readDNAStringSet(x, format = FORMAT) %>% Biostrings::writeXStringSet(y)
  })

  # next make each of the fasta entries occur on only one line
  future_walk2(genome_files$genome_capitalized, genome_files$genome_oneline, \(x, y) {
    system(paste0(
      "awk '/^[>;]/ { if (seq) { print seq }; ", 'seq=""; ', "print } /^[^>;]/ { seq = seq $0 } END { print seq }' ",
      x, " > ", y
    ))
  })


  # fafq filter the fasta files
  future_walk2(genome_files$genome_oneline, genome_files$genome_satc, \(x, y) system(
    paste0(
      file.path(opt$satc_util_bin, "fafq_filter"), " -i ", x, " -o ", y,
      " -d ", anchor_filter,
      " -n 2 --anchor_len ", opt$anchor_len, " --target_len ", opt$target_len
    )
  ))

  # dump the satc files
  future_walk2(genome_files$genome_satc, genome_files$genome_dump, \(x, y) system(
    paste0(file.path(opt$satc_util_bin, "satc_dump "), x, " ", y)
  ))

  # merge them into one file
  walk2(genome_files$genome_dump, genome_files$genome_name, \(x, y) system(
    paste0("cut -f2,3,4 ", x, " | sed 's/^/", y, "\t/g' >> ", all_genome_file)
  ))

  # remove any lines that start or end in [ACTG]
  system(paste(
    "grep -v '^[ACTG]' ", all_genome_file, " | grep -v '[ACTG]$' > ",
    file.path(opt$temp_dir, "all_satc_merged_no_anchor.txt")
  ))
  system(paste("mv", file.path(opt$temp_dir, "all_satc_merged_no_anchor.txt"), all_genome_file))

  system(paste("rm -r", file.path(temp_dir, "capitalized")))
  system(paste("rm -r", file.path(temp_dir, "oneline")))
  system(paste("rm -r", file.path(temp_dir, "satc")))
  system(paste("rm -r", file.path(temp_dir, "dumped")))
}

# read in the dumped satc file
satc_dt <- fread(all_genome_file,
  header = F,
  col.names = c("sample", "anchor", "target", "count"),
  colClasses = c("character", "character", "character", "integer")
)

# grab the top anchor per cluster as a representative anchor
representative_anchors <- anchor_clusters %>%
  group_by(cluster_id) %>%
  distinct(cluster_id, .keep_all = T) %>%
  ungroup() %>%
  pull(anchor)

# read in the satc and pivot it wider
wide_satc <- merge(satc_dt, anchor_clusters, by = "anchor", all.x = TRUE)
wide_satc <- as.data.table(wide_satc)
wide_satc <- wide_satc[order(cluster_id, rank)]
wide_satc <- unique(wide_satc, by = c("sample", "cluster_id"))

representative_anchor_dt <- anchor_clusters %>%
  arrange(cluster_id, rank) %>%
  distinct(cluster_id, .keep_all = T) %>%
  select(cluster_id, anchor, rank)

missing_sample <- representative_anchor_dt %>%
  mutate(sample = "NULLSAMPLE") %>%
  mutate(target = strrep("N", opt$target_len)) %>%
  mutate(count = 1)

requested_missing_samples <- tidyr::crossing(
  sample = unique(genome_files$genome_name),
  representative_anchor_dt
) %>%
  mutate(target = strrep("N", opt$target_len)) %>%
  mutate(count = 0)

wide_satc <- bind_rows(as_tibble(wide_satc), missing_sample, requested_missing_samples) %>%
  mutate(is_missing_placeholder = count == 0) %>%
  arrange(sample, cluster_id, is_missing_placeholder, rank) %>%
  distinct(sample, cluster_id, .keep_all = TRUE) %>%
  select(-is_missing_placeholder)

wide_satc <- wide_satc %>%
  mutate(seq = str_c(anchor, target, sep = "")) %>%
  select(sample, cluster_id, seq)

wide_satc <- as.data.table(wide_satc)
wide_satc <- dcast(wide_satc, sample ~ cluster_id, value.var = "seq")
wide_satc <- as.data.frame(wide_satc)

# write a function that operates on every cluster column of the wide satc
# and aligns the nonrepresentative anchors to the representative anchor
align_to_representative <- function(x, colname, representative_anchor) {
  # x is a one column data frame
  x <- data.frame(sequence = x)
  # grab the unique sequences in the column
  unique_seqs <- x %>%
    group_by(sequence) %>%
    summarise(n = n()) %>%
    ungroup()

  # if there is only one anchor return the original column
  n_anchors <- unique_seqs %>%
    mutate(sequence = substr(sequence, 1, 27)) %>%
    filter(!is.na(sequence)) %>%
    pull(sequence) %>%
    unique() %>%
    length()
  if (n_anchors == 1) {
    colnames(x) <- c(colname)
    return(x)
  }

  # grab the sequence to align to (the most abundant one containing the anchor)
  representative_seq <- unique_seqs %>%
    filter(str_detect(sequence, representative_anchor)) %>%
    filter(n == max(n)) %>%
    head(1) %>%
    pull(sequence)

  # format the unique sequences as a character vector
  unique_seqs <- unique_seqs %>% pull(sequence)

  # function to perform alignment on a unique set of sequences
  perform_alignment <- function(seq, ref_seq) {
    # Calculate the Levenshtein edit moves
    lev_moves <- attr(adist(seq, ref_seq, counts = TRUE), "trafos")[[1]]

    # Initialize empty aligned sequence
    aligned_seq <- ""

    # Process each move backwards
    # split lev moves into a vector one letter each
    lev_moves <- strsplit(lev_moves, "")[[1]]
    # perform the moves on the seq to get the aligned seq (skip S)
    # if there is an M keep the corresponding letter from the seq
    # if there is a D, delete the corresponding letter from the seq
    # if there is a I, insert an N in the aligned seq
    # if there is an S, do nothing
    for (move in lev_moves) {
      if (move == "M") {
        aligned_seq <- paste0(aligned_seq, substr(seq, 1, 1))
        seq <- substr(seq, 2, nchar(seq))
      } else if (move == "D") {
        seq <- substr(seq, 2, nchar(seq))
      } else if (move == "I") {
        aligned_seq <- paste0(aligned_seq, "N")
      }
    }

    # first trim the leading Ns from the aligned seq
    aligned_seq <- substr(aligned_seq, 4, nchar(aligned_seq))

    # trim or pad with Ns the aligned sequence to be the same length as the reference
    if (nchar(aligned_seq) < nchar(ref_seq)) {
      aligned_seq <- str_pad(aligned_seq, nchar(ref_seq), pad = "N", side = "right")
    } else if (nchar(aligned_seq) > nchar(ref_seq)) {
      aligned_seq <- substr(aligned_seq, 1, nchar(ref_seq))
    }

    # Return the final aligned sequence
    return(aligned_seq)
  }

  # Function to align all sequences against the reference using the above alignment logic
  align_sequences <- function(seqs, ref_seq) {
    aligned_seqs <- sapply(seqs, function(seq) perform_alignment(seq, ref_seq))
    return(aligned_seqs)
  }

  # align the unique sequences to the representative anchor
  aligned_seqs <- align_sequences(unique_seqs, representative_seq)

  # mutate the column in x mapping each unique sequence to its aligned sequence
  x <- x %>% mutate(sequence = map_chr(sequence, \(x) ifelse(is.na(x), NA, aligned_seqs[x])))
  colnames(x) <- c(colname)

  return(x)
}

# # apply the alignment function to each cluster column
# wide_satc <- cbind(wide_satc[1],
#                    future_pmap_dfc(list(wide_satc[,2:ncol(wide_satc)],
#                                         colnames(wide_satc[,2:ncol(wide_satc)]),
#                                         representative_anchors),
#                                    \(x,y,z) align_to_representative(x, y, z)))

# add the representative anchors with Ns to the wide satc where there are NAs
wide_satc <- base::cbind(
  wide_satc[1],
  map2_df(
    wide_satc[, 2:ncol(wide_satc)],
    1:length(representative_anchors),
    \(x, y) ifelse(is.na(x), str_c(representative_anchors[y], strrep("N", opt$target_len), sep = ""), x)
  )
)
wide_satc <- wide_satc %>% ungroup()

wide_satc <- wide_satc %>% filter(sample != "NULLSAMPLE")

# join the columns together into one sequence and write to a tsv
seqs <- wide_satc %>% unite(seq, -sample, sep = "")
seqs %>% write_tsv(paste0(opt$output_prefix, "_sample_sequences.tsv"))

# write the data to a fasta file
dna <- Biostrings::DNAStringSet(seqs$seq)
names(dna) <- seqs$sample
writeXStringSet(dna, paste0(opt$output_prefix, "_sample_sequences.fasta"))
