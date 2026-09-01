# Load necessary libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

# Define command line options
option_list <- list(
  make_option(c("-a", "--blast_annotations"),
    type = "character", default = NULL,
    help = "Path to the annotations CSV file", metavar = "character"
  ),
  make_option(c("-c", "--coefficients"),
    type = "character", default = NULL,
    help = "Path to the coefficients CSV file", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "Path to the output CSV file", metavar = "character"
  ),
  make_option(c("-t", "--translation_table"),
    type = "integer", default = 11,
    help = "Path to the translation table file. Default is 11 for bacteria.",
    metavar = "integer"
  )
)


# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$blast_annotations) || is.null(opt$coefficients) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All arguments must be supplied", call. = FALSE)
}

ensure_columns <- function(tbl, cols) {
  for (col in cols) {
    if (!col %in% colnames(tbl)) {
      tbl[[col]] <- rep(NA, nrow(tbl))
    }
  }
  tbl
}

# Read in the data
annotations <- fread(opt$blast_annotations, header = TRUE, sep = "\t", nThread = 60)

if (str_detect(opt$blast_annotations, "blastp|swissprot")) {
  blastp_cols <- c("query", "evalue", "identity", "qcovs", "qframe", "stitle",
                   "NCBI_protein_accession", "UniProt_accession", "method", "GO")
  annotations <- ensure_columns(annotations, blastp_cols)
  annotations <- annotations %>%
    select(any_of(blastp_cols)) %>%
    mutate(cluster = str_extract(query, "(^.*cluster_\\d+|\\w+_kmer_\\d+)_", group = 1))

  # we also add on the translated sequence for as a column. The query column contains cluster_X_{Sequence} and the qframe column contains the frame
  # we extract {sequence} and translate it using the qframe
  TRANSLATION_TABLE <- opt$translation_table

  sequence_dt <- annotations %>%
    select(query, qframe) %>%
    distinct() %>%
    filter(!is.na(query)) %>%
    mutate(sequence = str_remove(query, "^.*cluster_\\d+_|\\w+_kmer_\\d+_")) %>%
    mutate(
      qframe_numeric = suppressWarnings(as.numeric(qframe)),
      reverse = ifelse(!is.na(qframe_numeric) & qframe_numeric < 0, TRUE, FALSE),
      frame = abs(qframe_numeric) - 1
    ) %>%
    distinct()

  if (nrow(sequence_dt) == 0) {
    sequence_dt <- sequence_dt %>%
      mutate(translated_sequence = character(), aligned_sequence = character())
  } else {
  duplicates <- sequence_dt %>%
    group_by(query) %>%
    filter(n() > 1) %>%
    mutate(cluster = str_remove(query, "_[ACTGNactgn]+$")) %>%
    pull(cluster) %>%
    unique()

  translate_sequences <- function(sequences, frames, reverse, translation_table) {
    translated_seqs <- tolower(sequences)
    translated_seqs <- map(translated_seqs, seqinr::s2c)

    translated <- pmap_vec(list(x = translated_seqs, y = frames, z = reverse), function(x, y, z) {
      out <- ifelse(!is.na(y),
        seqinr::c2s(seqinr::translate(seq = x, frame = y, sens = ifelse(z, "R", "F"), numcode = as.integer(translation_table))),
        NA
      )

      return(out)
    })

    return(translated)
  }

  translated_sequences <- translate_sequences(sequence_dt$sequence, sequence_dt$frame, sequence_dt$reverse, TRANSLATION_TABLE)
  if (is.null(translated_sequences) || length(translated_sequences) != nrow(sequence_dt)) {
    translated_sequences <- rep(NA_character_, nrow(sequence_dt))
  }
  sequence_dt <- sequence_dt %>% mutate(translated_sequence = translated_sequences)

  # also add a column removing anything after a stop codon and aligning the translated_sequences
  aa_temps <- sequence_dt %>%
    mutate(cluster = str_remove(query, "_[ACTGNactgn]+$")) %>%
    filter(!(cluster %in% duplicates)) %>%
    select(query, cluster, translated_sequence) %>%
    group_by(cluster) %>%
    group_split(.keep = FALSE)

  align_cluster <- function(aa_dt) {
    aa_dt_filtered <- aa_dt %>%
      mutate(translated_sequence = str_remove(translated_sequence, "\\*.+$")) %>%
      filter(!is.na(translated_sequence), nchar(translated_sequence) > 5)

    if (nrow(aa_dt_filtered) == 0) {
      return(aa_dt %>%
              select(-translated_sequence) %>%
              mutate(aligned_sequence = NA_character_))
    } else if (nrow(aa_dt_filtered) == 1) {
      return(aa_dt_filtered %>%
              transmute(query = query, aligned_sequence = translated_sequence))
    }

    aas <- Biostrings::AAStringSet(aa_dt_filtered$translated_sequence)
    names(aas) <- aa_dt_filtered$query

    aligned <- tryCatch(
      msa::msa(aas, order = "input"),
      error = function(e) {
        message("Skipping amino-acid MSA for cluster ",
                aa_dt_filtered$cluster[[1]],
                " because msa failed: ", conditionMessage(e))
        return(NULL)
      }
    )

    if (is.null(aligned)) {
      return(aa_dt_filtered %>%
              transmute(query = query, aligned_sequence = translated_sequence))
    }

    aligned <- Biostrings::AAStringSet(aligned)
    data.frame(query = names(aligned), aligned_sequence = as.character(aligned))
  }

  aa_aligned <- map(aa_temps, align_cluster)
  aa_aligned <- bind_rows(aa_aligned)
  if (!"query" %in% colnames(aa_aligned)) {
    aa_aligned <- tibble(query=character(), aligned_sequence=character())
  }
  sequence_dt <- sequence_dt %>% left_join(aa_aligned, by = "query")
  sequence_dt <- sequence_dt %>% select(query, qframe, translated_sequence, aligned_sequence)
  }

  # now bind it all together
  annotations <- annotations %>% left_join(sequence_dt, by = c("query", "qframe"))
} else {
  blast_cols <- c("query", "evalue", "identity", "qcovs", "features", "features_10000_window")
  annotations <- ensure_columns(annotations, blast_cols)
  annotations <- annotations %>%
    select(any_of(blast_cols), contains("window")) %>%
    mutate(cluster = str_extract(query, "(^.*cluster_\\d+|\\w+_kmer_\\d+)_", group = 1))
}



coefficients <- fread(opt$coefficients, header = TRUE)
coefficients <- ensure_columns(coefficients, c("metadata_category", "feature", "coefficients"))
coefficients <- coefficients %>% mutate(cluster = str_extract(feature, "(^.*cluster_\\d+|\\w+_kmer_\\d+)_", group = 1))

# because we have run grouped elastic net there will be multiple rows per cluster. We want to keep the one with the highest absolute coefficient
get_max_coef <- function(coef_string) {
  coef_string <- as.character(coef_string)
  if (length(coef_string) == 0 || is.na(coef_string) || !nzchar(coef_string)) {
    return(NA_real_)
  }
  coef_string <- gsub("\\[|\\]|c\\(|\\)", "", coef_string)
  coefs <- suppressWarnings(as.numeric(unlist(strsplit(coef_string, "[,;[:space:]]+"))))
  coefs <- coefs[is.finite(coefs)]
  if (length(coefs) == 0) {
    return(NA_real_)
  }
  return(max(abs(coefs)))
}

coefficients <- coefficients %>%
  rowwise() %>%
  mutate(max_coef = get_max_coef(.data$coefficients)) %>%
  arrange(cluster, -max_coef) %>%
  group_by(cluster) %>%
  slice(1) %>%
  ungroup()

merged_data <- coefficients %>% full_join(annotations, by = "cluster", relationship = "many-to-many")

merged_data <- merged_data %>%
  arrange(metadata_category, -abs(max_coef), -identity) %>%
  select(-max_coef)

# relocate sequences in blastp output
if (str_detect(opt$blast_annotations, "blastp")) {
  if (!"translated_sequence" %in% colnames(merged_data)) {
    merged_data$translated_sequence <- rep(NA_character_, nrow(merged_data))
  }
  if (!"aligned_sequence" %in% colnames(merged_data)) {
    merged_data$aligned_sequence <- rep(NA_character_, nrow(merged_data))
  }
  merged_data <- merged_data %>%
    relocate(translated_sequence, aligned_sequence,
      .after = "NCBI_protein_accession"
    )
}


# Write the merged data to a new CSV file
write_tsv(merged_data, opt$output, col_names = T, quote = "needed")
