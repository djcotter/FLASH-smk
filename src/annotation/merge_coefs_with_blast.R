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

extract_blast_species <- function(x) {
  vapply(as.character(x), function(value) {
    if (is.na(value) || !nzchar(value)) {
      return(NA_character_)
    }

    parts <- strsplit(value, "[", fixed = TRUE)[[1]]
    if (length(parts) < 2) {
      return(NA_character_)
    }

    candidate <- parts[[length(parts)]]
    closing_bracket <- regexpr("]", candidate, fixed = TRUE)[[1]]
    if (closing_bracket <= 1) {
      return(NA_character_)
    }

    species <- stringr::str_squish(substr(candidate, 1, closing_bracket - 1))
    if (nzchar(species)) species else NA_character_
  }, character(1))
}

clean_blast_species_field <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x[x == "" | toupper(x) %in% c("NA", "N/A", "NONE", "NAN")] <- NA_character_
  x
}

ensure_columns <- function(tbl, cols) {
  for (col in cols) {
    if (!col %in% colnames(tbl)) {
      tbl[[col]] <- rep(NA, nrow(tbl))
    }
  }
  tbl
}

annotation_cols <- c(
  "query", "subject", "sacc", "identity", "alignment_length", "mismatches",
  "gap_opens", "q_start", "q_end", "s_start", "s_end", "sstrand", "evalue",
  "qcovs", "qframe", "sgi", "slen", "staxids", "sscinames", "stitle",
  "species_origin", "NCBI_protein_accession", "UniProt_accession", "method", "GO",
  "features", "features_10000_window", "features_all", "blast_mode", "blast_scope"
)

compactor_cols <- c(
  "compactor_annotation", "compactor_query", "compactor_sequence",
  "compactor_blast_query_sequence", "compactor_length", "compactor_exact_support",
  "compactor_support", "compactor_expected_read_count", "compactor_extender_specificity",
  "compactor_num_extended", "compactor_support_threshold", "compactor_selection_reason",
  "compactor_raw_annotation", "compactor_species", "compactor_staxids",
  "compactor_sscinames", "compactor_blast_subject_id", "compactor_blast_accession",
  "compactor_ncbi_protein_accession", "compactor_uniprot_accession",
  "blast_subject_id_origin", "blast_accession_origin", "blast_annotation_source_id"
)

# Read in the data
annotations <- fread(opt$blast_annotations, header = TRUE, sep = "\t", nThread = 60)

if (str_detect(opt$blast_annotations, "blastp|swissprot")) {
  annotations <- ensure_columns(annotations, c(annotation_cols, compactor_cols))
  annotations <- annotations %>%
    select(any_of(c(annotation_cols, compactor_cols)), everything()) %>%
    mutate(blast_mode = coalesce(as.character(blast_mode), "blastx")) %>%
    mutate(cluster = str_extract(query, "(^.*cluster_\\d+|\\w+_kmer_\\d+)_", group = 1))
  if (!"sscinames" %in% colnames(annotations)) {
    annotations$sscinames <- NA_character_
  }
  if (!"species_origin" %in% colnames(annotations)) {
    annotations <- annotations %>% mutate(species_origin = coalesce(clean_blast_species_field(sscinames), extract_blast_species(stitle)))
  } else {
    annotations <- annotations %>% mutate(species_origin = coalesce(clean_blast_species_field(species_origin), clean_blast_species_field(sscinames), extract_blast_species(stitle)))
  }

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
      filter(nchar(translated_sequence) > 5)
    if (nrow(aa_dt_filtered) == 0) {
      return(aa_dt %>% select(-translated_sequence) %>% mutate(aligned_sequence = NA))
    } else if (nrow(aa_dt_filtered) == 1) {
      return(aa_dt_filtered %>% rename(aligned_sequence = translated_sequence))
    }
    aas <- Biostrings::AAStringSet(aa_dt_filtered$translated_sequence %>% str_remove("\\*.+$"))
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
    aligned_seqs <- data.frame(query = names(aligned), aligned_sequence = aligned)
    return(aligned_seqs)
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
  annotations <- ensure_columns(annotations, c(annotation_cols, compactor_cols))
  annotations <- annotations %>%
    select(any_of(c(annotation_cols, compactor_cols)), everything()) %>%
    mutate(blast_mode = coalesce(as.character(blast_mode), "blastn")) %>%
    mutate(cluster = str_extract(query, "(^.*cluster_\\d+|\\w+_kmer_\\d+)_", group = 1))
  if (!"sscinames" %in% colnames(annotations)) {
    annotations$sscinames <- NA_character_
  }
  if (!"species_origin" %in% colnames(annotations) && "stitle" %in% colnames(annotations)) {
    annotations <- annotations %>% mutate(species_origin = coalesce(clean_blast_species_field(sscinames), extract_blast_species(stitle)))
  } else if ("species_origin" %in% colnames(annotations) && "stitle" %in% colnames(annotations)) {
    annotations <- annotations %>% mutate(species_origin = coalesce(clean_blast_species_field(species_origin), clean_blast_species_field(sscinames), extract_blast_species(stitle)))
  }
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
  for (delimiter in c("[", "]", "(", ")")) {
    coef_string <- gsub(delimiter, " ", coef_string, fixed = TRUE)
  }
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
