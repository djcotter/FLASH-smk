suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(msa))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))


option_list <- list(
  make_option(c("--nonzero_annotations"), type = "character", default = NULL, 
              help = "Path to the nonzero annotations tsv file", metavar = "character"),
  make_option(c("--clusters"), type = "character", default = NULL, 
              help = "Path to the clusters tsv file", metavar = "character"),
  make_option(c("--feather_file"), type = "character", default = NULL, 
              help = "Path to the X matrix feather file", metavar = "character"),
  make_option(c("--sample_seqs"), type = "character", default = NULL, 
              help = "Path to the sample sequences file", metavar = "character"),
  make_option(c("--metadata"), type = "character", default = NULL, 
              help = "Path to the metadata tsv file", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL, 
              help = "Path to set of output plots", metavar = "character"),
  make_option(c("--products"), type= "logical", default=FALSE, action="store_true",
              help = "default to using products for column names instead of genes"),
  make_option(c("--num_hits"), type="numeric", default=10,
              help = "num nonzero coefficients to plot", metavar = "numeric")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# Check if all required arguments are provided
if (is.null(opt$nonzero_annotations) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All arguments must be supplied", call. = FALSE)
}

# set known_causes to be empty (can be changed for interactive experimentation on specific datasets)
known_causes = "NNNNNNNNNNNNNNN"
max_features_per_heatmap = 40

# # testing
# setwd("/oak/stanford/groups/horence/dcotter1/projects/metaSPLASH_pipeline")
# opt$nonzero_annotations = "results/canTrop-AzoleResistance-PRJNA946688/filter1/shiftDist-levFilter/hyena/normalized/canTrop-AzoleResistance-PRJNA946688_hyena_adelie_results_top20000_k54_s54_nonzero_coefficients_blastp_annotated.tsv"
# opt$clusters = "results/canTrop-AzoleResistance-PRJNA946688/filter1/shiftDist-levFilter/canTrop-AzoleResistance-PRJNA946688_sequences_per_cluster_top20000-clusters_k54_s54.tsv"
# opt$feather = "/scratch/users/dcotter1/metaSPLASH_workflows_v2/canTrop-AzoleResistance-PRJNA946688/canTrop-AzoleResistance-PRJNA946688_hyena_top_variance_features_for_glmnet_filter1_shiftDist-levFilter_top20000_k54_s54_normalized.feather"
# opt$sample_seqs = "/scratch/users/dcotter1/metaSPLASH_workflows_v2/canTrop-AzoleResistance-PRJNA946688/canTrop-AzoleResistance-PRJNA946688_prepared_sequences_filter1_shiftDist-levFilter_top20000_sample_sequences.tsv"
# opt$metadata = "/oak/stanford/groups/horence/dcotter1/utility_files/metadata/metaSPLASH_metadata/candida_tropicalis_PRJNA946688_cleaned_metadata.tsv"
# opt$output = "/oak/stanford/groups/horence/dcotter1/share/250506/test_eFac_more_blast_hits_out.pdf"


filename = data.frame(path=opt$nonzero_annotations)

filename <- filename %>% 
  mutate(num_clusters = str_extract(path, "top(\\d+)", group=1)) %>%
  mutate(path=dirname(gsub("^results/", "", path))) %>%
  mutate(path=gsub("/","_",path)) %>%
  dplyr::rename(paramater_set=path) %>%
  mutate(paramater_set=str_replace(paramater_set, "_", "/")) %>% 
  separate(paramater_set, into=c("dataset", "paramater_set"), sep="/") %>%
  mutate(model=str_extract(paramater_set,
                           'hyenaHG38_normalized|hyenaHG38_unnormalized|hyenaMarlowe_normalized|hyenaMarlowe_unnormalized|esm_normalized|esm_unnormalized|hyena_normalized|hyena_unnormalized|ohe')) %>%
  mutate(filter = str_extract(paramater_set, "(filter\\d)_",group=1)) %>% 
  mutate(cluster_approach = str_extract(paramater_set, "filter\\d_([A-Za-z-2]+)_", group=1)) 

paramaters <- filename %>% pivot_longer(everything(), names_to="paramater", values_to="value") %>% deframe()


# Define Function 
get_max_abs_value <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    max(abs(nums), na.rm = TRUE)
  })
}


get_first_coef <- function(x) {
  sapply(x, function(str) {
    nums <- as.numeric(strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]])
    nums[1]
  })
}

get_first_class <- function(x) {
  sapply(x, function(str) {
    classes <- strsplit(gsub("^\\[|\\]$", "", str), ",")[[1]]
    classes[1]
  })
}

clean_blast_label <- function(x) {
  boilerplate <- regex("Gnomon|Derived by automated computational analysis|Supporting evidence includes|coverage of the annotated genomic feature|support for all annotated introns", ignore_case=TRUE)
  vapply(as.character(x), function(value) {
    if (is.na(value)) return(NA_character_)
    parts <- unlist(str_split(value, ";+"), use.names=FALSE)
    parts <- parts[!str_detect(parts, boilerplate)]
    value <- str_trim(str_replace_all(paste(parts, collapse=";"), "\\s+", " "))
    value <- str_replace_all(str_replace_all(value, "LOC\\d+[- ]*", ""), "\\s*[,;]\\s*$", "")
    if (nchar(value) == 0) NA_character_ else value
  }, character(1), USE.NAMES=FALSE)
}

extract_feature_qualifier <- function(features, qualifier) {
  pattern <- paste0("'", qualifier, "': \\['([^']+)'\\]")
  str_extract(features, pattern, group=1) %>% clean_blast_label()
}

choose_feature_label <- function(products, genes, prefer_products = FALSE) {
  products <- clean_blast_label(products)
  genes <- clean_blast_label(genes)
  genes_are_loc <- !is.na(genes) & str_detect(genes, "^LOC\\d+$")
  use_products <- prefer_products | genes_are_loc | is.na(genes) | nchar(genes) < 2
  label <- ifelse(use_products & !is.na(products) & nchar(products) > 1, products, genes)
  label <- ifelse((is.na(label) | nchar(label) < 2) & !is.na(products), products, label)
  clean_blast_label(label)
}

parse_vector_text <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "^\\s*\\[|\\]\\s*$", "")
  x <- str_replace_all(x, "^\\s*c\\(|\\)\\s*$", "")
  vals <- str_split(x, ",", simplify = FALSE)[[1]]
  vals <- str_trim(str_replace_all(vals, "^['\"]|['\"]$", ""))
  vals[vals != ""]
}

parse_coef_values <- function(x) {
  suppressWarnings(as.numeric(parse_vector_text(x)))
}

parse_class_values <- function(x) {
  parse_vector_text(x)
}

class_interest_score <- function(label) {
  normalized <- str_to_lower(as.character(label)) %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish()
  negative_pattern <- paste(
    c("uninfected", "untreated", "unexposed", "noninfected", "non infected",
      "not infected", "not treated", "not exposed", "control", "mock", "placebo",
      "healthy", "negative", "none", "absent", "baseline", "wild type",
      "wildtype", "susceptible", "not given", "notgiven", "no"),
    collapse="|"
  )
  positive_pattern <- paste(
    c("infected", "infection", "treated", "treatment", "exposed", "case", "diseased", "positive",
      "present", "resistant", "mutant", "yes"),
    collapse="|"
  )
  if (str_detect(normalized, paste0("(^| )(", negative_pattern, ")($| )"))) {
    return(-1)
  }
  if (str_detect(normalized, paste0("(^| )(", positive_pattern, ")($| )"))) {
    return(1)
  }
  0
}

choose_interesting_binary_class <- function(classes, metadata_values=character()) {
  classes <- unique(as.character(classes))
  classes <- classes[!is.na(classes) & nchar(classes) > 0]
  if (length(classes) != 2) {
    return(classes)
  }

  scores <- vapply(classes, class_interest_score, numeric(1))
  if (sum(scores == max(scores)) == 1) {
    return(classes[which.max(scores)])
  }

  counts <- table(factor(as.character(metadata_values), levels=classes))
  if (sum(counts == min(counts)) == 1) {
    return(classes[which.min(counts)])
  }

  sort(classes)[2]
}

reduce_binary_class_specs <- function(specs, metadata, confusion_log=NULL) {
  if (nrow(specs) == 0) {
    return(specs)
  }

  keys <- specs %>% distinct(metadata_category, matrix, kind)
  retained <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    key <- keys[i,]
    current <- specs %>%
      filter(metadata_category == key$metadata_category,
             matrix == key$matrix,
             kind == key$kind)
    classes <- unique(as.character(current$focus_class))
    classes <- classes[!is.na(classes) & nchar(classes) > 0]

    if (key$kind == "classification" && length(classes) == 2) {
      metadata_values <- character()
      if (!is.null(confusion_log) &&
          all(c("row_type", "metadata_category", "true_label") %in% colnames(confusion_log))) {
        metadata_values <- confusion_log %>%
          filter(row_type == "sample", metadata_category == key$metadata_category) %>%
          pull(true_label)
      } else if (key$metadata_category %in% colnames(metadata)) {
        metadata_values <- metadata[[key$metadata_category]]
      }
      selected_class <- choose_interesting_binary_class(classes, metadata_values)
      message(paste0("Binary target ", key$metadata_category,
                     ": generating heatmaps only for class '", selected_class, "'."))
      current <- current %>% filter(focus_class == selected_class)
    }
    retained[[i]] <- current
  }
  bind_rows(retained)
}

safe_file_label <- function(x) {
  x <- as.character(x)
  x <- ifelse(is.na(x) | nchar(x) == 0, "all", x)
  x <- str_replace_all(x, "[^A-Za-z0-9._-]+", "_")
  str_replace_all(x, "_+", "_")
}

get_confusion_log_path <- function(nonzero_path, metadata_path) {
  base <- nonzero_path
  base <- str_replace(base, "_nonzero_coefficients_blastp_annotated.*\\.tsv$", "_confusion_matrices")
  base <- str_replace(base, "_nonzero_coefficients_blast_annotated.*\\.tsv$", "_confusion_matrices")
  base <- str_replace(base, "_nonzero_coefficients_annotated.*\\.tsv$", "_confusion_matrices")
  ext <- tolower(tools::file_ext(metadata_path))
  if (ext == "csv") {
    paste0(base, ".csv")
  } else {
    paste0(base, ".tsv")
  }
}

read_confusion_log <- function(nonzero_path, metadata_path) {
  log_path <- get_confusion_log_path(nonzero_path, metadata_path)
  if (!file.exists(log_path)) {
    message("No confusion/scatter log found at ", log_path, "; falling back to coefficient-derived heatmaps.")
    return(NULL)
  }
  sep <- ifelse(tolower(tools::file_ext(log_path)) == "csv", ",", "\t")
  fread(log_path, sep = sep, fill = TRUE, colClasses = "character")
}

make_coef_long <- function(dt) {
  if (!"coefficients" %in% colnames(dt) && "first_coef" %in% colnames(dt)) {
    rows <- list()
    if (!"first_class" %in% colnames(dt)) {
      dt$first_class <- vapply(dt$classes, function(x) {
        vals <- parse_class_values(x)
        if (length(vals) == 0) "__continuous_target__" else vals[[1]]
      }, character(1))
    }
    rows[[length(rows) + 1]] <- dt %>%
      transmute(metadata_category,
                accuracy = if ("accuracy" %in% colnames(dt)) accuracy else NA_character_,
                feature, cluster,
                class_label = as.character(first_class),
                class_coef = suppressWarnings(as.numeric(first_coef)))
    if ("second_coef" %in% colnames(dt) && "second_class" %in% colnames(dt)) {
      rows[[length(rows) + 1]] <- dt %>%
        transmute(metadata_category,
                  accuracy = if ("accuracy" %in% colnames(dt)) accuracy else NA_character_,
                  feature, cluster,
                  class_label = as.character(second_class),
                  class_coef = suppressWarnings(as.numeric(second_coef)))
    }
    return(bind_rows(rows) %>%
             filter(!is.na(class_label), nchar(class_label) > 0, is.finite(class_coef)) %>%
             mutate(abs_coef = abs(class_coef)))
  }

  if (!"coefficients" %in% colnames(dt)) {
    stop("Coefficient table must contain either coefficients or first_coef columns.")
  }

  rows <- vector("list", nrow(dt))
  for (i in seq_len(nrow(dt))) {
    classes <- parse_class_values(dt$classes[[i]])
    coefs <- parse_coef_values(dt$coefficients[[i]])
    if (length(coefs) == 0 || all(is.na(coefs))) {
      next
    }
    if (length(classes) == 0) {
      classes <- rep("__continuous_target__", length(coefs))
    }
    n <- min(length(classes), length(coefs))
    rows[[i]] <- tibble(
      metadata_category = dt$metadata_category[[i]],
      accuracy = if ("accuracy" %in% colnames(dt)) dt$accuracy[[i]] else NA_character_,
      feature = dt$feature[[i]],
      cluster = dt$cluster[[i]],
      class_label = classes[seq_len(n)],
      class_coef = coefs[seq_len(n)]
    )
  }
  bind_rows(rows) %>%
    filter(is.finite(class_coef)) %>%
    mutate(abs_coef = abs(class_coef))
}

build_feature_labels <- function(dt, dt2) {
  labels1 <- tibble(feature = character(), cluster = character(), label = character())
  if ("stitle" %in% colnames(dt)) {
    labels1 <- dt %>%
      mutate(annotation = str_remove_all(as.character(stitle), "\\[.+\\]$|MULTISPECIES:\\s|, partial")) %>%
      mutate(label = clean_blast_label(annotation)) %>%
      select(feature, cluster, label)
  } else if ("label" %in% colnames(dt)) {
    labels1 <- dt %>% mutate(label = clean_blast_label(label)) %>% select(feature, cluster, label)
  }

  labels2 <- tibble(feature = character(), cluster = character(), label = character())
  if (!is.null(dt2) && "features" %in% colnames(dt2)) {
    labels2 <- dt2 %>%
      separate_longer_delim(features, delim = "},") %>%
      mutate(products = extract_feature_qualifier(features, "product"),
             genes = extract_feature_qualifier(features, "gene"),
             label = choose_feature_label(products, genes, opt$products)) %>%
      select(feature, cluster, label)
  } else if (!is.null(dt2) && "label" %in% colnames(dt2)) {
    labels2 <- dt2 %>% mutate(label = clean_blast_label(label)) %>% select(feature, cluster, label)
  }

  bind_rows(labels1, labels2) %>%
    filter(!is.na(label), nchar(label) > 0) %>%
    group_by(feature, cluster) %>%
    summarise(label = paste(unique(label), collapse = "; "), .groups = "drop")
}

make_panel_title <- function(panel_type, category, focus_class) {
  paste0(panel_type, "\n",
         "feature: ", category, "; class: ", focus_class)
}

draw_feature_heatmap <- function(matrix_data, annotation_df, annotation_colors, column_labels, title) {
  if (nrow(matrix_data) == 0 || ncol(matrix_data) == 0) {
    message("Skipping empty heatmap: ", title)
    return(invisible(NULL))
  }
  matrix_data[!is.finite(matrix_data)] <- 0
  if (length(annotation_colors) > 0) {
    ha <- rowAnnotation(df = annotation_df, col = annotation_colors)
  } else {
    ha <- rowAnnotation(df = annotation_df)
  }
  lim <- max(abs(matrix_data), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) {
    lim <- 1
  }
  heatmap_plot <- Heatmap(
    matrix_data,
    name = "Value",
    column_labels = column_labels[colnames(matrix_data)],
    left_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    heatmap_legend_param = list(
      title = expression("embedding" %*% ~ beta),
      at = c(-lim, 0, lim),
      labels = c("Negative", "Zero", "Positive")
    ),
    column_title = "Embedding features",
    row_title = "Samples"
  )
  draw(heatmap_plot, column_title = title, column_title_gp = grid::gpar(fontsize = 16),
       padding = unit(c(100, 2, 2, 2), "pt"))
}

draw_correlation_heatmap <- function(matrix_data, column_labels, title) {
  if (nrow(matrix_data) < 2 || ncol(matrix_data) < 2) {
    message("Skipping correlation heatmap with too few rows/columns: ", title)
    return(NULL)
  }
  matrix_data[!is.finite(matrix_data)] <- NA_real_
  cor_mat <- suppressWarnings(cor(matrix_data, use = "pairwise.complete.obs"))
  cor_mat[!is.finite(cor_mat)] <- NA_real_
  cor_heatmap <- Heatmap(
    cor_mat,
    name = "Correlation",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_labels = column_labels[colnames(cor_mat)],
    row_labels = column_labels[rownames(cor_mat)],
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (is.finite(cor_mat[i, j])) {
        grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 6))
      }
    },
    heatmap_legend_param = list(title = "Correlation", at = c(-1, 0, 1), labels = c("-1", "0", "1")),
    column_title = "Embedding features",
    row_title = "Embedding features"
  )
  draw(cor_heatmap, column_title = title, column_title_gp = grid::gpar(fontsize = 16),
       padding = unit(c(100, 2, 2, 2), "pt"))
  cor_mat
}

# read in input files
dt <- fread(opt$nonzero_annotations)
blast_annotated_path <- gsub("blastp_annotated", "blast_annotated", opt$nonzero_annotations)
dt2 <- if (file.exists(blast_annotated_path)) fread(blast_annotated_path) else NULL
feather_dt <- feather::read_feather(opt$feather_file)
all_metadata <- fread(opt$metadata)
confusion_log <- read_confusion_log(opt$nonzero_annotations, opt$metadata)

coef_long <- make_coef_long(dt)
feature_labels <- build_feature_labels(dt, dt2)

out_csvs_prefix <- file.path(dirname(opt$output), "raw_matrices", tools::file_path_sans_ext(basename(opt$output)))
dir.create(file.path(dirname(opt$output), "raw_matrices"), recursive = TRUE, showWarnings = FALSE)

pdf(opt$output, width=20, height=16)

# write a title page first
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, paramaters['dataset'])
text(5, 7, paramaters['filter'])
text(5, 6, paramaters['cluster_approach'])
text(5, 5, paramaters['model'])
text(5, 4, paste("At most", paramaters['num_clusters'], "clusters"))
text(5,3, paste(Sys.Date()))

if (!"sample_name" %in% colnames(feather_dt)) {
  stop("Feather matrix must contain a sample_name column.")
}

if (!is.null(confusion_log) &&
    all(c("row_type", "sample_name", "matrix") %in% colnames(confusion_log))) {
  categorical_specs <- confusion_log %>%
    filter(row_type == "sample") %>%
    mutate(true_label = as.character(true_label),
           predicted_label = as.character(predicted_label)) %>%
    group_by(metadata_category) %>%
    summarise(focus_class = list(sort(unique(c(true_label, predicted_label)))),
              .groups = "drop") %>%
    unnest(focus_class) %>%
    mutate(matrix = "all_samples", kind = "classification")
  specs <- categorical_specs
} else {
  specs <- coef_long %>%
    distinct(metadata_category, class_label) %>%
    filter(class_label != "__continuous_target__") %>%
    mutate(matrix = "all_samples",
           focus_class = class_label,
           kind = "classification") %>%
    select(metadata_category, matrix, focus_class, kind)
}

specs <- reduce_binary_class_specs(specs, all_metadata, confusion_log)

for (i in seq_len(nrow(specs))) {
  spec <- specs[i,]
  tryCatch({
    category <- spec$metadata_category
    matrix_name <- spec$matrix
    focus_class <- spec$focus_class
    kind <- spec$kind

    category_coef <- coef_long %>%
      filter(metadata_category == category)
    category_coef <- category_coef %>% filter(class_label == focus_class)

    important_features <- category_coef %>%
      arrange(desc(abs_coef)) %>%
      distinct(feature, .keep_all = TRUE) %>%
      head(max_features_per_heatmap) %>%
      select(feature, class_coef) %>%
      deframe()

    important_features <- important_features[names(important_features) %in% colnames(feather_dt)]
    if (length(important_features) == 0) {
      message("No finite coefficients/features for ", category, " ", matrix_name, " ", focus_class)
      next
    }
    if (!is.null(confusion_log) && "sample_name" %in% colnames(confusion_log)) {
        sample_rows <- confusion_log %>%
          filter(row_type == "sample",
                 metadata_category == category) %>%
          transmute(sample_name,
                    matrix = tolower(as.character(matrix)),
                    true_label = as.character(true_label),
                    predicted_label = as.character(predicted_label)) %>%
          distinct()
    } else {
      sample_rows <- all_metadata %>% select(sample_name)
        sample_rows <- all_metadata %>%
          select(sample_name, all_of(category)) %>%
          dplyr::rename(true_label = all_of(category)) %>%
          mutate(matrix = "all_samples",
                 predicted_label = NA_character_)
    }

    sub_feather_unscaled <- feather_dt %>%
      select(sample_name, all_of(names(important_features))) %>%
      inner_join(sample_rows, by = "sample_name")

    if (nrow(sub_feather_unscaled) == 0) {
      message("No samples for ", category, " ", matrix_name, " ", focus_class)
      next
    }

    sample_order <- sub_feather_unscaled$sample_name
    unscaled_heatmap_data <- sub_feather_unscaled %>%
      select(sample_name, all_of(names(important_features))) %>%
      column_to_rownames("sample_name") %>%
      as.matrix()
    scaled_heatmap_data <- sweep(unscaled_heatmap_data, 2, important_features[colnames(unscaled_heatmap_data)], `*`)

    label_map <- category_coef %>%
      filter(feature %in% names(important_features)) %>%
      arrange(match(feature, names(important_features))) %>%
      distinct(feature, cluster) %>%
      left_join(feature_labels, by = c("feature", "cluster")) %>%
      mutate(label = ifelse(is.na(label) | nchar(label) == 0, "NO MATCH", label),
             label = str_wrap(str_trunc(gsub(",", ", ", label), 80, side = "right"), 40)) %>%
      select(feature, label) %>%
      deframe()
    missing_labels <- setdiff(colnames(scaled_heatmap_data), names(label_map))
    if (length(missing_labels) > 0) {
      label_map[missing_labels] <- missing_labels
    }
    label_map <- label_map[colnames(scaled_heatmap_data)]
    ann <- sample_rows %>%
        filter(sample_name %in% sample_order) %>%
        arrange(match(sample_name, sample_order)) %>%
        select(sample_name, matrix, true_label, predicted_label) %>%
        column_to_rownames("sample_name")
      all_labels <- sort(unique(na.omit(c(ann$true_label, ann$predicted_label))))
      label_cols <- setNames(colorRampPalette(brewer.pal(max(3, min(8, length(all_labels))), "Dark2"))(length(all_labels)), all_labels)
      matrix_levels <- sort(unique(na.omit(ann$matrix)))
      matrix_cols <- c("train" = "#4D4D4D", "test" = "#B2182B", "all_samples" = "#4D4D4D")
      missing_matrix_levels <- setdiff(matrix_levels, names(matrix_cols))
      if (length(missing_matrix_levels) > 0) {
        matrix_cols <- c(matrix_cols,
                         setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(missing_matrix_levels)),
                                  missing_matrix_levels))
      }
      matrix_cols <- matrix_cols[names(matrix_cols) %in% matrix_levels]
      annotation_colors <- list(matrix = matrix_cols, true_label = label_cols, predicted_label = label_cols)

    draw_feature_heatmap(scaled_heatmap_data, ann, annotation_colors, label_map,
                         make_panel_title("sample by feature heatmap", category, focus_class))
    cor_mat <- draw_correlation_heatmap(scaled_heatmap_data, label_map,
                                        make_panel_title("feature correlation", category, focus_class))

    clustered_order <- rownames(scaled_heatmap_data)
    if (nrow(scaled_heatmap_data) > 1) {
      clustering_matrix <- scaled_heatmap_data
      clustering_matrix[!is.finite(clustering_matrix)] <- 0
      clustered_order <- rownames(scaled_heatmap_data)[hclust(dist(clustering_matrix))$order]
    }

    file_label <- paste(safe_file_label(category), safe_file_label(matrix_name), safe_file_label(focus_class), sep = "_")
    write_csv(
      scaled_heatmap_data[clustered_order, , drop = FALSE] %>%
        as.data.frame() %>% rownames_to_column("sample_name") %>%
        left_join(ann %>% rownames_to_column("sample_name"), by = "sample_name"),
      paste(out_csvs_prefix, file_label, "sample_by_feature_clustered_scaled_by_beta.csv", sep = "_"),
      col_names = TRUE, quote = "needed"
    )
    write_csv(
      unscaled_heatmap_data[clustered_order, , drop = FALSE] %>%
        as.data.frame() %>% rownames_to_column("sample_name") %>%
        left_join(ann %>% rownames_to_column("sample_name"), by = "sample_name"),
      paste(out_csvs_prefix, file_label, "sample_by_feature_clustered_unscaled.csv", sep = "_"),
      col_names = TRUE, quote = "needed"
    )
    if (!is.null(cor_mat)) {
      write_csv(
        cor_mat %>% as.data.frame() %>% rownames_to_column("feature"),
        paste(out_csvs_prefix, file_label, "feature_correlation_matrix.csv", sep = "_"),
        col_names = TRUE, quote = "needed"
      )
    }
  }, error = function(e) {
    message(paste("Error processing heatmap panel:", paste(spec, collapse = " "), "\n", e$message))
  })
}
dev.off()
