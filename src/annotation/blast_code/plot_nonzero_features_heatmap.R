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
  x <- replace_na(x, "")
  x <- str_replace_all(x, "LOC\\d+[- ]*", "")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_replace_all(x, "\\s*[,;]\\s*$", "")
  x <- str_trim(x)
  ifelse(nchar(x) == 0, NA_character_, x)
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

# read in input files
dt <- fread(opt$nonzero_annotations)
if (TRUE) {dt2 <- fread(gsub("blastp_annotated", "blast_annotated", opt$nonzero_annotations))}
all_clusters <- fread(opt$clusters) %>% select(-kmer)
feather_dt <- feather::read_feather(opt$feather)
all_metadata <- fread(opt$metadata)

categories <- dt %>% select(metadata_category, accuracy) %>% distinct() %>% arrange(-accuracy) %>% pull(metadata_category)

out_csvs_prefix <- file.path(dirname(opt$output), "raw_matrices", tools::file_path_sans_ext(basename(opt$output)))
system(paste("mkdir -p", file.path(dirname(opt$output), "raw_matrices")))

pdf(opt$output, width=20, height=16)

# write a title page first
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, paramaters['dataset'])
text(5, 7, paramaters['filter'])
text(5, 6, paramaters['cluster_approach'])
text(5, 5, paramaters['model'])
text(5, 4, paste("At most", paramaters['num_clusters'], "clusters"))
text(5,3, paste(Sys.Date()))

for (category in categories) {
  tryCatch({
    new_dt <- dt %>% filter(is.na(query)) %>% select(-query) %>% left_join(all_clusters %>% mutate(query = paste0(cluster, "_", seq)) %>% select(-seq), by="cluster", relationship="many-to-many") %>% 
      filter(!is.na(query))
    summ_dt <- bind_rows(dt, new_dt) %>% filter(!is.na(query)) %>% filter(metadata_category==category) %>%
      mutate(first_coef=get_first_coef(coefficients)) %>% mutate(max_coefficient=abs(first_coef)) %>% 
      arrange(-max_coefficient) %>% mutate(first_class=get_first_class(classes)) %>% mutate(annotation = str_remove_all(stitle, "\\[.+\\]$|MULTISPECIES:\\s|, partial")) %>%
      rowwise() %>%
      mutate(classes=list(str_split_1(gsub("\\[|\\]", "", classes),pattern=","))) %>%
      ungroup() %>%
      select(metadata_category, accuracy, classes, first_class, first_coef, max_coefficient, cluster, feature, query, identity, qcovs, annotation) %>%
      mutate(query = str_remove(query, "cluster_\\d+_")) %>%
      distinct(cluster,annotation,query,.keep_all = T) %>% group_by(cluster)  
    
    
    summ_dt <- summ_dt %>% group_by(cluster,query) %>%
      mutate(label=ifelse(!is_empty(unique(na.omit(annotation))), paste(unique(na.omit(annotation)),collapse=";"), NA)) %>%
      distinct(cluster, query, label, .keep_all=T) %>% ungroup()
    
    if (TRUE) {
      summ_dt2 <- dt2 %>% filter(metadata_category==category) %>%
        separate_longer_delim(features, delim = "},") %>% 
        mutate(products=extract_feature_qualifier(features, "product")) %>% 
        mutate(genes=extract_feature_qualifier(features, "gene")) %>% 
        select(-features) %>% mutate(first_coef=get_first_coef(coefficients)) %>% mutate(max_coefficient=abs(first_coef)) %>% 
        arrange(-max_coefficient) %>% mutate(first_class=get_first_class(classes)) %>%
        rowwise() %>%
        mutate(classes=list(str_split_1(gsub("\\[|\\]", "", classes),pattern=","))) %>%
        ungroup() %>%
        select(metadata_category, accuracy, classes, first_class, first_coef, max_coefficient, cluster, feature, query, identity, products, genes) %>%
        mutate(query = str_remove(query, "cluster_\\d+_")) %>%
        group_by(cluster) %>%
        ungroup() %>%
        distinct(cluster,products,query,genes,.keep_all = T) %>% group_by(cluster) 
      
      summ_dt2 <- summ_dt2 %>%
        mutate(label = choose_feature_label(products, genes, opt$products)) %>%
        group_by(cluster,query) %>%
        mutate(label=ifelse(!is_empty(unique(na.omit(label))), paste(unique(na.omit(label)),collapse=";"), NA)) %>%
        distinct(cluster, query, label, .keep_all=T) %>% ungroup()
      if (!"qcovs" %in% colnames(summ_dt2)) {
        summ_dt2$qcovs <- NA
      }
      summ_dt2 <- summ_dt2 %>% select(cluster, query, identity, qcovs, label) %>% dplyr::rename(label2=label)
      summ_dt <- summ_dt %>% left_join(summ_dt2 %>% 
                                         select(cluster, query, identity, qcovs, label2), by=c("cluster", "query")) %>% 
        dplyr::rename(identity=`identity.x`, qcovs=`qcovs.x`) %>%
        mutate(identity = ifelse(is.na(label) & !is.na(label2), `identity.y`, identity)) %>%
        mutate(qcovs = ifelse(is.na(label) & !is.na(label2), `qcovs.y`, qcovs)) %>%
        mutate(label = ifelse(is.na(label) & !is.na(label2), label2, label)) %>%
        mutate(label = ifelse(is.na(label) | nchar(label)<2, NA, label))
    }
    
    summ_dt <- summ_dt %>% group_by(feature) %>% mutate(label = ifelse(rep(sum(!is.na(label))==0, length(label)) & (is.na(label)) & (!is.na(identity) | !is.na(identity.y)), "NO PROTEIN/GENE HIT", label))
    
    summ_dt <- summ_dt %>% group_by(cluster) %>% 
      mutate(label=ifelse(label=="",annotation,label)) %>%
      mutate(label = ifelse(rep(sum(!is.na(label))==0, length(label)), "NO MATCH", label)) %>% 
      mutate(hypothetical=length(unique(label))>1 & sum(str_detect(label, "(?i)hypothetical|uncharacterized"))>0) %>%
      mutate(hypothetical=replace_na(hypothetical, FALSE)) %>% 
      rowwise() %>%
      mutate(label = map2_vec(label, hypothetical, \(x,y) if (y) {str_c(str_trim(unlist(str_split(x, ";"))[str_detect(unlist(str_split(x, ";")), "(?i)hypothetical|uncharact", negate=T)]),sep = ",", collapse=",")} else {x})) %>%
      ungroup()
    
    my_classes <- summ_dt[1,]$classes %>% unlist()
    
    if (sum(str_detect(my_classes, "0.0")) > 0) {
      my_classes <- as.numeric(my_classes)
    }
    
    summ_dt <- summ_dt %>% group_by(cluster, feature, max_coefficient, first_coef) %>% summarise(label=paste(unique(na.omit(label)), collapse=","), label2=paste(unique(na.omit(label2)), collapse=",")) %>% 
      arrange(-max_coefficient) %>% ungroup()
    
    important_features <- summ_dt %>% select(feature, first_coef) %>% head(max_features_per_heatmap) %>% deframe()
    
    sub_feather_dt <- feather_dt %>% select(sample_name, all_of(names(important_features)))
    
    sub_feather_unscaled <- sub_feather_dt
    
    sub_feather_dt <- sub_feather_dt %>% mutate(across(all_of(names(important_features)), \(x) x * important_features[cur_column()]))
    
    sub_feather_dt <- sub_feather_dt %>% left_join(all_metadata %>% select(sample_name, !!category) %>% dplyr::rename(class:=!!category)) %>% 
      relocate(class, .after=sample_name) 
    
    sub_feather_unscaled <- sub_feather_unscaled %>% left_join(all_metadata %>% select(sample_name, !!category) %>% dplyr::rename(class:=!!category)) %>% 
      relocate(class, .after=sample_name) 
    
    column_labels <- summ_dt %>% ungroup() %>% select(feature,cluster,label) %>%
      filter(feature %in% names(important_features)) %>%
      # mutate(label=str_c(cluster, label, sep=" ")) %>% 
      select(-cluster) %>%
      mutate(label=str_wrap(str_trunc(gsub(",",", ", label), 80, side="right"), 40)) %>% deframe()
    
    # Reshape the data to wide format for heatmap
    heatmap_data <- sub_feather_dt %>%
      filter(class %in% my_classes) %>%
      select(-class) %>%
      column_to_rownames("sample_name") %>%
      as.matrix()
    
    unscaled_heatmap_data <- sub_feather_unscaled %>%
      filter(class %in% my_classes) %>%
      select(-class) %>%
      column_to_rownames("sample_name") %>%
      as.matrix()
    
    # define dynamic colors
    n_colors <- min(length(my_classes), 8)  # Set2 has max 8 colors
    class_palette <- colorRampPalette(brewer.pal(max(n_colors,3), "Dark2"))(length(my_classes))
    class_colors <- setNames(class_palette, my_classes)
    
    
    # Create class annotation 
    class_annotation <- sub_feather_dt %>%
      filter(class %in% my_classes) %>%
      select(sample_name, class) %>%
      distinct() %>%
      column_to_rownames("sample_name")
    
    # Create a heatmap annotation for classes
    ha <- rowAnnotation(df = class_annotation, 
                        col = list(class = class_colors))
    
    # Create the heatmap with hierarchical clustering
    heatmap_plot <- Heatmap(heatmap_data, 
                            name = "Value",
                            column_labels = column_labels[colnames(heatmap_data)],
                            left_annotation = ha,
                            cluster_rows = TRUE,  # Enable hierarchical clustering for rows
                            cluster_columns = FALSE,  # Enable hierarchical clustering for columns
                            show_row_names = FALSE, 
                            show_column_names = TRUE,
                            heatmap_legend_param = list(title = expression("embedding" %*% ~ beta), at = c(min(heatmap_data), 0, max(heatmap_data)), labels = c("Negative", "Zero", "Positive")),
                            column_title = "Embedding Features",
                            row_title = "Samples",
                            
    )
    
    
    draw(heatmap_plot,column_title=category, column_title_gp=grid::gpar(fontsize=16),padding=unit(c(100,2,2,2), "pt"))
    
    cov_mat <- cov(heatmap_data)
    
    cov_heatmap <- Heatmap(cov_mat, 
                           name = "Covariance",
                           cluster_rows = TRUE,  # Enable hierarchical clustering for rows
                           cluster_columns = TRUE,  # Enable hierarchical clustering for columns
                           show_row_names = TRUE, 
                           show_column_names = TRUE,
                           column_labels = column_labels[colnames(cov_mat)],
                           row_labels = column_labels[rownames(cov_mat)],
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.text(sprintf("%.2f", cov_mat[i, j]), x, y, gp = gpar(fontsize = 6))},
                           heatmap_legend_param = list(title = "Covariance", at = c(min(cov_mat), 0, max(cov_mat)), labels = c("Negative", "Zero", "Positive")),
                           column_title = "Embedding Features",
                           row_title = "Embedding Features",
                           
    )
    
    draw(cov_heatmap)
    
    out_csv_beta_name <- paste(out_csvs_prefix, str_replace_all(category, "/", "_"), "nonzero_feature_matrix_scaled_by_beta.csv", sep = "_")
    write_csv(heatmap_data %>% as.data.frame() %>% rownames_to_column("sample_name"), out_csv_beta_name, col_names = T, quote="needed")
    
    out_csv_no_beta_name <- paste(out_csvs_prefix, str_replace_all(category, "/", "_"), "nonzero_feature_matrix_unsacled.csv", sep="_")
    write_csv(unscaled_heatmap_data %>% as.data.frame() %>% rownames_to_column("sample_name"), out_csv_no_beta_name, col_names = T, quote="needed")
    
  }, error = function(e) {
    message(paste("Error processing category:", category, "\n", e$message))
    # Optionally log the error or take other actions
  })
  
}
dev.off()
