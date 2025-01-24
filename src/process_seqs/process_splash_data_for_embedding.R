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
  make_option(c("-a", "--anchor_file"), "File with list of anchors", type="character"),
  make_option(c("-c", "--cluster_file"), "File with list of anchor clusters", type="character"),
  make_option(c("i", "--id_mapping"), "File with sample id mapping", type="character"),
  make_option(c("-s", "--satc_files"), "Path to SPLASH SATC file directory", type="character"),
  make_option(c("-o", "--output_prefix"), "Output prefix.", type="character"),
  make_option(c("--temp_dir"), "Temporary directory to store intermediate files", 
              type="character"),
  make_option(c("--satc_util_bin"), "Path to SATC Util binary folder",
              type="character", default="/oak/stanford/groups/horence/dcotter1/satc_utils/"),
  make_option(c("--num_cores"), "Number of cores to use", type="integer", default = 8)
)


# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# set up parallel processing
plan(multicore, workers = opt$num_cores)

# check that user specified all files 
if (!file.exists(opt$anchor_file) | !file.exists(opt$cluster_file) | is.null(opt$satc_files) | is.null(opt$output_prefix) | !file.exists(opt$id_mapping)) {
  stop("Must provide anchor file, cluster file, satc files, id mapping file, and output prefix")
}

# create a temporary directory to store intermediate files
if (!is.null(opt$temp_dir)) {
  temp_dir <- ifelse(grepl("/$", opt$temp_dir),
                     opt$temp_dir, 
                     paste0(opt$temp_dir, "/"))
  system(paste("mkdir -p", temp_dir))
} else {
  temp_dir <- file.path(dirname(opt$output_prefix), "tmp/")
  system(paste("mkdir -p", temp_dir))
}

# read in the anchor cluster file
anchor_clusters <- fread(opt$cluster_file, 
                    header = F, col.names = c("cluster_id", "anchor")) %>%
  group_by(cluster_id) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# list all of the .satc files in the result_satc folder
satc_files <- list.files(opt$satc_files, pattern = ".satc", full.names = T)
satc_files <- data.frame(satc_file=satc_files) %>% 
  mutate(satc_dump = gsub(".satc", ".satc.dump", 
                          file.path(opt$temp_dir, 'dumped', basename(satc_file))))
system(paste("mkdir -p", file.path(opt$temp_dir, "dumped")))

# declare a satc file for the output of all the dump files
all_satc_file <- file.path(opt$temp_dir, "all_satc_merged.txt")

if (!file.exists(all_satc_file)) {
  # dump the satc files
  future_walk2(satc_files$satc_file, satc_files$satc_dump, \(x,y) system(
    paste0(file.path(opt$satc_util_bin, "satc_dump"), " --anchor_list ", opt$anchor_file,
           " --sample_names ", opt$id_mapping, " ", x, " ", y)))
  
  # merge them into one file
  walk(satc_files$satc_dump, \(x) system(paste("cat", x, ">>", all_satc_file)))

  # remove any lines that start or end in [ACTG]
  system(paste("grep -v '^[ACTG]{3}' ", all_satc_file, " | grep -v '[ACTG]$' > ", 
               file.path(opt$temp_dir, "all_satc_merged_no_anchor.txt")))
  system(paste("mv", file.path(opt$temp_dir, "all_satc_merged_no_anchor.txt"), all_satc_file))
}

# undump the satc files
all_satc_undumped <- file.path(opt$temp_dir, "all_satc_merged.undumped.txt")
all_satc_temp_mapping <- file.path(opt$temp_dir, "all_satc_merged.temp_mapping.txt")
if (!file.exists(all_satc_undumped) | !file.exists(all_satc_temp_mapping)) {
  system(paste(file.path(opt$satc_util_bin, "satc_undump"), "-i", all_satc_file, 
               "-o", all_satc_undumped, "-m", all_satc_temp_mapping))
}

# filter the satc files 
all_satc_filtered <- file.path(opt$temp_dir, "all_satc.filtered.txt")
if (!file.exists(all_satc_filtered)) {
  system(paste(file.path(opt$satc_util_bin, "satc_filter"), 
               "-i", all_satc_undumped, "-o", all_satc_filtered,
               "-d", opt$anchor_file, "-n", 1))
}

# redump the satc file 
all_satc_filtered_dump <- file.path(opt$temp_dir, "all_satc.filtered.dump")
if (!file.exists(all_satc_filtered_dump)) {
  system(paste(file.path(opt$satc_util_bin, "satc_dump"),
               "--sample_names", all_satc_temp_mapping,
               all_satc_filtered, all_satc_filtered_dump))
}

# read in the dumped satc file
satc_dt <- fread(all_satc_filtered_dump, header=F,
                 col.names=c("sample", "anchor", "target", "count"))

# grab the top anchor per cluster as a representative anchor
representative_anchors <- anchor_clusters %>% group_by(cluster_id) %>% 
  distinct(cluster_id, .keep_all=T) %>% ungroup() %>% pull(anchor)

# read in the satc and pivot it wider 
wide_satc <- merge(satc_dt, anchor_clusters, by="anchor", all.x=TRUE)

# identify clusters that are not in the satc file
missing_clusters <- setdiff(anchor_clusters$cluster_id, wide_satc$cluster_id)

# add the missing cluster for one sample with the representative anchor and Ns
if (!is_empty(missing_clusters)) {
  missing_df <- data.frame(sample=unique(head(wide_satc$sample))[1],
                           anchor=representative_anchors[missing_clusters],
                           target=strrep("N", nchar(representative_anchors[missing_clusters])),
                           count=0,
                           cluster_id=missing_clusters,
                           rank=1)
  
  wide_satc <- rbind(wide_satc, missing_df)
}

wide_satc <- as.data.table(wide_satc)
wide_satc <- wide_satc[order(cluster_id, rank)]
wide_satc <- unique(wide_satc, by=c("sample", "cluster_id"))

wide_satc <- wide_satc %>% mutate(seq=str_c(anchor, target, sep="")) %>% 
  select(sample, cluster_id, seq)

wide_satc <- dcast(wide_satc, sample ~ cluster_id, value.var="seq")
wide_satc <- as.data.frame(wide_satc)

# add the representative anchors to the wide satc
wide_satc <- cbind(wide_satc[1],
                   map2_df(wide_satc[,2:ncol(wide_satc)], 
                           1:length(representative_anchors), 
                           \(x,y) ifelse(is.na(x), str_c(representative_anchors[y], strrep("N", 27), sep = ""), x)))
wide_satc <- wide_satc %>% ungroup()

# join the columns together into one sequence and write to a tsv
seqs <- wide_satc %>% unite(seq, -sample, sep="")
seqs %>% write_tsv(paste0(opt$output_prefix, "_sample_sequences.tsv"))

# write the data to a fasta file
dna <- Biostrings::DNAStringSet(seqs$seq)
names(dna) <- seqs$sample
writeXStringSet(dna, paste0(opt$output_prefix, "_sample_sequences.fasta"))

