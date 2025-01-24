# calculate_kmer_dinucleotide_freqs.R
# Date: 2025-01-24
# Description: Given a fasta file of kmers, and an ordering file of kmer ids for each sample, calculate a matrix of 
# dinucleotide frequencies for eahc sample. The rows will be the samples and the columns will be dinculeotide frequencies
# for each cluster of kmers. The first n columns will be the dinucleotide frequencies for the first cluster of kmers, the
# next n columns will be the dinucleotide frequencies for the second cluster of kmers, and so on. The output will be a
# matrix of dinucleotide frequencies for each sample.

## import packages --------
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(furrr))

## parse arguments ------
# define command line arguments
option_list <- list(
  make_option(c("-o", "--ordering"), help="Ordering file", type="character"),
  make_option(c("-k", "--kmers"), help="Kmers fasta file", type="character"),
  make_option(c("-p", "--output"), help="Output file.", type="character", default=""),
  make_option(c("-t", "--threads"), help="Number of threads to use.", type="integer", default=1),
  make_option(c("--normalized"), help="Normalize the dinucleotide frequencies.", type="logical", default=FALSE)
)

# parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# check if oredering, kmers, and output files exist
if (!file.exists(opt$order) | !file.exists(opt$kmers)) {
  stop("Ordering or kmers file does not exist.")
}

if (opt$output == "") {
  stop("Output file not specified.")
}

# set up future
plan(multisession, workers = opt$threads)

## read in kmers and ordering data
# read in kmers
cat("Loading the kmers file...\n")
kmers <- readDNAStringSet(opt$kmers)

# read in ordering
# load the ordering file
cat("Loading the ordering file...\n")
ordering <- fread(opt$ordering, header=F, sep="\t", 
                  col.names = c("sample_name", "seq", "kmer", "start", "end")) 
ordering <- ordering %>% select(sample_name, kmer)

# labeling clusters in the ordering file (they appear one per row per sample)
ordering <- ordering %>% group_by(sample_name) %>% mutate(cluster=row_number()) %>%
  mutate(cluster=cluster-1) %>% ungroup()

## calculate dinucleotide frequencies for each kmer
if (opt$normalized) {
  dinucleotide_freqs <- dinculeotideFrequency(kmers, as.prob=TRUE)
} else {
  dinucleotide_freqs <- dinculeotideFrequency(kmers, as.prob=FALSE)
}
rownames(dinucleotide_freqs) <- names(kmers)
dinucleotide_freqs <- as.data.frame(dinucleotide_freqs) %>% rownames_to_column("kmer")

## merge the dinucleotide frequencies with the ordering file
merged_dt <- ordering %>% left_join(freqs, by="kmer") %>% select(-kmer)

# pivot wider and rename the columns so they are cluster_num_{dinucleotide_freq}
wide_dt <- merged_dt %>% group_by(sample_name) %>% mutate(cluster=paste0("cluster_", cluster)) %>%
  pivot_wider(names_from=cluster, values_from=colnames(dinucleotide_freqs)[-1], names_glue="{cluster}_{.value}")

cat("Writing per sample per cluster dinucleotide freqs to ", opt$output, "\n")
feather::write_feather(wide_dt, opt$output)