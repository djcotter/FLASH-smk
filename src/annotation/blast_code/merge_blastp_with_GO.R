# Load necessary libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(furrr))
#library(Biostrings)

# gene refseq file: https://ftp.ncbi.nlm.nih.gov/refseq/uniprotkb/gene_refseq_uniprotkb_collab.gz
# go mapping file: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
# filtered with cut -f1,7 idmapping_selected.tab | awk '{if ($2) print $0;}' > id_mapping_selected_go_mapping_filtered.tab

# Define command line options
option_list <- list(
  make_option(c("-b", "--blast_folder"), type = "character", default = NULL, 
              help = "Path to the blastp output folder", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL, 
              help = "Path to the output tsv file", metavar = "character"),
  make_option(c("-u", "--uniprot_mapping_path"), type = "character", 
              help = "Path to the uniprot refseq mapping file", metavar = "character",
              default="/scratch/users/dcotter1/gene_refseq_uniprotkb_collab"),
  make_option(c("-g", "--go_mapping_path"), type = "character",  
              help = "Path to the gene2go mapping file", metavar = "character",
              default="/scratch/users/dcotter1/idmapping_selected_go_mapping_filtered.tab"),
  make_option(c("-t", "--max_workers"), type = "integer", default = 2, 
              help = "Number of threads to use", metavar = "integer")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$blast_folder) || is.null(opt$output_file)) {
  print_help(opt_parser)
  stop("All arguments must be supplied", call. = FALSE)
}

# set future plan
plan(multisession, workers=opt$max_workers)

# Read in the data
blast_files <- list.files(path=opt$blast_folder, pattern="blastout.tsv", full.names = T)

empty_blastp_df <- function() {
  tibble(query=character(), subject=character(), identity=numeric(),
         alignment_length=numeric(), mismatches=numeric(), gap_opens=numeric(),
         q_start=numeric(), q_end=numeric(), s_start=numeric(), s_end=numeric(),
         sstrand=character(), evalue=numeric(), qcovs=numeric(), qframe=character(),
         sgi=character(), sacc=character(), slen=numeric(), staxids=character(),
         stitle=character())
}

read_blastp_file <- function(file) {
  if (!file.exists(file) || file.info(file)$size == 0) {
    return(empty_blastp_df())
  }
  preview <- tryCatch(fread(file, sep="\t", nrow=5), error=function(e) NULL)
  if (is.null(preview) || ncol(preview) == 0) {
    return(empty_blastp_df())
  }
  if (ncol(preview) == 19) {
    return(fread(file, sep="\t", col.names = c("query", "subject", "identity", "alignment_length",
                                               "mismatches", "gap_opens", "q_start", "q_end",
                                               "s_start", "s_end", "sstrand", "evalue", "qcovs", "qframe",
                                               "sgi", "sacc", "slen", "staxids", "stitle")))
  }
  fread(file, sep="\t", col.names = c("query", "subject", "identity", "alignment_length",
                                      "mismatches", "gap_opens", "q_start", "q_end",
                                      "s_start", "s_end", "sstrand", "evalue", "qcovs",
                                      "sgi", "sacc", "slen", "staxids", "stitle")) %>%
    mutate(qframe="-")
}

blast_dfs <- future_map(blast_files, read_blastp_file)

# remove any data frames that had no data 
valid_blast_dfs <- map_vec(blast_dfs, \(x) nrow(x)>0)
blast_files <- blast_files[valid_blast_dfs]
blast_dfs <- blast_dfs[valid_blast_dfs]

df_lengths <- map(blast_dfs, \(x) nrow(x)) %>% unlist()

if (length(df_lengths) == 0) {
  empty_blastp_df() %>%
    mutate(NCBI_protein_accession=character(), UniProt_accession=character(),
           method=character(), GO=character()) %>%
    select(query, identity, evalue, qcovs, qframe, staxids, stitle,
           NCBI_protein_accession, UniProt_accession, method, GO) %>%
    write_tsv(opt$output_file, col_names = T, quote="needed")
  quit(save="no", status=0)
}

if (mean(as.numeric(df_lengths), na.rm=TRUE) > 350) {
  plan(multisession, workers=4)
}

merge_on_go_terms <- function(file, df, uniprot_mapping, go_mapping) {
  # grep the uniprot mapping file for presence of identical or similar genes
  temp_uniprot_file = str_replace(file, "out.tsv", "_uniprot_matches.tsv")
  system(paste0("cut -f15 ", file, " | sort | uniq | grep -Ff - ", uniprot_mapping, " | sort -k2 > ", temp_uniprot_file))
  
  # grep the idmapping file to match these uniprot ids to go terms
  temp_go_file = str_replace(file, "out.tsv", "_idmapping.tsv")
  system(paste0("cut -f2 ", temp_uniprot_file, " | sort | uniq | grep -Ff - ", go_mapping, " > ", temp_go_file))
  
  go_df <- fread(temp_go_file, col.names=c("UniProt_accession",
                                           "GO")) %>%
    select(UniProt_accession, GO)
  uniprot_df <- fread(temp_uniprot_file, col.names = c("NCBI_protein_accession", "UniProt_accession", 
                                                         "NCBI_tax_id",	"UniProtKB_tax_id",	"method")) %>%
    select(NCBI_protein_accession, UniProt_accession, method)
  
  df <- df %>% mutate(NCBI_protein_accession=str_extract(subject, "ref\\|(.+)\\|", group=1)) %>% 
    left_join(uniprot_df, by="NCBI_protein_accession", relationship = "many-to-many") %>% 
    left_join(go_df, by="UniProt_accession") %>% 
    mutate(staxids = as.character(staxids))
  return(df)
}

if (TRUE) {
  merged_df <- map(blast_dfs, \(x) x %>% mutate(staxids=as.character(staxids))) %>% bind_rows() %>% mutate(NCBI_protein_accession=str_extract(subject, "ref\\|(.+)\\|", group=1)) %>%
    mutate(UniProt_accession=str_extract(subject, "sp\\|(.+)\\|", group=1), method=NA, GO=NA) %>%
    select(query, identity, evalue, qcovs, qframe, staxids, stitle, NCBI_protein_accession, UniProt_accession, method, GO)
} else {
  merged_df <- future_map2(blast_files, blast_dfs, \(x,y) merge_on_go_terms(x,y,opt$uniprot_mapping_path,opt$go_mapping_path)) %>% bind_rows()
}

merged_df %>% select(query, identity, evalue, qcovs, qframe, staxids, stitle, NCBI_protein_accession, UniProt_accession, method, GO) %>%
  write_tsv(opt$output_file, col_names = T, quote="needed")
