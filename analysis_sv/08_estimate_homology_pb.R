# compute homology for recombination intervals
# using illumina and pacbio data
# KMS 2019

library("tidyverse")
library("parallel")
library("seqinr")

###############################################################
# PacBio Consensus FASTAs
###############################################################
# list the pacbio fastas
pb_files <- list.files("analysis_sv/data/homology/pacbio", pattern = "fasta", full.names = TRUE)
fasta_file <- pb_files[1]

# extract the individual ID from the file name
Indiv <- pb_files %>% gsub(".*/", "", .) %>% 
  sub("_recom.*", "", .) %>%
  paste0("S", .)

file_df <- data.frame(Indiv, pb_files)


# get the individual codes for the sample IDs
meta_df <- read.table("meta/fasta_samplelist_recode.txt", h = F)
names(meta_df) <- c("Indiv", "ID")

# join to the file data
file_df <- left_join(meta_df, file_df)

# prepare file list excluding mv2-25
pb_files <- file_df %>%
  filter(ID != "MV2-25") %>%
  filter(!is.na(pb_files)) %>%
  pull(pb_files)

# prepare MV2-25 (precomputed, is the same for each)
fasta_file <- file_df %>%
  filter(ID == "MV2-25") %>%
  pull(pb_files)

tmp_name <- gsub(".fasta", "_tmp.fasta", fasta_file)

sys_call <- paste0("sed 's/[atcg]/-/g' ", fasta_file, "> ", tmp_name)
system(sys_call)

# get the fasta sequence
fasta_raw <- read.fasta(tmp_name, forceDNAtolower = FALSE)
file.remove(tmp_name)

# get the list of sequence names
# repair the "-roup" name
seq_names <- names(fasta_raw) %>%
  gsub("-rou", "grou", .)

names(fasta_raw) <- seq_names

# function for converting seqs to data frames
seq_to_df <- function(seq_name){
  
  CHROM <- seq_name %>% gsub(":.*", "", .)
  pos1 <- seq_name %>% gsub(".*:", "", .) %>% gsub("-.*", "", .) %>% as.numeric()
  pos2 <- seq_name %>% gsub(".*:", "", .) %>% gsub(".*-", "", .) %>% as.numeric()
  pos2 <- pos2 - 1
  POS <- pos1:(pos2)
   
  gt <- fasta_raw[seq_name][[1]] %>% as.character

  data.frame(seq_name, CHROM, POS, gt)
}

# build the mv genotype data frame
mv_gt_df <- lapply(seq_names, seq_to_df) %>%
  bind_rows

names(mv_gt_df)[4] <- "mv_gt"
rm(fasta_raw)
gc()

fasta_file <- pb_files[1]

# function for computing homology for each individual
# windows are interated within this function 
compute_homology_pac_bio <- function(fasta_file){
  
  message(fasta_file)
  
  # shell call to replace lowercase DNA with missing data
  # genomicconsensus adds in UNCALLED reference alleles as lowercase DNA
  # we really don't want these, at all
  tmp_name <- gsub(".fasta", "_tmp.fasta", fasta_file)
  
  sys_call <- paste0("sed 's/[atcg]/-/g' ", fasta_file, "> ", tmp_name)
  system(sys_call)
  
  # get the fasta sequence
  fasta_raw <- read.fasta(tmp_name, forceDNAtolower = FALSE)
  file.remove(tmp_name)

  # get the list of sequence names
  # repair the "-roup" name error
  seq_names <- names(fasta_raw) %>%
    gsub("-rou", "grou", .)
  
  names(fasta_raw) <- seq_names
  
  # build the mv genotype data frame
  gt_df <- lapply(seq_names, seq_to_df) %>%
    bind_rows
  
  gt_df <- cbind(gt_df, mv_gt_df$mv_gt)
  names(gt_df)[5] <- "mv_gt"
  
  # count differences + n_sites
  gt_df <- gt_df %>% 
    mutate(gt_diff = ifelse(gt != "-" & mv_gt != "-", as.numeric(gt == mv_gt), NA))
  
  # get the individual id 
  ID <- file_df %>%
    filter(pb_files == fasta_file) %>%
    pull(ID)
  
  gt_df %>%
    group_by(seq_name) %>%
    summarise(sum_diff = sum(gt_diff, na.rm = TRUE), n_sites = sum(!is.na(gt_diff))) %>%
    mutate(ID = ID) %>%
    select(ID, seq_name, sum_diff, n_sites)
  
}

pb_hmlgy <- lapply(pb_files, compute_homology_pac_bio)

#pb_hmlgy <- mclapply(pb_files, compute_homology_pac_bio, mc.cores = 4)
bind_rows(pb_hmlgy) %>%
  saveRDS("analysis_sv/data/homology/pb_homology.rds")
