# compute sequence divergence from reference for 
# aligned recombination genes
# KMS Jan 2019

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biostrings", version = "3.8")

options(device = "quartz")

library("tidyverse")
library("seqinr")

# the haplotype summaries from the previous script
# (to get relative position data)

hap_df <- readRDS("data/fasta_recomb/summaries/cds_alt_haps.rds")




# function for translating sequences
# and computing the number of difference vs. the reference CDS

translate_cds_seqs <- function(aln_cds_file){
  
  # want: 
  # gene_id sample_id hap_num codon_num amino_acid 
  
  # read in the fastas for all the samples
  aln_tmp <- read.alignment(aln_cds_file, format = "fasta")
  
  # translate all the sequences
  aa_df <- lapply(1:length(aln_tmp$nam), function(x) data.frame(info_string = aln_tmp$nam[x], amino_acid =
    translate(aln_tmp$seq[x] %>% 
                unlist %>% 
                strsplit(split = "") %>% 
                unlist))) %>%
      bind_rows
  
  return(aa_df)
  
}

# aligned CDS files
aln_cds_files <- list.files("data/fasta_recomb/cds_aligned", full.names = TRUE)

# apply the translation function
aa_df <- lapply(aln_cds_files, translate_cds_seqs) %>% 
  bind_rows

# parse the sample string for gene/sample/hap
aa_df <- aa_df %>%
  mutate(gene_name_full = gsub("_samp.*", "", info_string)) %>%
  mutate(sample_id = gsub(".*=|_.*", "", info_string)) %>%
  mutate(hap_num = gsub(".*_", "", info_string)) %>%
  select(gene_name_full, sample_id, hap_num, amino_acid)

# assign codon positions
aa_df <- aa_df %>%
  group_by(gene_name_full, sample_id, hap_num) %>%
  mutate(residue_num = 1:length(amino_acid)) %>%
  ungroup


# pull out CDS seq
cds_aa_df <- aa_df %>%
  filter(sample_id == "CDS", hap_num == "hap1") %>%
  rename(cds_amino_acid = amino_acid) %>%
  select(gene_name_full, residue_num, cds_amino_acid)

# join back into main data frame
# count differences at site level
aa_diff_df <- aa_df %>%
  left_join(cds_aa_df) %>%
  group_by(gene_name_full, sample_id, hap_num) %>%
  spread(key = hap_num, value = amino_acid) %>%
  mutate(aa_diff1 = as.numeric(hap1 != cds_amino_acid)) %>%
  mutate(aa_diff2 = as.numeric(hap2 != cds_amino_acid)) %>%
  mutate(aa_diff_sum = aa_diff1 + aa_diff2) %>%
  ungroup

# get synonymous changes
# decompose haplotypes into codons

split_codons_gene <- function(gene_name){
  
  message(gene_name)
  
  hap_sub <- hap_df %>%
    filter(gene_name_full == gene_name)
  
  # get sample_ids 
  sample_ids <- hap_sub$sample_id %>% unique
  
  # split codons 
  
  split_codon_sample <- function(samp_id){
    
    hap_sub_samp <- hap_sub %>%
      filter(sample_id == samp_id)
    
    # split codons
    codon_split <- split(hap_sub_samp$hap1, ceiling(seq_along(hap_sub_samp$hap1)/3))
    codon_split <- lapply(names(codon_split), function(x) data.frame(residue_num = as.numeric(x), 
                                                                     hap1 = codon_split[[x]])) %>% bind_rows
     
    data.frame(hap_sub_samp, residue_num = codon_split$residue_num)
  }
  
  codon_df <- lapply(sample_ids, split_codon_sample) %>%
    bind_rows
  
  codon_df
  
}

# compute residue positions
codon_df <- lapply(hap_df$gene_name_full %>% unique, split_codons_gene) %>%
  bind_rows

# build haplotype residues and translate
residue_df <- codon_df %>%
  group_by(gene_name_full, sample_id, residue_num) %>%
  summarise(hap1_codon = ifelse(length(hap1) == 3, paste0(hap1, collapse = ""), NA_character_), 
            hap2_codon = ifelse(length(hap2) == 3, paste0(hap2, collapse = ""), NA_character_),
            hap1_aa = ifelse(length(hap1) == 3, translate(hap1), NA_character_),
            hap2_aa = ifelse(length(hap1) == 3, translate(hap2), NA_character_)) %>%
  ungroup

# assign synonymous vs. nonsynt

# pull out CDS seq and join into residue_df
cds_residue_df <- residue_df %>%
  filter(sample_id == "CDS") %>%
  rename(cds_codon = hap1_codon, cds_aa = hap1_aa) %>%
  select(gene_name_full, residue_num, cds_codon, cds_aa)

residue_df <- residue_df %>%
  left_join(cds_residue_df)

residue_df <- residue_df %>%
  mutate(hap1_aa_diff = case_when(hap1_aa == cds_aa & hap1_codon == cds_codon ~ "identical",
                                  hap1_aa == cds_aa & hap1_codon != cds_codon ~ "synonymous",
                                  hap1_aa != cds_aa & hap1_codon != cds_codon ~ "non-synonymous",
                                  TRUE ~ NA_character_)) %>%
  mutate(hap2_aa_diff = case_when(hap2_aa == cds_aa & hap2_codon == cds_codon ~ "identical",
                                  hap2_aa == cds_aa & hap2_codon != cds_codon ~ "synonymous",
                                  hap2_aa != cds_aa & hap2_codon != cds_codon ~ "non-synonymous",
                                  TRUE ~ NA_character_))

# write out the codon/residue info
saveRDS(residue_df, "data/fasta_recomb/summaries/residue_df.rds")

# 
