# compute homology for recombination intervals
# using illumina and pacbio data
# KMS 2019

library("tidyverse")
library("parallel")
library("vcfR")

###############################################################
# Illumina 'all sites' data
###############################################################

# list the illumina files
il_files <- list.files("analysis_sv/data/homology/illumina", full.names = TRUE)
vcf_file <- il_files[1]

# get the individual codes for the sample IDs
meta_df <- read.table("meta/fasta_samplelist_recode.txt", h = F)
names(meta_df) <- c("Indiv", "ID")

parse_illumina_allsites <- function(vcf_file){
  
  # read in the allsites vcf for one region
  message(vcf_file)
  
  # tidy and subset the vcf
  vcf_df <- read.vcfR(vcf_file) %>%
    vcfR2tidy(single_frame = TRUE)
  
  vcf_df <- vcf_df$dat %>% 
    select(-ID) %>%
    left_join(meta_df) %>% 
    select(-Indiv) %>%
    rename(Indiv = ID) %>%
    select(CHROM, POS, Indiv, gt_GQ, gt_RGQ, gt_DP, gt_GT_alleles)
  
  # apply genotype filters
  # DP > 15, RGQ|GQ > 80
  
  vcf_df <- vcf_df %>%
    filter(gt_DP > 10, gt_RGQ > 40 | gt_GQ > 40) %>%
    select(-gt_GQ, -gt_RGQ, -gt_DP)
  
  # divide alleles into two columns
  vcf_df <- vcf_df %>%
    separate(gt_GT_alleles, into = c("gt1", "gt2"))
  
  # extract the tester line genotypes
  # (all samples will be compared to this)
  
  mv_df <- vcf_df %>%
    filter(Indiv == "MV2-25") %>%
    rename(mv_gt1 = gt1, mv_gt2 = gt2) %>%
    select(-Indiv)
  
  # compute number of differences
  # separate snps and indels
  vcf_df <- vcf_df %>%
    left_join(mv_df) %>%
    arrange(Indiv, CHROM, POS) %>%
    drop_na %>%
    filter(Indiv != "MV2-25")
  
  # count snp differences 
  vcf_snps <- vcf_df %>%
    filter(nchar(gt1) == 1 & nchar(gt2) == 1 & nchar(mv_gt1) == 1 & nchar(mv_gt2) == 1) %>%
    mutate(gt1_diff = gt1 != mv_gt1) %>%
    mutate(gt2_diff = gt2 != mv_gt2) %>%
    mutate(gt1_gc_loss = ifelse(gt1 %in% c("A", "T") & mv_gt1 %in% c("G", "C"), 1, 0)) %>%
    mutate(gt2_gc_loss = ifelse(gt2 %in% c("A", "T") & mv_gt2 %in% c("G", "C"), 1, 0)) %>%
    mutate(gt1_gc_gain = ifelse(gt1 %in% c("G", "C") & mv_gt1 %in% c("A", "T"), 1, 0)) %>%
    mutate(gt2_gc_gain = ifelse(gt2 %in% c("G", "C") & mv_gt2 %in% c("A", "T"), 1, 0)) %>%
    mutate(diff_sum = (gt1_diff + gt2_diff)/2) %>%
    mutate(diff_sum_gc_loss = (gt1_gc_loss + gt2_gc_loss)/2) %>%
    mutate(diff_sum_gc_gain = (gt1_gc_gain + gt2_gc_gain)/2) %>%
    mutate(n_sites = 1) %>%
    mutate(var_class = "snp") %>%
    select(-gt1_diff, -gt2_diff, -gt1_gc_loss, -gt2_gc_loss, -gt1_gc_gain, -gt2_gc_gain)
  
  # count indel differences
  vcf_indels <- vcf_df %>%
    filter(nchar(gt1) != 1 | nchar(gt2) != 1 | nchar(mv_gt1) != 1 | nchar(mv_gt2) != 1) %>%
    rowwise %>%
    mutate(gt1_diff = adist(gt1, mv_gt1)) %>%
    mutate(gt2_diff = adist(gt2, mv_gt2)) %>%
    mutate(diff_sum = (gt1_diff + gt2_diff)/2) %>%
    mutate(n_sites = max(nchar(gt1), nchar(gt2), nchar(mv_gt1), nchar(mv_gt2))) %>%
    mutate(diff_sum_gc_loss = NA, diff_sum_gc_gain = NA) %>%
    ungroup %>%
    mutate(var_class = "indel") %>%
    select(-gt1_diff, -gt2_diff)
  
  diff_df <- bind_rows(vcf_snps, vcf_indels) %>%
    select(-gt1, -gt2, -mv_gt1, -mv_gt2) %>%
    group_by(Indiv, var_class) %>%
    summarise(difference_count = sum(diff_sum), gc_gain = sum(diff_sum_gc_gain, na.rm = TRUE),
              gc_loss = sum(diff_sum_gc_loss, na.rm = TRUE), n_sites = sum(n_sites))
  
  # output differences to file
  
  interval_name <-  vcf_file %>%
    strsplit(split = "_") %>%
    unlist
  
  data.frame(CHROM = interval_name[1], pos1 = interval_name[2], pos2 = interval_name[3], 
             diff_df)
    
}

il_homology <- mclapply(il_files, parse_illumina_allsites, mc.cores = 4)

il_homology %>%
  bind_rows %>%
  saveRDS(file = "analysis_sv/data/homology/illumina_homology.rds")
