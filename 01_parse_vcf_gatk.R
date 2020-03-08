# parse mpileup vcfs 

library("tidyverse")
library("vcfR")
library("parallel")

#list of the vcf files
vcf_files <- list.files("data/gatk_vcf", pattern = ".vcf$", full.names = TRUE)

parse_vcf_gatk <- function(vcf_file){
  
  cat(vcf_file)
  
  # read in the raw vcf for one line 
  # (as tidy data frame)
  gt_df <- read.vcfR(file = vcf_file) %>% vcfR2tidy(single_frame = TRUE) 
  
  gt_df <- gt_df$dat %>%
    select(-ID, -FILTER)
  
  # apply GATK best pratices filters for SNPs (biallelic)
  snps <- gt_df %>%
    filter(!grepl("[A-Z]{2,}",ALT)) %>% 
    filter(QD > 2.0, FS < 60.0, MQ > 40.0, MQRankSum > -12.5)
  
  # apply GATK best pratices filters for INDELs (biallelic)
  indels <- gt_df %>% 
    filter(grepl("[A-Z\\*]{2,}",ALT)) %>% 
    filter(QD > 2.0, FS < 200.0)
  
  # bind snps and indels back into a single frame
  gt_df <- bind_rows(snps, indels) %>%
    select(CHROM, POS, REF, ALT, Indiv, gt_AD, gt_DP, gt_GQ, gt_GT, gt_PGT, gt_GT_alleles)
  
  # remove data from unassigned contigs
  # only retain variable sites
  gt_df <- gt_df %>%
    filter(!(grepl("Unk", CHROM)))%>%
    filter(!is.na(ALT)) %>%
    filter(!is.na(gt_GT))
  
  # parse the indiv name into a sane sample name
  # and output
  gt_df  %>%
    mutate(plate = gsub("_[A-H]{1}[0-9]{2}$", "", Indiv)) %>%
    mutate(well = gsub(".*_", "", Indiv)) %>%
    select(CHROM, POS, Indiv, plate, well, REF, ALT, gt_AD, gt_DP, gt_GQ, gt_GT, gt_PGT, gt_GT_alleles)
  
}

#gt_df_list <- lapply(vcf_files, parse_vcf_gatk)
gt_df_list <- mclapply(vcf_files, parse_vcf_gatk, mc.cores = 2)

gt_df_out <- bind_rows(gt_df_list)

write_rds(gt_df_out, "data/gt_df/gt_seq_all_gt_df_gatk.rds", compress = "gz")
