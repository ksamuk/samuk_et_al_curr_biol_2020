# join all sv data types
# PacBio pbsv + Smoove SV + PopPoolation TE2
# KMS 2019

#install.packages("vcfR")
#install.packages("tidyverse")
library("vcfR")
library("tidyverse")

######################################
# smoove
######################################

# read in the smoove vcf
smoove_file <- list.files("analysis_sv/data/sv_raw", full.names = TRUE, pattern = "smoove")
vcf_df <- read.vcfR(smoove_file)
vcf_df <- vcfR2tidy(vcf_df, single_frame = TRUE)

# write smoove annotation descriptions to file
write.table(vcf_df$meta, "analysis_sv/meta/smoove_vcf_metadata.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

# retain main table remove unknown contigs
vcf_df <- vcf_df$dat %>%
  filter(!grepl("Unknown", CHROM)) 

# save a tidied df
saveRDS(vcf_df, "analysis_sv/data/sv_tidy/tidy_smoove.rds")

rm(vcf_df)


######################################
# pbsv
######################################

# read in the smoove vcf
pbsv_file <- list.files("analysis_sv/data/sv_raw", full.names = TRUE, pattern = "pacbio")
vcf_df <- read.vcfR(pbsv_file)
vcf_df <- vcfR2tidy(vcf_df, single_frame = TRUE)

# write smoove annotation descriptions to file
write.table(vcf_df$meta, "analysis_sv/meta/pbsv_vcf_metadata.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

# retain main table remove unknown contigs
vcf_df <- vcf_df$dat %>%
  filter(!grepl("Unknown", CHROM)) 

# save a tidied sm
saveRDS(vcf_df, "analysis_sv/data/sv_tidy/tidy_pbsv.rds")

rm(vcf_df)


######################################
# popte2
######################################

# read in the popte2 signatures file
isote_file <- list.files("analysis_sv/data/sv_raw", full.names = TRUE, pattern = "wider")
iso_te_df <- read_delim(isote_file, delim = "\t", 
                        col_types = "ncncccccn", col_names = c("Indiv", "CHROM", "pos1", "pos2", 
                                                               "family", "te_type", "direction", "strand", "freq"))

# # doesn't work -- popte2 doesn't actually allow recovery of individual TE info :(
# 
# # read in the te sizes
# te_lengths <- read.FASTA("repeat_masker/drosophila_tes_hill_et_al_2018.fasta")
# te_name <- names(te_lengths) %>% unlist
# 
# te_lengths <- data.frame(te_name) %>%
#   separate(te_name, sep = "\t", into = c("name", "name2", "family", "te_type", "source", "sv_length")) %>%
#   select(-name2) %>%
#   select(family, te_type, sv_length)
# 
# # join in te lengths
# iso_te_df %>%
#   left_join(te_lengths)

# remove unknown contigs
iso_te_df <- iso_te_df %>%
  filter(!grepl("Unknown", CHROM)) 


# save a tidied sm
saveRDS(iso_te_df, "analysis_sv/data/sv_tidy/tidy_popte2.rds")

rm(vcf_df)

######################################
# gatk indels
######################################

# read in the popte2 signatures file

create_indel_df <- function(gatk_file){
  
  vcf_df <- read.vcfR(gatk_file)
  vcf_df <- vcfR2tidy(vcf_df, single_frame = TRUE, format_fields = c("DP", "RGQ", "GQ", "GT"), info_fields = c("QD", "FS", "ReadPosRankSum"))
  
  # write smoove annotation descriptions to file
  #write.table(vcf_df$meta, "analysis_sv/meta/gatk_vcf_metadata.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # retain main table remove unknown contigs
  vcf_df <- vcf_df$dat %>%
    filter(!grepl("Unknown", CHROM)) 
  
  # GATK BP hard filters for indels
  # "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
  
  # only retain indels
  vcf_df <- vcf_df %>%
    filter(!is.na(ALT)) %>%
    filter(QD > 2.0, FS < 200, ReadPosRankSum > -20.0) %>% 
    filter(gt_DP > 20, gt_RGQ > 80 | gt_GQ > 80) %>%
    filter(!((nchar(ALT) == 1) & nchar(REF == 1))) %>% 
    filter(!(grepl(",", ALT))) %>%
    filter(!(grepl("[A-Z]{1},[A-Z]{1}", ALT))) %>%
    filter(!(grepl("[A-Z]{1},[A-Z]{1},[A-Z]{1}", ALT))) %>%
    filter(!(grepl("[A-Z]{1},[A-Z]{1},[A-Z]{1},[A-Z]{1}", ALT)))  %>%
    select(CHROM, POS, REF, ALT, Indiv, gt_GT, gt_GT_alleles)
  
  indel_df <- vcf_df %>%
    mutate(indel_type = case_when(
      nchar(REF) < nchar(ALT) ~ "INS",
      nchar(REF) > nchar(ALT) ~ "DEL",
      grepl("\\*", ALT) ~ "Complex",
      TRUE ~ as.character("Unknown"))) %>%
    separate(gt_GT_alleles, into = c("gt1", "gt2")) %>%
    mutate(indel_type_line = case_when(
      nchar(gt1) > nchar(REF) ~ "INS",
      nchar(gt1) < nchar(REF) ~ "DEL",
      nchar(gt2) > nchar(REF) ~ "INS",
      nchar(gt2) < nchar(REF) ~ "DEL",
      TRUE ~ as.character("REF")
    )) %>% 
    mutate(indel_geno = gt_GT) %>%
    mutate(indel_len_gt1 = abs(nchar(gt1) - nchar(REF))) %>%
    mutate(indel_len_gt2 = abs(nchar(gt2) - nchar(REF))) %>%
    mutate(indel_len = (indel_len_gt1 + indel_len_gt2)/2) %>%
    select(Indiv, CHROM, POS, indel_type, indel_type_line, indel_geno, indel_len) 
  
  file_slug <- gatk_file %>% gsub(".*/", "", .) %>% gsub("_recomb_regions.vcf.gz", "", .)
  
  # save a tidied sm
  saveRDS(indel_df, paste0("analysis_sv/data/gatk_indels/", file_slug, ".rds"))
  
}

gatk_files <- list.files("analysis_sv/data/homology/illumina/", full.names = TRUE)

lapply(gatk_files, create_indel_df)


gatk_indel_files <- list.files("analysis_sv/data/gatk_indels/", full.names = TRUE)

gatk_indel_df <- lapply(gatk_indel_files, readRDS) %>%
  bind_rows %>%
  mutate(pos1 = POS, pos2 = POS, sv_type = indel_type, sv_length = indel_len, genotype = indel_geno, copy_number = 1, depth = 10, method = "gatk", notes = NA) %>%
  rename(ind = Indiv, chrom = CHROM) %>%
  select(-POS, -indel_type, -indel_geno, -indel_len)

saveRDS(gatk_indel_df, "analysis_sv/data/sv_tidy/tidy_gatk.rds")
