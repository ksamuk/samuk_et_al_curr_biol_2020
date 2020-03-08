# determine genotypic sex for each individual 
# KMS Nov 2018

library("tidyverse")
select <- dplyr::select

#############################################################
# basic data filtration (from 02_ script)
#############################################################

# read in amplicon target region meta data
sites_df <- read.table("meta/snp_pair_meta_data.txt", h = T)

# read in snp identities
snp_df <- readRDS("meta/sampled_gt_seq_snps_all.rds") %>%
  select(CHROM, POS, unique_in, pair_id) %>%
  distinct

# read in list of actually used primers
prim_df <- read.csv("meta/ordered_primers.csv", h = F) 
names(prim_df) <- "slug"
prim_df <- prim_df %>%
  mutate(pair_id = gsub("_[0-9]*$", "", slug))

# determine allele identities for polarizing ancestry
# (sometimes the "derived" allele is actually the unique MV/Flag allele)
snp_ident <- readRDS("meta/sampled_gt_seq_snps_all.rds") %>%
  #select(CHROM, POS, unique_in, pair_id) %>%
  filter(Indiv %in% c("Flag14", "MV2-25")) %>%
  select(CHROM, POS, pair_id, Indiv, unique_in, gt_GT_alleles) %>%
  spread(key = Indiv, value = gt_GT_alleles)

names(snp_ident)[5:6] <- c("flag_alelle", "mv_allele") 

# compute first snp in pair (for plotting)
pair_df <- readRDS("meta/sampled_gt_seq_snps_all.rds") %>%
  group_by(CHROM, pair_id) %>%
  summarise(POS = mean(POS), first_pos = min(POS)) %>%
  ungroup %>%
  distinct %>%
  select(pair_id, first_pos)

snp_df <- left_join(snp_df, pair_df, by = c("pair_id"))
snp_df <- left_join(snp_df, snp_ident)

# only allow snps that were in the actual primer set (!)
snp_df <- snp_df %>% 
  filter(pair_id %in% prim_df$pair_id )

# read in the gt seq vcf
gt_df <- readRDS(file = "data/gt_df/gt_seq_all_gt_df_gatk.rds")

# join in unique_in info for each snp
# (note that there are lots of genotyped sites that are non-target, and thus will be NA here)
# re: the distinct, the upstream script generated duplicated rows somehow -- potentially a bcftools thing?

gt_df <- gt_df %>%
  left_join(snp_df, by = c("CHROM", "POS")) %>%
  distinct %>%
  arrange(CHROM, POS, Indiv) 

# mean representation of individuals
depth_df <- gt_df %>%
  group_by(plate, Indiv) %>%
  summarise(mean_dp = mean(gt_DP, na.rm = TRUE))

# join into gt_df for filtering
gt_df <- left_join(gt_df, depth_df)

# also pull out the flag sites
gt_df_flag <- gt_df %>%
  filter(!is.na(unique_in)) %>%
  filter(mean_dp > 30) %>%
  filter(ALT %in% c("A", "C", "T", "G")) %>%
  filter(REF %in% c("A", "C", "T", "G")) %>%
  filter(nchar(gt_GT_alleles) <=3) %>%
  #filter(!(unique_in == "flag" & (gt_GT == "0/0" |gt_GT == "1/1" ))) %>%
  #filter(!(unique_in == "mv" & gt_GT == "0/0")) %>%
  filter(unique_in == "flag") %>%
  mutate(gt_adj = ifelse(gt_GT_alleles == flag_alelle, "0/0", gt_GT))%>%
  mutate(gt_adj = ifelse((gt_GT_alleles != flag_alelle) & (gt_GT == "0/0"), "1/1", gt_adj)) # correct for shared non-ref alleles

#############################################################
# sexing individuals via homozygosity (using flag alleles)
#############################################################

sex_df <- gt_df_flag %>%
  select(CHROM, POS, Indiv, gt_adj) %>% 
  filter(grepl("X.*", CHROM)) %>%
  group_by(Indiv) %>%
  summarise(mean_hom = mean(!(gt_adj == "0/1"), na.rm = TRUE)) %>%
  mutate(sex = ifelse(mean_hom >= 0.90, "M", "F")) %>% 
  filter(mean_hom >= 0.95 | mean_hom <= 0.20)

write.csv(sex_df, "data/sex_df.csv", row.names = FALSE)



