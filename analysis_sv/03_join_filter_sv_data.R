# join and filter structural variant calls
# KMS 2019

library("tidyverse")

select <- dplyr::select

######################################
# read in tidy SV data
######################################

pb_df <- readRDS("analysis_sv/data/sv_tidy/tidy_pbsv.rds")
pte_df <- readRDS("analysis_sv/data/sv_tidy/tidy_popte2.rds")
smv_df <- readRDS("analysis_sv/data/sv_tidy/tidy_smoove.rds")
gatk_df <- readRDS("analysis_sv/data/sv_tidy/tidy_gatk.rds")

######################################
# merging SVs 
######################################

# basic ideas:
# - for pb and smoove, filtering on depth as a minimum (popte2 doesn't provide depth annots)
# - subanalysis of SVs that have support from at least two methods (e.g. malapbsv + smoove)
# - empirical cutoffs for bizzare svs (outliers in annotation values for their SV class)
# - idea for emp cutoffs: mahalanobis distance (for each class separately)

# <do addition pre filtering here if necessary>

# simplified representation of SVs
# ind chrom pos1 pos2 sv_type source

# subset and rename SV datasets for merging

pb_sub <- pb_df %>%
  filter(FILTER == "PASS") %>%
  mutate(SVLEN = ifelse(SVTYPE == "INV", abs(END - POS), SVLEN)) %>%
  mutate(sv_type_line = SVTYPE) %>%
  select(Indiv, CHROM, POS, END, SVTYPE, sv_type_line, SVLEN, gt_GT,  gt_CN,  gt_DP) %>%
  mutate(method = "pbsv", notes = NA)%>%
  mutate(Indiv = as.character(Indiv)) %>%
  mutate(SVLEN = as.numeric(SVLEN))

names(pb_sub) <- c("ind", "chrom", "pos1", "pos2", "sv_type", "sv_type_line", "sv_length", "genotype", "copy_number", "depth", "method", "notes")

# apply filters to smv data:
# there are a number of very large SVs called by smoove
# however these are not supported by the pacbio data
# manual inspection of the raw reads seems to support the existance of the break points, 
# it seems odd to "connect the dots" and assume the breakpoints define very large SVs (e.g. a 6MB deletion?)
# so, until we have better data here, I am dropping very large SVs called only by smoove

# also filtering on allele ratio, to mirror treatment of SNPs/indels
smv_sub <- smv_df %>%
  mutate(SU = as.numeric(SU)) %>%
  filter(gt_GQ > 50,  gt_DP > 10, SU > 100) %>%
  filter(SVTYPE != "BND") %>% 
  mutate(sv_type_line = SVTYPE) %>%
  mutate(SVLEN = as.numeric(SVLEN)) %>% 
  filter(abs(SVLEN) < 150000) %>%
  select(Indiv, CHROM, POS, END, SVTYPE, sv_type_line, SVLEN, gt_GT, gt_DP)%>%
  mutate(copy_number = NA, method = "smoove", notes = NA)%>%
  select(Indiv, CHROM, POS, END, SVTYPE, sv_type_line, SVLEN, gt_GT, copy_number, gt_DP, method, notes)%>%
  mutate(Indiv = as.character(Indiv)) 

names(smv_sub) <- c("ind", "chrom", "pos1", "pos2", "sv_type", "sv_type_line", "sv_length", "genotype", "copy_number", "depth", "method", "notes")

pte_sub <- pte_df %>%
  mutate(sv_type = "TE", sv_length = NA, genotype = "1/1", copy_number = NA, method = "popte2",  depth = NA) %>%
  mutate(sv_type_line = sv_type) %>%
  mutate(notes = paste(family, te_type, direction, sep = "_")) %>%
  mutate(pos2 = pos1 + 1 ) %>%
  select(Indiv, CHROM, pos1, pos2, sv_type, sv_type_line, sv_length, genotype, copy_number, depth, method, notes) %>%
  mutate(Indiv = as.character(Indiv))
names(pte_sub) <- c("ind", "chrom", "pos1", "pos2", "sv_type", "sv_type_line", "sv_length", "genotype", "copy_number", "depth", "method", "notes")


# one df to rule them all!

sv_df <- bind_rows(pb_sub, pte_sub, smv_sub, gatk_df)

######################################
# basic look at depth
######################################

#sv_df %>%
#  filter(!is.na(depth)) %>%
#  ggplot(aes(x = depth))+
#  geom_histogram()+
#  facet_wrap(~method, scales = "free")

# clearly some intense skew 
# opted to filter out the extreme depth outliers 
# on a per-method basis

sv_df %>%
  filter(!is.na(depth)) %>%
  filter(depth > 2) %>%
  filter((method == "pbsv" & depth < 25) | (method == "smoove" & depth < 100) | (method == "gatk" & depth < 100)) %>%
  ggplot(aes(x = depth))+
  geom_histogram()+
  facet_wrap(~method, scales = "free")
  
# this results in much more plausible depth distributions
sv_df <- sv_df %>%
  filter((method == "pbsv" & depth < 25 & depth > 4) | (method == "smoove" & depth < 100 & depth > 4) | method == "popte2" | method == "gatk") 

######################################
# repair individual ids 
######################################

# read in individual labels
ind_names_sr <- read.table("meta/fasta_samplelist_recode.txt", header = FALSE) 
names(ind_names_sr) <- c("ind", "ind2")
ind_names_sr $ind <- gsub("S", "", ind_names_sr $ind) %>% as.numeric
ind_names_sr$label_method <- "sr"

# pb sample names
ind_names_pb <- read.table("analysis_sv/meta/pac_bio_sample_ids.txt", header = TRUE) 
names(ind_names_pb) <- c("ind", "ind2")
ind_names_pb$label_method <- "pb"

ind_names <- bind_rows(ind_names_pb, ind_names_sr)

# clean and normalize individual labels
sv_df <- sv_df %>%
  mutate(ind = gsub("S", "", ind) %>% as.numeric) %>%
  mutate(label_method = ifelse(method == "pbsv", "pb", "sr")) %>%
  left_join(ind_names) %>%
  select(-ind) %>%
  rename(ind = ind2) %>%
  select(ind, everything()) %>%
  mutate(ind = gsub("FC|C", "", ind)) %>%
  mutate(ind = gsub("MV225", "MV2-25", ind)) %>%
  mutate(ind = gsub("VY-F16", "VY", ind)) %>%
  select(-label_method)

saveRDS(sv_df, "analysis_sv/data/sv_joined.rds")
  
  
  


