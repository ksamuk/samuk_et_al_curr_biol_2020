# format genetic map data for input into rqtl 
# KMS Jan 2019

library("tidyverse")
select <- dplyr::select

# read in filtered and polarized genotypes
gt_df <- readRDS("data/gt_df/gt_seq_filtered.rds")

# and sex_df (needed later)
sex_df <- read.csv("data/sex_df.csv") %>%
  select(Indiv, sex)

# subset the master data file
# remove two plates with barcode errors

rqtl_df <- gt_df %>%
  filter(!(plate %in% c("AFC_19_17", "AFC_60_19"))) %>%
  filter(!is.na(sex)) %>%
  filter(!is.na(unique_in)) %>%
  filter(unique_in == "mv") %>%
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  select(line, plate, sex, Indiv, CHROM_schf, POS_schf, gt_adj) %>%
  arrange(line, plate, Indiv, CHROM_schf, POS_schf) %>%
  distinct %>%
  filter(!is.na(gt_adj)) %>%
  rename(gt_GT = gt_adj)

# check for missingness due to polarization failure
gt_df_sites <- gt_df %>% 
  filter(unique_in == "mv") %>%
  mutate(marker = paste0(CHROM_schf, "_", POS_schf)) %>% 
  pull(marker) %>% 
  unique

rqtl_df_sites <- rqtl_df %>%
  mutate(marker = paste0(CHROM_schf, "_", POS_schf)) %>% 
  pull(marker) %>% 
  unique

# should be a big barrel of TRUE
gt_df_sites %in% rqtl_df_sites


# determine sex-linked markers
auto_x_df <- rqtl_df %>%
  group_by(sex, CHROM_schf, POS_schf) %>%
  summarise(p_hom = sum(gt_GT == "0/0" | gt_GT == "1/1") / length(gt_GT)) %>%
  spread(key = sex, value = p_hom) %>% 
  mutate(hom_diff = M - `F`) %>%
  mutate(chr = ifelse(M > 0.99 & hom_diff > 0.01, "X", "A")) %>%
  select(CHROM_schf, POS_schf, chr)

# convert to AB format
rqtl_df <- rqtl_df %>%  
  left_join(auto_x_df, by = c("CHROM_schf", "POS_schf")) %>%
  mutate(gt_GT = ifelse(gt_GT == "1/1", "BB", gt_GT))%>%
  mutate(gt_GT = ifelse(gt_GT == "0/1", "AB", gt_GT))%>%
  mutate(gt_GT = ifelse(gt_GT == "0/0", "AA", gt_GT))%>%
  mutate(gt_GT = ifelse(is.na(gt_GT), "-", gt_GT))%>%
  mutate(gt_GT = ifelse(gt_GT == "NA", "-", gt_GT))

# how many sites?
rqtl_df %>%
  select(CHROM_schf,  POS_schf) %>%
  distinct %>%
  tally

# function for creating the input files for rqtl

build_rqtl_input_files <- function(line_id) {
  
  message(line_id)
  
  # build the phenotype file for rqtl
  sex_df %>%
    select(sex, id = Indiv) %>%
    filter(grepl(paste0(line_id, "_"), id)) %>%
    data.frame(pgm = 1) %>%
    write.csv(paste0("data/rqtl/raw/", line_id, "_pheno_dat.csv"), row.names = FALSE)
  
  # build the gt table
  gt_tab <- rqtl_df %>%
    filter(!is.na(chr)) %>%
    filter(grepl(paste0(line_id, "_"), Indiv)) %>%
    mutate(marker = paste0(chr, "_chr_",CHROM_schf, "_pos_", POS_schf)) %>%
    select(Indiv, marker, gt_GT) %>%
    group_by(Indiv) %>%
    filter(!duplicated(marker)) %>%
    ungroup %>%
    spread(key = marker, value = gt_GT) %>%
    select(id = Indiv, everything())
  
  gt_tab[is.na(gt_tab)] <- "-"
  
  # filter out sites with low coverage
  #geno_counts <- lapply(gt_tab, function(x) ((x != "-") %>% sum) > 100) %>% unlist
  #geno_counts[1] <- TRUE
  #geno_counts <- gt_tab[,geno_counts]
  
  # the chromosome names line
  chr_names <- names(gt_tab) %>%
    gsub("_.*", "", .) 
    #gsub("XL|XR", "X", .) # re-enable to preserve X chr
  
  chr_names[1] <- ""
  names(chr_names) <- names(gt_tab)
  
  header_lines <- data.frame(chr_names) %>% t
  
  # build the genotype file for rqtl
  header_lines %>%
    write.table(paste0("data/rqtl/raw/", line_id, "_geno_dat.csv"), row.names = FALSE, sep = ",")
  
  gt_tab %>%
    filter(grepl(paste0(line_id, "_"), id)) %>%
    write.table(paste0("data/rqtl/raw/", line_id, "_geno_dat.csv"), append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  
  
}

# all the line_ids
line_ids <- rqtl_df$Indiv %>% unique %>% gsub("_[0-9]*_.{3}$", "", .) %>% unique
lapply(line_ids, build_rqtl_input_files)


### combined file for combined map

# build the phenotype file for rqtl
sex_df %>%
    select(sex, id = Indiv) %>%
    data.frame(pgm = 1) %>%
    write.csv(paste0("data/rqtl/raw/all_pheno_dat.csv"), row.names = FALSE)
  
# build the gt table
gt_tab <- rqtl_df %>%
  filter(!is.na(chr)) %>%
  mutate(marker = paste0(chr, "_chr_",CHROM_schf, "_pos_", POS_schf)) %>%
  select(Indiv, marker, gt_GT) %>%
  group_by(Indiv) %>%
  filter(!duplicated(marker)) %>%
  ungroup %>%
  spread(key = marker, value = gt_GT) %>%
  select(id = Indiv, everything())
  
gt_tab[is.na(gt_tab)] <- "-"
  
  # filter out sites with low coverage
  #geno_counts <- lapply(gt_tab, function(x) ((x != "-") %>% sum) > 100) %>% unlist
  #geno_counts[1] <- TRUE
  #geno_counts <- gt_tab[,geno_counts]
  
# the chromosome names line
chr_names <- names(gt_tab) %>%
  gsub("_.*", "", .) %>%
  gsub("A", "1", .)

chr_names[1] <- ""
names(chr_names) <- names(gt_tab)
  
header_lines <- data.frame(chr_names) %>% t
  
  # build the genotype file for rqtl
header_lines %>%
  write.table(paste0("data/rqtl/raw/all_geno_dat.csv"), row.names = FALSE, sep = ",")
  
gt_tab %>%
  write.table(paste0("data/rqtl/raw/all_geno_dat.csv"), append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
