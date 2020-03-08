# combine crossover data with miller reference positions
# to get corrected cM/MB
# KMS Mar 2019

library("tidyverse")

###########################################################
# read and format datasets 
###########################################################

# read the blast marker alignments
miller_pos <- read.table("analysis_cm_mb/out/gt_seq_miller_blast.txt", h = T)

# read the map distances from the rqtl "all" map
map_dat <- read.csv("data/rqtl/maps/all.csv", nrows = 3, h = F) %>% t 

map_dat <- map_dat[-c(1:3),] %>%
  data.frame

names(map_dat) <- c("marker", "CHROM_rqtl", "map_pos_cM")
row.names(map_dat) <- NULL

# read the coverted positions based on schaeffer 2008
remap_df <- readRDS("analysis_cm_mb/data/schaffer2008_remapped_positions.rds")

# read the amplicon sequence meta data
# (some amplicons ended up not being used)
amp_dat <- read.table("meta/amplicon_meta_master.txt", h = T)

###########################################################
# filtration of sites 
###########################################################

# harmonize the pair ids
# filter blast results for actively used markers
miller_pos <- miller_pos %>%
  mutate(pair_id = gsub(".*\\|", "", QueryID) %>% gsub("_[0-9]*:.*", "", .)) %>%
  filter(pair_id %in% as.character(amp_dat$pair_id))

# chose the best hit (max length seems to work) from the blast results for each marker
# drop unused columns
miller_pos <- miller_pos %>%
  group_by(pair_id) %>%
  filter(Alignment.Length == max(Alignment.Length)) %>%
  ungroup %>%
  select(pair_id, SubjectID, S.start, S.end) %>%
  rename(miller_contig = SubjectID, miller_pos1 = S.start, miller_pos2 = S.end)

miller_pos <- miller_pos %>%
  mutate(miller_contig = gsub("Consensus_Consensus_Consensus_|_pilon_pilon_pilon", "", miller_contig))

# harmonize map dat marker names

amp_dat <- amp_dat %>%
  mutate(marker = paste0("chr_", CHROM, "_pos_", POS)) %>%
  mutate(marker = ifelse(grepl("XR|XL", as.character(CHROM)), 
                         paste0("X_", marker), 
                         paste0("A_", marker)))

# create the unified position list
# whew
pos_df <- map_dat %>% 
  left_join(amp_dat) %>%
  left_join(remap_df) %>%
  mutate(map_pos_cM = as.numeric(as.character(map_pos_cM))) %>%
  left_join(miller_pos %>% select(-miller_pos2)) %>%
  select(pair_id, CHROM, POS, chrom_schaeffer, pos_schaeffer, miller_contig, miller_pos1, CHROM_rqtl, map_pos_cM)

# compute intermarker physical distances
pos_df <- pos_df %>%
  mutate(dist_map = ifelse(lag(CHROM_rqtl) == CHROM_rqtl, abs(map_pos_cM - lag(map_pos_cM)), NA)) %>%
  mutate(dist_ref = ifelse(lag(CHROM) == CHROM, abs(POS - lag(POS)), NA))%>%
  mutate(dist_sch = ifelse(lag(chrom_schaeffer) == chrom_schaeffer, abs(pos_schaeffer - lag(pos_schaeffer)), NA))%>%
  mutate(dist_miller = ifelse(lag(miller_contig) == miller_contig, abs(miller_pos1 - lag(miller_pos1)), NA)) %>%
  mutate(dev_sch = abs(dist_ref - dist_sch)) %>%
  mutate(dev_miller = abs(dist_ref - dist_miller))
  

# plot   
pos_df %>%
  mutate(cm_mb_schf = dist_map/(dist_sch/1000000)) %>%
  mutate(cm_mb_miller = dist_map/(dist_miller/1000000)) %>%
  filter(cm_mb_schf < 100) %>%
  ggplot(aes(x = pos_schaeffer, y = cm_mb_miller))+
  geom_point()+
  facet_wrap(~chrom_schaeffer)



  

