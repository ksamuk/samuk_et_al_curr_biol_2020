# get three-way positions for the reference, adjusted orders (schf 2008)
# and miller 2018 positions (the latter are approximately, to ~500bp)
# KMS Mar 2019

library("tidyverse")

###########################################################
# read and format datasets 
###########################################################

# read the blast marker alignments
miller_pos <- read.table("analysis_cm_mb/out/gt_seq_miller_blast.txt", h = T)

# read the coverted positions based on schaeffer 2008
remap_df <- readRDS("analysis_cm_mb/data/schaffer2008_remapped_positions.rds")

###########################################################
# filtration of sites 
###########################################################

# harmonize the pair ids
# filter blast results for actively used markers
miller_pos <- miller_pos %>%
  mutate(CHROM = gsub("_[0-9]*:[0-9]*|", "", QueryID) %>%
           gsub("all\\|", "", .)  %>% gsub("_[0-9]*$", "", .)) %>%
  mutate(POS = gsub(".*_", "", QueryID)) %>%
  mutate(ref_pos1 = gsub(".*:", "", POS))%>%
  mutate(ref_pos2 = gsub(":.*", "", POS)) 

# chose the best hit (max length seems to work) from the blast results for each marker
# drop unused columns
miller_pos <- miller_pos %>%
  group_by(CHROM, POS) %>%
  filter(Alignment.Length == max(Alignment.Length)) %>%
  ungroup %>%
  select(CHROM, ref_pos1, ref_pos2, SubjectID, S.start, S.end) %>%
  rename(miller_contig = SubjectID, miller_pos1 = S.start, miller_pos2 = S.end) %>%
  mutate(miller_contig = gsub("Consensus_Consensus_Consensus_|_pilon_pilon_pilon", "", miller_contig))


# expand the site ranges for joining
# the miller positions are approximate within 500bp
# converting the blast results to base-wise positions ended up being way too complicated
# next time, use bwa
miller_pos <- miller_pos %>%
  rowwise %>%
  do(CHROM_ref = rep(.$CHROM, length(.$ref_pos1:.$ref_pos2)), 
     POS_ref = .$ref_pos1:.$ref_pos2, 
     CHROM_miller = rep(.$miller_contig, length(.$ref_pos1:.$ref_pos2)),
     POS_miller = rep(floor(mean(c(.$miller_pos1, .$miller_pos2))), length(.$ref_pos1:.$ref_pos2))) %>%
  unnest()
  
# create the unified position list
# whew
pos_df <- remap_df %>% 
  rename(CHROM_ref = CHROM, POS_ref = POS, CHROM_schf = chrom_schaeffer, POS_schf = pos_schaeffer) %>%
  left_join(miller_pos) 

write_rds(pos_df, "analysis_cm_mb/out/remapped_positions.rds")


  

