# prepare cm/mb file for k korunes 
# kms apr 2019

library("tidyverse")

# read in map positon data
ord_df <- read.table("data/rqtl/marker_orders/all.txt", h = T) %>%
  mutate(CHROM_schf = gsub("*._chr_", "", marker) %>% 
           gsub("_pos.*", "", .)) %>%
  mutate(POS_schf = as.numeric(gsub(".*_pos_", "", marker))) %>%
  group_by(chr) %>%
  mutate(rqtl_rank = dense_rank(map_pos)) %>%
  ungroup %>%
  rename(rqtl_chr = chr, rqtl_map_pos = map_pos) %>%
  select(CHROM_schf, POS_schf, rqtl_chr, rqtl_map_pos, rqtl_rank) 

# read in remapped position data
remap_df <- readRDS("analysis_cm_mb/out/remapped_positions.rds")

map_df <- left_join(ord_df, remap_df)

head(map_df)

# construct recombination intervals

interval_list <- vector("list", length(unique(map_df$CHROM_schf)))

for (i in 1:length(unique(map_df$CHROM_schf))){
  
  map_chr_sub <- map_df %>%
    filter(CHROM_schf == unique(map_df$CHROM_schf)[i])
  
  phys_positions <- unique(map_chr_sub$POS_schf)

  # construct pairwise map and phyiscal distances
  for (j in 1:length(phys_positions)){
    
    marker_sub <- map_chr_sub %>%
      filter(POS_schf %in% c(phys_positions[j], phys_positions[j + 1]))
    
    interval_list[[i]][[j]] <- data.frame(
               CHROM_ref1 = marker_sub$CHROM_ref[1], CHROM_ref2 = marker_sub$CHROM_ref[2],
               POS_ref1 = marker_sub$POS_ref[1], POS_ref2 = marker_sub$POS_ref[2],
               CHROM_schf1 = marker_sub$CHROM_schf[1], CHROM_schf2 = marker_sub$CHROM_schf[2],
               POS_schf1 = marker_sub$POS_schf[1], POS_schf2 = marker_sub$POS_schf[2], 
               CHROM_miller1 = marker_sub$CHROM_miller[1], CHROM_miller2 = marker_sub$CHROM_miller[2],
               POS_miller1 = marker_sub$POS_miller[1], POS_miller2 = marker_sub$POS_miller[2],
               rqtl_chr1 = marker_sub$rqtl_chr[1], rqtl_chr2 = marker_sub$rqtl_chr[2],
               map_pos1 = marker_sub$rqtl_map_pos[1], map_pos2 = marker_sub$rqtl_map_pos[2])  
      
      
  }
  
  interval_list[[i]] <- bind_rows(interval_list[[i]])
    
}

interval_df <- interval_list %>%
  bind_rows 

interval_df <- interval_df %>%
  mutate(cm_dist_map = ifelse(rqtl_chr1 == rqtl_chr2, abs(map_pos1 - map_pos2))) %>%
  mutate(mb_dist_schf = ifelse(CHROM_schf1 == CHROM_schf2, abs(POS_schf1 - POS_schf2), NA))%>%
  mutate(mb_dist_ref = ifelse(CHROM_ref1 == CHROM_ref2, abs(POS_ref1 - POS_ref2), NA))%>%
  mutate(mb_dist_miller = ifelse(CHROM_miller1 == CHROM_miller2, abs(POS_miller1 - POS_miller2), NA))

options(scipen = 999)

interval_df %>%
  ggplot(aes(x = mb_dist_schf, y = mb_dist_miller))+
  geom_point() + 
  xlab("MB distance in Schaeffer et al. 2008")+
  ylab("MB distance in Miller et al. 2018")

interval_df %>%
  ggplot(aes(x = mb_dist_schf, y = mb_dist_ref))+
  geom_point()+ 
  xlab("MB distance in Schaeffer et al. 2008")+
  ylab("MB distance in Dpse Reference V3")

cmmb_df <- interval_df %>%
  mutate(cm_mb_ref = ifelse(cm_dist_map == 0, 0, (cm_dist_map / mb_dist_ref) * 1000000)) %>%
  mutate(cm_mb_schf = ifelse(cm_dist_map == 0, 0, (cm_dist_map / mb_dist_schf) * 1000000)) %>%
  mutate(cm_mb_miller = ifelse(cm_dist_map == 0, 0, (cm_dist_map / mb_dist_miller) * 1000000))

# write the full raw
cmmb_df %>%
  #filter(!is.na(cm_mb_miller) & !is.infinite(cm_mb_miller)& ! is.nan(cm_mb_miller)) %>%
  filter(mb_dist_schf > 100000) %>% 
  write.table("analysis_cm_mb/out/dpse_recomb_conservative.txt", quote = FALSE, row.names = FALSE)

cmmb_df %>%
  #filter(!is.na(cm_mb_miller) & !is.infinite(cm_mb_miller)& ! is.nan(cm_mb_miller)) %>%
  filter(mb_dist_schf > 100000) %>% 
  write.table("analysis_cm_mb/out/dpse_recomb_all.txt", quote = FALSE, row.names = FALSE)  
  ggplot(aes(x = (POS_schf1+POS_schf2)/2, y = cm_mb_schf))+
  geom_step()+
  facet_wrap(~CHROM_schf1)




