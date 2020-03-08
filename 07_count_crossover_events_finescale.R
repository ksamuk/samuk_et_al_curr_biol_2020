library("tidyverse")
library("parallel")

# AFC_14_15_B02	
# options(device = "quartz")

# the ancestry segments
seg_df <- readRDS("data/crossover_segments_min_len_2_rqtl.rds")

# determine which segments are the result of crossovers
# i.e. switches in genotype
co_df <- seg_df %>%
  filter(len_seg >= 2) %>%
  mutate(plate = gsub("_[A-H][0-9]{2}$", "", Indiv))%>%
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  mutate(pop = gsub("_.*$", "", Indiv)) %>%
  group_by(Indiv, CHROM_schf) %>%
  mutate(crossover = ifelse(gt_adj != lag(gt_adj), TRUE, FALSE)) %>%
  mutate(crossover = ifelse(is.na(crossover), FALSE, crossover)) %>%
  ungroup

#############################################################################
# create a per-interval crossover summary for each line
#############################################################################

# want to catch the specific interval (two markers) 
# between which each crossover occured

# in terms of the ancestry segments, that means 
# the last marker of 

# something like:
# indiv line chromosome pos1 pos2 co_count

# re-inflate the segmented data
# (this allows for more precise filtering of short ancestry blocks, i.e. to detect false double XOs)
inf_df <- seg_df %>% 
  group_by(Indiv, CHROM_schf, interval_id) %>%
  mutate(interval_length = end_rank - start_rank) %>%
  filter(interval_length >= 2) %>%
  do(rqtl_rank = .$start_rank:.$end_rank, 
     gt_adj = rep(.$gt_adj, length(.$start_rank:.$end_rank))) %>%
  unnest

head(inf_df)

# re-build the marker order data
ord_df <- read.table("data/rqtl/marker_orders/all.txt", h = T) %>%
  mutate(CHROM_schf = gsub("*._chr_", "", marker) %>% 
           gsub("_pos.*", "", .)) %>%
  mutate(POS_schf = as.numeric(gsub(".*_pos_", "", marker))) %>%
  group_by(chr) %>%
  mutate(rqtl_rank = dense_rank(map_pos)) %>%
  ungroup %>%
  rename(rqtl_chr = chr, rqtl_map_pos = map_pos) %>%
  select(CHROM_schf, POS_schf, rqtl_chr, rqtl_map_pos, rqtl_rank) 

# break rank ties using physical position
# accounts for contigs being in reverse order v rqtl positions (arbitrary)
ord_df <- ord_df %>%
  group_by(rqtl_chr) %>%
  mutate(POS_direction = ifelse(lm(POS_schf ~ rqtl_rank) %>% coefficients %>% .[2] < 0, "neg", "pos" )) %>%
  group_by(rqtl_chr, rqtl_rank) %>%
  mutate(rank_2 = ifelse(POS_direction == "neg", dense_rank(-POS_schf), dense_rank(POS_schf))) %>%
  ungroup %>%
  mutate(rqtl_rank = dense_rank(rqtl_rank + 0.1*rank_2)) %>%
  select(-rank_2, -POS_direction) 

# annotate the inflated (genotype) data with the crossover data
inf_df <- inf_df %>% 
  ungroup %>%
  left_join(ord_df) %>%
  filter(!is.na(rqtl_chr)) %>%
  left_join(co_df %>% select(Indiv, CHROM_schf, interval_id, len_seg, crossover))

# collapse each pair of markers into an interval, and recode crossovers
# at the interval level (TRUE = crossover occured somewhere in that interval)
# write to rds 

inf_df <- inf_df %>%
  select(-interval_id, -len_seg) %>%
  group_by(Indiv, CHROM_schf) %>%
  arrange(Indiv, CHROM_schf, rqtl_rank) %>%
  mutate(POS_schf1 = POS_schf, POS_schf2 = lead(POS_schf)) 
  
# first check: annotate and filter map intervals that are unusually long AND are assembly breakpoints
# these are intervals where the rqtl marker order disagrees (greatly) with the reference
# we can't deal with these in the homology analysis because it relies on reference-called SVs
  
ord_df <- ord_df %>%
  mutate(POS_schf1 = POS_schf, POS_schf2 = lead(POS_schf)) %>% 
  mutate(interval_length = abs(POS_schf1 - POS_schf2)) %>%
  filter(interval_length > 60000) %>%
  mutate(long_interval = interval_length > 1200000) %>%
  mutate(breakpoint = sign(lag(POS_schf1) - lag(POS_schf2)) != sign(POS_schf1 - POS_schf2)) %>% 
  mutate(bad_interval = breakpoint & long_interval) %>% 
  filter(!bad_interval) %>%
  select(-rqtl_rank, -POS_schf, -bad_interval, -breakpoint, -long_interval) 

# second check: intervals can never overlap each other in reference space
# default to keeping first interval defined for overlap

ord_list <- list()

for (j in 1:length(unique(ord_df$CHROM_schf))){
  
  message(unique(ord_df$CHROM_schf)[j])
  
  ord_sub <- ord_df %>%
    filter(CHROM_schf == unique(ord_df$CHROM_schf)[j])
  
  pos_list <- c()
  ord_sub$interval_overlaps_previous <- rep(NA, nrow(ord_sub))
  
  for (i in 1:nrow(ord_sub)){
    
    cat(paste0(ord_sub$POS_schf1[i], "\n"))
    
    row_pos <- ord_sub$POS_schf1[i]:ord_sub$POS_schf2[i]
    
    if(sum(row_pos %in% pos_list) > 10){
      
      ord_sub$interval_overlaps_previous[i] <- TRUE
      
    } else{
      
      ord_sub$interval_overlaps_previous[i] <- FALSE
      pos_list <- c(pos_list, row_pos)
      
    }
    
    
    
  }
  rm(pos_list)
  ord_list[[j]] <- ord_sub
  
}

rm(pos_list)

# apply the overlap filter
ord_df <- ord_list %>%
  bind_rows %>% 
  filter(!interval_overlaps_previous) %>%
  select(rqtl_chr, rqtl_map_pos, CHROM_schf, POS_schf1, POS_schf2)

# all the invervals in this list have passed filter
ord_df$valid_interval <- TRUE

# join the filtered intervals to the inflated genotypes
# remove anything that failed the breakpoint/overlap filters
# write out to rds
inf_df %>% 
  left_join(ord_df) %>%
  filter(!is.na(valid_interval)) %>% 
  mutate(map_pos1 = rqtl_map_pos, map_pos2 = lead(rqtl_map_pos)) %>% 
  mutate(crossover_int = ifelse(lead(crossover) & !(crossover), TRUE, FALSE)) %>% 
  filter(!is.na(POS_schf2) & !is.na(POS_schf1)) %>%
  select(-rqtl_rank, -rqtl_chr, -POS_schf, -crossover, -valid_interval) %>%
  saveRDS("data/crossovers_fine_scale.rds")

plot_df <- inf_df %>% 
  left_join(ord_df) %>%
  filter(!is.na(valid_interval)) %>% 
  mutate(map_pos1 = rqtl_map_pos, map_pos2 = lead(rqtl_map_pos)) %>% 
  mutate(crossover_int = ifelse(lead(crossover) & !(crossover), TRUE, FALSE)) %>% 
  filter(!is.na(POS_schf2) & !is.na(POS_schf1)) %>%
  select(-rqtl_rank, -rqtl_chr, -POS_schf, -crossover, -valid_interval) 

plot_df %>%
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  filter(grepl("AFC_60", line)) %>%
  filter(grepl("_B", Indiv)) %>%
  filter(CHROM_schf == "2") %>%
  #filter(Indiv == "AFC_19_19_G05") %>%
  ggplot(aes(x = as.factor(POS_schf1), y = Indiv, 
             fill = as.factor(gt_adj)))+
  geom_tile()+
  theme(strip.text.y = element_text(angle = 180),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_y_discrete(position = "right")+
  xlab("")+
  facet_grid(~CHROM_schf ,scales = "free")


inf_df <- readRDS("data/crossovers_fine_scale.rds")

# chr2
inf_df %>%
  select(Indiv, CHROM_schf, POS_schf1, crossover_int) %>%
  setNames(nm = c("Individual", "Chromosome", "Position", "Crossover")) %>%
  write.csv("figures/TableS2.csv", row.names = FALSE, quote = FALSE)

