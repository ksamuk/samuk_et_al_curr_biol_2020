# infer crossovers

library("tidyverse")
library("parallel")
library("parallelsugar")
library("plotly")

#options(device = "quartz")

# read in filtered data
gt_df <- readRDS("data/gt_df/gt_seq_filtered.rds")

# read in rQTL marker orders
ord_df <- read.table("data/rqtl/marker_orders/all.txt", h = T) %>%
  mutate(CHROM_schf = gsub("*._chr_", "", marker) %>% 
                 gsub("_pos.*", "", .)) %>%
  mutate(POS_schf = as.numeric(gsub(".*_pos_", "", marker))) %>%
  group_by(chr) %>%
  mutate(rqtl_rank = dense_rank(map_pos)) %>%
  ungroup %>%
  rename(rqtl_chr = chr, rqtl_map_pos = map_pos) %>%
  select(CHROM_schf, POS_schf, rqtl_chr, rqtl_map_pos, rqtl_rank) 

# join in rqtl marker ranks
gt_df <- gt_df %>%
  left_join(ord_df, by = c("CHROM_schf", "POS_schf")) %>%
  filter(!is.na(rqtl_chr))

# representative individual genotype
gt_df  %>% 
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  filter(grepl("AFC_60", line)) %>%
  filter(sex == "F") %>%
  filter(unique_in != "flag") %>%
  filter(rqtl_chr == "X") %>%
  #filter(Indiv == "AFC_19_19_G05") %>%
  ggplot(aes(x = as.factor(rqtl_rank), y = Indiv, 
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
  facet_grid(~rqtl_chr ,scales = "free")

##############################################
# compute ancestry tracts from filtered data
##############################################

ind_list <- gt_df %>% pull(Indiv) %>% unique

build_ancestry_intervals <- function(indiv, chromosome = NULL, gt_df = NULL, min_seg_length = 3){
  
  cat(indiv)
  cat("\n")
  
  gt_sub <- gt_df %>%
    filter(Indiv == indiv, CHROM_schf == chromosome) %>%
    arrange(rqtl_rank) %>%
    select(rqtl_rank, POS_schf, gt_adj)
    
  
  # if there are too few genotypes, abort mission
  if(nrow(gt_sub) < 10){
    
    return(NULL)
    
  } else{
    
    # compute intervals
    
    runs <- gt_sub$gt_adj %>% rle
    end <- cumsum(runs$lengths)
    start <- c(1, lag(end)[-1] + 1)
    
    # tabulate intervals
    
    
    out_df <- data.frame(Indiv = indiv, CHROM_schf = chromosome, interval_id = 1:length(start), 
               start_rank = gt_sub$rqtl_rank[start], end_rank = gt_sub$rqtl_rank[end], 
               start_POS = gt_sub$POS_schf[start], end_POS = gt_sub$POS_schf[end], 
               gt_adj = runs$values)
    
    # apply threshold for min interval size 
    # interpolate across small intervals based on surroundings
    # merge down any redundant segments
    
    out_df <- out_df %>%
      mutate(len_seg = end_rank - start_rank + 1) %>%
      filter(len_seg >= min_seg_length)
    
    # inject POS string for easing joining later
    
    out_df$pos_string <- NA
    
    for (i in 1:nrow(out_df)){
      
      pos_list <- gt_sub$POS_schf[gt_sub$rqtl_rank >= out_df$start_rank[i] & 
                             gt_sub$rqtl_rank <= out_df$end_rank[i]]
      
      out_df$pos_string[i] <- paste(pos_list, collapse = ",")
      
    }
    
    out_df
    
  }
  
}

#lapply(ind_list, build_ancestry_intervals, chromosome = "2", gt_df = gt_df)
#build_ancestry_intervals(indiv = "MC_6_16_A07", chromosome = "2")

chr2 <- mclapply(ind_list, build_ancestry_intervals, chromosome = "2", gt_df = gt_df, mc.cores = 3)
chr3 <- mclapply(ind_list, build_ancestry_intervals, chromosome = "3", gt_df = gt_df, mc.cores = 3)
chr4 <- mclapply(ind_list, build_ancestry_intervals, chromosome = "4", gt_df = gt_df, mc.cores = 3)
chrXR <- mclapply(ind_list, build_ancestry_intervals, chromosome = "XR", gt_df = gt_df, mc.cores = 3)
chrXL <- mclapply(ind_list, build_ancestry_intervals, chromosome = "XL", gt_df = gt_df, mc.cores = 3)

# chr2 <- lapply(ind_list, build_ancestry_intervals, chromosome = "2", gt_df = gt_df)
# chr3 <- lapply(ind_list, build_ancestry_intervals, chromosome = "3", gt_df = gt_df)
# chr4 <- lapply(ind_list, build_ancestry_intervals, chromosome = "4", gt_df = gt_df)
#chrX <- lapply(ind_list, build_ancestry_intervals, chromosome = "X", gt_df = gt_df)
#chrX[which(is.na(chrX))] <- NULL

seg_df <- bind_rows(chr2, chr3, chr4, chrXR, chrXL)

saveRDS(seg_df, "data/crossover_segments_min_len_2.rds")
