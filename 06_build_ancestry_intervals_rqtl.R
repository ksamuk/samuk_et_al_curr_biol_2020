# infer crossovers

library("tidyverse")
library("qtl")
library("parallel")
library("parallelsugar")
library("plotly")

#options(device = "quartz")

# read in rqtl cross
rqtl_cross <- read.cross(file = "data/rqtl/maps/all.csv", format = "csv", genotypes = c("AA", "AB"))

# tidy the cross genotypes into a dataframe
geno_df <- pull.pheno(rqtl_cross) %>% select(-pgm)
geno_tmp <- pull.geno(rqtl_cross) %>% data.frame
geno_df <- cbind(geno_df, geno_tmp) %>% 
  gather(key = marker, value = geno, -id, -sex) %>%
  arrange(id, marker) %>%
  mutate(CHROM_schf = gsub("*._chr_", "", marker) %>% 
           gsub("_pos.*", "", .)) %>%
  mutate(POS_schf = as.numeric(gsub(".*_pos_", "", marker))) %>%
  rename(Indiv = id)

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

# break rank ties using physical position
# accounts for contigs being in reverse order v rqtl positions (arbitrary)
ord_df <- ord_df %>%
  group_by(rqtl_chr) %>%
  mutate(POS_direction = ifelse(lm(POS_schf ~ rqtl_rank) %>% coefficients %>% .[2] < 0, "neg", "pos" )) %>%
  group_by(rqtl_chr, rqtl_rank) %>%
  mutate(rank_2 = ifelse(POS_direction == "neg", dense_rank(-POS_schf), dense_rank(POS_schf))) %>%
  ungroup %>%
  mutate(rqtl_rank = dense_rank(rqtl_rank + 0.1*rank_2)) %>%
  select(-rank_2, -POS_direction) %>%
  arrange(rqtl_chr, rqtl_rank)
  
# join in rqtl marker ranks
geno_df <- geno_df %>%
  left_join(ord_df, by = c("CHROM_schf", "POS_schf")) %>%
  filter(!is.na(rqtl_chr))

# representative individual genotype
geno_df  %>% 
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  filter(grepl("AFC_60", line)) %>%
  filter(sex == "F") %>%
  filter(rqtl_chr == "X") %>%
  #filter(Indiv == "AFC_19_19_G05") %>%
  ggplot(aes(x = as.factor(rqtl_rank), y = Indiv, 
             fill = as.factor(geno)))+
  geom_tile()+
  theme(strip.text.y = element_text(angle = 180),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_discrete(position = "right")+
  xlab("")+
  facet_grid(~rqtl_chr ,scales = "free")

##############################################
# compute ancestry tracts from filtered data
##############################################

ind_list <- geno_df %>% pull(Indiv) %>% as.character %>% unique

#indiv <- ind_list[10]
#chromosome <- "X"
#gt_df <- geno_df
#min_seg_length <- 3

build_ancestry_intervals <- function(indiv, chromosome = NULL, gt_df = NULL, min_seg_length = 3){
  
  cat(indiv)
  cat("\n")
  
  gt_sub <- gt_df %>%
    filter(Indiv == indiv, CHROM_schf == chromosome) %>%
    arrange(rqtl_rank) %>%
    select(rqtl_rank, POS_schf, geno)
    
  
  # if there are too few genotypes, abort mission
  if(nrow(gt_sub) < 10){
    
    return(NULL)
    
  } else{
    
    # compute intervals
    
    runs <- gt_sub$geno %>% rle
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

lapply(ind_list, build_ancestry_intervals, chromosome = "2", gt_df = geno_df)
build_ancestry_intervals(indiv = "MC_6_16_A07", chromosome = "2", gt_df = geno_df)

chr2 <- mclapply(ind_list, build_ancestry_intervals, chromosome = "2", gt_df = geno_df, mc.cores = 4)
chr3 <- mclapply(ind_list, build_ancestry_intervals, chromosome = "3", gt_df = geno_df, mc.cores = 4)
chr4 <- mclapply(ind_list, build_ancestry_intervals, chromosome = "4", gt_df = geno_df, mc.cores = 4)
chrXR <- mclapply(ind_list, build_ancestry_intervals, chromosome = "XR", gt_df = geno_df, mc.cores = 4)
chrXL <- mclapply(ind_list, build_ancestry_intervals, chromosome = "XL", gt_df = geno_df, mc.cores = 4)

chr2 <- lapply(ind_list, build_ancestry_intervals, chromosome = "2", gt_df = gt_df)
chr3 <- lapply(ind_list, build_ancestry_intervals, chromosome = "3", gt_df = gt_df)
chr4 <- lapply(ind_list, build_ancestry_intervals, chromosome = "4", gt_df = gt_df)
chrX <- lapply(ind_list, build_ancestry_intervals, chromosome = "X", gt_df = gt_df)
#chrX[which(is.na(chrX))] <- NULL

seg_df <- bind_rows(chr2, chr3, chr4, chrXR, chrXL)

saveRDS(seg_df, "data/crossover_segments_min_len_2_rqtl.rds")
