# join and filter structural variant calls
# KMS 2019

library("tidyverse")
#library("parallel")
library("parallelsugar")

######################################
# read in tidy SV data
######################################

sv_df <- readRDS("analysis_sv/data/sv_joined.rds")

# make a list of all positions (1 and 2 represented in the sv data)

gt_df <- bind_rows(data.frame(CHROM = sv_df$chrom, POS = sv_df$pos1), 
          data.frame(CHROM = sv_df$chrom, POS = sv_df$pos2)) %>% distinct
  
# the marker orders suggested by shaeffer et al.
scaff_df <- read.csv("meta/scaffold_order_corrections.csv") %>%
  mutate(seg_start = gsub(",", "", start) %>% as.numeric) %>%
  mutate(seg_end = gsub(",", "", end) %>% as.numeric) %>% 
  select(-start, -end) %>%
  mutate(seg_length = ifelse(seg_end > seg_start, seg_end - seg_start, seg_start - seg_end)) %>%
  mutate(direction = ifelse(seg_end < seg_start, "rev", "fwd"))

# create new coordinates in the schaeffer orderings by summing contig lengths

scaff_df <- scaff_df %>%
  group_by(CHROM) %>%
  mutate(new_end = cumsum(seg_length)) %>% 
  mutate(new_start = new_end - seg_length + 1) %>%
  select(CHROM, segment, rank, seg_start, seg_end, seg_length, direction, new_start, new_end) %>%
  ungroup

# function for translating beteween the two coordinate sets
# expects x to be dataframe row index with CHROM and POS

translate_positions <- function(x){
  
  x_df  <- gt_df[x,]
  
  # ignore chr 2 and 3 (these don't need reordering)
  if(x_df$CHROM %in% c(2,3)){
    
    return(data.frame(CHROM_orig = x_df$CHROM, POS_orig = x_df$POS,
                      CHROM_new = x_df$CHROM,  POS_new = x_df$POS))
    
  }
  
  # subset the scaffold data for the target chromosome
  scaff_sub <- scaff_df %>%
    mutate(segment = gsub("Ch", "", segment)) %>%
    filter(segment == x_df$CHROM) %>%
    filter((x_df$POS >= seg_start & x_df$POS <= seg_end) | 
             (x_df$POS <= seg_start & x_df$POS >= seg_end))
  
  # of the marker doesn't exist in the scaffold list, leave it alone
  if(nrow(scaff_sub) != 1){
    
    return(data.frame(CHROM_orig = x_df$CHROM, POS_orig = x_df$POS,
                      CHROM_new = x_df$CHROM,  POS_new = x_df$POS))
    
  }
  
  # if the scaffold is "forward" orientated w/r/t to the v.3 ref:
  if(scaff_sub$direction == "fwd"){
    
    # compute the offset (position from the start)
    # of the marker from the start of the segment
    dist_from_start <- (x_df$POS - scaff_sub$seg_start)
    
    # compute new position by offsetting from the amended startpoint
    POS_new <- dist_from_start + scaff_sub$new_start
    
    # add columns based on the updated positions
    scaff_sub$CHROM_orig <- x_df$CHROM
    scaff_sub$POS_orig <- x_df$POS
    scaff_sub$offset <- dist_from_start
    scaff_sub$CHROM_new <- scaff_sub$CHROM
    scaff_sub$POS_new <- POS_new
    
    return(scaff_sub %>% select(-CHROM))
    
    
  } else{
    
    # if the segment is reversed, then compute the new postion from the 
    # end of the segment 
    
    # compute the offset (position from the start)
    # of the marker from the start of the segment
    dist_from_start <- (x_df$POS - scaff_sub$seg_end)
    
    # compute new position by offsetting from the amended startpoint
    # still the startpoint for reversed segments
    POS_new <- scaff_sub$new_end - dist_from_start 
    
    # return(data.frame(CHROM_orig = x_df$CHROM, POS_orig = x_df$POS,
    #                  CHROM_new = scaff_sub$CHROM,  POS_new = POS_new))
    
    scaff_sub$CHROM_orig <- x_df$CHROM
    scaff_sub$POS_orig <- x_df$POS
    scaff_sub$offset <- dist_from_start
    scaff_sub$CHROM_new <- scaff_sub$CHROM
    scaff_sub$POS_new <- POS_new
    
    return(scaff_sub %>% select(-CHROM))
    
  }
  
}

remap_df <- mclapply(1:nrow(gt_df %>% distinct), translate_positions, mc.cores = 2)
remap_df <- bind_rows(remap_df) 

head(remap_df)

write_rds(remap_df, "analysis_sv/meta/remap.rds")
