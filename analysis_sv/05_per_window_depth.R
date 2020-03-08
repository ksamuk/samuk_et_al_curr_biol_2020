# compute per window depth

# extract per-window estimates of sequencing/mapping depth
# for both the PacBio and Illumina data
# KMS Aug 2019

library("tidyverse")
library("R.utils")
library("vcfR")

# the list of recombination intervals from shaeffer et al.
recomb_intervals <- read.table("analysis_sv/meta/recom_intervals_ref_schf.txt", h = T) %>%
  rowwise %>%
  mutate(pos1_ref = min(POS_ref1, POS_ref2),
         pos2_ref = max(POS_ref1, POS_ref2)) %>%
  mutate(POS_schf1 = pos1_schf, POS_schf2 = pos2_schf) %>%
  select(CHROM_ref, pos1_ref, pos2_ref, POS_schf1, POS_schf2, CHRO) %>% 
  ungroup 


###################################################
# pac bio
###################################################

pb_keep_files <- list.files("analysis_sv/data/coverage_pb/tmp", full.names = TRUE)

compress_depth_pb <- function(pb_keep_file){
  
  message(pb_keep_file)
  
  # get the ind id
  ind <- pb_keep_file %>% gsub(".*/", "", .) %>% gsub("_keep.txt", "", .)
  
  if(!file.exists(paste0("analysis_sv/data/coverage_pb/compressed/", ind, ".txt"))){
    
    # read in the pac bio depth file
    pb_df <- read_delim(pb_keep_file, delim = " ", col_types = c("fdd"), col_names = FALSE)
    pb_df <- setNames(pb_df, c("CHROM", "POS", "depth"))

    # compress to 1kb scale
    pb_df <- pb_df %>%
      mutate(pos1 = (floor(POS/1000)) * 1000,
             pos2 = ((floor(POS/1000) + 1) * 1000)-1) %>%
      group_by(CHROM, pos1, pos2) %>%
      summarise(sum_depth = sum(depth, na.rm = TRUE), n_sites = n()) %>%
      mutate(Indiv = ind) %>%
      select(Indiv, everything())
    
    write.table(pb_df, file = paste0("analysis_sv/data/coverage_pb/compressed/", ind, ".txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
    
  }
  
  

}

pb_keep_files <- list.files("analysis_sv/data/coverage_pb/tmp", full.names = TRUE)
lapply(pb_keep_files, compress_depth_pb)


###################################################
# illumina
###################################################

compress_depth_illumina <- function(il_keep_file){
  
  message(il_keep_file)
  
  coords <- il_keep_file %>% gsub(".*/", "", .) %>% gsub("_keep.txt.gz", "", .) 
  
  if(!file.exists(paste0("analysis_sv/data/coverage_illumina/compressed/", coords, ".txt"))){
  
  # read in the pac bio depth file
  il_df <- read_delim(il_keep_file, delim = " ", col_names = c("CHROM", "POS", "Indiv", "depth"), col_types = c("cdcd"))
  
  # compress to 1kb scale
  il_df <- il_df %>%
    mutate(pos1 = (floor(POS/1000)) * 1000,
           pos2 = ((floor(POS/1000) + 1) * 1000)-1) %>%
    group_by(Indiv, CHROM, pos1, pos2) %>%
    summarise(sum_depth = sum(depth, na.rm = TRUE), n_sites = n()) %>%
    select(Indiv, everything())
  
  write.table(il_df, file = paste0("analysis_sv/data/coverage_illumina/compressed/", coords, ".txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
}

il_keep_files <- list.files("analysis_sv/data/coverage_illumina",  pattern = ".txt", full.names = TRUE)

lapply(il_keep_files, compress_depth_illumina)
mclapply(il_keep_files, compress_depth_illumina, mc.cores = 2)

###################################################
# combine depth measures
###################################################

# illumina
il_df <- list.files("analysis_sv/data/coverage_illumina/compressed", full.names = TRUE) %>%
  lapply(read.table, h = T,  colClasses = c("factor", "factor", "numeric", "numeric", "numeric", "numeric")) %>%
  bind_rows() %>%
  mutate(method = "illumina")

# illumina  sample names
ind_names <- read.table("meta/fasta_samplelist_recode.txt", header = FALSE) 
names(ind_names) <- c("Indiv", "ind2")
ind_names$Indiv <- gsub("S", "", ind_names$Indiv) %>% as.numeric

il_df <- il_df %>% 
  mutate(Indiv = gsub("S", "", Indiv) %>% as.numeric) %>%
  left_join(ind_names) %>%
  select(-Indiv) %>%
  rename(Indiv = ind2) %>%
  select(Indiv, everything())

pb_df <- list.files("analysis_sv/data/coverage_pb/compressed", full.names = TRUE)%>%
  lapply(read.table, h = T) %>%
  bind_rows() %>%
  mutate(method = "pacbio")

# pb sample names
ind_names <- read.table("analysis_sv/meta/pac_bio_sample_ids.txt", header = TRUE) 
names(ind_names) <- c("Indiv", "ind2")

pb_df <- pb_df %>% 
  left_join(ind_names) %>%
  select(-Indiv) %>%
  rename(Indiv = ind2) %>%
  select(Indiv, everything())

depth_df <- bind_rows(pb_df, il_df)


recomb_join_df <- list()
#pos_df <- depth_df %>%
#  select(CHROM, pos1) %>%
#  distinct %>%
#  mutate(pos1_ref = NA, pos2_ref = NA, POS_schf1 = NA, POS_schf2 = NA)

for(i in 1:nrow(recomb_intervals)){
  
  
  recomb_sub <- recomb_intervals[i,]
  pos1 <- seq(from = floor(recomb_sub$pos1_ref/1000), to = floor(recomb_sub$pos2_ref/1000), by = 1)*1000
  
  recomb_join_df[[i]] <- data.frame(CHROM = recomb_sub$CHROM_ref, pos1, 
                                    pos1_ref = recomb_sub$pos1_ref, pos2_ref = recomb_sub$pos2_ref, 
                                    POS_schf1 = recomb_sub$POS_schf1, POS_schf2 = recomb_sub$POS_schf2)
  
  
}

recomb_join_df <- bind_rows(recomb_join_df)

depth_df <- depth_df %>%
  left_join(recomb_join_df) %>% 
  filter(pos1 >= pos1_ref, pos2 <= pos2_ref)

saveRDS(depth_df, "analysis_sv/data/depth.rds")

