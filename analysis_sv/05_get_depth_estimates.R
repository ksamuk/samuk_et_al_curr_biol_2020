# extract per-window estimates of sequencing/mapping depth
# for both the PacBio and Illumina data
# KMS Aug 2019

library("tidyverse")
library("R.utils")
library("vcfR")

# the list of recombination intervals from shaeffer et al.
recomb_intervals <- read.table("analysis_sv/meta/recom_intervals_ref_schf.txt", h = T)

# prepare bed file

recomb_intervals %>%
  rowwise %>%
  mutate(pos1_ref = min(POS_ref1, POS_ref2),
         pos2_ref = max(POS_ref1, POS_ref2)) %>%
  mutate(POS_schf1 = pos1_schf, POS_schf2 = pos2_schf) %>%
  select(CHROM_ref, pos1_ref, pos2_ref) %>% 
  ungroup %>%
  write.table("analysis_sv/meta/coverage_recomb_intervals.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

########################################################################
# per window depth estimates: PacBio data
########################################################################

# get list of chrom:pos pairs to extract from file

head(recomb_intervals)

expand_interval <- function(row_num){
  
  int_sub <- recomb_intervals[row_num,]
  
  POS <- int_sub$POS_ref1:int_sub$POS_ref2
  CHROM <- int_sub$CHROM_ref
  
  data.frame(CHROM, POS)
}

interval_list <- lapply(1:nrow(recomb_intervals), expand_interval) %>% bind_rows

# precomputed depth estimates with bedtools, see genomics pipeline

per_window_depth <- function(pb_depth_file){
  
  message(pb_depth_file)
  
  file_number <- gsub(".*/", "", pb_depth_file) %>% gsub(".pacbio.coverage.txt.gz", "", .)
  
  temp_file <- "analysis_sv/data/coverage_pb/tmp/depth_file.tmp"
  file.remove(temp_file)
  
  gunzip(pb_depth_file, destname = temp_file, remove = FALSE, overwrite = TRUE)
  
  file.remove(paste0("analysis_sv/data/coverage_pb/tmp/",file_number,"_keep.txt"))
  
  con = file(temp_file, "r")
  
  repeat{
    
    cat(".")
    
    depth_df <- read.table(con, nrows = 2000000, colClasses = c("character", "numeric", "numeric"), header = FALSE)
    names(depth_df) <- c("CHROM", "POS", "depth")
    
    keep_list <- list()
    
    chrom_list <- unique(interval_list$CHROM) %>% as.character
    chrom_list <- chrom_list[chrom_list %in% unique(depth_df$CHROM)]
    
    cat(paste0(chrom_list, "..."))
    
    if(length(chrom_list) > 0){
      
      for (i in 1:length(chrom_list)){
        
        chrom_targ <- chrom_list[i]
        
        pos_targ <- interval_list %>% filter(CHROM == chrom_targ) %>% pull(POS)
        
        keep_list[[i]] <- depth_df %>%
          filter(CHROM == chrom_targ) %>%
          filter(POS %in% pos_targ) 
        
      }
      
      keep_list <- bind_rows(keep_list)
      
      write.table(keep_list, paste0("analysis_sv/data/coverage_pb/tmp/", file_number, "_keep.txt"), append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      
    }
    
    if (nrow(depth_df) != 2000000){
      print('Processed all files!')
      break}
  }
  
  close(con)
  file.remove("analysis_sv/data/coverage_pb/tmp/depth_file.tmp")
  
}

pb_depth_files <- list.files("analysis_sv/data/coverage_pb", full.names = TRUE)
pb_depth_files <- pb_depth_files[20:21]

#pb_depth_file <- pb_depth_files[[1]]
#per_window_depth(pb_depth_files[[1]])

lapply(pb_depth_files, per_window_depth)

########################################################################
# per window depth estimates: illumina data
########################################################################


compute_depth_illumina <- function(il_depth_file){
  
  coordinates <- gsub(".*/", "", il_depth_file) %>% gsub("_recomb.*", "", .)
  
  if(file.exists(paste0("analysis_sv/data/coverage_illumina/", coordinates, "_keep.txt.gz"))){
    
    cat("file exists, skipping...")
    
  } else{
    
    il_df <- read.vcfR(il_depth_file, check_keys = FALSE, checkFile = FALSE) %>%
      vcfR2tidy(single_frame = TRUE, info_fields = c("CHROM", "POS"), format_fields = "DP", alleles = FALSE, format_types = c(DP = "i"))
    
    il_df <- il_df$dat %>%
      select(CHROM, POS, Indiv, gt_DP)
    
    keep_list <- list()
    
    chrom_list <- unique(interval_list$CHROM) %>% as.character
    chrom_list <- chrom_list[chrom_list %in% unique(il_df$CHROM)]
    
    cat(paste0(il_depth_file, "..."))
    
    if(length(chrom_list) > 0){
      
      for (i in 1:length(chrom_list)){
        
        chrom_targ <- chrom_list[i]
        
        pos_targ <- interval_list %>% filter(CHROM == chrom_targ) %>% pull(POS)
        
        keep_list[[i]] <- il_df %>%
          filter(CHROM == chrom_targ) %>%
          filter(POS %in% pos_targ) 
        
      }
    }
    keep_list <- bind_rows(keep_list)
    
    write.table(keep_list, gzfile(paste0("analysis_sv/data/coverage_illumina/", coordinates, "_keep.txt.gz")), append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    
  }
  
  
}

il_depth_files <- list.files("analysis_sv/data/homology/illumina", full.names = TRUE)
il_depth_files <- il_depth_files[-c(1:7)]

mclapply(il_depth_files, compute_depth_illumina, mc.cores = 2)
