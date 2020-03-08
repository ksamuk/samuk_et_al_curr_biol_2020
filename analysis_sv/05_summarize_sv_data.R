# summarize sv data for each window
# KMS 2019

library("tidyverse")
library("parallel")
#library("parallelsugar")

#options(device = "quartz")

######################################
# read in tidy SV data
######################################

# sv data
sv_df <- readRDS("analysis_sv/data/sv_joined.rds")

# coverted positions
sv_remap <- readRDS("analysis_sv/meta/remap.rds") %>%
  select(CHROM_orig, POS_orig, CHROM_new, POS_new)

# apply conversion to sv data
# handling pos1 and pos2 separately

pos1_sv <- sv_df %>%
  select(ind, chrom, pos1) %>%
  rename(CHROM_orig = chrom, POS_orig = pos1) %>%
  left_join(sv_remap)%>%
  rename(pos1_schf = POS_new)

pos2_sv <- sv_df %>%
  select(ind, chrom, pos2) %>%
  rename(CHROM_orig = chrom, POS_orig = pos2) %>%
  left_join(sv_remap) %>%
  rename(pos2_schf = POS_new)
  

sv_df <- cbind(sv_df, pos1_sv$pos1_schf, pos2_sv$pos2_schf) %>%
  rename(pos1_schf = `pos1_sv$pos1_schf`) %>%
  rename(pos2_schf = `pos2_sv$pos2_schf`)

sv_df <- sv_df %>%
  select(ind, chrom, pos1, pos2, pos1_schf, pos2_schf, everything())

######################################
# read in recombination data
######################################

# this is from count_crossover_events_finescale
recomb_df <- readRDS("data/crossovers_fine_scale.rds")

#######################################################
# build a list of recombination intervals for each line
#######################################################

interval_df <- recomb_df %>%
  ungroup %>%
  select(Indiv, CHROM_schf, POS_schf1, POS_schf2) %>%
  mutate(line = gsub("_[0-9]{1,2}_.{3}$", "",  Indiv) %>% 
           gsub("_", "", .) %>%
           gsub("AFC", "A", .) %>%
           gsub("MC", "M", .)) %>%
  select(line, CHROM_schf, POS_schf1, POS_schf2) %>%
  distinct
  
############################################################
# for each line and interval, 
# count the number of SVs of each type (?) in the interval
############################################################

indiv_list <- interval_df$line %>% unique
indiv <- indiv_list[1]

indiv <- "A14"

split_collapse_genotype <- function(x){
  
  strsplit(x, split = "/") %>% paste0(collapse = "") %>% as.character
  
}

tabulate_svs_per_interval <- function(indiv, verbose = TRUE){
  
  message(indiv)
  
  # subset recomb data for the target individual
  interval_df <- interval_df %>%
    filter(line == indiv)
  
  # get the line info from the individual's ID
  indiv_line <- indiv %>% 
    gsub("_[0-9]{1,2}_.{3}$", "", .) %>% 
    gsub("_", "", .) %>% 
    gsub("AFC", "A", .) %>%
    gsub("MC", "M", .)
  
  # subset the SV data for the individual's line
  sv_sub <- sv_df %>%
    filter(ind == indiv_line)
  
  # also grab the MV225 data
  # we want to compare SVs to this
  # (our MV225 has some non-ref variation)
  sv_sub_mv <- sv_df %>%
    filter(ind == "MV2-25")
  

  # count the number of each type of SV in each recombination interval
  
  #i <- 2
  chr_list <- vector("list", length(unique(interval_df$CHROM_schf)))
  
  for (i in 1:length(unique(interval_df$CHROM_schf))) {
    
    recomb_sub_chrom <- interval_df %>%
      filter(CHROM_schf == unique(interval_df$CHROM_schf)[i])
    
    interval_list <- vector("list", length(recomb_sub_chrom$POS_schf1))
    #j <- 10
    for (j in 1:length(recomb_sub_chrom$POS_schf1)){
      
      if(verbose){
        
        cat(paste0(unique(interval_df$CHROM_schf)[i], "\n"))
        cat(paste0(j, "\n"))
        
      }
      
      # define the start and end of the recombination interval
      # in some cases the chromosome is in reversed order, so need 
      # to rearrange intervals to match the SV data
      pos2 <- recomb_sub_chrom %>% 
        filter(POS_schf1 == recomb_sub_chrom$POS_schf1[j]) %>% 
        pull(POS_schf2) -1
      
      interval <- sort(c(recomb_sub_chrom$POS_schf1[j], pos2))
      
      interval_start <- interval[1]
      interval_end <- interval[2]
      
      # what i want here
      # chrom pos1 pos2 crossover no_TE_diffs len_TE_diffs no_INS_diffs len_TE_diffs
      
      # find the svs in the target interval
      svs_in_interval <- sv_sub %>%
        #filter(!(genotype %in% c("0/0","0/1"))) %>%
        filter(pos1_schf >= interval_start & pos2_schf <= interval_end) %>%
        select(ind, pos1_schf, pos2_schf, method, sv_type, sv_type_line, sv_length, genotype)
        
      
      if(nrow(svs_in_interval) == 0){
        
        svs_in_interval <- data.frame(ind = indiv, pos1_schf = NA, pos2_schf = NA, 
                                         method = NA, sv_type = NA, sv_type_line = NA, sv_length = 0, genotype = NA,
                                         length_class = NA, is_te = NA)
      }
      
      mv_svs_in_interval <- sv_sub_mv %>%
        #filter(!(genotype %in% c("0/0","0/1"))) %>%
        filter(pos1_schf >= interval_start & pos2_schf <= interval_end) %>%
        select(ind, pos1_schf, pos2_schf, method, sv_type, sv_type_line, sv_length, genotype) %>%
        rename(mv_genotype = genotype, sv_type_mv = sv_type_line, sv_length_mv = sv_length) %>%
        select(-ind)
      

      if(nrow(mv_svs_in_interval) == 0){
        
        mv_svs_in_interval <- data.frame(pos1_schf = 0, pos2_schf = 0, method = NA, 
                                         sv_type = NA, sv_type_mv = NA, sv_length_mv = NA, mv_genotype = NA)
      }
      
      sv_join <- full_join(svs_in_interval, mv_svs_in_interval, by = c("pos1_schf", "pos2_schf", "method", "sv_type"))
      
      sv_join_simple <- sv_join %>%
        filter(sv_type != "Complex") %>%
        mutate(sv_length = ifelse(sv_type == "TE", pos2_schf - pos1_schf, sv_length))%>%
        mutate(sv_length = ifelse(method == "gatk" & sv_type == "DEL", -(sv_length), sv_length))%>%
        mutate(sv_length = as.numeric(sv_length)) %>%
        mutate(length_class = ifelse(sv_length > 0, "ins", "del")) %>%
        mutate(is_te = (sv_type == "TE")) %>%
        rowwise %>%
        mutate(len_diff = case_when(
          genotype == mv_genotype ~ 0,
          mv_genotype != genotype ~ sv_length,
          TRUE ~ NA_real_)) %>%
        mutate(len_diff = ifelse(any(c(genotype, mv_genotype) %in% c("0/1", "0/2", "1/2")), len_diff/2, len_diff)) %>%
        rowwise %>%
        mutate(diff_count = ifelse(!is.na(genotype) & !is.na(mv_genotype), 
                                   adist(split_collapse_genotype(genotype), split_collapse_genotype(mv_genotype)), NA)) %>%
        ungroup
      
      sv_join_complex <- sv_join %>%
        filter(sv_type == "Complex") %>% 
        mutate(sv_length = as.numeric(sv_length)) %>%
        mutate(length_class = ifelse(sv_length > 0, "ins", "del")) %>%
        mutate(is_te = (sv_type == "TE")) %>%
        rowwise %>%
        mutate(len_diff = case_when(
          genotype == mv_genotype ~ 0,
          sv_type_line == "DEL" & sv_type_mv == "REF" ~ -sv_length,
          sv_type_line == "REF" & sv_type_mv == "DEL" ~ sv_length,
          sv_type_line == "INS" & sv_type_mv == "REF" ~ sv_length,
          sv_type_line == "REF" & sv_type_mv == "INS" ~ -sv_length,
          sv_type_line == "INS" & sv_type_mv == "DEL" ~ 2*sv_length,
          sv_type_line == "DEL" & sv_type_mv == "INS" ~ -(2*sv_length),
          mv_genotype == "1/1" & genotype == "0/0" ~ sv_length,
          mv_genotype == "0/0" & genotype == "1/1" ~ sv_length,
          TRUE ~ NA_real_)) %>%
        mutate(len_diff = ifelse(any(c(genotype, mv_genotype) %in% c("0/1", "0/2", "1/2")), len_diff/2, len_diff)) %>%
        mutate(diff_count = ifelse(!is.na(genotype) & !is.na(mv_genotype), 
                                   adist(split_collapse_genotype(genotype), split_collapse_genotype(mv_genotype)), NA)) %>%
        ungroup
      
      sv_join <- bind_rows(sv_join_simple, sv_join_complex)
      
      sv_summary <- sv_join %>%
        group_by(method, sv_type) %>%
        summarise(n_genotyped = sum(!is.na(diff_count)),
                  diff_count  = sum(diff_count, na.rm = TRUE),
                  len_diff  = sum(len_diff, na.rm = TRUE)) 
      

      
      sv_summary <- data.frame(Indiv = indiv, CHROM_schf = unique(interval_df$CHROM_schf)[i], 
                               pos1 = interval_start, pos2 = interval_end, sv_summary)
      
      interval_list[[j]] <- sv_summary
      
    }
    
    chr_list[[i]] <- bind_rows(interval_list) %>%
      filter(!is.na(method))
    
  }
  
  bind_rows(chr_list) %>%
    distinct()
  
}

sv_summary_df <- lapply(indiv_list, tabulate_svs_per_interval)
  
#sv_summary_df <- mclapply(indiv_list, tabulate_svs_per_interval, mc.cores = 6)

sv_summary_df <- bind_rows(sv_summary_df)

saveRDS(sv_summary_df, "analysis_sv/data/sv_window_summary.rds")


