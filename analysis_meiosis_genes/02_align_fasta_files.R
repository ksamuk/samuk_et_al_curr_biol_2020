# align recombination gene fasta files using MAFFT
# KMS Jan 2019

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biostrings", version = "3.8")

options(device = "quartz")

library("tidyverse")
library("seqinr")


# the anderston et al genes and annotations
and_df <- read.table("meta/anderson_gene_locations.txt", h = F)
names(and_df) <- c("chr", "start_pos_pad", "end_pos_pad", "header_string", "dmel_gene_name", "cds_seq", "strand")
and_df$cds_id <- gsub(".*_", "", and_df$header_string)

# list the files
fasta_files <- list.files("data/fasta_recomb", full.names = TRUE, pattern = ".fasta")
dir.create("data/fasta_recomb/aligned")

align_with_mafft <- function(fasta_file){
  
  file_slug <- gsub(".*/", "", fasta_file)
  mafft_flags <- " --thread 4 --maxiterate 1000 "
  
  # build a system call to mafft for the fasta file
  sys_call <- paste0("mafft", mafft_flags, fasta_file, " > data/fasta_recomb/aligned/aln_", file_slug)
  system(sys_call)
  
}

# align all the genes
lapply(fasta_files, align_with_mafft)

# some inconsistent behaviour with mafft like this --
# had to rerun some files manually
# align_with_mafft("data/fasta_recomb/recomb_genes_cav_FBpp0284020_all.fasta")

aligned_fastas <- list.files("data/fasta_recomb/aligned", full.names = TRUE)

# function to parse the aligned fasta files

parse_fasta_file <- function(fasta_file){
  
  message(fasta_file)
  
  fasta_tmp <- read.fasta(fasta_file)
  
  fasta_names <- names(fasta_tmp)
  
  out_list <- lapply(fasta_names, function(x) data.frame(info_string = x, geno = fasta_tmp[x] %>% 
                                               as.character() %>% 
                                               gsub("\\*|\\(|\\)|\\\"|\\,|\\n| ", "", .) %>% 
                                               strsplit(split = "")  %>% unlist))
  bind_rows(out_list)
  
}

fasta_df <- lapply(aligned_fastas, parse_fasta_file) %>%
  bind_rows()

fasta_df$geno %>% unique

# parse the gene name/sample id
fasta_df <- fasta_df %>%
  mutate(sample_id = gsub("_.*", "", info_string)) %>%
  mutate(gene_name = gsub("_[A-Za-z0-9]*$", "", info_string) %>% gsub(".*_", "", .)) %>%
  mutate(cds_id = gsub(".*_", "", info_string)) %>%
  group_by(gene_name, cds_id, sample_id) %>%
  mutate(rel_pos = 1:length(geno)) %>%
  ungroup %>%
  select(gene_name, cds_id, sample_id, rel_pos, geno)

# join in the strand data
fasta_df <- fasta_df %>%
  left_join(and_df %>% select(cds_id, strand))

fasta_df$gene_name %>% unique

# plot the haplotypes to check alignments
fasta_df %>%
  mutate(gene_name_full = paste0(gene_name, "_", cds_id)) %>%
  mutate(geno_color = case_when(geno == "a" ~ "a",
                                geno == "t" ~ "t",
                                geno == "c" ~ "c",
                                geno == "g" ~ "g",
                                geno == "-" ~ "-",
                                TRUE ~ "z_poly"
                                )) %>%
  ggplot(aes(x = rel_pos, y = as.factor(sample_id), fill = geno_color))+
  geom_raster(color = NA)+
  facet_wrap(~gene_name_full, scale = "free_x")+
  scale_fill_viridis_d()

# revcomp the -ve strand sequneces
# (CDS's were flipped for alignment)
# only the genotypes, not the (arbitrary) relative position

fasta_df <- fasta_df %>%
  mutate(gene_name_full = paste0(gene_name, "_", cds_id)) %>%
  group_by(gene_name_full, sample_id) %>%
  do(rel_pos = .$rel_pos,
     geno = ifelse(.$strand == "-", rev(.$geno) %>% chartr("ATGC", "TACG", .), .$geno)) %>% 
  ungroup %>%
  unnest

# flip between two plots to check revcomp was successful
fasta_df %>%  
  mutate(geno_color = case_when(geno == "a" ~ "a",
                                geno == "t" ~ "t",
                                geno == "c" ~ "c",
                                geno == "g" ~ "g",
                                geno == "-" ~ "-",
                                TRUE ~ "z_poly")) %>%
  ggplot(aes(x = rel_pos, y = as.factor(sample_id), fill = geno_color))+
  geom_raster(color = NA)+
  facet_wrap(~gene_name_full, scale = "free_x")+
  scale_fill_viridis_d()
  
# build the cds masks

cds_mask <- fasta_df %>%
  filter(sample_id == "CDS", geno != "-") %>%
  mutate(in_cds = TRUE) %>%
  select(-sample_id) %>%
  arrange(gene_name_full, rel_pos)

# validate the start codons in the CDSs
# find the first occurance of a start codon
cds_start_pos <- cds_mask %>%
  group_by(gene_name_full) %>%
  mutate(start_pos = min(rel_pos)) %>% 
  mutate(atg_motif_t = geno == "t" & 
           lag(geno) == "a" & 
           lead(geno) == "g" ) %>%
  mutate(atg_motif = 
           geno == "g" & lag(atg_motif_t) | 
           geno == "a" & lead(atg_motif_t) |
           atg_motif_t) %>%
  select(-atg_motif_t) %>%
  filter(!is.na(atg_motif)) %>%
  mutate(start_pos = min(rel_pos[atg_motif])) %>%
  select(gene_name_full, start_pos) %>%
  ungroup %>%
  distinct
  
# clip any bases before the first start codon
cds_mask <- cds_mask %>%
  left_join(cds_start_pos) %>%
  filter(rel_pos >= start_pos) %>%
  arrange(gene_name_full, rel_pos) %>%
  select(-start_pos, -geno)
  
# join mask back into sequence data
cds_df <- fasta_df %>%
  left_join(cds_mask) %>%
  filter(in_cds)

# visualize the atg motifs in the cds
# there should be at least one at the start of the sequence
cds_df %>%  
  mutate(atg_motif_t = geno == "t" & 
           lag(geno) == "a" & 
           lead(geno) == "g" ) %>%
  mutate(atg_motif = 
           geno == "g" & lag(atg_motif_t) | 
           geno == "a" & lead(atg_motif_t) |
           atg_motif_t) %>%
  group_by(gene_name_full, sample_id) %>%
  mutate(rel_pos_cds = 1:length(geno)) %>%
  ungroup %>%
  mutate(geno_color = case_when(geno == "a" ~ "a",
                                geno == "t" ~ "t",
                                geno == "c" ~ "c",
                                geno == "g" ~ "g",
                                geno == "-" ~ "-",
                                TRUE ~ "z_poly")) %>%
  ggplot(aes(x = rel_pos_cds, y = as.factor(sample_id), fill = atg_motif))+
  geom_raster()+
  facet_wrap(~gene_name_full, scale = "free_x")+
  scale_fill_manual(values = c("white", "black"))


# functions for converting polymorphisms to one of two types
build_hap1 <- function(x){
  
  if(x %in% c("a", "t", "c", "g", "-")){
    
    return(x)
    
  } else if(x %in% c("u","r","y","m","k","s","w","b","d","h","v","n")){
    
    seqinr::amb(x) %>% unlist %>% .[1]
    
  } else{
    
    message(paste0("INVALID BASE: ", x))
    
  }
  
}

build_hap2 <- function(x){
  
  if(x %in% c("a", "t", "c", "g", "-")){
    
    return(x)
    
  } else if(x %in% c("u","r","y","m","k","s","w","b","d","h","v","n")){
    
    seqinr::amb(x) %>% unlist %>% .[2]
    
  } else{
    
    message(paste0("INVALID BASE: ", x))
  
  }
  
}

# make two haplotypes for each sample, to account for ambiguities
build_alt_haps <- function(sample_name){
  
  message(sample_name)
  
  # subset the cds df
  cds_sub <- cds_df %>%
    filter(sample_id == sample_name)
  
  # check for polymorphisms
  if (any(!(cds_sub$geno %in% c("a", "t", "c", "g", "-")))){
    
    # build two alternative haplotypes
    cds_sub$hap1 <- cds_sub$geno %>% lapply(build_hap1) %>% unlist
    cds_sub$hap2 <- cds_sub$geno %>% lapply(build_hap2) %>% unlist
    
  } else{
    
    cds_sub$hap1 <- cds_sub$geno
    cds_sub$hap2 <- cds_sub$geno
    
  }
  
  return(cds_sub)
  
}

alt_haps <- lapply(unique(cds_df$sample_id), build_alt_haps) %>%
  bind_rows

# write out alt haps

saveRDS(alt_haps, "data/fasta_recomb/summaries/cds_alt_haps.rds")

# output the clipped cds to a (new) fasta file
cds_out <- alt_haps %>%
  group_by(gene_name_full, sample_id) %>%
  summarize(fasta_seq1 = paste0(hap1, collapse = ""),
            fasta_seq2 = paste0(hap2, collapse = ""),
            fasta_header = paste0(">", gene_name_full[1], "_sample=", sample_id[1])) %>%
  mutate(fasta_lines = paste0(fasta_header,"_hap1", "\n", fasta_seq1, "\n",
                               fasta_header,"_hap2", "\n", fasta_seq2))

write_cds_fastas <- function(gene){
  
  lines_sub <- cds_out %>% filter(gene_name_full == gene) %>% pull(fasta_lines)
  write_lines(lines_sub, paste0("data/fasta_recomb/cds_aligned/", gene, "_","cds_aligned.fasta"))
  
}

lapply(unique(cds_out$gene_name_full), write_cds_fastas)

