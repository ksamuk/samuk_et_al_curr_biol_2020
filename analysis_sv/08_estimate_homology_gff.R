# compute homology for recombination intervals
# via pacbio gff files ONLY
# KMS 2019

library("tidyverse")
library("rtracklayer")

# get the individual codes for the sample IDs
meta_df <- read.table("meta/fasta_samplelist_recode.txt", h = F)
names(meta_df) <- c("Indiv", "ID")

parse_gff <- function(ind_id){
  
  # get the individual id 
  id <- meta_df %>%
    filter(Indiv == paste0("S", ind_id)) %>%
    pull(ID)
  
  
  gff_df <- readGFF("analysis_sv/data/homology/pacbio_gff/1_variants.gff") %>%
    select(seqid, type, start, end, reference, variantSeq, coverage, confidence)
  
  gff_df %>%
    mutate(Indiv = id) %>%
    rename(CHROM = seqid, pos1 = start, pos2 = end) %>%
    select(Indiv, CHROM, pos1, pos2, type, everything())
  
}
