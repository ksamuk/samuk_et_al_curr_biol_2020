# get annotation information for meiosis-associated genes
# KMS Jan 2019

#if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#BiocManager::install("InterMineR", version = "3.8")

#devtools::install_github("ropensci/rflybase")

library("InterMineR")
library("tidyverse")
library("seqinr")

###########################################################
# anderson et al. 2009 gene analysis
###########################################################

# read in anderson et al gene ids
and_gene_ids <- read.csv("meta/anderson_et_al_gene_ids.csv")

# from Hunter et al. 2016
hunter_gene_ids <- c("CG10864", "CG33970", "Eip75B", "lola", "Ptp61F")

# load FlyMine and HumanMine 
im.fly <- initInterMine(listMines()["FlyMine"])

# load templates
templates.fly <- getTemplates(im.fly)

# load data models
model.fly <- getModel(im.fly)

# build a query function for dmel -> dpse orthologues 

query_gene_dros_ortho <- getTemplateQuery(
  im = im.fly, 
  name = "Gene_DrosophilaOrthologues"
)

query_gene_dros_ortho$where[[1]]$value <- "Drosophila pseudoobscura"

query_for_pseudoobscura_orthologs <- function(gene_symbol){
  
  query_tmp <- query_gene_dros_ortho
  query_tmp$where[[2]]$value <- gene_symbol
  
  runQuery(im.fly, query_tmp)
  
}

# get the dmel -> dpse orthologs for the anderson gene list
dpse_ortho_df <- lapply(as.character(and_gene_ids$Gene), query_for_pseudoobscura_orthologs) %>%
  bind_rows()

dpse_ortho_df2 <- lapply(as.character(hunter_gene_ids), query_for_pseudoobscura_orthologs) %>%
  bind_rows()

# combine hits
dpse_ortho_df <- bind_rows(dpse_ortho_df, dpse_ortho_df2)

# read in the dpse flybase gff file

# unzip data/dpse-all-r3.04.gff.gz into data/ and then ran:
# grep 'FlyBase' data/dpse-all-r3.04.gff | grep -v '##' > data/flybase_genes.txt

# fun: column type imputation failed and needed to be manually specified 
dpse_gff <- read_delim("data/flybase_genes.txt", delim = "\t", col_types = "cccnncccc", col_names = FALSE)
names(dpse_gff) <- c("chr", "db", "feature_type", "start_pos", "end_pos", "X6", "strand", "X8", "data_string")

# basic parse of the field string for genes
# "id", "gene_name", "alias", 
dpse_genes <- dpse_gff %>%
  filter(feature_type == "gene") %>%
  mutate(data_string = gsub(";Alias.*|;fullname.*", "", data_string))%>%
  mutate(data_string = gsub("ID=|Name=", "", data_string)) %>%
  separate(data_string, into = c("dpse_gene_id", "dpse_gene_name"), sep = ";")

# same for CDS/mMRA
# TBD
dpse_cds <- dpse_gff %>%
  filter(feature_type == "CDS")
#  mutate(data_string = gsub(".*=", "", data_string)) %>%
#  rename(dpse_gene_id = data_string)

target_genes <- dpse_genes %>%
  filter(dpse_gene_id %in% dpse_ortho_df$Gene.homologues.homologue.primaryIdentifier) %>%
  left_join(dpse_ortho_df %>% 
              rename(dpse_gene_id = Gene.homologues.homologue.primaryIdentifier) %>% 
              rename(dmel_gene_id = Gene.primaryIdentifier) %>% 
              rename(dmel_gene_name = Gene.symbol) %>% 
              select(dpse_gene_id, dmel_gene_name, dmel_gene_id))

control_genes <- dpse_genes %>%
  filter(!(dpse_gene_id %in% dpse_ortho_df$Gene.homologues.homologue.primaryIdentifier)) %>%
  sample_n(40) %>% 
  left_join(dpse_ortho_df %>% 
              rename(dpse_gene_id = Gene.homologues.homologue.primaryIdentifier) %>% 
              rename(dmel_gene_id = Gene.primaryIdentifier) %>% 
              rename(dmel_gene_name = Gene.symbol) %>% 
              select(dpse_gene_id, dmel_gene_name, dmel_gene_id))
  


# pull down CDS sequences via URL interface to FlyBase
# e.g:
# http://flybase.org/download/sequence/batch/?ids=FBgn0000490,FBgn0013765&type=CDS&output=file
# (had to build the link and paste it in the browser to get the fasta file)

fb_url_prefix <- "http://flybase.org/download/sequence/batch/?ids="
fb_url_suffix <- "&type=CDS&output=file"
fb_url <- paste0(fb_url_prefix, paste(target_genes$dpse_gene_id, collapse = ","), fb_url_suffix)

# pull out fasta sequences for each gene id
fb_cds <- read.fasta("data/flybase_cds.fasta", forceDNAtolower = FALSE)

# functio to simply cds sequences into a data frame
collapse_cds <- function(x){
  
  annot <- attr(x, "Annot") 
  dpse_gene_id <- annot %>% gsub(".*parent=", "", .) %>% gsub(",FBtr.*", "", .)
  dpse_cds_id <- annot %>% gsub(".*dbxref=.*FlyBase:", "", .) %>% gsub(",.*", "", .)
  
  cds_seq <- x %>% as.character %>% paste0(collapse = "")
  
  data.frame(dpse_gene_id, dpse_cds_id, cds_seq)
  
}

cds_df <- lapply(fb_cds, collapse_cds) %>% bind_rows

# join in cds sequences for each gene
target_genes <- target_genes %>%
  left_join(cds_df)

control_genes <- control_genes %>%
  left_join(cds_df)

# format the target genes df for output
# create bounds, a header, etc.
# also included 1000bp up and downstream of gene

target_genes <- target_genes %>%
  mutate(start_pos_pad = start_pos - 1000) %>%
  mutate(end_pos_pad = end_pos + 1000) %>%
  mutate(dmel_gene_name = paste0(dmel_gene_name, "_", dpse_cds_id)) %>%
  select(chr, start_pos_pad, end_pos_pad, dpse_gene_id, dmel_gene_id, dmel_gene_name, cds_seq, strand) %>%
  mutate(header_string = paste0(">chr_", chr, "_", start_pos_pad, "-", end_pos_pad, 
                                "_dpse-", dpse_gene_id, "_dmel-", dmel_gene_id, "_", dmel_gene_name)) %>%
  mutate(cds_seq = ifelse(strand == "-", stringi::stri_reverse(cds_seq) %>% chartr("ATGC", "TACG", .), cds_seq)) %>%
  select(chr, start_pos_pad, end_pos_pad, header_string, dmel_gene_name, cds_seq, strand) %>%
  distinct

control_genes <- control_genes %>%
  mutate(start_pos_pad = start_pos - 1000) %>%
  mutate(end_pos_pad = end_pos + 1000) %>%
  mutate(dmel_gene_name = paste0(dmel_gene_name, "_", dpse_cds_id)) %>%
  select(chr, start_pos_pad, end_pos_pad, dpse_gene_id, dmel_gene_id, dmel_gene_name, cds_seq, strand) %>%
  mutate(header_string = paste0(">chr_", chr, "_", start_pos_pad, "-", end_pos_pad, 
                                "_dpse-", dpse_gene_id, "_dmel-", dmel_gene_id, "_", dmel_gene_name)) %>%
  mutate(cds_seq = ifelse(strand == "-", stringi::stri_reverse(cds_seq) %>% chartr("ATGC", "TACG", .), cds_seq)) %>%
  select(chr, start_pos_pad, end_pos_pad, header_string, dmel_gene_name, cds_seq, strand) %>%
  distinct

write.table(target_genes, "meta/anderson_gene_locations.txt", 
            col.names = FALSE, quote = FALSE, row.names = FALSE)

col_names <- c("chr", "start_pos_pad", "end_pos_pad", "header_string", "dmel_gene_name", "cds_seq", "strand")

candidate_gene_suppmat <- read.table("meta/anderson_gene_locations.txt", col.names = col_names)


candidate_gene_suppmat %>%
  mutate(dpse_gene = gsub(".*dpse-|_dmel.*" ,"", header_string)) %>%
  mutate(dmel_gene = gsub(".*_dmel-|_.*" ,"", header_string)) %>%
  mutate(start_pos = start_pos_pad + 1000) %>%
  mutate(end_pos = end_pos_pad - 1000) %>%
  select(dmel_gene, dpse_gene, chr, start_pos, end_pos) %>%
  write.csv("figures/TableS3.csv", row.names = FALSE, quote = FALSE)

