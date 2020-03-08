# rblast

devtools::install_github("mhahsler/rBLAST")
library("rBLAST")
library("tidyverse")

# point to blast
# for mac os
#Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/Users/ksamuk/Dropbox (Duke Bio_Ea)/Projects/gt_seq_rev/analysis_cm_mb/bin/ncbi-blast-2.8.1/bin", 
#                        sep= .Platform$path.sep))
# check that blast exists!
Sys.which("blastn")
Sys.which("megablast")
Sys.which("makeblastdb")

##################################################################
# create a blast db for the miller reference
##################################################################

blast_help("makeblastdb")

## create a database for some example sequences
seq <- readDNAStringSet("analysis_cm_mb/data/Dpse.pass.minimap2.racon.x3.pilon.x3.fasta")

## 1. write the FASTA file to a temp dir
dir <- tempdir()
writeXStringSet(seq, filepath = file.path(dir, "miller_ref.fasta"))

## 2. make database
makeblastdb(file.path(dir, "miller_ref.fasta"), dbtype = "nucl")

## 3. open database
db <- blast(file.path(dir, "miller_ref.fasta"))
db

##################################################################
# blast the gt seq primer targets
##################################################################

gt_fa <- readDNAStringSet("analysis_cm_mb/data/pseudoobscura_target_regions_blast.fasta")

gt_hits <- predict(db, gt_fa)
View(gt_hits)


##################################################################
# blast the dpse genes
##################################################################

## 3. open database
gene_fa <- readDNAStringSet("analysis_cm_mb/data/dpse-all-gene-r3.04.fasta")

blast_gene <- function(x){
  
  cat(paste0(x, "\n"))
  
  predict(db, gene_fa[x,], 
          BLAST_args = "-reward 2 -penalty -3 -word_size 20 -gapopen 5 -gapextend 2 -num_threads 3")
  
  
}

# attempt to recapitulate megablast-ish settings for longer sequences
gene_hits <- lapply(1:length(gene_fa), blast_gene)
gene_hits_out <- bind_rows(gene_hits)

View(gene_hits)

##################################################################
# write the results to file
##################################################################

write.table(gt_hits, "analysis_cm_mb/out/gt_seq_miller_blast.txt", row.names = FALSE, quote = FALSE)
write.table(gene_hits_out, "analysis_cm_mb/out/dpse_genes_miller_blast.txt", row.names = FALSE, quote = FALSE)
