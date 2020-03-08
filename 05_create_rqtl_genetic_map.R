library("tidyverse")
library("qtl")
library("ASMap")

#options(device = "quartz")


# function for automatic estimation of a genetic map
cross_id <- "all"
data_dir <- "data/rqtl/raw"
missing_marker_cutoff = 0.05
missing_ind_cutoff = 0.05
marker_co_cutoff = 0.95
max_cross_overs = 15

create_genetic_map <- function(cross_id, data_dir,
                               missing_marker_cutoff = 0.05, 
                               missing_ind_cutoff = 0.05, 
                               marker_co_cutoff = 0.95,
                               max_cross_overs = 15){
  
  #############################################################
  # reading in raw data + filtration
  #############################################################
  
  geno_file <- paste0(cross_id, "_geno_dat.csv")
  pheno_file <- paste0(cross_id, "_pheno_dat.csv")
  
  message(paste0("Creating map for ", cross_id, "..."))
  
  # read in data
  # NOTE: investigate if pgm is being inferred properly/matters
  qtl_dat <- read.cross("csvs", dir = data_dir, genfile = geno_file, phefile = pheno_file,
                        genotypes = c("BB", "AB", "AA", "D", "C"), crosstype = "bc")
  
  # fill genotypes based on maximal marginal genotypes
  # vqtl_dat <- fill.geno(qtl_dat, method = "no_dbl_XO") 
  qtl_dat <- fill.geno(qtl_dat, method = "maxmarginal")
  
  # plot missing data
  #plotMissing(qtl_dat)
  
  # remove markers with evidence of seg dist
  # expected ratio for BC is 0.5
  gt <- geno.table(qtl_dat)
  gt$prop <- gt$AA/(gt$AB+gt$AA)
  
  gt <- data.frame(marker_id = row.names(gt), gt)
  
  #write.table(gt, file = "meta/rqtl_missingness.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  qtl_dat <- drop.markers(qtl_dat, row.names(gt[ (gt$prop < 0.4)|(gt$prop > 0.6) | is.na(gt$prop),]))
  
  # plot coverage / representation
  par(mfrow=c(1,2), las=1)
  plot(ntyped(qtl_dat), ylab="No. typed markers", main="No. genotypes by individual")
  plot(ntyped(qtl_dat, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")
  
  # drop sites in the nth quantiles of coverage
  nt.bymar <- ntyped(qtl_dat, "mar")
  n_typed_cutoff <- quantile(nt.bymar, missing_marker_cutoff) %>% as.numeric
  
  todrop <- names(nt.bymar[nt.bymar < n_typed_cutoff])
  qtl_dat <- drop.markers(qtl_dat, todrop)
  
  # drop individuals in the 0-25th quantiles of marker num
  n_markers_cutoff <- quantile(ntyped(qtl_dat), missing_ind_cutoff) %>% as.numeric
  qtl_dat <- subset(qtl_dat, ind = (ntyped(qtl_dat) > n_markers_cutoff))
  #qtl_dat <- fill.geno(qtl_dat, method = "no_dbl_XO") 
  
  # genotype classes for all markers
  # g <- pull.geno(qtl_dat)
  # gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
  # gfreq <- t(t(gfreq) / colSums(gfreq))
  # par(mfrow=c(1,3), las=1)
  # for(i in 1:3)
  #   plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))
  
  
 # check for switched alleles
  qtl_dat <- est.rf(qtl_dat)
  checkAlleles(qtl_dat, threshold = 3)
  
  # plot rf lod scores vs. recombination fractions
  # should look like a decay towards ~0.4-0.5 (not a u-shape)
  rf <- pull.rf(qtl_dat)
  lod <- pull.rf(qtl_dat, what="lod")
  plot(as.numeric(rf), as.numeric(lod), xlab = "Recombination fraction", ylab = "LOD score")
  
  #############################################################
  # map creation + marker ordering
  #############################################################
  
  # find the initial optimal marker orderwith ASmap
  qtl_reorder <- mstmap(qtl_dat, id = "id", bychr = FALSE, dist.fun = "kosambi", 
                        trace = FALSE, detectBadData = TRUE, mvest.bc = FALSE, p.value = 1e-15)
  qtl_reorder <- est.rf(qtl_reorder)
  plotRF(qtl_reorder, what = "rf")
  
  # plot the order map and rfs
  #plotRF(qtl_reorder, what = "rf")
  
  # find and drop problem individuals (very high CO rates)
  plot(countXO(qtl_reorder), ylab = "Number of crossovers")
  hist(countXO(qtl_reorder))
  qtl_reorder_clean <- subset(qtl_reorder, ind = (countXO(qtl_reorder) <= max_cross_overs))
  
  plot(countXO(qtl_reorder_clean), ylab = "Number of crossovers")
  hist(countXO(qtl_reorder_clean))
  
  # reorder with mst map
  qtl_reorder_clean <- mstmap(qtl_reorder_clean, id = "id", bychr = TRUE, dist.fun = "kosambi", 
                              trace = FALSE, detectBadData = TRUE, mvest.bc = TRUE)
  qtl_reorder_clean <- est.rf(qtl_reorder_clean)
  plotRF(qtl_reorder_clean, what = "rf")
  
  # find an drop markers with abnormally high (>95th percentile) CO rates
  # this doesn't have much of an effect probably because of the upstream filtering
  marker_co <- statMark(qtl_reorder_clean)$marker
  marker_co$marker_name <- row.names(marker_co)
  marker_co_cutoff <- quantile(marker_co$dxo, 0.95) %>% as.numeric
  bad_markers <-  marker_co %>% filter(dxo >= marker_co_cutoff) %>% pull(marker_name)
  qtl_reorder_clean <- drop.markers(qtl_reorder_clean, bad_markers)
  
  # re-estimate map (again)
  qtl_reorder_clean <- mstmap(qtl_reorder_clean, id = "id", bychr = TRUE, dist.fun = "kosambi", 
                              trace = FALSE, detectBadData = FALSE, mvest.bc = TRUE, p.value = 1e-10)
  qtl_reorder_clean <- est.rf(qtl_reorder_clean)
  plotRF(qtl_reorder_clean, what = "rf")
  
  rf_dat <- pull.rf(qtl_reorder_clean)
  # perform dropone to identify problematic markers
  # not included in final pipeline
  # difficult to justify imo
  # drop_df2 <- droponemarker(qtl_reorder_clean, map.function = "kosambi", maxit = 10)


  # rename chromosomes based on majority reference chromosome in marker names
  
  get_majority_chromosome_name <- function(chr_name){
    
    marker_list <- qtl_reorder_clean$geno[chr_name][[1]]$map %>% names
    marker_list <- marker_list %>% gsub(".*_chr_", "", .) %>% 
      gsub("_.*", "", .) %>% gsub("XR|XL", "X", .)
    
    data.frame(orig = chr_name, 
               actual = names(which.max(table(marker_list))),
               stringsAsFactors = FALSE)
    
  }
  
  chr_names <- names(qtl_reorder_clean$geno)
  chr_names <- lapply(chr_names, get_majority_chromosome_name) %>% bind_rows
  
  names(qtl_reorder_clean$geno) <- chr_names$actual
  
  countXO(qtl_reorder_clean, chr = "X") %>% hist
  
  # rename chromosomes
  
  #############################################################
  # output map and meta data
  #############################################################
  
  #  RF plot
  png(paste0("data/rqtl/plots/", cross_id, ".png"), width = 800, height = 800)
  plotRF(qtl_reorder_clean, what = "rf")
  dev.off()
  
  pdf(paste0("data/rqtl/plots/", cross_id, ".pdf"), width = 8.5, height = 8.5)
  plotRF(qtl_reorder_clean, what = "rf")
  dev.off()
  
  #mix_qtl <- 
  
  #par(bg = 0, col.lab = "white", col.axis = "white")
  #plotRF(qtl_dat, what = "rf")
  #plotRF(qtl_reorder, what = "rf")
  #plotRF(qtl_reorder_clean, what = "rf")
  
  # write map to find
  cross_id <- geno_file %>% gsub("_geno_dat.csv", "", .)
  dir.create(paste0("data/rqtl/maps/", cross_id))
  write.cross(qtl_reorder_clean, "csv", filestem = paste0("data/rqtl/maps/", cross_id))
  
  # tally crossovers and write to file
  count_cross_overs <- function(chr){
    
    id <- countXO(qtl_reorder_clean, chr = chr) %>% names
    crossovers <-  as.numeric(countXO(qtl_reorder_clean, chr = chr))
    data.frame(id, chr, crossovers)
    
  }
  
  co_df <- suppressWarnings(lapply(names(qtl_reorder_clean$geno), count_cross_overs) %>% bind_rows)
  
  write.table(co_df, paste0("data/rqtl/co_estimates/", cross_id, ".txt"), row.names = FALSE, quote = FALSE)
  
  # extract marker orders
  chr2 <- qtl_reorder_clean$geno$`2`$map %>% data.frame(chr = "2", marker = names(.), map_pos = as.numeric(.), row.names = NULL) %>% select(chr, marker, map_pos)
  chr3 <- qtl_reorder_clean$geno$`3`$map %>% data.frame(chr = "3", marker = names(.), map_pos = as.numeric(.), row.names = NULL) %>% select(chr, marker, map_pos)
  chr4 <- qtl_reorder_clean$geno$`4`$map %>% data.frame(chr = "4", marker = names(.), map_pos = as.numeric(.), row.names = NULL) %>% select(chr, marker, map_pos)
  chrX <- qtl_reorder_clean$geno$`X`$map %>% data.frame(chr = "X", marker = names(.), map_pos = as.numeric(.), row.names = NULL) %>% select(chr, marker, map_pos)

  marker_order <- bind_rows(chr2, chr3, chr4, chrX)
  marker_order <- data.frame(cross_id, marker_order)
  
  write.table(marker_order, paste0("data/rqtl/marker_orders/", cross_id, ".txt"), row.names = FALSE, quote = FALSE)
  
  # output list of individuals used to build the map
  # i.e. which were filtered due to data quality issues
  # similar data for markers is implicit (i.e. if they were filtered they dont have map positions)
  included_inds <- data.frame(countXO(qtl_reorder_clean))
  included_inds$id <- row.names(data.frame(countXO(qtl_reorder_clean))) 
  included_inds <- included_inds %>%
    rename(xo_count_rqtl = countXO.qtl_reorder_clean.) %>%
    select(id, xo_count_rqtl)
  
  write.table(included_inds, paste0("data/rqtl/marker_orders/", cross_id, "_indivs.txt"), row.names = FALSE, quote = FALSE)
  
}

cross_ids <- list.files("data/rqtl/raw", pattern = "geno") %>% gsub("_geno_dat.csv", "", .)

lapply(cross_ids, create_genetic_map, data_dir = "data/rqtl/raw", 
       missing_marker_cutoff = 0.05,
       missing_ind_cutoff = 0.05,
       marker_co_cutoff = 0.95,
       max_cross_overs = 15)

# analyze co data


# plot concordance between maps for chromosome XL

#############################################################
# Supplementary Figures 1-3
#############################################################

mark_df <- lapply(list.files("data/rqtl/marker_orders", full.names = TRUE), read.table, h = T) %>% bind_rows()

# assign marker numbers to the "all" set
marker_ids <- mark_df %>% filter(cross_id == "all") %>% arrange(chr, map_pos) %>% select(chr, marker, map_pos)
marker_ids <- marker_ids %>% group_by(chr) %>% do(marker = .$marker, mark_num = dense_rank(.$map_pos)) %>% unnest

# join in marker numbers to the varied set

mark_df %>%
  left_join(marker_ids) %>%
  filter(chr == 2) %>%
  ggplot(aes(x = cross_id, y = map_pos, group = mark_num))+
  geom_point()+
  geom_path()

library("GGally")

mark_df %>%
  left_join(marker_ids) %>%
  filter(!is.na(mark_num)) %>%
  filter(chr == "X") %>%
  group_by(cross_id) %>%
  mutate(map_pos_rel = map_pos/max(map_pos)) %>%
  ungroup %>%
  select(cross_id, chr, mark_num, map_pos_rel) %>%
  spread(key = cross_id, value = map_pos_rel) %>%
  ggpairs(columns = 3:20)
  ggplot(aes(x = cross_id, y = map_pos_rel, group = mark_num))+
  geom_point()+
  geom_path()
  
  
chrom_order_plot <- mark_df %>%
  left_join(marker_ids) %>%
  filter(!is.na(mark_num)) %>%
  filter(chr == "X") %>%
  filter(cross_id!="all") %>%
  group_by(cross_id) %>%
  mutate(map_pos_rel = map_pos/max(map_pos)) %>%
  mutate(map_pos_rank = dense_rank(map_pos)) %>%
  ungroup %>%
  ggplot(aes(x = cross_id, y = map_pos_rank, group = mark_num, fill = mark_num))+
  geom_path(color = "grey75")+
  geom_point(size = 3, pch = 22)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  ylab("Marker Rank on Chr X")+
  xlab("Cross")+
  scale_fill_viridis_c()

 ggsave("figures/FigureS1_new.pdf", plot = chrom_order_plot, height = 8.5, width = 8.5)
  

dotplot <- mark_df %>%
  left_join(marker_ids) %>%
  filter(!is.na(mark_num)) %>%
  #filter(chr == "X") %>%
  filter(cross_id!="all") %>%
  group_by(cross_id) %>%
  mutate(map_pos_rel = map_pos/max(map_pos)) %>%
  mutate(map_pos_rank = dense_rank(map_pos)) %>%
  ungroup %>%
  mutate(flip = ifelse(chr == 3 & cross_id == "AFC_19", 1, 0)) %>%
  mutate(flip = ifelse(chr == 2 & cross_id %in% c("AFC_14", "AFC_30", "AFC_49", "AFC_56", 
                                                  "MC_13", "MC_14", "MC_15"), 1, flip)) %>%
  group_by(cross_id, chr) %>%
  mutate(map_pos_rank = ifelse(flip == 1, rev(map_pos_rank), map_pos_rank)) %>%
  ggplot(aes(x = mark_num, y = map_pos_rank, group = mark_num))+
  geom_point(size = 1)+
  theme_bw()+
  facet_wrap(chr~cross_id, scales = "free", nrow = 4, ncol = 17)+
  xlab("Marker rank on combined map")+
  ylab("Marker rank on individual map") 
  ggsave("figures/FigureS2_raw.pdf", plot = dotplot, height = 6.5, width = 10.5)
  
  

# concordance of CO estimates
co_df <- lapply(list.files("data/rqtl/co_estimates", full.names = TRUE, pattern = "_.*"), read.table, h = T) %>% bind_rows()
co_df_all <- lapply(list.files("data/rqtl/co_estimates", full.names = TRUE, pattern = "all"), read.table, h = T) %>% bind_rows()

names(co_df_all)[3] <- "crossovers_all"
co_df_all <- left_join(co_df, co_df_all) 
co_df_all <- co_df_all %>% filter(complete.cases(co_df_all))
cor(x = co_df_all$crossovers, y= co_df_all$crossovers_all, method = "pearson")
# [1] 0.9174831

cor.test(x = co_df_all$crossovers, y= co_df_all$crossovers_all, method = "spearman")

