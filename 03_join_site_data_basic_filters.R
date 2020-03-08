# format and filter raw genotypes from GT-seq
# KMS 2018-2019

library("tidyverse")
library("viridis")
library("patchwork")

#############################################################
# read in raw genotypes and check which markers are informative
#############################################################

# read in the gt seq vcf
gt_df <- readRDS(file = "data/gt_df/gt_seq_all_gt_df_gatk.rds")

# get the identities of all known informative snps

snp_ident <- readRDS("meta/sampled_gt_seq_snps_all.rds") %>%
  #select(CHROM, POS, unique_in, pair_id) %>%
  filter(Indiv %in% c("Flag14", "MV2-25")) %>%
  select(CHROM, POS, pair_id, Indiv, unique_in, gt_GT_alleles) %>%
  spread(key = Indiv, value = gt_GT_alleles) %>%
  select(-pair_id)

names(snp_ident)[4:5] <- c("flag_alelle", "mv_allele") 

# join in unique_in (informativeness) info for each snp

gt_df <- gt_df %>%
  left_join(snp_ident)

# check number of informative sites
# should be ~700 (350 * 2)
gt_df %>%
  filter(!is.na(unique_in)) %>%
  select(CHROM, POS) %>%
  distinct %>%
  nrow

# filter out non-informative sites (!)
# also low depth samples and any non-biallelic sites
gt_df <- gt_df %>%
  distinct %>%
  arrange(CHROM, POS, Indiv) %>%
  filter(!is.na(unique_in)) %>%
  filter(gt_GT %in% c("0/0", "0/1", "1/1")) %>%
  filter(gt_DP > 7)

gt_df %>%
  filter(!is.na(unique_in)) %>%
  select(CHROM, POS) %>%
  distinct %>%
  nrow

# 687 sites (~1/2 that are mapping informative)
# so like 340 that are mapping informative

# get remapped positions
# see "analysis_cm_mb" folder
remapped_df <- readRDS("analysis_cm_mb/out/remapped_positions.rds") %>%
  rename(CHROM = CHROM_ref, POS = POS_ref)

gt_df <- gt_df %>%
  left_join(remapped_df) 

#############################################################
# genotyping basic filters
#############################################################

# and/or have depth 
#gt_df_sub <- gt_df %>%
#  filter(gt_DP > 20)

# mean representation of individuals
depth_df <- gt_df %>%
  group_by(plate, Indiv) %>%
  summarise(mean_dp = mean(gt_DP, na.rm = TRUE))

# join into gt_df for filtering
gt_df <- left_join(gt_df, depth_df)

#############################################################
# exploratory plots
#############################################################

# genotypes per individual
# expecting ~350*2 max
gt_df %>%
  group_by(Indiv) %>%
  tally

# mean depth per amplicon
# expecting ???
gt_df %>%
  mutate(pair_id = paste0(CHROM,"_",POS)) %>%
  group_by(pair_id) %>%
  summarise(mean_dp = mean(gt_DP)) %>%
  ggplot(aes(x = pair_id, y = mean_dp))+
  geom_point()

# mean representation of barcodes per well
gt_df %>%
  mutate(plate_row = gsub("[0-9]*", "", well) %>% factor(., levels = rev(sort(unique(.))))) %>%
  mutate(plate_col = gsub("[^0-9]*", "", well)) %>%
  group_by(plate_row, plate_col) %>%
  summarise(mean_depth = mean(gt_DP, na.rm = TRUE)) %>%
  ggplot(aes(x = plate_col, y = plate_row, fill = mean_depth))+
  geom_raster()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_viridis()+
  ylab("Row") + xlab("Column")

### how many individuals are acutally usable?
gt_df %>% 
  select(plate, Indiv, mean_dp) %>%
  distinct %>%
  filter(mean_dp >= 20) %>%
  group_by(plate) %>%
  tally %>%
  pull(n) 

# tonig info file
# built from the contig metadata in the VCF header
contig_df <- read.csv("meta/dpse_ref_contig_info.csv", h= T) %>%
  mutate(contig_mid = (contig_start + contig_end)/2)

names(contig_df)[1] <- "CHROM"

Figure_S1A <- gt_df %>%
  filter(!(unique_in == "flag")) %>%
  mutate(CHROM_group = gsub("_.*", "", CHROM)) %>%
  select(CHROM, CHROM_group, POS, unique_in) %>%
  distinct %>%
  left_join(contig_df) %>%
  mutate(POS_adj = POS + contig_start) %>%
  ggplot(aes(x = POS_adj/1000000, y = 0, xend = POS_adj/1000000, yend = 0.5))+ 
  geom_rect(data = contig_df, aes(x = 0, y = 0, xmin = contig_start/1000000, xend = 0, yend = 0.5,
                                  ymin= 0, xmax = contig_end/1000000, ymax = 0.50, fill = CHROM))+
  geom_rect(data = contig_df, aes(x = 0, y = 0, xmin = 0, xend = 0, yend = 0.5,
                                  ymin= 0, xmax = CHROM_group_length/1000000, ymax = 0.50), 
            color = "black", size = 0.5, alpha = 0.5, fill = NA)+
  #ggrepel::geom_text_repel(data = contig_df, aes(x = contig_mid/1000000,y = 0.5, label = CHROM_sub_id, xend = 0), 
  #          color = "black", size = 4, direction = "both", min.segment.length = 0, nudge_y = 1.0,
  #          segment.size = 0.25, ylim = 0.75)+
  geom_segment(size = 0.1)+
  facet_grid(CHROM_group~.)+
  theme(text = element_text(family = "Lato", size = 14),
        strip.background = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        panel.spacing.y = unit(2, "mm"),
        axis.ticks.y = element_blank(),
       # axis.line.x = element_line(),
        legend.position = "none")+
  scale_x_continuous(breaks = seq(0, 30, 5))+
  scale_y_discrete(expand = c(0.05,0))+
  scale_fill_manual(values = c(rep(c("white", "grey"), 7), "white"))+
  #scale_fill_manual(values = wesanderson::wes_palette("Zissou1", 15, type = c("continuous")))+
  xlab("Position (Mb)")+
  ylab("")+
  labs(subtitle = "A")

FigureS1B <- gt_df %>%
  filter(!(unique_in == "flag")) %>%
  mutate(CHROM_group = gsub("_.*", "", CHROM)) %>%
  select(CHROM, CHROM_group, POS, unique_in) %>%
  distinct %>%
  group_by(CHROM) %>%
  mutate(pos_diff = POS - lag(POS)) %>%
  ungroup %>% 
  filter(pos_diff > 1000) %>%
  ggplot(aes(x = pos_diff/1000))+
  geom_histogram(bins = 40)+
  geom_vline(xintercept = 267481/1000, size = 0.5, linetype = 3)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  coord_cartesian(expand = FALSE)+
  theme_bw() +
  theme(text = element_text(family = "Lato", size = 14),
        strip.background = element_blank())+
  xlab("Inter-amplicon Distance (kb)")+
  ylab("Frequency")+
  labs(subtitle = "B")

# depth per plate for each run
FigureS1C <- gt_df %>%
  group_by(plate, Indiv) %>%
  summarise(mean_dp = mean(gt_DP, na.rm = TRUE)) %>%
  ungroup %>%
  mutate(run = ifelse(grepl("_15$|_17$|_19$", plate), "1st Pool (NextSeq 500 High Output)", "2nd Pool (NextSeq 500 Mid Output)")) %>%
  ggplot(aes(x = plate, y = mean_dp)) +
  geom_boxplot(size = 0.5, width = 0.75, outlier.size = 0.05)+
  facet_wrap(~run, scales = "free_x")+
  theme_bw() +
  theme(text = element_text(family = "Lato", size = 14),
        strip.background = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        axis.text.x = element_text(angle = 90, size = 6,hjust=1))+
  expand_limits(x = 0)+
  ylab("Mean Depth Per Amplicon")+
  xlab("Plate ID")+
  labs(subtitle = "C")

Figure_S1 <- (Figure_S1A | FigureS1B)/FigureS1C

ggsave(Figure_S1, filename = "figures/raw/FigureS1_raw.pdf", width = 8.5, height = 6, device = cairo_pdf)

# evenness per amplicon
ampli_stats <- gt_df %>%
  group_by(plate, Indiv) %>%
  summarise(mean_dp = mean(gt_DP, na.rm = TRUE), 
            sd_dp = sd(gt_DP, na.rm = TRUE),
            coeff_var = sd_dp/mean_dp,
            n_amplicons = length(gt_DP)/2) %>%
  #filter(n_amplicons > 200) %>% # we filter these individuals out later
  ungroup 

ampli_stats %>%
  gather(key = stat, value = value, -plate, -Indiv) %>%
  ggplot(aes(x = value))+
  geom_histogram()+
  facet_wrap(~stat, scales = "free")+
  theme_bw()


#############################################################
# additional filtering
#############################################################

# representation filter
# i.e. remove individuals with a lot of missing data
# rationale is that if you are missing a lot of data
# your individual genotype calls might be unreliable as well
# (now do this as part of rqtl)

gt_df  %>% 
  filter(unique_in != "flag" & !is.na(unique_in)) %>%
  #left_join(prop_het) %>%
  #filter(gt_GT != "0/0") %>%
  #filter(grepl("X", CHROM)) %>%
  filter(grepl("2$", CHROM_schf)) %>%
  filter(grepl("AFC_14", plate)) %>%
  filter(grepl("H.*|A.*", well)) %>%
  ggplot(aes(x = as.factor(POS_schf), y = Indiv, 
             fill = as.factor(gt_GT)))+
  #ggplot(aes(x = as.factor(), y = Indiv, fill = as.factor(gt_GT)))+
  geom_tile()+
  theme(strip.text.y = element_text(angle = 180),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_discrete(position = "right")+
  xlab("")+
  facet_wrap(~CHROM_schf ,scales = "free_x")

# polarize ancestry
gt_pol <- gt_df %>%
  filter(gt_GT %in% c("0/0", "0/1", "1/1")) %>%
  mutate(gt_adj = case_when(gt_GT_alleles == mv_allele ~ "0/0",
                            unique_in == "mv" & gt_GT == "1/1" & gt_GT_alleles != mv_allele ~ "1/1",
                            unique_in == "mv" & gt_GT == "0/1" ~ "0/1",
                            TRUE ~ gt_GT))

# grab and join in sex data
sex_df <- read.csv("data/sex_df.csv") %>%
  select(Indiv, sex)

gt_pol <- gt_pol %>% 
  left_join(sex_df)

# pre filtering numbers of markers
gt_pol %>%
  filter(unique_in == "mv") %>%
  select(CHROM_schf, POS) %>%
  distinct %>%
  group_by(CHROM_schf) %>%
  tally

# CHECK: allele ratios make sense
# looking for a ~1:1 depth ratio for hets
# all good
gt_pol %>%
  filter(unique_in == "mv") %>%
  filter(POS == 488847, CHROM == "2") %>%
  filter(gt_GT %in% c("0/0", "0/1", "1/1")) %>% 
  separate(gt_AD, into = c("ref_depth", "alt_depth"), convert = TRUE) %>% 
  mutate(allele_ratio = alt_depth/(ref_depth + alt_depth)) %>% 
  ggplot(aes(x = ref_depth, y = alt_depth, color = gt_GT)) +
  geom_point()
  
# CHECK: genotype ratios make sense (per chromosome)  
# some lines have excess of 0/0 on chr4 -- tracked this down to
# mis-crossed individuals (i.e. inds with wrong cross direction)
gt_pol  %>%
  filter(unique_in == "mv") %>%
  filter(grepl("^4", CHROM_schf)) %>%
  mutate(line = gsub("_[0-9]*$", "", plate)) %>%
  group_by(line, CHROM_schf, POS, gt_GT) %>%
  tally %>%
  mutate(n = n/sum(n)) %>%
  ggplot(aes(x = n, fill = gt_GT))+
  geom_histogram() +
  facet_grid(line~gt_GT, scales = "free_y")

# throw out inds with evidence of cross failure (excess of 0/0 on autosomes)
# these are likely the result of mis-sexing individuals for crosses (e.g. accidently a MV male)
# allowed for 0/0 of up to 0.2 to account for differences in marker informativeness between lines

good_inds <- gt_pol  %>%
  filter(unique_in == "mv") %>%
  filter(grepl("^2$|^3$|^4$", CHROM_schf)) %>%
  group_by(Indiv, CHROM_schf) %>%
  summarise(mean_00 = mean(gt_GT == "0/0")) %>%
  summarise(high_chrom_excess = any(mean_00 >= 0.20)) %>%
  filter(!(high_chrom_excess)) %>%
  pull(Indiv)

length(good_inds)

gt_pol <- gt_pol %>%
  filter(Indiv %in% good_inds)

# CHECK: cross failure filter was successful
# expecting to see no individuals with long 0/0 segments
gt_pol  %>% 
  filter(unique_in != "flag" & !is.na(unique_in)) %>%
  filter(grepl("2$", CHROM_schf)) %>%
  filter(grepl("AFC_14", plate)) %>%
  filter(grepl("H.*|A.*", well)) %>%
  ggplot(aes(x = as.factor(POS_schf), y = Indiv, fill = as.factor(gt_adj)))+
    geom_tile()+
    theme(strip.text.y = element_text(angle = 180),
          axis.text = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    scale_y_discrete(position = "right")+
    xlab("")+
  facet_wrap(~CHROM_schf ,scales = "free_x")

# after chasing down some things downstream of this analysis
# it turns out that some anc alleles on chr 2 and 3 were misspecified
# this creates a number of issues in later analysis
# so, while hacky, i just flipped the identities here for specific markers (homs only)
# this probably results in a slight underestimate of recombination in these regions
gt_pol <- gt_pol  %>% 
  mutate(gt_adj = ifelse((POS > 10006213 & POS < 11007213 & grepl("2$", CHROM_schf) & gt_adj == "0/0"), "1/1", gt_adj)) %>%
  mutate(gt_adj = ifelse((POS > 16546545 & grepl("3$", CHROM_schf) & gt_adj == "0/0"), "1/1", gt_adj))

# find markers displaying segregation distortion
# well-behaved markers should have 50% hets, 50% homs (of exactly one type)
# two major classes of distorted sites: 
# 1. incorrectly polarized alleles (50:50 ratio, but hom is wrong genotype)
# 2. site where the marker failed to be informative (e.g. two types of homozygotes)
# the rationale here is to err on the side of knocking out uninformative alleles
# and imputing in corrected ancestry later (using neighbouring sites)

dist_markers_autosomes <- gt_pol  %>%
  filter(sex != "M") %>%
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  filter(unique_in == "mv") %>%
  group_by(line, CHROM_schf, POS_schf, gt_adj) %>%
  tally %>%
  mutate(n = n/sum(n)) %>%
  spread(gt_adj, n) %>%
  mutate(outlier_type = case_when(`0/1` > 0.7 | `0/1` < 0.3 ~ "out_het",
                                  `0/0` > 0.1 ~ "out_hom00",
                                  `1/1` < 0.3 | `1/1` > 0.7 ~ "out_hom11",
                                  TRUE ~ "normal")) %>%
  mutate(maj_hom_allele = case_when(`0/0` > 0.2 &  `1/1` < 0.2 ~ "0/0",
                                    `1/1` > 0.2 &  `0/0` < 0.2 ~ "1/1",
                                  TRUE ~ "1/1"))

                                      
# plot distorted markers in context
gt_pol  %>% 
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  filter(unique_in != "flag") %>%
  left_join(dist_markers_autosomes) %>%
  filter(outlier_type != "out_het" ) %>%
  filter(grepl("AFC_14", plate)) %>%
  filter(grepl("H.*", well)) %>%
  ggplot(aes(x = as.factor(POS_schf), y = Indiv, 
             fill = as.factor(gt_adj), color = outlier_type))+
  geom_tile()+
  theme(strip.text.y = element_text(angle = 180),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_discrete(position = "right")+
  xlab("")+
  facet_grid(sex~CHROM_schf ,scales = "free")+
  scale_color_manual(values = c("red", "yellow", NA))

# after looking at these plots, its become clear that the issue is 
# a result of markers with excess 0/0 genotypes
# probably a result of non-differentiation between lines etc. 
# simpliest and least biased way i can think of dealing with this is to 
# simply toss all the 0/0 genotypes (except for XL and XR for the males)
# and then fill the missing genos back in during the imputation step

gt_pol <- gt_pol %>%
  mutate(gt_adj = ifelse(gt_adj == "0/0" & !grepl("X", CHROM), NA, gt_adj)) %>%
  mutate(gt_adj = ifelse(gt_adj == "0/0" & sex == "F" & grepl("X", CHROM_schf), NA, gt_adj)) %>%
  mutate(gt_adj = ifelse(gt_adj == "0/1" & sex == "M" & grepl("X", CHROM_schf), NA, gt_adj)) %>%
  filter(!is.na(gt_adj)) 

dist_markers_autosomes <- gt_pol  %>%
  filter(sex != "M") %>%
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  filter(unique_in == "mv") %>%
  group_by(line, CHROM_schf, POS_schf, gt_adj) %>%
  tally %>%
  mutate(n = n/sum(n)) %>%
  spread(gt_adj, n) %>%
  mutate(outlier_type = case_when(`0/1` > 0.7 | `0/1` < 0.3 ~ "out_het",
                                  `1/1` < 0.3 | `1/1` > 0.7 ~ "out_hom11",
                                  TRUE ~ "normal"))

gt_pol  %>% 
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  filter(unique_in != "flag") %>%
  left_join(dist_markers_autosomes) %>%
  filter(outlier_type != "out_het" ) %>%
  filter(grepl("AFC_14", plate)) %>%
  filter(grepl("H.*", well)) %>%
  ggplot(aes(x = as.factor(POS_schf), y = Indiv, 
             fill = as.factor(gt_adj), color = outlier_type))+
  geom_tile()+
  theme(strip.text.y = element_text(angle = 180),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_discrete(position = "right")+
  xlab("")+
  facet_grid(sex~CHROM_schf ,scales = "free")+
  scale_color_manual(values = c("red", "yellow", NA))

saveRDS(gt_pol, "data/gt_df/gt_seq_filtered.rds")
