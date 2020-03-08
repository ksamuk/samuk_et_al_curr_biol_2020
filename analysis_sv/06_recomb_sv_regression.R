# summarize sv data for each window
# KMS 2019

library("tidyverse")
library("parallelsugar")
library("wesanderson")
library("scales")
#library("lme4")
#library("visreg")
#library("MASS")
#library("pscl")
library("patchwork")
#options(device = "quartz")

######################################
# read in summarized SV data
######################################

sv_df <- readRDS("analysis_sv/data/sv_window_summary.rds") %>%
  ungroup %>%
  rename(POS_schf1 = pos1, POS_schf2 = pos2, line = Indiv)

recomb_df <- readRDS("data/crossovers_fine_scale.rds") %>% 
  rowwise %>%
  mutate(pos1_schf = min(POS_schf1, POS_schf2),
         pos2_schf = max(POS_schf1, POS_schf2)) %>%
  mutate(POS_schf1 = pos1_schf, POS_schf2 = pos2_schf) %>%
  select(-pos1_schf, -pos2_schf) %>% 
  ungroup

# build a bed file for the recombination regions

# read in map positon data
ord_df <- read.table("data/rqtl/marker_orders/all.txt", h = T) %>%
  mutate(CHROM_schf = gsub("*._chr_", "", marker) %>% 
           gsub("_pos.*", "", .)) %>%
  mutate(POS_schf = as.numeric(gsub(".*_pos_", "", marker))) %>%
  group_by(chr) %>%
  mutate(rqtl_rank = dense_rank(map_pos)) %>%
  ungroup %>%
  rename(rqtl_chr = chr, rqtl_map_pos = map_pos) %>%
  select(CHROM_schf, POS_schf, rqtl_chr, rqtl_map_pos, rqtl_rank) 

# read in remapped position data
remap_df <- readRDS("analysis_cm_mb/out/remapped_positions.rds")

map_df <- left_join(ord_df, remap_df)

bed_out <- recomb_df %>%
  ungroup %>%
  select(CHROM_schf, POS_schf1, POS_schf2) %>%
  group_by(CHROM_schf, POS_schf1, POS_schf2) %>%
  tally %>% 
  arrange(CHROM_schf, POS_schf1) %>%
  filter(n > 100) %>%
  select(-n)

bed_out1 <- bed_out %>%
  ungroup %>%
  select(CHROM_schf, POS_schf1) %>%
  rename(POS_schf = POS_schf1) %>%
  left_join(map_df) %>%
  select(CHROM_schf, POS_schf, CHROM_ref, POS_ref)%>%
  rename(POS_schf1 = POS_schf, POS_ref1 = POS_ref)

bed_out2 <- bed_out %>%
  ungroup %>%
  select(CHROM_schf, POS_schf2) %>%
  rename(POS_schf = POS_schf2) %>%
  left_join(map_df) %>%
  select(POS_schf, POS_ref) %>%
  rename(POS_schf2 = POS_schf, POS_ref2 = POS_ref)

cbind(bed_out1, bed_out2) %>%
  rowwise %>%
  mutate(pos1 = min(POS_ref1, POS_ref2),
         pos2 = max(POS_ref1, POS_ref2)) %>%
  select(CHROM_ref , pos1, pos2) %>%
  mutate(pos2 = pos2 - 1) %>%
  filter((pos2 - pos1) > 10000) %>% 
  filter((pos2 - pos1) < 1500000) %>% 
  rename(CHROM = CHROM_ref, POS = pos1, POS_TO = pos2) %>%
  arrange(CHROM, POS) %>%
  distinct %>%
  write.table("analysis_sv/meta/recom_intervals.txt", sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)

cbind(bed_out1, bed_out2) %>%
  rowwise %>%
  mutate(pos1 = min(POS_ref1, POS_ref2),
         pos2 = max(POS_ref1, POS_ref2)) %>%
  mutate(pos1_schf = min(POS_schf1, POS_schf2),
         pos2_schf = max(POS_schf1, POS_schf2)) %>%
  mutate(pos2 = pos2 - 1) %>%
  filter((pos2 - pos1) > 10000) %>% 
  filter((pos2 - pos1) < 1500000) %>% 
  distinct %>%
  select(CHROM_ref, POS_ref1, POS_ref2, CHROM_schf, POS_schf1, POS_schf2, pos1, pos2, pos1_schf, pos2_schf) %>% 
  write.table("analysis_sv/meta/recom_intervals_ref_schf.txt", sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = TRUE)

# summarize recombination data for each line and window w/i in
recomb_sum_df <- recomb_df %>% 
  mutate(line = gsub("_[0-9]{1,2}_.{3}$", "",  Indiv) %>% 
           gsub("_", "", .) %>%
           gsub("AFC", "A", .) %>%
           gsub("MC", "M", .)) %>%
  mutate(POS_schf2 = POS_schf2 - 1) %>%
  group_by(line, CHROM_schf, POS_schf1, POS_schf2) %>%
  summarise(crossover_count = sum(crossover_int, na.rm = TRUE)) %>%
  arrange(line, CHROM_schf, POS_schf1)

# remove invalid windows (from bad interval conversion)
recomb_sum_df <- recomb_sum_df %>%
  filter(abs(POS_schf1 - POS_schf2) < 1500000)

saveRDS(recomb_sum_df, "analysis_sv/recomb_counts_fine_scale.rds")

# read in the fine scale counts and flip windows to match sv_df (listing smallest position as POS1)
recomb_sum_df <- readRDS("analysis_sv/recomb_counts_fine_scale.rds") %>%
  rowwise %>%
  mutate(pos1_schf = min(POS_schf1, POS_schf2),
         pos2_schf = max(POS_schf1, POS_schf2)) %>%
  mutate(POS_schf1 = pos1_schf, POS_schf2 = pos2_schf) %>%
  select(-pos1_schf, -pos2_schf) %>% 
  ungroup

# # # check for harmonization  of positions
# tmp1 <- filter(recomb_df, CHROM_schf == "XL") %>%
#   ungroup %>%
#   select(POS_schf1, POS_schf2) %>%
#   distinct %>%
#   arrange(POS_schf1)%>%
#   mutate(dataset = "recomb_sum")
# 
# tmp2 <- filter(sv_df, CHROM_schf == "XL") %>%
#   select(POS_schf1, POS_schf2) %>%
#   distinct %>%
#   arrange(POS_schf1) %>%
#   mutate(dataset = "sv_df")
# 
# full_join(tmp2, tmp1, by = "POS_schf1")

# filter out any windows that are suspiciously small or larger (all should be ~200kb)
# also repair off by one error for POS_schf1

recomb_sum_sv <- sv_df %>%
  mutate(CHROM_schf = as.character(CHROM_schf)) %>%
  mutate(POS_schf1 = POS_schf1 + 1) %>%
  mutate(POS_schf2 = POS_schf2 - 1) %>%
  rowwise %>%
  mutate(pos1 = min(POS_schf1, POS_schf2),
         pos2 = max(POS_schf1, POS_schf2)) %>%
  mutate(POS_schf1 = pos1, POS_schf2 = pos2) %>%
  select(-pos1, -pos2) %>% 
  ungroup %>%
  left_join(recomb_sum_df, by = c("line", "CHROM_schf", "POS_schf1", "POS_schf2")) %>%
  #filter(abs(POS_schf1 - POS_schf2) < 1500000)%>%
  #filter(abs(POS_schf1 - POS_schf2) > 50000) %>%
  filter(!is.na(method))


###################################################################
# DEPTH NORMALIZATION
###################################################################

recomb_intervals <- read.table("analysis_sv/meta/recom_intervals_ref_schf.txt", h = T) %>%
  rowwise %>%
  mutate(pos1_ref = min(POS_ref1, POS_ref2),
         pos2_ref = max(POS_ref1, POS_ref2)) %>%
  mutate(POS_schf1 = pos1_schf, POS_schf2 = pos2_schf) %>%
  select(CHROM_ref, pos1_ref, pos2_ref, POS_schf1, POS_schf2, CHROM_schf) %>% 
  ungroup 

depth_df <- readRDS("analysis_sv/data/depth.rds") %>% 
  select(-pos1, -pos2) %>% 
  ungroup %>%
  select(-pos1_ref, -pos2_ref) %>%
  left_join(recomb_intervals) %>%
  mutate(CHROM_schf = as.character(CHROM_schf)) %>%
  rowwise %>%
  mutate(pos1 = min(POS_schf1, POS_schf2),
         pos2 = max(POS_schf1, POS_schf2)) %>%
  mutate(POS_schf1 = pos1, POS_schf2 = pos2) %>%
  select(-pos1, -pos2) %>%
  ungroup

depth_df_pb <- depth_df %>%
  filter(method != "illumina") %>%
  mutate(method = "pbsv")

depth_df_gatk <- depth_df %>%
  filter(method != "pacbio") %>%
  mutate(method = "gatk")

depth_df_smoove <- depth_df %>%
  filter(method != "pacbio") %>%
  mutate(method = "smoove")

depth_df_popte <- depth_df %>%
  filter(method != "pacbio") %>%
  mutate(method = "popte2")

depth_df <- bind_rows(depth_df_pb, depth_df_gatk, depth_df_smoove, depth_df_popte) %>%
  rename(line = Indiv) %>%
  mutate(line = gsub("FC|C", "", line)) %>%
  group_by(line, method, CHROM_schf, POS_schf1) %>%
  summarise(depth = sum(sum_depth, na.rm = TRUE) / sum(n_sites, na.rm = TRUE), 
            n_sites = sum(n_sites, na.rm = TRUE), coverage_1k = sum(sum_depth > 0, na.rm = TRUE), 
            coverage_prop = (coverage_1k/n_sites)*1000) %>%
  ungroup

depth_df %>% filter(CHROM_schf == "XR") %>% pull(POS_schf1) %>% unique %>% sort
recomb_sum_sv %>% filter(CHROM_schf == "XR") %>% pull(POS_schf1) %>% unique %>% sort

rm(list = c("depth_df_pb", "depth_df_gatk", "depth_df_smoove", "depth_df_popte"))

# allow for a fuzzy left join
recomb_sum_sv <- recomb_sum_sv %>%
  mutate(pos_join = round(POS_schf1/1000, 0) * 1000)

recomb_sum_sv <- depth_df %>%
  mutate(pos_join = round(POS_schf1/1000, 0) * 1000) %>%
  select(-POS_schf1)  %>%
  left_join(recomb_sum_sv,., by = c("line", "method", "CHROM_schf", "pos_join")) %>%
  filter(!is.na(depth)) %>%
  select(-pos_join)

###################################################################
# FURTHER QC FOR SVS
###################################################################

# there is a giant deletion on XL, detected by all methods
# given the size of the deletion, i suspect this might be an assembly artifact
# excelsior!

recomb_sum_sv %>%
  filter(CHROM_schf == "XL") %>%
  ggplot(aes(x = POS_schf1, y = diff_count))+
  geom_point()+
  facet_grid(line~.)


###################################################################
# FIGURES
###################################################################

plot_pal <- c(wes_palette("Darjeeling1", 5) %>% as.character, (wes_palette("Darjeeling2", 2) %>% as.character)[2])

sv_sub <- readRDS("analysis_sv/data/sv_joined.rds") %>%
  select(method, pos1, pos2, sv_type, sv_length) %>%
  filter(!(sv_type == "REF"), !(sv_type == "BND"), !is.na(sv_type)) %>%
  filter(!(is.na(method) & is.na(sv_type))) %>%
  distinct %>%
  mutate(sv_type = gsub("Comp.*", "COMP", sv_type)) %>%
  mutate(sv_length = ifelse(method =="popte2", abs(pos2-pos1), sv_length)) %>%
  mutate(sv_type = factor(sv_type, levels = c("INS", "DEL", "COMP", "DUP", "INV", "TE")))

################
# FIGURE 5
################

figure_5a <- recomb_sum_sv %>%
  filter(!(sv_type == "REF")) %>%
  filter(!(is.na(method) & is.na(sv_type))) %>%
  mutate(sv_type = gsub("Comp.*", "COMP", sv_type)) %>%
  mutate(sv_type = factor(sv_type, levels = c("INS", "DEL", "COMP", "DUP", "INV", "TE"))) %>%
  #mutate(size_diff = abs(line_length - mv_length)) %>%
  #mutate(size_kb = size_diff / 1000) %>%
  mutate(read_length = ifelse(method == "pbsv", "long", "short")) %>%
  group_by(read_length, sv_type) %>%
  summarise(n = sum(n_genotyped, na.rm = TRUE)) %>%
  ggplot(aes(y = n, x = sv_type, fill = sv_type))+
  geom_bar(stat = "identity")+
  ylab("Count")+
  xlab("Variant Class")+
  facet_grid(~read_length , scales = "free_x", space = "free")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = plot_pal)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

figure_5b <- sv_sub %>%
  filter(method != "popte2") %>%
  filter(sv_length != 0) %>% 
  #sample_frac(0.01) %>%
  mutate(read_length = ifelse(method == "pbsv", "long", "short")) %>%
  distinct %>%
  ggplot(aes(y = ceiling(abs(sv_length)), x = sv_type, fill = sv_type, outlier.fill = sv_type))+
  geom_boxplot(outlier.size = 0.75, outlier.alpha = 0.1, outlier.shape = 16)+
  #geom_jitter(alpha = 0.5, shape = 16)+
  ylab("Size (bp)")+
  xlab("Variant Class")+
  facet_grid(~read_length, scales = "free_x", space = "free")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = plot_pal)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = plot_pal)

figure5 <- figure_5a/figure_5b

ggsave("figures/Figure5_raw.pdf", figure5, width = 6, height = 6, useDingbats = FALSE)

################
# FIGURE S3?
################

recomb_sum_sv %>%
  filter(method == "pbsv") %>%
  ggplot(aes(x = depth))+
  geom_histogram()+
  facet_wrap(~line)

recomb_sum_sv %>%
  #filter(!line %in% c("A60", "A48", "A49")) %>%
  filter(!is.na(sv_type)) %>%
  filter(!is.na(method)) %>%
  group_by(method, line, CHROM_schf, POS_schf1) %>%
  filter(depth > 5) %>%
  filter(!(method == "pbsv" & line %in% c("A48", "A49", "M14", "M17"))) %>% # remove lines with v. low pacbio depth
  mutate(diff_count_norm = ifelse(method == "pbsv", (diff_count/coverage_1k)/(depth^2), (diff_count/coverage_1k)/(depth))) %>%
  summarise(sum_length_diff = sum(len_diff, na.rm = TRUE) , 
            sum_coverage = sum(coverage_1k*1000, na.rm = TRUE),
            mean_count_diff = mean(diff_count_norm, na.rm = TRUE),
            mean_depth = mean(depth, na.rm = TRUE),
            mean_coverage = mean(coverage_prop, na.rm = TRUE)) %>%
  ggplot(aes(y = mean_count_diff, x = line , color = mean_depth))+
  #geom_boxplot()+
  #geom_point()+
  geom_jitter(alpha = 0.5, size = 2, shape = 16)+
  #facet_grid(CHROM_schf~method, scales = "free", space = "free")+
  facet_wrap(CHROM_schf~method, scales = "free", ncol = 4)+
  theme_bw()+
  theme(strip.background = element_blank(),
        text = element_text(size=12))+
  scale_color_viridis_c(option = "D")

################
# FIGURE 6
################

plot_pal <- c(wes_palette("Darjeeling1", 3) %>% as.character)

sv_zscore <- recomb_sum_sv %>%
  filter(!is.na(method), method != "popte2") %>%
  filter(!is.na(sv_type)) %>%
  filter(!is.na(method)) %>%
  group_by(method, line, CHROM_schf, POS_schf1, crossover_count) %>%
  filter(depth > 5) %>%
  filter(!(method == "pbsv" & line %in% c("A48", "A49", "A60", "M14", "M17"))) %>% # remove lines with v. low pacbio depth
  mutate(diff_count_norm = ifelse(method == "pbsv", (diff_count/coverage_1k)/(depth^2), (diff_count/coverage_1k)/(depth))) %>%
  mutate(len_diff_norm = ifelse(method == "pbsv", (len_diff/coverage_1k)/(depth^2), (len_diff/coverage_1k)/(depth))) %>%
  summarise(sum_length_diff = sum(len_diff, na.rm = TRUE) , 
            sum_length_diff_norm = sum(abs(len_diff_norm) , na.rm = TRUE) , 
            sum_coverage = sum(coverage_1k*1000, na.rm = TRUE),
            mean_count_diff = mean(diff_count_norm, na.rm = TRUE),
            mean_depth = mean(depth, na.rm = TRUE),
            mean_coverage = mean(coverage_prop, na.rm = TRUE),
            crossover_mean = mean(crossover_count, na.rm = TRUE)) %>%
  ungroup %>%
  group_by(method, CHROM_schf, POS_schf1) %>%
  mutate(line = line, crossover_z = scale(crossover_count), 
         length_diff_z = scale(sum_length_diff, center = FALSE), 
         count_diff_z = scale(mean_count_diff))



figure_6a <- sv_zscore %>%
  ggplot(aes(x = count_diff_z, y = crossover_z, color = method, fill = method))+
  geom_point(size = 0.5, alpha = 0.4, shape = 21)+
  #geom_smooth(size = 1, alpha = 0.3, level = 0.95, span = 0.8)+
  stat_smooth(size = 1, alpha = 0.3, level = 0.95, span = 1.0, n = 80, method = "loess")+
  facet_wrap(~method, scale = "free_x") +
  xlab("Difference in Window SV Count (Z-score)")+
  ylab("Recombination Rate (Z-score)")+
  theme_bw()+
  theme(legend.title = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    text = element_text(size=12))+
  scale_color_manual(values = plot_pal)+
  scale_fill_manual(values = plot_pal)

figure_6b <- sv_zscore  %>%
  ggplot(aes(x = abs(sum_length_diff)+1, y = crossover_z, color = method, fill = method))+
  geom_point(size = 0.5, alpha = 0.4, shape = 21)+
  stat_smooth(size = 1, alpha = 0.3, level = 0.95, span = 1.0, n = 80, method = "loess")+
  #geom_smooth(size = 1, alpha = 0.3, alpha = 0.3, level = 0.95)+
  #geom_smooth(method = "lm", size = 1, alpha = 0.3, alpha = 0.3, level = 0.95, data = filter(sv_zscore, length_diff_z > 0))+
  #geom_smooth(method = "lm", size = 1, alpha = 0.3, alpha = 0.3, level = 0.95, data = filter(sv_zscore, length_diff_z < 0))+
  #geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~method, scale = "free_x") +
  xlab("Difference in Sequence Length (bp)")+
  ylab("Recombination Rate (Z-score)")+
  theme_bw()+
  theme(legend.title = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    text = element_text(size=12))+
  scale_color_manual(values = plot_pal)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = plot_pal)
 
figure6 <- figure_6a / figure_6b

ggsave("figures/Figure5_raw.pdf", figure6, width = 8, height = 6.5, useDingbats = FALSE)

###################################################################
# LINEAR MODELS 
###################################################################

sv_zscore <- sv_zscore %>% 
  ungroup %>%
  mutate(window_unique = paste0(CHROM_schf, "_", POS_schf1))
library("lme4")

# TEST OF LENGTH EFFECT

sv_mod1 <- glmer(crossover_count ~ method + (1|window_unique)+ (1|line), family = "poisson", data = sv_zscore)
sv_mod2 <- glmer(crossover_count ~ method * abs(sum_length_diff_norm) + (1|window_unique)+ (1|line), family = "poisson", data = sv_zscore)

#visreg::visreg(sv_mod2)
anova(sv_mod1, sv_mod2)

# Data: sv_zscore
# Models:
# sv_mod1: crossover_count ~ method + (1 | window_unique) + (1 | line)
# sv_mod2: crossover_count ~ method * abs(sum_length_diff_norm) + (1 | window_unique) + 
# sv_mod2:     (1 | line)
# Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# sv_mod1  5 24267 24300 -12128    24257                         
# sv_mod2  8 24270 24323 -12127    24254 2.9533      3     0.3989

# TEST OF COUNT EFFECT

sv_mod1 <- glmer(crossover_count ~  method +  (1|window_unique)+ (1|line), family = "poisson", data = sv_zscore)
sv_mod2 <- glmer(crossover_count ~  method + mean_count_diff + (1|window_unique)+ (1|line), family = "poisson", data = sv_zscore)
sv_mod3 <- glmer(crossover_count ~  method * mean_count_diff+ (1|window_unique)+ (1|line), family = "poisson", data = sv_zscore)

#visreg::visreg(sv_mod2)
anova(sv_mod1, sv_mod2, sv_mod3)

# Data: sv_zscore
# Models:
# sv_mod1: crossover_count ~ method + (1 | window_unique) + (1 | line)
# sv_mod2: crossover_count ~ method * mean_count_diff + (1 | window_unique) + 
# sv_mod2:     (1 | line)
# Df   AIC   BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# sv_mod1  5 24267 24300 -12128    24257                         
# sv_mod2  8 24272 24325 -12128    24256 1.1637      3     0.7617





