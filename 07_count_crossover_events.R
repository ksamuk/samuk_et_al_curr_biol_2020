library("tidyverse")
library("parallel")
library("parallelsugar")
library("aod")
library("rvg")
library("officer")


# options(device = "quartz")

# the ancestry segments
seg_df <- readRDS("data/crossover_segments_min_len_2_rqtl.rds")

# representative individual genotype
seg_df %>% 
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  filter(grepl("AFC_60", line)) %>%
  #filter(sex == "F") %>%
  #filter(unique_in != "flag") %>%
  filter(CHROM_schf == "XL") %>%
  #filter(Indiv == "AFC_19_19_G05") %>%
  ggplot(aes(x = as.factor(rqtl_rank), y = Indiv, 
             fill = as.factor(gt_adj)))+
  geom_tile()+
  theme(strip.text.y = element_text(angle = 180),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_y_discrete(position = "right")+
  xlab("")+
  facet_grid(~CHROM_schf ,scales = "free")

###########################################
# tally crossovers
###########################################

# CHROMOSOME-WIDE

# count crossovers:
# 1. filter out short segments (< 3 contiguous markers in length)
# 2. count a crossover as a switch in segment genotype


# determine which segments are the result of cross overs
co_df <- seg_df %>%
  filter(len_seg >= 2) %>%
  mutate(plate = gsub("_[A-H][0-9]{2}$", "", Indiv))%>%
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  mutate(pop = gsub("_.*$", "", Indiv)) %>%
  group_by(Indiv, CHROM_schf) %>%
  mutate(crossover = ifelse(gt_adj != lag(gt_adj), TRUE, FALSE)) %>%
  mutate(crossover = ifelse(is.na(crossover), FALSE, crossover)) %>%
  ungroup

# CHROMOSOME-LEVEL
co_counts <- co_df %>%
  filter(!(plate %in% c("AFC_19_17", "AFC_60_19"))) %>%
  group_by(pop, line, plate, Indiv, CHROM_schf, crossover) %>%
  tally %>% 
  spread(key = crossover, value = n) %>%
  rename(crossovers = `TRUE`) %>%
  mutate(crossovers = ifelse(is.na(crossovers), 0, crossovers)) %>%
  select(-`FALSE`) %>%
  ungroup

write.table(co_counts, "data/co_counts_manual.txt", row.names = FALSE, quote = FALSE)

co_counts_mean <- co_counts %>%
  group_by(pop, line) %>%
  summarise(mean_co = mean(crossovers))

write.table(co_counts_mean, "data/fasta_recomb/summaries/co_rate_line.txt", col.names = T, quote = FALSE)

###########################################################
# Plots: Chromosome-level
###########################################################

theme_presentation<- function(base_size = 24, title_size = 28, legend_size = 18, 
                              bg_col = "grey25", panel_col = "white", text_col = "white",
                              grid_col = "grey50",
                              base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = base_size,colour=text_col),
      strip.text.y = element_text(size = base_size,colour=text_col),
      axis.text.x = element_text(size= base_size,colour=text_col),
      axis.text.y = element_text(size= base_size,colour=text_col),
      axis.ticks =  element_line(colour = text_col), 
      axis.title.x= element_text(size=title_size,colour=text_col),
      axis.title.y= element_text(size=title_size,angle=90,colour=text_col),
      #legend.position = "none", 
      panel.background = element_rect(fill=panel_col, color = NULL), 
      panel.border =element_blank(), 
      panel.grid = element_line(color = grid_col, size = 0.5),
      panel.grid.major.x = element_blank(), 
      panel.grid.minor.x = element_blank(), 
      panel.grid.minor.y = element_blank(), 
      panel.margin = unit(1.0, "lines"), 
      plot.background = element_rect(fill=bg_col, color = NULL), 
      plot.title = element_text(size=base_size,colour=text_col), 
      plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
      axis.line = element_line(colour = text_col),
      legend.background=element_rect(fill=bg_col),
      legend.title = element_text(size=legend_size,colour=text_col),
      legend.text = element_text(size=legend_size,colour=text_col),
      legend.key = element_rect( fill = bg_col),
      legend.key.size = unit(c(0.5, 0.5), "lines")
    )
}

# for ordering plots below
co_counts_mean <- co_counts %>%
  mutate(cross_chr = paste0(line, "_", CHROM_schf)) %>%
  mutate(cross_chr = factor(cross_chr)) %>%
  group_by(pop, line, CHROM_schf, cross_chr) %>%
  summarise(mean_co = mean_cl_boot(crossovers)$y, 
            upper_ci = mean_cl_boot(crossovers)$ymax,
            lower_ci = mean_cl_boot(crossovers)$ymin) %>%
  mutate(cross_chr = reorder(cross_chr, mean_co))

# MEAN COUNTS PER LINE 
# ALL CHROMOSOMES
# (with error bars)
(co_line <- co_counts %>%
  filter(crossovers <= 10) %>%
  filter(!(plate %in% c("AFC_19_17", "AFC_60_19"))) %>%
  ungroup %>%
  mutate(cross_chr = paste0(line, "_", CHROM_schf)) %>%
  left_join(co_counts_mean) %>%
  mutate(mean_co = ifelse(pop == "MC", mean_co + 100, mean_co)) %>% 
  mutate(cross_chr = factor(cross_chr)) %>%
  mutate(cross_chr = reorder(cross_chr, mean_co)) %>%
  ggplot(aes(x = as.factor(CHROM_schf), y = crossovers, color = pop, group = cross_chr, label = line))+
  #geom_boxplot(position = position_dodge(width = 0.75))+
  stat_summary(fun.data = "mean_cl_boot", size = 0.5, position = position_dodge(width = 0.825))+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Mean number of crossovers\n")+
  xlab("\nChromosome")+
  theme_presentation(base_size = 22, title_size = 24, legend_size = 18, 
                     panel_col = "grey25", bg_col = "black", grid_col = "grey50") +
  scale_color_brewer(palette = "Set1", name = "Population"))


read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(co_line), left = 0, top = 0, width = 10, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  print(target = "slides/co_rates_chromosome.pptx")


# MEAN COUNTS PER POP
(co_pop <- co_counts %>%
  filter(crossovers <= 10) %>%
  ggplot(aes(x = as.factor(CHROM_schf), y = crossovers, color = pop))+
  #geom_point()+
  stat_summary(fun.data = "mean_cl_boot", size = 1.25, position = position_dodge(width = 0.5))+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Mean number of crossovers\n")+
  xlab("\nChromosome")+
  theme_presentation(base_size = 22, title_size = 24, legend_size = 18, 
                     panel_col = "grey25", bg_col = "black", grid_col = "grey50") +
  scale_color_brewer(palette = "Set1", name = "Population"))


read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(co_pop), left = 0, top = 0, width = 10, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  print(target = "slides/co_rates_chromosome_pop.pptx")

###########################################
#  local variation in CO rate
###########################################

# re-inflate the segmented data
# (this allows for more precise filtering of short ancestry blocks, i.e. to detect false double XOs)
inf_df <- seg_df %>%
  mutate(interval_length = end_rank - start_rank) %>%
  filter(interval_length >= 2) %>%
  group_by(Indiv, CHROM_schf, interval_id) %>%
  do(rqtl_rank = .$start_rank:.$end_rank, 
     gt_adj = rep(.$gt_adj, length(.$start_rank:.$end_rank)))

inf_df <- inf_df %>% 
  unnest 

# extract the position data for all intervals
# read in rQTL marker orders
ord_df <- read.table("data/rqtl/marker_orders/all.txt", h = T) %>%
  mutate(CHROM_schf = gsub("*._chr_", "", marker) %>% 
           gsub("_pos.*", "", .)) %>%
  mutate(POS_schf = as.numeric(gsub(".*_pos_", "", marker))) %>%
  group_by(chr) %>%
  mutate(rqtl_rank = dense_rank(map_pos)) %>%
  ungroup %>%
  rename(rqtl_chr = chr, rqtl_map_pos = map_pos) %>%
  select(CHROM_schf, POS_schf, rqtl_chr, rqtl_map_pos, rqtl_rank) 

# compute inter-marker distnaces (for filtering out redundant markers)
ord_df <- ord_df %>%
  arrange(CHROM_schf, POS_schf) %>%
  group_by(CHROM_schf) %>%
  mutate(nn_dist_prev = abs(POS_schf - lag(POS_schf))) %>%
  mutate(nn_dist_next = abs(POS_schf - lead(POS_schf))) %>%
  mutate(nn_dist_adj = ifelse(is.na(nn_dist_prev)|nn_dist_prev > nn_dist_next, nn_dist_next, nn_dist_prev )) %>%
  select(-nn_dist_prev, -nn_dist_next)

ord_df %>%
  filter(nn_dist_adj < 1000) %>%

# join in chromosomal coordinates
inf_df <- inf_df %>%
  left_join(ord_df) 

seg_df %>%
filter(Indiv == "AFC_56_16_A11")


# a marker-level list of all the crossover events
co_df <- seg_df %>%
  filter(len_seg >= 2) %>%
  group_by(Indiv, CHROM_schf) %>%
  mutate(crossover = ifelse(gt_adj != lag(gt_adj), TRUE, FALSE)) %>%
  mutate(crossover = ifelse(is.na(crossover), FALSE, crossover)) %>%
  filter(crossover) %>%
  mutate(rqtl_rank = start_rank) %>% 
  select(Indiv, CHROM_schf, start_POS, end_POS, rqtl_rank, crossover)

co_df <- co_df %>% 
  unnest

# annotate the inf_df with the crossover information
inf_df <- left_join(inf_df, co_df) %>%
  mutate(crossover = ifelse(is.na(crossover), FALSE, crossover))

# for checking xo counts visually
xo_sum_df <- inf_df %>%
  group_by(Indiv, CHROM_schf) %>%
  summarise(xo_count = sum(crossover))

# CHECK
# plot the inflated genotypes

inf_df %>%
  left_join(xo_sum_df) %>%
  filter(grepl("AFC", Indiv)) %>%
  filter(CHROM_schf == "XR") %>%
  ggplot(aes(x = as.factor(rqtl_rank), y = Indiv, 
             fill = as.factor(gt_adj)))+
  geom_raster()+
  theme(strip.text.y = element_text(angle = 180),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_discrete(position = "right")+
  facet_wrap(~xo_count, scales = "free_y")+
  xlab("")+
  scale_fill_brewer(palette = "Set1")

# line-level summary
co_rates_line <- inf_df %>%
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  mutate(pop = gsub("_.*$", "", Indiv)) %>%
  group_by(pop, line, CHROM_schf, POS_schf, rqtl_rank) %>%
  summarise(co_rate = mean(crossover, na.rm = TRUE), 
            co_rate_n = sum(!is.na(crossover))) %>%
  mutate(co_rate_sd = sqrt(co_rate * (1 - co_rate) / co_rate_n)) %>%
  ungroup

co_rates_pop <- inf_df %>%
  mutate(line = gsub("_[0-9]{2}_[A-H][0-9]{2}$", "", Indiv))%>%
  mutate(pop = gsub("_.*$", "", Indiv)) %>%
  group_by(pop, CHROM_schf, POS_schf, rqtl_rank) %>%
  summarise(co_rate = mean(crossover, na.rm = TRUE), 
            co_rate_n = sum(!is.na(crossover))) %>%
  mutate(co_rate_sd = sqrt(co_rate * (1 - co_rate) / co_rate_n)) %>%
  ungroup

co_rate_diff <- co_rates_pop %>%
  select(-co_rate_n, -co_rate_sd) %>%
  spread(key = pop, value = co_rate) %>%
  mutate(co_diff = AFC - MC)

###########################################
# PLOTS: local variation in CO rate
###########################################

# by line
(local_rates_pop_marker <- co_rates_pop %>% 
  mutate(CHROM_schf = factor(CHROM_schf, levels = c("2", "3", "4", "XR", "XL"))) %>%
  mutate(POS_schf_Mb = POS_schf / 1000000) %>%
  mutate(POS_schf_Mb = ifelse(CHROM_schf == "XR", (POS_schf_Mb - max(POS_schf_Mb)) * -1, POS_schf_Mb)) %>%
  ggplot(aes(x = rqtl_rank, 
             y = as.numeric(co_rate), 
             group = pop,
             color = pop))+
  geom_step(size = 0.75)+
  geom_hline(yintercept = 0, color = "white")+
  facet_wrap(~CHROM_schf, scales = "free_x")+
  theme_presentation(base_size = 20, title_size = 22, legend_size = 20, 
                     panel_col = "grey15", bg_col = "black", grid_col = "grey40") +
  ylab("Mean Crossover Rate\n")+
  xlab("\nMarker Index")+
  scale_color_brewer(palette = "Set1", name = "Line"))

read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(local_rates_pop_marker), left = 0, top = 0, width = 10, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  print(target = "slides/co_rates_local_marker.pptx")

(local_rates_pop_cm <- co_rates_pop %>% 
    mutate(CHROM_schf = factor(CHROM_schf, levels = c("2", "3", "4", "XR", "XL"))) %>%
    mutate(POS_schf_Mb = POS_schf / 1000000) %>%
    mutate(POS_schf_Mb = ifelse(CHROM_schf == "XR", (POS_schf_Mb - max(POS_schf_Mb)) * -1, POS_schf_Mb)) %>%
    ggplot(aes(x = POS_schf_Mb, 
               y = as.numeric(co_rate), 
               group = pop,
               color = pop))+
    geom_step(size = 0.75)+
    geom_hline(yintercept = 0, color = "white")+
    facet_wrap(~CHROM_schf, scales = "free_x")+
    theme_presentation(base_size = 20, title_size = 22, legend_size = 20, 
                       panel_col = "grey15", bg_col = "black", grid_col = "grey40") +
    ylab("Mean Crossover Rate\n")+
    xlab("\nMarker Index")+
    scale_color_brewer(palette = "Set1", name = "Line"))

read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(local_rates_pop_cm), left = 0, top = 0, width = 10, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  print(target = "slides/co_rates_local_cm.pptx")

################################################
# PLOTS: correlation in genome-wide recomb rate
################################################
library("GGally")

my_bin <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point()+
    geom_smooth(method = "lm", se = TRUE, fill = "grey75", alpha = 0.25, color = "white")+
    scale_color_brewer(palette = "Set1", name = "Line")
}

(co_cor <- co_counts %>%
  mutate(CHROM_schf = paste0("chr", CHROM_schf)) %>%
  group_by(pop, line, CHROM_schf) %>%
  summarise(mean_xo_rate = mean(crossovers)) %>%
  spread(key = CHROM_schf, value = mean_xo_rate) %>%
  ggpairs(data = ., columns = c("chr2", "chr3", "chr4", "chrXL", "chrXR"),
          color = "pop",
          mapping=ggplot2::aes(colour = pop),
          diag = list(continuous = "blank"),
          lower = list(continuous = my_bin ),
          upper = list(continuous = "blank"))+
  theme_presentation(base_size = 14, title_size = 16, legend_size = 16, 
                     panel_col = "grey15", bg_col = "black", grid_col = "grey40"))

read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(co_cor), left = 0, top = 0, width = 10, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  print(target = "slides/co_rates_correlation.pptx")

(co_cor2 <- co_counts %>%
    mutate(CHROM_schf = paste0("chr", CHROM_schf)) %>%
    group_by(pop, line, CHROM_schf) %>%
    summarise(mean_xo_rate = mean(crossovers, na.rm = TRUE)) %>% 
    spread(key = CHROM_schf, value = mean_xo_rate) %>%
    ggplot(aes(x = chrXR, y = chr3, color = pop))+
    geom_point(size = 3)+
    geom_smooth(method = "lm", color = "white", fill = "grey75", alpha = 0.10, size = 2)+
    theme_presentation(base_size = 14, title_size = 16, legend_size = 16, 
                       panel_col = "grey15", bg_col = "black", grid_col = "grey40")+
      ylab("Chr3 Crossover Rate\n")+
      xlab("\nChrXR Crossover Rate")+
  scale_color_brewer(palette = "Set1", name = "Population"))

read_pptx() %>% 
  add_slide(layout = "Title and Content", master = "Office Theme") %>% 
  ph_with_vg_at(code = print(co_cor2), left = 0, top = 0, width = 10, height = 7.5) %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  print(target = "slides/co_rates_correlation2.pptx")


# pop values
co_rates_pop %>%
  ggplot(aes(x = rqtl_rank, 
             y = as.numeric(co_rate), 
             group = pop,
             color = pop))+
  geom_step(size = 0.75)+
  facet_wrap(~CHROM_schf, scales = "free_x")+
  theme_classic()+
  ylab("Mean Crossover Rate")+
  xlab("Marker Rank")+
  scale_color_brewer(palette = "Set1", name = "Population")

# pop values
co_rate_diff %>%
  ggplot(aes(x = AFC, 
             y = MC))+
  geom_point()+
  geom_smooth(method = "lm") +
  theme_classic()

# pop differences
co_rate_diff  %>%
  ggplot(aes(x = rqtl_rank, 
             y = as.numeric(co_diff)))+
  geom_step(size = 0.5)+
  facet_wrap(~CHROM_schf, scales = "free_x")+
  theme_classic()+
  ylab("Mean Crossover Rate (+ve = higher in AFC)")+
  xlab("Marker Rank")+
  geom_hline(yintercept  = 0, lty = 3) +
  scale_color_brewer(palette = "Set1", name = "Population")


co_counts %>%
  filter(!(plate %in% c("AFC_19_17", "AFC_60_19"))) %>%
  group_by(pop, line, CHROM_schf) %>%
  summarise(mean_crossovers = mean(crossovers, na.rm = TRUE)) %>%
  spread(key = CHROM_schf, value = mean_crossovers) %>%
  rowwise %>%
  mutate(genome_wide = mean(c(`2`, `3`, `4`, XL, XR))) %>%
  arrange(XL)
filter(line %in% c("AFC_57", "AFC_48", "MC_27", "AFC_30"))

AFC48 (h) = 0.662
AFC30 (h) = 0.697
MC27 (l) = 0.679 
AFC57 (l) = 0.703


# recomb correlation plot
tmp <- co_counts_mean %>%
  spread(key = CHROM_reorder, value = mean_co)

names(tmp)[3:4] <- c("chr2", "chr3")

lm(tmp$chr2 ~ tmp$chr3) %>% summary

recomb_cor <- tmp %>%
  ggplot(aes(x = chr2, y = chr3))+
  geom_point(color = "white", size = 4) +
  geom_smooth(method = "lm", size = 2) +
  theme_presentation()+
  ylab("Chromosome 3 Crossovers")+
  xlab("Chromosome 2 Crossovers")




