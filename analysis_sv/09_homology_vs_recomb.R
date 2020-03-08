# interset homology data with recombination data

# summarize sv data for each window
# KMS 2019

library("tidyverse")
library("parallelsugar")
library("broom")

#options(device = "quartz")

############################################
# read in summarized recombination data
############################################

recomb_sum_df <- readRDS("analysis_sv/data/sv_window_summary.rds") %>%
  ungroup %>%
  rename(POS_schf1 = pos1, POS_schf2 = pos2, line = Indiv)

recomb_df <- readRDS("data/crossovers_fine_scale.rds")

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

# compute the window-based summaries (counts of COs)
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
  filter(abs(POS_schf1 - POS_schf2) < 1500000) %>%
  ungroup %>%
  rename(Indiv = line)

str(recomb_sum_df )

############################################
# read in homology data
############################################

il_hm_df <- readRDS("analysis_sv/data/homology/illumina_homology.rds") %>%
  rename(pos1_ref = pos1, pos2_ref = pos2) %>%
  mutate(pos1_ref = ifelse(grepl("group", vcf_file), gsub(".*group[a-z0-9]{1,4}_", "", vcf_file) %>% gsub("_.*", "", .), pos1_ref)) %>%
  mutate(pos2_ref = ifelse(grepl("group", vcf_file), gsub(".*group[a-z0-9]{1,4}_", "", vcf_file) %>% gsub("_recomb_regions.*", "", .) %>% gsub(".*_", "", .), pos2_ref)) %>%
  mutate(pos1_ref = as.numeric(pos1_ref)) %>%
  mutate(pos2_ref = as.numeric(pos2_ref))

recomb_intervals <- read.table("analysis_sv/meta/recom_intervals_ref_schf.txt", h = T) %>%
  rename(pos1_ref = pos1, pos2_ref = pos2)

il_hm_df <- il_hm_df %>%
  left_join(recomb_intervals, by = c("pos1_ref", "pos2_ref")) 
  
il_hm_df <- il_hm_df %>%
  mutate(POS_schf2 = POS_schf2 - 1) %>%
  select(Indiv, CHROM_schf, POS_schf1, POS_schf2, var_class, difference_count, gc_gain, gc_loss, n_sites)

il_hm_df <- left_join(il_hm_df, recomb_sum_df) 


il_hm_df %>%
  filter(n_sites > 10000)  %>%
  mutate(homology_est  = difference_count / n_sites) %>%
  group_by(var_class, CHROM_schf, POS_schf1) %>%
  do(homology_lm = glm(.$crossover_count ~ .$homology_est, family = "poisson")) %>%
  tidy(homology_lm) %>%
  filter(term != "(Intercept)")

  
il_hm_df %>%
  filter(var_class == "indel") %>%
  filter(n_sites > 5000)  %>%
  mutate(homology_est  = difference_count / n_sites) %>%
  ggplot(aes(x = homology_est, y = crossover_count))+
  geom_point()

il_hm_df %>%
  filter(var_class == "snp") %>%
  filter(n_sites > 5000)  %>%
  mutate(homology_est  = difference_count / n_sites) %>%
  ggplot(aes(x = homology_est, y = crossover_count))+
  geom_point()


library("lme4")
library("visreg")

mod1 <- il_hm_df %>%
  filter(var_class == "indel") %>%
  filter(n_sites > 10000) %>%
  mutate(homology_est  = difference_count / n_sites) %>%
  mutate(window = paste(CHROM_schf, POS_schf1, POS_schf2, sep = "_")) %>% 
  glmer(crossover_count ~ homology_est + (1|Indiv) + (1 | window) ,
        data = ., family = "poisson")

mod2 <- il_hm_df %>%
  filter(var_class == "indel") %>%
  filter(n_sites > 10000) %>%
  mutate(homology_est  = difference_count / n_sites) %>%
  mutate(window = paste(CHROM_schf, POS_schf1, POS_schf2, sep = "_")) %>% 
  glmer(crossover_count ~1 + (1|Indiv) + (1 | window) ,
        data = ., family = "poisson")


anova(mod2, mod1)

summary(mod1)
visreg(mod1)
lmerTest::ranova(mod1)

car::Anova(mod1)

mod2 <- il_hm_df %>%
  filter(var_class == "snp") %>%
  filter(n_sites > 5000) %>%
  mutate(homology_est  = difference_count / n_sites) %>%
  mutate(window = paste(CHROM_schf, POS_schf1, POS_schf2, sep = "_")) %>% 
  glmer(crossover_count ~ homology_est + (1|Indiv) + (1 | window),
        data = ., family = "poisson")

summary(mod2)

car::Anova(mod2)
install.packages("rlm")
library("rlm")

il_hm_df %>%
  filter(n_sites > 10000) %>%
  mutate(gc_gain  = gc_gain / n_sites) %>%
  mutate(gc_loss  = gc_gain / n_sites) %>%
  #filter(gc_loss < 1.5e-07) %>%
  mutate(homology_est  = difference_count / n_sites) %>%
  ggplot(aes(x = homology_est, y = crossover_count))+
  geom_point(aes(color = gc_loss))+
  stat_smooth(method = "lm", size = 1, color = "red")+
  facet_wrap(~var_class, scales = "free_x")+
  scale_color_viridis_c()

il_hm_df %>%
  filter(n_sites > 5000) %>%
  filter(var_class == "indel") %>%
  mutate(gc_gain_prop  = gc_gain / difference_count) %>%
  mutate(gc_loss_prop  = gc_loss / difference_count) %>%
  #filter(gc_loss < 1.5e-07) %>%
  mutate(homology_est  = (n_sites - difference_count) / n_sites) %>%
  ggplot(aes(y = crossover_count, x = homology_est))+
  geom_point()+
  stat_smooth()+
  scale_color_viridis_c()


il_hm_df %>%
  filter(n_sites > 5000) %>%
  filter(var_class == "snp") %>%
  mutate(gc_gain  = gc_gain / n_sites) %>%
  mutate(gc_loss  = gc_gain / n_sites) %>%
  ggplot(aes(x = gc_gain, y = crossover_count))+
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, color = "red")+
  facet_wrap(~var_class, scales = "free_x")

il_hm_df %>%
  filter(n_sites > 100000) %>%
  filter(var_class == "snp") %>%
  mutate(gc_gain  = gc_gain / n_sites) %>%
  mutate(gc_loss  = gc_gain / n_sites) %>%
  filter(gc_loss < 1.5e-07) %>%
  mutate(window = paste(CHROM_schf, POS_schf1, POS_schf2, sep = "_")) %>% 
  ggplot(aes(x = gc_gain, y = crossover_count))+
  geom_point()+
  stat_smooth(method = "lm")+
  facet_wrap(~window, scales = "free")+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())


mod2 <- il_hm_df %>%
  filter(var_class == "snp") %>%
  filter(n_sites > 5000) %>%
  mutate(gc_gain  = gc_gain / n_sites) %>%
  mutate(window = paste(CHROM_schf, POS_schf1, POS_schf2, sep = "_")) %>% 
  glmer(crossover_count ~ gc_gain + (1 | window),
        data = ., family = "poisson")

summary(mod2)
car::Anova(mod2)


mod2 <- il_hm_df %>%
  filter(var_class == "snp") %>%
  filter(n_sites > 5000) %>%
  mutate(gc_loss  = gc_loss / n_sites) %>%
  mutate(window = paste(CHROM_schf, POS_schf1, POS_schf2, sep = "_")) %>% 
  glmer(crossover_count ~ gc_loss + (1 | window),
        data = ., family = "poisson")

summary(mod2)
car::Anova(mod2)
