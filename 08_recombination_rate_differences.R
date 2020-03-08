#install.packages("effects")

library("tidyverse")
library("parallel")
library("lme4")
library("wesanderson")
library("GGally")
library("optimx")
library("patchwork")
library("visreg")
library("broom")
library("rsq")
source("functions/theme_presentation.R")

#options(device = "quartz")
#options(device = "cairo_pdf")

select <- dplyr::select

###########################################################
# read in recombination rate data
###########################################################

pop_labels <- read.table("data/co_counts_manual.txt", h = T) %>%
  select(pop, line, plate, Indiv)

# read in an format crossover estimates from r/qtl
# also perform filtering
co_counts_rqtl <- read.table("data/rqtl/co_estimates/all.txt", h = T) %>%
  rename(Indiv = id) %>%
  left_join(pop_labels)%>%
  filter(!grepl("AFC_19_17|AFC_60_19", Indiv)) %>% # these plates accidently recieved two plate barcodes
  filter(crossovers <= 5) %>%
  distinct %>%
  mutate(pop = as.factor(pop), line = as.factor(line)) %>%
  group_by(pop, line, Indiv) %>%
  summarise(crossovers = sum(crossovers, na.rm = TRUE)) %>%
  ungroup

# hypo/hypergamous check
co_counts_chr_wide <- read.table("data/rqtl/co_estimates/all.txt", h = T) %>%
  rename(Indiv = id) %>%
  left_join(pop_labels)%>%
  filter(!grepl("AFC_19_17|AFC_60_19", Indiv)) %>%
  filter(crossovers <= 5) %>%
  distinct %>%
  select(pop, line, plate, Indiv, everything()) %>%
  mutate(chr = paste0("chr_", chr)) %>%
  spread(key = chr, value = crossovers)

# write out crossovers counts for suppmat
setNames(co_counts_chr_wide, 
         c("Population", "Line", "Plate", "Individual", 
           "Chr2_Crossovers", "Chr3_Crossovers", "Chr4_Crossovers", "ChrX_Crossovers")) %>%
  write.csv(file = "figures/TableS1.csv", row.names = FALSE, quote = FALSE)
  

co_counts_fine <- readRDS("data/crossovers_fine_scale.rds")

##########################################################################
# FIGURE 1A + 1B
##########################################################################

# perform statistical comparison of crossover rates

# poisson glmers 
co_lm_pois_rqtl1 <- lme4::glmer(crossovers ~ (1|line),  data = co_counts_rqtl, family = "poisson")
co_lm_pois_rqtl2 <- lme4::glmer(crossovers ~ pop + (1|line),  data = co_counts_rqtl, family = "poisson")

# type II wald chisquare test
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: crossovers
# Chisq Df Pr(>Chisq)  
# pop 5.5256  1    0.01874 *
car::Anova(co_lm_pois_rqtl2)

# likeihood ratio test (comparing random effect only vs. random effect + population)
# Data: co_counts_rqtl
# Models:
# co_lm_pois_rqtl1: crossovers ~ (1 | line)
# co_lm_pois_rqtl2: crossovers ~ pop + (1 | line)
# Df   AIC   BIC logLik deviance Chisq Chi Df Pr(>Chisq)  
# co_lm_pois_rqtl1  2 32849 32862 -16422    32845                          
# co_lm_pois_rqtl2  3 32846 32867 -16420    32840 4.794      1    0.02856 *
anova(co_lm_pois_rqtl1, co_lm_pois_rqtl2, test = "LRT")

#co_lm_pois_rqtl2 <- glm(crossovers ~ pop,  data = co_counts_rqtl, family = "poisson")

df <- augment(co_lm_pois_rqtl2) %>%
  mutate(pearson_resid = (.resid/sd(.resid)))

figA <- df %>%
  ggplot(aes(sample = pearson_resid))+  ## scale to variance=1
  stat_qq(distribution = stats::qpois, 
          dparams = list(lambda = mean(df$pearson_resid)), 
          size = 2, alpha = 0.1)+
  stat_qq_line(distribution = stats::qpois, 
               dparams = list(lambda = mean(df$pearson_resid)), lty = 2,
               detrend = TRUE)+
  xlab("Theoretical Poisson Quantiles")+
  ylab("Sample Poisson Quantiles")+
  theme_bw()

figB <- df %>%
  ggplot(aes(x = .fitted, y = pearson_resid - 10))+
  geom_point(size = 2, alpha = 0.1)+
  geom_hline(yintercept = 0, lty = 2)+
  xlab("Fitted Value")+
  ylab("Pearson Residual")+
  theme_bw()

sup_fig <- (figA|figB) + plot_annotation(tag_levels = 'A')
ggsave(filename = "figures/FigureS3.pdf", sup_fig, width = 8, height = 4, useDingbats = FALSE)


# compute confidence intervals (for plotting only)
# gaussian confint (rqtl)
co_lm_gaus_rqtl <- lme4::lmer(crossovers ~ pop + (1|line),  data = co_counts_rqtl)
co_lm_gaus_rqtl <- lme4::lmer(crossovers ~ pop + 0 + (1|line),  data = co_counts_rqtl)
sum_lm <- summary(co_lm_gaus_rqtl)$coefficients[,1]
means <- sum_lm %>% as.numeric
cis <- confint(co_lm_gaus_rqtl, level = 0.90)
ci_upper <- cis[,2][3:4]
ci_lower <- cis[,1][3:4]

# confidence intervals for population means
co_pop_rqtl_df <- data.frame(pop = c("AFC", "MC"), mean_co = means, ci_upper, ci_lower)

# palette 
plot_pal <- wes_palette("Darjeeling1", n = 2, type = "discrete")%>% as.character

# for ordering plots below
co_counts_mean_rqtl <- co_counts_rqtl %>%
  filter(crossovers <= 10) %>%
  group_by(pop, line) %>%
  summarise(mean_co = mean_cl_boot(crossovers)$y, 
            upper_ci = mean_cl_boot(crossovers)$ymax,
            lower_ci = mean_cl_boot(crossovers)$ymin)

mean_co_new <- read.table("data/rqtl/co_estimates/all.txt", h = T) %>%
  rename(Indiv = id) %>%
  left_join(pop_labels)%>%
  filter(!grepl("AFC_19_17|AFC_60_19", Indiv)) %>%
  filter(crossovers <= 5) %>%
  distinct %>% 
  group_by(pop, line) %>%
  summarise(mean_co = mean_cl_boot(crossovers)$y, 
            upper_ci = mean_cl_boot(crossovers)$ymax,
            lower_ci = mean_cl_boot(crossovers)$ymin)

# summary stats (from raw counts, not from models)
co_counts_mean_rqtl %>% pull(mean_co) %>% summary

afc_recomb <- mean_co_new %>%
  filter(pop == "AFC") %>%
  pull(mean_co)

mc_recomb <- mean_co_new %>%
  filter(pop == "MC") %>%
  pull(mean_co)

t.test(afc_recomb, mc_recomb)
wilcox.test(afc_recomb, mc_recomb)



#### plots

# MEAN COUNTS PER LINE 
# (with error bars)
fig_1a <- co_counts_rqtl %>%
    filter(crossovers <= 10) %>%
    left_join(co_counts_mean_rqtl) %>%
    mutate(line = gsub("_", " ", line)) %>%
    mutate(line = fct_reorder(line, mean_co)) %>%
    ggplot(aes(x = as.factor(line), y = crossovers, group = pop, color = pop))+
    stat_summary(fun.data = "mean_cl_boot", size = 0.75)+
    theme_bw()+
    theme(legend.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        text = element_text(size=14),
        axis.text.x = element_text(angle = 45, size = 10, hjust = 1))+
    ylab("Crossovers per genome")+
    xlab("Inbred Line")+
    scale_color_manual(values = plot_pal)+
    coord_cartesian( ylim= c(4.10, 6.0))

# global recombination rate between populations
fig_1b <- co_counts_mean_rqtl %>%
  ggplot(aes(x = as.factor(pop), y = mean_co, group = pop, color = pop))+
  geom_jitter(width = 0.08, alpha = 0.4, size = 2, shape = 16)+
  geom_pointrange(data = co_pop_rqtl_df, aes(x = pop, ymin = ci_lower, ymax = ci_upper), size = 0.8, lwd = 0)+
  xlab("Population")+
  ylab("")+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        text = element_text(size=14),
        axis.text.x = element_text(size = 12))+
  scale_color_manual(values = plot_pal)+
  coord_cartesian( ylim= c(4.10, 6.0))

fig_1b

(fig_1a | fig_1b)
##########################################################################
# FIGURE 1C
##########################################################################

# drop windows with too much missing data
geno_df <- co_counts_fine %>%
  mutate(pop = gsub("_.*", "", Indiv )) %>%
  mutate(line = gsub("_[0-9]{2}_[A-Z,0-9]{3}$", "", Indiv )) %>%
  ungroup %>%
  group_by(CHROM_schf, POS_schf1,POS_schf2) %>%
  summarise(n_genotyped = length(crossover_int)) %>%
  ungroup()

############# PER POPULATION

co_pop <- co_counts_fine %>%
  mutate(pop = gsub("_.*", "", Indiv )) %>%
  mutate(line = gsub("_[0-9]{2}_[A-Z,0-9]{3}$", "", Indiv )) %>%
  ungroup %>%
  left_join(geno_df) %>%
  filter(n_genotyped > 4000) %>%
  group_by(pop,  CHROM_schf, POS_schf1,POS_schf2) %>%
  summarise(co_count = sum(crossover_int, na.rm = TRUE), co_n = length(crossover_int)) %>%
  ungroup %>%
  mutate(cM = 50*log((1/(1-(2*co_count / co_n))))) %>%
  mutate(Mb = abs(POS_schf1 - POS_schf2) / 1000000) %>%
  mutate(cM_Mb = cM / Mb) %>%
  group_by(CHROM_schf) %>%
  mutate(window_rank = dense_rank(POS_schf1)) %>%
  ungroup %>%
  filter(cM_Mb < 50)



############## PER LINE

co_line <- co_counts_fine %>%
  mutate(pop = gsub("_.*", "", Indiv )) %>%
  mutate(line = gsub("_[0-9]{2}_[A-Z,0-9]{3}$", "", Indiv )) %>%
  ungroup %>%
  left_join(geno_df) %>%
  filter(n_genotyped > 4000) %>%
  group_by(pop, line, CHROM_schf, POS_schf1,POS_schf2) %>%
  summarise(co_count = sum(crossover_int, na.rm = TRUE), co_n = length(crossover_int)) %>%
  ungroup %>%
  mutate(cM = 50*log((1/(1-(2*(co_count / co_n)))))) %>%
  mutate(Mb = abs(POS_schf1 - POS_schf2) / 1000000) %>%
  mutate(cM_Mb = cM / Mb) %>%
  group_by(CHROM_schf) %>%
  mutate(window_rank = dense_rank(POS_schf1)) %>%
  ungroup %>%
  filter(cM_Mb < 50)

fig_1c <- co_line %>%
  ggplot(aes(x = POS_schf1/1000000, y = cM_Mb, color = line))+
  geom_step(alpha = 0.2) +
  geom_step(data = co_pop, aes(color = pop), size = 0.75) +
  facet_grid(pop~CHROM_schf, space = "free_x", scales = "free_x")+
  theme_bw()+
  xlab("Position (MB)")+
  ylab("Recombination Rate (cM/MB)") +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        text = element_text(size=14))+
  scale_color_manual(values = c(rep(plot_pal[1], 11), rep(plot_pal[2],8)))

# compose figure 1
figure1 <- (fig_1a | fig_1b) / fig_1c + plot_layout(nrow = 2, heights=c(2,3))
ggsave(filename = "figures/raw/Figure1_raw.pdf", figure1, width = 10, height = 8, useDingbats = FALSE)

##########################################################################
# FIGURE 2
##########################################################################

# Figure 2a: correlation between population maps

# test of correlation

mc_recomb <- co_pop %>%
  filter(pop == "MC") %>%
  select(CHROM_schf, POS_schf1, co_count, co_n) %>%
  mutate(MC_recomb = co_count/co_n)%>%
  rename(co_count_MC = co_count, co_n_MC = co_n)

afc_recomb <- co_pop %>%
  filter(pop == "AFC") %>%
  select(CHROM_schf, POS_schf1, co_count, co_n) %>%
  mutate(AFC_recomb = co_count/co_n) %>%
  rename(co_count_AFC = co_count, co_n_AFC = co_n)


# join and create a
# create a per-interval id variable
recomb_both <- left_join(afc_recomb, mc_recomb)%>%
  mutate(interval_id = paste0(CHROM_schf, "_", POS_schf1)) %>%
  mutate(recomb_diff = AFC_recomb - MC_recomb)

global_diff <- recomb_both %>%
  summarise(co_AFC = sum(co_count_AFC), n_AFC = sum(co_n_AFC),
            co_MC = sum(co_count_MC), n_MC = sum(co_n_MC)) %>%
  unlist()

# compare each interval between the two populations
interval_diff <- recomb_both %>%
  group_by(interval_id) %>%
  do(recomb_diff = prop.test(x = c(.$co_count_AFC, .$co_count_MC),
                             n = c(.$co_n_AFC, .$co_n_MC)) %>% tidy) %>%
  unnest %>%
  select(-statistic, -parameter, -method, -alternative) 

# same for lines
interval_diff_line <- co_line %>%
  mutate(interval_id = paste0(CHROM_schf, "_", POS_schf1)) %>%
  mutate(co_prop = co_count/co_n) %>%
  group_by(interval_id) %>%
  do(recomb_diff = glm(.$co_prop~ .$line, weights = .$co_n, family = "binomial") %>% car::Anova(.) %>% tidy) %>%
  unnest %>%
  filter(!is.na(df)) %>%
  select(interval_id, p.value)

interval_diff_var <- co_line %>%
  mutate(interval_id = paste0(CHROM_schf, "_", POS_schf1)) %>%
  mutate(co_prop = co_count/co_n) %>%
  group_by(interval_id) %>%
  #do(recomb_diff = glm(.$co_prop~.$line, weights = .$co_n, family = "binomial")  %>% tidy) %>%
  do(var_line = var(.$co_prop)) %>%
  unnest 

interval_diff_line <- left_join(interval_diff_line, interval_diff_var)

# add in Mb for each interval
mb_dat <- co_pop %>% select(CHROM_schf, POS_schf1, Mb) %>%
  mutate(interval_id = paste0(CHROM_schf, "_", POS_schf1))

interval_diff <- left_join(interval_diff, mb_dat) %>%
  mutate(co_diff = estimate1 - estimate2) %>%
  distinct()

interval_diff_line <- left_join(interval_diff_line, mb_dat)

# comparison of population-averaged CO rates
lm(data = recomb_both , AFC_recomb~MC_recomb)
cor.test(recomb_both$AFC_recomb, recomb_both$MC_recomb)

# figure 2a 
figure2a <- recomb_both %>%
  ggplot(aes(x = AFC_recomb, y = MC_recomb)) +
  geom_point(size = 2, alpha = 0.2, pch = 16)+
  geom_smooth(method = "lm", se = TRUE, size = 1.5)+
  theme_bw()+
  ylab("Crossover Frequency (MC)")+
  xlab("Crossover Frequency (AFC)")+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 14))

# local differences in recombination rate

Figure_S3 <- interval_diff %>%
  filter(!is.nan(p.value)) %>%
  mutate(p_adjust = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = p_adjust < 0.05) %>%
  ggplot(aes(x = POS_schf1/1000000, y = co_diff/Mb, 
             ymin = conf.low/Mb, ymax = conf.high/Mb, fill = p_adjust, color = significant)) +
  geom_pointrange(size = 0.5, pch = 21)+
  geom_hline(yintercept = 0, linetype = 3)+
  facet_wrap(~CHROM_schf, scales = "free_x")+
  theme_bw()+
  ylab("Difference in CO Rate (AFC - MC)")+
  xlab("Chromosomal position (Mb)")+
  theme(text = element_text(family = "Lato", size = 14),
        strip.background = element_blank(),
        panel.spacing.x = unit(4, "mm"))+
  scale_fill_viridis_c(option = "viridis", direction = -1, limits = c(0.0,1))+
  scale_color_manual(values = c("black", "red"))+
  scale_x_continuous(breaks = scales::pretty_breaks(3))+
  expand_limits(x = 0)

Figure_S3_hist <- interval_diff %>%
  filter(!is.nan(p.value)) %>%
  ggplot(aes(x = co_diff/Mb))+
  geom_histogram()+
  theme_bw()+
  ylab("Count")+
  xlab("Difference in CO Rate (AFC - MC)")

Figure_S4 <- co_line  %>%
  mutate(interval_id = paste0(CHROM_schf, "_", POS_schf1)) %>%
  left_join(interval_diff_line) %>%
  #interval_diff_line %>%
  filter(!is.nan(p.value)) %>%
  mutate(p_adjust = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = p_adjust < 0.05) %>%
  distinct %>%
  ggplot(aes(x = POS_schf1/1000000, y = var_line, fill = p_adjust, color = significant)) +
  geom_point(size = 2, pch = 21)+
  geom_hline(yintercept = 0, linetype = 3)+
  facet_wrap(~CHROM_schf, scales = "free_x")+
  theme_bw()+
  ylab("Between Line Variance")+
  xlab("Chromosomal position (Mb)")+
  theme(text = element_text(family = "Lato", size = 14),
        strip.background = element_blank(),
        panel.spacing.x = unit(4, "mm"))+
  scale_fill_viridis_c(option = "viridis", direction = -1, limits = c(0.0,1))+
  scale_color_manual(values = c("black", "red"))+
  scale_x_continuous(breaks = scales::pretty_breaks(3))+
  expand_limits(x = 0)

ggsave(Figure_S3, filename = "figures/raw/FigureS3_raw.pdf", width = 8, height = 5, device = cairo_pdf)
ggsave(Figure_S4, filename = "figures/raw/FigureS4_raw.pdf", width = 8, height = 5, device = cairo_pdf)
  

# "#FF0000" "#00A08A" "#F2AD00" "#F98400" "#5BBCD6" "#046C9A"
plot_pal <- c(wes_palette("Darjeeling1", 5) %>% as.character, (wes_palette("Darjeeling2", 2) %>% as.character)[2])


pair_wise_df <- co_line %>%
  mutate(CHROM_schf_join = ifelse(CHROM_schf %in% c("XR", "XL"), "X", CHROM_schf)) %>%
  group_by(pop, line, CHROM_schf_join) %>%
  summarise(xo_freq = sum(co_count)/sum(co_n)) %>% 
  mutate(CHROM_schf_join = paste0("chr", CHROM_schf_join)) %>%
  ungroup(CHROM_schf_join)


comparisons <- data.frame(chr1 = c("chr2", "chr2", "chr2", "chr3", "chr3", "chr4"), 
                          chr2 = c("chr3", "chr4", "chrX", "chr4", "chrX", "chrX"))

# pair_wise_df <- co_line %>%
#   group_by(pop, line, CHROM_schf) %>%
#   summarise(recomb_mean = sum(co_count)/sum(Mb)) %>% 
#   mutate(CHROM_schf = paste0("chr", CHROM_schf)) %>%
#   ungroup(CHROM_schf)
# 
# comparisons <- data.frame(chr1 = c("chr2", "chr2", "chr2", "chr2", "chr3", "chr3","chr3", "chr4", "chr4"), 
#                           chr2 = c("chr3", "chr4", "chrXL", "chrXR", "chr4", "chrXL", "chrXR", "chrXL", "chrXR"))

recomb_pairs_df <- list()
for (i in 1:nrow(comparisons)){
  
  chrs_target <- comparisons[i,] %>% unlist %>% as.character
  
  pair_chr_df <- pair_wise_df %>%
    filter(CHROM_schf_join %in% chrs_target) %>%
    spread(key = CHROM_schf_join, value = xo_freq)
  
  names(pair_chr_df)[3:4] <- c("recomb_chr1", "recomb_chr2")
  recomb_pairs_df[[i]] <- data.frame(pair_chr_df, comparison = rep(paste(chrs_target, collapse = "_"), nrow(pair_chr_df)))
  
}

recomb_pairs_df <- bind_rows(recomb_pairs_df)

plot_pal <- c(wes_palette("Darjeeling1", 5) %>% as.character, (wes_palette("Darjeeling2", 2) %>% as.character)[2])

scaleFUN <- function(x) sprintf("%.4f", x)

figure2b<- recomb_pairs_df %>%
  mutate(comparison = gsub("_", " vs. ", comparison)) %>%
  mutate(comparison = gsub("chr", "", comparison)) %>%
  ggplot(aes(x = recomb_chr1, y = recomb_chr2, color = comparison))+
  geom_point(size = 2, alpha = 0.4, pch = 16)+
  geom_smooth(method = "lm", se = TRUE, size = 1.5)+
  theme_bw()+
  xlab("Crossover Frequency")+
  ylab("Crossover Frequency")+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "none")+
  scale_color_manual(values = plot_pal)+
  facet_wrap(~comparison, scales = "free")+ 
  scale_x_continuous(labels=scaleFUN)

figure2 <- figure2a / figure2b  + plot_layout(ncol = 2, widths = c(1, 2))

ggsave(filename = "figures/raw/Figure2_raw.pdf", figure2, width = 8, height = 5, useDingbats = FALSE)

# stats on average correlation
recomb_cor_stats <- recomb_pairs_df %>%
  group_by(comparison) %>%
  do(cor_recomb = cor.test(.$recomb_chr1, .$recomb_chr2) %>% tidy) %>%
  unnest 

recomb_cor_stats$estimate^2 %>% mean

####################################################
# not used
####################################################

# flag windows with too much missing data
# and/or unusually high recomb rates (for inspecting re: assembly errors)

geno_df <- co_counts_fine %>%
  mutate(pop = gsub("_.*", "", Indiv )) %>%
  mutate(line = gsub("_[0-9]{2}_[A-Z,0-9]{3}$", "", Indiv )) %>%
  ungroup %>%
  group_by(CHROM_schf, POS_schf1,POS_schf2) %>%
  summarise(n_genotyped = length(crossover_int)) %>%
  ungroup()

################## HEAT MAP STYLE
recomb_heat <- co_counts_fine %>%
  mutate(pop = gsub("_.*", "", Indiv )) %>%
  mutate(line = gsub("_[0-9]{2}_[A-Z,0-9]{3}$", "", Indiv )) %>%
  ungroup %>%
  left_join(geno_df) %>%
  #filter(n_genotyped > 4000) %>%
  group_by(pop, line, CHROM_schf, POS_schf1,POS_schf2) %>%
  summarise(co_count = sum(crossover_int, na.rm = TRUE), co_n = length(crossover_int), n_geno = mean(n_genotyped, na.rm = TRUE)) %>%
  ungroup %>%
  mutate(cM = (co_count / co_n) * 100) %>%
  mutate(Mb = abs(POS_schf1 - POS_schf2) / 1000000) %>%
  mutate(cM_Mb = cM / Mb) %>%
  group_by(CHROM_schf) %>%
  mutate(window_rank = dense_rank(POS_schf1)) %>%
  ungroup %>%
  mutate(high_recomb_flag = cM_Mb > 20) 

# flag suspect windows for closer examination
suspect_windows <- recomb_heat %>%
  group_by(CHROM_schf, POS_schf1, POS_schf2) %>%
  summarise(prop_high_recomb = sum(high_recomb_flag, na.rm = TRUE)/length(high_recomb_flag),
            n_geno = mean(n_geno))

saveRDS(suspect_windows, "meta/suspect_windows.rds")  
  
recomb_heat %>%  
  ggplot(aes(x = window_rank, fill = cM_Mb, y = line, color = high_recomb_flag))+
  geom_tile() +
  facet_wrap(.~CHROM_schf, scales = "free")+
  scale_fill_viridis_c()+
  scale_color_manual(values = c(NA, "red"))


my_bin <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(size = 2)+
    #geom_point(shape = 21, fill = "white", size = 6)+
    #geom_text(size = 3)+
    geom_smooth(method = "lm", se = TRUE, fill = "grey75", alpha = 0.4, color = "black")+
    scale_color_manual(values = plot_pal)
}

figure_2b <- co_line %>%
  group_by(pop, line, CHROM_schf) %>%
  summarise(recomb_mean = sum(co_count)/sum(Mb)) %>% 
  group_by(CHROM_schf) %>%
  #mutate(recomb_rank = dense_rank(recomb_mean)) %>% 
  #select(-recomb_mean) %>%
  ungroup() %>%
  #summarise(recomb_median = median(cM_Mb, na.rm = TRUE)) %>% 
  #summarise(recomb_rank = dense_rank(cM_Mb)) %>% 
  mutate(CHROM_schf = paste0("chr", CHROM_schf)) %>%
  spread(key = CHROM_schf, value = recomb_mean) %>%
  ungroup %>%
  mutate(line = gsub("[A-Z_]+", "", line)) %>%
  ggpairs(columns = c("chr2", "chr3", "chr4", "chrXL", "chrXR"),
          color = "pop",
          mapping=ggplot2::aes(color = pop, label = line),
          diag = list(continuous = "blank"),
          lower = list(continuous = my_bin ),
          upper = list(continuous = "blank"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 16))
