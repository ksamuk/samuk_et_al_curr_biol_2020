library("tidyverse")
library("parallel")
library("parallelsugar")
library("aod")
library("rvg")
library("officer")
library("SNPRelate")
library("lmerTest")
library("lme4")
library("boot")
library("wesanderson")
library("optimx")
source("functions/theme_presentation.R")

#options(device = "quartz")

pop_labels <- read.table("data/co_counts_manual.txt", h = T) %>%
  select(pop, line, plate, Indiv)

# the crossover counts from rqtl
co_counts_rqtl <- read.table("data/rqtl/co_estimates/all.txt", h = T) %>%
  rename(Indiv = id) %>%
  left_join(pop_labels)%>%
  filter(crossovers <= 4) %>%
  distinct %>%
  mutate(pop = as.factor(pop), line = as.factor(line)) %>%
  group_by(pop, line, Indiv) %>%
  summarise(crossovers = sum(crossovers, na.rm = TRUE)) %>%
  ungroup

###########################################
# QST
###########################################

# QST CALCULATION

# QST = sigma^2(pop) / [ sigma^2(pop))+ (sigma^2(family(pop)))]
# "VB (pop) is the among population component of variance for the trait, 
# and VAW (family(pop)) is the additive genetic variance within populations."
# THE FACTOR OF 2 IN FROM OF VW/VAW IS DROPPED FOR COMPLETELY INBRED LINES
# PAGE 453 WALSH AND LYNCH VOL 2, BONIN et al. 1996

compute_qst <- function(filtered_df, bootstrap = FALSE){
  
  filtered_df <- filtered_df %>%
    mutate(pop = as.factor(pop), line = as.factor(line))
  
  var_glm <- lmer(crossovers ~ (1|pop/line), data = filtered_df, 
                   control = lmerControl(
                     optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
  #var_glm <- lmer(crossovers ~ (1|pop/line), data = filtered_df,  control = lmerControl(optimizer ="Nelder_Mead"))
  #var_glm <- glmer(crossovers ~ (1|pop/line), family = "poisson", data = filtered_df)
  
  var_comp <- VarCorr(var_glm) %>% data.frame
  var_pop <- var_comp$vcov[2]
  var_cross <- var_comp$vcov[1]
  
  qst <- var_pop / (var_pop + var_cross)
  
  return(list(qst = qst, var_within = var_cross, var_between = var_pop, lm = var_glm))

} 

# compute observed qst and variance components
qst_vals <- compute_qst(co_counts_rqtl)
# variance components
# $qst
# [1] 0.2269951
# 
# $var_within
# [1] 0.08817176
# 
# $var_between
# [1] 0.02589189
summary(qst_vals$lm)

###########################################
# FST
###########################################

# sample meta data
meta_dat <- read.csv("meta/plate_guide.csv")

meta_dat <- meta_dat %>% 
  mutate(id = gsub("P", "L", plate) %>% paste0(., "_", flip)) %>%
  select(id, pop, line_id, type, unique_id)

# sample depth data
dp_dat <- read.table("meta/mean_depth.txt", header = TRUE)
names(dp_dat)[1] <- "id"

meta_dat <- left_join(dp_dat, meta_dat)

# sample homozygostity data (for confirming wild/inbred status)
hom_dat <- read.table("meta/h_obs_filter_scores.txt", header = TRUE)
meta_dat <- left_join(meta_dat, hom_dat)

# indentification of individuals with highest depth in F1 pairs + Mothers
# (we only want to keep one good individual from each F1 family, including mothers)

meta_dat <- meta_dat %>%
  mutate(type2 = ifelse(type %in% c("F", "F1"), "F-F1", "M-Lab")) %>%
  group_by(pop, type2, line_id) %>%
  arrange(pop, type2, line_id) %>%
  mutate(dp_rank = rank(-mean_DP)) %>% 
  mutate(dp_rank = ifelse(type2 == "F-F1", dp_rank, 1)) 

# pick out the wild individuals with rank == 1 (see above)
# also filter out low DP individuals (hamfisted, consider site-level filter)
wild_dat <- meta_dat %>%
  filter(mean_DP > 20) %>%
  #filter(p_homozygous < 0.92) %>%
  filter(dp_rank == 1) %>%
  filter(pop %in% c("AFC", "MC")) %>% 
  filter(type != "L") 

mc_ind <- wild_dat %>%
  filter(pop == "MC") %>%
  pull(id) %>%
  as.character

afc_inds <- wild_dat %>%
  filter(pop == "AFC") %>%
  pull(id) %>%
  #sample(size = length(mc_ind))%>%
  as.character

fst_inds <- c(mc_ind, afc_inds)
fst_pops <- c(rep("MC", length(mc_ind)), rep("AFC", length(afc_inds)))

############################################
# load raw data and prune snps for LD + HWE
############################################

genofile <- snpgdsOpen("data/snp_relate/dpse_ddrad.gds", allow.duplicate = TRUE)

# hwe check
# now handled upstream, but the code is shown below

#hwe_snps <- snpgdsHWE(genofile, with.id = TRUE, sample.id = fst_inds)

#snpset_vec <- data.frame(hwe_snps$snp.id, hwe_snps$pvalue) %>%
#  filter(hwe_snps.pvalue > 0.05) %>%
#  pull(hwe_snps.snp.id) %>%
#  as.character

# LD pruning per Lotterhos reccomendations
snpset <- snpgdsLDpruning(genofile, autosome.only = FALSE, maf = 0.1, 
                          missing.rate = 0.5, method = "r", slide.max.bp = 10000,
                          ld.threshold = 0.2, sample.id = fst_inds)

snpset_vec <- unlist(snpset)
names(snpset_vec) <- NULL

################
# FST Bootstrap
################


# compiute a single observed value of fst
fst_obs <- function(maf_cut = 0.1, missing_rate_cut = 0.2){
  
  snpset_vec_filt <- snpset_vec[!grepl("3_", snpset_vec)]
  
  snp_tab <- snpgdsGetGeno(genofile, snp.id = snpset_vec_filt, sample.id = fst_inds, with.id = TRUE, verbose = FALSE)
  snp_boot <- snp_tab$genotype
  
  #snpgdsCreateGeno("data/snp_relate/dpse_ddrad_boot1.gds", genmat = snp_boot, sample.id = snp_tab$sample.id, snp.id = NULL, snpfirstdim = FALSE)
  genoboot <- snpgdsOpen("data/snp_relate/dpse_ddrad.gds", allow.duplicate = TRUE)
  
  fst_dat <- snpgdsFst(genoboot, sample.id = fst_inds, method=c("W&C84"),
                       population = factor(fst_pops), maf = maf_cut, 
                       autosome.only = FALSE, missing.rate = missing_rate_cut, verbose = FALSE)
  
  snpgdsClose(genoboot)
  #file.remove("C:/Users/Kieran/Dropbox (Duke Bio_Ea)/Projects/gt_seq_rev/data/snp_relate/dpse_ddrad_boot1.gds")
  #file.remove("data/snp_relate/dpse_ddrad_boot1.gds")
  
  return(fst_dat)
  
}

fst_dat <- fst_obs()

# generate a non-parametric bootstrap sample of FST by resampling loci
# part of the QST-FST bootstrap pipeline 
# (note that null QST is estimated via *parametric* bootstap later)

fst_boot <- function(maf_cut = 0.1, missing_rate_cut = 0.2){
  
  #genoboot_file <- "C:/Users/Kieran/Dropbox (Duke Bio_Ea)/Projects/gt_seq_rev/data/snp_relate/dpse_ddrad_boot1.gds"
  genoboot_file <- "data/snp_relate/dpse_ddrad_boot1.gds"
  
  snpset_vec_filt <- snpset_vec[!grepl("3_", snpset_vec)]
  
  snp_tab <- snpgdsGetGeno(genofile, snp.id = snpset_vec_filt, sample.id = fst_inds, with.id = TRUE, verbose = FALSE)
  snp_boot <- snp_tab$genotype[sample(1:nrow(snp_tab$genotype), replace = TRUE),]
  
  snpgdsCreateGeno(genoboot_file, genmat = snp_boot, sample.id = snp_tab$sample.id, snp.id = NULL, snpfirstdim = FALSE)
  genoboot <- snpgdsOpen(genoboot_file, allow.duplicate = TRUE)
  
  fst_dat <- snpgdsFst(genoboot, sample.id = fst_inds, method=c("W&C84"),
                       population = factor(fst_pops), maf = maf_cut, 
                       autosome.only = FALSE, missing.rate = missing_rate_cut, verbose = FALSE)
  
  snpgdsClose(genoboot)
  file.remove(genoboot_file)
  #file.remove("data/snp_relate/dpse_ddrad_boot1.gds")
  
  fst_dat$MeanFst
  
}


######################################################################
# estimate the distributon of QST - FST via non parametric bootstrap
######################################################################

# QST parametric bootstrap (after whitlock + guillame, gilbert + whitlock, adapted for inbred lines)

qst_parboot <- function(var_within, n_fam, n_pop) {
  
  cat(".")
  
  # bootstrap estimate of variation within
  # modified from whitlock and guillame
  # dfs are all 1 here, so 
  var_within_boot <- rchisq(1, n_fam) * var_within^2
  
  # fst boostrap
  fst_boot_mean <- fst_boot()
  
  # bootstrap estimate of variation between
  var_between_boot =  ((fst_boot_mean * var_within_boot) / ((1 - fst_boot_mean))) * rchisq(1, n_pop -1)
  
  # return bootstrap estimate of qst
  qst_boot <- var_between_boot / (var_between_boot + var_within_boot)
  
  data.frame(fst_boot_mean, var_within_boot, var_between_boot, qst_boot)
  
}

qst_vals_rqtl <- compute_qst(co_counts_rqtl)

# perform boostrap replicates for both variance estimates
qst_boot_rqtl <- replicate(10000, qst_parboot(qst_vals_rqtl$var_within, n_pop = 2, n_fam = 17), simplify = FALSE) 
qst_boot_rqtl <- qst_boot_rqtl %>% bind_rows

saveRDS(qst_boot_rqtl , "data/qst_boot/qst_boot_rqtl.rds")
qst_boot_rqtl <- readRDS("data/qst_boot/qst_boot_rqtl.rds") %>% bind_rows()

qst_boot_rqtl  <- qst_boot_rqtl  %>%
  bind_rows %>%
  mutate(var_between_boot2 = ((2*fst_boot_mean * var_within_boot) / ((1 - fst_boot_mean))) * rchisq(1, 1)) %>%
  mutate(qst_boot2 = var_between_boot2 / (var_between_boot2 + var_within_boot)) %>%
  mutate(qst_fst_boot2 = qst_boot2 - fst_boot_mean)

qst_vals

var_within <- 0.08817176 # REML assuming normality
var_between <- 0.02589189  # REML 
qst_obs <- var_between / (var_within + var_between)

# fst
mean_fst <- 0.006928394

qst_fst_obs <- qst_obs - mean_fst

# qst_n pvalue
# 9.999e-05
qst_boot_rqtl$qst_fst_boot <- qst_boot_rqtl$qst_boot - qst_boot_rqtl$fst_boot_mean
(sum(qst_boot_rqtl$qst_fst_boot > qst_fst_obs ) + 1)/ (length(qst_boot_rqtl$qst_fst_boot) + 1)

# plots

# wes palatte ftw
plot_pal<- wes_palette("Darjeeling1", n = 3, type = "discrete")%>% as.character

# build annotation dataframe
# for displating annotations on top of main plot
# 95% CI for obs QST-FST
qst_fst_obs_boot <- qst_boot_rqtl %>%
  mutate(qst_fst_obs_boot = qst_fst_obs - qst_fst_boot) %>%
  pull(qst_fst_obs_boot)

qst_fst_obs_mean <- median(qst_fst_obs_boot)
qst_fst_obs_95 <- quantile(qst_fst_obs_boot, probs = 0.95) %>% as.numeric
qst_fst_obs_05 <- quantile(qst_fst_obs_boot, probs = 0.05) %>% as.numeric

qst_fst_exp_boot <- qst_boot_rqtl %>%
  pull(qst_fst_boot)

qst_fst_exp_mean <- median(qst_fst_exp_boot)
qst_fst_exp_95 <- quantile(qst_fst_exp_boot, probs = 0.95) %>% as.numeric
qst_fst_exp_05 <- quantile(qst_fst_exp_boot, probs = 0.05) %>% as.numeric


is_fst <- !(c(F,F,T))
qst_fst_mean <- c(qst_fst_exp_mean, qst_fst_obs_mean, NA)
qst_fst_95 <- c(qst_fst_exp_95, qst_fst_obs_95, qst_fst_obs_mean)
qst_fst_05 <- c(qst_fst_exp_05, qst_fst_obs_05, qst_fst_obs_mean)
qst_y <- c(100, 100, 0)
qst_yend <- c(100, 100, 100)
qst_text_height <- c(300, 300, 120)
qst_text_pos <- c(qst_fst_exp_mean, qst_fst_obs_mean, qst_fst_obs_mean)
qst_text <- c("QST-FST\n Expected", "QST-FST\n Observed", "QST Observed")
statistic <- c("qst_fst_boot", "qst_fst_obs_boot", "fst_obs")

annot_df <- data.frame(statistic, is_fst, qst_fst_mean, qst_text_height,
                       qst_text_pos, qst_fst_95, qst_fst_05, qst_text, qst_y, qst_yend)

# observed values of fst
fst_obs <- data.frame(statistic = "fst_obs", value= fst_dat$FstSNP)

figure3<- qst_boot_rqtl %>%
  mutate(qst_fst_obs_boot = qst_fst_obs - qst_fst_boot) %>%
  select(qst_fst_obs_boot, qst_fst_boot) %>%
  gather(key = statistic, value = value) %>% 
  bind_rows(fst_obs) %>%
  mutate(is_fst = !(statistic == "fst_obs")) %>%
  ggplot(aes(x = value, fill = statistic)) +
  geom_histogram(bins = 400, position="identity", alpha = 0.95)+
  geom_point(data = annot_df, aes(x = qst_fst_mean, y = qst_y ), fill = "black", size = 3)+
  geom_segment(data = annot_df, aes(x = qst_fst_05, xend = qst_fst_95, 
                                    y = qst_y, yend = qst_yend), color = "black", fill = "black", size = 1)+
  geom_text(data = annot_df, aes(x = qst_text_pos , y = qst_text_height, label = qst_text), size = 5)+
  #facet_grid(is_fst~.,scales = "free_y") +
  facet_wrap(~is_fst,nrow = 2,scales = "free")+
  ylab("Frequency")+
  xlab("QST-FST or FST")+
  xlim(-0.05, 0.25)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = plot_pal, labels = c("FST", "QST-FST (expected)", "QST-FST (observed)"))+
  theme_bw()+
  theme(strip.text.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=16))

ggsave(filename = "figures/Figure3.pdf", figure3, width = 8, height = 8, useDingbats = FALSE)


