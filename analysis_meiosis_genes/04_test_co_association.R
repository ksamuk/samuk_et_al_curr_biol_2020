# test association between CO rate and divergence in recomb genes

library("tidyverse")
library("broom")
library("rvg")
library("officer")
library("wesanderson")
library("patchwork")

# read in residue genotypes
residue_df <- readRDS("data/fasta_recomb/summaries/residue_df.rds")

# determine residues with at least one non-synonymous change
ns_residues <- residue_df %>%
  group_by(gene_name_full, residue_num) %>%
  summarise(has_ns_genos = any(hap1_aa_diff == "non-synonymous" | hap2_aa_diff == "non-synonymous"),
            count_ns_genos = sum(hap1_aa_diff == "non-synonymous" | hap2_aa_diff == "non-synonymous")) %>%
  filter(has_ns_genos)

# total n genes:
residue_df$gene_name_full %>% unique %>% length

# filter for sites with at least one  change
ns_df <- residue_df %>%
  left_join(ns_residues) %>%
  filter(has_ns_genos) %>%
  filter(count_ns_genos > 4) %>%
  select(-has_ns_genos, -count_ns_genos)

# toss singletons

# collapse the aa to bialleic genotypes
# (ignoring non-synonymous changes)

ns_df <- ns_df %>%
  mutate(aa_geno = (hap1_aa_diff == "non-synonymous") + (hap2_aa_diff == "non-synonymous")) %>%
  select(gene_name_full, sample_id, residue_num, aa_geno)

# join in co data
co_df <- read.table("data/fasta_recomb/summaries/co_rate_line.txt") %>%
  rename(sample_id = line) %>%
  mutate(sample_id = as.character(sample_id)) %>%
  select(-pop)

# harmonize names
ns_df <- ns_df %>%
  filter(grepl("M[0-9]*|A[0-9]*", sample_id)) %>%
  filter(!grepl("MV.*", sample_id)) %>%
  filter(!grepl("A6|A12", sample_id)) %>%
  mutate(sample_id = gsub("A", "AFC_", sample_id))%>%
  mutate(sample_id = gsub("M", "MC_", sample_id))

ns_co <- ns_df %>%
  left_join(co_df) 

# confirm filters
ns_freqs <- ns_co %>%
  group_by(gene_name_full, residue_num) %>%
  summarise(ns_freq = sum(aa_geno)/(length(aa_geno)*2))

# fit a linear model to each gene
co_fits <- ns_co %>%
  left_join(ns_freqs) %>%
  filter(ns_freq > 0.1) %>% # maf filter
  group_by(gene_name_full, residue_num) %>%
  do(co_geno_fit = lm(.$mean_co ~ .$aa_geno)) %>%
  tidy(co_geno_fit) 

co_fits$p.adjust <- p.adjust(co_fits$p.value, method = "BY")

# 357 residues
co_fits %>%
  ungroup %>%
  filter(term != "(Intercept)") %>% View
  filter(p.adjust < 0.05)

# 33 genes
ns_genes <- co_fits$gene_name_full %>% unique
length(ns_genes)

residue_df %>%
  filter(gene_name_full %in% ns_genes) %>%
  select(gene_name_full, residue_num) %>%
  distinct %>% 
  nrow


# Figure 4a

plot_pal <- wes_palette("Darjeeling1", n = 2, type = "discrete")%>% as.character

co_plot <- co_fits %>%
  ungroup %>%
  filter(term != "(Intercept)") %>%
  mutate(effect_size = abs(estimate)) %>%
  mutate(effect_size_cohen = abs(estimate) / (std.error*sqrt(17))) %>%
  mutate(gene_name_full = gsub("_.*", "", gene_name_full)) %>%
  mutate(signif = p.adjust < 0.05)

figure_4a <- co_plot %>%
  ggplot(aes(x = gene_name_full, y = estimate, color = signif, fill = std.error))+
  geom_hline(yintercept = 0, lty = 2) +
  geom_jitter(width = 0.0, pch = 21, size = 2.5, alpha = 0.9)+
  geom_jitter(data = co_plot %>% filter(signif), width = 0.0, pch = 21, size = 2.5, alpha = 1.0, fill = "red", color = "black")+
  #facet_grid(.~gene_name_full, scale = "free_x")+
  xlab("Gene Name")+
  ylab("Effect Size\n(Regression Coefficient)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(family = "Lato",size=14))+
  scale_color_manual(values = c("black","red"))+
  scale_fill_viridis_c()

co_fits %>%
  ungroup %>%
  filter(term != "(Intercept)") %>%
  mutate(effect_size = abs(estimate)) %>%
  mutate(effect_size_cohen = abs(estimate) / (std.error*sqrt(17))) %>%
  mutate(gene_name_full = gsub("_.*", "", gene_name_full)) %>%
  mutate(signif = p.adjust < 0.1) %>%
  filter(signif)

figure_4b <- ns_co %>%
  mutate(gene_residue = paste0(gene_name_full, "_", residue_num)) %>%
  filter(gene_name_full == "asp_FBpp0282877" & residue_num == 591 | 
         gene_name_full == "mei-41_FBpp0282161" & residue_num == 1791 | 
         gene_name_full == "mei-41_FBpp0282161" & residue_num == 1642) %>% 
  filter(aa_geno != 1) %>%
  mutate(gene_residue = gsub("_FB.*_", ":", gene_residue)) %>%
  mutate(pop = gsub("_.*", "", sample_id)) %>%
  mutate(aa_geno = ifelse(aa_geno == 0, "0/0", "1/1")) %>%
  ggplot(aes(x = as.factor(aa_geno), y = mean_co, color = as.factor(aa_geno))) +
  #geom_boxplot()+
  stat_summary(fun.y = mean, geom = "point", size = 2) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, size = 0.75)+
  geom_jitter(width = 0.1, alpha = 0.7, shape = 16, size = 2)+
  facet_wrap(~gene_residue)+
  xlab("Genotype")+
  ylab("Mean COs per\nChromosome Arm")+
  theme_bw()+
  theme(text = element_text(family = "Lato", size=14),
        strip.background = element_blank(),
        legend.position = "none")+
  scale_color_manual(values = rev(plot_pal))+
  scale_alpha_continuous(trans = "reciprocal")

figure4 <- figure_4a / figure_4b  


ggsave(filename = "figures/raw/Figure5_raw.pdf", figure4, width = 10, height = 6,  device = cairo_pdf)

################### total substitution count analysis
# decided this doesn't really make sense, so sticking with the NS sites analysis above

ns_genes <- ns_df %>%
  pull(gene_name_full) %>%
  unique

ns_per_gene <- residue_df %>%
  left_join(ns_residues) %>%
  filter(sample_id != "CDS") %>%
  filter(gene_name_full %in% ns_genes) %>%
  select(gene_name_full, sample_id, residue_num, hap1_aa, hap2_aa, cds_aa, hap1_aa_diff, hap2_aa_diff) %>%
  mutate(count_ns = (hap1_aa_diff == "non-synonymous") + (hap2_aa_diff == "non-synonymous")) %>%
  group_by(gene_name_full, sample_id) %>%
  summarise(count_ns = sum(count_ns, na.rm = TRUE), n_sites = max(residue_num)) %>%
  mutate(prop_ns = count_ns / n_sites) %>%
  mutate(prop_ns_zscore = scale(prop_ns))

ns_per_sample <- residue_df %>%
  left_join(ns_residues) %>%
  filter(sample_id != "CDS") %>%
  select(gene_name_full, sample_id, residue_num, hap1_aa, hap2_aa, cds_aa, hap1_aa_diff, hap2_aa_diff) %>%
  mutate(count_ns = (hap1_aa_diff == "non-synonymous") + (hap2_aa_diff == "non-synonymous")) %>%
  group_by(gene_name_full, sample_id) %>%
  summarise(count_ns = sum(count_ns, na.rm = TRUE), n_sites = max(residue_num)) %>%
  ungroup %>%
  group_by(sample_id) %>%
  summarise(total_ns = sum(count_ns, na.rm = TRUE), total_sites = sum(n_sites, na.rm = TRUE)) %>%
  mutate(prop_ns = total_ns / total_sites)


# checking that summaries appear sane
ns_per_sample %>%
  ggplot(aes(x = prop_ns))+
  geom_histogram() 

ns_per_gene %>%
  ggplot(aes(x = prop_ns))+
  geom_histogram() +
  facet_wrap(~gene_name_full)

# join in co data
co_df <- read.table("data/fasta_recomb/summaries/co_rate_line.txt") %>%
  rename(sample_id = line) %>%
  mutate(sample_id = as.character(sample_id)) %>%
  select(-pop)

# harmonize names
ns_per_gene <- ns_per_gene %>%
  filter(grepl("M[0-9]*|A[0-9]*", sample_id)) %>%
  filter(!grepl("MV.*", sample_id)) %>%
  filter(!grepl("A6|A12", sample_id)) %>%
  mutate(sample_id = gsub("A", "AFC_", sample_id))%>%
  mutate(sample_id = gsub("M", "MC_", sample_id))

ns_per_sample <- ns_per_sample %>%
  filter(grepl("M[0-9]*|A[0-9]*", sample_id)) %>%
  filter(!grepl("MV.*", sample_id)) %>%
  filter(!grepl("A6|A12", sample_id)) %>%
  mutate(sample_id = gsub("A", "AFC_", sample_id))%>%
  mutate(sample_id = gsub("M", "MC_", sample_id))

ns_per_gene <- ns_per_gene %>%
  left_join(co_df) 

ns_per_sample <- ns_per_sample %>%
  left_join(co_df)

ns_per_gene %>%
  ggplot(aes(x = prop_ns_zscore, y = mean_co))+
  geom_point()+
  facet_wrap(~gene_name_full)

ns_per_sample %>%
  ggplot(aes(x = prop_ns, y = mean_co))+
  geom_point()
