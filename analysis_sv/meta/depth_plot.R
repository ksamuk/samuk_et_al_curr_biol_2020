# sv histogram

library("tidyverse")

cov_df <- read.table("analysis_sv/data/14.coverage.txt")
names(cov_df) <- c("chrom", "depth", "count", "total_len", "density")

cov_df %>%
  filter(depth < 100) %>%
  filter(!grepl("Unk", chrom)) %>%
  ggplot(aes(x = depth, y = count))+
  geom_histogram(stat = "identity", pad = 0)+
  facet_wrap(~chrom, scales = "free_y")
