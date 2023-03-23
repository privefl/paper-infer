# setwd("~/florian/ldpred2-inference/")

library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
  transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2,
            af_UKBB = freq)


#### Type 1 diabetes (T1D) ####

# Download at https://datadryad.org/stash/dataset/doi:10.5061/dryad.ns8q3
# unzip("../sumstats/doi_10.5061_dryad.ns8q3__v1.zip", exdir = "../sumstats/T1D")

# Affymetrix
sumstats <- fread2(
  paste0("../sumstats/T1D/meta_chr_", 1:22),
  select = c("chromosome", "position", "a0", "a1",
             "info_score.A", "info_score.I", "beta.meta", "se.meta"),
  col.names = c("chr", "pos", "a0", "a1", "info1", "info2", "beta", "beta_se")) %>%
  filter(info1 > 0.4, info2 > 0.4, beta_se > 0) %>%
  select(-info1, -info2)


info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 8,531,891 variants to be matched.
# 79 ambiguous SNPs have been removed.
# Some duplicates were removed.
# 1,127,489 variants have been matched; 519 were flipped and 1,769 were reversed.

info_snp$n_eff <- 4 / (1 / 1930 + 1 / (1455 + 1490 + 1884)) + 4 / (1 / 3983 + 1 / 3999)


sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))

is_bad <-
  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05

ggplot(slice_sample(data.frame(sd_ldref, sd_ss, is_bad), n = 50e3)) +
  geom_point(aes(sd_ldref, sd_ss, color = is_bad), alpha = 0.5) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

mean(is_bad)  # 0.005969903

bigassertr::assert_dir("data/sumstats")
saveRDS(filter(info_snp, !is_bad), "data/sumstats/T1D.rds")
