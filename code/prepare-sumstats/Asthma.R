# setwd("~/florian/ldpred2-inference/")

library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
  transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2,
            af_UKBB = freq)


#### Asthma ####

zip <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006862/TAGC_meta-analyses_results_for_asthma_risk.zip",
  dir = "../sumstats")
unzip(zip, exdir = "../sumstats", overwrite = FALSE)


sumstats <- fread2(
  "~/florian/sumstats/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv",
  select = c("chr", "position", "alternate_allele", "reference_allele",
             "European_ancestry_beta_fix", "European_ancestry_se_fix"),
  col.names = c("chr", "pos", "a1", "a0", "beta", "beta_se"))


info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 2,001,280 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 946,092 variants have been matched; 0 were flipped and 0 were reversed.

info_snp$n_eff <- 4 / (1 / 19954 + 1 / 107715)

sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))

is_bad <-
  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05

ggplot(slice_sample(data.frame(sd_ldref, sd_ss, is_bad), n = 50e3)) +
  geom_point(aes(sd_ldref, sd_ss, color = is_bad), alpha = 0.5) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

mean(is_bad)  # 0.000988276

bigassertr::assert_dir("data/sumstats")
saveRDS(filter(info_snp, !is_bad), "data/sumstats/Asthma.rds")
