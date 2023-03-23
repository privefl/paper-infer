library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
  transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2,
            af_UKBB = freq)


#### Breast cancer (BRCA) ####

tgz <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
  dir = "tmp-data", fname = "sumstats_BRCA.txt.gz")
R.utils::gunzip(tgz, overwrite = FALSE, remove = FALSE)

sumstats <- fread2("tmp-data/sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_onco_icogs_gwas_beta", "bcac_onco_icogs_gwas_se",
                              "bcac_onco_icogs_gwas_eaf_controls"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se", "freq"))
sumstats$n_eff <- 4 / (1 / 137045 + 1 / 119078)

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 11,792,542 variants to be matched.
# 110 ambiguous SNPs have been removed.
# Some duplicates were removed.
# 1,411,710 variants have been matched; 0 were flipped and 6 were reversed.

# comparing allele frequencies
freq2 <- ifelse(info_snp$beta * sumstats$beta[info_snp$`_NUM_ID_.ss`] < 0,
                1 - info_snp$freq, info_snp$freq)
hist(diff <- freq2 - info_snp$af_UKBB, "FD", xlim = c(-0.08, 0.08))

# better to use af of GWAS and INFO scores as well (then can use 0.7 instead of 0.5 in L35)
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))

is_bad <-
  abs(diff) > 0.05 | sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05

library(ggplot2)
ggplot(slice_sample(data.frame(sd_ldref, sd_ss, is_bad), n = 50e3)) +
  geom_point(aes(sd_ldref, sd_ss, color = is_bad), alpha = 0.5) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

mean(is_bad)  # 0.006124487

bigassertr::assert_dir("data/sumstats")
saveRDS(filter(info_snp, !is_bad), "data/sumstats/BrCa.rds")
