library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
  transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2,
            af_UKBB = freq)


#### Vitamin D ####

# See google drive link in https://doi.org/10.1038/s41467-017-02662-2
# Should have used Allele1 -> a1 and Allele2 -> a0
sumstats <- fread2("../sumstats/upload_25HydroxyVitaminD_QC.METAANALYSIS1.txt",
                   select = c("Chr", "Pos", "MarkerName", "Allele1", "Allele2",
                              "Effect", "StdErr", "SAMPLESIZE", "HetPVal"),
                   col.names = c("chr", "pos", "rsid", "a1", "a0",
                                 "beta", "beta_se", "n_eff", "het_pval")) %>%
  mutate_at(c("a0", "a1"), toupper)

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 2,543,566 variants to be matched.
# 35 ambiguous SNPs have been removed.
# 1,028,171 variants have been matched; 6,370 were flipped and 498,902 were reversed.

hist(-log10(info_snp$het_pval), "FD")

info_snp2 <- info_snp %>%
  mutate(sd_ldref = sqrt(2 * af_UKBB * (1 - af_UKBB)),
         sd_ss = 1 / sqrt(n_eff * beta_se^2 + beta^2),
         sd_ss = sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5))

info_snp2$is_bad <- with(info_snp2, het_pval < 1e-5 |
                           sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) |
                           sd_ss < 0.05 | sd_ldref < 0.05)

ggplot(slice_sample(info_snp2, n = 50e3)) +
  geom_point(aes(sd_ldref, sd_ss, color = is_bad), alpha = 0.5) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

mean(info_snp2$is_bad)  # 0.00399739

bigassertr::assert_dir("data/sumstats")
saveRDS(filter(info_snp2, !is_bad), "data/sumstats/VitaminD.rds")
