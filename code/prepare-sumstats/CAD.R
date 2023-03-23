library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
  transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2,
            af_UKBB = freq)


#### Coronary artery disease (CAD) ####

txt <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003116/cad.add.160614.website.txt",
  dir = "tmp-data", fname = "sumstats_CAD.txt")
sumstats <- fread2(txt,
                   select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "se_dgc", "p_dgc",
                              "effect_allele_freq", "median_info", "n_studies", "het_pvalue"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se",
                                 "p", "freq", "info", "n_studies", "het_pvalue"))

study_size <- c(
  # from Supp Table 1
  "5719 6545","206 259","278 312","505 1021","392 410","1010 3998",
  "1628 368","2083 2048","1216 653","658 5841","1802 466","2099 2690",
  "634 1608","1207 1288","1061 1467","1089 1147","877 2187",
  "2700 2758","361 2778","487 1381","758 3337","2791 3757","2095 503",
  "933 468","2905 2998","947 1008","1294 1529","843 318","933 468",
  "119 830","631 334","836 761","426 594","814 5999","322 857",
  "1926 2938","4651 4452","4380 3929","1535 772","1007 22286","402 448",
  "745 1389","397 2474","506 5335","259 4202","334 3446","2034 3210","454 8443")
study_size2 <- fread2(text = study_size, col.names = c("Nca", "Nco"))
nrow(study_size2)   # 48
colSums(study_size2)   # Nca: 61289 -- Nco: 126310  (slightly different from 60801 -- 123504)
(sumstats$n_eff <- sum(4 / rowSums(1 / study_size2)))  # 129014.3
4 / (1 / 60801 + 1 / 123504)  # 162972.6

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 9,455,778 variants to be matched.
# 80 ambiguous SNPs have been removed.
# 1,325,052 variants have been matched; 0 were flipped and 959,263 were reversed.

# comparing allele frequencies
freq2 <- ifelse(info_snp$beta * sumstats$beta[info_snp$`_NUM_ID_.ss`] < 0,
                1 - info_snp$freq, info_snp$freq)
hist(diff <- freq2 - info_snp$af_UKBB, "FD", xlim = c(-0.1, 0.1))

# better to use af of GWAS and INFO scores as well (then can use 0.7 instead of 0.5 in L35)
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))

is_bad <-
  abs(diff) > 0.1 | sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05

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

mean(is_bad)  # 0.01604843

bigassertr::assert_dir("data/sumstats")
saveRDS(filter(info_snp, !is_bad), "data/sumstats/CAD.rds")
