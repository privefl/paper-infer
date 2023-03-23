library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
  transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2,
            af_UKBB = freq)


#### Prostate cancer (PRCA) ####

zip <- runonce::download_file(
  "http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip",
  dir = "tmp-data", fname = "sumstats_PRCA.zip")
unzip(zip, exdir = "tmp-data", overwrite = FALSE)

sumstats <- fread2(
  "tmp-data/meta_v3_onco_euro_overall_ChrAll_1_release.txt",
  select = c("Chr", "position", "Allele1", "Allele2", "Effect", "StdErr",
             "Pvalue", "Freq1", "OncoArray_imputation_r2"),
  col.names = c("chr", "pos", "a1", "a0", "beta", "beta_se", "p", "freq", "info")
) %>%
  mutate(a0 = toupper(a0), a1 = toupper(a1)) %>%
  filter(info > 0.4, beta_se > 0)

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 11,792,542 variants to be matched.
# 110 ambiguous SNPs have been removed.
# Some duplicates were removed.
# 1,411,710 variants have been matched; 0 were flipped and 6 were reversed.

sizes <- bigreadr::fread2(text = c("44825\t27904","20219\t20440","1854\t1894","3650\t3940",
                                   "474\t482","1458\t512","2068\t2993", "4600\t2941"))
info_snp$n_eff <- print(sum(4 / rowSums(1 / sizes)))  # 135316.1
4 / (1 / 79148 + 1 / 61106)  # 137933.1

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

mean(is_bad)  # 0.01076782

bigassertr::assert_dir("data/sumstats")
saveRDS(filter(info_snp, !is_bad), "data/sumstats/PrCa.rds")
