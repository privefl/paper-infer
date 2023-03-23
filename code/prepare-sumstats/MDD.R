library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
  transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2,
            af_UKBB = freq)


#### MDD ####

gz <- runonce::download_file(
  "https://figshare.com/ndownloader/files/34427408",
  dir = "../sumstats", fname = "daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz")
sumstats <- fread2(gz,
                   select = c("CHR", "BP", "A1", "A2", "OR", "SE",
                              "INFO", "FRQ_A_45396", "FRQ_U_97250", "Neff_half"),
                   col.names = c("chr", "pos", "a1", "a0", "or", "beta_se",
                                 "info", "freq1", "freq2", "Neff_half")) %>%
  as_tibble() %>%
  mutate(beta = log(or), or = NULL, chr = as.integer(chr),
         freq = (freq1 * 45396 + freq2 * 97250) / (45396 + 97250),
         freq1 = NULL, freq2 = NULL,
         n_eff = 2 * Neff_half, Neff_half = NULL,
         info = pmin(info, 1)) %>%
  filter(n_eff > (0.5 * max(n_eff)),
         info > 0.4) %>%
  print()


info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 9,343,125 variants to be matched.
# 76 ambiguous SNPs have been removed.
# 1,314,499 variants have been matched; 0 were flipped and 612,915 were reversed.


info_snp2 <- info_snp %>%
  mutate(sd_af = sqrt(2 * freq * (1 - freq)),
         sd_ss = 2 / sqrt(n_eff * beta_se^2 + beta^2),
         sd_ss2 = sd_ss / sqrt(info))

info_snp2$freq2 <- ifelse(info_snp2$beta * sumstats$beta[info_snp2$`_NUM_ID_.ss`] < 0,
                          1 - info_snp2$freq, info_snp2$freq)

qplot(af_UKBB, freq2, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr() +
  coord_equal() +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Allele frequencies from the validation set",
       y = "Allele frequencies from the summary statistics")
# Some are reversed?

hist(diff <- with(info_snp2, abs(af_UKBB - freq2)), "FD", xlim = c(0, 0.1))

info_snp2$is_bad <- with(info_snp2, diff > 0.05 |
                           sd_ss2 < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
                           sd_ss2 < 0.1 | sd_af < 0.05)

qplot(sd_af, sd_ss2, color = is_bad, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 50e3)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

mean(info_snp2$is_bad)  # 0.002820086

bigassertr::assert_dir("data/sumstats")
saveRDS(filter(info_snp2, !is_bad), "data/sumstats/MDD.rds")
