library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
  transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2,
            af_UKBB = freq)


#### Type 2 diabetes (T2D) ####

# Download from http://diagram-consortium.org/downloads.html
# DIAGRAM 1000G GWAS meta-analysis Stage 1 Summary statistics
# Published in Scott et al (2017)
# unzip("../sumstats//METAANALYSIS_DIAGRAM_SE1.zip", exdir = "../sumstats")
sumstats <- fread2("../sumstats//METAANALYSIS_DIAGRAM_SE1.txt")
sumstats <- tidyr::separate(sumstats, "Chr:Position", c("chr", "pos"), convert = TRUE)
names(sumstats) <- c("chr", "pos", "a1", "a0", "beta", "beta_se", "p", "N")
sumstats$n_eff <- 72143

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 12,056,346 variants to be matched.
# 89 ambiguous SNPs have been removed.
# 1,408,283 variants have been matched; 154 were flipped and 650,482 were reversed.

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

mean(is_bad)  # 0.02798159

bigassertr::assert_dir("data/sumstats")
saveRDS(filter(info_snp, !is_bad), "data/sumstats/T2D.rds")
