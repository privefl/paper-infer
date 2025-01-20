library(bigsnpr)
library(bigreadr)

NCORES <- 15  # TO MODIFY

## Information for the variants provided in the LD reference
map_ldref <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = "tmp-data", fname = "map_hm3_plus.rds"))

## The one before is for HapMap3+; if using HapMap3 instead, get
readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/36360900",
  dir = "tmp-data", fname = "map_hm3.rds"))


## Breast cancer summary statistics
# options(timeout = 1000)
tgz <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
  dir = "tmp-data", fname = "sumstats_BRCA.txt.gz")
R.utils::gunzip(tgz, overwrite = FALSE, remove = FALSE)

sumstats <- fread2("tmp-data/sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_onco_icogs_gwas_beta",
                              "bcac_onco_icogs_gwas_se"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se"))
sumstats$n_eff <- 4 / (1 / 137045 + 1 / 119078)

info_snp <- snp_match(sumstats, map_ldref)
# 11,792,542 variants to be matched.
# 110 ambiguous SNPs have been removed.
# 1,411,710 variants have been matched; 0 were flipped and 6 were reversed.
(info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))

# better to use af of GWAS and INFO scores as well (then can use 0.7 instead of 0.5 in L43)
# (cf. https://privefl.github.io/bigsnpr/articles/LDpred2.html#quality-control-of-gwas-summary-statistics)
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))

is_bad <-
  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05

library(ggplot2)
qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

df_beta <- info_snp[!is_bad, ]

# Here, you also want to restrict to the variants present
# in your test data as well. For this, you can use something like
in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], map_test[, c("chr", "pos")])
df_beta <- df_beta[in_test, ]

tmp <- tempfile(tmpdir = "tmp-data")

for (chr in 1:22) {

  cat(chr, ".. ", sep = "")

  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

  corr_chr <- readRDS(paste0("data/corr_hm3_plus/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# SNP-heritability estimation from LD score regression,
# to be used as a starting value in LDpred2-auto.
# Here `ld` is a pre-computed column of `map_ldref` and therefore `df_beta`.
# It corresponds to the pre-computed LD scores for the full set of HM3+ variants 
# (therefore the need for `ld_size = nrow(map_ldref)` instead of `length(ld)`).
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))
#    int  int_se     h2   h2_se
# 1.0377 0.00767 0.1693 0.01196
h2_est <- ldsc[["h2"]]


# LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
                               allow_jump_sign = FALSE, shrink_corr = 0.95,
                               ncores = NCORES) # 5 min

# Filter for best chains and average remaining ones
# -> the effects sizes of your polygenic score
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

# Use your favorite tool (e.g. PLINK) to get the polygenic score 'pred_auto'
# corresponding to the vectors of effects in 'beta_auto', e.g. here I can do
library(bigsnpr)
ukb <- snp_attach("data/UKBB_hm3_plus.rds")
G <- ukb$genotypes
map <- dplyr::transmute(ukb$map,
                        chr = as.integer(chromosome), pos = physical.pos,
                        a0 = allele1, a1 = allele2)  # reversed somehow..
map_pgs <- df_beta[1:4]; map_pgs$beta <- 1
map_pgs2 <- snp_match(map_pgs, map)

pred_auto <- big_prodVec(G, beta_auto * map_pgs2$beta,
                         ind.col = map_pgs2[["_NUM_ID_"]],
                         ncores = NCORES)

# Testing of final PGS
AUCBoot_no_NA <- function(pred, target) {
  is_not_na <- which(!is.na(target))
  AUCBoot(pred[is_not_na], target[is_not_na])
}
y <- readRDS("data/ukbb-binary-pheno/174.1.rds")[ukb$fam$ind_csv]
AUCBoot_no_NA(pred_auto, y + 0)
#   Mean   2.5%  97.5%      Sd
# 0.6579 0.6531 0.6627 0.00244

# cleanup
file.remove(paste0(tmp, ".sbk"))
