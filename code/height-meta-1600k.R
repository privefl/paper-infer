library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)


map <- transmute(snp_attach("data/UKBB_hm3_plus.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, rsid, af_val = freq,
                 a0 = allele2, a1 = allele1)


gz <- runonce::download_file(
  "https://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz",
  dir = "tmp-data")
sumstats <- fread2(
  gz,
  select = c("CHR", "POS", "RSID", "EFFECT_ALLELE", "OTHER_ALLELE", "EFFECT_ALLELE_FREQ", "BETA", "SE", "N"),
  col.names = c("chr", "pos", "rsid", "a0", "a1", "freq", "beta", "beta_se", "n_eff"))


# Ancestry inference
all_freq <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31620968",
                         dir = "tmp-data", fname = "ref_freqs.csv.gz"))
projection <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31620953",
                         dir = "tmp-data", fname = "projection.csv.gz"))
matched <- snp_match(sumstats, all_freq[1:5], return_flip_and_rev = TRUE) %>%
  mutate(freq = ifelse(`_REV_`, 1 - freq, freq))

res <- snp_ancestry_summary(
  freq = 1 - matched$freq,  # need to switch ref alleles
  info_freq_ref = all_freq[matched$`_NUM_ID_`, -(1:5)],
  projection = projection[matched$`_NUM_ID_`, -(1:5)],
  correction = c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                 1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)
)
# Note that some ancestry groups from the reference are very close to one another,
# and should be merged a posteriori.
group <- colnames(all_freq)[-(1:5)]
group[group %in% c("Scandinavia", "United Kingdom", "Ireland")]   <- "Europe (North West)"
group[group %in% c("Europe (South East)", "Europe (North East)")] <- "Europe (East)"
grp_fct <- factor(group, unique(group))
final_res <- tapply(res, grp_fct, sum)
round(100 * final_res[final_res > 0.001], 1)
# Africa (West)   Ashkenazi    Europe (East)    Finland  Europe (North West)  Europe (South West)
#           0.2         1.5              9.5        6.5                 81.9                  0.3


# Matching and QC
info_snp <- as_tibble(snp_match(sumstats, map, return_flip_and_rev = TRUE))
# 1,373,020 variants to be matched.
# 11 ambiguous SNPs have been removed.
# 1,052,908 variants have been matched; 1 were flipped and 5 were reversed.

hist(info_snp$n_eff)
max(info_snp$n_eff)  # 1,597,374

info_snp2 <- info_snp %>%
  filter(n_eff > (0.8 * max(n_eff)), !`_FLIP_`, !`_REV_`) %>%
  mutate(sd_af = sqrt(2 * freq * (1 - freq)),
         sd_ss = 1 / sqrt(n_eff * beta_se^2 + beta^2),
         sd_ss = sd_ss * print(sqrt(0.5) / quantile(sd_ss, 0.99)))


info_snp2$is_bad <- with(info_snp2,
                         sd_ss < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
                           sd_ss < 0.05 | sd_af < 0.05)
mean(info_snp2$is_bad) # 0.0001489666


qplot(sd_af, sd_ss, color = ifelse(is_bad, "Yes", "No"), alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr(0.8) +
  coord_equal() +
  theme(legend.position = c(0.2, 0.8)) +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

df_beta <- filter(info_snp2, !is_bad)
nrow(df_beta)  # 1,013,499


corr <- runonce::save_run({

  for (chr in 1:22) {

    cat(chr, ".. ", sep = "")

    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ## indices in 'corr_chr'
    ind.chr3 <- match(ind.chr2, which(map$chr == chr))

    corr_chr <- readRDS(paste0("data/corr_hm3_plus/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]

    if (chr == 1) {
      corr <- as_SFBM(corr_chr, "tmp-data/corr_height", compact = TRUE)
    } else {
      corr$add_columns(corr_chr, nrow(corr))
    }
  }
  corr
}, file = "tmp-data/corr_height.rds")


NCORES <- 13

# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
(ldsc <- snp_ldsc2(corr, df_beta, blocks = 200, intercept = NULL, ncores = NCORES))
#        int     int_se         h2      h2_se
# 2.31384970 0.06809280 0.39210820 0.01719189

# LDpred2-auto
coef_shrink <- 0.95
multi_auto <- runonce::save_run(
  snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                   vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                   burn_in = 500, num_iter = 500, report_step = 20,
                   allow_jump_sign = TRUE, shrink_corr = coef_shrink,
                   ncores = NCORES),
  file = "tmp-data/mod_height_1600K.rds"
)

# derive other results
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)) & range < 3))

hist(all_h2    <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est,    500)))
unname(round(100 * quantile(all_h2, c(0.5, 0.025, 0.975)), 1))  # 54.2 [53.9, 54.5]
hist(all_p     <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est,     500)))
unname(round(100 * quantile(all_p, c(0.5, 0.025, 0.975)), 1))  # 5.9 [5.6, 6.3]
hist(all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500)))
unname(round(quantile(all_alpha, c(0.5, 0.025, 0.975)), 2))  # -0.78 [-0.82, -0.76]

# inferred r2
bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
  print(ic)
  b1 <- bsamp[[ic]]
  Rb1 <- apply(b1, 2, function(x)
    coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
  b2 <- do.call("cbind", bsamp[-ic])
  b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
}))
hist(all_r2)
unname(round(100 * quantile(all_r2, c(0.5, 0.025, 0.975)), 1)) # 47.0 [46.8, 47.1]


postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp))
sum(postp > 0.95)  # 1753

qplot(y = postp, color = as.factor(df_beta$chr), alpha = I(0.1)) +
  theme_bw(13) +
  scale_color_manual(values = rep_len(c("#000000", "#999999"), 22)) +
  scale_x_continuous(breaks = seq(0, 1.1e6, by = 100e3)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(legend.position = "none") +
  labs(x = "Index of variant",
       y = "Posterior probability of being causal (for height)")
# ggsave("figures/postp_height.png", width = 12, height = 6)


untar(runonce::download_file(
  "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz",
  dir = "tmp-data"), exdir = "tmp-data/baseline_annot")
annot <- bigreadr::fread2(paste0("tmp-data/baseline_annot/baselineLD.", 1:22, ".annot.gz"))

ind <- vctrs::vec_match(map[df_beta$`_NUM_ID_`, 1:2], transmute(annot, chr = CHR, pos = BP))

var <- annot %>%
  sapply(function(x) identical(sort(unique(x)), 0:1)) %>%
  which() %>%
  names() %>%
  grep(pattern = "flanking", value = TRUE, invert = TRUE) %>%
  print()
# 50

library(future.apply)
all_ind1 <- sapply(annot[var][ind, ], function(one_annot) which(one_annot == 1))
plan("multisession", workers = 15)
options(future.globals.maxSize = 1000 * 1024^2)

bsamp <- do.call("cbind", bsamp)
all_enrich <- future_sapply(cols_along(bsamp), function(j) {

  x_all <- bsamp[, j]
  h2_all <- coef_shrink * crossprod(x_all, bigsparser::sp_prodVec(corr, x_all)) +
    (1 - coef_shrink) * crossprod(x_all)

  zero <- x_all * 0

  enrich <- sapply(all_ind1, function(ind1) {
    x <- zero; x[ind1] <- x_all[ind1]
    h2_1 <- coef_shrink * crossprod(x, bigsparser::sp_prodVec(corr, x)) +
      (1 - coef_shrink) * crossprod(x)
    (h2_1 / length(ind1)) / (h2_all / length(x_all))
  })
})

names <- paste0(names(all_ind1), " [", round(100 * lengths(all_ind1) / length(ind), 1), "%]")
res_ldpred2 <- all_enrich %>%
  apply(1, quantile, probs = c(0.5, 0.025, 0.975)) %>%
  as.data.frame() %>%
  setNames(names) %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols = -rowname) %>%
  tidyr::pivot_wider(names_from = "rowname", values_from = "value") %>%
  mutate(name = factor(name, levels = names[order(colMeans(all_enrich))])) %>%
  print()

# saveRDS(res_ldpred2, "tmp-data/mod_height_enrich_meta.rds")
