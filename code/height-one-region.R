library(bigsnpr)
library(dplyr)

NCORES <- 13


fam <- snp_attach(paste0("data/UKBB_hm3_plus.rds"))$fam %>%
  filter(set == "train")

# Get top hit of GWAS in HM3+
gwas0 <- purrr::map_dfr(paste0("GWAS/height_hm3_plus_part", 1:30, ".rds"), readRDS)
ind_top <- which.max(abs(gwas0$score))
map_hm3_plus <- readRDS("map_hm3_plus.rds")
stopifnot(nrow(map_hm3_plus) == nrow(gwas0))

map_hm3_plus[ind_top, ]
# chr       pos a0    a1    rsid      af_UKBB    ld block_id  pos_hg18  pos_hg38
#   3 141121814 A     C     rs2871960   0.448  27.1       50 142604504 141402972

# Use 1 Mb region around the top variant
pos <- map_hm3_plus$pos[ind_top]
chr <- map_hm3_plus$chr[ind_top]
mfi <- bigreadr::fread2(paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt"))
snp_id <- mfi %>%
  filter(V6 > 0.005, V8 > 0.5, between(V3, pos - 500e3, pos + 500e3)) %>%
  with(paste(chr, V3, V4, V5, sep = "_"))
length(snp_id)  # 3881

# Read dosage data
snp_readBGEN(
  bgenfiles = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = chr),
  list_snp_id = list(snp_id),
  backingfile = "tmp-data/height_one_region",
  ind_row = fam$ind_bgen,
  ncores = NCORES
)


ukbb <- snp_attach("tmp-data/height_one_region.rds")
G <- ukbb$genotypes
y <- readRDS("data/ukbb-quant-pheno/height.rds")[fam$ind_csv]
covar <- fam %>% select(PC1:date) %>% as.matrix()
ind.train <- which(!is.na(y) & complete.cases(covar))
length(ind.train)  # 305338

gwas <- big_univLinReg(G, y[ind.train], ind.train = ind.train,
                       covar.train = covar[ind.train, ], ncores = NCORES)

lpval <- predict(gwas)
library(ggplot2)
qplot(x = ukbb$map$physical.pos / 1e6, y = -lpval, alpha = I(0.5)) +
  labs(title = "Manhattan plot", x = glue::glue("Position (Mb, chr #{chr})"),
       y = expression(-log[10](italic("p-value")))) +
  theme_bigstatsr()
# ggsave("figures/manhattan_height_one_region.pdf", width = 13, height = 7)


bigparallelr::set_blas_ncores(NCORES)
corr0 <- big_cor(G, ind.row = ind.train)[]
corr <- corr0 %>% as("dgCMatrix") %>% as_SFBM(compact = TRUE)
df_beta <- transmute(gwas, beta = estim, beta_se = std.err, n_eff = length(ind.train))

# SuSiE-RSS
res <- with(df_beta, susieR::susie_rss(
  bhat = beta, shat = beta_se, R = corr0, n = n_eff[1]))  # L = 10 by default
res2 <- with(df_beta, susieR::susie_rss(
  bhat = beta, shat = beta_se, R = corr0, n = n_eff[1], L = 100))


# LDpred2-auto
# Heritability estimation of LD score regression to be used as a starting value in LDpred2-auto
(ldsc_h2 <- snp_ldsc2(corr, df_beta, ncores = NCORES)[["h2"]])  # 0.004631656

ldpred2_pip <- function(max_p = 1) {

  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc_h2, 0.001),
                                 vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                                 burn_in = 500, num_iter = 500, p_bounds = c(1e-5, max_p),
                                 allow_jump_sign = FALSE, shrink_corr = 0.98,
                                 use_MLE = FALSE, ncores = NCORES)

  # post-QC
  print(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
  (keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

  postp <- rowMeans(purrr::map_dfc(multi_auto[keep], "postp_est"))
}

postp1 <- ldpred2_pip(max_p = 1)
postp2 <- ldpred2_pip(max_p = 0.01)
postp3 <- ldpred2_pip(max_p = 0.001)
plot(data.frame(postp1, postp2, postp3))


res_all <-
  tibble(`1` = postp1, `0.01` = postp2, `0.001` = postp3, `10` = res$pip, `100` = res2$pip,
         pval = predict(gwas, log10 = FALSE)) %>%
  tidyr::pivot_longer(c(`10`, `100`), names_to = "L", values_to = "pip") %>%
  tidyr::pivot_longer(c(`1`, `0.01`, `0.001`), names_to = "max(p)", values_to = "postp") %>%
  arrange(desc(`max(p)`), desc(L))

p1 <- ggplot(res_all, aes(postp, pip, color = L)) +
  facet_grid(~ `max(p)`, labeller = label_both) +
  geom_abline(color = "red", linetype = 2) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  theme_bigstatsr(0.9) +
  coord_equal() +
  theme(legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = "PIPs with LDpred2-auto", y = "PIPs with SuSiE-RSS")

p2 <- ggplot(res_all, aes(postp, pip, color = `max(p)`)) +
  facet_grid(~ L, labeller = label_both) +
  geom_abline(color = "red", linetype = 2) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("#009E73", "#F0E442", "#CC79A7")) +
  theme_bigstatsr(0.9) +
  coord_equal() +
  theme(legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = "PIPs with LDpred2-auto", y = "PIPs with SuSiE-RSS")

cowplot::plot_grid(p1, p2, ncol = 1, labels = LETTERS[1:2], label_size = 16, scale = 0.9)
# ggsave("figures/pip_height_one_region.pdf", width = 9, height = 8)


top_pip <- bind_cols(ukbb$map, postp = postp2, pip = res$pip, gwas) %>%
  mutate(id = row_number()) %>%
  filter(postp > 0.3 | pip > 0.3) %>%
  select(-marker.ID) %>%
  print(n = Inf)
#  hromosome rsid             physical.pos allele1 allele2   freq  info postp   pip   estim std.err score    id
# 1 03         rs2871960           141121814 A       C       0.448  0.999 0.233 0.336  0.549   0.0162 33.9   1903
# 2 03         rs16851441          141173521 A       G       0.0905 0.991 0.259 0.347 -0.0555  0.0282 -1.97  2084
# 3 03         3:141302726_AC_A    141302726 AC      A       0.0128 0.816 0.753 0.550 -0.539   0.0796 -6.77  2398


## Top hits from GWAS

length(ind <- which(gwas$score > 33))  # 19
min(corr0[ind, ind]^2)  # 0.96 -> 19 top variants are highly correlated

sum(print(postp1[ind]))  # 5.77
sum(print(postp2[ind]))  # 2.39
sum(print(postp3[ind]))  # 1.31

sum(res$pip[ind])  # 1.02
sum(res2$pip[ind])  # 1.41


## Build credible sets from LDpred2's PIPs

postp <- postp2
all_cs <- list()
keep <- rep(TRUE, length(postp))
while(any(keep)) {
  ind_maxpostp <- which.max(postp * keep)
  if (postp[ind_maxpostp] < 0.1) break
  ind_inLD <- Matrix::which(corr[, ind_maxpostp]^2 > 0.25 & keep)
  ord <- order(postp[ind_inLD], decreasing = TRUE)
  cumprob <- cumsum(postp[ind_inLD][ord])
  n_enough <- which(cumprob > 0.95)
  if (length(n_enough) > 0) {
    ind_cs <- sort(ind_inLD[ord[seq_len(n_enough[1])]])
    keep[ind_cs] <- FALSE
    all_cs[[length(all_cs) + 1]] <- ind_cs
    print(ind_cs)
    for (j in ind_cs) {
      ind_inhighLD <- Matrix::which(corr[, j]^2 > 0.95)
      keep[ind_inhighLD] <- FALSE
    }
  } else {
    keep[ind_maxpostp] <- FALSE
  }
}
# [1] 1899 1905 2398
# [1] 1850 1871 1903 1911 1914
# [1] 2056 2057 2058 2061 2062 2073 2084
# [1] 1864 1917 1953 1955 1962 1967 1968 1970 1975 1977 1981 1983 1985 2006 2020 2021 2022 2037 2039 2040

res$sets$cs
# [1] 1843 1850 1854 1861 1870 1871 1874 1887 1891 1903 1911 1914
# [1] 2057 2058 2061 2084
# [1] 1899 1905 2297 2341 2398
