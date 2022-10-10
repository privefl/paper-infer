library(bigsnpr)
library(dplyr)

pheno_file <- "data/ukbb-quant-pheno/height.rds"
NCORES <- 13


# Get GWAS and its sample size
pheno <- sub("\\.rds$", "", basename(pheno_file))
gwas <- purrr::map_dfr(paste0("GWAS/", pheno, "_hm3_plus_part", 1:30, ".rds"), readRDS)
N <- as.integer(sub(".+, df = ([0-9]+), .+", "\\1", body(attr(gwas, "predict"))[2]))

df_beta <- gwas %>%
  transmute(beta = estim, beta_se = std.err, n_eff = N)


corr <- readRDS("data/corr_hm3_plus.rds")


# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
(ldsc <- snp_ldsc2(corr, df_beta, blocks = 200, intercept = NULL, ncores = NCORES))
#        int     int_se         h2      h2_se
# 1.11317774 0.01532732 0.59707067 0.02219329


# LDpred2-auto
coef_shrink <- 0.95
multi_auto <- runonce::save_run(
  snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                   vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                   burn_in = 500, num_iter = 500, report_step = 20,
                   allow_jump_sign = TRUE, shrink_corr = coef_shrink,
                   ncores = NCORES),
  file = "tmp-data/mod_height.rds"
)

# derive other results
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)) & range < 3))

hist(all_h2    <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est,    500)))
unname(round(100 * quantile(all_h2, c(0.5, 0.025, 0.975)), 1))  # 63.2 [62.0, 64.4]
hist(all_p     <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est,     500)))
unname(round(100 * quantile(all_p, c(0.5, 0.025, 0.975)), 1))  # 2.3 [2.0, 2.5]
hist(all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500)))
unname(round(quantile(all_alpha, c(0.5, 0.025, 0.975)), 2))  # -0.74 [-0.76, -0.72]

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
unname(round(100 * quantile(all_r2, c(0.5, 0.025, 0.975)), 1)) # 42.7 [42.2, 43.1]


untar(runonce::download_file(
  "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz",
  dir = "tmp-data"), exdir = "tmp-data/baseline_annot")
annot <- bigreadr::fread2(paste0("tmp-data/baseline_annot/baselineLD.", 1:22, ".annot.gz"))

map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
  transmute(chr = as.integer(chromosome), pos = physical.pos)
ind <- vctrs::vec_match(map, transmute(annot, chr = CHR, pos = BP))


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
options(future.globals.maxSize = 600 * 1024^2)

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

same_names <- function(name) {
  res_ldpred2$name[match(sub("(.+) \\[.+\\]", "\\1", name),
                         sub("(.+) \\[.+\\]", "\\1", res_ldpred2$name))]
}

both_res_ldpred2 <- bind_rows(
  bind_cols(N = "305K", res_ldpred2),
  # from 'enrichment-height-100k.R':
  mutate(readRDS("tmp-data/mod_height_enrich_100k.rds"),
         N = "100K", name = same_names(name)),
  # from 'enrichment-height-meta.R':
  mutate(readRDS("tmp-data/mod_height_enrich_meta.rds"),
         N = "1.6M", name = same_names(name))
)

library(ggplot2)
both_res_ldpred2 %>%
  mutate(name = factor(name, levels = with(res_ldpred2, name[order(`50%`)])),
         N = factor(N, levels = c("100K", "305K", "1.6M"))) %>%
  ggplot(aes(name, `50%`, color = N)) +
  theme_bigstatsr(0.8) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(position = position_dodge(width = 0.7)) +
  coord_flip() +
  scale_color_manual(values = c("#E69F00", "#0072B2", "#CC79A7"),
                     guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(minor_breaks = seq(0, 10, by = 0.1)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.5,
                position = position_dodge(width = 0.7)) +
  theme(legend.position = c(0.75, 0.15)) +
  labs(y = "Heritability enrichment", x = NULL, color = "GWAS sample size")
# ggsave("figures/height_enrichment.pdf", width = 10, height = 11)
