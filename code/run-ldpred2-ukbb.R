library(bigsnpr)
library(dplyr)

SETS <- c(
  "hm3",        # HM3 set from LDpred2
  "hm3_plus",   # extending prev HM3 to fill gaps of tagging
  "hm3_small",  # HM3 but with corr based on 2000 indiv only
  "hm3_altpop"  # HM3 but with corr based on 10K indiv from around South Europe
)

# build SFBM (same for all pheno)
for (set in SETS) {

  corr <- runonce::save_run({

    for (chr in 1:22) {

      cat(chr, ".. ", sep = "")

      if (set == "hm3") {
        corr_chr <- readRDS(paste0("../misspec/ldref/LD_with_blocks_chr", chr, ".rds"))
      } else {
        corr_chr <- readRDS(paste0("data/corr_", set, "/LD_with_blocks_chr", chr, ".rds"))
      }

      if (chr == 1) {
        corr <- as_SFBM(corr_chr, paste0("data/corr_", set), compact = TRUE)
      } else {
        corr$add_columns(corr_chr, nrow(corr))
      }
    }
    corr
  }, file = paste0("data/corr_", set, ".rds"))
  print(dim(corr))
}
# 1054330 x 1054330
# 1444196 x 1444196
# 1054330 x 1054330
# 1054330 x 1054330
rm(corr)

pheno_files <- c(list.files("data/ukbb-quant-pheno", full.names = TRUE),
                 list.files("data/ukbb-binary-pheno", full.names = TRUE))

bigassertr::assert_dir("ldpred2")

grid <- tibble(pheno_file = pheno_files) %>%
  # slice_sample(n = 10) %>%
  tidyr::expand_grid(set = SETS, use_mle = c(TRUE, FALSE), jump_sign = c(TRUE, FALSE)) %>%
  filter(set == "hm3" | !jump_sign) %>%
  mutate(res_file = paste0("ldpred2/", set, ifelse(use_mle, "_MLE", "_noMLE"),
                           ifelse(jump_sign, "_jump", "_nojump"),
                           "_", basename(pheno_file))) %>%
  # filter(!file.exists(res_file)) %>%
  print(n = Inf)


NCORES <- 13
library(future.batchtools)
plan(workers = nrow(grid) + 10, batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))
# future::plan("multisession", workers = 14, gc = TRUE)

res_grid <- furrr::future_pmap_dfr(grid, function(pheno_file, set, use_mle, jump_sign, res_file) {

  # pheno_file <- "data/ukbb-quant-pheno/log_lipoA.rds"
  # set <- "hm3_small"
  # res_file <- paste0("ldpred2/", set, "_", basename(pheno_file))

  res <- runonce::save_run(file = res_file, {

    # Get GWAS and its sample size
    pheno <- sub("\\.rds$", "", basename(pheno_file))
    gwas <- if (set == "hm3_plus") {
      purrr::map_dfr(paste0("GWAS/", pheno, "_hm3_plus_part", 1:30, ".rds"), readRDS)
    } else {
      purrr::map_dfr(paste0("GWAS/", pheno, "_hm3_part",      1:30, ".rds"), readRDS)
    }
    gwas[is.na(gwas$score), ] <- list(0, 1, 0)  # not ideal, but only 9 SNPs for 1 pheno
    N <- as.integer(sub(".+, df = ([0-9]+), .+", "\\1", body(attr(gwas, "predict"))[2]))

    df_beta <- gwas %>%
      transmute(beta = estim, beta_se = std.err, n_eff = N)


    corr <- readRDS(paste0("data/corr_", set, ".rds"))


    # Heritability estimation of LD score regression
    # to be used as a starting value in LDpred2-auto
    (ldsc <- snp_ldsc2(corr, df_beta, blocks = 200, intercept = NULL, ncores = NCORES))


    # LDpred2-auto
    coef_shrink <- 0.95
    time <- system.time(
      multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                                     vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                                     burn_in = 500, num_iter = 500, report_step = 20,
                                     use_MLE = use_mle,
                                     allow_jump_sign = jump_sign,
                                     shrink_corr = coef_shrink,
                                     ncores = NCORES)
    )

    # derive other results
    (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
    (keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)) & range < 3))

    if (length(keep) < 2) {

      list(NULL)

    } else {

      (all_h2    <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est,    500)))
      (all_p     <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est,     500)))
      (all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500)))

      # inferred r2
      bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
      all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
        # print(ic)
        b1 <- bsamp[[ic]]
        Rb1 <- apply(b1, 2, function(x)
          coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
        b2 <- do.call("cbind", bsamp[-ic])
        b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
      }))
      # hist(all_r2)

      beta <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
      ukbb <- snp_attach(`if`(set == "hm3_plus", paste0("data/UKBB_hm3_plus.rds"),
                              paste0("data/UKBB_hm3.rds")))
      G <- ukbb$genotypes
      ind_test <- which(ukbb$fam$set == "test")
      pred <- big_prodVec(G, beta, ind.row = ind_test, ncores = NCORES)
      y_test <- readRDS(pheno_file)[ukbb$fam$ind_csv[ind_test]]
      covar_test <- dplyr::select(ukbb$fam, -(1:4))[ind_test, ]
      cor <- pcor(pred, y_test, covar_test)
      boot_cor <- replicate(1000, {
        # cat(".")
        ind <- sample(length(pred), replace = TRUE)
        pcor(pred[ind], y_test[ind], covar_test[ind, ])[1]
      })
      cor[2:3] <- quantile(boot_cor, c(0.025, 0.975))
      (r2 <- sign(cor) * cor^2)

      # local h2 (per-block)
      blocks <- if (set == "hm3") {
        readRDS("tmp-data/map_hm3_ldpred2.rds")$group_id
      } else {
        all_size <- unlist(readRDS(paste0("data/corr_", set, "/all_final_grp.rds"))$all_size)
        rep(seq_along(all_size), all_size)
      }
      list_ind <- split(seq_along(blocks), blocks)
      bsamp <- do.call("cbind", bsamp)
      all_local_h2 <- sapply(list_ind, function(ind) {
        corr_sub <- corr[ind, ind]
        bsamp_sub <- bsamp[ind, ]
        Rb <- coef_shrink * corr_sub %*% bsamp_sub + (1 - coef_shrink) * bsamp_sub
        Matrix::colSums(bsamp_sub * Rb)
      })

      # per-variant prob causal
      postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp_est))

      list2 <- function(...) setNames(list(...), sapply(substitute(list(...))[-1], deparse))
      list2(ldsc, postp, all_h2, all_p, all_alpha, all_r2, r2, all_local_h2, time)
    }
  })

  if (identical(res, list(NULL))) {
    tibble(n_keep = 0)
  } else {
    list(res) %>%
      tibble(res = .) %>%
      tidyr::unnest_wider(res) %>%
      rowwise() %>%
      mutate(
        n_keep = ncol(all_h2),
        local_h2 = list(colMeans(all_local_h2)),
        h2_ldsc_est = list(ldsc[["h2"]] + c(0, -1.96, 1.96) * ldsc[["h2_se"]]),
        across(c(h2_est = all_h2, p_est = all_p, alpha_est = all_alpha, r2_est = all_r2),
               ~ list(unname(quantile(., probs = c(0.5, 0.025, 0.975), na.rm = TRUE))))
      ) %>%
      select(-ldsc, -(all_h2:all_r2), -all_local_h2)
  }
})

gc(reset = TRUE)
c(nrow(grid), nrow(res_grid))


library(dplyr)
library(ggplot2)
theme_set(theme_bw(14))

grid2 <- bind_cols(grid, res_grid) %>%
  filter(!grepl("rint_", pheno_file), !grepl("raw_", pheno_file), !jump_sign)

with(grid2, table(set, use_mle))

qplot(data = grid2, n_keep, geom = "density", color = set,
      linetype = ifelse(use_mle, "Yes", "No")) +
  labs(x = "# chains kepts", y = "Density",
       color = "Set of variants", linetype = "Use MLE?") +
  theme(legend.position = c(0.35, 0.75), legend.box = "horizontal")
# ggsave("figures/ukbb_chains.pdf", width = 7, height = 5)

qplot(data = filter(grid2, set %in% c("hm3", "hm3_plus")),
      purrr::map_dbl(time, "elapsed") / 60, geom = "density",
      color = set, linetype = ifelse(use_mle, "Yes", "No")) +
  xlim(0, NA) +
  labs(x = "Runtime  (in minutes)", y = "Density", color = "Set of variants",
       linetype = "Use MLE?") +
  theme(legend.position = c(0.75, 0.7))
# ggsave("figures/ukbb_runtimes.pdf", width = 7, height = 5)

median(purrr::map_dbl(filter(grid2, set == "hm3_plus")$time, 3)) /
  median(purrr::map_dbl(filter(grid2, set == "hm3")$time, 3))
# 1.45

grid3 <- grid2 %>%
  filter(n_keep > 0) %>%
  tidyr::unnest_wider("r2", names_sep = "_") %>%
  tidyr::unnest_wider("r2_est", names_sep = "_") %>%
  tidyr::unnest_wider("h2_ldsc_est", names_sep = "_") %>%
  tidyr::unnest_wider("h2_est", names_sep = "_") %>%
  tidyr::unnest_wider("p_est", names_sep = "_") %>%
  tidyr::unnest_wider("alpha_est", names_sep = "_")


grid3_hm3 <- filter(grid3, set == "hm3")

grid3_hm3 %>%
  mutate(h2 = cut(h2_est_1, c(0, 0.05, 0.1, 0.2, 1))) %>%
  ggplot(aes(r2_1, r2_est_1, color = n_keep)) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_abline(color = "chartreuse3", linetype = 2) +
  geom_errorbar(aes(ymin = r2_est_2, ymax = r2_est_3),
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(xmin = r2_2, xmax = r2_3),
                position = position_dodge(width = 0.9), width = 0) +
  ggrepel::geom_label_repel(aes(label = sub("\\.rds$", "", basename(pheno_file))),
                            color = "purple", min.segment.length = 0, force = 500, size = 3, seed = 4,
                            data = mutate(filter(grid3_hm3, grepl("height", pheno_file)),
                                          h2 = cut(h2_est_1, c(0, 0.05, 0.1, 0.2, 1)))) +
  facet_wrap(~ use_mle + h2, scales = "free", labeller = label_both, nrow = 2) +
  geom_point(size = 1.5) +
  theme(legend.position = "top", legend.key.width = unit(40, "pt")) +
  guides(color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred r2", x = "r2 in test set", color = "# chains")
# ggsave("figures/ukbb_r2.pdf", width = 8.5, height = 6)

grid3_hm3 %>%
  filter(grepl("height", pheno_file) | grepl("/years", pheno_file)) %>%
  select(res_file, starts_with("r2_")) %>%
  print(width = Inf)
#   res_file                                     r2_1  r2_2  r2_3 r2_est_1 r2_est_2 r2_est_3
#  1 ldpred2/hm3_MLE_nojump_F_height.rds         0.317  0.308  0.326    0.352    0.345    0.358
#  2 ldpred2/hm3_noMLE_nojump_F_height.rds       0.314  0.306  0.324    0.360    0.355    0.365
#  3 ldpred2/hm3_MLE_nojump_height.rds           0.369  0.363  0.377    0.399    0.394    0.404
#  4 ldpred2/hm3_noMLE_nojump_height.rds         0.367  0.360  0.373    0.406    0.402    0.410
#  5 ldpred2/hm3_MLE_nojump_M_height.rds         0.313  0.302  0.323    0.322    0.315    0.329
#  6 ldpred2/hm3_noMLE_nojump_M_height.rds       0.309  0.299  0.319    0.330    0.324    0.336
#  7 ldpred2/hm3_MLE_nojump_sitting_height.rds   0.238  0.231  0.245    0.247    0.243    0.251
#  8 ldpred2/hm3_noMLE_nojump_sitting_height.rds 0.236  0.229  0.242    0.253    0.250    0.256
#  9 ldpred2/hm3_MLE_nojump_years_of_edu.rds     0.0336 0.0305 0.0371   0.0334   0.0313   0.0356
# 10 ldpred2/hm3_noMLE_nojump_years_of_edu.rds   0.0331 0.0297 0.0365   0.0343   0.0329   0.0356

grid3_hm3 %>%
  mutate(p = cut(p_est_1, c(1e-5, 1e-3, 0.01, 0.04))) %>%
  ggplot(aes(h2_ldsc_est_1, h2_est_1, color = n_keep)) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_abline(color = "chartreuse3", linetype = 2) +
  geom_errorbar(aes(ymin = h2_est_2, ymax = h2_est_3),
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(xmin = h2_ldsc_est_2, xmax = h2_ldsc_est_3),
                position = position_dodge(width = 0.9), width = 0) +
  facet_wrap(~ p + use_mle, labeller = label_both, ncol = 2) +
  coord_equal() +
  geom_point(size = 1.5) +
  theme(legend.key.height = unit(25, "pt")) +
  scale_x_continuous(breaks = 0:10 / 10) + scale_y_continuous(breaks = 0:10 / 10) +
  labs(y = "Inferred h2 with LDpred2-auto",
       x = "Inferred h2 with LDSc regression", color = "# chains")
# ggsave("figures/ukbb_h2.pdf", width = 8, height = 6.5)


p1 <- grid3_hm3 %>%
  filter(use_mle) %>%
  ggplot(aes(p_est_1, h2_est_1, color = n_keep)) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_vline(xintercept = 1e-5, color = "grey60", linetype = 3) +
  geom_hline(yintercept = 1e-3, color = "grey60", linetype = 3) +
  geom_errorbar(aes(ymin = h2_est_2, ymax = h2_est_3),
                position = position_dodge(width = 0.9), width = 0, alpha = 0.6) +
  geom_errorbar(aes(xmin = p_est_2, xmax = p_est_3),
                position = position_dodge(width = 0.9), width = 0, alpha = 0.6) +
  geom_point(size = 1.5) +
  scale_x_log10(breaks = signif(seq_log(1e-5, 1, 11), 1),
                minor_breaks = unique(signif(seq_log(1e-5, 1, 100), 1))) +
  scale_y_log10(breaks = signif(seq_log(1e-5, 1, 11), 1),
                minor_breaks = unique(signif(seq_log(1e-5, 1, 100), 1))) +
  labs(y = "Inferred h2  (log-scale)", x = "Inferred p  (log-scale)", color = "# chains") +
  theme(legend.position = c(0.13, 0.78), legend.key.height = unit(20, "pt"))
(p1.2 <- ggExtra::ggMarginal(p1, type = "histogram", margins = "both"))
# ggsave("figures/ukbb_h2_p.pdf", p1.2, width = 8, height = 6)

p2 <- grid3_hm3 %>%
  filter(use_mle) %>%
  filter(n_keep > 25) %>%
  ggplot(aes(alpha_est_1, h2_est_1, color = n_keep)) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_vline(xintercept = c(-1.5, 0.5), color = "grey60", linetype = 3) +
  geom_errorbar(aes(ymin = h2_est_2, ymax = h2_est_3),
                position = position_dodge(width = 0.9), width = 0, alpha = 0.6) +
  geom_errorbar(aes(xmin = alpha_est_2, xmax = alpha_est_3),
                position = position_dodge(width = 0.9), width = 0, alpha = 0.6) +
  geom_point(size = 1.5) +
  scale_y_log10(breaks = signif(seq_log(1e-5, 1, 11), 1), limits = c(NA, 1),
                minor_breaks = unique(signif(seq_log(1e-5, 1, 100), 1))) +
  scale_x_continuous(minor_breaks = seq(-1.5, 0.5, by = 0.1)) +
  theme(legend.position = c(0.75, 0.9), legend.key.width = unit(20, "pt"),
        legend.direction = "horizontal") +
  labs(y = "Inferred h2  (log-scale)", x = "Inferred alpha", color = "# chains")
(p2.2 <- ggExtra::ggMarginal(p2, type = "histogram", margins = "both"))
# ggsave("figures/ukbb_h2_alpha.pdf", p2.2, width = 8, height = 6)


#### Comparing hm3 vs hm3_plus and alternative LD ####

grid4 <- grid2 %>%
  rowwise() %>%
  filter(set != "hm3" | n_keep > 25,
         !sub("\\.rds$", "", basename(pheno_file)) %in% c("F_height", "M_height")) %>%
  mutate(across(c(r2, r2_est), ~ list(`if`(n_keep == 0, rep(0, 3), .))),
         across(c(alpha_est, p_est, h2_est, h2_ldsc_est),
                ~ list(`if`(n_keep == 0, rep(NA_real_, 3), .)))) %>%
  ungroup() %>%
  tidyr::pivot_longer(c(r2, r2_est, alpha_est, p_est, h2_est, h2_ldsc_est)) %>%
  tidyr::pivot_wider(c(pheno_file, use_mle, name), names_from = set, values_from = value) %>%
  tidyr::unnest_wider("hm3", names_sep = "_") %>%
  tidyr::unnest_wider("hm3_plus", names_sep = "_") %>%
  tidyr::unnest_wider("hm3_small", names_sep = "_") %>%
  tidyr::unnest_wider("hm3_altpop", names_sep = "_")


grid4 %>%
  filter(use_mle) %>%
  ggplot(aes(hm3_1, hm3_plus_1)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_errorbar(aes(ymin = hm3_plus_2, ymax = hm3_plus_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(xmin = hm3_2, xmax = hm3_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_abline(color = "red", linetype = 2) +
  geom_point(size = 1.2) +
  facet_wrap(~ name, scales = "free") +
  labs(x = "Estimates with HapMap3 variants", y = "Estimates with HapMap3+ variants")
# ggsave("figures/ukbb_compare_hm3_plus.pdf", width = 9, height = 6)

grid4 %>% filter(grepl("lipoA", pheno_file)) %>% select(2:9) %>% mutate_if(is.numeric, round, digits = 3)
#    use_mle name         hm3_1  hm3_2  hm3_3 hm3_plus_1 hm3_plus_2 hm3_plus_3
#  1 TRUE    r2           0.344  0.336  0.352      0.516      0.509      0.524
#  2 TRUE    r2_est       0.32   0.316  0.323      0.502      0.497      0.506
#  3 TRUE    alpha_est   -0.931 -1.24  -0.503     -0.868     -1.19      -0.547
#  4 TRUE    p_est        0      0      0          0          0          0
#  5 TRUE    h2_est       0.324  0.32   0.329      0.508      0.502      0.514
#  6 TRUE    h2_ldsc_est  0.203 -0.131  0.536      0.287     -0.214      0.788
#  7 FALSE   r2           0.344  0.336  0.353      0.516      0.508      0.524
#  8 FALSE   r2_est       0.319  0.316  0.323      0.501      0.497      0.506
#  9 FALSE   alpha_est   NA     NA     NA         NA         NA         NA
# 10 FALSE   p_est        0      0      0          0          0          0
# 11 FALSE   h2_est       0.324  0.32   0.329      0.508      0.503      0.515
# 12 FALSE   h2_ldsc_est  0.203 -0.131  0.536      0.287     -0.214      0.788

grid4 %>%
  filter(name == "r2", !grepl("lipoA", pheno_file), use_mle) %>%
  deming::deming(hm3_plus_1 ~ hm3_1 + 0, data = .,
                 xstd = hm3_3 - hm3_2, ystd = hm3_plus_3 - hm3_plus_2)
# Slope     1.061172 0.01049487   1.040602   1.081741

grid4 %>%
  filter(use_mle) %>%
  ggplot(aes(hm3_1, hm3_small_1)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_errorbar(aes(ymin = hm3_small_2, ymax = hm3_small_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(xmin = hm3_2, xmax = hm3_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_abline(color = "red", linetype = 2) +
  geom_point(size = 1.2) +
  facet_wrap(~ name, scales = "free") +
  labs(x = "Estimates with HapMap3 variants",
       y = "Estimates with HapMap3 variants (and small LD)")
# ggsave("figures/ukbb_compare_hm3_small.pdf", width = 9, height = 6)

grid4 %>% filter(name == "r2", hm3_small_1 == 0) %>% pull(pheno_file)
# log_ALP -- log_lipoA -- MCH


grid4 %>%
  filter(use_mle) %>%
  ggplot(aes(hm3_1, hm3_altpop_1)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_errorbar(aes(ymin = hm3_altpop_2, ymax = hm3_altpop_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(xmin = hm3_2, xmax = hm3_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_abline(color = "red", linetype = 2) +
  geom_point(size = 1.2) +
  facet_wrap(~ name, scales = "free") +
  labs(x = "Estimates with HapMap3 variants",
       y = "Estimates with HapMap3 variants (and S.Eur. LD)")
# ggsave("figures/ukbb_compare_hm3_altpop.pdf", width = 9, height = 6)

grid4 %>% filter(name == "r2", hm3_altpop_1 == 0) %>% pull(pheno_file)
# log_leukocyte -- log_lymphocyte -- log_monocyte -- 557.1


#### Local h2 in HapMap3+ ####

res_local_h2 <- grid2 %>%
  filter(set == "hm3", use_mle) %>%
  rowwise() %>%
  mutate(h2 = sum(local_h2), max_local_h2 = max(local_h2)) %>%
  ungroup() %>%
  mutate(pheno = sub("\\.rds$", "", basename(pheno_file)),
         pheno = factor(pheno, levels = pheno[order(h2)]),
         h2 > 0.1) %>%
  filter((max_local_h2 / h2) > 0.1) %>%
  tidyr::pivot_longer(c(h2, max_local_h2))

library(ggplot2)
ggplot(res_local_h2, aes(pheno, value, fill = name, color = name)) +
  theme_bw(13) +
  scale_color_manual(values = c("#0072B2FF", "#00000033")) +  # "#E69F0000"
  scale_fill_manual(values = c("#0072B200", "#E69F0099")) +
  geom_col(position = "identity") +
  facet_wrap(~ `h2 > 0.1`, scales = "free") +
  theme(legend.position = c(0.15, 0.75), strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Phenotype", y = NULL, fill = NULL, color = NULL)
# ggsave("figures/ukbb_local_h2.pdf", width = 11, height = 6)

res_median_local_h2 <- grid2 %>%
  filter(set == "hm3_plus", n_keep > 25, use_mle) %>%
  print() %>%
  pull(local_h2) %>%
  simplify2array() %>%
  matrixStats::rowMedians() %>%
  unname()

all_size <- unlist(readRDS(paste0("data/corr_hm3_plus/all_final_grp.rds"))$all_size)

map_hm3_plus <- snp_attach("data/UKBB_hm3_plus.rds")$map
last <- cumsum(all_size)
first <- head(c(1L, last + 1L), -1)
label <- glue::glue("chr{chr} [{beg}-{end}]",
                    chr = as.integer(map_hm3_plus$chromosome[first]),
                    beg = round(map_hm3_plus$physical.pos[first] / 1e6, 1),
                    end = round(map_hm3_plus$physical.pos[last] / 1e6, 1))
label[which.max(res_median_local_h2)]  # chr6 [22.1-41.4] -> HLA

qplot(all_size, res_median_local_h2) + #, label = label) + +
  # ggrepel::geom_label_repel(color = "purple", min.segment.length = 0) +
  labs(x = "Number of variants in block",
       y = "Median heritability of block (across 169 phenotypes)")
# ggsave("figures/median_local_h2.pdf", width = 8, height = 6)


# posterior probs of being causal
res_median_postp <- grid2 %>%
  filter(set == "hm3_plus", n_keep > 25) %>%
  ungroup() %>%
  pull(postp) %>%
  simplify2array() %>%
  matrixStats::rowMedians() %>%
  print()

(ind_keep <- which(res_median_postp > 0.01))
info <- map_hm3_plus %>%
  slice(ind_keep) %>%
  print(n = Inf) %>%
  pull(rsid) %>%
  rsnps::ncbi_snp_query() %>%
  select(1:6) %>%
  print(n = Inf)

df_for_label <- info %>%
  mutate(id = ind_keep, postp = res_median_postp[ind_keep]) %>%
  filter(gene != "") %>%
  group_by(chromosome) %>%
  arrange(desc(postp)) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(id) %>%
  print(n = Inf)
#    query     chromosome        bp class rsid      gene             id  postp
#  1 rs1260326 2           27508073 snv   rs1260326 GCKR         133360 0.0132
#  2 rs724016  3          141386728 snv   rs724016  ZBTB38       310045 0.0110
#  3 rs256903  5           56513311 snv   rs256903  C5orf67      460613 0.0103
#  4 rs1800562 6           26092913 snv   rs1800562 HFE/HFE-AS1  541003 0.0112
#  5 rs1051921 7           73593613 snv   rs1051921 MLXIPL       657705 0.0106
#  6 rs1461729 8            9329732 snv   rs1461729 LOC157273    707827 0.0100
#  7 rs505922  9          133273813 snv   rs505922  ABO          840626 0.0121
#  8 rs7073746 10          63144311 snv   rs7073746 NRBF2        878517 0.0113
#  9 rs174537  11          61785208 snv   rs174537  MYRF         950291 0.0132
# 10 rs653178  12         111569952 snv   rs653178  ATXN2       1048228 0.0139
# 11 rs3751812 16          53784548 snv   rs3751812 FTO         1227728 0.0115
# 12 rs4420638 19          44919689 snv   rs4420638 APOC1       1357252 0.0121

# Numbers of traits associated with gene (according to the GWAS Catalog):
# GCKR -> 293 | ZBTB38 -> 87 | C5orf67 -> 156 | HFE -> 72 | MLXIPL -> 149 | LOC157273 -> ? |
# ABO -> 301 | NRBF2 -> 50 | MYRF -> 119 | ATXN2 -> 236 | FTO -> 156 | APOC1 -> 207

qplot(y = res_median_postp, color = map_hm3_plus$chromosome) +
  ggrepel::geom_label_repel(
    with(df_for_label, aes(id, postp, color = chromosome, label = gene)),
    color = "purple", min.segment.length = 0) +
  scale_color_manual(values = rep_len(c("#000000", "#999999"), 22)) +
  scale_x_continuous(breaks = seq(0, 1.5e6, by = 200e3)) +
  ylim(0, NA) +
  theme(legend.position = "none") +
  labs(x = "Index of variant",
       y = "Median probability of being causal (across 169 phenotypes)")
# ggsave("figures/median_postp.png", width = 12, height = 7)

# verif it is not population structure
hm3_plus <- snp_attach("data/UKBB_hm3_plus.rds")
obj.pcadapt <- runonce::save_run({
  G <- hm3_plus$genotypes
  ind_train <- which(hm3_plus$fam$set == "train")
  U <- svd(select(hm3_plus$fam, PC1:PC4)[ind_train, ], nu = 4, nv = 0)$u
  snp_pcadapt(G, U, ind.row = ind_train, ncores = nb_cores())
}, file = "tmp-data/pcadapt.rds")
pcor(log(obj.pcadapt$score), res_median_postp, NULL)  # -5.5%

# verif it is not small ld
ld <- bigsnpr:::ld_scores_sfbm(readRDS("data/corr_hm3_plus.rds"), compact = TRUE, ncores = nb_cores())
pcor(log(ld), res_median_postp, NULL)  # -6.8%
pcor(ld, res_median_postp, NULL)  # 11.6%


#### RINT vs log-transformation ####

rint_pheno_files <- grep("rint_", grid$pheno_file, value = TRUE)

grid5 <- bind_cols(grid, res_grid) %>%
  filter(set == "hm3", use_mle, !jump_sign,
         pheno_file %in% c(rint_pheno_files, sub("rint_", "log_", rint_pheno_files))) %>%
  mutate(transfo = ifelse(grepl("rint_", pheno_file), "RIN", "log"),
         pheno = sub("^rint_|log_", "", basename(pheno_file))) %>%
  tidyr::pivot_longer(c(r2, r2_est, alpha_est, p_est, h2_est)) %>%
  tidyr::pivot_wider(c(pheno, name), names_from = transfo, values_from = value) %>%
  tidyr::unnest_wider("log", names_sep = "_") %>%
  tidyr::unnest_wider("RIN", names_sep = "_")

TO_SHOW <- grid5 %>%
  filter(name == "h2_est", (log_2 > RIN_3 | RIN_2 > log_3)) %>%
  pull(pheno)

grid5 %>%
  filter(name != "r2_est") %>%
  ggplot(aes(RIN_1, log_1, label = sub("\\.rds$", "", pheno))) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_errorbar(aes(ymin = log_2, ymax = log_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(xmin = RIN_2, xmax = RIN_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_abline(color = "red", linetype = 2) +
  geom_point(size = 1.2) +
  facet_wrap(~ name, scales = "free") +
  # ggrepel::geom_label_repel(color = "purple", min.segment.length = 0,
  #                           data = filter(grid5, pheno %in% TO_SHOW, name %in% c("r2", "h2_est"))) +
  labs(x = "Estimates for RIN-transformed phenotypes",
       y = "Estimates for log-transformed phenotypes")
# ggsave("figures/ukbb_compare_transfo.pdf", width = 6, height = 6)


raw_pheno_files <- grep("raw_", grid$pheno_file, value = TRUE)

grid6 <- bind_cols(grid, res_grid) %>%
  filter(set == "hm3", use_mle, !jump_sign,
         pheno_file %in% c(raw_pheno_files, sub("raw_", "log_", raw_pheno_files))) %>%
  mutate(transfo = ifelse(grepl("raw_", pheno_file), "none", "log"),
         pheno = sub("^raw_|log_", "", basename(pheno_file))) %>%
  tidyr::pivot_longer(c(r2, r2_est, alpha_est, p_est, h2_est)) %>%
  tidyr::pivot_wider(c(pheno, name), names_from = transfo, values_from = value) %>%
  tidyr::unnest_wider("log", names_sep = "_") %>%
  tidyr::unnest_wider("none", names_sep = "_")

grid6 %>%
  filter(name != "r2_est") %>%
  ggplot(aes(none_1, log_1, label = sub("\\.rds$", "", pheno))) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_errorbar(aes(ymin = log_2, ymax = log_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(xmin = none_2, xmax = none_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_abline(color = "red", linetype = 2) +
  geom_point(size = 1.2) +
  facet_wrap(~ name, scales = "free") +
  labs(x = "Estimates for raw phenotypes",
       y = "Estimates for log-transformed phenotypes")
# ggsave("figures/ukbb_compare_notransfo.pdf", width = 6, height = 6)


#### Option 'jump_sign' ####

grid7 <- bind_cols(grid, res_grid) %>%
  filter(set == "hm3", use_mle,
         !grepl("rint_", pheno_file), !grepl("raw_", pheno_file)) %>%
  tidyr::pivot_longer(c(r2, r2_est, alpha_est, p_est, h2_est)) %>%
  tidyr::pivot_wider(c(pheno_file, name), names_from = jump_sign, values_from = value) %>%
  tidyr::unnest_wider("TRUE", names_sep = "_") %>%
  tidyr::unnest_wider("FALSE", names_sep = "_")

anyNA(grid7) # FALSE
grid7 %>%
  filter(name != "r2_est") %>%
  ggplot(aes(TRUE_1, FALSE_1)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_errorbar(aes(ymin = FALSE_2, ymax = FALSE_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(xmin = TRUE_2, xmax = TRUE_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_abline(color = "red", linetype = 2) +
  geom_point(size = 1.2) +
  facet_wrap(~ name, scales = "free") +
  labs(x = "Estimates with normal sampling",
       y = "Estimates with extra robustness (no jump sign)")
# ggsave("figures/ukbb_compare_jumpsign.pdf", width = 6, height = 6)


#### Option 'use_MLE' ####

grid8 <- bind_cols(grid, res_grid) %>%
  filter(!jump_sign, !grepl("rint_", pheno_file), !grepl("raw_", pheno_file)) %>%
  rowwise() %>%
  mutate(across(r2, ~ list(`if`(n_keep == 0, rep(0, 3), .)))) %>%
  ungroup() %>%
  tidyr::pivot_longer(c(r2, r2_est, alpha_est, p_est, h2_est)) %>%
  filter(name == "r2") %>%
  tidyr::pivot_wider(c(pheno_file, set), names_from = use_mle, values_from = value) %>%
  tidyr::unnest_wider("TRUE", names_sep = "_") %>%
  tidyr::unnest_wider("FALSE", names_sep = "_")

anyNA(grid8) # FALSE
grid8 %>%
  ggplot(aes(TRUE_1, FALSE_1)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_errorbar(aes(ymin = FALSE_2, ymax = FALSE_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_errorbar(aes(xmin = TRUE_2, xmax = TRUE_3), color = "chartreuse3",
                position = position_dodge(width = 0.9), width = 0) +
  geom_abline(color = "red", linetype = 2) +
  geom_point(size = 1.2) +
  facet_wrap(~ set, scales = "free") +
  labs(x = "Estimates with new sampling (MLE)",
       y = "Estimates with previous sampling (no MLE)")
# ggsave("figures/ukbb_compare_mle.pdf", width = 6, height = 6)
