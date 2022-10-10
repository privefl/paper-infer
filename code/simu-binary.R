library(bigsnpr)

#### Run simulations ####

grid <- tidyr::expand_grid(
  p     = signif(seq_log(0.1, 5e-4, 4), 1),
  h2    = c(1, 3, 10, 30) / 100,
  alpha = c(0, -0.5, -1),
  # N     = c(20, 200) * 1000,
  K     = c(0.02, 0.2),
  num   = 1,
  logistic = c(TRUE, FALSE)
)

library(future.batchtools)
NCORES <- 13
plan(workers = print(nrow(grid)) + 10, batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "60g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results_simu")

grid$res <- furrr::future_pmap(grid, function(p, h2, alpha, K, num, logistic) {

  print(params <- paste(c(p, h2, alpha, K, num, logistic), collapse = "_"))
  runonce::save_run(file = paste0("results_simu/res_binary_", params, ".rds"), {

    corr <- readRDS("data/corr_simu.rds")

    simu <- snp_attach("data/ukbb4simu.rds")
    G <- simu$genotypes

    # simu quantitative pheno
    simu_pheno <- snp_simuPheno(G, h2 = h2, M = round(ncol(G) * p), K = K,
                                alpha = alpha, ncores = NCORES)

    # GWAS to get sumstats
    load("data/ukbb4simu_ind.RData")
    y <- simu_pheno$pheno[ind.gwas]
    if (logistic) {
      gwas <- big_univLogReg(G, y, ind.train = ind.gwas, ncores = NCORES)
      Neff <- 4 / (1 / sum(y == 1) + 1 / sum(y == 0))
    } else {
      gwas <- big_univLinReg(G, y, ind.train = ind.gwas, ncores = NCORES)
      Neff <- length(ind.gwas)
    }
    df_beta <- dplyr::transmute(gwas, beta = estim, beta_se = std.err, n_eff = Neff)

    # ldsc estimate
    (ldsc <- snp_ldsc2(corr, df_beta, blocks = 100, intercept = NULL, ncores = NCORES))

    # LDpred2-auto
    coef_shrink <- 0.95
    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                                   vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                                   burn_in = 500, num_iter = 500, report_step = 20,
                                   allow_jump_sign = FALSE, shrink_corr = coef_shrink,
                                   ncores = NCORES)

    # derive other results
    (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
    (keep <- (range > (0.95 * quantile(range, 0.95))))

    hist(all_h2    <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est,    500)))
    hist(all_p     <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est,     500)))
    hist(all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500)))

    bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
    all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
      # print(ic)
      b1 <- bsamp[[ic]]
      Rb1 <- apply(b1, 2, function(x)
        coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
      b2 <- do.call("cbind", bsamp[-ic])
      b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
    }))
    hist(all_r2)

    beta <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    pred <- big_prodVec(G, beta, ind.row = ind.test, ncores = NCORES)
    cor <- pcor(pred, simu_pheno$pheno[ind.test], z = NULL)
    r2 <- sign(cor) * cor^2  # TODO: use bootstrap again

    postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp_est))

    list2 <- function(...) setNames(list(...), sapply(substitute(list(...))[-1], deparse))
    list2(simu_pheno, ldsc, postp, all_h2, all_p, all_alpha, all_r2, r2, Neff)
  })
})


vec_to_liab <- function(x, K, logistic) {
  x * purrr::map2_dbl(K, logistic, ~ coef_to_liab(.x, `if`(.y, 0.5, .x)))
}

#### Just one simulation with CIs of the estimates ####

library(dplyr)
library(ggplot2)

grid3 <- grid %>%
  filter(num == 1) %>%
  tidyr::unnest_wider("res") %>%
  mutate(n_keep = purrr::map_dbl(all_h2, ncol))

grid3$r2_est <- purrr::map(grid3$all_r2, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))
grid3 %>%
  mutate(p = as.factor(p)) %>%
  tidyr::unnest_wider("r2", names_sep = "_") %>%
  tidyr::unnest_wider("r2_est", names_sep = "_") %>%

  mutate(r2_1 = vec_to_liab(r2_1, K, FALSE)) %>%
  mutate(across(starts_with("r2_est_"), ~ vec_to_liab(., K, logistic))) %>%

  ggplot(aes(p, r2_est_1, color = n_keep, shape = ifelse(logistic, "Yes", "No"))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.65) +
  ylim(min(unlist(grid3$r2_est)), NA) +
  geom_point(alpha = 0) +
  facet_grid(h2 ~ K + alpha, scales = "free_y", labeller = label_both) +
  geom_point(size = 1.5, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = r2_est_2, ymax = r2_est_3),
                position = position_dodge(width = 0.7), width = 0.25) +
  geom_segment(aes(x = as.numeric(p) - ifelse(logistic, 0, 0.35),
                   xend = as.numeric(p) + ifelse(logistic, 0.35, 0),
                   y = r2_1, yend = r2_1), color = "chartreuse3", linetype = 1,
               position = position_dodge(width = 0.7)) +
  theme(legend.position = "top", legend.key.width = unit(30, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt"))) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred r2  (+ 95% CI of the estimate)",
       x = "Simulated polygenicity (p)",
       color = "# chains", shape = "Use Neff??")
# ggsave("figures/est_binary_r2_one.pdf", width = 10.5, height = 7)

grid3$h2_est <- purrr::map(grid3$all_h2, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))
grid3$h2_ldsc_est <- purrr::map(grid3$ldsc, ~ .[["h2"]] + c(0, -1.96, 1.96) * .[["h2_se"]])
grid3 %>%
  tidyr::pivot_longer(c(h2_est, h2_ldsc_est), values_to = "h2_est") %>%
  mutate(name = c("h2_est" = "LDpred2-auto", "h2_ldsc_est" = "LDSc reg")[name],
         n_keep = ifelse(name == "LDSc reg", NA, n_keep), use_Neff = logistic) %>%
  tidyr::unnest_wider("h2_est", names_sep = "_") %>%

  # transformation to the liability scale
  mutate(across(starts_with("h2_est_"), ~ vec_to_liab(., K, logistic))) %>%

  filter(h2_est_1 > -1) %>%
  ggplot(aes(as.factor(p), h2_est_1, shape = name, color = n_keep)) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "chartreuse3") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = h2), color = "grey40", linetype = 2) +
  geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = h2_est_2, ymax = h2_est_3),
                position = position_dodge(width = 0.5), width = 0.5) +
  facet_grid(K + h2 ~ use_Neff + alpha, scales = "free_y", labeller = label_both) +
  scale_shape_manual(values = c(1, 4)) +
  labs(y = "Inferred h2  (+ 95% CI of the estimate)", x = "Simulated polygenicity (p)",
       shape = "Method", color = "# chains") +
  theme(legend.key.height = unit(20, "pt")) +
  guides(shape = guide_legend(override.aes = list(color = c("black", "chartreuse3"))))
# ggsave("figures/est_binary_h2_one.pdf", width = 12, height = 8)

grid3$p_est <- purrr::map(grid3$all_p, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))
grid3 %>%
  # filter(!logistic) %>%
  tidyr::unnest_wider("p_est", names_sep = "_") %>%
  ggplot(aes(as.factor(h2), p_est_1, color = n_keep, shape = ifelse(logistic, "Yes", "No"))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_hline(aes(yintercept = p), color = "grey40", linetype = 2) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = p_est_2, ymax = p_est_3),
                position = position_dodge(width = 0.5), width = 0.5) +
  facet_grid(p ~ K + alpha, scales = "free_y", labeller = label_both) +
  scale_y_log10() +
  theme(legend.position = "top", legend.key.width = unit(30, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt"))) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred polygenicity  (+ 95% CI of the estimate)",
       x = "Simulated heritability (h2)",
       color = "# chains", shape = "Use Neff??")
# ggsave("figures/est_binary_p_one.pdf", width = 11, height = 7)


grid3$alpha_est <- purrr::map(grid3$all_alpha, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))
grid3 %>%
  tidyr::unnest_wider("alpha_est", names_sep = "_") %>%
  ggplot(aes(as.factor(h2), alpha_est_1, color = n_keep, shape = ifelse(logistic, "Yes", "No"))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_hline(aes(yintercept = alpha), color = "grey40", linetype = 2) +
  geom_hline(yintercept = c(-1.5, 0.5), color = "black", linetype = 3) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = alpha_est_2, ymax = alpha_est_3),
                position = position_dodge(width = 0.5), width = 0.5) +
  facet_grid(alpha ~ K + p, labeller = label_both) +
  theme(legend.position = "top", legend.key.width = unit(30, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt"))) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred alpha  (+ 95% CI of the estimate)",
       x = "Simulated heritability (h2)",
       color = "# chains", shape = "Use Neff??")
# ggsave("figures/est_binary_alpha_one.pdf", width = 11, height = 7)
