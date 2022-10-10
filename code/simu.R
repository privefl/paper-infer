library(bigsnpr)

#### Prepare LD matrix ####

dim(runonce::save_run({

  for (chr in seq(3, 22, by = 3)) {

    print(chr)

    corr0 <- readRDS(paste0("../misspec/ldref/LD_with_blocks_chr", chr, ".rds"))

    if (chr == 3) {
      corr <- as_SFBM(corr0, "data/corr_simu", compact = TRUE)
    } else {
      corr$add_columns(corr0, nrow(corr))
    }
  }

  corr
}, file = "data/corr_simu.rds"))
# 322805 x 322805


#### Run simulations ####

grid <- tidyr::expand_grid(
  p     = signif(seq_log(0.1, 5e-4, 4), 1),
  h2    = c(1, 3, 10, 30) / 100,
  alpha = c(0, -0.5, -1),
  N     = c(20, 200) * 1000,
  num   = 1,
  jump_sign = c(FALSE, TRUE)
)

library(future.batchtools)
NCORES <- 13
plan(workers = print(nrow(grid)) + 10, batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "60g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results_simu")

grid$res <- furrr::future_pmap(grid, function(p, h2, alpha, N, num, jump_sign) {

  print(params <- paste(c(p, h2, alpha, N, num, jump_sign), collapse = "_"))
  runonce::save_run(file = paste0("results_simu/res_", params, ".rds"), {

    corr <- readRDS("data/corr_simu.rds")

    simu <- snp_attach("data/ukbb4simu.rds")
    G <- simu$genotypes

    load("data/ukbb4simu_ind.RData")
    ind.gwas0 <- ind.gwas

    # simu quantitative pheno
    simu_pheno <- snp_simuPheno(G, h2 = h2, M = round(ncol(G) * p), alpha = alpha, ncores = NCORES)

    # GWAS to get sumstats
    ind.gwas <- sort(sample(ind.gwas0, N))
    gwas <- big_univLinReg(G, simu_pheno$pheno[ind.gwas], ind.train = ind.gwas, ncores = NCORES)
    df_beta <- dplyr::transmute(gwas, beta = estim, beta_se = std.err, n_eff = length(ind.gwas))

    # ldsc estimate
    (ldsc <- snp_ldsc2(corr, df_beta, blocks = 100, intercept = NULL, ncores = NCORES))

    # LDpred2-auto
    coef_shrink <- 0.95
    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                                   vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                                   burn_in = 500, num_iter = 500, report_step = 20,
                                   allow_jump_sign = jump_sign, shrink_corr = coef_shrink,
                                   ncores = NCORES)

    # derive other results
    (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
    (keep <- (range > (0.95 * quantile(range, 0.95))))

    hist(all_h2    <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est,    500)))
    hist(all_p     <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est,     500)))
    hist(all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500)))

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
    hist(all_r2)

    # r2 in test set
    beta <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    pred <- big_prodVec(G, beta, ind.row = ind.test, ncores = NCORES)
    cor <- pcor(pred, simu_pheno$pheno[ind.test], z = NULL)
    r2 <- sign(cor) * cor^2

    # local h2 (per-block)
    blocks <- subset(readRDS("../misspec/ldref/map.rds"), chr %% 3 == 0)$group_id
    list_ind <- split(seq_along(blocks), blocks)
    bsamp <- do.call("cbind", bsamp)
    all_local_h2 <- sapply(list_ind, function(ind) {
      corr_sub <- corr[ind, ind]
      bsamp_sub <- bsamp[ind, ]
      Rb <- coef_shrink * corr_sub %*% bsamp_sub + (1 - coef_shrink) * bsamp_sub
      Matrix::colSums(bsamp_sub * Rb)
    }) # 32 sec

    all_true_h2 <- sapply(list_ind, function(ind) {
      ind2 <- which(simu_pheno$set %in% ind)
      b <- simu_pheno$effects[ind2]
      ind3 <- simu_pheno$set[ind2]
      crossprod(b, corr[ind3, ind3] %*% b)[1]
    })

    # posterior probabilities of being causal (per-variant)
    postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp_est))

    list2 <- function(...) setNames(list(...), sapply(substitute(list(...))[-1], deparse))
    list2(simu_pheno, ldsc, postp, all_h2, all_p, all_alpha, all_r2, r2,
          all_local_h2, all_true_h2)
  })
})


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
  ggplot(aes(p, r2_est_1, color = n_keep, shape = ifelse(jump_sign, "Yes", "No"))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.65) +
  ylim(min(unlist(grid3$r2_est)), NA) +
  geom_point(alpha = 0) +
  facet_grid(h2 ~ N + alpha, scales = "free_y", labeller = label_both) +
  geom_point(size = 1.5, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = r2_est_2, ymax = r2_est_3),
                position = position_dodge(width = 0.6), width = 0.3) +
  geom_segment(aes(x = as.numeric(p) - ifelse(jump_sign, 0, 0.3),
                   xend = as.numeric(p) + ifelse(jump_sign, 0.3, 0),
                   y = r2_1, yend = r2_1), color = "chartreuse3", linetype = 1,
               position = position_dodge(width = 0.6)) +
  theme(legend.position = "top", legend.key.width = unit(30, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt"))) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred r2  (+ 95% CI of the estimate)",
       x = "Simulated polygenicity (p)",
       color = "# chains", shape = "Can jump sign?")
# ggsave("figures/est_r2_one.pdf", width = 10.5, height = 7)

grid3$h2_est <- purrr::map(grid3$all_h2, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))
grid3$h2_ldsc_est <- purrr::map(grid3$ldsc, ~ .[["h2"]] + c(0, -1.96, 1.96) * .[["h2_se"]])
grid3 %>%
  tidyr::pivot_longer(c(h2_est, h2_ldsc_est), values_to = "h2_est") %>%
  mutate(name = c("h2_est" = "LDpred2-auto", "h2_ldsc_est" = "LDSc reg")[name],
         n_keep = ifelse(name == "LDSc reg", NA, n_keep)) %>%
  tidyr::unnest_wider("h2_est", names_sep = "_") %>%
  ggplot(aes(as.factor(p), h2_est_1, shape = name, color = n_keep)) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "chartreuse3") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = h2), color = "grey40", linetype = 2) +
  geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = h2_est_2, ymax = h2_est_3),
                position = position_dodge(width = 0.5), width = 0.5) +
  facet_grid(N + h2 ~ jump_sign + alpha, scales = "free_y", labeller = label_both) +
  scale_shape_manual(values = c(1, 4)) +
  theme(legend.key.height = unit(20, "pt")) +
  labs(y = "Inferred h2  (+ 95% CI of the estimate)", x = "Simulated polygenicity (p)",
       shape = "Method", color = "# chains") +
  guides(shape = guide_legend(override.aes = list(color = c("black", "chartreuse3"))))
# ggsave("figures/est_h2_one.pdf", width = 12, height = 8)

grid3$p_est <- purrr::map(grid3$all_p, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))
grid3 %>%
  # filter(!jump_sign) %>%
  tidyr::unnest_wider("p_est", names_sep = "_") %>%
  ggplot(aes(as.factor(h2), p_est_1, color = n_keep, shape = ifelse(jump_sign, "Yes", "No"))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_hline(aes(yintercept = p), color = "grey40", linetype = 2) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = p_est_2, ymax = p_est_3),
                position = position_dodge(width = 0.5), width = 0.5) +
  facet_grid(p ~ N + alpha, scales = "free_y", labeller = label_both) +
  scale_y_log10() +
  theme(legend.position = "top", legend.key.width = unit(30, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt"))) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred polygenicity  (+ 95% CI of the estimate)",
       x = "Simulated heritability (h2)",
       color = "# chains", shape = "Can jump sign?")
# ggsave("figures/est_p_one.pdf", width = 11, height = 7)


grid3$alpha_est <- purrr::map(grid3$all_alpha, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))
grid3 %>%
  tidyr::unnest_wider("alpha_est", names_sep = "_") %>%
  ggplot(aes(as.factor(h2), alpha_est_1, color = n_keep, shape = ifelse(jump_sign, "Yes", "No"))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_hline(aes(yintercept = alpha), color = "grey40", linetype = 2) +
  geom_hline(yintercept = c(-1.5, 0.5), color = "black", linetype = 3) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = alpha_est_2, ymax = alpha_est_3),
                position = position_dodge(width = 0.5), width = 0.5) +
  facet_grid(alpha ~ N + p, labeller = label_both) +
  theme(legend.position = "top", legend.key.width = unit(30, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt"))) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred alpha  (+ 95% CI of the estimate)",
       x = "Simulated heritability (h2)",
       color = "# chains", shape = "Can jump sign?")
# ggsave("figures/est_alpha_one.pdf", width = 12, height = 7)


grid3$prop <- purrr::pmap(grid3[c("simu_pheno", "postp")], function(simu_pheno, postp) {

  ind <- split(seq_along(postp), cut(postp, seq_log(1e-5, 1, 11)))

  lapply(ind, function(.) {
    is_causal <- . %in% simu_pheno$set
    n <- length(is_causal)
    if (n < 2) return(NULL)
    unlist(prop.test(sum(is_causal), n)[c("estimate", "conf.int")], use.names = FALSE)
  }) %>%
    setNames(sapply(ind, function(.) mean(postp[.]))) %>%
    .[lengths(.) > 0]
})

grid3 %>%
  select(p:jump_sign, n_keep:prop) %>%
  filter(alpha == -0.5, !jump_sign) %>%
  tidyr::unnest_longer("prop", indices_to = "mid", transform = list(mid = as.double)) %>%
  tidyr::unnest_wider("prop", names_sep = "_") %>%
  ggplot(aes(mid, prop_1, color = n_keep)) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  scale_y_log10(breaks = 10^-(5:0), minor_breaks = seq_log(1e-5, 1, 11)) +
  scale_x_log10(breaks = c(1e-4, 0.1), minor_breaks = seq_log(1e-5, 1, 11)) +
  coord_equal() +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_point() +
  geom_abline(lty = 2, lwd = 1, color = "chartreuse2") +
  geom_errorbar(aes(ymin = prop_2, ymax = prop_3), width = 0.2) +
  facet_grid(h2 ~ N + p, labeller = label_both) +
  theme(legend.position = "top", legend.key.width = unit(30, "pt")) +
  guides(color = guide_colorbar(title.vjust = 0.9)) +
  labs(x = "Mean of posterior probabilities (in bin)",
       y = "Proportion of causal variants (in bin)", color = "# chains")
# ggsave("figures/postp_calib.pdf", width = 9.5, height = 7)

grid3 %>%
  select(p:jump_sign, n_keep, all_local_h2, all_true_h2) %>%
  filter(alpha == -0.5, !jump_sign) %>%
  mutate(across(all_local_h2, ~ lapply(., colMeans))) %>%
  tidyr::unnest(c(all_local_h2, all_true_h2)) %>%
  ggplot(aes(all_local_h2 / h2, all_true_h2 / h2, color = n_keep)) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  # coord_equal() +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_point() +
  geom_abline(lty = 2, lwd = 1, color = "chartreuse2") +
  facet_grid(N + p ~ h2, labeller = label_both, scales = "free") +
  theme(legend.key.height = unit(25, "pt")) +
  guides(color = guide_colorbar(title.vjust = 0.9)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.02)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.02)) +
  labs(x = "Estimated per-block heritability  (as percentage of h2)",
       y = "Simulated per-block heritability  (as percentage of h2)",
       color = "# chains")
# ggsave("figures/local_h2_calib.pdf", width = 11, height = 9)
