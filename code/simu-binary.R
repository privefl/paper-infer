library(bigsnpr)

#### Run simulations ####

grid <- tidyr::expand_grid(
  p     = signif(seq_log(0.1, 5e-4, 4), 1),
  h2    = c(1, 3, 10, 30) / 100,
  alpha = c(0, -0.5, -1),
  K     = c(0.02, 0.2),
  num   = 1
)

NCORES <- 13
library(future.batchtools)
plan(workers = print(nrow(grid)) + 10, batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "60g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))
# future::plan("multisession", workers = 14)

bigassertr::assert_dir("results_simu_binary")

grid$res <- furrr::future_pmap(grid, function(p, h2, alpha, K, num) {

  # p <- 0.02
  # h2 <- 0.03
  # alpha <- -0.5
  # K <- 0.2
  # num <- 1
  print(params <- paste(c(p, h2, alpha, K, num), collapse = "_"))

  corr <- readRDS("data/corr_simu.rds")
  simu <- snp_attach("data/ukbb4simu.rds")
  G <- simu$genotypes
  load("data/ukbb4simu_ind.RData")

  # simu binary pheno
  simu_pheno <- runonce::save_run(file = paste0("results_simu_binary/simu_", params, ".rds"), {
    snp_simuPheno(G, h2 = h2, M = round(ncol(G) * p), K = K, alpha = alpha, ncores = NCORES)
  })

  # GWAS to get sumstats (linear reg)
  df_beta <- runonce::save_run(file = paste0("results_simu_binary/gwas_lin_", params, ".rds"), {
    gwas <- big_univLinReg(G, simu_pheno$pheno[ind.gwas], ind.train = ind.gwas, ncores = NCORES)
    dplyr::transmute(gwas, beta = estim, beta_se = std.err, n_eff = length(ind.gwas))
  })

  # ldsc estimate
  ldsc <- runonce::save_run(file = paste0("results_simu_binary/ldsc_lin_", params, ".rds"), {
    snp_ldsc2(corr, df_beta, blocks = 100, intercept = NULL, ncores = NCORES)
  })

  # Run methods
  source('code/simu-methods.R', local = TRUE)
  res_lin <- list(
    LDpred2_noMLE = runonce::save_run(
      run_ldpred2(jump_sign = FALSE, use_mle = FALSE, coef_shrink = 0.95),
      file = paste0("results_simu_binary/LDpred2_noMLE_lin_", params, ".rds")),
    LDpred2_MLE = runonce::save_run(
      run_ldpred2(jump_sign = FALSE, use_mle = TRUE,  coef_shrink = 0.95),
      file = paste0("results_simu_binary/LDpred2_nojump_lin_", params, ".rds")),
    SBayesS = runonce::save_run(
      run_sbayess(),
      file = paste0("results_simu_binary/SBayesS_lin_", params, ".rds")),
    LDSc = list(all_h2 = ldsc)
  )


  # GWAS to get sumstats (log reg)
  df_beta <- runonce::save_run(file = paste0("results_simu_binary/gwas_log_", params, ".rds"), {
    y <- simu_pheno$pheno[ind.gwas]
    gwas <- big_univLogReg(G, y, ind.train = ind.gwas, ncores = NCORES)
    Neff <- 4 / (1 / sum(y == 1) + 1 / sum(y == 0))
    dplyr::transmute(gwas, beta = estim, beta_se = std.err, n_eff = Neff)
  })

  # ldsc estimate
  ldsc <- runonce::save_run(file = paste0("results_simu_binary/ldsc_log_", params, ".rds"), {
    snp_ldsc2(corr, df_beta, blocks = 100, intercept = NULL, ncores = NCORES)
  })

  # Run methods
  res_log <- list(
    LDpred2_noMLE = runonce::save_run(
      run_ldpred2(jump_sign = FALSE, use_mle = FALSE, coef_shrink = 0.95),
      file = paste0("results_simu_binary/LDpred2_noMLE_log_", params, ".rds")),
    LDpred2_MLE = runonce::save_run(
      run_ldpred2(jump_sign = FALSE, use_mle = TRUE,  coef_shrink = 0.95),
      file = paste0("results_simu_binary/LDpred2_nojump_log_", params, ".rds")),
    SBayesS = runonce::save_run(
      run_sbayess(),
      file = paste0("results_simu_binary/SBayesS_log_", params, ".rds")),
    LDSc = list(all_h2 = ldsc)
  )


  # return both
  tibble::tibble(logistic = c(FALSE, TRUE), res = list(res_lin, res_log))
})


library(dplyr)

grid %>%
  tidyr::unnest_wider(res) %>%
  tidyr::unnest(cols = c(logistic, res)) %>%
  tidyr::unnest_wider(res) %>%
  tidyr::pivot_longer(-(1:6), names_to = "Method") %>%
  tidyr::unnest_wider("value") %>%
  # head(20) %>%
  rowwise() %>%
  mutate(
    r2 = list(unname(r2)),
    n_keep = `if`(grepl("LDpred2", Method), `if`(is.null(all_h2), 0, NCOL(all_h2)), NA),
    r2_est = list(`if`(Method == "LDSc", NULL,
                       unname(quantile(all_r2, probs = c(0.5, 0.025, 0.975))))),
    h2_est = list(
      `if`(Method == "LDSc", all_h2[["h2"]] + c(0, -1.96, 1.96) * all_h2[["h2_se"]],
           unname(quantile(all_h2, probs = c(0.5, 0.025, 0.975))))
    ),
    p_est = list(`if`(Method == "LDSc", NULL,
                      unname(quantile(all_p, probs = c(0.5, 0.025, 0.975))))),
    alpha_est = list(`if`(Method %in% c("LDSc", "LDpred2_noMLE"), NULL,
                          unname(quantile(all_alpha, probs = c(0.5, 0.025, 0.975))))),
  ) %>%
  relocate(n_keep, .after = use_mle) %>%
  ungroup() %>%
  mutate(Method = factor(Method, levels = Method[1:4],
                         labels = sub("LDpred2_", "LDpred2-auto_", Method[1:4]))) %>%
  select(-starts_with("all"), -postp) %>%
  tidyr::unnest_wider(where(is.list), names_sep = "_") %>%
  mutate(across(r2_1:alpha_est_3, ~ signif(., digits = 4))) %>%
  print(n = 20) %>%
  bigreadr::fwrite2("results_infer_simu_binary.csv")


#### Just one simulation with CIs of the estimates ####

library(ggplot2)

grid2 <- grid %>%
  filter(num == 1) %>%
  tidyr::unnest_wider(res) %>%
  tidyr::unnest(cols = c(logistic, res)) %>%
  tidyr::unnest_wider(res) %>%
  tidyr::pivot_longer(-(1:6), names_to = "Method") %>%
  tidyr::unnest_wider("value") %>%
  mutate(n_keep = purrr::map_dbl(all_h2, ~ `if`(is.null(.), 0, NCOL(.))),
         Method = factor(Method, levels = Method[1:4],
                         labels = sub("LDpred2_", "LDpred2-auto_", Method[1:4])),
         p = as.factor(p)) %>%
  relocate(n_keep, .after = Method) %>%
  print(n = 100)

vec_to_liab <- function(x, K, logistic) {
  x * purrr::map2_dbl(K, logistic, ~ coef_to_liab(.x, `if`(.y, 0.5, .x)))
}

grid2_r2 <- grid2 %>%
  filter(Method != "LDSc") %>%
  mutate(r2 = lapply(r2, unname),
         r2_est = purrr::map(all_r2, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("r2",     names_sep = "_") %>%
  tidyr::unnest_wider("r2_est", names_sep = "_") %>%
  # transformation to liability scale
  mutate(across(matches("r2_[1-3]"),     ~ vec_to_liab(., K, FALSE))) %>%
  mutate(across(matches("r2_est_[1-3]"), ~ vec_to_liab(., K, logistic)))

grid2_r2 %>%
  ggplot(aes(Method, r2_est_1, color = ifelse(Method == "SBayesS", NA, n_keep),
             shape = ifelse(logistic, "Yes", "No"))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  facet_grid(h2 + p ~ K + alpha, scales = "free_y", labeller = label_both) +
  geom_point(size = 1.5, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = r2_est_2, ymax = r2_est_3),
                position = position_dodge(width = 0.8), width = 0.25) +
  geom_errorbarh(aes(y = r2_1, xmin = as.numeric(Method) - ifelse(logistic, 0, 0.4),
                     xmax = as.numeric(Method) + ifelse(logistic, 0.4, 0), height = r2_3 - r2_2),
                 color = "chartreuse3") +
  theme(legend.position = "top", legend.key.width = unit(40, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt")),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred r2  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains", shape = "Use Neff?")
# ggsave("figures/est_binary_r2_one.pdf", width = 10, height = 13)


grid2 %>%
  mutate(h2_est = purrr::map2(all_h2, Method, ~ {
    `if`(.y == "LDSc", .x[["h2"]] + c(0, -1.96, 1.96) * .x[["h2_se"]],
         unname(quantile(.x, probs = c(0.5, 0.025, 0.975))))
  })) %>%
  tidyr::unnest_wider("h2_est", names_sep = "_") %>%
  # transformation to the liability scale
  mutate(across(starts_with("h2_est_"), ~ vec_to_liab(., K, logistic))) %>%
  ggplot(aes(Method, h2_est_1, shape = ifelse(logistic, "Yes", "No"),
         color = ifelse(grepl("LDpred2", Method), n_keep, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = h2), color = "chartreuse3", linetype = 2) +
  geom_point(size = 1.5, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = h2_est_2, ymax = h2_est_3),
                position = position_dodge(width = 0.8), width = 0.5) +
  facet_grid(h2 + p ~ K + alpha, scales = "free_y", labeller = label_both) +
  theme(legend.position = "top", legend.key.width = unit(40, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt")),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred h2  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains", shape = "Use Neff?")
# ggsave("figures/est_binary_h2_one.pdf", width = 10, height = 13)


grid2 %>%
  filter(Method != "LDSc") %>%
  mutate(p_est = purrr::map(all_p, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("p_est", names_sep = "_") %>%
  ggplot(aes(Method, p_est_1, shape = ifelse(logistic, "Yes", "No"),
             color = ifelse(grepl("LDpred2", Method), n_keep, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = as.double(as.character(p))),
             color = "chartreuse3", linetype = 2) +
  geom_point(size = 1.5, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = p_est_2, ymax = p_est_3),
                position = position_dodge(width = 0.8), width = 0.5) +
  scale_y_log10() +
  facet_grid(h2 + p ~ K + alpha, scales = "free_y", labeller = label_both) +
  theme(legend.position = "top", legend.key.width = unit(40, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt")),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred polygenicity  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains", shape = "Use Neff?")
# ggsave("figures/est_binary_p_one.pdf", width = 10, height = 13)


grid2 %>%
  filter(!Method %in% c("LDSc", "LDpred2-auto_noMLE")) %>%
  mutate(alpha_est = purrr::map(all_alpha, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("alpha_est", names_sep = "_") %>%
  ggplot(aes(Method, alpha_est_1, shape = ifelse(logistic, "Yes", "No"),
             color = ifelse(grepl("LDpred2", Method), n_keep, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = alpha), color = "chartreuse3", linetype = 2) +
  geom_hline(yintercept = c(-1.5, 0.5), color = "black", linetype = 3) +
  geom_point(size = 1.5, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = alpha_est_2, ymax = alpha_est_3),
                position = position_dodge(width = 0.8), width = 0.5) +
  facet_grid(h2 + alpha ~ K + p, scales = "free_y", labeller = label_both) +
  theme(legend.position = "top", legend.key.width = unit(40, "pt"),
        legend.margin = margin(0, 1, 0, 1, unit = "cm"),
        legend.text = element_text(margin = margin(l = -10, r = 5, unit = "pt")),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  guides(shape = guide_legend(override.aes = list(size = 2.5), label.hjust = 1),
         color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred alpha  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains", shape = "Use Neff?")
# ggsave("figures/est_binary_alpha_one.pdf", width = 9, height = 13)
