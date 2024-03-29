library(bigsnpr)
library(dplyr)

#### Prepare LD matrix ####

dim(corr <- runonce::save_run({

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

# Same data, but format for GCTB
binfile <- "tmp-data/corr_simu.ldm.sparse.bin"
infofile <- sub("\\.bin$", ".info", binfile)
runonce::skip_run_if(files = c(binfile, infofile), {

  con <- file(binfile, open = "wb")
  system.time(
    info <- sapply(1:ncol(corr), function(j) {
      if (j %% 1000 == 0) print(j)
      one_col <- corr[, j]
      I <- one_col@i
      X <- one_col@x
      L <- length(I)
      writeBin(as.integer(I), con, size = 4)
      writeBin(as.double(X),  con, size = 4)
      c(I[1], I[L], L, sum(X))
    })
  ) # 35 min (for 322,805 variants)
  close(con)

  map <- subset(readRDS("../misspec/ldref/map.rds"), chr %in% seq(3, 22, by = 3))
  stopifnot(nrow(map) == ncol(corr))
  info2 <- dplyr::transmute(map, Chrom = chr, ID = rsid, GenPos = 0, PhysPos = pos,
                            A1 = a1, A2 = a0, A2Freq = af_UKBB,
                            Index = 0:(ncol(info) - 1),
                            WindStart = info[1, ],
                            WindEnd = info[2, ],
                            WindSize = info[3, ],
                            WindWidth = -1,
                            N = 360e3,
                            SamplVar = -1,
                            LDsum = info[4, ])

  bigreadr::fwrite2(info2, infofile, sep = " ", scipen = 50)
})


#### Run simulations ####

grid <- tidyr::expand_grid(
  p     = signif(seq_log(0.1, 5e-4, 4), 1),
  h2    = c(1, 3, 10, 30) / 100,
  alpha = c(0, -0.5, -1),
  N     = c(20, 200) * 1000,
  num   = 1:10
) %>%
  filter(num == 1 | alpha == -0.5)

NCORES <- 13
library(future.batchtools)
plan(workers = print(nrow(grid)) + 10, batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "60g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))
# future::plan("multisession", workers = 14)

bigassertr::assert_dir("results_simu")

grid$res <- furrr::future_pmap(grid[1:5], function(p, h2, alpha, N, num) {

  # p <- 0.1
  # h2 <- 0.01
  # alpha <- -0.5
  # N <- 20e3
  # num <- 1
  print(params <- paste(c(p, h2, alpha, N, num), collapse = "_"))

  simu <- snp_attach("data/ukbb4simu.rds")
  G <- simu$genotypes

  load("data/ukbb4simu_ind.RData")
  ind.gwas0 <- ind.gwas

  # simu quantitative pheno
  simu_pheno <- runonce::save_run(file = paste0("results_simu/simu_", params, ".rds"), {
    snp_simuPheno(G, h2 = h2, M = round(ncol(G) * p), alpha = alpha, ncores = NCORES)
  })

  # GWAS to get sumstats
  df_beta <- runonce::save_run(file = paste0("results_simu/gwas_", params, ".rds"), {
    ind.gwas <- sort(sample(ind.gwas0, N))
    gwas <- big_univLinReg(G, simu_pheno$pheno[ind.gwas], ind.train = ind.gwas, ncores = NCORES)
    dplyr::transmute(gwas, beta = estim, beta_se = std.err, n_eff = length(ind.gwas))
  })

  # ldsc estimate
  ldsc <- runonce::save_run(file = paste0("results_simu/ldsc_", params, ".rds"), {
    snp_ldsc2(corr, df_beta, blocks = 100, intercept = NULL, ncores = NCORES)
  })


  # Run methods
  source('code/simu-methods.R', local = TRUE)
  list(
    `LDpred2-auto_altfilter` = runonce::save_run(
      run_ldpred2_new(jump_sign = FALSE, use_mle = TRUE, coef_shrink = 0.95),
      file = paste0("results_simu/LDpred2_newfilter_", params, ".rds")),
    `LDpred2-auto_noMLE` = runonce::save_run(
      run_ldpred2(jump_sign = FALSE, use_mle = FALSE, coef_shrink = 0.95),
      file = paste0("results_simu/LDpred2_noMLE_", params, ".rds")),
    `LDpred2-auto_nojump` = runonce::save_run(
      run_ldpred2(jump_sign = FALSE, use_mle = TRUE,  coef_shrink = 0.95),
      file = paste0("results_simu/LDpred2_nojump_", params, ".rds")),
    `LDpred2-auto_jump` = runonce::save_run(
      run_ldpred2(jump_sign = TRUE,  use_mle = TRUE,  coef_shrink = 0.95),
      file = paste0("results_simu/LDpred2_", params, ".rds")),
    SBayesS = runonce::save_run(
      run_sbayess(),
      file = paste0("results_simu/SBayesS_", params, ".rds")),
    LDSC = list(all_h2 = ldsc)
  )
})


library(dplyr)

grid %>%
  # filter(num == 1) %>%
  tidyr::unnest_wider("res") %>%
  tidyr::pivot_longer(-(1:5), names_to = "Method") %>%
  tidyr::unnest_wider("value") %>%
  # head(20) %>%
  rowwise() %>%
  mutate(
    r2 = list(unname(r2)),
    n_keep_infer = `if`(grepl("LDpred2", Method), `if`(is.null(all_h2), 0, NCOL(all_h2)), NA),
    n_keep_pred  = `if`(grepl("LDpred2", Method), `if`(is.null(all_r2), 0, NCOL(all_r2) / 500 * 20), NA),
    r2_est = list(`if`(Method == "LDSC", NULL,
                       unname(quantile(all_r2, probs = c(0.5, 0.025, 0.975))))),
    h2_est = list(
      `if`(Method == "LDSC", all_h2[["h2"]] + c(0, -1.96, 1.96) * all_h2[["h2_se"]],
           unname(quantile(all_h2, probs = c(0.5, 0.025, 0.975))))
    ),
    p_est = list(`if`(Method == "LDSC", NULL,
                      unname(quantile(all_p, probs = c(0.5, 0.025, 0.975))))),
    alpha_est = list(`if`(Method %in% c("LDSC", "LDpred2-auto_noMLE"), NULL,
                          unname(quantile(all_alpha, probs = c(0.5, 0.025, 0.975))))),
  ) %>%
  relocate(n_keep_infer, n_keep_pred, .after = use_mle) %>%
  ungroup() %>%
  select(-starts_with("all"), -postp) %>%
  tidyr::unnest_wider(where(is.list), names_sep = "_") %>%
  mutate(across(r2_1:alpha_est_3, ~ signif(., digits = 4))) %>%
  print(n = 20) %>%
  bigreadr::fwrite2("results_infer_simu.csv")


#### Just one simulation with CIs of the estimates ####

library(ggplot2)
library(latex2exp)

grid2 <- grid %>%
  filter(num == 1) %>%
  tidyr::unnest_wider("res") %>%
  tidyr::pivot_longer(-(1:5), names_to = "Method") %>%
  tidyr::unnest_wider("value") %>%
  mutate(r2 = lapply(r2, unname),
         n_keep_infer = purrr::map_dbl(all_h2, ~ `if`(is.null(.), 0, NCOL(.))),
         n_keep_pred  = purrr::map_dbl(all_r2, ~ `if`(is.null(.), 0, NCOL(.) / 500 * 20)),
         Method = factor(Method, levels = Method[1:6]),
         p = as.factor(p)) %>%
  relocate(n_keep_infer, n_keep_pred, .after = Method) %>%
  print(n = 20)


ggplot(filter(grid2, grepl("LDpred2", Method))) +
  geom_histogram(aes(n_keep_infer, fill = Method), color = "black", linewidth = 0.1,
                 breaks = seq(0, 50, by = 2)) +
  bigstatsr::theme_bigstatsr(0.8) +
  facet_grid(p ~ N + h2, labeller = label_both) +
  labs(x = "Number of chains kepts") +
  scale_y_continuous(breaks = 0:10 * 2) +
  theme(legend.position = "top")


grid2_r2 <- grid2 %>%
  filter(Method != "LDSC") %>%
  mutate(r2_est = purrr::map(all_r2, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("r2",     names_sep = "_") %>%
  tidyr::unnest_wider("r2_est", names_sep = "_")

grid2_r2 %>%
  filter(N == 200e3) %>%
  ggplot(aes(Method, r2_est_1, color = ifelse(Method == "SBayesS", NA, n_keep_pred))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_point(alpha = 0) +
  facet_grid(h2 + p ~ N + alpha, scales = "free_y", labeller = label_both) +
  geom_errorbarh(aes(y = r2_1, xmin = as.numeric(Method) - 0.42,
                     xmax = as.numeric(Method) + 0.42, height = r2_3 - r2_2),
                 color = "chartreuse3") +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = r2_est_2, ymax = r2_est_3), width = 0.3) +
  theme(legend.key.height = unit(50, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  labs(y = TeX("Inferred $r^2$  (+ 95% CI of the estimate)"),
       x = "Method", color = "# chains")
# ggsave("figures/est_r2_one_20K.pdf",  width = 9, height = 13)
# ggsave("figures/est_r2_one_200K.pdf", width = 9, height = 13)


grid2 %>%
  filter(N == 20e3, alpha == -0.5) %>%
  mutate(h2_est = purrr::map2(all_h2, Method, ~ {
    `if`(.y == "LDSC", .x[["h2"]] + c(0, -1.96, 1.96) * .x[["h2_se"]],
         unname(quantile(.x, probs = c(0.5, 0.025, 0.975))))
  })) %>%
  tidyr::unnest_wider("h2_est", names_sep = "_") %>%
  ggplot(aes(Method, h2_est_1, color = ifelse(grepl("LDpred2", Method), n_keep_infer, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = h2), color = "chartreuse3", linetype = 2) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = h2_est_2, ymax = h2_est_3), width = 0.5) +
  facet_grid(h2 ~ N + alpha + p, scales = "free_y", labeller = label_both) +
  theme(legend.key.height = unit(40, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.6), "cm")) +
  labs(y = TeX("Inferred $h^2$  (+ 95% CI of the estimate)"),
       x = "Method", color = "# chains")
# ggsave("figures/est_h2_one_main.pdf", width = 10, height = 6)

grid2 %>%
  filter(N == 20e3) %>%
  mutate(h2_est = purrr::map2(all_h2, Method, ~ {
    `if`(.y == "LDSC", .x[["h2"]] + c(0, -1.96, 1.96) * .x[["h2_se"]],
         unname(quantile(.x, probs = c(0.5, 0.025, 0.975))))
  })) %>%
  tidyr::unnest_wider("h2_est", names_sep = "_") %>%
  ggplot(aes(Method, h2_est_1, color = ifelse(grepl("LDpred2", Method), n_keep_infer, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = h2), color = "chartreuse3", linetype = 2) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = h2_est_2, ymax = h2_est_3), width = 0.5) +
  facet_grid(h2 + p ~ N + alpha, scales = "free_y", labeller = label_both) +
  theme(legend.key.height = unit(50, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  labs(y = TeX("Inferred $h^2$  (+ 95% CI of the estimate)"),
       x = "Method", color = "# chains")
# ggsave("figures/est_h2_one_20K.pdf",  width = 10, height = 13)
# ggsave("figures/est_h2_one_200K.pdf", width = 10, height = 13)


grid2 %>%
  filter(N == 20e3, Method != "LDSC", alpha == -0.5) %>%
  mutate(p_est = purrr::map(all_p, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("p_est", names_sep = "_") %>%
  ggplot(aes(Method, p_est_1, color = ifelse(grepl("LDpred2", Method), n_keep_infer, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = as.double(as.character(p))),
             color = "chartreuse3", linetype = 2) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = p_est_2, ymax = p_est_3), width = 0.5) +
  facet_grid(p ~ N + alpha + h2, scales = "free_y", labeller = label_both) +
  scale_y_log10() +
  theme(legend.key.height = unit(40, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.6), "cm")) +
  labs(y = "Inferred polygenicity  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains")
# ggsave("figures/est_p_one_main.pdf", width = 10, height = 6)

grid2 %>%
  filter(N == 200e3, Method != "LDSC") %>%
  mutate(p_est = purrr::map(all_p, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("p_est", names_sep = "_") %>%
  ggplot(aes(Method, p_est_1, color = ifelse(grepl("LDpred2", Method), n_keep_infer, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = as.double(as.character(p))),
             color = "chartreuse3", linetype = 2) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = p_est_2, ymax = p_est_3), width = 0.5) +
  facet_grid(h2 + p ~ N + alpha, scales = "free_y", labeller = label_both) +
  scale_y_log10() +
  theme(legend.key.height = unit(50, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  labs(y = "Inferred polygenicity  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains")
# ggsave("figures/est_p_one_20K.pdf",  width = 9, height = 13)
# ggsave("figures/est_p_one_200K.pdf", width = 9, height = 13)


grid2 %>%
  filter(N == 20e3, !Method %in% c("LDSC", "LDpred2-auto_noMLE"), alpha == -0.5) %>%
  mutate(alpha_est = purrr::map(all_alpha, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("alpha_est", names_sep = "_") %>%
  ggplot(aes(Method, alpha_est_1, color = ifelse(grepl("LDpred2", Method), n_keep_infer, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = alpha), color = "chartreuse3", linetype = 2) +
  geom_hline(yintercept = c(-1.5, 0.5), color = "black", linetype = 3) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = alpha_est_2, ymax = alpha_est_3), width = 0.5) +
  facet_grid(h2 ~ N + alpha + p, scales = "free_y", labeller = label_both) +
  theme(legend.key.height = unit(40, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.7), "cm")) +
  labs(y = "Inferred alpha  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains")
# ggsave("figures/est_alpha_one_main.pdf", width = 9, height = 6)

grid2 %>%
  filter(N == 200e3, !Method %in% c("LDSC", "LDpred2-auto_noMLE")) %>%
  mutate(alpha_est = purrr::map(all_alpha, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("alpha_est", names_sep = "_") %>%
  ggplot(aes(Method, alpha_est_1, color = ifelse(grepl("LDpred2", Method), n_keep_infer, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = alpha), color = "chartreuse3", linetype = 2) +
  geom_hline(yintercept = c(-1.5, 0.5), color = "black", linetype = 3) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = alpha_est_2, ymax = alpha_est_3), width = 0.5) +
  facet_grid(h2 + alpha ~ N + p, scales = "free_y", labeller = label_both) +
  theme(legend.key.height = unit(50, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  labs(y = "Inferred alpha  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains")
# ggsave("figures/est_alpha_one_20K.pdf",  width = 9, height = 12)
# ggsave("figures/est_alpha_one_200K.pdf", width = 9, height = 12)


grid$res2 <- purrr::pmap(grid[1:5], function(p, h2, alpha, N, num) {
  print(params <- paste(c(p, h2, alpha, N, num), collapse = "_"))
  list(
    simu_pheno = readRDS(paste0("results_simu/simu_", params, ".rds")),
    `LDpred2-auto_nojump` = readRDS(paste0("results_simu/LDpred2_nojump_", params, ".rds"))
  )
})

grid3 <- grid %>%
  filter(num == 1) %>%
  select(-res) %>%
  tidyr::unnest_wider(res2) %>%
  tidyr::unnest_wider(`LDpred2-auto_nojump`) %>%
  mutate(n_keep = purrr::map_dbl(all_h2, ~ `if`(is.null(.), 0, NCOL(.))))

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
  select(p:N, n_keep:prop) %>%
  filter(alpha == -0.5) %>%
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
  theme(legend.position = "top", legend.key.width = unit(40, "pt")) +
  guides(color = guide_colorbar(title.vjust = 0.9)) +
  labs(x = "Mean of posterior inclusion probabilities (in bin)",
       y = "Proportion of causal variants (in bin)", color = "# chains")
# ggsave("figures/postp_calib.pdf", width = 10, height = 7)

grid3 %>%
  select(p:N, n_keep, all_local_h2, all_true_h2) %>%
  filter(alpha == -0.5) %>%
  mutate(across(all_local_h2, ~ lapply(., colMeans))) %>%
  tidyr::unnest(c(all_local_h2, all_true_h2)) %>%
  ggplot(aes(all_local_h2 / h2, all_true_h2 / h2, color = n_keep)) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50)) +
  coord_equal() +
  bigstatsr::theme_bigstatsr(0.7) +
  geom_point() +
  geom_abline(lty = 2, lwd = 1, color = "chartreuse2") +
  facet_grid(h2 ~ N + p, labeller = label_both) +
  theme(legend.position = "top", legend.key.width = unit(40, "pt")) +
  guides(color = guide_colorbar(title.vjust = 0.9)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.04),
                     minor_breaks = seq(0, 1, by = 0.01)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.02)) +
  labs(x = TeX("Estimated per-block heritability  (as percentage of $h^2$)"),
       y = TeX("Simulated per-block heritability  (as percentage of $h^2$)"),
       color = "# chains")
# ggsave("figures/local_h2_calib.pdf", width = 10, height = 7)


#### Multiple simulations to validate CIs ####

grid4 <- grid %>%
  select(-any_of("res2")) %>%
  filter(alpha == -0.5) %>%
  tidyr::unnest_wider("res") %>%
  tidyr::pivot_longer(-(1:5), names_to = "Method") %>%
  filter(Method != "LDpred2-auto_altfilter") %>%
  tidyr::unnest_wider("value") %>%
  mutate(r2 = lapply(r2, unname),
         n_keep = purrr::map_dbl(all_h2, ~ `if`(is.null(.), 0, NCOL(.))),
         Method = factor(Method, levels = Method[1:5]),
         p = as.factor(p)) %>%
  relocate(n_keep, .after = Method) %>%
  print(n = 20)

grid4_r2 <- grid4 %>%
  filter(Method != "LDSC") %>%
  mutate(r2_est = purrr::map(all_r2, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("r2",     names_sep = "_") %>%
  tidyr::unnest_wider("r2_est", names_sep = "_")

grid4_r2 %>%
  filter(N == 20e3) %>%
  ggplot(aes(Method, r2_est_1, color = ifelse(Method == "SBayesS", NA, n_keep))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_point(alpha = 0) +
  facet_grid(h2 + p ~ N + alpha + num, scales = "free_y", labeller = label_both) +
  geom_errorbarh(aes(y = r2_1, xmin = as.numeric(Method) - 0.42,
                     xmax = as.numeric(Method) + 0.42, height = r2_3 - r2_2),
                 color = "chartreuse3") +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = r2_est_2, ymax = r2_est_3), width = 0.3) +
  theme(legend.position = "top", legend.key.width = unit(50, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  guides(color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = TeX("Inferred $r^2$  (+ 95% CI of the estimate)"),
       x = "Method", color = "# chains")
# ggsave("figures/est_r2_ten_20K.pdf", width = 11, height = 14)
# ggsave("figures/est_r2_ten_200K.pdf", width = 11, height = 14)


grid4 %>%
  filter(N == 20e3) %>%
  mutate(h2_est = purrr::map2(all_h2, Method, ~ {
    `if`(.y == "LDSC", .x[["h2"]] + c(0, -1.96, 1.96) * .x[["h2_se"]],
         unname(quantile(.x, probs = c(0.5, 0.025, 0.975))))
  })) %>%
  tidyr::unnest_wider("h2_est", names_sep = "_") %>%
  ggplot(aes(Method, h2_est_1, color = ifelse(grepl("LDpred2", Method), n_keep, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = h2), color = "chartreuse3", linetype = 2) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = h2_est_2, ymax = h2_est_3), width = 0.5) +
  facet_grid(h2 + p ~ N + alpha + num, scales = "free_y", labeller = label_both) +
  theme(legend.position = "top", legend.key.width = unit(50, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  guides(color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = TeX("Inferred $h^2$  (+ 95% CI of the estimate)"),
       x = "Method", color = "# chains")
# ggsave("figures/est_h2_ten_20K.pdf", width = 11, height = 14)
# ggsave("figures/est_h2_ten_200K.pdf", width = 11, height = 14)


grid4 %>%
  filter(N == 20e3, Method != "LDSC") %>%
  mutate(p_est = purrr::map(all_p, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("p_est", names_sep = "_") %>%
  ggplot(aes(Method, p_est_1, color = ifelse(grepl("LDpred2", Method), n_keep, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = as.double(as.character(p))),
             color = "chartreuse3", linetype = 2) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = p_est_2, ymax = p_est_3), width = 0.5) +
  facet_grid(h2 + p ~ N + alpha + num, scales = "free_y", labeller = label_both) +
  scale_y_log10() +
  theme(legend.position = "top", legend.key.width = unit(50, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  guides(color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred polygenicity  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains")
# ggsave("figures/est_p_ten_20K.pdf",  width = 11, height = 14)
# ggsave("figures/est_p_ten_200K.pdf", width = 11, height = 14)


grid4 %>%
  filter(N == 20e3, !Method %in% c("LDSC", "LDpred2-auto_noMLE")) %>%
  mutate(alpha_est = purrr::map(all_alpha, ~ unname(quantile(., probs = c(0.5, 0.025, 0.975))))) %>%
  tidyr::unnest_wider("alpha_est", names_sep = "_") %>%
  ggplot(aes(Method, alpha_est_1, color = ifelse(grepl("LDpred2", Method), n_keep, NA))) +
  scale_color_gradient(high = "#0072B2", low = "#D55E00", limits = c(0, 50),
                       na.value = "black") +
  bigstatsr::theme_bigstatsr(0.65) +
  geom_hline(aes(yintercept = alpha), color = "chartreuse3", linetype = 2) +
  geom_hline(yintercept = c(-1.5, 0.5), color = "black", linetype = 3) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = alpha_est_2, ymax = alpha_est_3), width = 0.5) +
  facet_grid(h2 + p ~ N + alpha + num, scales = "free_y", labeller = label_both) +
  theme(legend.position = "top", legend.key.width = unit(50, "pt"),
        axis.text.x = element_text(angle = 40, vjust = 0.95, hjust = 0.95)) +
  guides(color = guide_colorbar(title.vjust = 0.9)) +
  labs(y = "Inferred alpha  (+ 95% CI of the estimate)",
       x = "Method", color = "# chains")
# ggsave("figures/est_alpha_ten_20K.pdf",  width = 10, height = 14)
# ggsave("figures/est_alpha_ten_200K.pdf", width = 10, height = 14)
