library(bigsnpr)
library(dplyr)

grid <- tibble::tribble(
  ~ pheno,    ~ pheno_file,
  "CAD",      "data/ukbb-binary-pheno/411.4.rds",
  "PrCa",     "data/ukbb-binary-pheno/185.rds",
  "BrCa",     "data/ukbb-binary-pheno/174.1.rds",
  "T1D",      "data/ukbb-binary-pheno/250.1.rds",
  "T2D",      "data/ukbb-binary-pheno/250.2.rds",
  "MDD",      "data/ukbb-binary-pheno/296.2.rds",
  "Asthma",   "data/ukbb-binary-pheno/495.rds",
  "VitaminD", "data/ukbb-quant-pheno/log_vitaminD.rds"
) %>%
  mutate(K = purrr::map_dbl(pheno_file, ~ {
    y <- readRDS(.)
    `if`(is.logical(y), mean(y, na.rm = TRUE), NA) })) %>%
  print() %>%
  #   pheno    pheno_file                                    K
  # 1 CAD      data/ukbb-binary-pheno/411.4.rds        0.0609
  # 2 PrCa     data/ukbb-binary-pheno/185.rds          0.0508
  # 3 BrCa     data/ukbb-binary-pheno/174.1.rds        0.0731
  # 4 T1D      data/ukbb-binary-pheno/250.1.rds        0.00867
  # 5 T2D      data/ukbb-binary-pheno/250.2.rds        0.0609
  # 6 MDD      data/ukbb-binary-pheno/296.2.rds        0.0756
  # 7 Asthma   data/ukbb-binary-pheno/495.rds          0.0763
  # 8 VitaminD data/ukbb-quant-pheno/log_vitaminD.rds NA
  tidyr::expand_grid(set = c("hm3", "hm3_plus", "hm3_plus_regul"),
                     mle = c(TRUE, FALSE),
                     coef_shrink = seq(0.3, 1, by = 0.1)) %>%
  mutate(res_file = paste0("ldpred2-ext/", pheno, "_", set, ifelse(mle, "", "_noMLE"),
                           "_", coef_shrink, ".rds")) %>%
  # filter(!file.exists(res_file)) %>%
  relocate(K, .after = last_col())

bigassertr::assert_dir("ldpred2-ext")

NCORES <- 13
library(future.batchtools)
plan(workers = nrow(grid) + 10, batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

# future::plan("multisession", workers = 14)

grid$res <- furrr::future_pmap(grid[1:6], function(pheno, pheno_file, set, mle,
                                                   coef_shrink, res_file) {

  # pheno <- "T1D"
  # bigassertr::assert_exist(pheno_file <- "data/ukbb-binary-pheno/250.1.rds")
  # set <- "hm3_plus_regul"
  # mle <- TRUE
  # coef_shrink <- 0.9

  runonce::save_run(file = res_file, {

    ukbb <- snp_attach(`if`(grepl("hm3_plus", set), "data/UKBB_hm3_plus.rds",
                            "data/UKBB_hm3.rds"))
    map_ldref <- ukbb$map %>%
      transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele1, a1 = allele2)
    df_beta <- paste0("data/sumstats/", pheno, ".rds") %>%
      readRDS() %>%
      select(chr, pos, a0, a1, beta, beta_se, n_eff) %>%
      snp_match(map_ldref)

    corr <- readRDS(paste0("data/corr_", set, ".rds"))
    ind.sub <- df_beta$`_NUM_ID_`

    blocks <- {
      if (set == "hm3") {
        readRDS("tmp-data/map_hm3_ldpred2.rds")$group_id
      } else {
        all_size <- unlist(readRDS(paste0("data/corr_", set, "/all_final_grp.rds"))$all_size)
        rep(seq_along(all_size), all_size)
      }
    }[ind.sub]

    # Heritability estimation of LD score regression
    # to be used as a starting value in LDpred2-auto
    (ldsc <- snp_ldsc2(corr, df_beta, ind.beta = ind.sub,
                       blocks = blocks, intercept = NULL, ncores = NCORES))
    #        int     int_se         h2      h2_se
    # 0.99550915 0.01242761 0.17885210 0.01171032

    # LDpred2-auto
    time <- system.time(
      multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]],
                                     vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                                     burn_in = 500, num_iter = 500, report_step = 20,
                                     allow_jump_sign = FALSE, shrink_corr = coef_shrink,
                                     use_MLE = mle, ind.corr = ind.sub, ncores = NCORES)
    )

    # derive other results
    (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
    (keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

    if (length(keep) < 2) list(NULL) else {

      hist(all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est,   500)))
      hist(all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est,     500)))
      (all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500)))

      # inferred r2
      bsamp <- lapply(multi_auto[keep], function(auto) {
        b <- auto$sample_beta
        Matrix::sparseMatrix(i = ind.sub[b@i + 1L] - 1L, p = b@p, x = b@x,
                             dims = c(ncol(corr), ncol(b)), index1 = FALSE)
      })
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
      pred <- big_prodVec(ukbb$genotypes, beta, ind.col = ind.sub, ncores = NCORES)

      y_test <- readRDS(pheno_file)[ukbb$fam$ind_csv]
      if (pheno == "VitaminD") y_test <- exp(y_test)
      covar_test <- dplyr::select(ukbb$fam, -(1:4))
      cor <- pcor(pred, y_test, covar_test)
      boot_cor <- replicate(1000, {
        # cat(".")
        ind <- sample(length(pred), replace = TRUE)
        pcor(pred[ind], y_test[ind], covar_test[ind, ])[1]
      })
      cor[2:3] <- quantile(boot_cor, c(0.025, 0.975))
      (r2 <- sign(cor) * cor^2)

      # local h2 (per-block)
      list_ind <- split(seq_along(blocks), blocks)
      bsamp <- do.call("cbind", bsamp)
      all_local_h2 <- sapply(list_ind, function(ind) {
        corr_sub <- corr[ind, ind]
        bsamp_sub <- bsamp[ind, ]
        Rb <- coef_shrink * corr_sub %*% bsamp_sub + (1 - coef_shrink) * bsamp_sub
        Matrix::colSums(bsamp_sub * Rb)
      })
      plot(colMeans(all_local_h2))

      # per-variant prob causal
      postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp_est))

      list2 <- function(...) setNames(list(...), sapply(substitute(list(...))[-1], deparse))
      list2(ldsc, postp, all_h2, all_p, all_alpha, all_r2, r2, all_local_h2, time)

    }
  })
})


library(dplyr)
library(ggplot2)
theme_set(theme_bw(14))
library(latex2exp)

grid2 <- grid %>%
  select(-res_file) %>%
  tidyr::unnest_wider("res", names_sep = "_") %>%
  rename_with(~ sub("res_", "", .)) %>%
  rowwise() %>%
  mutate(n_keep = `if`(is.null(all_h2), 0, ncol(all_h2)))

table(grid2$set)
table(grid2$n_keep)
# 0  33  34  36  39  40  41  43  44  45  46  47  48  49  50
# 9   1   1   2   1   3   2   4   3   2   1   5   2   7 341

scale_to_liab <- function(K, logistic) {
  ifelse(is.na(K), 1,
         purrr::map2_dbl(K, logistic, ~ coef_to_liab(.x, `if`(.y, 0.5, .x))))
}

grid3 <- grid2 %>%
  filter(n_keep > 0) %>%
  mutate(across(c(h2_est = all_h2, p_est = all_p, alpha_est = all_alpha, r2_est = all_r2),
                ~ list(unname(quantile(., probs = c(0.5, 0.025, 0.975), na.rm = TRUE))))) %>%
  ungroup() %>%
  tidyr::pivot_longer(c(r2, r2_est, alpha_est, p_est, h2_est)) %>%
  select(pheno, set, mle, coef_shrink, K, name, value) %>%
  tidyr::unnest_wider("value", names_sep = "_") %>%
  # transformation to the liability scale
  mutate(scaler = case_when(
    name == "r2"     ~ scale_to_liab(K, FALSE),
    name == "r2_est" ~ scale_to_liab(K, TRUE),
    name == "h2_est" ~ scale_to_liab(K, TRUE),
    TRUE ~ 1
  ),
  across(starts_with("value_"), ~ . * scaler))


for (PHENO in unique(grid3$pheno)) {

  grid3 %>%
    filter(pheno == PHENO) %>%
    mutate(mle = ifelse(mle, "Yes", "No"),
           name = factor(name, levels = c(
             "alpha_est", "h2_est", "h2_ldsc_est", "p_est", "r2", "r2_est"),
             labels = sapply(c("$alpha$ (LDpred2-auto)", "$h^2$ (LDpred2-auto)",
                               "$h^2$ (LDSC)", "$p$ (LDpred2-auto)",
                               "$r^2$ (in test set)", "$r^2$ (LDpred2-auto)"), TeX))) %>%
    ggplot(aes(coef_shrink, value_1, color = set, shape = mle, linetype = mle)) +
    bigstatsr::theme_bigstatsr(0.65) +
    facet_wrap(~ pheno + name, scales = "free_y", labeller = label_parsed) +
    geom_point(size = 2) +
    geom_line() +
    theme(legend.position = c(0.85, 0.2), legend.key.width = unit(2, "line")) +
    scale_shape_manual(values = c(3, 16, 17)) +
    scale_color_manual(values = c('#E69F00', '#56B4E9', '#CC79A7')) +
    labs(y = "Estimate", x = "LD shrinkage factor",
         color = "Set", shape = "Use MLE?", linetype = "Use MLE?")

  ggsave(paste0("figures/ldpred2_", PHENO, ".pdf"), width = 9, height = 6)
}
