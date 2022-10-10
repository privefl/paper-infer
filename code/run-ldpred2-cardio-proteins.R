library(dplyr)
urls <- readLines("https://zenodo.org/record/2615265#.YvuAdBxBy01") %>%
  grep(value = TRUE, pattern = "type=\"application/gzip\" href=\"https://zenodo.org/record/2615265/files/.+\\.txt\\.gz\">") %>%
  sub(pattern = ".+href=\"(https://zenodo.org/record/2615265/files/.+\\.txt\\.gz)\">", replacement = "\\1")

bigassertr::assert_dir("cardio_proteins")

NCORES <- 13
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

# future::plan("sequential")

res <- furrr::future_map(urls, function(url) {

  protein <- sub("\\.txt\\.gz$", "", basename(url))
  res_file <- paste0("cardio_proteins/", protein, ".rds")
  runonce::save_run(file = res_file, {

    sumstats <- bigreadr::fread2(runonce::download_file(url, dir = "tmp-data/SCALLOP_CVD1")) %>%
      transmute(chr = as.integer(sub("(.+):.+:.+_.+", "\\1", MarkerName)),
                pos = as.integer(sub(".+:(.+):.+_.+", "\\1", MarkerName)),
                across(c(a1 = Allele1, a0 = Allele2), toupper),
                freq = Freq1, beta = Effect, beta_se = StdErr, n_eff = TotalSampleSize)

    library(bigsnpr)
    map <- snp_attach("data/UKBB_hm3_plus.rds")$map %>%
      transmute(chr = as.integer(chromosome), pos = physical.pos, a0 = allele2, a1 = allele1, freq)

    (info_snp <- tibble::as_tibble(snp_match(sumstats, map)))

    maf <- map$freq[info_snp$`_NUM_ID_`]
    sd_val <- sqrt(2 * maf * (1 - maf))
    sd_ss <- with(info_snp, 1 / sqrt(n_eff * beta_se^2 + beta^2))
    is_bad <-
      sd_ss < (0.6 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.05 | sd_val < 0.05
    mean(is_bad)

    df_beta <- info_snp[!is_bad, ] %>%
      filter(n_eff > 0.7 * max(n_eff))

    tmp <- tempfile(tmpdir = "tmp-data")
    on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

    for (chr in 1:22) {

      cat(chr, ".. ", sep = "")

      ## indices in 'df_beta'
      ind.chr <- which(df_beta$chr == chr)
      ## indices in 'map'
      ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
      ## indices in 'corr_chr'
      ind.chr3 <- match(ind.chr2, which(map$chr == chr))

      corr_chr <- readRDS(paste0("data/corr_hm3_plus/LD_with_blocks_chr", chr, ".rds"))
      corr_sub <- corr_chr[ind.chr3, ind.chr3]

      if (chr == 1) {
        corr <- as_SFBM(corr_sub, tmp, compact = TRUE)
      } else {
        corr$add_columns(corr_sub, nrow(corr))
      }
    }


    # Heritability estimation of LD score regression
    # to be used as a starting value in LDpred2-auto
    (ldsc <- snp_ldsc2(corr, df_beta, blocks = 200, intercept = NULL, ncores = NCORES))


    # LDpred2-auto
    coef_shrink <- 0.95
    time <- system.time(
      multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                                     vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                                     burn_in = 500, num_iter = 500, report_step = 20,
                                     allow_jump_sign = FALSE, shrink_corr = coef_shrink,
                                     ncores = NCORES)
    )

    # derive other results
    (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
    (keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

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

    # local h2 (per-block)
    all_size <- unlist(readRDS(paste0("data/corr_hm3_plus/all_final_grp.rds"))$all_size)
    blocks <- rep(seq_along(all_size), all_size)[df_beta$`_NUM_ID_`]
    list_ind <- split(seq_along(blocks), blocks)
    bsamp <- do.call("cbind", bsamp)
    all_local_h2 <- sapply(list_ind, function(ind) {
      # cat(".")
      corr_sub <- corr[ind, ind]
      bsamp_sub <- bsamp[ind, ]
      Rb <- coef_shrink * corr_sub %*% bsamp_sub + (1 - coef_shrink) * bsamp_sub
      Matrix::colSums(bsamp_sub * Rb)
    })
    local_h2 <- unname(colMeans(all_local_h2))
    plot(local_h2 / sum(local_h2))

    # per-variant prob causal
    postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp_est))

    list2 <- function(...) setNames(list(...), sapply(substitute(list(...))[-1], deparse))
    list2(ldsc, postp, all_h2, all_p, all_alpha, all_r2, all_local_h2, time, is_bad)
  })
})


library(dplyr)

res2 <- tibble(protein = sub("\\.txt\\.gz$", "", basename(urls)), res) %>%
  tidyr::unnest_wider("res")

summary(sapply(res2$is_bad, mean))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.001200 0.005005 0.005110 0.007247 0.008492 0.059975

res_local_h2 <- res2 %>%
  rowwise() %>%
  mutate(local_h2 = list(colMeans(all_local_h2)),
         h2 = sum(local_h2),
         max_local_h2 = max(local_h2)) %>%
  ungroup() %>%
  mutate(protein = factor(protein, levels = protein[order(h2)]))

with(res_local_h2, sum((max_local_h2 / h2) > 0.5))  # 22
with(res_local_h2, sum((max_local_h2 / h2) > 0.8))  # 8

library(ggplot2)
ggplot(tidyr::pivot_longer(res_local_h2, c(h2, max_local_h2)),
       aes(protein, value, fill = name, color = name)) +
  theme_bw(13) +
  scale_color_manual(values = c("#0072B2FF", "#00000033")) +  # "#E69F0000"
  scale_fill_manual(values = c("#0072B200", "#E69F0099")) +
  geom_col(position = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05), minor_breaks = seq(0, 1, by = 0.01)) +
  labs(y = NULL, fill = NULL, color = NULL)
# ggsave("figures/protein_local_h2.pdf", width = 13, height = 7)
