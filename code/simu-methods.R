################################################################################

list2 <- function(...) rlang::dots_list(..., .named = TRUE)

r2_boot <- function(x, y) {
  all_cor <- replicate(1000, {
    ind <- sample(length(x), replace = TRUE)
    cor(x[ind], y[ind])
  })
  cor <- c(cor(x, y), quantile(all_cor, c(0.025, 0.975)))
  sign(cor) * cor^2
}

conv_stat <- function(multi_auto, which, Niter) {
  # Two-Sample Cramer-von Mises Statistic -- Eq (6) of 10.1214/aoms/1177704477
  sapply(multi_auto, function(auto) {
    x <- tail(auto[[paste0("path_", which, "_est")]], Niter)
    if (anyNA(x)) return(NA_real_)
    r <- rank(x, na.last = "keep", ties.method = "average")
    N <- floor(length(x) / 2)
    M <- length(x) - N
    t <- (N / M * crossprod(sort(r[1:N]) - 1:N * (N + M) / N) +
            M / N * crossprod(sort(r[1:M + N]) - 1:M * (N + M) / M)) / (N + M)^2
    t[[1]]
  })
}

conv_stat_all <- function(multi_auto, Niter) {
  rowMeans(
    sapply(c("h2", "p", "alpha"), conv_stat, multi_auto = multi_auto, Niter = Niter),
    na.rm = TRUE)
}

################################################################################

run_ldpred2 <- function(jump_sign, use_mle, coef_shrink) {

  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                                 vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                                 burn_in = 500, num_iter = 500, report_step = 20,
                                 use_MLE = use_mle,
                                 allow_jump_sign = jump_sign,
                                 shrink_corr = coef_shrink,
                                 ncores = NCORES)

  # derive other results
  (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
  (keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))
  if (length(keep) < 2) return(list(NULL))

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

  # r2 in test set
  beta <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
  pred <- big_prodVec(G, beta, ind.row = ind.test, ncores = NCORES)
  r2 <- r2_boot(pred, simu_pheno$pheno[ind.test])

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

  list2(jump_sign, use_mle, postp, all_h2, all_p, all_alpha, all_r2, r2,
        all_local_h2, all_true_h2)
}

run_ldpred2_new <- function(jump_sign, use_mle, coef_shrink) {

  res0 <- list2(jump_sign, use_mle)

  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                                 vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                                 burn_in = 500, num_iter = 500, report_step = 20,
                                 use_MLE = use_mle,
                                 allow_jump_sign = jump_sign,
                                 shrink_corr = coef_shrink,
                                 ncores = NCORES)

  ## derive other results

  # inferred r2 per chain, to filter
  library(magrittr)
  all_r2_est <- sapply(multi_auto, function(auto) {
    cat(".")
    if (anyNA(auto$beta_est)) return(NA_real_)
    bsamp <- auto$sample_beta
    apply(bsamp, 2, function(x) bigsparser::sp_prodVec(corr, x)) %>%
      Matrix::crossprod(bsamp, .) %>%
      Matrix::tril(k = -5) %>%
      Matrix::drop0() %>%
      .@x %>%
      mean()
  })
  (keep_pred <- which(all_r2_est > (0.95 * quantile(all_r2_est, 0.95, na.rm = TRUE))))

  res_pred <- if (length(keep_pred) < 2) NULL else {

    # r2 in test set
    beta <- rowMeans(sapply(multi_auto[keep_pred], function(auto) auto$beta_est))
    pred <- big_prodVec(G, beta, ind.row = ind.test, ncores = NCORES)
    r2 <- r2_boot(pred, simu_pheno$pheno[ind.test])

    # inferred r2
    bsamp <- lapply(multi_auto[keep_pred], function(auto) auto$sample_beta)
    all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
      # print(ic)
      b1 <- bsamp[[ic]]
      Rb1 <- apply(b1, 2, function(x)
        coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
      b2 <- do.call("cbind", bsamp[-ic])
      b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
    }))
    # hist(all_r2); abline(v = r2, col = "red")

    list2(r2, all_r2)
  }


  (keep_infer <- which(conv_stat_all(multi_auto, 500) < 4))

  res_infer <- if (length(keep_infer) < 2) NULL else {

    (all_h2    <- sapply(multi_auto[keep_infer], function(auto) tail(auto$path_h2_est,    500)))
    (all_p     <- sapply(multi_auto[keep_infer], function(auto) tail(auto$path_p_est,     500)))
    (all_alpha <- sapply(multi_auto[keep_infer], function(auto) tail(auto$path_alpha_est, 500)))

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
    postp <- rowMeans(sapply(multi_auto[keep_infer], function(auto) auto$postp_est))

    list2(all_h2, all_p, all_alpha, all_local_h2, all_true_h2, postp)
  }


  res0 %>% append(res_pred) %>% append(res_infer)
}

################################################################################

# unzip(runonce::download_file(
#   "https://cnsgenomics.com/software/gctb/download/gctb_2.04.3_Linux.zip",
#   dir = "tmp-data"), exdir = "tmp-data")
gctb <- bigsnpr:::make_executable("tmp-data/gctb_2.04.3_Linux/gctb")

run_sbayess <- function() {

  prefix <- "tmp-data/corr_simu.ldm.sparse"
  burnin <- 1000

  tmp <- tempfile(tmpdir = "tmp-data")
  on.exit(unlink(paste0(tmp, "*")), add = TRUE)

  library(dplyr)
  gwas_file <- df_beta %>%
    bind_cols(simu$map) %>%
    transmute(SNP = rsid, A1 = allele1, A2 = allele2, freq, b = beta, se = beta_se,
              p = pchisq((b / se)^2, df = 1, lower.tail = FALSE), N = n_eff) %>%
    bigreadr::fwrite2(paste0(tmp, ".ma"), sep = " ")

  # maybe would need  --constraint="s03|s04|s05"
  system(glue::glue(
    gctb,
    " --ldm {prefix}",
    " --sbayes S",
    " --hsq {pmax(ldsc[['h2']], 0.001)}",
    " --gwas-summary {gwas_file}",
    " --chain-length 2000 --burn-in {burnin}",
    " --out {tmp} --out-freq 100"
  ))

  all_par <- bigreadr::fread2(paste0(tmp, ".mcmcsamples.Par"))
  (all_h2    <- tail(all_par$hsq, -burnin))
  (all_p     <- tail(all_par$Pi,  -burnin))
  (all_alpha <- tail(all_par$S,   -burnin))

  # inferred r2
  file_mcmc_eff <- paste0(tmp, ".mcmcsamples.SnpEffects")
  size <- file.size(file_mcmc_eff) / 4 / 3
  file.rename(file_mcmc_eff, paste0(file_mcmc_eff, ".bk"))

  IJ <- FBM(nrow = 3, ncol = size, type = "integer", backingfile = file_mcmc_eff,
            create_bk = FALSE, is_read_only = TRUE)[1:2, ]
  X <- FBM(nrow = 3, ncol = size, type = "float", backingfile = file_mcmc_eff,
           create_bk = FALSE, is_read_only = TRUE)[3, -1]
  bsamp <- Matrix::sparseMatrix(i = IJ[2, -1], j = IJ[1, -1], x = X,
                                dims = IJ[2:1, 1], index1 = FALSE)[, 101:200] %>%
    sweep(1, with(df_beta, sqrt(n_eff * beta_se^2 + beta^2)), '/')

  all_r2 <- apply(bsamp, 2, function(x) bigsparser::sp_prodVec(corr, x)) %>%
    Matrix::crossprod(bsamp, .) %>%
    Matrix::tril(k = -5) %>%
    Matrix::drop0() %>%
    .@x
  hist(all_r2)
  length(all_r2) # 4560

  # r2 in test set
  res_sbayess <- bigreadr::fread2(paste0(tmp, ".snpRes"))
  ind <- match(res_sbayess$Name, simu$map$rsid)
  stopifnot(identical(ind, seq_along(ind)))
  beta <- res_sbayess$A1Effect
  pred <- big_prodVec(G, beta, ind.row = ind.test, ncores = NCORES)
  r2 <- r2_boot(pred, simu_pheno$pheno[ind.test])
  abline(v = r2, col = "red")

  # posterior inclusion probabilities (i.e. of being causal, per-variant)
  postp <- res_sbayess$PIP

  list2(postp, all_h2, all_p, all_alpha, all_r2, r2)
}

################################################################################
