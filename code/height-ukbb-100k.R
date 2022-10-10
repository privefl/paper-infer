library(bigsnpr)
library(dplyr)

pheno_file <- "data/ukbb-quant-pheno/height.rds"
set <- "hm3_plus"

# Run the GWAS with less people (100K only)
ukbb <- snp_attach(paste0("data/UKBB_", set, ".rds"))
G <- ukbb$genotypes
ind_csv <- ukbb$fam$ind_csv
covar <- ukbb$fam %>% select(PC1:date) %>% as.matrix()
y <- readRDS(pheno_file)[ind_csv] + 0
set.seed(1); ind.train <- sort(sample(which(complete.cases(covar) & !is.na(y)), 100e3))


NCORES <- 14

gwas <- runonce::save_run({

  library(future.batchtools)
  plan(batchtools_slurm(resources = list(
    t = "12:00:00", c = NCORES + 2, mem = "50g",
    name = basename(rstudioapi::getSourceEditorContext()$path))))
  intervals <- bigparallelr::split_len(ncol(G), nb_split = 30)

  furrr::future_map_dfr(rows_along(intervals), function(block) {
    big_univLinReg(G, y[ind.train], ind.train = ind.train, covar.train = covar[ind.train, ],
                   ind.col = seq(intervals[block, "lower"], intervals[block, "upper"]),
                   ncores = NCORES)
  })
}, file = "tmp-data/GWAS-height-100K.rds")

df_beta <- gwas %>%
  transmute(beta = estim, beta_se = std.err, n_eff = length(ind.train))

corr <- readRDS(paste0("data/corr_", set, ".rds"))


# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
(ldsc <- snp_ldsc2(corr, df_beta, blocks = 200, intercept = NULL, ncores = NCORES))
#         int      int_se          h2       h2_se
# 1.018185266 0.007567497 0.645516599 0.027348155


# LDpred2-auto
coef_shrink <- 0.95
multi_auto <- runonce::save_run(
  snp_ldpred2_auto(corr, df_beta, h2_init = pmax(ldsc[["h2"]], 0.001),
                   vec_p_init = seq_log(1e-4, 0.2, length.out = 50),
                   burn_in = 500, num_iter = 500, report_step = 20,
                   allow_jump_sign = TRUE, shrink_corr = coef_shrink,
                   ncores = NCORES),
  file = "tmp-data/mod_height_100K.rds"
)

# derive other results
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)) & range < 3))

hist(all_h2    <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est,    500)))
unname(round(100 * quantile(all_h2, c(0.5, 0.025, 0.975)), 1))  # 60.2 [57.2, 63.2]
hist(all_p     <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est,     500)))
unname(round(100 * quantile(all_p, c(0.5, 0.025, 0.975)), 1))  # 1.1 [1.0, 1.4]
hist(all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500)))
unname(round(quantile(all_alpha, c(0.5, 0.025, 0.975)), 2))  # -0.71 [-0.75, -0.67]

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
unname(round(100 * quantile(all_r2, c(0.5, 0.025, 0.975)), 1)) # 29.6 [28.7, 30.5]


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

# saveRDS(res_ldpred2, "tmp-data/mod_height_enrich_100k.rds")
