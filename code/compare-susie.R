library(bigsnpr)
NCORES <- nb_cores()

load("data/ukbb4simu_ind.RData")
ind.gwas0 <- ind.gwas

simu <- snp_attach("data/ukbb4simu.rds")
G <- simu$genotypes

p <- 0.002
# p <- 0.02
h2 <- 0.1
alpha <- -0.5
N <- 100e3

#### Run GWAS and prepare sumstats ####

simu_pheno <- snp_simuPheno(G, h2 = h2, M = round(ncol(G) * p), alpha = alpha,
                            ncores = NCORES)

# GWAS to get sumstats
ind.gwas <- sort(sample(ind.gwas0, N))
gwas <- big_univLinReg(G, simu_pheno$pheno[ind.gwas], ind.train = ind.gwas, ncores = NCORES)
df_beta <- dplyr::transmute(gwas, beta = estim, beta_se = std.err, n_eff = length(ind.gwas))


# Run SuSiE-RSS

blocks <- subset(readRDS("../misspec/ldref/map.rds"), chr %% 3 == 0)$group_id
list_ind <- split(seq_along(blocks), blocks)
hist(lengths(list_ind))

bigassertr::assert_dir("tmp-data/corr_susie")
bigparallelr::set_blas_ncores(NCORES)

purrr::iwalk(list_ind, ~ {
  print(.y)
  res_file <- paste0("tmp-data/corr_susie/block_", .y, ".rds")
  runonce::skip_run_if({
    R <- big_cor(G, ind.col = .x, ind.row = sample(nrow(G), 20e3))[]
    saveRDS(R, res_file)
  }, files = res_file)
})

Z <- with(df_beta, beta / beta_se)
N_gwas <- length(ind.gwas)
system.time(
  postp2 <- unlist(purrr::imap(list_ind, ~ {
    print(.y)
    R <- readRDS(paste0("tmp-data/corr_susie/block_", .y, ".rds"))
    # print(isSymmetric(R))
    susieR::susie_rss(Z[.x], R, n = N_gwas, L = 10)$pip
  }))
)
#     user   system  elapsed
# 3741.050 3300.362  897.720
# 3479.697 3059.230  770.534


# Run LDpred2-auto

corr <- readRDS("data/corr_simu.rds")
ldsc <- snp_ldsc2(corr, df_beta)
system.time(
  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc[["h2"]], ncores = NCORES,
                                 vec_p_init = seq_log(1e-4, 0.2, length.out = NCORES),
                                 burn_in = 100, num_iter = 100,
                                 allow_jump_sign = FALSE, shrink_corr = 0.95)
)
#  user  system elapsed
# 2.657   1.730  59.752
# 2.625   2.084  47.646
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- (range > (0.95 * quantile(range, 0.95))))

postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp_est))


# Compare results

plot(postp, postp2, col = seq_along(postp) %in% simu_pheno$set + 1,
     log = "xy", pch = 20); abline(0, 1, col = "red")

ind <- split(seq_along(postp), cut(postp, seq_log(1e-4, 1, 17)))
mid <- sapply(ind, function(.) mean(postp[.]))
prop <- sapply(ind, function(.) {
  is_causal <- . %in% simu_pheno$set
  n <- length(is_causal)
  if (n < 2) return(rep(NA_real_, 3))
  unlist(prop.test(sum(is_causal), n)[c("estimate", "conf.int")], use.names = FALSE)
})

ind2 <- split(seq_along(postp2), cut(postp2, seq_log(1e-4, 1, 17)))
mid2 <- sapply(ind2, function(.) mean(postp2[.]))
prop2 <- sapply(ind2, function(.) {
  is_causal <- . %in% simu_pheno$set
  n <- length(is_causal)
  if (n < 2) return(rep(NA_real_, 3))
  unlist(prop.test(sum(is_causal), n)[c("estimate", "conf.int")], use.names = FALSE)
})

df <- rbind(
  cbind.data.frame(Method = "LDpred2-auto", mid, prop = prop[1, ], lo = prop[2, ], up = prop[3, ]),
  cbind.data.frame(Method = "SuSiE-RSS", mid = mid2, prop = prop2[1, ], lo = prop2[2, ], up = prop2[3, ])
)

library(ggplot2)
ggplot(df, aes(mid, prop, color = Method)) +
  geom_point() +
  geom_abline(lty = 2, alpha = 0.3) +
  theme_bw(14) +
  scale_y_log10(breaks = signif(seq_log(1e-4, 1, 9), 1)) +
  scale_x_log10(breaks = signif(seq_log(1e-4, 1, 9), 1)) +
  coord_equal() +
  labs(x = "Mean of posterior inclusion probabilities (in bin)",
       y = "Proportion of causal variants (in bin)") +
  geom_errorbar(aes(ymin = lo, ymax = up), width = 0.1) +
  theme(legend.position = c(0.25, 0.8)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9"))

# ggsave("figures/susie_02.pdf",  width = 10, height = 6)
# ggsave("figures/susie_002.pdf", width = 8,  height = 7)
