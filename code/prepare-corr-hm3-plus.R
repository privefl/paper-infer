library(bigsnpr)
library(dplyr)
ukb <- snp_attach("data/UKBB_hm3_plus.rds")
G     <- ukb$genotypes
CHR   <- as.integer(ukb$map$chromosome)
POS   <- ukb$map$physical.pos
POS2  <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = nb_cores())

library(future.batchtools)
NCORES <- 14
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("data/corr_hm3_plus")

library(furrr)
all_final_grp <- future_map_dfr(1:22, function(chr) {

  ind.chr <- which(CHR == chr)

  # LD for LDpred2
  corr0 <- runonce::save_run(
    snp_cor(G, ind.row = sort(sample(nrow(G), 20e3)), ind.col = ind.chr,
            infos.pos = POS2[ind.chr], size = 3 / 1000, ncores = NCORES),
    file = paste0("data/corr_hm3_plus/chr", chr, ".rds")
  )

  # find nearly independent LD blocks
  m <- ncol(corr0)
  (SEQ <- round(seq_log(m / 30, m / 5, length.out = 20)))
  splits <- snp_ldsplit(corr0, thr_r2 = 0.02, min_size = 200, max_size = SEQ, max_r2 = 0.2)

  best_split <- splits %>%
    arrange(cost2 * sqrt(5 + cost)) %>%
    print() %>%
    slice(1) %>%
    print()

  library(ggplot2)
  plot_grid(
    qplot(data = splits, perc_kept, cost, color = as.factor(max_size)) +
      geom_point(data = best_split, size = 2, color = "black") +
      theme_bw(12) +
      theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
      scale_x_continuous(limits = c(0, NA), breaks = seq(0, 1, by = 0.1)) +
      labs(y = "Sum of squared correlations outside blocks",
           x = "% of non-zero values kept", color = "Maximum block size"),
    qplot(data = splits, cost2, cost, color = as.factor(max_size)) +
      geom_point(data = best_split, size = 2, color = "black") +
      theme_bw(12) +
      theme(legend.position = "none") +
      scale_y_log10() +
      xlim(0, NA) +
      labs(y = "Sum of squared correlations outside blocks",
           x = "Sum of squared blocks", color = "Maximum block size"),
    scale = 0.95, ncol = 1, rel_heights = c(1.18, 1)
  )

  (all_size <- best_split$all_size[[1]])
  best_grp <- rep(seq_along(all_size), all_size)

  runonce::save_run({
    corr0T <- as(corr0, "dgTMatrix")
    corr0T@x <- ifelse(best_grp[corr0T@i + 1L] == best_grp[corr0T@j + 1L], corr0T@x, 0)
    as(Matrix::drop0(corr0T), "symmetricMatrix")
  }, file = paste0("data/corr_hm3_plus/LD_with_blocks_chr", chr, ".rds"))

  # return
  best_split
}, .options = furrr_options(scheduling = FALSE))


# verif
plot(all_final_grp$n_block)
sum(all_final_grp$n_block) # 431
plot(all_final_grp$cost)
plot(all_final_grp$cost2)
# saveRDS(all_final_grp, "data/corr_hm3_plus/all_final_grp.rds")

sum(file.size(paste0("data/corr_hm3_plus/LD_with_blocks_chr", 1:22, ".rds"))) /
  sum(file.size(paste0("data/corr_hm3_plus/chr", 1:22, ".rds")))
# 89%
