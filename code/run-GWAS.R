bigassertr::assert_dir("GWAS")

pheno_files <- c(list.files("data/ukbb-quant-pheno", full.names = TRUE),
                 list.files("data/ukbb-binary-pheno", full.names = TRUE))

grid <- tidyr::expand_grid(
  set = c("hm3", "hm3_plus"),
  block = 1:30
)

NCORES <- 10
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "50g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pwalk(grid, function(set, block) {

  library(bigsnpr)
  library(dplyr)
  ukbb <- snp_attach(paste0("data/UKBB_", set, ".rds"))
  G <- ukbb$genotypes
  ind_csv <- ukbb$fam$ind_csv
  covar <- ukbb$fam %>% select(PC1:date) %>% as.matrix()
  is_train <- (ukbb$fam$set == "train") & complete.cases(covar)

  intervals <- bigparallelr::split_len(ncol(G), nb_split = 30)
  ind <- seq(intervals[block, "lower"], intervals[block, "upper"])

  lapply(pheno_files, function(pheno_file) {

    pheno <- sub("\\.rds$", "", basename(pheno_file))
    res_file <- paste0("GWAS/", pheno, "_", set, "_part", block, ".rds")

    if (!file.exists(res_file)) {

      y <- readRDS(pheno_file)[ind_csv] + 0
      ind.train <- which(is_train & !is.na(y))

      if (length(table(covar[ind.train, 1])) == 1) {
        COVAR <- covar[, -1]  # rm sex as covariate
      } else {
        COVAR <- covar
      }

      gwas <- big_univLinReg(G, y[ind.train], ind.train = ind.train,
                             covar.train = COVAR[ind.train, ],
                             ind.col = ind, ncores = NCORES)

      saveRDS(gwas, res_file)
    }
  })
})


library(ggplot2)
library(bigsnpr)
ukbb <- snp_attach("data/UKBB_hm3_plus.rds")
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

pheno <- "log_lipoA"
gwas <- purrr::map_dfr(paste0("GWAS/", pheno, "_hm3_plus_part", 1:30, ".rds"), readRDS)

snp_manhattan(gwas, CHR, POS, npoints = 100e3) +
  scale_y_log10() +
  geom_hline(yintercept = -log10(5e-8), color = "red") +
  ggtitle(NULL, NULL)
