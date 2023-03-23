library(bigsnpr)
library(dplyr)
map <- snp_attach("data/UKBB_hm3_plus.rds")$map
CHR <- as.integer(map$chromosome)
POS <- map$physical.pos

all_final_grp <- readRDS("data/corr_hm3_plus/all_final_grp.rds")
grid <- all_final_grp %>%
  select(all_size) %>%
  mutate(chr = 1:22) %>%
  tidyr::unnest_longer(all_size) %>%
  mutate(
    ind_group = {
      id_group <- rep(seq_along(all_size), all_size)
      split(seq_along(id_group), id_group)
    },
    id_group = row_number(),
    all_size = NULL
  ) %>%
  # filter(!file.exists(paste0("data/corr_hm3_plus_regul/LD_part", id_group, ".rds"))) %>%
  print()


bigassertr::assert_dir("data/corr_hm3_plus_regul")

library(furrr)
plan("multisession", workers = 5)
furrr::future_pwalk(grid, function(chr, ind_group, id_group) {

  # id_group <- 1
  # chr <- grid$chr[id_group]
  # ind_group <- grid$ind_group[[id_group]]

  library(Matrix)
  corr_block <- runonce::save_run({
    ind.chr <- which(CHR == chr)
    ind.block <- match(ind_group, ind.chr)
    readRDS(paste0("data/corr_hm3_plus/LD_with_blocks_chr", chr, ".rds"))[ind.block, ind.block]
  }, file = paste0("data/corr_hm3_plus_regul/LD_part", id_group, ".rds"))

  new_corr <- runonce::save_run({
    corrT <- as(corr_block, "dgTMatrix")
    POS2 <- snp_asGeneticPos(CHR[ind_group], POS[ind_group], dir = "tmp-data")
    corrT@x <- corrT@x * exp(-0.5 * abs(POS2[corrT@i + 1L] - POS2[corrT@j + 1L]))
    corrT %>%
      as("dgCMatrix") %>%
      as("symmetricMatrix")
  }, file = paste0("data/corr_hm3_plus_regul/LDregul_part", id_group, ".rds"))

}, .options = furrr_options(scheduling = FALSE))


res_files <- paste0("data/corr_hm3_plus_regul/LDregul_part", 1:sum(all_final_grp$n_block), ".rds")
sum(file.size(res_files)) /
  sum(file.size(paste0("data/corr_hm3_plus/LD_with_blocks_chr", 1:22, ".rds")))
# 1.005

saveRDS(all_final_grp, "data/corr_hm3_plus_regul/all_final_grp.rds")

rm(corr)
corr <- runonce::save_run({

  for (file in res_files) {

    corr_block <- readRDS(file)

    if (!exists("corr")) {
      corr <- as_SFBM(corr_block, "data/corr_hm3_plus_regul", compact = TRUE)
    } else {
      corr$add_columns(corr_block, nrow(corr))
    }
  }
  corr
}, file = "data/corr_hm3_plus_regul.rds")
print(dim(corr))
