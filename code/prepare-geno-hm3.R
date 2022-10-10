# Get variant info
info <- readRDS("tmp-data/map_hm3_ldpred2.rds")
list_snp_id <- with(info, split(paste(chr, pos, a0, a1, sep = "_"), chr))

# same individuals as in simu
library(bigsnpr)
fam <- snp_attach("data/ukbb4simu.rds")$fam

system.time(
  rds <- snp_readBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_hm3",
    ind_row     = fam$ind_bgen,
    ncores      = nb_cores()
  )
)

ukbb <- snp_attach("data/UKBB_hm3.rds")

obj.svd <- runonce::save_run(
  snp_autoSVD(ukbb$genotypes, infos.chr = as.integer(ukbb$map$chromosome),
              infos.pos = ukbb$map$physical.pos, ncores = nb_cores()),
  file = "tmp-data/svd_subset.rds")

plot(obj.svd)  # 4
plot(obj.svd, type = "scores", scores = 1:6, coeff = 0.6)  # ok for 4

library(dplyr)
PC <- predict(obj.svd)[, 1:4] %>%
  as.data.frame() %>%
  setNames(paste0("PC", 1:4))

other_covar <- bigreadr::fread2(
  "UKBB/ukb41181.csv",
  select = c("34-0.0", "52-0.0", "22001-0.0", "21022-0.0", "189-0.0"),
  col.names = c("year", "month", "sex", "age", "deprivation_index")
) %>%
  mutate(date = (year - 1900) + (month - 0.5) / 12,
         year = NULL, month = NULL) %>%
  slice(fam$ind_csv)

set <- rep("train", nrow(fam))
set.seed(1); set[sample(length(set), 50e3)] <- "test"

ukbb$fam <- cbind(fam, set, PC, other_covar)

snp_save(ukbb)
