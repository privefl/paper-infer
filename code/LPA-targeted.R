library(bigsnpr)
library(dplyr)

fam <- snp_attach("data/UKBB_hm3.rds")$fam

# The region used (GRCh37): 6:160952515-161087407 +/- 1 Mb = 6:159952515-162087407
# + filter out variants with MAF < 0.001 and INFO < 0.5
info <- bigreadr::fread2("UKBB/mfi/ukb_mfi_chr6_v3.txt")
snp_id <- info %>%
  filter(V6 > 1e-3, V8 > 0.5, V3 >= 159952515, V3 <= 162087407) %>%
  with(paste(6, V3, V4, V5, sep = "_"))

snp_readBGEN(
  bgenfiles = "UKBB/bgen/ukb_imp_chr6_v3.bgen",
  backingfile = "tmp-data/UKBB_LPA",
  ind_row = fam$ind_bgen,
  list_snp_id = list(snp_id),
  ncores = nb_cores()
)

lpa <- snp_attach("tmp-data/UKBB_LPA.rds")
G <- lpa$genotypes
covar <- as.matrix(select(fam, -eid, -ind_csv, -ind_bgen, -set))
y <- readRDS("data/ukbb-quant-pheno/log_lipoA.rds")[fam$ind_csv]
length(ind_train <- which(!is.na(y) & fam$set == "train" & complete.cases(covar)))  # 232,233
length(ind_col <- with(lpa$map, which(pmin(freq, 1 - freq) > 1e-3, info > 0.5)))  # 13,572

mod <- runonce::save_run(
  big_spLinReg(G, y[ind_train], ind.train = ind_train, ind.col = ind_col, K = 10,
               covar.train = covar[ind_train, ],
               pf.covar = rep(0, ncol(covar)),
               power_scale = c(0, 0.5, 1),
               power_adaptive = c(0, 0.5, 1.5),
               lambda.min.ratio = 1e-6, nlam.min = 30,
               n.abort = 3, ncores = nb_cores()),
  file = "tmp-data/PLR_LPA.rds"
) # 6H

summary(mod)
summary(mod)$message
plot(mod)

ind_test <- which(fam$set == "test")

pred <- predict(mod, G, ind.row = ind_test, covar.row = matrix(0, length(ind_test), ncol(covar)))
pcor(pred, y[ind_test], covar[ind_test, ])^2  # 0.6765771 0.6711782 0.6819082
