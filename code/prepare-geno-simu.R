library(dplyr)

# Get variant info
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/36360900",
  dir = "tmp-data", fname = "map_hm3_ldpred2.rds")) %>%
  filter(chr %% 3 == 0)

list_snp_id <- with(info, split(paste(chr, pos, a0, a1, sep = "_"), chr))


# Get indiv info

df0 <- bigreadr::fread2(
  "UKBB/ukb41181.csv",
  select    = c("eid", "22020-0.0",   paste0("22009-0.", 1:16)),
  col.names = c("eid", "used_in_pca", paste0("PC", 1:16))
)

# Individuals still in data
sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
ind.indiv <- match(df0$eid, sample$ID_2)
id_to_rm <- scan("UKBB/w58024_20220222.csv")
ind.indiv[df0$eid %in% id_to_rm] <- NA
sub <- which(!is.na(ind.indiv) & df0$used_in_pca)

# Genetically homogeneous
dist <- bigutilsr::dist_ogk(as.matrix(df0[sub, -(1:2)]))
hist(log(dist), "FD")
sub2 <- sub[log(dist) < 4.5]
length(sub2) # 356,409


# Read dosage data
library(bigsnpr)
system.time(
  snp_readBGEN(
    bgenfiles = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = names(list_snp_id)),
    list_snp_id = list_snp_id,
    backingfile = "data/ukbb4simu",
    ind_row = ind.indiv[sub2],
    ncores = nb_cores()
  )
) # 14 min

simu <- snp_attach("data/ukbb4simu.rds")
simu$map$chromosome <- as.integer(simu$map$chromosome)
simu$map$ld <- info$ld
simu$fam <- data.frame(eid = df0$eid[sub2],
                       ind_csv = sub2,
                       ind_bgen = ind.indiv[sub2])
snp_save(simu)

set.seed(1)
ind.gwas <- sort(sample(length(sub2), 200e3))
ind.test <- setdiff(seq_along(sub2), ind.gwas)
save(ind.gwas, ind.test, file = "data/ukbb4simu_ind.RData")
