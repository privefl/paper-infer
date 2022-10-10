NCORES <- 50
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 10, mem = "200g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

# plan("multisession", workers = 12)

new_map <- furrr::future_map_dfr(1:22, function(chr) {

  print(chr)

  library(bigsnpr)
  bkfile <- paste0("tmp-data/hm3_plus/ukbb_hm3_plus_chr", chr, ".bk")
  rdsfile <- sub_bk(bkfile, ".rds")

  if (!file.exists(rdsfile)) {

    mfi <- bigreadr::fread2(paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt"))

    library(dplyr)
    snp_id <- mfi %>%
      filter(V6 > 0.005) %>%
      with(paste(chr, V3, V4, V5, sep = "_"))

    df0 <- bigreadr::fread2(
      "UKBB/ukb41181.csv",
      select = c("eid", "22020-0.0", "22006-0.0"),
      col.names = c("eid", "used_in_pca", "white_british")
    )

    # Individuals still in data
    sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
    ind.indiv <- match(df0$eid, sample$ID_2)
    id_to_rm <- scan("UKBB/w58024_20220222.csv")
    ind.indiv[df0$eid %in% id_to_rm] <- NA
    sub <- which(!is.na(ind.indiv) & df0$used_in_pca & is.na(df0$white_british))
    length(sub) # 69,630

    # Read dosage data
    snp_readBGEN(
      bgenfiles = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = chr),
      list_snp_id = list(snp_id),
      backingfile = sub_bk(bkfile),
      ind_row = ind.indiv[sample(sub, 10e3)],
      ncores = NCORES
    )
  }

  ukbb <- snp_attach(rdsfile)
  map <- dplyr::transmute(ukbb$map, chr = as.integer(chromosome), pos = physical.pos,
                          a1 = allele1, a0 = allele2, freq, info)

  corr0 <- runonce::save_run(
    snp_cor(ukbb$genotypes, infos.pos = ukbb$map$physical.pos,
            size = 1000, thr_r2 = 0.3, ncores = NCORES),
    file = paste0("tmp-data/hm3_plus/corr_hm3_plus_chr", chr, ".rds"))

  res <- runonce::save_run({

    corr2 <- as(corr0, "dgCMatrix")

    map_hm3 <- readRDS("tmp-data/map_hm3_ldpred2.rds")
    ind_hm3 <- snp_match(cbind(map_hm3, beta = 1), map, match.min.prop = 0)[["_NUM_ID_"]]
    select <- rep(FALSE, ncol(corr2)); select[ind_hm3] <- TRUE
    exclude <- with(map, nchar(a0) > 1 | nchar(a1) > 1 | pmin(freq, 1 - freq) < 0.005 |
                      info < 0.3 | paste(a0, a1) %in% c("A T", "T A", "C G", "G C"))

    Rcpp::sourceCpp('code/greedy-maxtag.cpp')
    greedy_maxtag(corr2@p, corr2@i, corr2@x, min_add = 2, remove_diag = FALSE,
                  select = select, exclude = exclude, sqrt_info = sqrt(map$info))

  }, file = paste0("tmp-data/hm3_plus/set_hm3_plus_chr", chr, ".rds"))

  map[which(!is.na(res[[1]])), ]
})


nrow(new_map)  # 1,444,196
table(new_map$chr)
#      1      2      3      4      5      6      7      8      9     10     11
# 117534 120469 100666  94010  91117  94264  81254  78428  65806  75354  72406
#    12     13     14     15     16     17     18     19     20     21     22
# 70201  53449  47714  43091  46228  40436  42976  30552  36426  20582  21233

list_snp_id <- with(new_map, split(paste(chr, pos, a1, a0, sep = "_"), chr))

library(bigsnpr)

# same individuals as in simu
fam <- snp_attach("data/UKBB_hm3.rds")$fam

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_hm3_plus",
    ind_row     = fam$ind_bgen,
    ncores      = nb_cores()
  )
) # 65 min with 15 cores

ukbb <- snp_attach("data/UKBB_hm3_plus.rds")
ukbb$fam <- fam

snp_save(ukbb)
