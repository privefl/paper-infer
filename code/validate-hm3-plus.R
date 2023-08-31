library(bigsnpr)
library(dplyr)

df0 <- bigreadr::fread2(
  "UKBB/ukb41181.csv",
  select    = c("eid", "22020-0.0",   paste0("22009-0.", 1:16)),
  col.names = c("eid", "used_in_pca", paste0("PC", 1:16))
)
PC_UKBB <- select(df0, PC1:PC16) %>%
  as.matrix()

POPS <- c("United Kingdom", "Iran", "India", "China", "Nigeria")
all_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE) %>%
  filter(Ancestry %in% POPS)
all_sq_dist <- apply(all_centers[-1], 1, function(one_center) {
  rowSums(sweep(PC_UKBB, 2, one_center, '-')^2)
})
thr_sq_dist <- max(dist(all_centers[-1])^2) * 0.002 / 0.16
group <- apply(all_sq_dist, 1, function(x) {
  ind <- which.min(x)
  if (isTRUE(x[ind] < thr_sq_dist)) all_centers$Ancestry[ind] else NA
})
group[-which(df0$used_in_pca == 1)] <- NA  # unrelated and QCed
table(group, exclude = NULL)
# China      India     Iran    Nigeria   United Kingdom       <NA>
#  1769       6019     1246       5017           374515     113939


snpid_hm3_plus <- with(
  snp_attach("data/UKBB_hm3_plus.rds")$map,
  paste(chromosome, physical.pos, allele1, allele2, sep = "_"))

eid_bgen <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")$ID_2[-1]

chr <- 22
NCORES <- nb_cores()

all_tagging <- lapply(POPS, function(pop) {

  print(pop)

  # pop <- POPS[1]

  res_file <- paste0("tmp-data/tagging_chr", chr, "_", substr(pop, 1, 3), ".rds")

  tagging <- runonce::save_run(file = res_file, {

    ind <- which(group == pop)
    if (length(ind) > 10e3) ind <- sample(ind, 10e3)


    mfi <- bigreadr::fread2(paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt"))
    snp_id <- mfi %>%
      filter(V6 > 0.005) %>%
      with(paste(chr, V3, V4, V5, sep = "_"))

    tmp <- tempfile(tmpdir = "tmp-data")
    rds <- bigsnpr::snp_readBGEN(
      bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen"),
      list_snp_id = list(snp_id),
      backingfile = tmp,
      ind_row     = match(df0$eid[ind], eid_bgen),
      ncores      = NCORES
    )
    ukbb_pop_chr <- snp_attach(rds)

    corr0 <- snp_cor(ukbb_pop_chr$genotypes, infos.pos = ukbb_pop_chr$map$physical.pos,
                     size = 1000, thr_r2 = 0.3, ncores = NCORES)

    file.remove(paste0(tmp, c(".bk", ".rds")))

    colMax_sp <- function(X) {
      res <- numeric(ncol(X))
      X2 <- as(X, "TsparseMatrix")
      tmp <- tapply(X2@x, X2@j, max, na.rm = TRUE)
      res[as.integer(names(tmp)) + 1L] <- tmp
      res[res < 0] <- 0  # replace the -Inf from max of nothing
      res
    }

    colMax_sp(corr0[snp_id %in% snpid_hm3_plus, ]^2)
  })

})
names(all_tagging) <- POPS
hist(all_tagging[[4]])

sapply(all_tagging, function(tagging) length(which(tagging > 0.5)) / length(tagging))
round(100 * .Last.value, 1)
# United Kingdom           Iran          India          China        Nigeria
#           82.1           80.0           79.3           75.4           66.6
sapply(all_tagging, function(tagging) length(which(tagging > 0.8)) / length(tagging))
round(100 * .Last.value, 1)
# United Kingdom           Iran          India          China        Nigeria
#           69.1           66.7           66.1           65.3           48.6
