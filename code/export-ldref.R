library(bigsnpr)
# library(bigreadr)
# library(dplyr)
# library(ggplot2)

#### Filter data from allele frequency errors ####

ukb <- snp_attach("data/UKBB_hm3_plus.rds")
G <- ukb$genotypes
map <- transmute(ukb$map,
                 chr = as.integer(chromosome), pos = physical.pos,
                 a0 = allele1, a1 = allele2, rsid)  # reversed somehow..
vctrs::vec_duplicate_any(map[, 1:2])  # FALSE

NCORES <- nb_cores()
map$af_UKBB <- runonce::save_run(
  big_colstats(G, ncores = NCORES)$sum / (2 * nrow(G)),
  file = "tmp-data/af-ldref.rds")


# Compute LD scores
map$ld <- do.call('c', lapply(1:22, function(chr) {
  cat(chr, ".. ", sep = "")
  corr_chr <- readRDS(paste0("data/corr_hm3_plus/chr", chr, ".rds"))
  Matrix::colSums(corr_chr^2)
}))


# Add block IDs
block_size <- unlist(readRDS("data/corr_hm3_plus/all_final_grp.rds")$all_size)
map$block_id <- rep(seq_along(block_size), block_size)


# add positions in different builds
liftOver <- runonce::download_file(
  "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver", "tmp-data")
map$pos_hg18 <- snp_modifyBuild(map, liftOver, from = "hg19", to = "hg18")$pos
map$pos_hg38 <- snp_modifyBuild(map, liftOver, from = "hg19", to = "hg38")$pos

saveRDS(map, "map_hm3_plus.rds", version = 2)

setwd("data/corr_hm3_plus")
for(chr in 1:22) {
  ld_file <- paste0("LD_with_blocks_chr", chr, ".rds")
  system(paste("zip -0 ../../ldref_hm3_plus.zip", ld_file))
}
# 14.1 GB
