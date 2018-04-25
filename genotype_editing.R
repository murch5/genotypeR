library(dplyr)
library(data.table)

missingness_bysite <- data.table::fread("/home/murch/DATA/project/pheno_geno_correlation/geno_data/all_fam_missing_site_v2.lmiss", header=TRUE)

missingness_bysite_over10 <- missingness_bysite %>%
  dplyr::filter(F_MISS < 0.9)

hist(missingness_bysite_over10)