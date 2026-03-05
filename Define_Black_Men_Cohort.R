## Create dataframe of Black men only ##

## Note: Requires running "Logistic_Regressions_testing_Conti_GRS" down to creation of "PCa_iv_covariates_GRS_severity"

BlackMen <- PCa_iv_covariates_GRS_severity %>%
  dplyr::filter(
    Sex == "Male",
    ethnicity_group_narrow == "Black",
  ) %>%
  dplyr::mutate(PrCa = if_else(!is.na(pre_diagnosed), 1L, 0L, missing = 0L)) %>%
  dplyr::mutate(FID = eid) %>%
  dplyr::select(c("eid", "FID", "Sex", "ethnicity_group_narrow", "PrCa")) %>%
  dplyr::rename("IID" = "eid")


write.table(BlackMen, "BlackMen.tsv",quote=FALSE,sep='\t',row.names = FALSE)

## write to gzipped tsv

write.table(
  BlackMen,
  file = gzfile("BlackMen.tsv.gz"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

system(paste("dx upload", "BlackMen.tsv"))
system(paste("dx upload", "BlackMen.tsv.gz")) # Main phenotype file used in Fine Mapping









