#==========================================#
#=== Create dataframe of Black men only ===#
#==========================================#

## Note: Requires running "Logistic_Regressions_testing_Conti_GRS" down to creation of "PCa_iv_covariates_GRS_severity"

## This is stored as "BlackMen.tsv.gz" in Callum/Derived_datasets

## To create "AllMen.tsv.gz" as found in Callum/Derived_datasets, comment out "ethnicity_group_narrow == "Black"" in the filter below and run the code again. Change any dataframe names if needed.

BlackMen <- PCa_iv_covariates_GRS_severity %>%
  dplyr::filter(
    Sex == "Male",
    ethnicity_group_narrow == "Black",
  ) %>%
  dplyr::mutate(PrCa = if_else(!is.na(pre_diagnosed), 1L, 0L, missing = 0L)) %>%
  dplyr::mutate(FID = eid) %>%
  dplyr::select(c("eid", "FID", "Sex", "ethnicity_group_narrow", "PrCa")) %>%
  dplyr::rename("IID" = "eid")


#############################################
# Important - remove withdrawn participants #
#############################################

exclude_withdrawn=function(df){
  system('dx download Callum/Withdrawals/withdrawn_20260310.csv --overwrite') ## This file is a list of participants who withdrew from the Biobank up to the date 10th March 2026. This was sent from the UK Biobank team via email to members of approved applications
  df2 = df %>% left_join(
    read_csv("withdrawn_20260310.csv", col_names = FALSE, show_col_types = FALSE) %>%
      mutate(w=1) %>%
      rename(IID=X1),
    by='IID'
  ) %>%
    filter(is.na(w))
  return(df2)
}

BlackMen <- exclude_withdrawn(BlackMen)

## write to tsv

write.table(BlackMen, "BlackMen.tsv",quote=FALSE,sep='\t',row.names = FALSE)

## write to gzipped tsv

write.table(
  BlackMen,
  file = gzfile("BlackMen.tsv.gz"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

system(paste("dx upload", "BlackMen.tsv"))
system(paste("dx upload", "BlackMen.tsv.gz")) # Main phenotype file used in Fine Mapping