## Extracting rs72725854 carriers

remotes::install_github("lcpilling/ukbrapR")
varlist <- data.frame(rsid=c("rs72725854_T", "rs72725854_G"), pos=c("128074815"), chr=c(8))
imputed_genotypes <- ukbrapR:::extract_variants(varlist, use_imp_pos=TRUE, overwrite=TRUE, 
                                                progress=TRUE, verbose=TRUE, very_verbose=TRUE)

#############################################
# Important - remove withdrawn participants #
#############################################

exclude_withdrawn=function(df){
  system('dx download Callum/Withdrawals/withdrawn_20260310.csv --overwrite') ## This file is a list of participants who withdrew from the Biobank up to the date 10th March 2026. This was sent from the UK Biobank team via email to members of approved applications
  df2 = df %>% left_join(
    read_csv("withdrawn_20260310.csv", col_names = FALSE, show_col_types = FALSE) %>%
      mutate(w=1) %>%
      rename(eid=X1),
    by='eid'
  ) %>%
    filter(is.na(w))
  return(df2)
}

imputed_genotypes <- exclude_withdrawn(imputed_genotypes)

write.csv(imputed_genotypes, file = "imputed_rs72725854.csv")
system(paste("dx upload", "imputed_rs72725854.csv"))

