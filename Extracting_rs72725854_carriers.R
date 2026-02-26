## Extracting rs72725854 carriers

remotes::install_github("lcpilling/ukbrapR")
varlist <- data.frame(rsid=c("rs72725854"), chr=c(8))
imputed_genotypes <- ukbrapR:::extract_variants(varlist)

write.csv(imputed_genotypes, file = "imputed_rs72725854.csv")
system(paste("dx upload", "imputed_rs72725854.csv"))

