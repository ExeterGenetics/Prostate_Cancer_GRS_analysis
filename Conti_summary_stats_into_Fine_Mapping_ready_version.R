## Converting Conti summary stats MASTER/AFR excel file into rs72725854 Fine Mapping-ready version

install.packages('readxl')
library(readxl)
library(dplyr)

system(paste("dx download", "Callum/ContiGWAS/ContiGWASsummaryStatsMASTER.xlsx")) ## This is Conti's Multiancestry GWAS summary stats as downloadable here: https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs001120/analyses/phs001120.pha005082.txt, and converted to Susie-compatible format
system(paste("dx download", "Callum/ContiGWAS/ContiGWASsummaryStatsAFR.xlsx")) ## This is Conti's AFR-specific GWAS summary stats as downloadable here: https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs001120/analyses/phs001120.pha005078.txt, and converted to Susie-compatible format

ContiGWASsummaryStatsMaster <- read_excel("ContiGWASsummaryStatsMASTER.xlsx")
ContiGWASsummaryStatsAFR <- read_excel("ContiGWASsummaryStatsAFR.xlsx")

## Step 1: isolate only Chromosome 8 entries

ContiGWASsummaryStatsChr8 <- ContiGWASsummaryStatsAFR %>%    ## Choose which summary stats version (Master or AFR) you want to use here
  dplyr::filter(CHROM == "8") %>%
  dplyr::rename("rsid" = "ID")

## Step 2: Remove all insertions and deletions

ContiGWASsummaryStatsChr8 <- ContiGWASsummaryStatsChr8 %>%
  dplyr::filter(
    (ALLELE1 == "T" | ALLELE1 == "A" | ALLELE1 == "G" | ALLELE1 == "C") &
      (ALLELE0 == "T" | ALLELE0 == "A" | ALLELE0 == "G" | ALLELE0 == "C")
  )

## Step 3: Create new ID column that is CHR:POS:ALT:REF and goes both ways (i.e., for each SNP, do CHR:POS:ALLELE1:ALLELE0 and CHR:POS:ALLELE0:ALLELE1)

library(stringr)

# Orientation 1: as-is
fwd <- ContiGWASsummaryStatsChr8 %>%
  transmute(
    CHROM, GENPOS, rsid,
    ALLELE1, ALLELE0,
    BETA, SE, P,
    ID = str_c(CHROM, GENPOS, ALLELE1, ALLELE0, sep = ":"),
    orientation = "A1_as_effect"
  )

# Orientation 2: reversed (swap alleles using temps, flip BETA)
rev <- ContiGWASsummaryStatsChr8 %>%
  transmute(
    CHROM, GENPOS, rsid,
    ALLELE1_rev = ALLELE0,   # use temps to avoid masking
    ALLELE0_rev = ALLELE1,
    BETA = -BETA,            # flip sign in reversed orientation
    SE, P
  ) %>%
  rename(ALLELE1 = ALLELE1_rev, ALLELE0 = ALLELE0_rev) %>%
  mutate(
    ID = str_c(CHROM, GENPOS, ALLELE1, ALLELE0, sep = ":"),
    orientation = "A0_as_effect"
  )

ContiGWASsummaryStatsChr8_long <- bind_rows(fwd, rev) %>%
  arrange(CHROM, GENPOS, ID, desc(orientation))

ContiGWASsummaryStatsChr8_noindel_clean <- ContiGWASsummaryStatsChr8_long %>%
  dplyr::select("CHROM", "GENPOS", "ID", "ALLELE1", "ALLELE0", "BETA", "SE", "P") %>%
  dplyr::mutate(
    OR = exp(BETA),
    OR_L95 = exp(BETA - 1.96 * SE),
    OR_U95 = exp(BETA + 1.96 * SE)
  )

# Write to tsv

write.table(ContiGWASsummaryStatsChr8_noindel_clean,'ContiGWASsummaryStatsChr8_noindel_clean.tsv',quote=FALSE,sep='\t',row.names = FALSE)

# Write to gzipped tsv
write.table(
  ContiGWASsummaryStatsChr8_noindel_clean,
  file = gzfile("ContiGWASsummaryStatsChr8_noindel_clean.tsv.gz"),
  quote = FALSE, sep = "\t", row.names = FALSE
)


## Upload to project

system(paste("dx upload", "ContiGWASsummaryStatsChr8_noindel_clean.tsv.gz"))   ### Primary file used in fine mapping
system(paste("dx upload", "ContiGWASsummaryStatsChr8_noindel_clean.tsv"))