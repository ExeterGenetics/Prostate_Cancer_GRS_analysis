###########################################################################################
#=========================================================================================#
# Converting Wang summary stats AFR excel file into rs72725854 Fine Mapping-ready version #
#=========================================================================================#
###########################################################################################

library(readxl)
library(dplyr)

system(paste("dx download", "Callum/WangGWAS/Wang2023African_harmonised.tsv")) ## This is Wang's AFR-specific GWAS summary stats, haramonised to GRCh38, as downloadable here: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/GCST90274715/. Go to /harmonised and download GCST90274715.h.tsv.gz
system(paste("dx download", "Callum/WangGWAS/Wang2023Asian_harmonised.tsv")) ## This is Wang's Asian-specific GWAS summary stats, haramonised to GRCh38, as downloadable here: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/GCST90274716/. Go to /harmonised and download GCST90274716.h.tsv.gz

WangGWASsummaryStats <- read.delim("Wang2023African_harmonised.tsv", sep = "\t") %>% # Change "African" to "Asian" in the filename to load the Asian-specific summary stats instead
  dplyr::rename(
    "CHROM" = "chromosome",
    "GENPOS" = "base_pair_location",
    "rsid" = "rsid",
    "ALLELE1" = "effect_allele",
    "ALLELE0" = "other_allele",
    "BETA" = "beta",
    "SE" = "standard_error",
    "P" = "p_value",
  )

###############################################
## Step 1: isolate only Chromosome 8 entries ##
###############################################

WangGWASsummaryStatsChr8 <- WangGWASsummaryStats %>%     
  dplyr::filter(CHROM == "8")


##################################################
### Step 2: Perform liftover (GRCh37 -> GRCh38) ## (commented out as no longer required - summary stats are already harmonised)
##################################################
#
## Install/load required Bioconductor packages
#
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}
#for (pkg in c("rtracklayer", "GenomicRanges")) {
#  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
#  library(pkg, character.only = TRUE)
#}
#
### Download UCSC chain file if not already present
#
#system(paste("dx download", "Callum/tools/hg19ToHg38.over.chain")) ## This chain file can be downloaded at https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz 
#
#
#chain_file <- "hg19ToHg38.over.chain"
#
#chain <- rtracklayer::import.chain(chain_file)
#
### Build GRanges from current (GRCh37/hg19) coordinates
#
#gr_hg19 <- GenomicRanges::GRanges(
#  seqnames = paste0("chr", WangGWASsummaryStatsChr8$CHROM),
#  ranges   = IRanges::IRanges(start = WangGWASsummaryStatsChr8$GENPOS,
#                              end   = WangGWASsummaryStatsChr8$GENPOS)
#)
#
### Run liftover
#
#lv <- rtracklayer::liftOver(gr_hg19, chain)
#
### Keep only uniquely mapped sites (exactly one target interval)
#
#is_unique <- elementNROWS(lv) == 1L
#gr_hg38   <- unlist(lv[is_unique])
#
### Indices for bookkeeping
#
#idx_unique  <- which(is_unique)
#idx_drop    <- which(elementNROWS(lv) == 0L)              # unmapped
#idx_multi   <- which(elementNROWS(lv)  > 1L)              # multi-mapped
#
### Optional: write reports so you can see what was lost/duplicated
#
#if (length(idx_drop)) {
#  write.table(WangGWASsummaryStatsChr8[idx_drop, ],
#              file = "liftover_chr8_unmapped_GRCh37.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#}
#if (length(idx_multi)) {
#  write.table(WangGWASsummaryStatsChr8[idx_multi, ],
#              file = "liftover_chr8_multimapped_GRCh37.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#}
#
### Subset your data to the uniquely lifted variants and update coordinates to GRCh38
#
#WangGWASsummaryStatsChr8 <- WangGWASsummaryStatsChr8[idx_unique, ]
#WangGWASsummaryStatsChr8$CHROM  <- as.character(GenomeInfoDb::seqnames(gr_hg38))
#WangGWASsummaryStatsChr8$GENPOS <- GenomicRanges::start(gr_hg38)
#
### Convert back from UCSC style ("chr8") to plain "8" to match the rest of your pipeline
#
#WangGWASsummaryStatsChr8$CHROM <- sub("^chr", "", WangGWASsummaryStatsChr8$CHROM)
#
###############################################################################
### Step 2b: Harmonize alleles to GRCh38 (ensure ALLELE0 == REF, ALLELE1 == ALT)
###############################################################################
#
### This section ensures allele/reference harmonization on GRCh38:
###  - Ensures ALLELE0 equals the GRCh38 reference base at GENPOS and ALLELE1 is the alt.
###  - Swaps alleles + flip BETA when reference sits on ALLELE1.
###  - Handles complement (strand) scenarios.
###  - Drops unresolved palindromic or mismatching variants and writes a report.
#
#for (pkg in c("BSgenome", "BSgenome.Hsapiens.UCSC.hg38", "Biostrings")) {
#  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
#  library(pkg, character.only = TRUE)
#}
#
#ref <- BSgenome.Hsapiens.UCSC.hg38
#
### Upper-case alleles for safety
#WangGWASsummaryStatsChr8$ALLELE0 <- toupper(WangGWASsummaryStatsChr8$ALLELE0)
#WangGWASsummaryStatsChr8$ALLELE1 <- toupper(WangGWASsummaryStatsChr8$ALLELE1)
#
### Pull GRCh38 reference base at each position
#ref_base <- as.character(
#  Biostrings::getSeq(
#    ref,
#    paste0("chr", WangGWASsummaryStatsChr8$CHROM),
#    WangGWASsummaryStatsChr8$GENPOS,
#    WangGWASsummaryStatsChr8$GENPOS
#  )
#)
#ref_base <- toupper(ref_base)
#
### Convenience functions/flags
#comp <- function(a) chartr("ACGT", "TGCA", a)
#A0  <- WangGWASsummaryStatsChr8$ALLELE0
#A1  <- WangGWASsummaryStatsChr8$ALLELE1
#is_pal <- (A0 %in% c("A","T") & A1 %in% c("A","T")) |
#  (A0 %in% c("C","G") & A1 %in% c("C","G"))
#
#match_ref_A0      <- (A0 == ref_base)
#match_ref_A1      <- (A1 == ref_base)
#match_ref_A0_comp <- (comp(A0) == ref_base)
#match_ref_A1_comp <- (comp(A1) == ref_base)
#
#swap      <- rep(FALSE, nrow(WangGWASsummaryStatsChr8))
#flip_beta <- rep(FALSE, nrow(WangGWASsummaryStatsChr8))
#drop      <- rep(FALSE, nrow(WangGWASsummaryStatsChr8))
#note      <- rep(NA_character_, nrow(WangGWASsummaryStatsChr8))
#
## Case 1: ALLELE0 already equals GRCh38 REF
#keep1 <- match_ref_A0
#note[keep1] <- "aligned_ref"
#
### Case 2: ALLELE1 equals GRCh38 REF -> swap + flip
#keep2 <- !keep1 & match_ref_A1
#swap[keep2]      <- TRUE
#flip_beta[keep2] <- TRUE
#note[keep2]      <- "swap_flip_ref_on_A1"
#
### Case 3: Complement scenario, ALLELE0 complement equals GRCh38 REF
#keep3 <- !(keep1 | keep2) & match_ref_A0_comp
#A0[keep3]   <- comp(A0[keep3])
#A1[keep3]   <- comp(A1[keep3])
#note[keep3] <- "complement_A0"
#
### Case 4: Complement + swap + flip (ALLELE1 complement equals GRCh38 REF)
#keep4 <- !(keep1 | keep2 | keep3) & match_ref_A1_comp
#A0[keep4]      <- comp(A0[keep4])
#A1[keep4]      <- comp(A1[keep4])
#swap[keep4]    <- TRUE
#flip_beta[keep4] <- TRUE
#note[keep4]    <- "complement_then_swap_flip"
#
### Unresolved: either true mismatch or palindromic ambiguity -> drop
#unresolved <- !(keep1 | keep2 | keep3 | keep4)
#drop[unresolved] <- TRUE
#note[unresolved] <- ifelse(is_pal[unresolved], "drop_palindromic_unresolved", "drop_ref_mismatch")
#
### Apply swap
#A0_final <- ifelse(swap, A1, A0)
#A1_final <- ifelse(swap, A0, A1)
#
#WangGWASsummaryStatsChr8$ALLELE0 <- A0_final
#WangGWASsummaryStatsChr8$ALLELE1 <- A1_final
#
### Flip beta where needed
#WangGWASsummaryStatsChr8$BETA[flip_beta] <- -WangGWASsummaryStatsChr8$BETA[flip_beta]
#
### Attach reference and notes
#WangGWASsummaryStatsChr8$REF38               <- ref_base
#WangGWASsummaryStatsChr8$harmonization_note  <- note
#
### Drop unresolved and write a report
#harmonized <- WangGWASsummaryStatsChr8[!drop, ]
#dropped    <- WangGWASsummaryStatsChr8[drop, ]
#if (nrow(dropped) > 0) {
#  write.table(dropped, "chr8_hg38_harmonization_dropped.tsv",
#              sep = "\t", quote = FALSE, row.names = FALSE)
#}
#WangGWASsummaryStatsChr8 <- harmonized

#################################################
## Step 3: Remove all insertions and deletions ##
#################################################

WangGWASsummaryStatsChr8 <- WangGWASsummaryStatsChr8 %>%
  dplyr::filter(
    (ALLELE1 == "T" | ALLELE1 == "A" | ALLELE1 == "G" | ALLELE1 == "C") &
      (ALLELE0 == "T" | ALLELE0 == "A" | ALLELE0 == "G" | ALLELE0 == "C")
  )

#################################################
## Step 4: Create new ID column that is        ##
## CHR:POS:ALT:REF and goes both ways          ##
## (i.e., for each SNP, do                     ##
## CHR:POS:ALLELE1:ALLELE0 and                 ##
## CHR:POS:ALLELE0:ALLELE1)                    ##
################################################# 

library(stringr)

# Orientation 1: as-is
fwd <- WangGWASsummaryStatsChr8 %>%
  dplyr::transmute(
    CHROM, GENPOS, rsid,
    ALLELE1, ALLELE0,
    BETA, SE, P,
    ID = str_c(CHROM, GENPOS, ALLELE1, ALLELE0, sep = ":"),
    orientation = "A1_as_effect"
  )

# Orientation 2: reversed (swap alleles using temps, flip BETA)
rev <- WangGWASsummaryStatsChr8 %>%
  dplyr::transmute(
    CHROM, GENPOS, rsid,
    ALLELE1_rev = ALLELE0,   # use temps to avoid masking
    ALLELE0_rev = ALLELE1,
    BETA = -BETA,            # flip sign in reversed orientation
    SE, P
  ) %>%
  dplyr::rename(ALLELE1 = ALLELE1_rev, ALLELE0 = ALLELE0_rev) %>%
  dplyr::mutate(
    ID = str_c(CHROM, GENPOS, ALLELE1, ALLELE0, sep = ":"),
    orientation = "A0_as_effect"
  )

WangGWASsummaryStatsChr8_long <- bind_rows(fwd, rev) %>%
  arrange(CHROM, GENPOS, ID, desc(orientation))

WangGWASsummaryStatsChr8_noindel_clean <- WangGWASsummaryStatsChr8_long %>%
  dplyr::select("CHROM", "GENPOS", "ID", "ALLELE1", "ALLELE0", "BETA", "SE", "P") %>%
  dplyr::mutate(
    OR = exp(BETA),
    OR_L95 = exp(BETA - 1.96 * SE),
    OR_U95 = exp(BETA + 1.96 * SE)
  )

# Write to tsv

write.table(WangGWASsummaryStatsChr8_noindel_clean,'WangGWASsummaryStatsChr8_noindel_clean.tsv',quote=FALSE,sep='\t',row.names = FALSE)

# Write to gzipped tsv
write.table(
  WangGWASsummaryStatsChr8_noindel_clean,
  file = gzfile("WangGWASsummaryStatsChr8_noindel_clean.tsv.gz"),
  quote = FALSE, sep = "\t", row.names = FALSE
)


## Upload to project

system(paste("dx upload", "WangGWASsummaryStatsChr8_noindel_clean.tsv.gz"))   ### Primary file used in fine mapping
system(paste("dx upload", "WangGWASsummaryStatsChr8_noindel_clean.tsv"))
