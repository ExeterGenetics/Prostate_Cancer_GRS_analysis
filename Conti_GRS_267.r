##########################
# Compute Conti/Wang GRS # 
##########################

## Note: this is how the following dataset files in "Callum/Derived_datasets" are made: 
## - Conti_multiethnicGRS_267.tsv
## - Conti_EuropeanGRS_265.tsv, 
## - Conti_AfricanGRS_246.tsv, 
## - Conti_East_AsianGRS_222.tsv, 
## - Conti_HispanicGRS_253.tsv 
## - Wang_multiethnicGRS_450.tsv
## - Wang_EuropeanGRS_445.tsv
## - Wang_AfricanGRS_444.tsv
## - Wang_East_AsianGRS_379.tsv
## - Wang_HispanicGRS_446.tsv

install.packages("remotes")
remotes::install_github("lcpilling/ukbrapR@v0.3.10",
                        force = TRUE, clean = TRUE, dependencies = TRUE)
install.packages("readxl")

library(dplyr)
library(readxl)
library(readr)
library(stringr)
library(tidyr)
library(tibble)

## Choose which OR column to use

# Options: "OR_MULTI", "OR_EUR", "OR_AFR", "OR_EAS", "OR_HIS"
OR_column <- "OR_MULTI"

## Choose the source GWAS: "Conti" or "Wang" (note: if changing source, restart the R session to clear old data as bgenix reuses temporary files between runs)
source <- "Wang" 


##########
# Inputs #
##########

system('dx download Callum/ContiGWAS/Conti2021supplementarytables.xlsx') # This is the supplementary table file from Conti et al. (2021), downloadable at: https://www.nature.com/articles/s41588-020-00748-0
system('dx download Callum/WangGWAS/Wang2023supplementarytables.xlsx') # This is the supplementary table file from Wang et al. (2023), downloadable at: https://pmc.ncbi.nlm.nih.gov/articles/PMC10841479/

#########################################################################
# Step 1: Load sheet 4 from Supplementary tables and extract OR columns #
#########################################################################

if (source == "Conti") {
  raw_data <- read_excel("Conti2021supplementarytables.xlsx", sheet = "S4", skip = 3, na = "NA")
  raw_data <- raw_data[1:269, ] %>% arrange(Chromosome, Position)
  
  get_or_col <- function(df, pattern) {
    nm <- grep(pattern, names(df), value = TRUE)
    if (length(nm) == 0) stop(paste("Could not find OR column matching:", pattern))
    as.numeric(df[[nm[1]]])
  }
  
  base_data <- raw_data %>%
    rename(
      rsID = `rs*`,
      CHR  = Chromosome,
      POS  = Position,
      effect_allele = `Risk Allele`,
      other_allele  = `Reference Allele`
    ) %>%
    mutate(
      OR_MULTI = get_or_col(., "^Multiethnic"),
      OR_EUR   = get_or_col(., "^European...16"),
      OR_AFR   = get_or_col(., "^African...19"),
      OR_EAS   = get_or_col(., "^East\\s*Asian...22"),
      OR_HIS   = get_or_col(., "^Hispanic...25")
    ) %>%
    mutate(CHR = ifelse(CHR %in% c("X","x"), "X", as.character(as.integer(CHR))))
  
} else if (source == "Wang") {
  raw_data <- read_excel("Wang2023supplementarytables.xlsx", sheet = "S4", skip = 3, na = "NA")
  
  # Rename columns as in Conti_script.R
  raw_data <- raw_data %>%
    dplyr::rename(
      EUR_Rsquared = European,
      AFR_Rsquared = African,
      EAS_Rsquared = Asian,          
      HIS_Rsquared = Hispanic,       
      OR_Multiethnic_Marginal = `OR...15`,
      CI95_Multiethnic_Marginal = `95%CI...16`,
      Pvalue_Multiethnic_Marginal = `P-value...17`,
      OR_Multiethnic_Conditional = `OR...18`,
      CI95_Multiethnic_Conditional = `95%CI...19`,
      Pvalue_Multiethnic_Conditional = `P-value...20`,
      RAF_EUR = `RAF...21`,
      OR_EUR = `OR...22`,
      CI95_EUR = `95%CI...23`,
      Pvalue_EUR = `P-value...24`,
      RAF_AFR = `RAF...25`,
      OR_AFR = `OR...26`,
      CI95_AFR = `95%CI...27`,
      Pvalue_AFR = `P-value...28`,
      RAF_EAS = `RAF...29`,
      OR_EAS = `OR...30`,
      CI95_EAS = `95%CI...31`,
      Pvalue_EAS = `P-value...32`,
      RAF_HIS = `RAF...33`,
      OR_HIS = `OR...34`,
      CI95_HIS = `95%CI...35`,
      Pvalue_HIS = `P-value...36`
    ) %>%
    # Clean P-values
    mutate(
      across(
        c(Pvalue_Multiethnic_Marginal, Pvalue_Multiethnic_Conditional,
          Pvalue_EUR, Pvalue_AFR, Pvalue_EAS, Pvalue_HIS),
        ~ .x |> as.character() |> str_replace("^\\s*<\\s*", "") |> as.numeric()
      ),
      across(
        c(EUR_Rsquared, AFR_Rsquared, EAS_Rsquared, HIS_Rsquared,
          OR_Multiethnic_Marginal, OR_Multiethnic_Conditional,
          RAF_EUR, OR_EUR, RAF_AFR, OR_AFR, RAF_EAS, OR_EAS, RAF_HIS, OR_HIS),
        ~ suppressWarnings(as.numeric(.x))
      )
    )
  
  # Filter to Wang 451-SNP GRS-specific variants
  raw_data <- raw_data %>%
    dplyr::filter(!is.na(`Wang et al., 451 SNPs`))
  
  base_data <- raw_data %>%
    rename(
      rsID = rsID,
      CHR  = Chromosome,
      POS  = `Position (GRCh37)`,
      effect_allele = `Risk Allele`,
      other_allele  = `Reference Allele`
    ) %>%
    mutate(
      OR_MULTI = OR_Multiethnic_Marginal,
      # OR_EUR, OR_AFR, OR_EAS, OR_HIS already present from Wang data
      CHR = ifelse(CHR %in% c("X","x"), "X", as.character(as.integer(CHR)))
    )
  
} else {
  stop("Invalid source. Choose 'Conti' or 'Wang'.")
}

#######################################################
# Step 2: Extract genotypes for these SNPs -> BED set #
#######################################################

# This varlist is just to make a BED; beta is a placeholder

varlist_for_bed <- base_data %>%
  transmute(
    rsid = rsID,
    chr  = CHR,
    pos  = POS,
    effect_allele,
    other_allele,
    beta = 0
  )

ukbrapR::make_imputed_bed(
  in_file     = varlist_for_bed,
  out_bed     = "Conti_subset",
  use_pos     = TRUE,    # use build-37 positions
  progress    = FALSE,
  verbose     = TRUE,
  very_verbose = TRUE
)

stopifnot(file.exists("Conti_subset.bed"),
          file.exists("Conti_subset.bim"),
          file.exists("Conti_subset.fam"))

###########################################
# Step 3: Export additive dosages to .raw #
###########################################

plink2 <- "/home/rstudio-server/_ukbrapr_tools/plink2"
plink1 <- "/home/rstudio-server/_ukbrapr_tools/plink"

if (file.exists(plink2)) {
  cmd <- sprintf('%s --bfile %s --export A --out %s',
                 shQuote(plink2), shQuote("Conti_subset"), shQuote("Conti_subset"))
} else if (file.exists(plink1)) {
  cmd <- sprintf('%s --bfile %s --recode A --out %s',
                 shQuote(plink1), shQuote("Conti_subset"), shQuote("Conti_subset"))
} else {
  stop("PLINK binaries not found. Run ukbrapR:::prep_tools() or adjust paths.")
}

status <- system(sprintf('%s 2>&1 | tee plink_export.log', cmd))
if (status != 0L) stop("PLINK failed (exit ", status, "). See plink_export.log.")
stopifnot(file.exists("Conti_subset.raw"))

############################################################
# Step 4: Read BIM + RAW and align to Conti effect alleles #
############################################################

bim <- readr::read_tsv("Conti_subset.bim",
                       col_names = c("chr","id","null","pos","a1","a2"),
                       show_col_types = FALSE)

raw <- read.table("Conti_subset.raw", header = TRUE, check.names = FALSE)

normalize_chr_to_int <- function(x) {
  x <- as.character(x); x <- trimws(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("X","x")] <- "23"
  suppressWarnings(as.integer(x))
}

Conti_key <- base_data %>%
  transmute(
    chr_join = normalize_chr_to_int(CHR),
    pos_join = as.integer(POS),
    effect_allele, other_allele,
    OR_MULTI, OR_EUR, OR_AFR, OR_EAS, OR_HIS
  ) %>%
  distinct(chr_join, pos_join, .keep_all = TRUE)

bim2 <- bim %>%
  mutate(chr_join = normalize_chr_to_int(chr),
         pos_join = as.integer(pos))

join1 <- inner_join(bim2, Conti_key, by = c("chr_join","pos_join"))

# Parse .raw headers
geno_cols  <- setdiff(names(raw), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))
raw_parsed <- tibble(
  raw_col = geno_cols,
  base    = sub("_([^_]*)$", "", geno_cols),
  counted = sub("^.*_", "", geno_cols)
)

pos_prefix  <- paste0(join1$chr_join, ":", join1$pos_join)
raw_bases   <- raw_parsed$base
raw_counted <- raw_parsed$counted

find_match_index <- function(k) {
  hits <- which(startsWith(raw_bases, pos_prefix[k]))
  if (length(hits) == 0) {
    id_pref <- paste0(join1$id[k])
    hits <- which(startsWith(raw_bases, id_pref))
  }
  if (length(hits) == 0) return(NA_integer_)
  alle_ok <- toupper(raw_counted[hits]) %in%
    toupper(c(join1$effect_allele[k], join1$other_allele[k]))
  if (any(alle_ok, na.rm = TRUE)) hits <- hits[which(alle_ok)]
  hits[1]
}

match_idx <- vapply(seq_len(nrow(join1)), find_match_index, integer(1))
keep <- !is.na(match_idx)
if (!any(keep)) {
  stop("Could not match any .raw columns; check headers and BIM IDs.")
}

panel <- dplyr::bind_cols(
  join1[keep, , drop = FALSE],
  raw_parsed[match_idx[keep], , drop = FALSE]
)

raw_sub   <- raw %>% dplyr::select(IID, dplyr::all_of(panel$raw_col))
dosage_mat <- as.matrix(raw_sub[,-1]); mode(dosage_mat) <- "numeric"
stopifnot(ncol(dosage_mat) == nrow(panel))

# Flip to Conti effect allele
counted_up <- toupper(panel$counted)
eff_up     <- toupper(panel$effect_allele)
oth_up     <- toupper(panel$other_allele)

keep_ok <- counted_up %in% eff_up | counted_up %in% oth_up
if (!all(keep_ok)) {
  message("Dropping ", sum(!keep_ok),
          " variants where counted allele isn't one of the Conti alleles.")
  panel      <- panel[keep_ok, , drop = FALSE]
  dosage_mat <- dosage_mat[, keep_ok, drop = FALSE]
  counted_up <- counted_up[keep_ok]; eff_up <- eff_up[keep_ok]; oth_up <- oth_up[keep_ok]
}

flip_to_eff <- counted_up == oth_up
if (any(flip_to_eff, na.rm = TRUE)) {
  dosage_mat[, flip_to_eff] <- 2 - dosage_mat[, flip_to_eff]
}

###############################################################
# Step 5: Build betas from a single OR column and compute GRS #
###############################################################

if (!OR_column %in% c("OR_MULTI","OR_EUR","OR_AFR","OR_EAS","OR_HIS")) {
  stop("OR_column must be one of: OR_MULTI, OR_EUR, OR_AFR, OR_EAS, OR_HIS")
}

# 1) Take the chosen ORs and keep only valid, positive values
or_vec <- panel[[OR_column]]
ok_or  <- is.finite(or_vec) & (or_vec > 0)

# If none are usable, write NAs and exit gracefully
if (!any(ok_or)) {
  warning("No usable SNPs (all OR missing/non-positive) for ", OR_column,
          ". Writing NA scores.")
  raw_ids <- raw_sub$IID
  sum_scores <- rep(NA_real_, length(raw_ids))
  avg_scores <- rep(NA_real_, length(raw_ids))
  tag       <- gsub("^OR_", "", OR_column)
  col_sum   <- paste0("Conti_GRS_", tag, "_sum")
  col_avg   <- paste0("Conti_GRS_", tag, "_avg")
  out_file  <- paste0("GRS_", tag, ".tsv")
  out <- tibble(eid = raw_ids, sum = sum_scores, avg = avg_scores)
  names(out)[2:3] <- c(col_sum, col_avg)
  write.table(out, out_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  message("Finished (no usable SNPs). Wrote: ", out_file)
  quit(save = "no")
}

# 2) Subset to SNPs with valid ORs (keep panel & dosage columns in sync)
if (sum(!ok_or) > 0) {
  message("Dropping ", sum(!ok_or), " SNPs with NA/non-positive ORs in ", OR_column, ".")
}
panel      <- panel[ok_or, , drop = FALSE]
dosage_mat <- dosage_mat[, ok_or, drop = FALSE]
beta_vec   <- log(or_vec[ok_or])

# 3) Handle missing dosages — keep NA in denominator, 0 for the sum
non_missing_counts <- rowSums(!is.na(dosage_mat))
tmp_G <- dosage_mat
tmp_G[is.na(tmp_G)] <- 0

# 4) If zero columns remain (unlikely after ok_or, but safe-guard anyway)
if (ncol(tmp_G) == 0) {
  sum_scores <- rep(NA_real_, nrow(tmp_G))
  avg_scores <- rep(NA_real_, nrow(tmp_G))
} else {
  sum_scores <- as.vector(tmp_G %*% beta_vec)       # sum_j G * beta
  alleles_counted <- 2 * non_missing_counts         # per-participant denominator
  avg_scores <- ifelse(alleles_counted > 0, sum_scores / alleles_counted, NA_real_)
}

########################
# Step 6: Save outputs #
########################

raw_ids   <- raw_sub$IID
tag       <- gsub("^OR_", "", OR_column)  # e.g., MULTI, EUR, AFR...
col_sum   <- paste0("Conti_GRS_", tag, "_sum")
col_avg   <- paste0("Conti_GRS_", tag, "_avg")
out_file  <- paste0("GRS_", tag, ".tsv")

out <- tibble(
  eid = raw_ids,
  sum = sum_scores,
  avg = avg_scores
)
names(out)[2:3] <- c(col_sum, col_avg)

#############################################
# Important - remove withdrawn participants #
#############################################

exclude_withdrawn=function(df){
  system('dx download Callum/Withdrawals/withdrawn_20260310.csv --overwrite') ## This file is a list of participants who withdrew from the Biobank up to the date 10th March 2026. This was sent from the UK Biobank team via email to members of approved applications
  df2 = df %>% left_join(
    read_csv("withdrawn_20260310.csv", col_names = FALSE, show_col_types = FALSE) %>%
      mutate(w=1) %>%
      rename(eid = X1),
    by='eid'
  ) %>%
    filter(is.na(w))
  return(df2)
}

out <- exclude_withdrawn(out)

write.table(out, out_file,
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)




# Optional: upload to RAP project

system(paste("dx upload", shQuote(out_file)))

message("Finished. Wrote: ", out_file)

