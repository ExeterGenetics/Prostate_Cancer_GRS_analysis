##########################
# Compute Conti/Wang GRS # 
##########################

## Note: this is how the following dataset files in "Callum/GRSs" are made: 
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
## - SchumacherGRS_145.tsv
## - BARCODE1GRS_129.tsv
## - SeibertGRS_52.tsv
## - PagadalaGRS_285.tsv

remotes::install_github("lcpilling/ukbrapR@v0.3.10",
                        force = TRUE, clean = TRUE, dependencies = TRUE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

library(readxl)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)

## Choose which OR column to use, then you can run the rest of the script to generate your chosen GRS

# Options: "OR_MULTI", "OR_EUR", "OR_AFR", "OR_EAS", "OR_HIS"
OR_column <- "OR_MULTI"

## Choose the source GWAS: "Conti", "Wang", "Schumacher", "BARCODE1", "Seibert", or "Pagadala"
source <- "Conti" 


##########
# Inputs #
##########

## GRSs trained to predict general prostate cancer diagnosis

system('dx download Callum/ContiGWAS/Conti2021supplementarytables.xlsx') # This is the supplementary table file from Conti et al. (2021), downloadable at: https://www.nature.com/articles/s41588-020-00748-0
system('dx download Callum/WangGWAS/Wang2023supplementarytables.xlsx') # This is the supplementary table file from Wang et al. (2023), downloadable at: https://pmc.ncbi.nlm.nih.gov/articles/PMC10841479/
system('dx download Callum/SchumacherGWAS/Schumacher.txt') # This is the list of SNPs and weights from Schumacher et al. (2018), downloadable at: https://www.pgscatalog.org/publication/PGP000019/ 
system('dx download Callum/SchumacherGWAS/BARCODE1.txt') # This is the list of SNPs and weights from BARCODE1 (2021), downloadable at: https://www.pgscatalog.org/publication/PGP000726/ 

## GRSs trained to predict aggressive prostate cancer diagnosis

system('dx download Callum/SeibertGWAS/Seibert.txt') # This is the list of SNPs and weights from Seibert et al. (2018), downloadable at: https://www.pgscatalog.org/publication/PGP000047/ 
system('dx download Callum/SeibertGWAS/Pagadala.txt') # This is the list of SNPs and weights from Pagadala et al. (2022), downloadable at: https://www.pgscatalog.org/publication/PGP000400/

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
      OR_HIS   = get_or_col(., "^Hispanic...25"),
      effect_weight = NA_real_
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
      CHR = ifelse(CHR %in% c("X","x"), "X", as.character(as.integer(CHR))),
      effect_weight = NA_real_
    )
  
} else if (source == "Schumacher") {
  raw_data <- read.delim("Schumacher.txt", comment.char = "#")
  
  base_data <- raw_data %>%
    arrange(chr_name, chr_position) %>%
    rename(
      rsID = rsID,
      CHR  = chr_name,
      POS  = chr_position,
      effect_allele = effect_allele,
      effect_weight = effect_weight
    ) %>%
    mutate(
      other_allele = NA_character_,
      OR_MULTI = NA_real_,
      OR_EUR   = NA_real_,
      OR_AFR   = NA_real_,
      OR_EAS   = NA_real_,
      OR_HIS   = NA_real_,
      CHR = ifelse(CHR %in% c("X","x"), "X", as.character(as.integer(CHR)))
    )
  
} else if (source == "BARCODE1" | source == "Seibert" | source == "Pagadala") {
  raw_data <- read.delim(paste0(source, ".txt"), comment.char = "#")
  
  base_data <- raw_data %>%
    arrange(chr_name, chr_position) %>%
    rename(
      rsID = rsID,
      CHR  = chr_name,
      POS  = chr_position,
      effect_allele = effect_allele,
      other_allele = other_allele,
      effect_weight = effect_weight
    ) %>%
    mutate(
      OR_MULTI = NA_real_,
      OR_EUR   = NA_real_,
      OR_AFR   = NA_real_,
      OR_EAS   = NA_real_,
      OR_HIS   = NA_real_,
      CHR = ifelse(CHR %in% c("X","x"), "X", as.character(as.integer(CHR)))
    )
  
} else {
  stop("Invalid source. Choose 'Conti', 'Wang', 'Schumacher', 'BARCODE1', 'Seibert', or 'Pagadala'.")
}

#######################################################
# Step 2: Extract genotypes for these SNPs -> BED set #
#######################################################

# This varlist is just to make a BED; beta is a placeholder

if (source %in% c("Conti", "Wang")) {
  run_tag <- paste0(source, "_", gsub("^OR_", "", OR_column))
} else {
  run_tag <- source
}
out_prefix <- paste0("Conti_subset_", run_tag)
raw_export_file <- paste0(out_prefix, ".raw")
bim_file <- paste0(out_prefix, ".bim")

# Remove stale artifacts from a prior run with the same tag.
unlink(
  c(
    paste0(out_prefix, ".bed"),
    paste0(out_prefix, ".bim"),
    paste0(out_prefix, ".fam"),
    paste0(out_prefix, ".log"),
    paste0(out_prefix, ".nosex"),
    raw_export_file,
    "plink_export.log"
  ),
  force = TRUE
)

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
  out_bed     = out_prefix,
  use_pos     = TRUE,    # use build-37 positions
  progress    = FALSE,
  verbose     = TRUE,
  very_verbose = TRUE
)

stopifnot(file.exists(paste0(out_prefix, ".bed")),
          file.exists(paste0(out_prefix, ".bim")),
          file.exists(paste0(out_prefix, ".fam")))

###########################################
# Step 3: Export additive dosages to .raw #
###########################################

plink2 <- "/home/rstudio-server/_ukbrapr_tools/plink2"
plink1 <- "/home/rstudio-server/_ukbrapr_tools/plink"

if (file.exists(plink2)) {
  cmd <- sprintf('%s --bfile %s --export A --out %s',
                 shQuote(plink2), shQuote(out_prefix), shQuote(out_prefix))
} else if (file.exists(plink1)) {
  cmd <- sprintf('%s --bfile %s --recode A --out %s',
                 shQuote(plink1), shQuote(out_prefix), shQuote(out_prefix))
} else {
  stop("PLINK binaries not found. Run ukbrapR:::prep_tools() or adjust paths.")
}

status <- system(sprintf('%s 2>&1 | tee plink_export.log', cmd))
if (status != 0L) stop("PLINK failed (exit ", status, "). See plink_export.log.")
stopifnot(file.exists(raw_export_file))

############################################################
# Step 4: Read BIM + RAW and align to Conti effect alleles #
############################################################

bim <- readr::read_tsv(bim_file,
                       col_names = c("chr","id","null","pos","a1","a2"),
                       show_col_types = FALSE)

raw <- read.table(raw_export_file, header = TRUE, check.names = FALSE)

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
    OR_MULTI, OR_EUR, OR_AFR, OR_EAS, OR_HIS,
    effect_weight
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

########################################################################
# Step 5: Build SNP weights (OR-derived or direct effect weight) and GRS #
########################################################################

if (source %in% c("Conti", "Wang")) {
  if (!OR_column %in% c("OR_MULTI","OR_EUR","OR_AFR","OR_EAS","OR_HIS")) {
    stop("OR_column must be one of: OR_MULTI, OR_EUR, OR_AFR, OR_EAS, OR_HIS")
  }
  
  source_vec <- panel[[OR_column]]
  ok_weight  <- is.finite(source_vec) & (source_vec > 0)
  weight_desc <- paste0("OR column ", OR_column)
} else {
  source_vec <- panel$effect_weight
  ok_weight  <- is.finite(source_vec)
  weight_desc <- paste0(source, " effect_weight")
}

# If none are usable, write NAs and exit gracefully
if (!any(ok_weight)) {
  warning("No usable SNP weights for ", weight_desc, ". Writing NA scores.")
  raw_ids <- raw_sub$IID
  sum_scores <- rep(NA_real_, length(raw_ids))
  avg_scores <- rep(NA_real_, length(raw_ids))
  
  if (source %in% c("Conti", "Wang")) {
    tag       <- gsub("^OR_", "", OR_column)
    col_sum   <- paste0("Conti_GRS_", tag, "_sum")
    col_avg   <- paste0("Conti_GRS_", tag, "_avg")
    out_file  <- paste0("GRS_", tag, ".tsv")
  } else {
    tag       <- source
    col_sum   <- paste0(source, "_GRS_sum")
    col_avg   <- paste0(source, "_GRS_avg")
    out_file  <- paste0("GRS_", source, ".tsv")
  }
  
  out <- tibble(eid = raw_ids, sum = sum_scores, avg = avg_scores)
  names(out)[2:3] <- c(col_sum, col_avg)
  write.table(out, out_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  message("Finished (no usable SNP weights). Wrote: ", out_file)
  quit(save = "no")
}

# Subset to SNPs with valid weights (keep panel & dosage columns in sync)
if (sum(!ok_weight) > 0) {
  message("Dropping ", sum(!ok_weight), " SNPs with unusable weights in ", weight_desc, ".")
}
panel      <- panel[ok_weight, , drop = FALSE]
dosage_mat <- dosage_mat[, ok_weight, drop = FALSE]

if (source %in% c("Conti", "Wang")) {
  beta_vec <- log(source_vec[ok_weight])
} else {
  beta_vec <- source_vec[ok_weight]
}

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

raw_ids <- raw_sub$IID

if (source %in% c("Conti", "Wang")) {
  tag       <- gsub("^OR_", "", OR_column)   # MULTI / EUR / AFR / EAS / HIS
  col_sum   <- paste0(source, "_GRS_", tag, "_sum")
  col_avg   <- paste0(source, "_GRS_", tag, "_avg")
  out_file  <- paste0("GRS_", source, "_", tag, ".tsv")
} else {
  tag       <- source
  col_sum   <- paste0(source, "_GRS_sum")
  col_avg   <- paste0(source, "_GRS_avg")
  out_file  <- paste0("GRS_", source, ".tsv")
}

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

