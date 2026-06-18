#############################################################################
# Create adjustedGRS by weighting Conti odds ratios by Ancestry Probability #
#############################################################################


#--------------------------------------------------------------------#
# Setup (run these installs separately, then the rest of the script) #
#--------------------------------------------------------------------#

remotes::install_github("lcpilling/ukbrapR@v0.3.10",
                        force = TRUE, clean = TRUE, dependencies = TRUE)


library(dplyr)
library(readxl)
library(stringr)
library(readr)
library(tidyr)

# -----------------------------
# Inputs
# -----------------------------

system('dx download Conti.xlsx') # This is the supplementary table file from Conti et al. (2021), downloadable at: https://www.nature.com/articles/s41588-020-00748-0
system("dx download Callum/Derived_datasets/AncestryProbability2.csv") # Genetic similarity probabilities for Ancestry group calculated using AncestryProbability2_plink folder script

AP <- read.csv("AncestryProbability2.csv") %>%
  dplyr::select(IID, AFR, AMR, CSA, EAS, EUR, MID) %>%
  rename(eid = IID)

# -----------------------------
# 1) Load Conti S4 and extract OR columns
# -----------------------------

Conti_raw <- read_excel("Conti.xlsx", sheet = "S4", skip = 3, na = "NA")
Conti_raw <- Conti_raw[1:269, ] %>% arrange(Chromosome, Position)

get_or_col <- function(df, pattern) {
  nm <- grep(pattern, names(df), value = TRUE)
  if (length(nm) == 0) stop(paste("Could not find OR column matching:", pattern))
  as.numeric(df[[nm[1]]])
}

Conti_base0 <- Conti_raw %>%
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
  # Keep rows with valid positive ORs across all needed ancestries
  #filter(OR_MULTI > 0, OR_EUR > 0, OR_AFR > 0, OR_EAS > 0, OR_HIS > 0) %>%
  mutate(CHR = ifelse(CHR %in% c("X","x"), "X", as.character(as.integer(CHR))))

Conti_base <- Conti_base0 %>%
  # Ensure CHR formatting 
  mutate(CHR = ifelse(CHR %in% c("X","x"), "X", as.character(as.integer(CHR)))) %>%
  # Replace invalid/missing ancestry ORs with Multiethnic
  mutate(
    OR_EUR = ifelse(is.finite(OR_EUR) & OR_EUR > 0, OR_EUR, OR_MULTI),
    OR_AFR = ifelse(is.finite(OR_AFR) & OR_AFR > 0, OR_AFR, OR_MULTI),
    OR_EAS = ifelse(is.finite(OR_EAS) & OR_EAS > 0, OR_EAS, OR_MULTI),
    OR_HIS = ifelse(is.finite(OR_HIS) & OR_HIS > 0, OR_HIS, OR_MULTI)
  )

# This varlist is just to make a BED; beta column not needed for bed creation
varlist_for_bed <- Conti_base %>%
  transmute(
    rsid = rsID,
    chr  = CHR,
    pos  = POS,
    effect_allele,
    other_allele,
    beta = 0  # placeholder; we won't use plink scoring
  )

# -----------------------------
# 3) Build a BED with only these SNPs (imputed POS, build 37)
# -----------------------------
ukbrapR::make_imputed_bed(
  in_file   = varlist_for_bed,
  out_bed   = "Conti_subset",
  use_pos   = TRUE,
  progress  = FALSE,
  verbose   = TRUE,
  very_verbose = TRUE
)


stopifnot(file.exists("Conti_subset.bed"))

## Create .raw subset file


plink2 <- "/home/rstudio-server/_ukbrapr_tools/plink2"
plink1 <- "/home/rstudio-server/_ukbrapr_tools/plink"



cmd <- NULL
if (file.exists(plink2)) {
  cmd <- sprintf('%s --bfile %s --export A --out %s',
                 shQuote(plink2), shQuote("Conti_subset"), shQuote("Conti_subset"))
} else if (file.exists(plink1)) {
  cmd <- sprintf('%s --bfile %s --recode A --out %s',
                 shQuote(plink1), shQuote("Conti_subset"), shQuote("Conti_subset"))
} else {
  stop("PLINK tools not found at expected paths. Try ukbrapR:::prep_tools() again.")
}


cmd_with_log <- sprintf('%s 2>&1 | tee plink_export.log', cmd)
status <- system(cmd_with_log)

if (status != 0L) {
  stop("PLINK command failed (exit code ", status, "). See plink_export.log.")
}
if (!file.exists("Conti_subset.raw")) {
  stop("PLINK finished but Conti_subset.raw was not created.\n",
       "Inspect plink_export.log; check that the .bim contains >0 variants.")
}


stopifnot(file.exists("Conti_subset.raw"))


# -----------------------------
# 4) Read BIM to align alleles, and RAW dosages (no ID renaming)
# -----------------------------

bim <- readr::read_tsv("Conti_subset.bim",
                       col_names = c("chr","id","null","pos","a1","a2"),
                       show_col_types = FALSE)

raw <- read.table("Conti_subset.raw", header = TRUE, check.names = FALSE)

# Normalise chromosome for joins (map "X"->23 if present)
normalize_chr_to_int <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x[x %in% c("X", "x")] <- "23"
  suppressWarnings(as.integer(x))
}

# Build Conti SNP key; join on chr:pos to bring ORs/alleles alongside BIM
Conti_key <- Conti_base %>%
  dplyr::transmute(
    chr_join = normalize_chr_to_int(CHR),
    pos_join = as.integer(POS),
    effect_allele, other_allele,
    OR_MULTI, OR_EUR, OR_AFR, OR_EAS, OR_HIS
  ) %>%
  dplyr::distinct(chr_join, pos_join, .keep_all = TRUE)

bim2 <- bim %>%
  dplyr::mutate(
    chr_join = normalize_chr_to_int(chr),
    pos_join = as.integer(pos)
  )

# Keep only variants present in both BIM and Conti by chr:pos
join1 <- dplyr::inner_join(bim2, Conti_key, by = c("chr_join","pos_join"))

#############################################
# Important - remove withdrawn participants #
#############################################

exclude_withdrawn=function(df){
  system('dx download Callum/Withdrawals/withdrawn_20260310.csv --overwrite') ## This file is a list of participants who withdrew from the Biobank up to the date 10th March 2026. This was sent from the UK Biobank team via email to members of approved applications
  df2 = df %>% left_join(
    read_csv("withdrawn_20260310.csv", col_names = FALSE, show_col_types = FALSE) %>%
      mutate(w=1) %>%
      rename(IID = X1),
    by='IID'
  ) %>%
    filter(is.na(w))
  return(df2)
}

raw <- exclude_withdrawn(raw)

# -----------------------------
# Robust prefix-based mapping from .raw columns to BIM rows
# -----------------------------

# 0) Parse .raw headers: split at LAST underscore -> base + counted allele
geno_cols <- setdiff(names(raw), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))
raw_parsed <- tibble::tibble(
  raw_col = geno_cols,
  base    = sub("_([^_]*)$", "", geno_cols),      # everything before last underscore
  counted = sub("^.*_", "", geno_cols)            # token after last underscore
)

# 1) Build the "<chr>:<pos>" prefix for join1 rows
pos_prefix <- paste0(join1$chr_join, ":", join1$pos_join)

# 2) For each variant, find raw columns whose 'base' starts with the <chr>:<pos> prefix
raw_bases <- raw_parsed$base
raw_counted <- raw_parsed$counted

find_match_index <- function(k) {
  # candidates starting with "<chr>:<pos>"
  hits <- which(startsWith(raw_bases, pos_prefix[k]))
  # if none, fallback to candidates starting with exact BIM id
  if (length(hits) == 0) {
    id_pref <- paste0(join1$id[k])
    hits <- which(startsWith(raw_bases, id_pref))
  }
  if (length(hits) == 0) return(NA_integer_)
  
  # Prefer hits whose COUNTED allele equals either Conti effect or other allele (case-insensitive)
  alle_ok <- toupper(raw_counted[hits]) %in% toupper(c(join1$effect_allele[k], join1$other_allele[k]))
  if (any(alle_ok, na.rm = TRUE)) {
    hits <- hits[which(alle_ok)]
  }
  # If multiple remain, take the first
  hits[1]
}

match_idx <- vapply(seq_len(nrow(join1)), find_match_index, integer(1))

# 3) Keep matched rows only
keep <- !is.na(match_idx)
if (!any(keep)) {
  message("No .raw columns matched by position or BIM id. Let's print a few clues:")
  message("Example expected prefixes: ", paste(unique(pos_prefix[1:min(5, length(pos_prefix))]), collapse = ", "))
  message("First 10 RAW bases: ", paste(head(raw_bases, 10), collapse = " | "))
  stop("Mapping failed; share the two lines above so we can extend the matcher.")
}

panel <- dplyr::bind_cols(
  join1[keep, , drop = FALSE],
  raw_parsed[match_idx[keep], , drop = FALSE]
)

message("Mapped ", nrow(panel), " of ", nrow(join1),
        " variants to .raw columns (", nrow(join1) - nrow(panel), " missing).")

# 4) Subset RAW columns in this exact order
raw_sub <- raw %>% dplyr::select(IID, dplyr::all_of(panel$raw_col))

# -----------------------------
# 5. Convert dosages to be coded on the Conti effect allele
# -----------------------------
dosage_mat <- as.matrix(raw_sub[,-1]); mode(dosage_mat) <- "numeric"
stopifnot(ncol(dosage_mat) == nrow(panel))

# Flip based on the .raw COUNTED allele vs Conti effect/other allele
counted_up <- toupper(panel$counted)
eff_up     <- toupper(panel$effect_allele)
oth_up     <- toupper(panel$other_allele)

keep_ok <- counted_up %in% eff_up | counted_up %in% oth_up
if (!all(keep_ok)) {
  message("Dropping ", sum(!keep_ok), " variants where counted allele in RAW isn't one of the Conti alleles.")
  panel <- panel[keep_ok, , drop = FALSE]
  dosage_mat <- dosage_mat[, keep_ok, drop = FALSE]
  counted_up <- counted_up[keep_ok]; eff_up <- eff_up[keep_ok]; oth_up <- oth_up[keep_ok]
}

# If counted allele equals the Conti 'other', flip to get dosage of effect allele
flip_to_eff <- counted_up == oth_up
if (any(flip_to_eff, na.rm = TRUE)) {
  dosage_mat[, flip_to_eff] <- 2 - dosage_mat[, flip_to_eff]
}

# -----------------------------
# 6) Build OR matrices and compute per-participant adjusted betas
# -----------------------------
# Your mapping: EUR*European + AFR*African + EAS*East Asian + AMR*Hispanic
#              + CSA*Multiethnic + MID*African
# (Feel free to change CSA/MID mappings here.)

# Participant weights (no renormalisation; use exactly your columns)
W <- AP %>%
  select(eid, AFR, AMR, CSA, EAS, EUR, MID)

# Align participants to RAW order; drop any without probabilities
W <- W %>% filter(eid %in% raw$IID)

raw_ids <- raw$IID

W <- W %>% slice(match(raw_ids, eid))
stopifnot(all(W$eid == raw_ids))

# 6 x m OR "reference" matrix (rows: EUR, AFR, EAS, AMR, CSA->MULTI, MID->AFR)
OR_ref <- rbind(
  EUR = panel$OR_EUR,
  AFR = panel$OR_AFR,
  EAS = panel$OR_EAS,
  AMR = panel$OR_HIS,
  CSA = panel$OR_MULTI,
  MID = panel$OR_AFR
)  # dimensions: 6 x m

# n x 6 participant weights
W_mat <- as.matrix(W %>% select(EUR, AFR, EAS, AMR, CSA, MID))
mode(W_mat) <- "numeric"

# Compute per-participant, per-SNP OR_mix: (n x 6) %*% (6 x m) = (n x m)
OR_mix <- W_mat %*% OR_ref

# Safety: enforce positivity
if (any(OR_mix <= 0, na.rm = TRUE)) {
  stop("Found non-positive OR after mixing; check inputs.")
}

# Betas are log of OR_mix
beta_mat <- log(OR_mix)

# -----------------------------
# 7) Compute GRS (sum and average-per-allele)
# -----------------------------
# Handle missing dosages: keep NA, we’ll compute denominators per person
non_missing_counts <- rowSums(!is.na(dosage_mat))

# Sum score: sum_j G * beta
# (Use rowSums with NA handling by replacing NA with 0, but then
# the denominator will reflect how many SNPs were non-missing.)
tmp_G <- dosage_mat
tmp_G[is.na(tmp_G)] <- 0

sum_scores <- rowSums(tmp_G * beta_mat)

# Denominator akin to "alleles counted" (2 per non-missing SNP)
alleles_counted <- 2 * non_missing_counts
avg_scores <- ifelse(alleles_counted > 0, sum_scores / alleles_counted, NA_real_)

# -----------------------------
# 8) Save outputs
# -----------------------------
out <- tibble(
  eid = raw_ids,
  Conti_ORfirst_sum = sum_scores,
  Conti_ORfirst_avg = avg_scores
)



write.table(out, "OR_adjustedGRS267.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Optional upload
system("dx upload OR_adjustedGRS267.tsv")
