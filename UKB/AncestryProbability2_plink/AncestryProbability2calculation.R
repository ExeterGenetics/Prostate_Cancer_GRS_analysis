#===============================================================#
# Part 1 — Fetch HGDP+1KG projection files and load in UKB files #
# ===============================================================#

# Session & helpers

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  library(data.table)
})

library(data.table)

# Check PLINK availability (we need PLINK 1.9+ here)
plink1 <- Sys.which("plink")
if (!nzchar(plink1)) {
  message("PLINK not found in PATH; downloading a local copy via bigsnpr...")
  if (!requireNamespace("bigsnpr", quietly = TRUE)) install.packages("bigsnpr")
  plink1_path <- bigsnpr::download_plink(dir = "~")  # installs plink (1.x)
  Sys.chmod(plink1_path, "0755")
  Sys.setenv(PATH = paste(dirname(plink1_path), Sys.getenv("PATH"), sep = ":"))
  plink1 <- Sys.which("plink"); stopifnot(nzchar(plink1))
}

# Fetch HGDP+1KG projection assets into the current working directory
# These Principal Components loadings and scores are publicly available and can be downloaded from:
# https://console.cloud.google.com/storage/browser/covid19-hg-public/pca_projection;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

dxdownload("Callum/HGDP_1KG/hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.tsv")  
dxdownload("Callum/HGDP_1KG/hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.afreq")
dxdownload("Callum/HGDP_1KG/hgdp_tgp_pca_covid19hgi_snps_scores.txt")                 
dxdownload("Callum/HGDP_1KG/hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.tsv")

loadings_tsv <- "hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.tsv"
afreq_file   <- "hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.afreq"
scores_ref   <- "hgdp_tgp_pca_covid19hgi_snps_scores.txt"


rsid_loadings <- "hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.tsv"

# Make a one-column RSID list

system(sprintf("awk 'NR>1{print $1}' %s > variants.rsid", rsid_loadings))

dxdownload("Callum/tools/plink2"); Sys.chmod("plink2", "0755")  # can be downloaded from https://www.cog-genomics.org/plink/2.0/
plink2_bin <- normalizePath("./plink2")



#===========================================================================#
# Part 2A — Build subsetted PGENs with SNPs used to calculate reference PCs #
# (bgenix subset → index → pgen)                                            #
#===========================================================================#

stopifnot(file.exists(loadings_tsv), file.exists(afreq_file), file.exists("variants.rsid"))

bgen_path <- '/mnt/project/Bulk/Imputation/UKB imputation from genotype'




### Load in bgenix 

system("dx download 'Callum/tools/bgen.tgz' -o bgen.tgz") # can be downloaded from http://code.enkre.net/bgen/tarball/release/bgen.tgz

## Install/Expose bgenix + cat-bgen from a .tgz (binary or source)


# Clean staging folders
unlink(c("tools/bgenix_src", "tools/bgenix_local"), recursive = TRUE, force = TRUE)
dir.create("tools/bgenix_src", recursive = TRUE, showWarnings = FALSE)
dir.create("tools/bgenix_local", recursive = TRUE, showWarnings = FALSE)
local_prefix <- normalizePath("tools/bgenix_local", mustWork = TRUE)

stopifnot(file.exists("bgen.tgz"))

# Use base R untar to peek at contents and to extract

contents <- utils::untar("bgen.tgz", list = TRUE)
stopifnot(length(contents) > 0)

# Determine top-level directory name inside the tarball

top_dirs <- unique(sub("/.*$", "", contents))
top_dir  <- top_dirs[1]

# Extract into a known location

utils::untar("bgen.tgz", exdir = "tools/bgenix_src")
src_dir <- normalizePath(file.path("tools/bgenix_src", top_dir), mustWork = TRUE)

# Case 1: Prebuilt binaries are present (bin/bgenix or a lone `bgenix`)

bin_candidates <- c(
  file.path(src_dir, "bin", "bgenix"),
  file.path(src_dir, "bgenix")
)
bin_exists <- bin_candidates[file.exists(bin_candidates)]

if (length(bin_exists)) {
  # Treat as prebuilt release
  bgenix_bin <- bin_exists[1]
  cat_bgen_bin <- file.path(dirname(bgenix_bin), "cat-bgen")
  
  # Install into a local prefix
  dir.create(file.path(local_prefix, "bin"), recursive = TRUE, showWarnings = FALSE)
  ok1 <- file.copy(bgenix_bin, file.path(local_prefix, "bin", "bgenix"), overwrite = TRUE)
  if (file.exists(cat_bgen_bin)) {
    ok2 <- file.copy(cat_bgen_bin, file.path(local_prefix, "bin", "cat-bgen"), overwrite = TRUE)
  } else {
    ok2 <- TRUE
  }
  Sys.chmod(file.path(local_prefix, "bin", c("bgenix", "cat-bgen")), "0755")
  
} else {
  # 3) Case 2: Source tree — build with waf (Python 3)
  # Sanity checks for expected files
  waf <- file.path(src_dir, "waf")
  view_cpp <- file.path(src_dir, "src", "View.cpp")
  stopifnot(file.exists(waf) || file.exists(file.path(src_dir, "wscript")))
  stopifnot(file.exists(view_cpp))  # indicates proper source tree
  
  # Optional patch (older GCC): std::ios::streampos -> std::streampos
  suppressWarnings(
    system(sprintf("grep -n 'std::ios::streampos' %s || true", shQuote(view_cpp)))
  )
  invisible(system(sprintf("sed -i 's/std::ios::streampos/std::streampos/g' %s", shQuote(view_cpp))))
  
  # Build & install
  py <- Sys.which("python3"); if (!nzchar(py)) stop("python3 not found on PATH.")
  cmd <- sprintf(
    "bash -lc 'cd %s && %s waf configure --prefix=%s && %s waf && %s waf install'",
    shQuote(src_dir), shQuote(py), shQuote(local_prefix), shQuote(py), shQuote(py)
  )
  status <- system(cmd)
  if (status != 0L) stop("waf build failed — scroll up in the log for the first error.")
}

# Expose to PATH and verify
Sys.setenv(PATH = paste(file.path(local_prefix, "bin"), Sys.getenv("PATH"), sep = ":"))
message("bgenix on PATH? ", Sys.which("bgenix"))
system("bgenix -help")   # should print usage
if (!nzchar(Sys.which("bgenix"))) stop("bgenix is not available on PATH after install.")

# Optional: record exact tool versions for reproducibility
system("bgenix -help | head -n 3")
system("which bgenix; ls -l $(which bgenix)")


# Helper: choose a writable work directory with an actual write test
choose_writable_dir <- function(candidates) {
  # return first path for which we can create+delete a file
  for (d in candidates) {
    if (!nzchar(d)) next
    ok <- try(dir.create(d, recursive = TRUE, showWarnings = FALSE), silent = TRUE)
    if (!dir.exists(d)) next
    tf <- file.path(d, sprintf(".writetest_%s", as.integer(runif(1, 1e6, 9e6))))
    can <- tryCatch({ file.create(tf) }, warning = function(e) FALSE, error = function(e) FALSE)
    if (isTRUE(can)) { file.remove(tf); return(normalizePath(d, mustWork = TRUE)) }
  }
  return(NA_character_)
}

# Prefer project subfolders if writable; otherwise fall back to CWD/home/tempdir.
job_tag <- Sys.getenv("DX_JOB_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
work_dir <- choose_writable_dir(c(
  file.path("/mnt/project", "work", paste0("pca_subset_", job_tag)),
  file.path("/mnt/project", paste0("pca_subset_", job_tag)),
  file.path(getwd(),       paste0("pca_subset_", job_tag)),
  file.path("~",           paste0("pca_subset_", job_tag)),
  tempdir()  # last resort
))
if (is.na(work_dir)) {
  stop("Could not find a writable work directory. Check project permissions.")
}
message("Using work_dir: ", work_dir)

# For downstream parts, reuse 'work_dir' consistently instead of Sys.getenv("TMPDIR")

tmp_dir <- work_dir

message("Using tmp_dir: ", tmp_dir)

# Ensure tmp_dir exists and is writable; if not, fall back to a fresh temp path

dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

.testfile <- file.path(tmp_dir, sprintf(".writetest_%s", as.integer(Sys.time())))
.canwrite <- tryCatch(isTRUE(file.create(.testfile)),
                      warning = function(e) FALSE, error = function(e) FALSE)
if (.canwrite) unlink(.testfile, force = TRUE) else {
  job_tag <- Sys.getenv("DX_JOB_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
  tmp_dir <- file.path(tempdir(), paste0("pca_subset_", job_tag))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  .testfile <- file.path(tmp_dir, sprintf(".writetest_%s", as.integer(Sys.time())))
  .canwrite <- tryCatch(isTRUE(file.create(.testfile)),
                        warning = function(e) FALSE, error = function(e) FALSE)
  if (.canwrite) unlink(.testfile, force = TRUE) else {
    stop("No writable tmp_dir available. Check filesystem permissions.")
  }
}
message("Final tmp_dir (writable): ", tmp_dir)

# Define per‑loop helpers here so they’re in scope

stopifnot(file.exists("variants.rsid"))
extract_file    <- normalizePath("variants.rsid")
subset_prefixes <- character(22)

##

stopifnot(file.exists("variants.rsid"))
extract_file    <- normalizePath("variants.rsid")  # absolute path, one RSID per line
subset_prefixes <- character(22)                    # will hold per‑chr output prefixes


threads <- "16"
mem_gb <- "64"

for (chr in 1:22) {
  src_bgen   <- sprintf("%s/ukb22828_c%d_b0_v3.bgen",  bgen_path, chr)
  src_sample <- sprintf("%s/ukb22828_c%d_b0_v3.sample", bgen_path, chr)
  stopifnot(file.exists(src_bgen), file.exists(src_sample))
  
  # Copy the tiny .sample next to our temp outputs (space-free name)
  local_sample <- src_sample
  stopifnot(file.exists(local_sample))
  
  # Subset BGEN in tmp_dir (same as before)
  subset_bgen <- file.path(tmp_dir, sprintf("subset_chr%d.bgen", chr))
  if (file.exists(subset_bgen)) file.remove(subset_bgen)
  
  cmd_subset <- sprintf("bgenix -g %s -incl-rsids %s > %s",
                        shQuote(src_bgen), shQuote(extract_file), shQuote(subset_bgen))
  message(sprintf("[chr %d] bgenix subset → %s", chr, subset_bgen))
  st1 <- system(cmd_subset)
  if (st1 != 0L || !file.exists(subset_bgen) || file.size(subset_bgen) == 0) {
    stop(sprintf("bgenix subset failed for chr %d; check RSIDs/index.", chr))
  }
  
  message(sprintf("[chr %d] bgenix index (clobber)", chr))
  st2 <- system2("bgenix", c("-g", subset_bgen, "-index", "-clobber"))
  if (st2 != 0L || !file.exists(paste0(subset_bgen, ".bgi"))) {
    stop(sprintf("bgenix -index failed for chr %d.", chr))
  }
  
  # Convert the small BGEN to PGEN (use the staged .sample)
  pfx <- file.path(tmp_dir, sprintf("ukb_subset_chr%d", chr))
  subset_prefixes[chr] <- pfx
  
  args_make_pgen <- c(
    "--threads", threads,
    "--memory",  as.character(as.integer(mem_gb) * 1024),
    "--bgen",    subset_bgen, "ref-first",
    "--sample",  local_sample,
    "--make-pgen",
    "--out",     pfx
  )
  
  message(sprintf("[chr %d] subset.bgen → PGEN: %s", chr, pfx))
  out_txt <- system2(plink2_bin, args_make_pgen, stdout = TRUE, stderr = TRUE)
  
  if (!file.exists(paste0(pfx, ".pgen"))) {
    logf <- paste0(pfx, ".log")
    if (file.exists(logf)) {
      cat("---- PLINK log tail ----\n",
          paste(tail(readLines(logf), 50), collapse = "\n"), "\n", sep = "")
    }
    stop(sprintf("PLINK --make-pgen failed for chr %d.", chr))
  }
}

# Optional: you can delete subset_bgen + .bgi now to reclaim space
# unlink(c(subset_bgen, paste0(subset_bgen, ".bgi")), force = TRUE)

system("df -h")
system(sprintf("df -h %s", shQuote(work_dir)))

# ===================================== #
# Part 2B — Merge to a single `.pfile`  #
# ===================================== #

merge_list <- file.path(tmp_dir, "pmerge_list.txt")

# Which per-chr pfiles exist?

existing <- which(file.exists(paste0(subset_prefixes, ".pgen")))
if (!length(existing)) stop("No per-chromosome .pgen files were produced. Check Part 2A logs.")

# Choose the first existing prefix as the merge base

base_idx    <- existing[1]
base_prefix <- subset_prefixes[base_idx]

# Everything else (excluding the base) goes into --pmerge-list

merge_targets <- subset_prefixes[setdiff(existing, base_idx)]
writeLines(merge_targets, con = merge_list)

merged_pfx <- file.path(tmp_dir, "ukb_subset_merged")

args_merge <- c(
  "--threads", threads,
  "--memory",  as.character(as.integer(mem_gb) * 1024),
  "--pfile",   base_prefix,        # <- now guaranteed to exist
  "--pmerge-list", merge_list,
  "--out",     merged_pfx
)

message("Merging per-chromosome subset PGENs → genome-wide .pfile")
system2(plink2_bin, args_merge, stdout = TRUE, stderr = TRUE)

# This `--pmerge-list` route is the documented way to merge PLINK2 pfiles. [2](https://zenodo.org/records/15420125)

## Sanity check - how many variants carried over from HGDP+1KG loadings to UKB ?

# Overlap audit between loadings and merged .pfile (

library(data.table)

merged_pfx <- file.path("ukb_subset_merged")

audit_overlap_all <- function(loadings_tsv, rsid_loadings, merged_pfx, plink2_bin, tmp_dir) {
  stopifnot(file.exists(paste0(merged_pfx, ".pgen")),
            file.exists(paste0(merged_pfx, ".pvar")),
            file.exists(paste0(merged_pfx, ".psam")))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # RSID-based sanity check 
  
  rsid_vec <- fread(rsid_loadings, nThread = 1)[[1]]
  rsid_vec <- unique(rsid_vec[!is.na(rsid_vec) & rsid_vec != "rsid" & !grepl("^#", rsid_vec)])
  
  pvar <- fread(paste0(merged_pfx, ".pvar"), nThread = 1, data.table = TRUE)
  # usual columns: '#CHROM','POS','ID','REF','ALT'
  stopifnot("ID" %in% names(pvar))
  pvar_ids <- unique(pvar$ID)
  
  rsid_overlap <- length(intersect(rsid_vec, pvar_ids))
  
  # Position-keyed checks with different patterns ---
  
  pos_ids <- unique(fread(loadings_tsv, nThread = 1)[[1]][-1])  # skip header
  pos_ids <- pos_ids[!is.na(pos_ids)]
  writeLines(pos_ids, file.path(tmp_dir, "loadings.posid"))
  
  mk_snplist_and_count <- function(pattern, tag) {
    out <- file.path(tmp_dir, paste0("merged_posid_", tag))
    st <- system2(plink2_bin, c("--pfile", merged_pfx,
                                "--set-all-var-ids", pattern,
                                "--write-snplist",
                                "--out", out),
                  stdout = TRUE, stderr = TRUE)
    snpfile <- paste0(out, ".snplist")
    if (!file.exists(snpfile)) return(list(tag=tag, n=NA_integer_, file=NA_character_))
    merged_ids <- unique(readLines(snpfile, warn = FALSE))
    n <- length(intersect(pos_ids, merged_ids))
    list(tag=tag, n=n, file=snpfile)
  }
  
  # Try four patterns: no 'chr'/with 'chr' × REF/ALT and ALT/REF
  tries <- list(
    mk_snplist_and_count("@:#:$r:$a", "nochr_ra"),
    mk_snplist_and_count("@:#:$a:$r", "nochr_ar"),
    mk_snplist_and_count("chr@:#:$r:$a", "chr_ra"),
    mk_snplist_and_count("chr@:#:$a:$r", "chr_ar")
  )
  
  # Pick the best positional pattern
  best <- tries[[which.max(sapply(tries, `[[`, "n"))]]
  
  cat("==== Overlap audit ====\n")
  cat(sprintf("RSID overlap (loadings.rsid vs merged IDs): %d\n", rsid_overlap))
  cat(sprintf("PosID overlap (best of 4 patterns): %s = %s\n",
              best$tag, ifelse(is.na(best$n), "NA", best$n)))
  
  invisible(list(
    rsid_overlap = rsid_overlap,
    pos_overlap  = best$n,
    pos_pattern  = best$tag
  ))
}

# Run
audit <- audit_overlap_all(loadings_tsv, rsid_loadings, merged_pfx, plink2_bin, tmp_dir)



## Upload to project for future reference (will only work if you have made these during the session, otherwise just dx download them in Part C)

system(paste("dx upload", "pca_subset_job-J6Bqxz0J6jP38G4ZBF9bjYqK/ukb_subset_merged.pgen"))
system(paste("dx upload", "pca_subset_job-J6Bqxz0J6jP38G4ZBF9bjYqK/ukb_subset_merged.psam"))
system(paste("dx upload", "pca_subset_job-J6Bqxz0J6jP38G4ZBF9bjYqK/ukb_subset_merged.pvar"))
system(paste("dx upload", "pca_subset_job-J6Bqxz0J6jP38G4ZBF9bjYqK/ukb_subset_merged.log"))

# ================================================== #
# Part 2C — Single score pass from the merged .pfile #  From here onwards, this is self-contained
# ================================================== #

## Load in merged UKB subset from part 2B

dxdownload("Callum/HGDP_1KG/ukb_subset_merged.pgen")
dxdownload("Callum/HGDP_1KG/ukb_subset_merged.psam")
dxdownload("Callum/HGDP_1KG/ukb_subset_merged.pvar")
dxdownload("Callum/HGDP_1KG/ukb_subset_merged.log")

# Ensure RSID AF file is present
if (!file.exists("hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.afreq")) {
  dxdownload("Callum/HGDP_1KG/hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.afreq")
}
afreq_rsid    <- "hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.afreq"
loadings_rsid <- rsid_loadings  # "hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.tsv"

merged_pfx <- file.path("ukb_subset_merged")

score_out <- file.path("UKB_on_HGDP1KG_merged_rsid")

threads <- "16"
mem_gb <- "64"

args_score <- c(
  "--threads", threads,
  "--memory",  as.character(as.integer(mem_gb) * 1024),
  "--pfile",   merged_pfx,
  
  "--score",   loadings_rsid, "variance-standardize",
  "cols=-scoreavgs,+scoresums", "header-read", "list-variants",
  "--score-col-nums", "3-22",
  "--read-freq", afreq_rsid,
  "--out",       score_out
)


message("Scoring PCs 1–20 from merged .pfile (single pass)")
system2(plink2_bin, args_score, stdout = TRUE, stderr = TRUE)  # linear scoring behavior per PLINK2 docs [4](https://pan.ukbb.broadinstitute.org/)


## Scale to HGDP+1KG PCs

library(data.table)

ss <- fread(paste0(score_out, ".sscore"))           
nv <- nrow(fread(paste0(score_out, ".sscore.vars"),  
                 header = FALSE))                  

pc_sum <- intersect(paste0("PC", 1:20, "_SUM"), names(ss))
stopifnot(length(pc_sum) >= 2)

ukb_scaled <- ss[, .(IID)]
for (i in seq_along(pc_sum)) {
  nm <- sub("_SUM$", "", pc_sum[i])                 # "PC1".."PC20"
  ukb_scaled[[nm]] <- ss[[pc_sum[i]]] / sqrt(nv)    # scale by sqrt(#variants)
}

#############################################
# Important - remove withdrawn participants #
#############################################

library(dplyr)

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

ukb_scaled <- exclude_withdrawn(ukb_scaled)

fwrite(ukb_scaled, "HGDP_1KG_PCs_UKB2_scaled.csv")




# =================================================== #
# Part 2D — Parse the single .sscore → UKB PCs (CSV)  #
# =================================================== #

sscore <- paste0(score_out, ".sscore")
stopifnot(file.exists(sscore))
sc <- data.table::fread(paste0(score_out, ".sscore"))

# With cols=-scoreavgs,+scoresums we should have PC1_SUM..PC20_SUM present. [4](https://pan.ukbb.broadinstitute.org/)

sum_cols <- intersect(paste0("PC", 1:20, "_SUM"), names(sc))
if (!length(sum_cols)) {
  # fallback: SCORE1_SUM..SCORE20_SUM
  sum_cols <- intersect(paste0("SCORE", 1:20, "_SUM"), names(sc))
}

keep     <- c("IID", intersect(sum_cols, names(sc)))
ukb_dt   <- sc[, ..keep]
data.table::setnames(ukb_dt, old = intersect(sum_cols, names(sc)),
                     new = sub("_SUM$", "", intersect(sum_cols, names(sc))))

#############################################
# Important - remove withdrawn participants #
#############################################

ukb_dt <- exclude_withdrawn(ukb_dt)

data.table::fwrite(ukb_dt, "HGDP_1KG_PCs_UKB2_unscaled.csv")  

# Reference PCs unchanged

ref_scores <- data.table::fread(scores_ref)
data.table::setnames(ref_scores, "s", "IID", skip_absent = TRUE)
ref_out <- data.frame(IID = ref_scores$IID, ref_scores[, paste0("PC", 1:20), with = FALSE])
data.table::fwrite(ref_out, "HGDP_1KG_PCs2.csv")

## Save to project

system(paste("dx upload", "HGDP_1KG_PCs2.csv"))
system(paste("dx upload", "HGDP_1KG_PCs_UKB2_unscaled.csv"))
system(paste("dx upload", "HGDP_1KG_PCs_UKB2_scaled.csv"))




# You may keep or later delete the merged pfile if you plan to reuse it for other tasks.



#=====================================================#
# Part 3 - Use PCs 1–6 and scikit‑learn Random Forest #
#=====================================================#

library(data.table); library(dplyr)

# Load PCs (from Part 2) + metadata
HGDP_1KG_PCs <- read.csv("HGDP_1KG_PCs2.csv", check.names = FALSE)
UKB_PCs      <- read.csv("HGDP_1KG_PCs_UKB2_scaled.csv", check.names = FALSE)

dxdownload("Callum/HGDP_1KG/release_3.1_secondary_analyses_hgdp_1kg_v2_metadata_and_qc_gnomad_meta_updated.tsv") # This can be downloaded from https://console.cloud.google.com/storage/browser/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc

meta <- fread("release_3.1_secondary_analyses_hgdp_1kg_v2_metadata_and_qc_gnomad_meta_updated.tsv")
lab  <- meta[, .(IID = s, Region = `hgdp_tgp_meta.Genetic.region`)]

HGDP_1KG_PCs_w_labels <- merge(HGDP_1KG_PCs, lab, by = "IID") %>%
  dplyr::select(c("IID", "PC1", "PC2", "PC3", "PC4", "PC5" , "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "Region"))

# Merge labels, keep PC1..PC6, drop OCE

keepPCs <- paste0("PC", 1:6)
train <- HGDP_1KG_PCs %>%
  left_join(lab, by="IID") %>%
  dplyr::filter(!is.na(Region), Region != "OCE") %>%
  dplyr::select(c("IID", all_of(keepPCs), "Region"))

ukb_pred <- UKB_PCs %>%
  dplyr::filter(!is.na(PC1)) %>%
  dplyr::select(c("IID", all_of(keepPCs)))

# Use scikit-learn via reticulate to reproduce Pan-UKB helper defaults

if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
library(reticulate)
py_install(packages = c("scikit-learn==1.4.2"), pip = TRUE)

sk <- import("sklearn", delay_load = FALSE)
RFC <- import("sklearn.ensemble")$RandomForestClassifier

X_train <- as.matrix(train[, keepPCs])
y_train <- as.character(train$Region)

set.seed(42)
rf <- RFC(n_estimators = as.integer(10000), max_features = "sqrt", random_state = as.integer(42))
rf$fit(X_train, y_train)

# Predict probability for each class on UKB
X_test <- as.matrix(ukb_pred[, keepPCs])
proba  <- rf$predict_proba(X_test)
classes <- rf$classes_

probDF <- data.frame(IID = ukb_pred$IID, proba)
colnames(probDF) <- c("IID", as.character(classes))

# Reorder to canonical column order if present
canon <- c("AMR","AFR","CSA","EAS","EUR","MID")
probDF <- probDF[, c("IID", intersect(canon, colnames(probDF)))]

# Add most probable ancestry, margin, and entropy columns

prob_cols <- c(
  "AMR",
  "AFR",
  "CSA",
  "EAS",
  "EUR",
  "MID"
)

eps <- 1e-12
P <- as.matrix(probDF[, prob_cols])

max_prob  <- apply(P, 1, max)
top_idx   <- max.col(P, ties.method = "first")
top_group <- prob_cols[top_idx]

P2 <- P
P2[cbind(seq_len(nrow(P2)), top_idx)] <- -Inf
second_prob <- apply(P2, 1, max)

margin  <- max_prob - second_prob
entropy <- -rowSums(pmax(P, eps) * log(pmax(P, eps)))

probDF$max_prob    <- max_prob
probDF$top_group   <- top_group
probDF$second_prob <- second_prob
probDF$margin      <- margin
probDF$entropy     <- entropy


## Fetch Pan-UKB ancestry labels

dxdownload("Callum/Derived_datasets/Age_Sex_PRS_GA.csv") # Dataset of Age at Recruitment (p21022), Sex (p31), Genetically-inferred Sex (p22001), Standard PRS for Prostate cancer (p26267), Enhanced PRS for Prostate cancer (p26268), Genomic Ancestry (p30079), created using the UKB-RAP cohort browser

Genomic_ancestry <- read.csv("Age_Sex_PRS_GA.csv") %>%
  dplyr::select("eid", "p30079") %>%
  dplyr::rename("IID" = "eid", "Genomic_ancestry" = "p30079") %>%
  dplyr::filter(Genomic_ancestry == "European ancestry (EUR)" |
                  Genomic_ancestry == "African ancestry (AFR)" |
                  Genomic_ancestry == "East Asian ancestry (EAS)" |
                  Genomic_ancestry == "Central/South Asian ancestry (CSA)" |
                  Genomic_ancestry == "Middle Eastern ancestry (MID)" |
                  Genomic_ancestry == "Admixed American ancestry (AMR)"
  )

## Attach labels and PCs to compare

probDF_w_results <- merge(probDF, Genomic_ancestry, by = "IID", all.x = T)
probDF_w_results <- merge(probDF_w_results, UKB_PCs, by = "IID", all.x = T)

## Export to project

write.csv(probDF, "AncestryProbability2.csv", row.names = FALSE)
system("dx upload AncestryProbability2.csv")



#============================================================================================================#
# Sanity check 1 - does Pan-UKB-defined "Genomic ancestry" field match 50% likelihood of fitting in cluster? #
#============================================================================================================#

probDF_w_results %>%
  filter(
    Genomic_ancestry == "European ancestry (EUR)" & EUR < 0.5 |
      Genomic_ancestry == "African ancestry (AFR)" & AFR < 0.5 |
      Genomic_ancestry == "East Asian ancestry (EAS)" & EAS < 0.5 |
      Genomic_ancestry == "Admixed American ancestry (AMR)" & AMR < 0.5 |
      Genomic_ancestry == "Middle Eastern ancestry (MID)" & MID < 0.5 |
      Genomic_ancestry == "Central/South Asian ancestry (CSA)" & CSA < 0.5
  ) %>%
  summarise(n = n())


#====================================#
# Sanity check 2 - Visual separation #
#====================================#


library(ggplot2)

HGDP_1KG_PCs_w_labels <- HGDP_1KG_PCs_w_labels %>%
  dplyr::mutate(Genomic_ancestry = case_when(
    Region %in% c("EUR") ~ "European ancestry (EUR)",
    Region %in% c("AFR") ~ "African ancestry (AFR)",
    Region %in% c("EAS") ~ "East Asian ancestry (EAS)",
    Region %in% c("CSA") ~ "Central/South Asian ancestry (CSA)",
    Region %in% c("MID") ~ "Middle Eastern ancestry (MID)",
    Region %in% c("AMR") ~ "Admixed American ancestry (AMR)",
  ))

## First take a look at HGDP+1KG Principal Components

ggplot(HGDP_1KG_PCs_w_labels, aes(x = PC1, y = PC2, color = Genomic_ancestry)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  # coord_cartesian(xlim = c(-50, 80), ylim = c(-50, 60)) +
  labs(
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Genomic ancestry"
  )

## Second take a look at UKB Principal Components

test_population <- probDF_w_results #%>%
#  dplyr::filter(
#    ethnicity_group_narrow == "Mixed White and Black" | ethnicity_group_narrow == "Black" | ethnicity_group_narrow == "White"
#  )

ggplot(test_population, aes(x = PC1, y = PC2, color = Genomic_ancestry)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  # coord_cartesian(xlim = c(-50, 80), ylim = c(-50, 60)) +
  labs(
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Genomic_ancestry"
  )



#=========================================================#
# Sanity check 3  - assign own Pan-UKB groups and compare #
#=========================================================#

library(MASS); library(data.table)

k <- length(keepPCs)

# Reference stats (means & inverse covariances) per Region from 'train'

ref_stats <- as.data.table(train)[, {
  X  <- as.matrix(.SD)
  mu <- colMeans(X)
  S  <- cov(X, use = "pairwise.complete.obs")
  S  <- S + diag(1e-8, k)              # small ridge for stability
  iS <- chol2inv(chol(S))              # inverse via Cholesky
  list(mu = list(mu), iS = list(iS), S = list(S))
}, by = Region, .SDcols = keepPCs]
setkey(ref_stats, Region)

# Helper: fast vectorized MD^2 given X, mu, iS

md_block <- function(X, mu, iS) {
  C <- sweep(X, 2, mu, "-")
  M <- C %*% iS
  rowSums(M * C)                       # MD^2
}

# (Optional) empirical per-group MD^2 thresholds from REFERENCE

#    This adapts the cutoff to each group's real dispersion
emp_cut <- ref_stats[, {
  X  <- as.matrix(train[Region==.BY$Region, ..keepPCs])
  md <- md_block(X, mu[[1]], iS[[1]])
  list(thr_emp95 = as.numeric(quantile(md, 0.95, na.rm = TRUE)))
}, by = Region]

# Choose thresholding strategy

use_empirical <- TRUE
if (use_empirical) {
  ref_stats <- merge(ref_stats, emp_cut, by="Region", all.x=TRUE)
} else {
  ref_stats[, thr_emp95 := qchisq(0.95, df = k)]
}

# Attach predicted top_group and max_prob to the UKB PCs (from your RF)

ukb <- as.data.table(UKB_PCs)[, c("IID", keepPCs), with = FALSE]
probDT <- as.data.table(probDF)[, .(IID, top_group, max_prob)]
ukb <- merge(ukb, probDT, by = "IID", all.x = TRUE)

# probability > 0.5 (Pan-UKBB rule)

ukb[, stage1_label := ifelse(!is.na(top_group) & max_prob >= 0.5, top_group, "Other")]

# Compute MD^2 to the stage1_label centroid

ukb[stage1_label != "Other", MD2 := {
  st <- ref_stats[.BY$stage1_label]
  X  <- as.matrix(.SD)
  md_block(X, st$mu[[1]], st$iS[[1]])
}, by = stage1_label, .SDcols = keepPCs]

# Apply per-group cutoff

ukb <- merge(ukb, ref_stats[, .(Region, thr_emp95)], 
             by.x = "stage1_label", by.y = "Region", all.x = TRUE)

ukb[, stage2_label := stage1_label]
ukb[stage1_label != "Other" & !is.na(MD2) & MD2 > thr_emp95, stage2_label := "Other"]

### Compare to Pan-UKB labels

library(caret)

# Bring in Pan-UKBB labels (already in probDF_w_results data)

pan <- as.data.table(probDF_w_results)[, .(IID, Genomic_ancestry)]

# Map Pan-UKBB long labels -> short codes
map_lab <- c(
  "European ancestry (EUR)"             = "EUR",
  "African ancestry (AFR)"              = "AFR",
  "East Asian ancestry (EAS)"           = "EAS",
  "Admixed American ancestry (AMR)"     = "AMR",
  "Middle Eastern ancestry (MID)"       = "MID",
  "Central/South Asian ancestry (CSA)"  = "CSA",
  "Other/Unassigned/NA"                 = "Other"  # just in case
)
pan[, pan_label := unname(map_lab[Genomic_ancestry])]
pan[is.na(pan_label), pan_label := "Other"]

# Merge with your stage-2 labels

cmp <- merge(ukb[, .(IID, stage1_label, stage2_label, max_prob, MD2)], 
             pan[, .(IID, pan_label)], by = "IID", all.x = TRUE)

# Confusion matrices

cm1 <- confusionMatrix(factor(cmp$stage1_label), factor(cmp$pan_label))
cm2 <- confusionMatrix(factor(cmp$stage2_label), factor(cmp$pan_label))

cat("\n=== Stage 1 (prob >= 0.5) vs Pan-UKBB ===\n")
print(cm1$overall["Accuracy"]); print(cm1$byClass[, c("Precision","Recall","F1")])

cat("\n=== Stage 2 (prob >= 0.5 + centroid distance) vs Pan-UKBB ===\n")
print(cm2$overall["Accuracy"]); print(cm2$byClass[, c("Precision","Recall","F1")])

# Coverage / drop counts

coverage1 <- mean(cmp$stage1_label != "Other", na.rm=TRUE)
coverage2 <- mean(cmp$stage2_label != "Other", na.rm=TRUE)
cat(sprintf("\nCoverage Stage 1: %.2f%%   Coverage Stage 2: %.2f%%\n", 100*coverage1, 100*coverage2))

# Inspect the set that failed the distance filter

failed <- cmp[stage1_label != "Other" & stage2_label == "Other"]
summary(failed$MD2)

