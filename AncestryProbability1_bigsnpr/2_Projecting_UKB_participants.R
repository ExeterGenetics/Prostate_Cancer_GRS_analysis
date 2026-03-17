#====================================================================#
# Projecting UKB participants to HGDP+1KG Principal Components Space #
#====================================================================#

#############################
# ----- Step 0: Setup ----- #            Run at start of every session, Takes ~ 5 min
#############################

## Avoid nested parallelism

Sys.setenv(
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  OPENBLAS_NUM_THREADS = "1",
  BLAS_NUM_THREADS     = "1",
  NUMEXPR_NUM_THREADS  = "1"
)
if (!requireNamespace("RhpcBLASctl", quietly = TRUE)) install.packages("RhpcBLASctl")
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

## Load in ukbextractR ## 

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

## Load in HGDP + 1KG PLINK files created in Step 1

dxdownload("Callum/HGDP_1KG/HGDP_1KGG_autosomes_pruned_masked_norelatives_p1.bed")
dxdownload("Callum/HGDP_1KG/HGDP_1KGG_autosomes_pruned_masked_norelatives_p1.bim")
dxdownload("Callum/HGDP_1KG/HGDP_1KGG_autosomes_pruned_masked_norelatives_p1.fam")

## Load in bigsnpr (package for projecting PCs)

install.packages("bigsnpr")         
library(bigsnpr)

## Install PLINK

library(bigsnpr)
plink1_path <- bigsnpr::download_plink(dir = "~")
Sys.chmod(plink1_path, "0755")
Sys.setenv(PATH = paste(dirname(plink1_path), Sys.getenv("PATH"), sep = ":"))
plink1 <- Sys.which("plink"); stopifnot(nzchar(plink1))

## Note - each Step 1-4 is self-contained and writes a file to the project directory that can be used in the next step

###############################################################################
# ----- Step 1: Calculate HGDP+1KG PCAs and retrieve the variant subset ----- #      Takes ~ 40 minutes
###############################################################################

library(data.table)

ref_base <- "/home/rstudio-server/HGDP_1KGG_autosomes_pruned_masked_norelatives_p1"
obj.ref  <- bed(paste0(ref_base, ".bed"))

k <- 20
svd.ref <- bed_autoSVD(
  obj.bed   = obj.ref,
  k         = k,
  thr.r2    = 0.2,
  roll.size = 50,
  ncores    = 1   # keep single-core here
)

ref_subset <- attr(svd.ref, "subset")
ref_bim    <- fread(paste0(ref_base, ".bim"),
                    col.names = c("CHR","SNP","CM","POS","A1","A2"))
ref_bim_sub <- ref_bim[ref_subset, .(CHR, POS, A1, A2)]

## Save ref_bim_sub (the SNPs actually used to calculate HGDP+1KG PCs)

write.csv(ref_bim_sub, "SNPs_in_HGDP_1KG_PCs.csv")
system(paste("dx upload", "SNPs_in_HGDP_1KG_PCs.csv"))

###############################################################################
# ---- Step 2: Lift over GRCh38 (HGDP + 1KG type) into GRCh37 (UKB type) ---- #    Takes ~10 minutes
###############################################################################

dxdownload("Callum/HGDP_1KG/SNPs_in_HGDP_1KG_PCs.csv")
ref_bim_sub <- read_csv("SNPs_in_HGDP_1KG_PCs.csv")

BiocManager::install(c("GenomicRanges", "rtracklayer"), ask = FALSE, update = FALSE) # Takes ~7 minutes
library(GenomicRanges)
library(rtracklayer)
library(data.table)

dxdownload("Callum/tools/hg38ToHg19.over.chain")	# This chain can be installed from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
chain <- import.chain("hg38ToHg19.over.chain")

gr38 <- GRanges(seqnames = paste0("chr", ref_bim_sub$CHR),
                ranges   = IRanges(start = ref_bim_sub$POS, end = ref_bim_sub$POS))

hits   <- liftOver(gr38, chain)
one2one <- lengths(hits) == 1L
gr19    <- unlist(hits[one2one])
src_idx <- which(one2one)  # indices in ref_bim_sub that mapped uniquely

ref19 <- data.table(
  CHR = as.integer(sub("^chr","", as.character(seqnames(gr19)))),
  POS = as.integer(start(gr19)),
  SRC = src_idx
)
setkey(ref19, CHR, POS)
ref19 <- unique(ref19, by = c("CHR","POS"))   # to guard rare duplicates


# Convert to data.table once

if (!data.table::is.data.table(ref_bim_sub)) data.table::setDT(ref_bim_sub)

# Drop ambiguous A/T and C/G early

is_ambig <- function(a1, a2) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  (a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C")
}
ref_bim_sub[, AMB := is_ambig(A1, A2)]

# Filter ref19 to non-ambiguous SNPs via the source indices from liftOver

ref19 <- ref19[SRC %in% which(!ref_bim_sub$AMB)]

write.csv(ref19, "ref19.csv")
system(paste("dx upload", "ref19.csv"))

#######################################################################################
# ---- Step 3: Build per‑chr keep‑lists from UKB .bim, then extract reduced sets ---- #    Takes ~11 minutes
#######################################################################################

dxdownload("Callum/HGDP_1KG/ref19.csv")
ref19 <- read_csv("ref19.csv")

## Load in UKB PLINK files (takes ~ 10 minutes) (requires a project with Bulk data dispensed)

for (chr in 1:22) {
  base <- sprintf("ukb22418_c%d_b0_v2", chr)
  dxdownload(paste0('"/Bulk/Genotype Results/Genotype calls/', base, '.bed"'))
  dxdownload(paste0('"/Bulk/Genotype Results/Genotype calls/', base, '.bim"'))
  dxdownload(paste0('"/Bulk/Genotype Results/Genotype calls/', base, '.fam"'))
}


for (chr in 1:22) {
  # Read UKB bim for this chr
  ukb_bim <- fread(sprintf("ukb22418_c%d_b0_v2.bim", chr),
                   col.names = c("CHR","SNP","CM","POS","A1","A2"))
  setkey(ukb_bim, CHR, POS)
  # Join by (CHR, POS) in hg19
  keep <- ukb_bim[ref19, nomatch = 0L, on = .(CHR, POS)]$SNP
  # Write keep-list and extract a small BED
  keepfile <- sprintf("keep_chr%d.snplist", chr)
  fwrite(data.table(SNP = keep), keepfile, col.names = FALSE)
  system2(plink1, c(
    "--bfile", sprintf("ukb22418_c%d_b0_v2", chr),
    "--extract", keepfile,
    "--make-bed",
    "--out", sprintf("ukb_refSNPs_c%d", chr)
  ))
}

# ---- Step 3.5: Merge the reduced per‑chr sets into one small UKB BED ---- #      Takes ~ 2 minutes


bases_small <- sprintf("ukb_refSNPs_c%d", 1:22)
writeLines(bases_small[-1], "_ukb_small_merge_list.txt")

system2(plink1, c(
  "--bfile", bases_small[1],
  "--merge-list", "_ukb_small_merge_list.txt",
  "--make-bed",
  "--out", "UKB_refSNPs_merged"
))

system(paste("dx upload", "UKB_refSNPs_merged.bed"))
system(paste("dx upload", "UKB_refSNPs_merged.bim"))
system(paste("dx upload", "UKB_refSNPs_merged.fam"))

############################################################
# ---- Step 4: Project with bigsnpr::bed_projectPCA() ---- #      Takes ~ 8 minutes
############################################################

dxdownload("Callum/HGDP_1KG/UKB_refSNPs_merged.bed")
dxdownload("Callum/HGDP_1KG/UKB_refSNPs_merged.bim")
dxdownload("Callum/HGDP_1KG/UKB_refSNPs_merged.fam")

dxdownload("Callum/tools/liftOver")
liftOver_bin <- "liftOver"
Sys.chmod(liftOver_bin, "0755")

ref_base <- "HGDP_1KGG_autosomes_pruned_masked_norelatives_p1"
obj.ref <- bed(paste0(ref_base, ".bed"))	# GRCh38
obj.ukb <- bed("UKB_refSNPs_merged.bed")  	# GRCh37

k <- 20
nc <- min(16, bigstatsr::nb_cores())  
res <- bed_projectPCA(
  obj.bed.ref = obj.ref,
  obj.bed.new = obj.ukb,
  k           = k,
  join_by_pos = TRUE,
  strand_flip = TRUE,
  build.ref   = "hg38",
  build.new   = "hg19",
  liftOver    = liftOver_bin,     
  verbose     = TRUE,
  ncores      = nc
)

ref_pcs <- predict(res$obj.svd.ref)
ukb_pcs <- res$OADP_proj


# ---- Step 4.5: turn into dataframes with ID column ---- #      Takes ~ 2 minutes


library(data.table)

## Read IDs in the exact order bigsnpr used

ref_ids <- fread(paste0(ref_base, ".fam"))$V2
ukb_ids <- fread("UKB_refSNPs_merged.fam")$V2

## Convert matrices and name columns

ref_mat <- as.matrix(ref_pcs)       # predict(res$obj.svd.ref)
ukb_mat <- as.matrix(ukb_pcs)       # res$OADP_proj
colnames(ref_mat) <- paste0("PC", seq_len(ncol(ref_mat)))
colnames(ukb_mat) <- paste0("PC", seq_len(ncol(ukb_mat)))

## Attach IIDs as rownames (locks the mapping into the object)

stopifnot(nrow(ref_mat) == length(ref_ids),
          nrow(ukb_mat) == length(ukb_ids))
rownames(ref_mat) <- ref_ids
rownames(ukb_mat) <- ukb_ids

## Build data.frames with an IID column (and optional cohort)

ref_df <- data.frame(IID = ref_ids, cohort = "HGDP+1KG", ref_mat, check.names = FALSE)
ukb_df <- data.frame(IID = ukb_ids, cohort = "UKB",      ukb_mat, check.names = FALSE)

## Sanity checks: order and uniqueness
stopifnot(all(ref_df$IID == ref_ids),
          all(ukb_df$IID == ukb_ids),
          length(unique(ref_ids)) == length(ref_ids),
          length(unique(ukb_ids)) == length(ukb_ids))




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

ukb_df <- exclude_withdrawn(ukb_df)

ref_df <- exclude_withdrawn(ukb_df)

## Write to project directory ##

write.csv(ukb_df, "HGDP_1KG_PCs_UKB.csv")
system(paste("dx upload", "HGDP_1KG_PCs_UKB.csv"))

write.csv(ref_df, "HGDP_1KG_PCs.csv")
system(paste("dx upload", "HGDP_1KG_PCs.csv"))





#########################################################################################
#---- Sanity check - how many of the HGDP+1KG SNPs are actually used in projection? ----#
#########################################################################################

library(data.table)

# Load the variants used in the HGDP+1KG PCA (Step 1 output) ---
#ref_bim_sub <- fread("SNPs_in_HGDP_1KG_PCs.csv")  # columns: CHR, POS, A1, A2
ref_bim_sub[, IDX := .I]                          # index rows to join via SRC

# Load GRCh37-lifted, non-ambiguous positions (Step 2 output) ---
#ref19 <- fread("ref19.csv")                       # columns: CHR, POS, SRC (1-based index into ref_bim_sub)

# Sanity: how many PCA SNPs (hg38) and how many lifted, non-ambig (hg19)?
N_ref_PCA              <- nrow(ref_bim_sub)
N_ref_PCA_hg19_nonamb <- nrow(ref19)

# Add ref alleles (from hg38 list) back to lifted positions via SRC index
ref19a <- merge(ref19, ref_bim_sub[, .(IDX, A1_ref = toupper(A1), A2_ref = toupper(A2))],
                by.x = "SRC", by.y = "IDX", all.x = TRUE)
ref19a[, c("CHR","POS") := .(as.integer(CHR), as.integer(POS))]
setkey(ref19a, CHR, POS)

# Load the merged reduced UKB BIM (Step 3.5) ---
ukb_bim <- fread("UKB_refSNPs_merged.bim",
                 col.names = c("CHR","SNP","CM","POS","A1","A2"))
ukb_bim[, `:=`(CHR = as.integer(CHR), POS = as.integer(POS),
               A1 = toupper(A1), A2 = toupper(A2))]
setkey(ukb_bim, CHR, POS)

# Join by position (what bed_projectPCA(join_by_pos=TRUE) does) ---
M <- merge(ukb_bim, ref19a[, .(CHR, POS, A1_ref, A2_ref)],
           by = c("CHR","POS"), all = FALSE)

# Count raw position overlap
N_overlap_pos <- nrow(M)

# Allele-orientation logic: direct match, or safe reverse-strand complement ---
comp_base <- c(A="T", T="A", C="G", G="C")
comp_vec  <- function(v) unname(comp_base[toupper(v)])

is_ambig <- function(a,b) {
  (a %in% c("A","T") & b %in% c("A","T")) | (a %in% c("C","G") & b %in% c("C","G"))
}

M[, `:=`(
  ambig        = is_ambig(A1, A2),
  match_direct = ( (A1==A1_ref & A2==A2_ref) | (A1==A2_ref & A2==A1_ref) ),
  match_comp   = {
    c1 <- comp_vec(A1); c2 <- comp_vec(A2)
    ( (c1==A1_ref & c2==A2_ref) | (c1==A2_ref & c2==A1_ref) )
  }
)]

# Keep if direct/swap OR (complement and non-ambiguous)
M[, keep := match_direct | (!ambig & match_comp)]
N_used_final <- M[, sum(keep)]

# Summaries ---
cat("\n=== Projection overlap summary (OADP pipeline) ===\n")
cat("Variants used to compute HGDP+1KG PCA (hg38):             ", N_ref_PCA, "\n")
cat("Lifted to hg19 and non-ambiguous (A/T,C/G removed):       ", N_ref_PCA_hg19_nonamb, "\n")
cat("UKB reduced set: position-overlap with lifted set:        ", N_overlap_pos, "\n")
cat("UKB projection: allele-consistent (kept) among overlaps:  ", N_used_final, "\n")
cat(sprintf("Proportion of PCA SNPs (hg38) that contributed to projection: %.2f%%\n",
            100 * N_used_final / N_ref_PCA))

# Optional: per-chromosome counts
by_chr <- M[, .(
  n_pos_overlap = .N,
  n_kept        = sum(keep)
), by = CHR][order(CHR)]
print(by_chr)




