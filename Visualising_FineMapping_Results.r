## Taking a look at fine mapping results

# Note: these fine mapping results were created using Robin Beaumont's "Finemapping_WGS” applet. Preparation of GWAS summary statistics for use with this tool can be found in the script "Conti_summary_stats_into_Fine_Mapping_ready_version"

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AllMen.cs")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AllMen.snp")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AllMen.summary")

dxdownload("Callum/FineMappingResults/UsingContiGWAS/Conti_FineMap_rs72725854_AllMen.cs")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/Conti_FineMap_rs72725854_AllMen.snp")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/Conti_FineMap_rs72725854_AllMen.summary")


dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_BlackOnly.cs")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_BlackOnly.snp")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_BlackOnly.summary")

dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AFR.cs")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AFR.snp")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AFR.summary")

dxdownload("Callum/FineMappingResults/UsingWangGWAS/Wang_FineMap_rs72725854_BlackOnly.cs")
dxdownload("Callum/FineMappingResults/UsingWangGWAS/Wang_FineMap_rs72725854_BlackOnly.snp")
dxdownload("Callum/FineMappingResults/UsingWangGWAS/Wang_FineMap_rs72725854_BlackOnly.summary")

dxdownload("Callum/FineMappingResults/UsingWangGWAS/Wang_FineMap_rs72725854_BlackOnly_MAF0.001.cs")
dxdownload("Callum/FineMappingResults/UsingWangGWAS/Wang_FineMap_rs72725854_BlackOnly_MAF0.001.snp")
dxdownload("Callum/FineMappingResults/UsingWangGWAS/Wang_FineMap_rs72725854_BlackOnly_MAF0.001.summary")


cs_generalGWAS <- read.delim("FineMap_rs72725854_AllMen.cs")
snp_generalGWAS <- read.delim("FineMap_rs72725854_AllMen.snp")
summary_generalGWAS <- read.delim("~/FineMap_rs72725854_AllMen.summary", comment.char="#")

cs_generalGWAS_repeated <- read.delim("Conti_FineMap_rs72725854_AllMen.cs")
snp_generalGWAS_repeated <- read.delim("Conti_FineMap_rs72725854_AllMen.snp")
summary_generalGWAS_repeated <- read.delim("~/Conti_FineMap_rs72725854_AllMen.summary", comment.char="#")

cs_BlackOnlyGWAS <- read.delim("FineMap_rs72725854_BlackOnly.cs")
snp_BlackOnlyGWAS <- read.delim("FineMap_rs72725854_BlackOnly.snp")
summary_BlackOnlyGWAS <- read.delim("~/FineMap_rs72725854_BlackOnly.summary", comment.char="#")

cs_AFRGWAS <- read.delim("FineMap_rs72725854_AFR.cs")
snp_AFRGWAS <- read.delim("FineMap_rs72725854_AFR.snp")
summary_AFRGWAS <- read.delim("~/FineMap_rs72725854_AFR.summary", comment.char="#")

#cs_Wang_BlackOnlyGWAS <- read.delim("Wang_FineMap_rs72725854_BlackOnly.cs")
#snp_Wang_BlackOnlyGWAS <- read.delim("Wang_FineMap_rs72725854_BlackOnly.snp")
#summary_Wang_BlackOnlyGWAS <- read.delim("~/Wang_FineMap_rs72725854_BlackOnly.summary", comment.char="#")

cs_Wang_BlackOnlyGWAS <- read.delim("Wang_FineMap_rs72725854_BlackOnly_MAF0.001.cs")
snp_Wang_BlackOnlyGWAS <- read.delim("Wang_FineMap_rs72725854_BlackOnly_MAF0.001.snp")
summary_Wang_BlackOnlyGWAS <- read.delim("~/Wang_FineMap_rs72725854_BlackOnly_MAF0.001.summary", comment.char="#")


## Add MAX PIP to snp_Wang_BlackOnlyGWAS

library(dplyr)

snp_Wang_BlackOnlyGWAS_edited <- snp_Wang_BlackOnlyGWAS %>%
  dplyr::mutate(
    MAX_PIP = pmax(PIP.CS1., PIP.CS2., PIP.CS3., PIP.CS4., PIP.CS5.)
  )

## Plot variants on x axis, MAX_PIP on y axis

library(ggplot2)
library(readr)

snp_Wang_BlackOnlyGWAS_edited <- snp_Wang_BlackOnlyGWAS_edited %>% mutate(BP = reorder(BP, MAX_PIP))


# Ensure numeric ordering of IDs (works for "1", "ID_2", "sample-10", etc.)
snp_Wang_BlackOnlyGWAS_edited <- snp_Wang_BlackOnlyGWAS_edited %>%
  mutate(
    BP = as.character(BP),
    BP_num = parse_number(BP)          # pulls out the numeric part safely
  ) %>%
  arrange(BP_num) %>%
  mutate(BP = factor(BP, levels = unique(BP)))  # lock in numeric order on the x-axis

wang_bp_levels <- as.character(snp_Wang_BlackOnlyGWAS_edited$BP)
wang_x_axis_breaks <- unique(wang_bp_levels[
  round(seq(1, length(wang_bp_levels), length.out = 5))
])

snp_Wang_BlackOnlyGWAS_peaks <- snp_Wang_BlackOnlyGWAS_edited %>%
  mutate(
    is_peak = MAX_PIP > dplyr::lag(MAX_PIP, default = -Inf) &
      MAX_PIP > dplyr::lead(MAX_PIP, default = -Inf)
  ) %>%
  filter(is_peak, MAX_PIP > 0.1)

ggplot(snp_Wang_BlackOnlyGWAS_edited, aes(x = BP, y = MAX_PIP)) +
  geom_col(width = 2, fill = "#04dca4") +
  geom_text(
    data = snp_Wang_BlackOnlyGWAS_peaks,
    aes(y = MAX_PIP * 0.5, label = BP),
    angle = 90,
    size = 2.2,
    vjust = 0.5,
    check_overlap = TRUE
  ) +
  scale_x_discrete(breaks = wang_x_axis_breaks) +
  labs(x = "SNP position on chromosome 8", y = "MAX PIP", title = "") +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    # Vertical labels 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)
  )


## Add MAX PIP to snp_BlackOnlyGWAS (Conti)

library(dplyr)

snp_BlackOnlyGWAS_edited <- snp_BlackOnlyGWAS %>%
  dplyr::mutate(
    MAX_PIP = pmax(PIP.CS1., PIP.CS2., PIP.CS3.)
  )

## Plot variants on x axis, MAX_PIP on y axis

library(ggplot2)
library(readr)

snp_BlackOnlyGWAS_edited <- snp_BlackOnlyGWAS_edited %>% mutate(BP = reorder(BP, MAX_PIP))


# Ensure numeric ordering of IDs (works for "1", "ID_2", "sample-10", etc.)
snp_BlackOnlyGWAS_edited <- snp_BlackOnlyGWAS_edited %>%
  mutate(
    BP = as.character(BP),
    BP_num = parse_number(BP)          # pulls out the numeric part safely
  ) %>%
  arrange(BP_num) %>%
  mutate(BP = factor(BP, levels = unique(BP)))  # lock in numeric order on the x-axis

conti_bp_levels <- as.character(snp_BlackOnlyGWAS_edited$BP)
conti_x_axis_breaks <- unique(conti_bp_levels[
  round(seq(1, length(conti_bp_levels), length.out = 5))
])

snp_BlackOnlyGWAS_peaks <- snp_BlackOnlyGWAS_edited %>%
  mutate(
    is_peak = MAX_PIP > dplyr::lag(MAX_PIP, default = -Inf) &
      MAX_PIP > dplyr::lead(MAX_PIP, default = -Inf)
  ) %>%
  filter(is_peak, MAX_PIP > 0.1)

ggplot(snp_BlackOnlyGWAS_edited, aes(x = BP, y = MAX_PIP)) +
  geom_col(width = 2, fill = "#04dca4") +
  geom_text(
    data = snp_BlackOnlyGWAS_peaks,
    aes(y = MAX_PIP * 0.5, label = BP),
    angle = 90,
    size = 2.2,
    vjust = 0.5,
    check_overlap = TRUE
  ) +
  scale_x_discrete(breaks = conti_x_axis_breaks) +
  labs(x = "SNP position on chromosome 8", y = "MAX PIP", title = "") +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    # Vertical labels 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)
  )


## View initial summary stats document, to find real betas

system(paste("dx download", "Callum/WangGWAS/Wang2023African_harmonised.tsv")) ## This is Wang's AFR-specific GWAS summary stats, haramonised to GRCh38, as downloadable here: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90274001-GCST90275000/GCST90274715/. Go to /harmonised and download GCST90274715.h.tsv.gz

WangGWASsummaryStatsAFR <- read.delim("Wang2023African_harmonised.tsv", sep = "\t") %>%
  dplyr::rename(
    "CHROM" = "chromosome",
    "GENPOS" = "base_pair_location",
    "rsid" = "rsid",
    "ALLELE1" = "effect_allele",
    "ALLELE0" = "other_allele",
    "BETA" = "beta",
    "SE" = "standard_error",
    "P" = "p_value",
  ) %>%
  dplyr::filter(CHROM == 8)

betafinder <- WangGWASsummaryStatsAFR %>%
  dplyr::filter(GENPOS == 127062570 | GENPOS == 127091724 | GENPOS == 127012808 | GENPOS == 127102860 | GENPOS == 127146838)


