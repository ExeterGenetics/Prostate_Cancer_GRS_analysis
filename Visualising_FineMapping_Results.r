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

cs_Wang_BlackOnlyGWAS <- read.delim("Wang_FineMap_rs72725854_BlackOnly.cs")
snp_Wang_BlackOnlyGWAS <- read.delim("Wang_FineMap_rs72725854_BlackOnly.snp")
summary_Wang_BlackOnlyGWAS <- read.delim("~/Wang_FineMap_rs72725854_BlackOnly.summary", comment.char="#")




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

snp_Wang_BlackOnlyGWAS_peaks <- snp_Wang_BlackOnlyGWAS_edited %>%
  mutate(
    is_peak = MAX_PIP > dplyr::lag(MAX_PIP, default = -Inf) &
      MAX_PIP > dplyr::lead(MAX_PIP, default = -Inf)
  ) %>%
  filter(is_peak, MAX_PIP > 0.1)

ggplot(snp_Wang_BlackOnlyGWAS_edited, aes(x = BP, y = MAX_PIP)) +
  geom_col(width = 0.9, fill = "#04dca4") +
  geom_text(
    data = snp_Wang_BlackOnlyGWAS_peaks,
    aes(y = MAX_PIP * 0.5, label = BP),
    angle = 90,
    size = 2.2,
    vjust = 0.5,
    check_overlap = TRUE
  ) +
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

snp_BlackOnlyGWAS_peaks <- snp_BlackOnlyGWAS_edited %>%
  mutate(
    is_peak = MAX_PIP > dplyr::lag(MAX_PIP, default = -Inf) &
      MAX_PIP > dplyr::lead(MAX_PIP, default = -Inf)
  ) %>%
  filter(is_peak, MAX_PIP > 0.1)

ggplot(snp_BlackOnlyGWAS_edited, aes(x = BP, y = MAX_PIP)) +
  geom_col(width = 0.9, fill = "#04dca4") +
  geom_text(
    data = snp_BlackOnlyGWAS_peaks,
    aes(y = MAX_PIP * 0.5, label = BP),
    angle = 90,
    size = 2.2,
    vjust = 0.5,
    check_overlap = TRUE
  ) +
  labs(x = "SNP position on chromosome 8", y = "MAX PIP", title = "") +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    # Vertical labels 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)
  )
