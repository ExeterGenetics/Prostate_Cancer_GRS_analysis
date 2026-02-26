## Scatter plot HGDP+1KG Principal Components

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

library(ggplot2)

dxdownload("Callum/HGDP_1KG/HGDP_1KG_PCs.csv") # HGDP+1KG Principal components as calculated using AncestryProbability1_bigsnpr folder scripts
dxdownload("Callum/HGDP_1KG/release_3.1_secondary_analyses_hgdp_1kg_v2_metadata_and_qc_gnomad_meta_updated.tsv") # This can be downloaded from https://console.cloud.google.com/storage/browser/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc
dxdownload("Callum/HGDP_1KG/hgdp_tgp_pca_covid19hgi_snps_scores.txt") # can be downloaded from:
# https://console.cloud.google.com/storage/browser/covid19-hg-public/pca_projection;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

## Working out which column of the HGDP+1KG metadata file to use as ancestry labels

HGDP_1KG_metadata <- read.delim("release_3.1_secondary_analyses_hgdp_1kg_v2_metadata_and_qc_gnomad_meta_updated.tsv")

## Decided on using "Genetic Region" as it is closest to what we see in UK Biobank "Genomic Ancestry" variable

HGDP_1KG_metadata_ancestrylabel <- HGDP_1KG_metadata %>%
  dplyr::select(c("s", "hgdp_tgp_meta.Genetic.region")) %>%
  dplyr::rename("IID" = "s", "Genomic_ancestry" = "hgdp_tgp_meta.Genetic.region")

HGDP_1KG_PCs <- read.csv("HGDP_1KG_PCs.csv") # PCs as I calculated them
HGDP_1KG_PCs_w_labels <- merge(HGDP_1KG_PCs, HGDP_1KG_metadata_ancestrylabel, by = "IID")

HGDP_1KG_PCs_og <- read.delim("hgdp_tgp_pca_covid19hgi_snps_scores.txt") %>% # Publicly available PCs
  dplyr::rename("IID" = "s")
HGDP_1KG_PCs_og_w_labels <- merge(HGDP_1KG_PCs_og, HGDP_1KG_metadata_ancestrylabel, by = "IID")

################################################################################
#-------------------------------- Setup done ----------------------------------#
################################################################################

#test_population <- HGDP_1KG_PCs_og_w_labels              # Publicly available PCs
#test_population <- HGDP_1KG_PCs_w_labels                 # PCs as I calculated them
test_population <- PCa_iv_covariates_GRS_clean           # Projected PCs of UKB participants - requires running "Logistic_Regresions_testing_Conti_GRS.R" to create

ggplot(test_population, aes(x = PC1, y = PC2, color = Genomic_ancestry)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  # coord_cartesian(xlim = c(-50, 80), ylim = c(-50, 60)) +
  labs(
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Genomic ancestry"
  )

