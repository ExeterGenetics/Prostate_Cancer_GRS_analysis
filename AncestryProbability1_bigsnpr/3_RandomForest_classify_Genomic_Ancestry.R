############# Random Forest to classify Genomic Ancestry (HGDP+1KG edition) #############

# Install Packages and download HGDP+1KG PCs, plus projected Biobank PCs

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

for (p in c("data.table","randomForest","MASS","ggplot2")) {
  if (!require(p, character.only = TRUE)) install.packages(p)
}
library(data.table); library(randomForest); library(dplyr)

dxdownload("Callum/HGDP_1KG/HGDP_1KG_PCs.csv") # From Step 2
dxdownload("Callum/HGDP_1KG/release_3.1_secondary_analyses_hgdp_1kg_v2_metadata_and_qc_gnomad_meta_updated.tsv")  # This can be downloaded from https://console.cloud.google.com/storage/browser/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc
dxdownload("Callum/HGDP_1KG/HGDP_1KG_PCs_UKB.csv") # From Step 2
dxdownload("Callum/Derived_datasets/Age_Sex_PRS_GA.csv") # Dataset of Age at Recruitment (p21022), Sex (p31), Genetically-inferred Sex (p22001), Standard PRS for Prostate cancer (p26267), Enhanced PRS for Prostate cancer (p26268), Genomic Ancestry (p30079), created using the UKB-RAP cohort browser
dxdownload("Callum/Derived_datasets/ethnicity.csv") # Dataset of self-reported Ethnic Background (p21000_i0), created using the UKB-RAP cohort browser

## Working out which column of the HGDP+1KG metadata file to use as ancestry labels

HGDP_1KG_metadata <- read.delim("release_3.1_secondary_analyses_hgdp_1kg_v2_metadata_and_qc_gnomad_meta_updated.tsv")
HGDP_1KG_metadata_ancestrygroups <- HGDP_1KG_metadata %>%
  dplyr::select(c("s", "project_meta.sample_id", "project_meta.project_pop", "project_meta.project_subpop", "project_meta.subpop_description", "population_inference.training_pop", "population_inference.pop", "population_inference.training_pop_all", "bergstrom.region", "hgdp_tgp_meta.Study.region", "hgdp_tgp_meta.Population", "hgdp_tgp_meta.Genetic.region", "population"))

table(HGDP_1KG_metadata$project_meta.project_pop)
table(HGDP_1KG_metadata$hgdp_tgp_meta.Genetic.region)

## Decided on using "Genetic Region" as it is closest to what we see in UK Biobank "Genomic Ancestry" variable

HGDP_1KG_metadata_ancestrylabel <- HGDP_1KG_metadata_ancestrygroups %>%
  dplyr::select(c("s", "hgdp_tgp_meta.Genetic.region")) %>%
  dplyr::rename("IID" = "s", "Region" = "hgdp_tgp_meta.Genetic.region")

HGDP_1KG_PCs <- read.csv("HGDP_1KG_PCs.csv")
HGDP_1KG_PCs_w_labels <- merge(HGDP_1KG_PCs, HGDP_1KG_metadata_ancestrylabel, by = "IID") %>%
  dplyr::select(c("IID", "PC1", "PC2", "PC3", "PC4", "PC5" , "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "Region"))

# Setup training data

train <- HGDP_1KG_PCs_w_labels %>%
  dplyr::filter(!is.na(Region) & !is.na(PC1) & Region != "OCE")





## Train Random Forest on HGDP+1000G PCs to predict ancestry

keepPCs <- paste0("PC", 1:20)  

set.seed(1)
rf <- randomForest(x = train[keepPCs],
                   y = factor(train$Region),
                   ntree = 10000, importance = TRUE)




# Load in HGDP+1KG projected UKB PCs

rf_data <- read.csv("HGDP_1KG_PCs_UKB.csv") %>%
  dplyr::select(c("IID", 
                  "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20"))

Genomic_ancestry <- read_csv("Age_Sex_PRS_GA.csv") %>%
  dplyr::select(c("eid", "p30079")) %>%
  dplyr::rename(c("IID" = "eid", "Genomic_ancestry" = "p30079"))

ethnicity <- read_csv("ethnicity.csv")
ethnicity <- ethnicity %>%
  dplyr::rename(c("IID" = "eid", 'Ethnicity' = 'p21000_i0'))

ethnicity <- ethnicity %>% # stricter ethnicity grouping
  dplyr::mutate(ethnicity_group_narrow = case_when(
    Ethnicity %in% c("British", "White", "Any other white background", "Irish") ~ "White",
    Ethnicity %in% c("African", "Caribbean", "Any other Black background", "Black or Black British") ~ "Black",
    Ethnicity %in% c("Pakistani", "Indian", "Bangladeshi") ~ "South Asian",
    Ethnicity %in% c("Any other Asian background", "White and Asian", "Chinese", "Asian or Asian British") ~ "East Asian",
    Ethnicity %in% c("White and Black African", "White and Black Caribbean") ~ "Mixed White and Black",
    Ethnicity %in% c("Other ethnic group", "Any other mixed background", "Prefer not to answer", "Do not know", "Mixed", "NA") ~ "Other",
    TRUE ~ NA_character_
  ))

rf_data <- merge(rf_data, Genomic_ancestry, by = "IID", all.x=T)
rf_data <- merge(rf_data, ethnicity, by = "IID", all.x=T)

# Remove anyone who has NA for PCAs

rf_data <- rf_data %>%
  dplyr::filter(
    !is.na(PC1),
  )


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

rf_data <- exclude_withdrawn(rf_data)



# Predict for all participants with Principal_Components
probs <- predict(rf, rf_data[, keepPCs], type = "prob")

probDF <- cbind(
  IID = rf_data$IID,
  as.data.frame(probs, check.names = FALSE)
)



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



probDF_w_results <- merge(probDF, rf_data, by = "IID") %>%
  dplyr::select(c("IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20",  "AMR", "AFR", "CSA", "EAS", "EUR", "MID", "ethnicity_group_narrow", "Genomic_ancestry", "max_prob", "top_group", "second_prob"))


## Save

write.csv(probDF, "AncestryProbability1.csv")
system(paste("dx upload", "AncestryProbability1.csv"))


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
    Region %in% c("CSA") ~ "Central/South Asian ancestry (AFR)",
    Region %in% c("MID") ~ "Middle Eastern ancestry (MID)",
    Region %in% c("AMR") ~ "Admixed American ancestry (AFR)",
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

test_population <- probDF_w_results %>%
  dplyr::filter(
    ethnicity_group_narrow == "Mixed White and Black" | ethnicity_group_narrow == "Black" | ethnicity_group_narrow == "White"
  )

ggplot(test_population, aes(x = PC1, y = PC2, color = ethnicity_group_narrow)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  # coord_cartesian(xlim = c(-50, 80), ylim = c(-50, 60)) +
  labs(
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Self-reported ethnicity"
  )


#========================================================#
# Sanity check 3 - assign own Pan-UKB groups and compare #
#========================================================#

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

# Empirical per-group MD^2 thresholds from REFERENCE

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

# Attach predicted top_group and max_prob to the UKB PCs (from RF)

ukb <- as.data.table(ukb_df)[, c("IID", keepPCs), with = FALSE]
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

# Compare to Pan-UKB labels

library(caret)

# Bring in Pan-UKBB labels (already in your probDF_w_results data)
#    Keep only IIDs present in UKB PCs and harmonize
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

# Merge with stage-2 labels
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

