#############################################################################
# Symptomatic Angle - how well can GRS predict Prostate cancer diagnosis in #
# ------------------- symptomatic participants? ----------------------------#
#############################################################################

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

install.packages("RMySQL")
library(RMySQL)
library(dplyr)
install.packages("readr")
library(readr)
install.packages("readstata13")
library(readstata13)
library(ggplot2)
install.packages("pROC")
library(pROC)
library(matrixStats)
install.packages("survminer")
library(survminer)
library(survival)
install.packages("tidyverse")
library(tidyverse) 
install.packages("DiagrammeR")
library(DiagrammeR)
library(epiR)

## Download files from project workspace into RStudio session

dxdownload("Callum/Derived_datasets/KLK3.csv")			# Dataset of KLK3 (Olink Instance 0 NPX Result proteomics) readings, created using the UKB-RAP cohort browser
dxdownload("Callum/Derived_datasets/imputed_rs72725854.csv")	# Dataset of rs72725854 SNP carrier status, created using "Extracting_rs72725854_carriers.R"
dxdownload("Callum/Derived_datasets/Age_Sex_PRS_GA.csv")	# Dataset of Age at Recruitment (p21022), Sex (p31), Genetically-inferred Sex (p22001), Standard PRS for Prostate cancer (p26267), Enhanced PRS for Prostate cancer (p26268), Genomic Ancestry (p30079), created using the UKB-RAP cohort browser
dxdownload("Callum/Derived_datasets/ethnicity.csv")		# Dataset of self-reported Ethnic Background (p21000_i0), created using the UKB-RAP cohort browser
dxdownload("Callum/HGDP_1KG/HGDP_1KG_PCs_UKB2_scaled.csv")	# Projected Principal Components of UKB participants in the principal components space of the HGDP+1000 Genomes reference panel, calculated using "AncestryProbability2calculation.R" in AncestryProbability2_plink folder
dxdownload("Callum/Derived_datasets/AncestryProbability.csv")	# Genetic similarity probabilities for Ancestry group calculated using AncestryProbability1_bigsnpr folder scripts
dxdownload("Callum/Derived_datasets/AncestryProbability2.csv") 	# Genetic similarity probabilities for Ancestry group calculated using AncestryProbability2_plink folder script

dxdownload("Callum/Derived_datasets/FH_PrCa_BrCa.csv")		# Dataset of Family History of Prostate Cancer and Breast Cancer, created using "Family_history.R" script
dxdownload("Callum/GRSs/Conti_multi_ethnic.pgs.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "Multiethnic Analysis"
dxdownload("Callum/GRSs/Conti_European.pgs.tsv") 		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "European...16"
dxdownload("Callum/GRSs/Conti_African.pgs.tsv")			# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "African...19"
dxdownload("Callum/GRSs/Conti_East_Asian.pgs.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "East Asian...22"
dxdownload("Callum/GRSs/Conti_Hispanic.pgs.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "Hispanic...25"

dxdownload("Callum/Derived_datasets/OR_adjustedGRS.tsv")	# Calculated using "Create_adjustedGRS_weighting_Conti_odds_ratios_by_AncestryProbability2.R" script


system("dx download Callum/LUTS/Green2022supplementarytable1.csv") # This is the supplementary table 1 from Harry's LUTS paper, converted to .csv format using excel. Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC9553867/
Green2022supplementarytable1 <- read.csv("Green2022supplementarytable1.csv")


#################################################################
# Step 1: Fetch Prostate Cancer symptom records from GP records #
#################################################################

PrCa_symptom_codes <- Green2022supplementarytable1$read_3 |> as.character()
print(PrCa_symptom_codes)

## Using these codes, select corresponding GP records

# Method 1) Exactly as codes are written, all full stops and ellipses included

PrCaSymptom_gp_records1 <- read_GP(c(
  "1A27.", "K16y8", "XaNFc", "X30Ni", "1A1Z.", "1A11.", "1A1..", "R084.", "1A12.", "R084z", "R0840",
  "1A1..", "1A1..", "1A34.", "1A34.", "R083z", "R083.", "1AZ6.", "1AZ60", "1AZ61", "1AZ62", "1A2..",
  "1A2Z.", "R086z", "8D7..", "R08..", "R08zz", "66K3.", "1A4..", "1A...", "1A...", "1AZ..", "1AZZ.",
  "1AH1.", "R086.", "R08z.", "Kz...", "Ryu4.", "XaB9O", "XaXHi", "XaXHj", "XaXHk", "Xa96j", "1A13.",
  "R0842", "1A33.", "1A31.", "1A3..", "1A3Z.", "1A3..", "R0861", "R0863", "R0860", "317C.", "1A37.",
  "1A36.", "XaD2w", "X77SF", "X76Y0", "R082.", "R0824", "1A32.", "K196.", "R0820", "1A32.", "R0822",
  "1A25.", "R0862", "1A25.", "R15y0", "B7C20", "14270", "ZV104", "B834.", "1J08.", "B58y5", "B8340",
  "B46..", "XaC0j", "Xa3fu", "XaXGk", "XaFwo", "XaKyV"
))

PrCaSymptom_gp_records1 <- PrCaSymptom_gp_records1 %>%                                          ## Change all blanks to NA
  mutate(across(all_of(c("read_2", "read_3")),
                ~ na_if(str_trim(as.character(.)), "")))


PrCaSymptom_gp_records1$read_code <- coalesce(PrCaSymptom_gp_records1$read_2, PrCaSymptom_gp_records1$read_3)



PrCaSymptom_gp_records1_over40 <- PrCaSymptom_gp_records1 %>%
  filter(event_age >= 40) %>%
  mutate(category = case_when(
    read_code %in% c("1A27.", "K16y8", "XaNFc", "X30Ni") ~ "DoubleVoiding",
    read_code %in% c("1A1Z.", "1A11.", "1A1..", "R084.", "1A12.", "R084z", "R0840", "1A1..", "1A1..") ~ "Frequency",
    read_code %in% c("1A34.", "1A34.") ~ "Hesitancy",
    read_code %in% c("R083z", "R083.") ~ "Incontinence",
    read_code %in% c("1AZ6.", "1AZ60", "1AZ61", "1AZ62", "1A2..", "1A2Z.", "R086z", "8D7..", "R08..", "R08zz", "66K3.", "1A4..", "1A...", "1A...", "1AZ..", "1AZZ.", "1AH1.", "R086.", "R08z.", "Kz...", "Ryu4.", "XaB9O", "XaXHi", "XaXHj", "XaXHk", "Xa96j") ~ "LUTS",
    read_code %in% c("1A13.", "R0842") ~ "Nocturia",
    read_code %in% c("1A33.", "1A31.", "1A3..", "1A3Z.", "1A3..", "R0861", "R0863", "R0860", "317C.", "1A37.", "1A36.", "XaD2w", "X77SF", "X76Y0") ~ "PoorStream",
    read_code %in% c("R082.", "R0824", "1A32.", "K196.", "R0820", "1A32.", "R0822") ~ "Retention",
    read_code %in% c("1A25.", "R0862", "1A25.") ~ "Urgency",
    read_code %in% c("R15y0", "B7C20", "14270", "ZV104", "B834.", "1J08.", "B58y5", "B8340", "B46..", "XaC0j", "Xa3fu", "XaXGk", "XaFwo", "XaKyV") ~ "ProstateCancer",
    TRUE ~ NA_character_
  ))




# Method 2) Codes as written, but without full stops and ellipses

PrCaSymptom_gp_records2 <- read_GP(c(
  "1A27", "K16y8", "XaNFc", "X30Ni", "1A1Z", "1A11", "1A1", "R084", "1A12", "R084z", "R0840",
  "1A1", "1A1", "1A34", "1A34", "R083z", "R083", "1AZ6", "1AZ60", "1AZ61", "1AZ62", "1A2",
  "1A2Z", "R086z", "8D7", "R08", "R08zz", "66K3", "1A4", "1A", "1A", "1AZ", "1AZZ",
  "1AH1", "R086", "R08z", "Kz", "Ryu4", "XaB9O", "XaXHi", "XaXHj", "XaXHk", "Xa96j", "1A13",
  "R0842", "1A33", "1A31", "1A3", "1A3Z", "1A3", "R0861", "R0863", "R0860", "317C", "1A37",
  "1A36", "XaD2w", "X77SF", "X76Y0", "R082", "R0824", "1A32", "K196", "R0820", "1A32", "R0822",
  "1A25", "R0862", "1A25", "R15y0", "B7C20", "14270", "ZV104", "B834", "1J08", "B58y5", "B8340",
  "B46", "XaC0j", "Xa3fu", "XaXGk", "XaFwo", "XaKyV"
))

PrCaSymptom_gp_records2 <- PrCaSymptom_gp_records2 %>%                                          ## Change all blanks to NA
  mutate(across(all_of(c("read_2", "read_3")),
                ~ na_if(str_trim(as.character(.)), "")))


PrCaSymptom_gp_records2$read_code <- coalesce(PrCaSymptom_gp_records2$read_2, PrCaSymptom_gp_records2$read_3)



starts_with_any <- function(x, prefixes) {
  # returns a logical vector same length as x
  Reduce(`|`, lapply(prefixes, function(p) startsWith(x, p)))
}


PrCaSymptom_gp_records2_over40 <- PrCaSymptom_gp_records2 %>%
  filter(event_age >= 40) %>%
  mutate(category = case_when(
    starts_with_any(read_code, c("1A27", "K16y8", "XaNFc", "X30Ni")) ~ "DoubleVoiding",
    starts_with_any(read_code, c("1A1Z", "1A11", "1A1", "R084", "1A12", "R084z", "R0840", "1A1", "1A1")) ~ "Frequency",
    starts_with_any(read_code, c("1A34", "1A34")) ~ "Hesitancy",
    starts_with_any(read_code, c("R083z", "R083")) ~ "Incontinence",
    starts_with_any(read_code, c("1AZ6", "1AZ60", "1AZ61", "1AZ62", "1A2", "1A2Z", "R086z", "8D7",
                                 "R08", "R08zz", "66K3", "1A4", "1A", "1A", "1AZ", "1AZZ", "1AH1",
                                 "R086", "R08z", "Kz…", "Ryu4", "XaB9O", "XaXHi", "XaXHj", "XaXHk", "Xa96j")) ~ "LUTS",
    starts_with_any(read_code, c("1A13", "R0842")) ~ "Nocturia",
    starts_with_any(read_code, c("1A33", "1A31", "1A3", "1A3Z", "1A3", "R0861", "R0863", "R0860",
                                 "317C", "1A37", "1A36", "XaD2w", "X77SF", "X76Y0")) ~ "PoorStream",
    starts_with_any(read_code, c("R082", "R0824", "1A32", "K196", "R0820", "1A32", "R0822")) ~ "Retention",
    starts_with_any(read_code, c("1A25", "R0862", "1A25")) ~ "Urgency",
    starts_with_any(read_code, c("R15y0", "B7C20", "14270", "ZV104", "B834", "1J08", "B58y5", "B8340",
                                 "B46", "XaC0j", "Xa3fu", "XaXGk", "XaFwo", "XaKyV")) ~ "ProstateCancer",
    TRUE ~ NA_character_
  ))



## Chose which one you want to use here

chosen_dataframe <- PrCaSymptom_gp_records1_over40 %>%
  dplyr::filter(!is.na(category)) %>%
  dplyr::select(c("eid", "event_dt", "read_code", "event_age", "category"))

## Create new binary symptom columns - All symptoms, Frequency, Hesitancy, Incontinenece, LUTS, Nocturia, PoorStream, Retention, Urgency

PrCa_symptoms <- chosen_dataframe %>%
  dplyr::mutate(
    All_symptoms = if_else(!is.na(category) & category != "ProstateCancer", 1L, 0L, missing = 0L),
    Frequency = if_else(category == "Frequency", 1L, 0L, missing = 0L),
    Hesitancy = if_else(category == "Hesitancy", 1L, 0L, missing = 0L),
    Incontinence = if_else(category == "Incontinence", 1L, 0L, missing = 0L),
    LUTS = if_else(category == "LUTS", 1L, 0L, missing = 0L),
    Nocturia = if_else(category == "Nocturia", 1L, 0L, missing = 0L),
    PoorStream = if_else(category == "PoorStream", 1L, 0L, missing = 0L),
    Retention = if_else(category == "Retention", 1L, 0L, missing = 0L),
    Urgency = if_else(category == "Urgency", 1L, 0L, missing = 0L)
  )

## Create individual dataframes that collapse into earliest symptom presentation

PrCa_all_symptoms_earliest <- PrCa_symptoms %>%
  dplyr::filter(All_symptoms == 1) %>%
  dplyr::select(c("eid", "All_symptoms", "event_dt")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_All_symptoms = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    All_symptoms = first(All_symptoms),
    .groups = "drop"
)
  

PrCa_Frequency_earliest <- PrCa_symptoms %>%
  dplyr::filter(Frequency == 1) %>%
  dplyr::select(c("eid", "Frequency", "event_dt")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_Frequency = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    Frequency = first(Frequency),
    .groups = "drop"
  )

PrCa_Hesitancy_earliest <- PrCa_symptoms %>%
  dplyr::filter(Hesitancy == 1) %>%
  dplyr::select(c("eid", "Hesitancy", "event_dt")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_Hesitancy = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    Hesitancy = first(Hesitancy),
    .groups = "drop"
  )

PrCa_Incontinence_earliest <- PrCa_symptoms %>%
  dplyr::filter(Incontinence == 1) %>%
  dplyr::select(c("eid", "Incontinence", "event_dt")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_Incontinence = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    Incontinence = first(Incontinence),
    .groups = "drop"
  )

PrCa_LUTS_earliest <- PrCa_symptoms %>%
  dplyr::filter(LUTS == 1) %>%
  dplyr::select(c("eid", "LUTS", "event_dt")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_LUTS = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    LUTS = first(LUTS),
    .groups = "drop"
  )

PrCa_Nocturia_earliest <- PrCa_symptoms %>%
  dplyr::filter(Nocturia == 1) %>%
  dplyr::select(c("eid", "Nocturia", "event_dt")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_Nocturia = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    Nocturia = first(Nocturia),
    .groups = "drop"
  )

## To be continued



################################################################################
# Step 2: Fetch Prostate Cancer Diagnosis records from HES and Cancer Registry #
################################################################################


PCaCases_HES<-read_ICD10('C61') # Creates a dataframe of *almost* all recorded prostate cancer diagnoses in HES records

PCaCases_ICD9<-read_ICD9(185)
PCaCases_ICD9 <- PCaCases_ICD9 %>%
  mutate(
    diag_icd10 = case_when(
      diag_icd9 %in% c("1859") ~ "C61",
      TRUE ~ NA_character_
    )
  )
PCaCases_ICD9 <- PCaCases_ICD9 %>%
  select(c("eid", "diag_icd10", "assess_date_initial", "epistart", "epiend", "event_age"))

PCaCases_HES <- bind_rows(PCaCases_HES, PCaCases_ICD9) # Now we have all, after adding old ICD9 diagnoses

PCaCases_cancerregistry<-read_cancer('C61') # Creates a dataframe of all recorded prostate cancer diagnoses in cancer registry records
PCaCases_death<-read_death('C61') # Creates a dataframe of all recorded prostate cancer deaths in death records

PCaCases_HES_earliest <- PCaCases_HES %>%
  dplyr::mutate(epistart = as.Date(epistart)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    epistart = if (all(is.na(epistart))) NA else min(epistart, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    assess_date_initial = first(assess_date_initial),
    diag_icd10 = first(diag_icd10),
    .groups = "drop"
  )

PCaCases_cancerregistry_earliest <- PCaCases_cancerregistry %>%
  dplyr::mutate(date = as.Date(date)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    date = if (all(is.na(date))) NA else min(date, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    assess_date_initial = first(assess_date_initial),
    ICD10 = first(ICD10),
    .groups = "drop"
  )

PCaCases_death_earliest <- PCaCases_death %>%
  dplyr::mutate(date_of_death = as.Date(date_of_death)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    date_of_death = if (all(is.na(date_of_death))) NA else min(date_of_death, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    assess_date_initial = first(assess_date_initial),
    cause_icd10 = first(cause_icd10),
    .groups = "drop"
  )


# Sanity check - how many rows per patient in PCaCases vs. PCaCases_earliest? (should have same n number)

PCaCases_HES %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

PCaCases_HES_earliest %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

PCaCases_cancerregistry %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

PCaCases_cancerregistry_earliest %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

PCaCases_death %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

PCaCases_death_earliest %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

# Merge all "earliest" datasets into one

PCaCases_earliest <- merge(PCaCases_HES_earliest, PCaCases_cancerregistry_earliest, by="eid", all = T)
PCaCases_earliest <- merge(PCaCases_earliest, PCaCases_death_earliest, by="eid", all = T)
PCaCases_earliest <- PCaCases_earliest %>%
  dplyr::rename("assess_date_initial_ICD10" = "assess_date_initial.x", "assess_date_initial_cr" = "assess_date_initial.y", "assess_date_initial_death" = "assess_date_initial",
                "event_age_ICD10" = "event_age.x", "event_age_cr" = "event_age.y", "event_age_death" = "event_age")


# Sanity check - are the assessment centre dates the same in both HES and cancer registry, unless NA on either side?

mismatched_assessdates <- with(PCaCases_earliest,
                               !is.na(assess_date_initial_ICD10) &
                                 !is.na(assess_date_initial_cr) &
                                 assess_date_initial_ICD10 != assess_date_initial_cr
)

any(mismatched_assessdates)

PCaCases_earliest$assess_date_initial <- coalesce(PCaCases_earliest$assess_date_initial_ICD10, PCaCases_earliest$assess_date_initial_cr, PCaCases_earliest$assess_date_initial_death)

# Tag patients who already had prostate cancer at assessment centre as "pre-diagnosed". Also make "earliest PrCa diagnosis date" which is the earliest of epistart and date

PCaCases_prediagnosis <- PCaCases_earliest %>%
  dplyr::mutate(
    pre_diagnosed =
      (is.na(epistart) & is.na(date)) |
      (!is.na(epistart) & epistart <= assess_date_initial) |
      (!is.na(date) & date <= assess_date_initial),
    earliest_PrCa_date = pmin(epistart, date, na.rm = TRUE)
  )

Earliest_PrCa_diagnosis <- PCaCases_prediagnosis %>%
  select('eid', 'earliest_PrCa_date')

