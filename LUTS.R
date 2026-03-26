#############################################################################
# Symptomatic Angle - how well can GRS predict Prostate cancer diagnosis in #
# ------------------- symptomatic participants? ----------------------------#
#############################################################################

## Note: this script requires first running "functions.r"

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

install.packages("RMySQL")
library(RMySQL)
library(dplyr)
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
dxdownload("Callum/Derived_datasets/AncestryProbability1.csv")	# Genetic similarity probabilities for Ancestry group calculated using AncestryProbability1_bigsnpr folder scripts
dxdownload("Callum/Derived_datasets/AncestryProbability2.csv") 	# Genetic similarity probabilities for Ancestry group calculated using AncestryProbability2_plink folder script
dxdownload("Callum/Derived_datasets/FH_PrCa_BrCa.csv")		# Dataset of Family History of Prostate Cancer and Breast Cancer, created using "Family_history.R" script

dxdownload("Callum/GRSs/Conti_multi_ethnic.pgs.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "Multiethnic Analysis"
dxdownload("Callum/GRSs/Conti_European.pgs.tsv") 		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "European...16"
dxdownload("Callum/GRSs/Conti_African.pgs.tsv")			# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "African...19"
dxdownload("Callum/GRSs/Conti_East_Asian.pgs.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "East Asian...22"
dxdownload("Callum/GRSs/Conti_Hispanic.pgs.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "Hispanic...25"

dxdownload("Callum/GRSs/Conti_multiethnicGRS_267.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_MULTI" and "source" set to "Conti"
dxdownload("Callum/GRSs/Conti_EuropeanGRS_265.tsv") 		# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_EUR" and "source" set to "Conti"
dxdownload("Callum/GRSs/Conti_AfricanGRS_246.tsv")			# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_AFR" and "source" set to "Conti"
dxdownload("Callum/GRSs/Conti_East_AsianGRS_222.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_EAS" and "source" set to "Conti"
dxdownload("Callum/GRSs/Conti_HispanicGRS_253.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_HIS" and "source" set to "Conti"

dxdownload("Callum/GRSs/Wang_multi_ethnic.pgs.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/Wang_European.pgs.tsv") 		# Calculated Wang GRS for all eligible participants using "Conti_script" 
dxdownload("Callum/GRSs/Wang_African.pgs.tsv")			# Calculated Wang GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/Wang_East_Asian.pgs.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/Wang_Hispanic.pgs.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_script"

dxdownload("Callum/GRSs/Wang_multiethnicGRS_450.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_MULTI" and "source" set to "Wang"
dxdownload("Callum/GRSs/Wang_EuropeanGRS_445.tsv") 		# Calculated Wang GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_EUR" and "source" set to "Wang"
dxdownload("Callum/GRSs/Wang_AfricanGRS_444.tsv")			# Calculated Wang GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_AFR" and "source" set to "Wang"
dxdownload("Callum/GRSs/Wang_East_AsianGRS_379.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_EAS" and "source" set to "Wang"
dxdownload("Callum/GRSs/Wang_HispanicGRS_446.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_HIS" and "source" set to "Wang"

dxdownload("Callum/GRSs/OR_adjustedGRS_267.tsv")	# Calculated using "Create_adjustedGRS_weighting_Conti_odds_ratios_by_AncestryProbability2.R" script

dxdownload("Callum/GRSs/Schumacher.pgs.tsv")    # Calculated Schumacher (2018) GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/SchumacherGRS_145.tsv")    # Calculated Schumacher (2018) GRS for all eligible participants using "Conti_GRS_267.R" with "source" set to "Schumacher"

dxdownload("Callum/GRSs/BARCODE1.pgs.tsv")    # Calculated BARCODE1 (2021) GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/BARCODE1GRS_129.tsv")    # Calculated BARCODE1 (2021) GRS for all eligible participants using "Conti_GRS_267.R" with "source" set to "BARCODE1"

dxdownload("Callum/GRSs/Seibert.pgs.tsv")    # Calculated Seibert (2018) GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/SeibertGRS_52.tsv")    # Calculated Seibert (2018) GRS for all eligible participants using "Conti_GRS_267.R" with "source" set to "Seibert"

dxdownload("Callum/GRSs/Pagadala.pgs.tsv")    # Calculated Pagadala (2022) GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/PagadalaGRS_285.tsv")    # Calculated Pagadala (2022) GRS for all eligible participants using "Conti_GRS_267.R" with "source" set to "Pagadala"

system("dx download Callum/LUTS/Green2022supplementarytable1.csv") # This is the supplementary table 1 from Harry Green's LUTS paper, converted to .csv format using excel. Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC9553867/
Green2022supplementarytable1 <- read.csv("Green2022supplementarytable1.csv")


#################################################################
# Step 1: Fetch Prostate Cancer symptom records from GP records #
#################################################################

PrCa_symptom_codes <- Green2022supplementarytable1$read_3 |> as.character()
print(PrCa_symptom_codes)

## Using these codes, select corresponding GP records

# Method 1) Exactly as codes are written, all full stops and ellipses (written as three ASCII dots) included, without Prostate Cancer codes

PrCaSymptom_gp_records1 <- read_GP(c(
  "1A27.", "K16y8", "XaNFc", "X30Ni", "1A1Z.", "1A11.", "1A1..", "R084.", "1A12.", "R084z", "R0840",
  "1A1..", "1A1..", "1A34.", "1A34.", "R083z", "R083.", "1AZ6.", "1AZ60", "1AZ61", "1AZ62", "1A2..",
  "1A2Z.", "R086z", "8D7..", "R08..", "R08zz", "66K3.", "1A4..", "1A...", "1A...", "1AZ..", "1AZZ.",
  "1AH1.", "R086.", "R08z.", "Kz...", "Ryu4.", "XaB9O", "XaXHi", "XaXHj", "XaXHk", "Xa96j", "1A13.",
  "R0842", "1A33.", "1A31.", "1A3..", "1A3Z.", "1A3..", "R0861", "R0863", "R0860", "317C.", "1A37.",
  "1A36.", "XaD2w", "X77SF", "X76Y0", "R082.", "R0824", "1A32.", "K196.", "R0820", "1A32.", "R0822",
  "1A25.", "R0862", "1A25.")
)

##### Prostate Cancer Codes
##  "R15y0", "B7C20", "14270", "ZV104", "B834.", "1J08.", "B58y5", "B8340",
##  "B46..", "XaC0j", "Xa3fu", "XaXGk", "XaFwo", "XaKyV"


PrCaSymptom_gp_records1 <- PrCaSymptom_gp_records1 %>%                                          ## Change all blanks to NA
  dplyr::mutate(across(all_of(c("read_2", "read_3")),
                       ~ na_if(str_trim(as.character(.)), "")))


PrCaSymptom_gp_records1$read_code <- coalesce(PrCaSymptom_gp_records1$read_2, PrCaSymptom_gp_records1$read_3) # Combine read2 and read3 into one column called "read_code"

## Assign symptom categories (note that this is strict - only exact matches will be categorized)

PrCaSymptom_gp_records1_over40_strict <- PrCaSymptom_gp_records1 %>%
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
    ##    read_code %in% c("R15y0", "B7C20", "14270", "ZV104", "B834.", "1J08.", "B58y5", "B8340", "B46..", "XaC0j", "Xa3fu", "XaXGk", "XaFwo", "XaKyV") ~ "ProstateCancer",
    TRUE ~ NA_character_
  ))


## Write uncategorized dataframe to investigate non-exact codes - see if we should drop or include them

uncategorized <- PrCaSymptom_gp_records1_over40_strict %>%
  dplyr::filter(is.na(category))


unique_uncategorized_codes <- uncategorized %>%
  mutate(read_code = str_trim(as.character(read_code))) %>%
  filter(!is.na(read_code), read_code != "") %>%
  distinct(read_code) %>%
  arrange(read_code)

unique_uncategorized_codes

## Read in read2 and read3 dictionaries to see what these codes mean

dxdownload("Callum/read2_read3_codes/read2_lkp3.csv")
dxdownload("Callum/read2_read3_codes/read3_lkp3.csv")

uncategorized_codes_table <- uncategorized %>%
  transmute(code = str_trim(as.character(read_code))) %>%
  filter(!is.na(code), code != "") %>%
  distinct()


read2_dictionary <- read.csv("read2_lkp3.csv") %>%
  transmute(code = as.character(read_2), term_r2 = as.character(term_description)) %>%
  distinct()


read3_dictionary <- read.csv("read3_lkp3.csv") %>%
  transmute(code = as.character(read_3), term_r3 = as.character(term_description)) %>%
  distinct()



mapped_uncategorized_codes_table <- uncategorized_codes_table %>%
  left_join(read2_dictionary, by = "code") %>%
  left_join(read3_dictionary,  by = "code") %>%
  mutate(
    scheme = case_when(
      !is.na(term_r2)   &  is.na(term_r3) ~ "Read v2",
      is.na(term_r2)    & !is.na(term_r3) ~ "Read v3",
      !is.na(term_r2)   & !is.na(term_r3) ~ "Both (check)",
      TRUE                                 ~ "Not found"
    ),
    term = coalesce(term_r2, term_r3)
  ) %>%
  select(code, scheme, term) %>%
  arrange(scheme, code)

mapped_uncategorized_codes_table


## Create a less strict categorization based on mapped_uncategorized_codes_table

PrCaSymptom_gp_records1_over40_liberal <- PrCaSymptom_gp_records1 %>%
  filter(event_age >= 40) %>%
  mutate(category = case_when(
    read_code %in% c("1A27.", "K16y8", "XaNFc", "X30Ni") ~ "DoubleVoiding",
    read_code %in% c("1A1Z.", "1A11.", "1A1..", "R084.", "1A12.", "R084z", "R0840", "1A1..", "1A1..") ~ "Frequency",
    read_code %in% c("1A34.", "1A34.") ~ "Hesitancy",
    read_code %in% c("R083z", "R083.", "1A23.", "1A24.", "1A26.", "1A22.", "1A220", "R0831", "R0832", "R0830") ~ "Incontinence",
    read_code %in% c("1AZ6.", "1AZ60", "1AZ61", "1AZ62", "1A2..", "1A2Z.", "R086z", "8D7..", "R08..", "R08zz", "66K3.", "1A4..", "1A...", "1A...", "1AZ..", "1AZZ.", "1AH1.", "R086.", "R08z.", "Kz...", "Ryu4.", "XaB9O", "XaXHi", "XaXHj", "XaXHk", "Xa96j") ~ "LUTS",
    read_code %in% c("1A13.", "R0842") ~ "Nocturia",
    read_code %in% c("1A33.", "1A31.", "1A3..", "1A3Z.", "1A3..", "R0861", "R0863", "R0860", "317C.", "1A37.", "1A36.", "XaD2w", "X77SF", "X76Y0", "1AZ3.", "Ryu40") ~ "PoorStream",
    read_code %in% c("R082.", "R0824", "1A32.", "K196.", "R0820", "1A32.", "R0822", "1AH0.", "R0823") ~ "Retention",
    read_code %in% c("1A25.", "R0862", "1A25.", "1A35.", "R08z2") ~ "Urgency",
    ##    read_code %in% c("R15y0", "B7C20", "14270", "ZV104", "B834.", "1J08.", "B58y5", "B8340", "B46..", "XaC0j", "Xa3fu", "XaXGk", "XaFwo", "XaKyV") ~ "ProstateCancer",
    TRUE ~ NA_character_
  ))

## Choose which one you want to use here

chosen_dataframe <- PrCaSymptom_gp_records1_over40_liberal %>%  ## Change this to toggle between strict inclusion and liberal inclusion
  dplyr::filter(!is.na(category)) %>%    
  dplyr::select(c("eid", "event_dt", "read_code", "event_age", "date_of_birth", "category")) %>%
  inner_join(baseline) %>%
  dplyr::filter(sex == "Male")

## Create new binary symptom columns - All symptoms, Frequency, Hesitancy, Incontinenece, LUTS, Nocturia, PoorStream, Retention, Urgency

PrCa_symptoms <- chosen_dataframe %>%
  dplyr::mutate(
    All_symptoms = 1,
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
  dplyr::select(c("eid", "All_symptoms", "event_dt", "event_age", "date_of_birth")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_symptom = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    All_symptoms = first(All_symptoms),
    date_of_birth = first(date_of_birth),
    .groups = "drop"
  )


PrCa_Frequency_earliest <- PrCa_symptoms %>%
  dplyr::filter(Frequency == 1) %>%
  dplyr::select(c("eid", "Frequency", "event_dt", "event_age", "date_of_birth")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_symptom = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    Frequency = first(Frequency),
    date_of_birth = first(date_of_birth),
    .groups = "drop"
  )

PrCa_Hesitancy_earliest <- PrCa_symptoms %>%
  dplyr::filter(Hesitancy == 1) %>%
  dplyr::select(c("eid", "Hesitancy", "event_dt", "event_age", "date_of_birth")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_symptom = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    Hesitancy = first(Hesitancy),
    date_of_birth = first(date_of_birth),
    .groups = "drop"
  )

PrCa_Incontinence_earliest <- PrCa_symptoms %>%
  dplyr::filter(Incontinence == 1) %>%
  dplyr::select(c("eid", "Incontinence", "event_dt", "event_age", "date_of_birth")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_symptom = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    Incontinence = first(Incontinence),
    date_of_birth = first(date_of_birth),
    .groups = "drop"
  )

PrCa_LUTS_earliest <- PrCa_symptoms %>%
  dplyr::filter(LUTS == 1) %>%
  dplyr::select(c("eid", "LUTS", "event_dt", "event_age", "date_of_birth")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_symptom = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    LUTS = first(LUTS),
    date_of_birth = first(date_of_birth),
    .groups = "drop"
  )

PrCa_Nocturia_earliest <- PrCa_symptoms %>%
  dplyr::filter(Nocturia == 1) %>%
  dplyr::select(c("eid", "Nocturia", "event_dt", "event_age", "date_of_birth")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_symptom = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    Nocturia = first(Nocturia),
    date_of_birth = first(date_of_birth),
    .groups = "drop"
  )

PrCa_PoorStream_earliest <- PrCa_symptoms %>%
  dplyr::filter(PoorStream == 1) %>%
  dplyr::select(c("eid", "PoorStream", "event_dt", "event_age", "date_of_birth")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_symptom = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    PoorStream = first(PoorStream),
    date_of_birth = first(date_of_birth),
    .groups = "drop"
  )

PrCa_Urgency_earliest <- PrCa_symptoms %>%
  dplyr::filter(Urgency == 1) %>%
  dplyr::select(c("eid", "Urgency", "event_dt", "event_age", "date_of_birth")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_symptom = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    Urgency = first(Urgency),
    date_of_birth = first(date_of_birth),
    .groups = "drop"
  )

PrCa_Retention_earliest <- PrCa_symptoms %>%
  dplyr::filter(Retention == 1) %>%
  dplyr::select(c("eid", "Retention", "event_dt", "event_age", "date_of_birth")) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    EarliestDate_symptom = if (all(is.na(event_dt))) NA else min(event_dt, na.rm = TRUE),
    event_age = if (all(is.na(event_age))) NA else min(event_age, na.rm = TRUE),
    Retention = first(Retention),
    date_of_birth = first(date_of_birth),
    .groups = "drop"
  )

chosen_symptom_earliest_dataframe <- PrCa_all_symptoms_earliest # change this to look at different symptom group in rest of analysis

## Sanity check 1 - is there truly 1 row per participant now?

test_df <- chosen_symptom_earliest_dataframe 

test_df %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )                                   # Should be same as n now, with 1 for every other number

## Sanity check 2 - does "event_age" always match EarliestDate - date_of_birth?

test_df <- chosen_symptom_earliest_dataframe %>%
  dplyr::mutate(test_age = as.numeric(EarliestDate_symptom - date_of_birth)/365.25)

test_df %>%
  dplyr::filter(
    test_age == event_age
  ) %>%
  summarise(n = n())  ## should be same as number of obs. in test.df


################################################################################
# Step 2: Fetch Prostate Cancer Diagnosis records from HES and Cancer Registry #
################################################################################


PCaCases_HES<-read_ICD10(c('C61', 'Z8546', 'R9721')) # Creates a dataframe of *almost* all recorded prostate cancer diagnoses in HES records
# Z85.46 is "Personal history of malignant neoplasm of prostate" and R97.21 is "Elevated prostate specific antigen (PSA)". Adding these to read_ICD10 does not add any cases

PCaCases_ICD9<-read_ICD9(c(185, 'V1046'))
# V10.46 is "Personal history of malignant neoplasm of prostate". Adding this to read_ICD9 does not add any cases
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
                "event_age_ICD10" = "event_age.x", "event_age_cr" = "event_age.y", "event_age_death" = "event_age",
                "icd10_HES" = "diag_icd10", "icd10_cr" = "ICD10", "icd10_death" = "cause_icd10")


# Sanity check - are the assessment centre dates the same in both HES and cancer registry, unless NA on either side?

mismatched_assessdates <- with(PCaCases_earliest,
                               !is.na(assess_date_initial_ICD10) &
                                 !is.na(assess_date_initial_cr) &
                                 assess_date_initial_ICD10 != assess_date_initial_cr
)

any(mismatched_assessdates)

PCaCases_earliest$assess_date_initial <- coalesce(PCaCases_earliest$assess_date_initial_ICD10, PCaCases_earliest$assess_date_initial_cr, PCaCases_earliest$assess_date_initial_death)

## Combine PrCa diagnosis dataframe with chosen symptom dataframe

PrCa_symptoms_diagnosis <- merge(chosen_symptom_earliest_dataframe, PCaCases_earliest, by = "eid", all.x = T)

# Tag participants who already had prostate cancer diagnosis before symptom event date as "pre-diagnosed". Also make "earliest PrCa diagnosis date" which is the earliest of epistart and date

PrCa_symptoms_diagnosis <- PrCa_symptoms_diagnosis %>%
  dplyr::mutate(
    pre_diagnosed =
      (!is.na(epistart) | !is.na(date) | !is.na(date_of_death)) &
      ((!is.na(epistart) & epistart < EarliestDate_symptom) | (!is.na(date) & date < EarliestDate_symptom)),
    earliest_PrCa_date = pmin(epistart, date, na.rm = TRUE)
  )

Earliest_PrCa_diagnosis <- PrCa_symptoms_diagnosis %>%
  select('eid', 'earliest_PrCa_date')

# Tag participants who had a prostate cancer record in HES but not in cancer registry or death records as "HES-only" (these will be excluded later)

PrCa_symptoms_diagnosis <- PrCa_symptoms_diagnosis %>%
  dplyr::mutate(HES_only = !is.na(icd10_HES) & is.na(icd10_cr) & is.na(icd10_death))






##################################################################################################
# Step 3 - Fetch cancer death, chemotherapy, surgery, radiotherapy, and androgen therapy records #
##################################################################################################

###### Severe = chemotherapy or death from cancer within 2 years of diagnosis

###### Actionable = radiotherapy, prostate surgery, or anti-androgen treatments 
###### within 2 years of diagnosis, or any of the Severe criteria

Cancer_death<-read_death('C')
Cancer_death=Cancer_death[grepl("^C", Cancer_death$cause_icd10),] # Creates a dataframe of all recorded cancer (general) deaths in death records
Cancer_death_earliest <- Cancer_death %>%
  dplyr::mutate(date_of_death = as.Date(date_of_death)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    date_of_death = if (all(is.na(date_of_death))) NA else min(date_of_death, na.rm = TRUE),
    assess_date_initial = first(assess_date_initial),
    .groups = "drop"
  )

# Dataframe for chemotherapy (IV, IM, unspecified),

chemotherapy  <- read_OPCS(c('X70', 'X71', 'X72'))
chemotherapy <- merge(chemotherapy, Earliest_PrCa_diagnosis, by = 'eid') 
chemotherapy <- chemotherapy %>%
  dplyr::filter(opdate >= earliest_PrCa_date)
chemotherapy_earliest <- chemotherapy %>%
  dplyr::mutate(opdate = as.Date(opdate)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    opdate = if (all(is.na(opdate))) NA else min(opdate, na.rm = TRUE),
    assess_date_initial = first(assess_date_initial),
    .groups = "drop"
  )

death_chemo <- merge(chemotherapy_earliest, Cancer_death_earliest, by = "eid", all = T)

# Check that all assessment centre dates match, then merge

mismatched_assessdates <- with(death_chemo,
                               !is.na(assess_date_initial.x) &
                                 !is.na(assess_date_initial.y) &
                                 assess_date_initial.x != assess_date_initial.y
)
any(mismatched_assessdates)

death_chemo$assess_date_initial <- coalesce(death_chemo$assess_date_initial.x, death_chemo$assess_date_initial.y)

death_chemo <- death_chemo %>%
  select('eid', 'chemo_opdate' = 'opdate', 'date_of_cancer_death' = 'date_of_death')




# Dataframe for prostatectomy

surgery <- read_OPCS(c('M61', 'M611', 'M612', 'M613'))
surgery <- merge(surgery, Earliest_PrCa_diagnosis, by = 'eid')
surgery <- surgery %>%
  dplyr::filter(opdate >= earliest_PrCa_date)
surgery_earliest <- surgery %>%
  dplyr::mutate(opdate = as.Date(opdate)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    opdate = if (all(is.na(opdate))) NA else min(opdate, na.rm = TRUE),
    .groups = "drop"
  )

# Dataframe for radiotherapy (external, brachytherapy, planning, unspecified)

radiotherapy <- read_OPCS(c('X65', 'X66', 'X67', 'X69'))
radiotherapy <- merge(radiotherapy, Earliest_PrCa_diagnosis, by = 'eid')
radiotherapy <- radiotherapy %>%
  dplyr::filter(opdate >= earliest_PrCa_date)
radiotherapy_earliest <- radiotherapy %>%
  dplyr::mutate(opdate = as.Date(opdate)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    opdate = if (all(is.na(opdate))) NA else min(opdate, na.rm = TRUE),
    .groups = "drop"
  )

# Dataframe for androgen (therapy?)

androgen <- read_OPCS(c('X741', 'X383', 'S525', 'S526'))
androgen <- merge(androgen, Earliest_PrCa_diagnosis, by = 'eid')
androgen <- androgen %>%
  dplyr::filter(opdate >= earliest_PrCa_date)
androgen_earliest <- androgen %>%
  dplyr::mutate(opdate = as.Date(opdate)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    opdate = if (all(is.na(opdate))) NA else min(opdate, na.rm = TRUE),
    .groups = "drop"
  )


# Sanity check - how many rows per patient in each dataframe vs. dataframe_earliest?

Cancer_death %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

Cancer_death_earliest %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

chemotherapy %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

chemotherapy_earliest %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

surgery %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

surgery_earliest %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

radiotherapy %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

radiotherapy_earliest %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

androgen %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

androgen_earliest %>%
  dplyr::count(eid) %>%
  dplyr::summarise(
    avg_rows_per_id = mean(n),
    median_rows_per_id = median(n),
    min_rows = min(n),
    max_rows = max(n),
    sample_size = dplyr::n()
  )

# Merge all actionable criteria (remember - "actionable" includes all severe criteria also)

surgery_radio_androgen <- merge(surgery_earliest, radiotherapy_earliest, by = "eid", all = T)
surgery_radio_androgen <- merge(surgery_radio_androgen, androgen_earliest, by = "eid", all = T)

surgery_radio_androgen <- surgery_radio_androgen %>%
  rename('surgery_opdate' = 'opdate.x', 'radio_opdate' = 'opdate.y', 'androgen_opdate' = 'opdate')

actionable_criteria <- merge(death_chemo, surgery_radio_androgen, by = "eid", all = T)

# Collect severity criteria

severity_criteria <- death_chemo




######################################################################################
# Step 4 - Load in independent variables, GRSs, and covariates from DNAnexus project #
######################################################################################

# Load in independent variable, covariates, and GRSs

iv <- read_csv("imputed_rs72725854.csv") # loads independent variable dataset
iv <- iv %>%
  select(c("eid", "rs72725854_G", "rs72725854_T")) 

covariates <- read_csv("Age_Sex_PRS_GA.csv")
covariates <- covariates %>% 
  rename('Age_at_assess_centre' = 'p21022', 'Sex' = 'p31', 'PRS' = 'p26267', 'Enhanced_PRS' = 'p26268', 'Genetic_sex' = 'p22001', 'Genomic_ancestry' = 'p30079')


ethnicity <- read_csv("ethnicity.csv")
ethnicity <- ethnicity %>%
  rename('Ethnicity' = 'p21000_i0')

ethnicity <- ethnicity %>% # "wide net" ethnicity grouping
  mutate(ethnicity_group_wide = case_when(
    Ethnicity %in% c("British", "White", "Any other white background", "Irish") ~ "White",
    Ethnicity %in% c("African", "Caribbean", "Any other Black background", "White and Black African", "White and Black Caribbean", "Black or Black British") ~ "Black",
    Ethnicity %in% c("Pakistani", "Indian", "Bangladeshi") ~ "South Asian",
    Ethnicity %in% c("Any other Asian background", "White and Asian", "Chinese", "Asian or Asian British") ~ "East Asian",
    Ethnicity %in% c("Other ethnic group", "Any other mixed background", "Prefer not to answer", "Do not know", "Mixed", "NA") ~ "Other",
    TRUE ~ NA_character_
  ))

ethnicity <- ethnicity %>% # stricter ethnicity grouping
  mutate(ethnicity_group_narrow = case_when(
    Ethnicity %in% c("British", "White", "Any other white background", "Irish") ~ "White",
    Ethnicity %in% c("African", "Caribbean", "Any other Black background", "Black or Black British") ~ "Black",
    Ethnicity %in% c("Pakistani", "Indian", "Bangladeshi") ~ "South Asian",
    Ethnicity %in% c("Chinese") ~ "East Asian",
    Ethnicity %in% c("White and Black African", "White and Black Caribbean") ~ "Mixed White and Black",
    Ethnicity %in% c("Other ethnic group", "Any other mixed background", "Prefer not to answer", "Do not know", "Mixed", "NA", "Any other Asian background", "White and Asian", "Asian or Asian British") ~ "Other",
    TRUE ~ NA_character_
  ))

principal_components<-read_csv("HGDP_1KG_PCs_UKB2_scaled.csv") %>%
  dplyr::rename("eid" = "IID") %>%
  dplyr::select("eid", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

family_history <- read.csv("FH_PrCa_BrCa.csv") %>%
  dplyr::select(c("eid", "FH_Prostate_cancer", "FH_Breast_cancer", "FH_PrCa_BrCa"))

covariates <- merge(covariates, ethnicity, by = "eid", all.x = T)
covariates <- merge(covariates, principal_components, by = "eid", all = T)
covariates <- merge(covariates, family_history, by = "eid", all.x = T)


## Load in Conti GRSs, rename columns, and create "top 10%" variables for each GRS (where "top 10%" is defined as being above the 90th percentile of the GRS distribution in the entire eligible UKB population, not just within each ancestry group)

ContimultiethnicGRS <- read_delim("Conti_multi_ethnic.pgs.tsv")
ContimultiethnicGRS <- ContimultiethnicGRS %>%
  dplyr::rename('ContimultiethnicGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_ContimultiethnicGRS = ContimultiethnicGRS >= quantile(ContimultiethnicGRS, probs = 0.9))

ContiEuropeanGRS <- read_delim("Conti_European.pgs.tsv")
ContiEuropeanGRS <- ContiEuropeanGRS %>%
  dplyr::rename('ContiEuropeanGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_ContiEuropeanGRS = ContiEuropeanGRS >= quantile(ContiEuropeanGRS, probs = 0.9))

ContiAfricanGRS <- read_delim("Conti_African.pgs.tsv")
ContiAfricanGRS <- ContiAfricanGRS %>%
  dplyr::rename('ContiAfricanGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_ContiAfricanGRS = ContiAfricanGRS >= quantile(ContiAfricanGRS, probs = 0.9))

ContiEast_AsianGRS <- read_delim("Conti_East_Asian.pgs.tsv")
ContiEast_AsianGRS <- ContiEast_AsianGRS %>%
  dplyr::rename('ContiEast_AsianGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_ContiEast_AsianGRS = ContiEast_AsianGRS >= quantile(ContiEast_AsianGRS, probs = 0.9))

ContiHispanicGRS <- read_delim("Conti_Hispanic.pgs.tsv")
ContiHispanicGRS <- ContiHispanicGRS %>%
  dplyr::rename('ContiHispanicGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_ContiHispanicGRS = ContiHispanicGRS >= quantile(ContiHispanicGRS, probs = 0.9))

ContimultiethnicGRS267 <- read_delim("Conti_multiethnicGRS_267.tsv")
ContimultiethnicGRS267 <- ContimultiethnicGRS267 %>%
  dplyr::select(c("eid", "Conti_GRS_MULTI_avg")) %>%
  dplyr::rename('ContimultiethnicGRS267' = 'Conti_GRS_MULTI_avg') %>%
  dplyr::mutate(top10_all_ContimultiethnicGRS267 = ContimultiethnicGRS267 >= quantile(ContimultiethnicGRS267, probs = 0.9)) 

ContiEuropeanGRS265 <- read_delim("Conti_EuropeanGRS_265.tsv")
ContiEuropeanGRS265 <- ContiEuropeanGRS265 %>%
  dplyr::select(c("eid", "Conti_GRS_EUR_avg")) %>%
  dplyr::rename('ContiEuropeanGRS265' = 'Conti_GRS_EUR_avg') %>%
  dplyr::mutate(top10_all_ContiEuropeanGRS265 = ContiEuropeanGRS265 >= quantile(ContiEuropeanGRS265, probs = 0.9))

ContiAfricanGRS246 <- read_delim("Conti_AfricanGRS_246.tsv")
ContiAfricanGRS246 <- ContiAfricanGRS246 %>%
  dplyr::select(c("eid", "Conti_GRS_AFR_avg")) %>%
  dplyr::rename('ContiAfricanGRS246' = 'Conti_GRS_AFR_avg') %>%
  dplyr::mutate(top10_all_ContiAfricanGRS246 = ContiAfricanGRS246 >= quantile(ContiAfricanGRS246, probs = 0.9))

ContiEast_AsianGRS222 <- read_delim("Conti_East_AsianGRS_222.tsv")
ContiEast_AsianGRS222 <- ContiEast_AsianGRS222 %>%
  dplyr::select(c("eid", "Conti_GRS_EAS_avg")) %>%
  dplyr::rename('ContiEast_AsianGRS222' = 'Conti_GRS_EAS_avg') %>%
  dplyr::mutate(top10_all_ContiEast_AsianGRS222 = ContiEast_AsianGRS222 >= quantile(ContiEast_AsianGRS222, probs = 0.9))

ContiHispanicGRS253 <- read_delim("Conti_HispanicGRS_253.tsv")
ContiHispanicGRS253 <- ContiHispanicGRS253 %>%
  dplyr::select(c("eid", "Conti_GRS_HIS_avg")) %>%
  dplyr::rename('ContiHispanicGRS253' = 'Conti_GRS_HIS_avg') %>%
  dplyr::mutate(top10_all_ContiHispanicGRS253 = ContiHispanicGRS253 >= quantile(ContiHispanicGRS253, probs = 0.9))

# Merge all data

All_Conti_GRS <- merge(ContimultiethnicGRS, ContiEuropeanGRS, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, ContiAfricanGRS, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, ContiEast_AsianGRS, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, ContiHispanicGRS, by = "eid")

All_Conti_GRS <- merge(All_Conti_GRS, ContimultiethnicGRS267, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, ContiEuropeanGRS265, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, ContiAfricanGRS246, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, ContiEast_AsianGRS222, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, ContiHispanicGRS253, by = "eid")

### Creating an adjusted GRS ###

AncestryProbability1 <- read.csv("AncestryProbability1.csv") %>%          
  dplyr::select(c("IID", "AMR", "AFR", "CSA", "EAS", "EUR", "MID")) %>%
  dplyr::rename("eid" = "IID")

AncestryProbability2 <- read.csv("AncestryProbability2.csv") %>%         
  dplyr::select(c("IID", "AMR", "AFR", "CSA", "EAS", "EUR", "MID")) %>%
  dplyr::rename("eid" = "IID")

GRS_plus_ancestry <- merge(All_Conti_GRS, AncestryProbability2, all=T, by="eid") # Change y to either AncestryProbability1 or AncestryProbability 2 depending on preference
GRS_plus_ancestry <- merge(GRS_plus_ancestry, covariates, by = "eid", all.x=T) %>%
  dplyr::select(c("eid", "ContimultiethnicGRS", "top10_all_ContimultiethnicGRS", "ContiEuropeanGRS", "top10_all_ContiEuropeanGRS", "ContiAfricanGRS", "top10_all_ContiAfricanGRS", "ContiEast_AsianGRS", "top10_all_ContiEast_AsianGRS", "ContiHispanicGRS", "top10_all_ContiHispanicGRS", "ContimultiethnicGRS267", "top10_all_ContimultiethnicGRS267", "ContiEuropeanGRS265", "top10_all_ContiEuropeanGRS265", "ContiAfricanGRS246", "top10_all_ContiAfricanGRS246", "ContiEast_AsianGRS222", "top10_all_ContiEast_AsianGRS222", "ContiHispanicGRS253", "top10_all_ContiHispanicGRS253", "AMR", "AFR", "CSA", "EAS", "EUR", "MID", "Genomic_ancestry"))

GRS_plus_ancestry <- GRS_plus_ancestry %>%
  dplyr::mutate(
    ContiadjustedGRS = (EUR*ContiEuropeanGRS)+(AFR*ContiAfricanGRS)+(EAS*ContiEast_AsianGRS)+(AMR*ContiHispanicGRS)+(CSA*ContimultiethnicGRS)+(MID*ContiAfricanGRS)
  )

All_Conti_GRS <- GRS_plus_ancestry %>%
  dplyr::select(c("eid", "EUR", "AFR", "EAS", "CSA", "MID", "AMR", "ContimultiethnicGRS", "top10_all_ContimultiethnicGRS", "ContiEuropeanGRS", "top10_all_ContiEuropeanGRS", "ContiAfricanGRS", "top10_all_ContiAfricanGRS", "ContiEast_AsianGRS", "top10_all_ContiEast_AsianGRS", "ContiHispanicGRS", "top10_all_ContiHispanicGRS", "ContimultiethnicGRS267", "top10_all_ContimultiethnicGRS267", "ContiEuropeanGRS265", "top10_all_ContiEuropeanGRS265", "ContiAfricanGRS246", "top10_all_ContiAfricanGRS246", "ContiEast_AsianGRS222", "top10_all_ContiEast_AsianGRS222", "ContiHispanicGRS253", "top10_all_ContiHispanicGRS253", "ContiadjustedGRS"))

## Add Odds Ratio-adjusted GRS (where ancestry probability is applied to the OR, not the final scores)

ContiORadjustedGRS <- read.table("OR_adjustedGRS_267.tsv", header = T) %>%
  dplyr::select(c("eid", "Conti_ORfirst_avg")) %>%
  dplyr::rename("ContiORadjustedGRS" = "Conti_ORfirst_avg")

All_Conti_GRS <- merge(All_Conti_GRS, ContiORadjustedGRS, by = "eid", all.x = T)

##################### Now repeat for Wang GRSs

WangmultiethnicGRS <- read_delim("Wang_multi_ethnic.pgs.tsv")
WangmultiethnicGRS <- WangmultiethnicGRS %>%
  dplyr::rename('WangmultiethnicGRS' = 'Wang') %>%
  dplyr::mutate(top10_all_WangmultiethnicGRS = WangmultiethnicGRS >= quantile(WangmultiethnicGRS, probs = 0.9))

WangEuropeanGRS <- read_delim("Wang_European.pgs.tsv")
WangEuropeanGRS <- WangEuropeanGRS %>%
  dplyr::rename('WangEuropeanGRS' = 'Wang') %>%
  dplyr::mutate(top10_all_WangEuropeanGRS = WangEuropeanGRS >= quantile(WangEuropeanGRS, probs = 0.9))

WangAfricanGRS <- read_delim("Wang_African.pgs.tsv")
WangAfricanGRS <- WangAfricanGRS %>%
  dplyr::rename('WangAfricanGRS' = 'Wang') %>%
  dplyr::mutate(top10_all_WangAfricanGRS = WangAfricanGRS >= quantile(WangAfricanGRS, probs = 0.9))

WangEast_AsianGRS <- read_delim("Wang_East_Asian.pgs.tsv")
WangEast_AsianGRS <- WangEast_AsianGRS %>%
  dplyr::rename('WangEast_AsianGRS' = 'Wang') %>%
  dplyr::mutate(top10_all_WangEast_AsianGRS = WangEast_AsianGRS >= quantile(WangEast_AsianGRS, probs = 0.9))

WangHispanicGRS <- read_delim("Wang_Hispanic.pgs.tsv")
WangHispanicGRS <- WangHispanicGRS %>%
  dplyr::rename('WangHispanicGRS' = 'Wang') %>%
  dplyr::mutate(top10_all_WangHispanicGRS = WangHispanicGRS >= quantile(WangHispanicGRS, probs = 0.9))

WangmultiethnicGRS450 <- read_delim("Wang_multiethnicGRS_450.tsv")
WangmultiethnicGRS450 <- WangmultiethnicGRS450 %>%
  dplyr::select(c("eid", "Conti_GRS_MULTI_avg")) %>%
  dplyr::rename('WangmultiethnicGRS450' = 'Conti_GRS_MULTI_avg') %>%
  dplyr::mutate(top10_all_WangmultiethnicGRS450 = WangmultiethnicGRS450 >= quantile(WangmultiethnicGRS450, probs = 0.9))

WangEuropeanGRS445 <- read_delim("Wang_EuropeanGRS_445.tsv")
WangEuropeanGRS445 <- WangEuropeanGRS445 %>%
  dplyr::select(c("eid", "Conti_GRS_EUR_avg")) %>%
  dplyr::rename('WangEuropeanGRS445' = 'Conti_GRS_EUR_avg') %>%
  dplyr::mutate(top10_all_WangEuropeanGRS445 = WangEuropeanGRS445 >= quantile(WangEuropeanGRS445, probs = 0.9))

WangAfricanGRS444 <- read_delim("Wang_AfricanGRS_444.tsv")
WangAfricanGRS444 <- WangAfricanGRS444 %>%
  dplyr::select(c("eid", "Conti_GRS_AFR_avg")) %>%
  dplyr::rename('WangAfricanGRS444' = 'Conti_GRS_AFR_avg') %>%
  dplyr::mutate(top10_all_WangAfricanGRS444 = WangAfricanGRS444 >= quantile(WangAfricanGRS444, probs = 0.9))

WangEast_AsianGRS379 <- read_delim("Wang_East_AsianGRS_379.tsv")
WangEast_AsianGRS379 <- WangEast_AsianGRS379 %>%
  dplyr::select(c("eid", "Conti_GRS_EAS_avg")) %>%
  dplyr::rename('WangEast_AsianGRS379' = 'Conti_GRS_EAS_avg') %>%
  dplyr::mutate(top10_all_WangEast_AsianGRS379 = WangEast_AsianGRS379 >= quantile(WangEast_AsianGRS379, probs = 0.9))

WangHispanicGRS446 <- read_delim("Wang_HispanicGRS_446.tsv")
WangHispanicGRS446 <- WangHispanicGRS446 %>%
  dplyr::select(c("eid", "Conti_GRS_HIS_avg")) %>%
  dplyr::rename('WangHispanicGRS446' = 'Conti_GRS_HIS_avg') %>%
  dplyr::mutate(top10_all_WangHispanicGRS446 = WangHispanicGRS446 >= quantile(WangHispanicGRS446, probs = 0.9))


All_Wang_GRS <- merge(WangmultiethnicGRS, WangEuropeanGRS, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangAfricanGRS, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangEast_AsianGRS, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangHispanicGRS, by = "eid")

All_Wang_GRS <- merge(All_Wang_GRS, WangmultiethnicGRS450, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangEuropeanGRS445, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangAfricanGRS444, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangEast_AsianGRS379, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangHispanicGRS446, by = "eid")

## Add Schumacher and BARCODE1 GRSs

SchumacherGRS <- read_delim("Schumacher.pgs.tsv")
SchumacherGRS <- SchumacherGRS %>%
  dplyr::rename('SchumacherGRS' = 'Schumacher') %>%
  dplyr::mutate(top10_all_SchumacherGRS = SchumacherGRS >= quantile(SchumacherGRS, probs = 0.9))

BARCODE1GRS <- read_delim("BARCODE1.pgs.tsv")
BARCODE1GRS <- BARCODE1GRS %>%
  dplyr::rename('BARCODE1GRS' = 'BARCODE1') %>%
  dplyr::mutate(top10_all_BARCODE1GRS = BARCODE1GRS >= quantile(BARCODE1GRS, probs = 0.9))

SchumacherGRS145 <- read.delim("SchumacherGRS_145.tsv")
SchumacherGRS145 <- SchumacherGRS145 %>%
  dplyr::select(c("eid", "Schumacher_GRS_avg")) %>%
  dplyr::rename('SchumacherGRS145' = 'Schumacher_GRS_avg') %>%
  dplyr::mutate(top10_all_SchumacherGRS145 = SchumacherGRS145 >= quantile(SchumacherGRS145, probs = 0.9))

BARCODE1GRS129 <- read.delim("BARCODE1GRS_129.tsv")
BARCODE1GRS129 <- BARCODE1GRS129 %>%
  dplyr::select(c("eid", "BARCODE1_GRS_avg")) %>%
  dplyr::rename('BARCODE1GRS129' = 'BARCODE1_GRS_avg') %>%
  dplyr::mutate(top10_all_BARCODE1GRS129 = BARCODE1GRS129 >= quantile(BARCODE1GRS129, probs = 0.9))

## Add Seibert and Pagadala GRSs

SeibertGRS <- read_delim("Seibert.pgs.tsv")
SeibertGRS <- SeibertGRS %>%
  dplyr::rename('SeibertGRS' = 'Seibert') %>%
  dplyr::mutate(top10_all_SeibertGRS = SeibertGRS >= quantile(SeibertGRS, probs = 0.9))

PagadalaGRS <- read_delim("Pagadala.pgs.tsv")
PagadalaGRS <- PagadalaGRS %>%
  dplyr::rename('PagadalaGRS' = 'Pagadala') %>%
  dplyr::mutate(top10_all_PagadalaGRS = PagadalaGRS >= quantile(PagadalaGRS, probs = 0.9))

SeibertGRS52 <- read_delim("SeibertGRS_52.tsv")
SeibertGRS52 <- SeibertGRS52 %>%
  dplyr::select(c("eid", "Seibert_GRS_avg")) %>%
  dplyr::rename('SeibertGRS52' = 'Seibert_GRS_avg') %>%
  dplyr::mutate(top10_all_SeibertGRS52 = SeibertGRS52 >= quantile(SeibertGRS52, probs = 0.9))

PagadalaGRS285 <- read_delim("PagadalaGRS_285.tsv")
PagadalaGRS285 <- PagadalaGRS285 %>%
  dplyr::select(c("eid", "Pagadala_GRS_avg")) %>%
  dplyr::rename('PagadalaGRS285' = 'Pagadala_GRS_avg') %>%
  dplyr::mutate(top10_all_PagadalaGRS285 = PagadalaGRS285 >= quantile(PagadalaGRS285, probs = 0.9))

## Merge all GRSs into one

All_GRS <- merge(All_Conti_GRS, All_Wang_GRS, by = "eid", all = T)
All_GRS <- merge(All_GRS, SchumacherGRS, by = "eid", all = T)
All_GRS <- merge(All_GRS, BARCODE1GRS, by = "eid", all = T)
All_GRS <- merge(All_GRS, SchumacherGRS145, by = "eid", all = T)
All_GRS <- merge(All_GRS, BARCODE1GRS129, by = "eid", all = T)
All_GRS <- merge(All_GRS, SeibertGRS, by = "eid", all = T)
All_GRS <- merge(All_GRS, PagadalaGRS, by = "eid", all = T)
All_GRS <- merge(All_GRS, SeibertGRS52, by = "eid", all = T)
All_GRS <- merge(All_GRS, PagadalaGRS285, by = "eid", all = T)

#########################################################################################
# Step 5 - Merge all variables/covariates into one dataframe with Prostate Cancer cases #
#########################################################################################

iv_covariates <- merge(iv, covariates, by = "eid", all.y = T)

PCa_iv_covariates <- merge(iv_covariates, PrCa_symptoms_diagnosis, by = "eid", all.y = TRUE)

PCa_iv_covariates_GRS <- merge(PCa_iv_covariates, All_GRS, by = "eid", all.x = TRUE)

PCa_iv_covariates_GRS_severity <- merge(PCa_iv_covariates_GRS, actionable_criteria, by = 'eid', all = T)

## Extra exclusion criteria for controls 
## (derived from https://phekb.org/phenotype/prostate-cancer-0 see tables: https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fphekb.org%2Fsites%2Fphenotype%2Ffiles%2FPrCa%2520Phenotyping%2520Algorithm%2520codes.xlsx&wdOrigin=BROWSELINK
## with alterations to remove ICD9: V84.03 and ICD10: Z15.03 since these describe genetic susceptibility to prostate cancer)
## AND with some creative interpretation for OPCS codes, since they are given as CPT codes in the spreadsheet

exclusions_ICD9 <- read_ICD9(c(185,         # 185: Malignant neoplasm of prostate
                              'V104',       # V10.4: Personal history of malignant neoplasm of genital organs
                               2334,        # 233.4: Carcinoma in situ of the prostate
                               2365,        # 236.5: Neoplasm of uncertain behavior of the prostate
                               6023,        # 602.3: Dysplasia of the prostate
                               6021,        # 60.21: Transurethral (ultrasound) guided laser induced prostatectomy (TULIP)
                               6029,        # 60.29: Other transurethral prostatectomy
                               603,         # 60.3: Suprapubic prostatectomy
                               604,         # 60.4: Retropubic prostatectomy
                               605,         # 60.5: Radical prostatectomy
                               6061,        # 60.61: Local excision of lesion of prostate
                               6062,        # 60.62: Perineal prostatectomy
                               6069)        # 60.69: Other prostatectomy
                               ) %>%
  dplyr::select("eid", "diag_icd9")

exclusions_ICD10 <- read_ICD10(c('C61',     # C61: Malignant neoplasm of prostate
                                'Z854',     # Z85.4: Personal History of malignant neoplasm of genital organs
                                'R972',     # R97.2: Elevated prostate specific antigen [PSA]
                                'D075',     # D07.5: Carcinoma in situ of prostate
                                'D400',     # D40.0: Neoplasm of uncertain behavior of prostate
                                'N423')     # N42.3: Dysplasia of prostate
                                )%>%
  dplyr::select("eid", "diag_icd10")

exclusions_OPCS <- read_OPCS(c('M61',       # M61: Prostatectomy
                               'M611',      # M61.1: Radical prostatectomy
                               'M612',      # M61.2: Retropubic Prostatectomy
                               'M613',      # M61.3: Transvesical Prostatectomy
                               'M614',     # M61.4: Perineal Prostatectomy
                               'X65',       # X65: Radiotherapy Delivery
                               'X67',       # X67: Preparation of radiotherapy
                               'X68',       # X68: Brachytherapy preparation
                               'M706',      # M70.6 Radioactive seed implantation into prostate
                               'Y35',       # Y35: Introduction Material Radioactive Removable NOC
                               'Y36',       # Y36: Introduction Material Non-removable NOC
                               'T856',      # T85.6: Block dissection of pelvic lymph nodes
                               'M702',      # M70.2: Perineal needle biopsy of prostate
                               'M703',      # M70.3: Rectal needle biopsy of prostate
                               'N04',       # N04: Orchidectomy
                               'M65',       # M65: Endoscopic resection of prostate
                               'M68',       # M68: Endoscopic insertion of prosthesis into prostate
                               'M671',      # M67.1: Endoscopic cryotherapy to lesion of prostate
                               'M711',      # M71.1: High intensity focused ultrasound of prostate
                               'M712')      # M71.2: Implantation of radioactive substance into prostate
                               ) %>%
  dplyr::select("eid", "oper4")

exclusions_cancerregistry <- read_cancer(c('C61',    # C61: Malignant neoplasm of prostate 
                                           'Z854',    # Z85.4: Personal History of malignant neoplasm of genital organs
                                           'R972',    # R97.2: Elevated prostate specific antigen [PSA] (note: no results returned)
                                           'D075',    # D07.5: Carcinoma in situ of prostate
                                           'D400',    # D40.0: Neoplasm of uncertain behavior of prostate
                                           'N423')    # N42.3: Dysplasia of prostate (note: no results returned)
)%>
  dplyr::select("eid", "ICD10")

possible_PrCa_cases <- merge(exclusions_ICD9, exclusions_ICD10, by = "eid", all = T)
possible_PrCa_cases <- merge(possible_PrCa_cases, exclusions_OPCS, by = "eid", all = T)
possible_PrCa_cases <- merge(possible_PrCa_cases, exclusions_cancerregistry, by = "eid", all = T)

# Collapse exclusions to one row per participant, leaving only eids
possible_PrCa_cases <- possible_PrCa_cases %>%
  group_by(eid) %>%
  dplyr::summarise(possible_PrCa_case = 1)

# Join onto main dataframe and create "exclude" variable to exclude controls who meet any of the exclusion criteria for controls (i.e. those with possible PrCa diagnoses but no evidence of prostate cancer diagnosis)

PCa_iv_covariates_GRS_severity <- merge(PCa_iv_covariates_GRS_severity, possible_PrCa_cases, by = "eid", all.x = TRUE) %>%
  dplyr::mutate(exclude = if_else(possible_PrCa_case == 1 & is.na(pre_diagnosed), 1L, 0L, missing = 0L))

## Remove anybody who is female, lacks GRS data, anyone pre-diagnosed, anyone with HES-only PrCa diagnosis

PCa_iv_covariates_GRS_clean <- PCa_iv_covariates_GRS_severity %>%
  dplyr::filter(
    Sex == 'Male',
    !is.na(ContimultiethnicGRS),
    pre_diagnosed == FALSE | is.na(pre_diagnosed), ## to remove pre-diagnosed == TRUE patients
    HES_only == FALSE | is.na(HES_only), ## to remove patients with HES-only PrCa diagnoses
    exclude == 0L | is.na(exclude) ## to remove controls who meet any of the exclusion criteria for controls
  )

# Sanity Check - how many patients in PCa_iv_covariates were female or lacked GRS data? Does it match the difference in n between PCa_iv_covariates and PCa_iv_covariates_clean ?

PCa_iv_covariates_GRS_severity %>%
  filter(
    Sex == "Female" |
      is.na(ContimultiethnicGRS) |
      pre_diagnosed == TRUE |
      HES_only == TRUE |
      exclude == 1L
  ) %>%
  summarise(n = n())

# Set Prostate Cancer cases to 1, controls to 0. Also rename event_age to simply Age

PCa_iv_covariates_GRS_clean <- PCa_iv_covariates_GRS_clean %>%
  dplyr::mutate(
    PrCa = if_else(!is.na(pre_diagnosed) & HES_only == FALSE, 1L, 0L, missing = 0L),
  ) %>%
  dplyr::rename(Age = event_age)

############################################# (these should have already been filtered out,
# Important - remove withdrawn participants #  since the withdrawn participants no longer have 
#############################################  a GRS calculated, but still worth making sure)

exclude_withdrawn=function(df){
  system('dx download Callum/Withdrawals/withdrawn_20260310.csv --overwrite')  ## This file is a list of participants who withdrew from the Biobank up to the date 10th March 2026. This was sent from the UK Biobank team via email to members of approved applications
  df2 = df %>% left_join(
    read_csv("withdrawn_20260310.csv", col_names = FALSE, show_col_types = FALSE) %>%
      mutate(w=1) %>%
      rename(eid=X1),
    by='eid'
  ) %>%
    filter(is.na(w))
  return(df2)
}

PCa_iv_covariates_GRS_clean <- exclude_withdrawn(PCa_iv_covariates_GRS_clean)

##################################################################################
# Step 6 - Set prediction horizons, including for general/actionable/severe PrCa #
##################################################################################

library(lubridate)

PCa_iv_covariates_GRS_predhorizon <- PCa_iv_covariates_GRS_clean %>%
  dplyr::mutate(
    PrCa_post_symptoms  = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom
    ),
    PrCa_2yrs  = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom &
        earliest_PrCa_date <= EarliestDate_symptom %m+% lubridate::years(2)
    ),
    PrCa_5yrs  = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom &
        earliest_PrCa_date <= EarliestDate_symptom %m+% lubridate::years(5)
    ),
    PrCa_10yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom &
        earliest_PrCa_date <= EarliestDate_symptom %m+% lubridate::years(10)
    ),
    PrCa_severe = as.integer(
      !is.na(earliest_PrCa_date) & (
        (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
          (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2))
      )
    ),
    PrCa_severe_2yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom &
        earliest_PrCa_date <= EarliestDate_symptom %m+% lubridate::years(2) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    ),
    PrCa_severe_5yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom &
        earliest_PrCa_date <= EarliestDate_symptom %m+% lubridate::years(5) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    ),
    PrCa_severe_10yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom &
        earliest_PrCa_date <= EarliestDate_symptom %m+% lubridate::years(10) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    ),
    PrCa_actionable = as.integer(
      !is.na(earliest_PrCa_date) & (
        (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
          (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2)) |
          (!is.na(surgery_opdate) & surgery_opdate >= earliest_PrCa_date & surgery_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
          (!is.na(radio_opdate) & radio_opdate >= earliest_PrCa_date & radio_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
          (!is.na(androgen_opdate) & androgen_opdate >= earliest_PrCa_date & androgen_opdate <= earliest_PrCa_date %m+% lubridate::years(2))
      )
    ),
    PrCa_actionable_2yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom &
        earliest_PrCa_date <= EarliestDate_symptom %m+% lubridate::years(2) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(surgery_opdate) & surgery_opdate >= earliest_PrCa_date & surgery_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(radio_opdate) & radio_opdate >= earliest_PrCa_date & radio_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(androgen_opdate) & androgen_opdate >= earliest_PrCa_date & androgen_opdate <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    ),
    PrCa_actionable_5yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom &
        earliest_PrCa_date <= EarliestDate_symptom %m+% lubridate::years(5) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(surgery_opdate) & surgery_opdate >= earliest_PrCa_date & surgery_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(radio_opdate) & radio_opdate >= earliest_PrCa_date & radio_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(androgen_opdate) & androgen_opdate >= earliest_PrCa_date & androgen_opdate <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    ),
    PrCa_actionable_10yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(EarliestDate_symptom) &
        earliest_PrCa_date >= EarliestDate_symptom &
        earliest_PrCa_date <= EarliestDate_symptom %m+% lubridate::years(10) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(surgery_opdate) & surgery_opdate >= earliest_PrCa_date & surgery_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(radio_opdate) & radio_opdate >= earliest_PrCa_date & radio_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(androgen_opdate) & androgen_opdate >= earliest_PrCa_date & androgen_opdate <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    )
  )

################################################################
# Step 7 - Subset data for different ethnicity/ancestry groups #
################################################################

# Prediction horizons for White patients only

PCa_iv_covariates_GRS_predhorizon_WhiteOnly <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(ethnicity_group_narrow == "White")

#PCa_iv_covariates_GRS_predhorizon_EurOnly_subset <- PCa_iv_covariates_GRS_predhorizon_EurOnly %>%


# Prediction horizons for Black patients only

PCa_iv_covariates_GRS_predhorizon_BlackOnly <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(ethnicity_group_narrow == "Black")

# Prediction horizons for Mixed patients only

PCa_iv_covariates_GRS_predhorizon_Mixed <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(ethnicity_group_narrow == "Mixed White and Black")

# Prediction horizons for Black+Mixed participants

PCa_iv_covariates_GRS_predhorizon_BlackMixed <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(ethnicity_group_wide == "Black")

# Prediction horizons for European participants

PCa_iv_covariates_GRS_predhorizon_EUROnly <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(Genomic_ancestry == "European ancestry (EUR)")

# Prediction horizons for African participants

PCa_iv_covariates_GRS_predhorizon_AFROnly <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(Genomic_ancestry == "African ancestry (AFR)")

# Prediction horizons for East Asian participants

PCa_iv_covariates_GRS_predhorizon_EASOnly <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(Genomic_ancestry == "East Asian ancestry (EAS)")

# Prediction horizons for Central/South Asian participants

PCa_iv_covariates_GRS_predhorizon_CSAOnly <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(Genomic_ancestry == "Central/South Asian ancestry (CSA)")
#dplyr::filter(ethnicity_group_narrow == "South Asian")

# Prediction horizons for Middle Eastern participants

PCa_iv_covariates_GRS_predhorizon_MIDOnly <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(Genomic_ancestry == "Middle Eastern ancestry (MID)")

# Prediction horizons for Admixed American participants

PCa_iv_covariates_GRS_predhorizon_AMROnly <- PCa_iv_covariates_GRS_predhorizon %>%
  dplyr::filter(Genomic_ancestry == "Admixed American ancestry (AMR)")



################################################################################
# Step 8 (optional) - Visually inspect GRS distribution for cases vs. controls #
################################################################################

ggplot(data=PCa_iv_covariates_GRS_predhorizon, aes(x=ContimultiethnicGRS,colour=as.factor(PrCa)))+ # PrCa in general
  geom_density()+
  theme_bw()

ggplot(data=PCa_iv_covariates_GRS_predhorizon, aes(x=ContimultiethnicGRS,colour=as.factor(PrCa_2yrs)))+ # PrCa within 2 years
  geom_density()+
  theme_bw()

ggplot(data=PCa_iv_covariates_GRS_predhorizon, aes(x=ContimultiethnicGRS,colour=as.factor(PrCa_5yrs)))+ # PrCa within 5 years
  geom_density()+
  theme_bw()

ggplot(data=PCa_iv_covariates_GRS_predhorizon, aes(x=ContimultiethnicGRS,colour=as.factor(PrCa_10yrs)))+ # PrCa within 10 years
  geom_density()+
  theme_bw()

###########################################################################
# Step 9 - Model Logistic Regression, Confusion Matrix, and OR/RR tables  #
###########################################################################

# Set "data" to either: 
#   - PCa_iv_covariates_GRS_predhorizon (all participants)
#   - PCa_iv_covariates_GRS_predhorizon_EurOnly (White participants)
#   - PCa_iv_covariates_GRS_predhorizon_BlackOnly (Black participants)
#   - PCa_iv_covariates_GRS_predhorizon_Mixed (Mixed White and Black participants)
#   - PCa_iv_covariates_GRS_predhorizon_BlackMixed (Black + Mixed White and Black participants)
#   - PCa_iv_covariates_GRS_predhorizon_EUROnly (EUR-like participants)
#   - PCa_iv_covariates_GRS_predhorizon_AFROnly (AFR-like participants)
#   - PCa_iv_covariates_GRS_predhorizon_EASOnly (EAS-like participants)
#   - PCa_iv_covariates_GRS_predhorizon_CSAOnly (CSA-like participants)
#   - PCa_iv_covariates_GRS_predhorizon_MIDOnly (MID-like participants)
#   - PCa_iv_covariates_GRS_predhorizon_AMROnly (AMR-like participants)
# 
# 
# Set "outcome" to either: 
#   - PrCa (Prostate Cancer diagnosis after symptom presentation)
#   - PrCa_2yrs (Prostate Cancer diagnosis within 2 years after symptom presentation)
#   - PrCa_5yrs (Prostate Cancer diagnosis within 5 years after symptom presentation)
#   - PrCa_10yrs (Prostate Cancer diagnosis within 10 years after symptom presentation)
#
#   - PrCa_actionable (Prostate Cancer diagnosis after symptom presentation that satisfies "Actionable" criteria within 2 years of diagnosis)
#   - PrCa_actionable_2yrs (Prostate Cancer diagnosis within 2 years after symptom presentation that satisfies "Actionable" criteria within 2 years of diagnosis)
#   - PrCa_actionable_5yrs (Prostate Cancer diagnosis within 5 years after symptom presentation that satisfies "Actionable" criteria within 2 years of diagnosis)
#   - PrCa_actionable_10yrs (Prostate Cancer diagnosis within 10 years after symptom presentation that satisfies "Actionable" criteria within 2 years of diagnosis)
#
#   - PrCa_severe (Prostate Cancer diagnosis after symptom presentation that satisfies "Severe" criteria within 2 years of diagnosis)
#   - PrCa_severe_2yrs (Prostate Cancer diagnosis within 2 years after symptom presentation that satisfies "Severe" criteria within 2 years of diagnosis)
#   - PrCa_severe_5yrs (Prostate Cancer diagnosis within 5 years after symptom presentation that satisfies "Severe" criteria within 2 years of diagnosis)
#   - PrCa_severe_10yrs (Prostate Cancer diagnosis within 10 years after symptom presentation that satisfies "Severe" criteria within 2 years of diagnosis)
#
# Set "predictor" to either: 
#
#    Note: these top 6 GRSs use "Conti_script.R", which drops 4 SNPs by default. The numbered GRSs below only drop 2 SNPs
#
#   - ContimultiethnicGRS (Conti's GRS with pan-Ancestry weights)
#   - ContiEuropeanGRS (Conti's GRS with European-specific weights)
#   - ContiAfricanGRS (Conti's GRS with African-specific weights)
#   - ContiEast_AsianGRS (Conti's GRS with East Asian-specific weights)
#   - ContiHispanicGRS (Conti's GRS with Hispanic-specific weights)
#   - ContiadjustedGRS (Conti's GRS adjusted for ancestry probability at the Beta level)
#
#   - ContimultiethnicGRS267 (Conti's GRS with pan-Ancestry weights, all 267 available SNPs)
#   - ContiEuropeanGRS265 (Conti's GRS with European-specific weights, all 265 available SNPs)
#   - ContiAfricanGRS246 (Conti's GRS with African-specific weights, all 246 available SNPs)
#   - ContiEast_AsianGRS222 (Conti's GRS with East Asian-specific weights)
#   - ContiHispanicGRS253 (Conti's GRS with Hispanic-specific weights)
#
#   - ContiORadjustedGRS (Conti's GRS adjusted for ancestry probability at the Odds Ratio level)
#
#   - WangmultiethnicGRS (Wang's GRS with pan-Ancestry weights)
#   - WangEuropeanGRS (Wang's GRS with European-specific weights)
#   - WangAfricanGRS (Wang's GRS with African-specific weights)
#   - WangEast_AsianGRS (Wang's GRS with East Asian-specific weights)
#   - WangHispanicGRS (Wang's GRS with Hispanic-specific weights)
#
#   - WangmultiethnicGRS450 (Wang's GRS with pan-Ancestry weights, all 450 available SNPs)
#   - WangEuropeanGRS445 (Wang's GRS with European-specific weights, all 445 available SNPs)
#   - WangAfricanGRS444 (Wang's GRS with African-specific weights, all 444 available SNPs)
#   - WangEast_AsianGRS379 (Wang's GRS with East Asian-specific weights, all 379 available SNPs)
#   - WangHispanicGRS446 (Wang's GRS with Hispanic-specific weights, all 446 available SNPs)
#
#   - SchumacherGRS (Schumacher's GRS (2018)))
#   - BARCODE1GRS (BARCODE1 GRS (2021)))
#   - SchumacherGRS145 (Schumacher's GRS with 145 available SNPs)
#   - BARCODE1GRS129 (BARCODE1 GRS with 129 available SNPs)
#
#    - SeibertGRS (Seibert's GRS (2018))
#    - PagadalaGRS (Pagadala's GRS (2022))
#    - SeibertGRS52 (Seibert's GRS with 52 available SNPs)
#    - PagadalaGRS285 (Pagadala's GRS with 285 available SNPs)
#
# Set "covariates" to either: (or add multiple using + between covariates)
#
#   - (without quote marks) NULL
#   - Age (Age at symptom presentation)
#   - rs72725854_T (carrier status of rs72725854 risk allele)
#


# Logistic Regression

model <- run_logreg(data = PCa_iv_covariates_GRS_predhorizon_BlackOnly,
                    outcome = "PrCa_10yrs",
                    predictor = "Age",
                    covariates = NULL,
                    plot_roc = TRUE)

model2 <- run_logreg(data = PCa_iv_covariates_GRS_predhorizon_BlackOnly, ## model2 is used for NRI comparison with model1
                     outcome = "PrCa_10yrs",
                     predictor = "ContimultiethnicGRS",
                     covariates = "Age",
                     plot_roc = TRUE)

# Confusion Matrix

matrix <- confusion_matrix(data = model$data, 
                           outcome = model$outcome, 
                           cutoff_value = 0.20, 
                           positive_level = 1, negative_level = 0)

# Odds Ratio Table

OR_table <- ORtable(
  data = model$data,
  outcome = model$outcome,
  group_col = "ethnicity_group_narrow",       # <- your 6-level grouping variable
  positive_level = 1,        # 1 denotes positive outcome
  use_existing_cols = FALSE  # If in doubt, leave as FALSE. Set to TRUE if you have already created "predtopXX" columns for the desired bins and want to reuse them (must be global bins, not within-group bins)
)

print(OR_table$wide_formatted)

# Risk Ratio Table

RR_table <- RRtable(
  data = model$data,
  outcome = model$outcome,
  group_col = "ethnicity_group_narrow",
  positive_level = 1,
  use_existing_cols = TRUE
)

print(RR_table$wide_formatted)


# NRI Calculation to compare two models

nri_result <- nri(
  data = model$data %>% dplyr::mutate(pred2 = model2$data$pred),
  outcome = model$outcome,
)

print(nri_result)

######################################################################################
# Step 10 - Generate comprehensive logreg summary table for all GRSs and populations # (the lazy way)
######################################################################################

## By default, logreg_table() will compute all combinations of population, outcome, GRS, and covariates

bulk_LUTS <- logreg_table()

## The below block adds, for each row that represents a GRS + Age model, the
## equivalent Age-only model, and compares confidence intervals between them

age_only_reference <- bulk_LUTS %>%
  dplyr::filter(Predictor == "Age", Covariates == "None") %>%
  dplyr::select(
    Outcome,
    Population,
    Age_only_ROC_AUC_CI_Upper = ROC_AUC_CI_Upper
  )

formatted_LUTS <- bulk_LUTS %>%   ## To present ROC AUC and 95% CIs to 4 decimal places
  dplyr::left_join(age_only_reference, by = c("Outcome", "Population")) %>%
  dplyr::mutate(
    ROC_AUC_4dp = dplyr::if_else(
      is.na(ROC_AUC),
      NA_character_,
      sprintf("%.4f", ROC_AUC)
    ),
    ROC_AUC_CI_95_4dp = dplyr::if_else(
      is.na(ROC_AUC_CI_Lower) | is.na(ROC_AUC_CI_Upper),
      NA_character_,
      sprintf("%.4f [%.4f-%.4f]", ROC_AUC, ROC_AUC_CI_Lower, ROC_AUC_CI_Upper)
    ),
    `GRS+Age > Age?` = dplyr::case_when(
      Covariates != "Age" | Predictor == "Age" ~ NA_character_,
      is.na(ROC_AUC_CI_Lower) | is.na(Age_only_ROC_AUC_CI_Upper) ~ NA_character_,
      ROC_AUC_CI_Lower > Age_only_ROC_AUC_CI_Upper ~ "YES",
      ROC_AUC_CI_Lower <= Age_only_ROC_AUC_CI_Upper ~ "NO"
    )
  ) %>%
  dplyr::select(c("Outcome", "Population", "Predictor", "Covariates", "N_Cases", "N_Controls", "ROC_AUC_CI_95_4dp", "Age_only_ROC_AUC_CI_Upper", "GRS+Age > Age?"))

## This block is to view a subset of the bulk logistic regression table. Change
## the filter to investigate a specific Population, Predictor, Outcome, or Covariate

subset_LUTS <- formatted_LUTS %>%                        
  dplyr::filter(                          
    Population == "Black",
    #Predictor == "PagadalaGRS",
    #Outcome == "PrCa_severe_10yrs",
    #Covariates == "Age" #| Covariates == "None"
  )    


######################################################################################
# Step 11 - Generate comprehensive NRI table for all GRSs, populations, and outcomes # 
######################################################################################

## By default, nri_table() will compute all combinations of population, outcome, and GRS for the NRI comparison between a GRS+Age model vs. an Age-only model

nri_bulk_LUTS <- nri_table()

nri_bulk_LUTS <- nri_bulk_LUTS %>%
  dplyr::mutate(
    NRI_significant = dplyr::case_when(
      is.na(p_NRI) ~ NA_character_,
      p_NRI < 0.05 ~ "YES",
      p_NRI >= 0.05 ~ "NO"
    )
  )


## This block is to view a subset of the bulk NRI table. Change the filter to investigate a specific Population, Predictor, or Outcome

nri_subset_LUTS <- nri_bulk_LUTS %>%
  dplyr::filter(
    Population == "Black",
    #GRS == "PagadalaGRS",
    #Outcome == "PrCa_severe_10yrs"
  )