#############################################################################
# Symptomatic Angle - how well can GRS predict Prostate cancer diagnosis in #
# ------------------- symptomatic participants? ----------------------------#
#############################################################################

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
dxdownload("Callum/Derived_datasets/AncestryProbability.csv")	# Genetic similarity probabilities for Ancestry group calculated using AncestryProbability1_bigsnpr folder scripts
dxdownload("Callum/Derived_datasets/AncestryProbability2.csv") 	# Genetic similarity probabilities for Ancestry group calculated using AncestryProbability2_plink folder script

dxdownload("Callum/Derived_datasets/FH_PrCa_BrCa.csv")		# Dataset of Family History of Prostate Cancer and Breast Cancer, created using "Family_history.R" script
dxdownload("Callum/GRSs/Conti_multi_ethnic.pgs.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "Multiethnic Analysis"
dxdownload("Callum/GRSs/Conti_European.pgs.tsv") 		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "European...16"
dxdownload("Callum/GRSs/Conti_African.pgs.tsv")			# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "African...19"
dxdownload("Callum/GRSs/Conti_East_Asian.pgs.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "East Asian...22"
dxdownload("Callum/GRSs/Conti_Hispanic.pgs.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_script" with "effect_weight" set to "Hispanic...25"

dxdownload("Callum/GRSs/Wang_multi_ethnic.pgs.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/Wang_European.pgs.tsv") 		# Calculated Wang GRS for all eligible participants using "Conti_script" 
dxdownload("Callum/GRSs/Wang_African.pgs.tsv")			# Calculated Wang GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/Wang_East_Asian.pgs.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/Wang_Hispanic.pgs.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_script"

dxdownload("Callum/GRSs/OR_adjustedGRS_267.tsv")	# Calculated using "Create_adjustedGRS_weighting_Conti_odds_ratios_by_AncestryProbability2.R" script


system("dx download Callum/LUTS/Green2022supplementarytable1.csv") # This is the supplementary table 1 from Harry's LUTS paper, converted to .csv format using excel. Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC9553867/
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




# Dataframe for radical prostatectomy

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
  select(c("eid", "rs72725854_G"))

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
  dplyr::rename("eid" = "IID")

family_history <- read.csv("FH_PrCa_BrCa.csv") %>%
  dplyr::select(c("eid", "FH_Prostate_cancer", "FH_Breast_cancer", "FH_PrCa_BrCa"))

covariates <- merge(covariates, ethnicity, by = "eid", all.x = T)
covariates <- merge(covariates, principal_components, by = "eid", all = T)
covariates <- merge(covariates, family_history, by = "eid", all.x = T)



multiethnicGRS <- read_delim("Conti_multi_ethnic.pgs.tsv")
multiethnicGRS <- multiethnicGRS %>%
  dplyr::rename('multiethnicGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_multiethnicGRS = multiethnicGRS >= quantile(multiethnicGRS, probs = 0.9))

EuropeanGRS <- read_delim("Conti_European.pgs.tsv")
EuropeanGRS <- EuropeanGRS %>%
  dplyr::rename('EuropeanGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_EuropeanGRS = EuropeanGRS >= quantile(EuropeanGRS, probs = 0.9))

AfricanGRS <- read_delim("Conti_African.pgs.tsv")
AfricanGRS <- AfricanGRS %>%
  dplyr::rename('AfricanGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_AfricanGRS = AfricanGRS >= quantile(AfricanGRS, probs = 0.9))

East_AsianGRS <- read_delim("Conti_East_Asian.pgs.tsv")
East_AsianGRS <- East_AsianGRS %>%
  dplyr::rename('East_AsianGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_East_AsianGRS = East_AsianGRS >= quantile(East_AsianGRS, probs = 0.9))

HispanicGRS <- read_delim("Conti_Hispanic.pgs.tsv")
HispanicGRS <- HispanicGRS %>%
  dplyr::rename('HispanicGRS' = 'Conti') %>%
  dplyr::mutate(top10_all_HispanicGRS = HispanicGRS >= quantile(HispanicGRS, probs = 0.9))

# Merge all data

All_Conti_GRS <- merge(multiethnicGRS, EuropeanGRS, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, AfricanGRS, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, East_AsianGRS, by = "eid")
All_Conti_GRS <- merge(All_Conti_GRS, HispanicGRS, by = "eid")

### Creating an adjusted GRS ###

AncestryProbability <- read.csv("AncestryProbability.csv") %>%          
  dplyr::select(c("IID", "AMR", "AFR", "CSA", "EAS", "EUR", "MID")) %>%
  dplyr::rename("eid" = "IID")

AncestryProbability2 <- read.csv("AncestryProbability2.csv") %>%         
  dplyr::select(c("IID", "AMR", "AFR", "CSA", "EAS", "EUR", "MID")) %>%
  dplyr::rename("eid" = "IID")

GRS_plus_ancestry <- merge(All_Conti_GRS, AncestryProbability2, all=T, by="eid") # Change y to either AncestryProbability or AncestryProbability 2 depending on preference
GRS_plus_ancestry <- merge(GRS_plus_ancestry, covariates, by = "eid", all.x=T) %>%
  dplyr::select(c("eid", "multiethnicGRS", "top10_all_multiethnicGRS", "EuropeanGRS", "top10_all_EuropeanGRS", "AfricanGRS", "top10_all_AfricanGRS", "East_AsianGRS", "top10_all_East_AsianGRS", "HispanicGRS", "top10_all_HispanicGRS", "AMR", "AFR", "CSA", "EAS", "EUR", "MID", "Genomic_ancestry"))

GRS_plus_ancestry <- GRS_plus_ancestry %>%
  dplyr::mutate(
    adjustedGRS = (EUR*EuropeanGRS)+(AFR*AfricanGRS)+(EAS*East_AsianGRS)+(AMR*HispanicGRS)+(CSA*multiethnicGRS)+(MID*AfricanGRS)
  )

All_Conti_GRS <- GRS_plus_ancestry %>%
  dplyr::select(c("eid", "EUR", "AFR", "EAS", "CSA", "MID", "AMR", "multiethnicGRS", "top10_all_multiethnicGRS", "EuropeanGRS", "top10_all_EuropeanGRS", "AfricanGRS", "top10_all_AfricanGRS", "East_AsianGRS", "top10_all_East_AsianGRS", "HispanicGRS", "top10_all_HispanicGRS", "adjustedGRS"))

## Add Odds Ratio-adjusted GRS (where ancestry probability is applied to the OR, not the final scores)

ORadjustedGRS <- read.table("OR_adjustedGRS_267.tsv", header = T) %>%
  dplyr::select(c("eid", "Conti_ORfirst_avg")) %>%
  dplyr::rename("ORadjustedGRS" = "Conti_ORfirst_avg")

All_Conti_GRS <- merge(All_Conti_GRS, ORadjustedGRS, by = "eid", all.x = T)

### Now repeat for Wang GRSs

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

All_Wang_GRS <- merge(WangmultiethnicGRS, WangEuropeanGRS, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangAfricanGRS, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangEast_AsianGRS, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangHispanicGRS, by = "eid")

All_GRS <- merge(All_Conti_GRS, All_Wang_GRS, by = "eid", all = T)



#########################################################################################
# Step 5 - Merge all variables/covariates into one dataframe with Prostate Cancer cases #
#########################################################################################

iv_covariates <- merge(iv, covariates, by = "eid", all.y = T)

PCa_iv_covariates <- merge(iv_covariates, PrCa_symptoms_diagnosis, by = "eid", all.y = TRUE)

PCa_iv_covariates_GRS <- merge(PCa_iv_covariates, All_GRS, by = "eid", all.x = TRUE)

PCa_iv_covariates_GRS_severity <- merge(PCa_iv_covariates_GRS, actionable_criteria, by = 'eid', all.x = T)

# Remove anybody who is female or who lacks GRS data (optionally, also remove anyone pre-diagnosed)

PCa_iv_covariates_GRS_clean <- PCa_iv_covariates_GRS_severity %>%
  dplyr::filter(
    Sex == 'Male',
    !is.na(multiethnicGRS),
    pre_diagnosed == FALSE | is.na(pre_diagnosed) ## to remove pre-diagnosed patients
  )

# Sanity Check - how many patients in PCa_iv_covariates were female or lacked GRS data? Does it match the difference in n between PCa_iv_covariates and PCa_iv_covariates_clean ?

PCa_iv_covariates_GRS_severity %>%
  filter(
    Sex == "Female" |
      is.na(multiethnicGRS) |
      pre_diagnosed == TRUE
  ) %>%
  summarise(n = n())

# Set Prostate Cancer cases to 1, controls to 0. 

PCa_iv_covariates_GRS_clean <- PCa_iv_covariates_GRS_clean %>%
  dplyr::mutate(
    PrCa = if_else(!is.na(epistart) | !is.na(date) | !is.na(date_of_death), 1L, 0L, missing = 0L),
  )

#############################################
# Important - remove withdrawn participants #
#############################################

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

ggplot(data=PCa_iv_covariates_GRS_predhorizon, aes(x=multiethnicGRS,colour=as.factor(PrCa)))+ # PrCa in general
  geom_density()+
  theme_bw()

ggplot(data=PCa_iv_covariates_GRS_predhorizon, aes(x=multiethnicGRS,colour=as.factor(PrCa_2yrs)))+ # PrCa within 2 years
  geom_density()+
  theme_bw()

ggplot(data=PCa_iv_covariates_GRS_predhorizon, aes(x=multiethnicGRS,colour=as.factor(PrCa_5yrs)))+ # PrCa within 5 years
  geom_density()+
  theme_bw()

ggplot(data=PCa_iv_covariates_GRS_predhorizon, aes(x=multiethnicGRS,colour=as.factor(PrCa_10yrs)))+ # PrCa within 10 years
  geom_density()+
  theme_bw()

#########################################################################################
# Step 9 - Set up logistic regression, confusion matrix, and Odds Ratio table functions #
#########################################################################################

run_logreg <- function(data,
                       outcome,
                       predictor,
                       covariates = NULL,
                       plot_roc = TRUE) {
  
  if (is.null(covariates)) {
    formula <- as.formula(paste(outcome, "~", predictor))
  } else {
    formula <- as.formula(
      paste(outcome, "~", predictor, "+", paste(covariates, collapse = " + "))
    )
  }
  
  logreg <- glm(formula, data = data, family = binomial)
  print(summary(logreg))
  data$pred <- predict(logreg, data, type = "response")
  data$predtop10 <- data$pred >= quantile(data$pred, probs = 0.9, na.rm = TRUE, names = FALSE)
  data$predtop20 <- data$pred >= quantile(data$pred, probs = 0.8, na.rm = TRUE, names = FALSE)
  data$predtop30 <- data$pred >= quantile(data$pred, probs = 0.7, na.rm = TRUE, names = FALSE)
  data$predtop40 <- data$pred >= quantile(data$pred, probs = 0.6, na.rm = TRUE, names = FALSE)
  data$predtop50 <- data$pred >= quantile(data$pred, probs = 0.5, na.rm = TRUE, names = FALSE)
  data$predtop60 <- data$pred >= quantile(data$pred, probs = 0.4, na.rm = TRUE, names = FALSE)
  data$predtop70 <- data$pred >= quantile(data$pred, probs = 0.3, na.rm = TRUE, names = FALSE)
  data$predtop80 <- data$pred >= quantile(data$pred, probs = 0.2, na.rm = TRUE, names = FALSE)
  data$predtop90 <- data$pred >= quantile(data$pred, probs = 0.1, na.rm = TRUE, names = FALSE)
  
  if (plot_roc) {
    roc_obj <- roc(data[[outcome]] ~ data$pred,
                   plot = TRUE,
                   print.auc = TRUE,
                   ci = TRUE)
  } else {
    roc_obj <- roc(data[[outcome]] ~ data$pred)
  }
  
  print(roc_obj)
  
  return(list(
    model = logreg,
    data = data,
    roc = roc_obj
  ))
}

confusion_matrix <- function(data,
                             outcome,
                             cutoff_value = 0.5,
                             positive_level = 1,
                             negative_level = 0,
                             print_epi = TRUE) {
  
  data$predbin <- factor(
    data$pred > cutoff_value,
    levels = c(TRUE, FALSE)
  )
  
  data$outcome_actual <- factor(
    data[[outcome]],
    levels = c(positive_level, negative_level)
  )
  
  # Build confusion matrix
  cm <- table(
    pred_outcome = data$predbin,
    outcome_actual = data$outcome_actual
  )
  
  # Print confusion matrix
  print(cm)
  
  # Run epi.tests
  if (print_epi) {
    epi <- epi.tests(cm)
    print(epi)
  } else {
    epi <- epi.tests(cm)
  }
  
  # Return useful objects
  return(list(
    confusion_matrix = cm,
    epi_tests = epi,
    data = data
  ))
}


ORtable <- function(data,
                    outcome,
                    group_col = "Group",
                    positive_level = 1,
                    pred_var = "pred",
                    bins = c(10, 20, 30, 40, 50, 60, 70, 80, 90),
                    use_existing_cols = TRUE,
                    digits = 2) {
  # Basic checks
  stopifnot(outcome %in% names(data))
  stopifnot(group_col %in% names(data))
  if (!pred_var %in% names(data) && !use_existing_cols) {
    stop("`pred_var` not found and `use_existing_cols = FALSE`.")
  }
  
  # Make Group a factor (preserve order if already set)
  data[[group_col]] <- as.factor(data[[group_col]])
  
  # Helper: retrieve (or compute) bin membership vector
  get_bin_vector <- function(bin_pct) {
    colname <- paste0("predtop", bin_pct)
    if (use_existing_cols && colname %in% names(data)) {
      v <- data[[colname]]
      # Coerce to logical if it came as numeric/integer
      if (!is.logical(v)) v <- as.logical(v)
      return(v)
    } else {
      # Compute threshold from pred quantile if not present
      if (!pred_var %in% names(data)) {
        stop("Missing `pred` column to compute quantile thresholds; either supply it or set `use_existing_cols = TRUE` with existing predtopXX columns.")
      }
      thr <- stats::quantile(data[[pred_var]], probs = 1 - bin_pct/100, na.rm = TRUE, names = FALSE)
      return(data[[pred_var]] >= thr)
    }
  }
  
  # Containers
  groups <- levels(data[[group_col]])
  bin_labels <- paste0("Top ", bins, "%")
  
  # Long results accumulator
  res_long <- list()
  
  # Iterate over bins and groups
  for (i in seq_along(bins)) {
    bin_pct <- bins[i]
    bin_label <- bin_labels[i]
    inbin_vec <- get_bin_vector(bin_pct)
    
    for (g in groups) {
      # Subset to group g and non-missing outcomes and bin flags
      keep <- !is.na(data[[group_col]]) & !is.na(data[[outcome]]) & !is.na(inbin_vec)
      keep <- keep & (data[[group_col]] == g)
      
      if (!any(keep)) {
        # No data for this group
        res_long[[length(res_long) + 1]] <- data.frame(
          bin = bin_label, group = g,
          n_in = NA_integer_, events_in = NA_integer_, rate_in = NA_real_,
          n_out = NA_integer_, events_out = NA_integer_, rate_out = NA_real_,
          OR_pos = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }
      
      y <- data[[outcome]][keep]
      inbin <- inbin_vec[keep]
      y_pos <- y == positive_level
      
      # 2x2 counts inside the group:
      #          Outcome+
      # inbin=1     a       (in & pos)
      # inbin=1     b       (in & neg)
      # inbin=0     c       (out & pos)
      # inbin=0     d       (out & neg)
      a <- sum(inbin & y_pos, na.rm = TRUE)
      b <- sum(inbin & !y_pos, na.rm = TRUE)
      c <- sum(!inbin & y_pos, na.rm = TRUE)
      d <- sum(!inbin & !y_pos, na.rm = TRUE)
      
      n_in  <- a + b
      n_out <- c + d
      
      # If no one falls in or out of the bin for this group, OR is undefined
      if (n_in == 0L || n_out == 0L) {
        res_long[[length(res_long) + 1]] <- data.frame(
          bin = bin_label, group = g,
          n_in = n_in, events_in = a, rate_in = ifelse(n_in > 0, a / n_in, NA_real_),
          n_out = n_out, events_out = c, rate_out = ifelse(n_out > 0, c / n_out, NA_real_),
          OR_pos = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }
      
      # Haldane–Anscombe correction if any zero cells
      zero_cell <- any(c(a, b, c, d) == 0L)
      a2 <- if (zero_cell) a + 0.5 else a
      b2 <- if (zero_cell) b + 0.5 else b
      c2 <- if (zero_cell) c + 0.5 else c
      d2 <- if (zero_cell) d + 0.5 else d
      
      # OR (positive outcome): (a/b) / (c/d) = (a*d) / (b*c)
      OR <- (a2 * d2) / (b2 * c2)
      
      # Wald CI on log OR
      se_logOR <- sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
      CI_low <- exp(log(OR) - 1.96 * se_logOR)
      CI_high <- exp(log(OR) + 1.96 * se_logOR)
      
      # Fisher's exact p-value (uses uncorrected counts)
      pval <- tryCatch(
        stats::fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
        error = function(e) NA_real_
      )
      
      res_long[[length(res_long) + 1]] <- data.frame(
        bin = bin_label, group = g,
        n_in = n_in, events_in = a, rate_in = a / n_in,
        n_out = n_out, events_out = c, rate_out = c / n_out,
        OR_pos = OR, CI_low = CI_low, CI_high = CI_high,
        p_value = pval,
        stringsAsFactors = FALSE
      )
    }
  }
  
  res_long <- do.call(rbind, res_long)
  
  # Nicely formatted "OR [L, U]; p=" string for display
  fmt_num <- function(x, d = digits) {
    ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = d))
  }
  res_long$OR_CI_p <- ifelse(
    is.na(res_long$OR_pos),
    NA_character_,
    paste0(
      fmt_num(res_long$OR_pos), " [",
      fmt_num(res_long$CI_low), ", ",
      fmt_num(res_long$CI_high), "]; p=",
      formatC(res_long$p_value, format = "g", digits = max(3, digits))
    )
  )
  
  # Build a wide table with one column per group using the formatted string
  # (Rows ordered by bin sequence)
  # We'll keep just one row per (bin, group) with the formatted metric.
  bins_order <- unique(res_long$bin)
  groups_order <- levels(data[[group_col]])
  
  wide_format <- do.call(
    cbind,
    c(list(data.frame(Bin = bins_order, stringsAsFactors = FALSE)),
      lapply(groups_order, function(g) {
        v <- res_long$OR_CI_p[match(
          paste(bins_order, g),
          paste(res_long$bin, res_long$group)
        )]
        data.frame(v, stringsAsFactors = FALSE)
      }))
  )
  names(wide_format) <- c("Bin", groups_order)
  
  # Also return a numeric wide table with raw ORs (optional)
  wide_OR <- do.call(
    cbind,
    c(list(data.frame(Bin = bins_order, stringsAsFactors = FALSE)),
      lapply(groups_order, function(g) {
        v <- res_long$OR_pos[match(
          paste(bins_order, g),
          paste(res_long$bin, res_long$group)
        )]
        data.frame(v, stringsAsFactors = FALSE)
      }))
  )
  names(wide_OR) <- c("Bin", groups_order)
  
  # Add helpful attributes
  attr(res_long, "note") <- paste0(
    "Within each Group, OR compares odds of ", outcome, " == ", positive_level,
    " for IN-bin vs OUT-of-bin. 95% CI via Wald on log-OR; p-value via Fisher's exact test. ",
    "Haldane–Anscombe +0.5 applied when a zero cell occurs."
  )
  
  list(
    long = res_long,          # one row per (bin, group) with counts, rates, OR, CI, p
    wide_formatted = wide_format,  # display table: OR [L, U]; p=
    wide_OR = wide_OR         # numeric ORs only
  )
}



RRtable <- function(data,
                    outcome,
                    group_col = "Group",
                    positive_level = 1,
                    pred_var = "pred",
                    bins = c(10, 20, 30, 40, 50, 60, 70, 80, 90),
                    use_existing_cols = TRUE,
                    within_group_bins = FALSE,   # NEW: compute thresholds within each group
                    compare_to = c("rest", "bottom"),  # NEW: compare to the rest (default) or exact bottom X%
                    digits = 2) {
  compare_to <- match.arg(compare_to)
  
  # Checks
  stopifnot(outcome %in% names(data))
  stopifnot(group_col %in% names(data))
  if (!pred_var %in% names(data) && !use_existing_cols) {
    stop("`pred_var` not found and `use_existing_cols = FALSE`.")
  }
  
  # Ensure Group is factor (preserve order if already a factor)
  data[[group_col]] <- as.factor(data[[group_col]])
  groups <- levels(data[[group_col]])
  bin_labels <- paste0("Top ", bins, "%")
  
  # Helper to build a robust logical for "outcome == positive_level"
  is_positive <- function(y) {
    if (is.factor(y) || is.character(y)) {
      as.character(y) == as.character(positive_level)
    } else {
      # numeric/logical/integer
      y == suppressWarnings(as.numeric(positive_level))
    }
  }
  
  # Helper: if reusing existing predtopXX columns
  get_existing_bin <- function(bin_pct) {
    colname <- paste0("predtop", bin_pct)
    if (!colname %in% names(data)) {
      stop("Requested to use existing '", colname,
           "' but it was not found in `data`.")
    }
    v <- data[[colname]]
    if (!is.logical(v)) v <- as.logical(v)
    v
  }
  
  res_long <- list()
  
  for (i in seq_along(bins)) {
    bin_pct <- bins[i]
    bin_label <- bin_labels[i]
    
    for (g in groups) {
      # Subset mask for the current group
      in_group <- !is.na(data[[group_col]]) & data[[group_col]] == g
      
      # Build in-bin logical vector
      if (use_existing_cols && !within_group_bins && compare_to == "rest") {
        # Fast path: reuse existing global predtopXX as IN-bin
        inbin_full <- get_existing_bin(bin_pct)
      } else {
        # Need to compute from `pred_var`
        if (!pred_var %in% names(data)) {
          stop("Missing `pred_var` to compute thresholds.")
        }
        
        pred <- data[[pred_var]]
        
        if (within_group_bins) {
          # thresholds computed within group g
          pred_g <- pred[in_group]
          thr_top <- stats::quantile(pred_g, probs = 1 - bin_pct/100, na.rm = TRUE, names = FALSE)
          
          if (compare_to == "rest") {
            inbin_full <- pred >= thr_top
          } else { # compare_to == "bottom"
            thr_bot <- stats::quantile(pred_g, probs = bin_pct/100, na.rm = TRUE, names = FALSE)
            # top X% vs bottom X% (others excluded)
            inbin_full <- pred >= thr_top
            outbin_full <- pred <= thr_bot
          }
          
        } else {
          # global thresholds across all data
          thr_top <- stats::quantile(pred, probs = 1 - bin_pct/100, na.rm = TRUE, names = FALSE)
          
          if (compare_to == "rest") {
            inbin_full <- pred >= thr_top
          } else {
            thr_bot <- stats::quantile(pred, probs = bin_pct/100, na.rm = TRUE, names = FALSE)
            inbin_full <- pred >= thr_top
            outbin_full <- pred <= thr_bot
          }
        }
      }
      
      # Now restrict to rows that are in this group and have all required info
      y_all <- data[[outcome]]
      keep <- in_group & !is.na(y_all) & !is.na(inbin_full)
      
      # If we're comparing top X% vs bottom X%, we also need outbin_full
      if (compare_to == "bottom") {
        keep <- keep & !is.na(outbin_full)
      }
      
      if (!any(keep)) {
        res_long[[length(res_long) + 1]] <- data.frame(
          bin = bin_label, group = g,
          n_in = NA_integer_, events_in = NA_integer_, risk_in = NA_real_,
          n_out = NA_integer_, events_out = NA_integer_, risk_out = NA_real_,
          RR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }
      
      y <- y_all[keep]
      y_pos <- is_positive(y)
      
      inbin <- inbin_full[keep]
      
      if (compare_to == "rest") {
        outbin <- !inbin
      } else {
        # exact bottom X% (exclude the middle)
        outbin <- outbin_full[keep]
      }
      
      # If either side has no obs, RR undefined
      if (sum(inbin, na.rm = TRUE) == 0L || sum(outbin, na.rm = TRUE) == 0L) {
        n_in  <- sum(inbin,  na.rm = TRUE)
        n_out <- sum(outbin, na.rm = TRUE)
        a <- sum(inbin & y_pos,  na.rm = TRUE)
        c <- sum(outbin & y_pos, na.rm = TRUE)
        res_long[[length(res_long) + 1]] <- data.frame(
          bin = bin_label, group = g,
          n_in = n_in, events_in = a, risk_in = ifelse(n_in > 0, a / n_in, NA_real_),
          n_out = n_out, events_out = c, risk_out = ifelse(n_out > 0, c / n_out, NA_real_),
          RR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }
      
      # 2x2 counts for this group/bin
      a <- sum(inbin  &  y_pos, na.rm = TRUE)  # in & pos
      b <- sum(inbin  & !y_pos, na.rm = TRUE)  # in & neg
      c <- sum(outbin &  y_pos, na.rm = TRUE)  # out & pos
      d <- sum(outbin & !y_pos, na.rm = TRUE)  # out & neg
      
      n_in  <- a + b
      n_out <- c + d
      
      # Haldane–Anscombe if any zero cell
      zero_cell <- any(c(a, b, c, d) == 0L)
      a2 <- if (zero_cell) a + 0.5 else a
      b2 <- if (zero_cell) b + 0.5 else b
      c2 <- if (zero_cell) c + 0.5 else c
      d2 <- if (zero_cell) d + 0.5 else d
      
      risk_in  <- a / n_in
      risk_out <- c / n_out
      
      # RR using corrected counts when necessary
      RR <- (a2 / (a2 + b2)) / (c2 / (c2 + d2))
      
      # Katz log-method SE and 95% CI
      se_logRR <- sqrt( (1 / a2) - (1 / (a2 + b2)) + (1 / c2) - (1 / (c2 + d2)) )
      CI_low <- exp(log(RR) - 1.96 * se_logRR)
      CI_high <- exp(log(RR) + 1.96 * se_logRR)
      
      # Fisher's exact p on the original counts (layout by row)
      pval <- tryCatch(
        stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value,
        error = function(e) NA_real_
      )
      
      res_long[[length(res_long) + 1]] <- data.frame(
        bin = bin_label, group = g,
        n_in = n_in, events_in = a, risk_in = risk_in,
        n_out = n_out, events_out = c, risk_out = risk_out,
        RR = RR, CI_low = CI_low, CI_high = CI_high,
        p_value = pval,
        stringsAsFactors = FALSE
      )
    }
  }
  
  res_long <- do.call(rbind, res_long)
  
  # Formatting helpers
  fmt_num <- function(x, d = digits) {
    ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = d))
  }
  
  res_long$RR_CI_p <- ifelse(
    is.na(res_long$RR),
    NA_character_,
    paste0(
      fmt_num(res_long$RR), " [",
      fmt_num(res_long$CI_low), ", ",
      fmt_num(res_long$CI_high), "]; p=",
      formatC(res_long$p_value, format = "g", digits = max(3, digits))
    )
  )
  
  # Wide display table (formatted)
  bins_order <- unique(res_long$bin)
  groups_order <- levels(data[[group_col]])
  
  wide_formatted <- do.call(
    cbind,
    c(list(data.frame(Bin = bins_order, stringsAsFactors = FALSE)),
      lapply(groups_order, function(g) {
        v <- res_long$RR_CI_p[match(
          paste(bins_order, g),
          paste(res_long$bin, res_long$group)
        )]
        data.frame(v, stringsAsFactors = FALSE)
      }))
  )
  names(wide_formatted) <- c("Bin", groups_order)
  
  # Wide numeric RR only
  wide_RR <- do.call(
    cbind,
    c(list(data.frame(Bin = bins_order, stringsAsFactors = FALSE)),
      lapply(groups_order, function(g) {
        v <- res_long$RR[match(
          paste(bins_order, g),
          paste(res_long$bin, res_long$group)
        )]
        data.frame(v, stringsAsFactors = FALSE)
      }))
  )
  names(wide_RR) <- c("Bin", groups_order)
  
  attr(res_long, "note") <- paste0(
    "Within each ", group_col, ", RR compares risk of ", outcome, " == ", positive_level,
    if (compare_to == "rest")
      " for IN-bin (Top X%) vs OUT-of-bin (the rest). "
    else
      " for Top X% vs Bottom X% (middle excluded). ",
    "95% CI via Katz log method; p-value via Fisher's exact test. ",
    "Haldane–Anscombe +0.5 applied if any zero cell occurs.",
    if (within_group_bins)
      " Thresholds computed within each group."
    else
      " Thresholds computed globally."
  )
  
  list(
    long = res_long,             # one row per (bin, group) with counts, risks, RR, CI, p
    wide_formatted = wide_formatted,  # display table: RR [L, U]; p=
    wide_RR = wide_RR            # numeric RRs only
  )
}



##########################################################################
# Step 9 - Model Logistic Regression, Confusion Matrix, and OR/RR tables #
##########################################################################

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
#   - PrCa (Prostate Cancer diagnosis after assessment centre)
#   - PrCa_2yrs (Prostate Cancer diagnosis within 2 years after assessment centre)
#   - PrCa_5yrs (Prostate Cancer diagnosis within 2 years after assessment centre)
#   - PrCa_10yrs (Prostate Cancer diagnosis within 2 years after assessment centre)
#
#   - PrCa_actionable (Prostate Cancer diagnosis after assessment centre that satisfies "Actionable" criteria within 2 years of diagnosis)
#   - PrCa_actionable_2yrs (Prostate Cancer diagnosis within 2 years after assessment centre that satisfies "Actionable" criteria within 2 years of diagnosis)
#   - PrCa_actionable_5yrs (Prostate Cancer diagnosis within 2 years after assessment centre that satisfies "Actionable" criteria within 2 years of diagnosis)
#   - PrCa_actionable_10yrs (Prostate Cancer diagnosis within 2 years after assessment centre that satisfies "Actionable" criteria within 2 years of diagnosis)
#
#   - PrCa_severe (Prostate Cancer diagnosis after assessment centre that satisfies "Actionable" criteria within 2 years of diagnosis)
#   - PrCa_severe_2yrs (Prostate Cancer diagnosis within 2 years after assessment centre that satisfies "Actionable" criteria within 2 years of diagnosis)
#   - PrCa_severe_5yrs (Prostate Cancer diagnosis within 2 years after assessment centre that satisfies "Actionable" criteria within 2 years of diagnosis)
#   - PrCa_severe_10yrs (Prostate Cancer diagnosis within 2 years after assessment centre that satisfies "Actionable" criteria within 2 years of diagnosis)
#
# Set "predictor" to either: 
#   - multiethnicGRS (Conti's GRS with pan-Ancestry weights)
#   - EuropeanGRS (Conti's GRS with European-specific weights)
#   - AfricanGRS (Conti's GRS with African-specific weights)
#   - East_AsianGRS (Conti's GRS with East Asian-specific weights)
#   - HispanicGRS (Conti's GRS with Hispanic-specific weights)
#   - adjustedGRS (Conti's GRS adjusted for ancestry probability)
#   - WangmultiethnicGRS (Wang's GRS with multi-ethnic weights)
#   - WangEuropeanGRS (Wang's GRS with European-specific weights)
#   - WangAfricanGRS (Wang's GRS with African-specific weights)
#   - WangEast_AsianGRS (Wang's GRS with East Asian-specific weights)
#   - WangHispanicGRS (Wang's GRS with Hispanic-specific weights)
#
#
# Set "covariates" to either: (or add multiple using + between covariates)
#   - event_age (Age at symptom presentation)
#   - rs72725854_T (carrier status of rs72725854 risk allele)
#


model <- run_logreg(data = PCa_iv_covariates_GRS_predhorizon_WhiteOnly,
                    outcome = "PrCa",
                    predictor = "multiethnicGRS",
                    covariates = "event_age",
                    plot_roc = TRUE)



matrix <- confusion_matrix(data = model$data, 
                           outcome = "PrCa_10yrs", 
                           cutoff_value = 0.20, 
                           positive_level = 1, negative_level = 0)

OR_table <- ORtable(
  data = model$data,
  outcome = "PrCa_10yrs",
  group_col = "ethnicity_group_narrow",       # <- your 6-level grouping variable
  positive_level = 1,        # 1 denotes positive outcome
  use_existing_cols = TRUE   # rely on predtop10..predtop90 already in data
)

print(OR_table$wide_formatted)


RR_table <- RRtable(
  data = model$data,
  outcome = "PrCa_10yrs",
  group_col = "ethnicity_group_narrow",
  positive_level = 1,
  use_existing_cols = TRUE
)

# View formatted RR table
print(RR_table$wide_formatted)
