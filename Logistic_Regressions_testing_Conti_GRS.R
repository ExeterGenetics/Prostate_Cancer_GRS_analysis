######### Logistic Regression Application - testing Conti's GRSs #########

## Setup ##

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

library(ggplot2)
install.packages("pROC")
library(pROC)
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

dxdownload("Callum/GRSs/Conti_multiethnicGRS_267.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_MULTI"
dxdownload("Callum/GRSs/Conti_EuropeanGRS_265.tsv") 		# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_EUR"
dxdownload("Callum/GRSs/Conti_AfricanGRS_246.tsv")			# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_AFR"
dxdownload("Callum/GRSs/Conti_East_AsianGRS_222.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_EAS"
dxdownload("Callum/GRSs/Conti_HispanicGRS_253.tsv")		# Calculated Conti GRS for all eligible participants using "Conti_GRS_267.R" with "OR_column" set to "OR_HIS"

dxdownload("Callum/GRSs/Wang_multi_ethnic.pgs.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/Wang_European.pgs.tsv") 		# Calculated Wang GRS for all eligible participants using "Conti_script" 
dxdownload("Callum/GRSs/Wang_African.pgs.tsv")			# Calculated Wang GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/Wang_East_Asian.pgs.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_script"
dxdownload("Callum/GRSs/Wang_Hispanic.pgs.tsv")		# Calculated Wang GRS for all eligible participants using "Conti_script"

dxdownload("Callum/GRSs/OR_adjustedGRS_267.tsv")	# Calculated using "Create_adjustedGRS_weighting_Conti_odds_ratios_by_AncestryProbability2.R" script

######################################################################
# Step 1 - Fetch PrCa Cases and collapse into earliest epistart/date #
######################################################################

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
  select(c("eid", "diag_icd10", "assess_date_initial", "epistart", "epiend"))

PCaCases_HES <- bind_rows(PCaCases_HES, PCaCases_ICD9) # Now we have all, after adding old ICD9 diagnoses

PCaCases_cancerregistry<-read_cancer('C61') # Creates a dataframe of all recorded prostate cancer diagnoses in cancer registry records
PCaCases_death<-read_death('C61') # Creates a dataframe of all recorded prostate cancer deaths in death records

PCaCases_HES_earliest <- PCaCases_HES %>%
  dplyr::mutate(epistart = as.Date(epistart)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    epistart = if (all(is.na(epistart))) NA else min(epistart, na.rm = TRUE),
    assess_date_initial = first(assess_date_initial),
    diag_icd10 = first(diag_icd10),
    .groups = "drop"
  )

PCaCases_cancerregistry_earliest <- PCaCases_cancerregistry %>%
  dplyr::mutate(date = as.Date(date)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    date = if (all(is.na(date))) NA else min(date, na.rm = TRUE),
    assess_date_initial = first(assess_date_initial),
    ICD10 = first(ICD10),
    .groups = "drop"
  )

PCaCases_death_earliest <- PCaCases_death %>%
  dplyr::mutate(date_of_death = as.Date(date_of_death)) %>%
  dplyr::group_by(eid) %>%
  dplyr::summarise(
    date_of_death = if (all(is.na(date_of_death))) NA else min(date_of_death, na.rm = TRUE),
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
  rename("assess_date_initial_ICD10" = "assess_date_initial.x", "assess_date_initial_cr" = "assess_date_initial.y", "assess_date_initial_death" = "assess_date_initial")


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


##################################################################################################
# Step 2 - Fetch cancer death, chemotherapy, surgery, radiotherapy, and androgen therapy records #
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
# Step 3 - Load in independent variables, GRSs, and covariates from DNAnexus project #
######################################################################################

# Load in independent variable, covariates, and GRSs

iv <- read_csv("imputed_rs72725854.csv") # loads independent variable dataset
iv <- iv %>%
  select(c("eid", "rs72725854_G"))

covariates <- read_csv("Age_Sex_PRS_GA.csv")
covariates <- covariates %>% 
  rename('Age' = 'p21022', 'Sex' = 'p31', 'PRS' = 'p26267', 'Enhanced_PRS' = 'p26268', 'Genetic_sex' = 'p22001', 'Genomic_ancestry' = 'p30079')


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

All_Wang_GRS <- merge(WangmultiethnicGRS, WangEuropeanGRS, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangAfricanGRS, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangEast_AsianGRS, by = "eid")
All_Wang_GRS <- merge(All_Wang_GRS, WangHispanicGRS, by = "eid")

## Merge all GRSs into one

All_GRS <- merge(All_Conti_GRS, All_Wang_GRS, by = "eid", all = T)

#########################################################################################
# Step 4 - Merge all variables/covariates into one dataframe with Prostate Cancer cases #
#########################################################################################

iv_covariates <- merge(iv, covariates, by = "eid", all.y = T)

PCa_iv_covariates <- merge(iv_covariates, PCaCases_prediagnosis, by = "eid", all = TRUE)

PCa_iv_covariates_GRS <- merge(PCa_iv_covariates, All_GRS, by = "eid", all.x = TRUE)

PCa_iv_covariates_GRS_severity <- merge(PCa_iv_covariates_GRS, actionable_criteria, by = 'eid', all = T)

# Remove anybody who is female or who lacks GRS data (optionally, also remove anyone pre-diagnosed)

PCa_iv_covariates_GRS_clean <- PCa_iv_covariates_GRS_severity %>%
  dplyr::filter(
    Sex == 'Male',
    !is.na(ContimultiethnicGRS),
    pre_diagnosed == FALSE | is.na(pre_diagnosed) ## to remove pre-diagnosed patients
  )

# Sanity Check - how many patients in PCa_iv_covariates were female or lacked GRS data? Does it match the difference in n between PCa_iv_covariates and PCa_iv_covariates_clean ?

PCa_iv_covariates_GRS_severity %>%
  filter(
    Sex == "Female" |
      is.na(ContimultiethnicGRS) |
      pre_diagnosed == TRUE
  ) %>%
  summarise(n = n())

# Set Prostate Cancer cases to 1, controls to 0. 

PCa_iv_covariates_GRS_clean <- PCa_iv_covariates_GRS_clean %>%
  dplyr::mutate(
    PrCa = if_else(!is.na(pre_diagnosed), 1L, 0L, missing = 0L),
  )

############################################# (these should have already been filtered out, 
# Important - remove withdrawn participants #  since the withdrawn participants no longer have 
#############################################  a GRS calculated, but still worth making sure)

exclude_withdrawn=function(df){
  system('dx download Callum/Withdrawals/withdrawn_20260310.csv --overwrite') ## This file is a list of participants who withdrew from the Biobank up to the date 10th March 2026. This was sent from the UK Biobank team via email to members of approved applications
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
# Step 5 - Set prediction horizons, including for general/actionable/severe PrCa #
##################################################################################

library(lubridate)

PCa_iv_covariates_GRS_predhorizon <- PCa_iv_covariates_GRS_clean %>%
  dplyr::mutate(
    PrCa_post_assessment  = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial
    ),
    PrCa_2yrs  = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial &
        earliest_PrCa_date <= assess_date_initial %m+% lubridate::years(2)
    ),
    PrCa_5yrs  = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial &
        earliest_PrCa_date <= assess_date_initial %m+% lubridate::years(5)
    ),
    PrCa_10yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial &
        earliest_PrCa_date <= assess_date_initial %m+% lubridate::years(10)
    ),
    PrCa_severe = as.integer(
      !is.na(earliest_PrCa_date) & (
        (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
          (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2))
      )
    ),
    PrCa_severe_2yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial &
        earliest_PrCa_date <= assess_date_initial %m+% lubridate::years(2) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    ),
    PrCa_severe_5yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial &
        earliest_PrCa_date <= assess_date_initial %m+% lubridate::years(5) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    ),
    PrCa_severe_10yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial &
        earliest_PrCa_date <= assess_date_initial %m+% lubridate::years(10) & (
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
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial &
        earliest_PrCa_date <= assess_date_initial %m+% lubridate::years(2) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(surgery_opdate) & surgery_opdate >= earliest_PrCa_date & surgery_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(radio_opdate) & radio_opdate >= earliest_PrCa_date & radio_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(androgen_opdate) & androgen_opdate >= earliest_PrCa_date & androgen_opdate <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    ),
    PrCa_actionable_5yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial &
        earliest_PrCa_date <= assess_date_initial %m+% lubridate::years(5) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(surgery_opdate) & surgery_opdate >= earliest_PrCa_date & surgery_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(radio_opdate) & radio_opdate >= earliest_PrCa_date & radio_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(androgen_opdate) & androgen_opdate >= earliest_PrCa_date & androgen_opdate <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    ),
    PrCa_actionable_10yrs = as.integer(
      !is.na(earliest_PrCa_date) & !is.na(assess_date_initial) &
        earliest_PrCa_date >= assess_date_initial &
        earliest_PrCa_date <= assess_date_initial %m+% lubridate::years(10) & (
          (!is.na(chemo_opdate) & chemo_opdate >= earliest_PrCa_date & chemo_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(date_of_cancer_death) & date_of_cancer_death >= earliest_PrCa_date & date_of_cancer_death <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(surgery_opdate) & surgery_opdate >= earliest_PrCa_date & surgery_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(radio_opdate) & radio_opdate >= earliest_PrCa_date & radio_opdate <= earliest_PrCa_date %m+% lubridate::years(2)) |
            (!is.na(androgen_opdate) & androgen_opdate >= earliest_PrCa_date & androgen_opdate <= earliest_PrCa_date %m+% lubridate::years(2))
        )
    )
  )

################################################################
# Step 6 - Subset data for different ethnicity/ancestry groups #
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
# Step 7 (optional) - Visually inspect GRS distribution for cases vs. controls #
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

#########################################################################################
# Step 8 - Set up logistic regression, confusion matrix, and Odds Ratio table functions #
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
#
#    Note: these top 6 GRSs use "Conti_script.R", which drops 4 SNPs by default. The numbered GRSs below only drop 2 SNPs
#
#   - ContimultiethnicGRS (Conti's GRS with pan-Ancestry weights)
#   - ContiEuropeanGRS (Conti's GRS with European-specific weights)
#   - ContiAfricanGRS (Conti's GRS with African-specific weights)
#   - ContiEast_AsianGRS (Conti's GRS with East Asian-specific weights)
#   - ContiHispanicGRS (Conti's GRS with Hispanic-specific weights)
#   - ContiadjustedGRS (Conti's GRS adjusted for ancestry probability)
#
#   - ContimultiethnicGRS267 (Conti's GRS with pan-Ancestry weights, all 267 available SNPs)
#   - ContiEuropeanGRS (Conti's GRS with European-specific weights)
#   - ContiAfricanGRS (Conti's GRS with African-specific weights)
#   - ContiEast_AsianGRS (Conti's GRS with East Asian-specific weights)
#   - ContiHispanicGRS (Conti's GRS with Hispanic-specific weights)
#   - ContiadjustedGRS (Conti's GRS adjusted for ancestry probability)
#
#
# Set "covariates" to either: (or add multiple using + between covariates)
#   - Age (Age at assessment centre visit)
#   - rs72725854_T (carrier status of rs72725854 risk allele)
#


model <- run_logreg(data = PCa_iv_covariates_GRS_predhorizon_BlackOnly,
                    outcome = "PrCa_10yrs",
                    predictor = "WangAfricanGRS",
                    covariates = "Age",
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
