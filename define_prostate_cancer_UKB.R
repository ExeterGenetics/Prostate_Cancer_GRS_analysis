#####################################################################################################################################
# This script defines prostate cancer case inclusion, actionable/severe treatment criteria, and control exclusion in the UK Biobank #
#####################################################################################################################################

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

## Authors: Callum Webb, Wawa Chen, and Harry Green

## Configuration - set phenotype codes here

# Case inclusion

ICD9_case_inclusion_codes <- c('185')                                            ## 185: Malignant neoplasm of prostate
ICD10_case_inclusion_codes <- c('C61')                                           ## C61: Malignant neoplasm of prostate

# Control exclusion

ICD9_control_exclusion_codes <- c('185',                                         ## 185: Malignant neoplasm of prostate
                                  'V104',                                        ## V10.4: Personal history of malignant neoplasm of genital organs
                                  '2334',                                        ## 233.4: Carcinoma in situ of the prostate 
                                  '2365',                                        ## 236.5: Neoplasm of uncertain behavior of the prostate 
                                  '6023')                                        ## 602.3: Dysplasia of the prostate 

ICD10_control_exclusion_codes <- c('C61',                                        ## C61: Malignant neoplasm of prostate
                                   'Z854',                                       ## Z85.4: Personal History of malignant neoplasm of genital organs
                                  #'R972',                                       ## R97.2: Elevated prostate specific antigen [PSA] (note: no results returned)
                                   'D075',                                       ## D07.5: Carcinoma in situ of prostate
                                   'D400',                                       ## D40.0: Neoplasm of uncertain behavior of prostate
                                   'N423')                                       ## N42.3: Dysplasia of prostate   

OPCS_control_exclusion_codes <- c('M61',                                         ## M61: Prostatectomy
                                  'M611',                                        ## M61.1: Radical prostatectomy
                                  'M612',                                        ## M61.2: Retropubic Prostatectomy
                                  'M613',                                        ## M61.3: Transvesical Prostatectomy
                                  'M614',                                        ## M61.4: Perineal Prostatectomy
                                 #'X65',                                         ## X65: Radiotherapy Delivery
                                 #'X67',                                         ## X67: Preparation of radiotherapy
                                 #'X68',                                         ## X68: Brachytherapy preparation
                                 #'X69',                                         ## X69: Other radiotherapy
                                 #'Y91',                                         ## Y91: External beam radiotherapy
                                  'M706',                                        ## M70.6 Radioactive seed implantation into prostate
                                 #'Y35',                                         ## Y35: Introduction Material Radioactive Removable NOC
                                 #'Y36',                                         ## Y36: Introduction Material Non-removable NOC
                                 #'T856',                                        ## T85.6: Block dissection of pelvic lymph nodes
                                 #'M702',                                        ## M70.2: Perineal needle biopsy of prostate
                                 #'M703',                                        ## M70.3: Rectal needle biopsy of prostate
                                 #'N04',                                         ## N04: Orchidectomy (note: no results returned)
                                  'M65',                                         ## M65: Endoscopic resection of prostate
                                  'M68',                                         ## M68: Endoscopic insertion of prosthesis into prostate
                                  'M671',                                        ## M67.1: Endoscopic cryotherapy to lesion of prostate
                                  'M711',                                        ## M71.1: High intensity focused ultrasound of prostate
                                  'M712')                                        ## M71.2: Implantation of radioactive substance into prostate

# Procedures (to define actionable/severe cases)

OPCS_chemotherapy_codes <- c('X70',                                              ## X70: Procurement of drugs for chemotherapy for neoplasm in Bands 1-5
                             'X71',                                              ## X71: Procurement of drugs for chemotherapy for neoplasm in Bands 6-10
                             'X72',                                              ## X72: Delivery of chemotherapy for neoplasm
                             'X73',                                              ## X73: Delivery of oral chemotherapy for neoplasm
                             'X74')                                              ## X74: Other chemotherapy drugs

OPCS_surgery_codes <- c('M61',                                                   ## M61: Open excision of prostate
                        'M611',                                                  ## M61.1: Total excision of prostate and capsule of prostate (Radical prostatectomy)
                        'M612',                                                  ## M61.2: Retropubic Prostatectomy
                        'M613',                                                  ## M61.3: Transvesical Prostatectomy
                        'M614')                                                  ## M61.4: Perineal Prostatectomy

OPCS_radiotherapy_codes <- c('X65',                                              ## X65: Radiotherapy Delivery
                          'X67',                                                 ## X67: Preparation for external beam radiotherapy
                          'X68',                                                 ## X68: Preparation for brachytherapy
                          'X69',                                                 ## X69: Other radiotherapy
                          'Y91',                                                 ## Y91: External beam radiotherapy
                          'M706')                                                ## M70.6 Radioactive seed implantation into prostate

OPCS_androgen_therapy_codes <- c('X741',                                         ## X74.1: Cancer hormonal treatment drugs Band 1
                                 'X383',                                         ## X38.3: Injection of hormone for local action NEC
                                 'S525',                                         ## S52.5: Insertion of hormone into subcutaneous tissue
                                 'S526',                                         ## S52.6: Replacement of hormone in subcutaneous tissue
                                 'X376')                                         ## X37.6: Intramuscular hormone therapy

####################################################
#------ Step 1: Define prostate cancer cases ------#
####################################################

PCaCases_HES<-read_ICD10(ICD10_case_inclusion_codes) # Creates a dataframe of *almost* all recorded prostate cancer diagnoses in HES records
# Z85.46 (Z8546 in UKB format) is "Personal history of malignant neoplasm of prostate" and R97.21 (R9721 in UKB format) is 
# "Elevated prostate specific antigen (PSA)". Adding these to read_ICD10 does not add any cases

PCaCases_ICD9<-read_ICD9(ICD9_case_inclusion_codes) 
# V10.46 (V1046 in UKB format) is "Personal history of malignant neoplasm of prostate". Adding this to read_ICD9 does not add any cases
PCaCases_ICD9 <- PCaCases_ICD9 %>%
  dplyr::mutate(
    diag_icd10 = case_when(
      diag_icd9 %in% c("185", "1851", "1852", "1853", "1854", "1855", "1856", "1857", "1858", "1859") ~ "C61",
      diag_icd9 %in% c("V1046") ~ "Z8546",
      TRUE ~ NA_character_
    )
  )

PCaCases_ICD9 <- PCaCases_ICD9 %>%
  select(c("eid", "diag_icd10", "assess_date_initial", "epistart", "epiend"))

PCaCases_HES <- bind_rows(PCaCases_HES, PCaCases_ICD9) # Now we have all, after adding old ICD9 diagnoses

PCaCases_cancerregistry<-read_cancer(icd9 = ICD9_case_inclusion_codes, icd10 = ICD10_case_inclusion_codes) %>% # Creates a dataframe of all recorded prostate cancer diagnoses in cancer registry records
  dplyr::mutate(
    ICD10 = case_when(
      ICD9 %in% c("185", "1851", "1852", "1853", "1854", "1855", "1856", "1857", "1858", "1859") ~ "C61",
      ICD9 %in% c("V1046") ~ "Z8546",
      ICD10 %in% c("C61") ~ "C61",
      TRUE ~ NA_character_
    )
  )

PCaCases_death<-read_death(ICD10_case_inclusion_codes) # Creates a dataframe of all recorded prostate cancer deaths in death records

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
  rename("assess_date_initial_ICD10" = "assess_date_initial.x", 
         "assess_date_initial_cr" = "assess_date_initial.y", 
         "assess_date_initial_death" = "assess_date_initial",
         "icd10_HES" = "diag_icd10",
         "icd10_cr" = "ICD10",
         "icd10_death" = "cause_icd10"
  )


# Sanity check - are the assessment centre dates the same in both HES and cancer registry, unless NA on either side?

mismatched_assessdates <- with(PCaCases_earliest,
                               (!is.na(assess_date_initial_ICD10) & !is.na(assess_date_initial_cr) & assess_date_initial_ICD10 != assess_date_initial_cr) |
                                 (!is.na(assess_date_initial_ICD10) & !is.na(assess_date_initial_death) & assess_date_initial_ICD10 != assess_date_initial_death) |
                                 (!is.na(assess_date_initial_cr) & !is.na(assess_date_initial_death) & assess_date_initial_cr != assess_date_initial_death)
)

any(mismatched_assessdates)

PCaCases_earliest$assess_date_initial <- coalesce(PCaCases_earliest$assess_date_initial_ICD10, PCaCases_earliest$assess_date_initial_cr, PCaCases_earliest$assess_date_initial_death)

# Tag patients who already had prostate cancer at assessment centre as "pre-diagnosed". These will be filtered out later 
# (any death-only PrCa cases will be treated as pre_diagnosed, as we don't have HES or Cancer Registry dates to prove that they were diagnosed post-assessment centre)
# (furthermore, death records are more up-to-date than cancer registry records, so in most cases there is simply a delay in release of cancer registry dates)

# Also make "earliest PrCa diagnosis date" which is the earliest of epistart, date, and date_of_death

PCaCases_prediagnosis <- PCaCases_earliest %>%
  dplyr::mutate(
    pre_diagnosed =
      (is.na(epistart) & is.na(date)) |
      (!is.na(epistart) & epistart <= assess_date_initial) |
      (!is.na(date) & date <= assess_date_initial),
    PrCa_case = if_else(icd10_HES == "C61" | icd10_cr == "C61" | icd10_death == "C61", 1L, 0L, missing = 0L),
    earliest_PrCa_date = pmin(epistart, date, date_of_death, na.rm = TRUE)
  )

Earliest_PrCa_diagnosis <- PCaCases_prediagnosis %>%
  select('eid', 'earliest_PrCa_date')

# Tag participants who had a prostate cancer record in HES but not in cancer registry or death records as "HES-only" (these will be excluded later) 
# also tag participants with prostate cancer record in death records but not in HES or cancer registry as "death-only" 

PCaCases_prediagnosis <- PCaCases_prediagnosis %>%
  dplyr::mutate(HES_only = !is.na(icd10_HES) & is.na(icd10_cr) & is.na(icd10_death),
                death_only = !is.na(icd10_death) & is.na(icd10_cr) & is.na(icd10_HES))










########################################################################
# ------ Step 2: Define actionable/severe prostate cancer cases ------ #
########################################################################

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

chemotherapy  <- read_OPCS(OPCS_chemotherapy_codes)
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

surgery <- read_OPCS(OPCS_surgery_codes)
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

radiotherapy <- read_OPCS(OPCS_radiotherapy_codes)
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

androgen <- read_OPCS(OPCS_androgen_therapy_codes)
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













######################################################
# ------ Step 3: Define (who aren't) controls ------ #
######################################################

## Extra exclusion criteria for controls 
## (derived from https://phekb.org/phenotype/prostate-cancer-0 see tables: https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fphekb.org%2Fsites%2Fphenotype%2Ffiles%2FPrCa%2520Phenotyping%2520Algorithm%2520codes.xlsx&wdOrigin=BROWSELINK
## with alterations to remove ICD9: V84.03 and ICD10: Z15.03 since these describe genetic susceptibility to prostate cancer)
## AND with some creative interpretation for OPCS codes, since they are given as CPT codes in the spreadsheet

icd9_exclusion_prefixes <- c(ICD9_control_exclusion_codes)

icd9_exclusion_regex <- paste0(
  "^(",
  paste(gsub("[^A-Za-z0-9]", "", toupper(icd9_exclusion_prefixes)), collapse = "|"),
  ")"
)

exclusions_ICD9 <- read_ICD9(icd9_exclusion_prefixes) %>%
  dplyr::select("eid", "diag_icd9") %>%
  dplyr::filter(
    !is.na(diag_icd9),
    grepl(icd9_exclusion_regex, gsub("[^A-Za-z0-9]", "", toupper(diag_icd9)))
  )

exclusions_ICD10 <- read_ICD10(ICD10_control_exclusion_codes) %>%
  dplyr::select("eid", "diag_icd10")

exclusions_cancerregistry_ICD10 <- read_cancer(icd10 = ICD10_control_exclusion_codes)%>%
  dplyr::select("eid", "ICD10")

exclusions_cancerregistry_ICD9 <- read_cancer(icd9 = ICD9_control_exclusion_codes)%>%
  dplyr::select("eid", "ICD9")

exclusions_death <- read_death(ICD10_control_exclusion_codes)%>%
  dplyr::select("eid", "cause_icd10")

exclusions_OPCS <- read_OPCS(OPCS_control_exclusion_codes) %>%
  dplyr::select("eid", "oper4")

possible_PrCa_cases <- merge(exclusions_ICD9, exclusions_ICD10, by = "eid", all = T)
possible_PrCa_cases <- merge(possible_PrCa_cases, exclusions_OPCS, by = "eid", all = T)
possible_PrCa_cases <- merge(possible_PrCa_cases, exclusions_cancerregistry_ICD10, by = "eid", all = T)
possible_PrCa_cases <- merge(possible_PrCa_cases, exclusions_cancerregistry_ICD9, by = "eid", all = T)
possible_PrCa_cases <- merge(possible_PrCa_cases, exclusions_death, by = "eid", all = T)

### Sanity check - did we isolate the intended codes correctly?

possible_PrCa_unique_codes <- possible_PrCa_cases %>%
  dplyr::select(diag_icd9, diag_icd10, oper4, ICD9, ICD10, cause_icd10) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "source_column",
    values_to = "code"
  ) %>%
  dplyr::filter(!is.na(code), code != "") %>%
  dplyr::distinct(source_column, code) %>%
  dplyr::arrange(source_column, code)

possible_PrCa_unique_codes

possible_PrCa_unique_codes %>%
  dplyr::count(source_column, name = "n_unique_codes")

# Collapse exclusions to one row per participant, leaving only eids. Add "possible_PrCa_case" column 
possible_PrCa_cases <- possible_PrCa_cases %>%
  group_by(eid) %>%
  dplyr::summarise(possible_PrCa_case = 1)






########################################################
#----------- Step 4: Merge all dataframes -------------#
########################################################

PrCa_cases_and_exclusions <- merge(PCaCases_prediagnosis, actionable_criteria, by = "eid", all = T)
PrCa_cases_and_exclusions <- merge(PrCa_cases_and_exclusions, possible_PrCa_cases, by = "eid", all = T)

# Define "exclude" as any control who has possible_PrCa_case == 1

PrCa_cases_and_exclusions <- PrCa_cases_and_exclusions %>%
  dplyr::mutate(
    exclude = if_else(possible_PrCa_case == 1 & PrCa_case == 0, 1L, 0L, missing = 0L) ## Tagged for exclusion, so that exclusion can occur after joining with main dataframe
  )



#######################################################################
# Step 5: Set prediction horizons to define prostate cancer endpoints #
#######################################################################

library(lubridate)

endpoints <- PrCa_cases_and_exclusions %>%
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

########################################################
# ---------- Step 6: Join with main dataframe ---------#
#-----------              (example)           ---------#
########################################################


dxdownload("Callum/Derived_datasets/Age_Sex_PRS_GA.csv")	# Dataset of Age at Recruitment (p21022), Sex (p31), Genetically-inferred Sex (p22001), Standard PRS for Prostate cancer (p26267), Enhanced PRS for Prostate cancer (p26268), Genomic Ancestry (p30079), created using the UKB-RAP cohort browser

main <- read_csv("Age_Sex_PRS_GA.csv")
main <- main %>% 
  rename('Age' = 'p21022', 'Sex' = 'p31', 'GenomicsPLC_PRS' = 'p26267', 'Enhanced_PRS' = 'p26268', 'Genetic_sex' = 'p22001', 'Genetic_similarity' = 'p30079')

PrCa_dataframe <- merge(main, endpoints, by = "eid", all.x = T)

## Remove excluded participants, females, and those with missing PRS, pre-diagnosed prostate cancer, or HES-only prostate cancer diagnoses

PrCa_dataframe <- PrCa_dataframe %>%
  dplyr::filter(
    Sex == 'Male',
    !is.na(GenomicsPLC_PRS),
    pre_diagnosed == FALSE | is.na(pre_diagnosed), ## to remove pre-diagnosed == TRUE participants
    HES_only == FALSE | is.na(HES_only), ## to remove participants with HES-only PrCa diagnoses
    exclude == 0 | is.na(exclude) ## to remove controls who meet any of the exclusion criteria for controls
  )