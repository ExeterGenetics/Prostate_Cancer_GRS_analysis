library(tidyverse)
system('dx download Callum/Derived_datasets/family_history_participant.csv') # Dataset of p20107_i0, p20107_i1, p20107_i2, p20107_i3, p20110_i0, p20110_i1, p20110_i2, p20110_i3, p20111_i0, p20111_i1, p20111_i2, p20111_i3, created using Cohort browser
FH=read.csv('family_history_participant.csv')

all_ids <- FH %>% distinct(eid)

# Step 2: Extract diseases into binary format
FH_diseases <- FH %>%
  pivot_longer(cols = starts_with("p"), names_to = "p", values_to = "disease") %>%
  filter(!is.na(disease), disease != "null", disease != "") %>%
  separate_rows(disease, sep = "\\|") %>%
  filter(disease != "", !is.na(disease)) %>%
  distinct(eid, disease) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = disease, values_from = value, values_fill = list(value = 0)) %>%
  rename_with(~paste0("FH_", .x), .cols = -eid)  # add FH_ prefix to disease columns only

# Step 3: Join back to retain all eids
FH_reformat <- all_ids %>% left_join(FH_diseases, by = "eid")

FH_reformat <- FH_reformat %>%
  rename_with(~gsub(" ", "_", .x))

FH_reformat$n_FH_answers=(FH[,2]=='')+(FH[,3]=='')+(FH[,3]=='')+(FH[,4]=='')

# Step 4: select only eid, breast cancer FH, prostate cancer FH, then create a combined column

FH_PrCa_BrCa <- FH_reformat %>%
  dplyr::select(c("eid", "FH_Prostate_cancer", "FH_Breast_cancer"))

FH_PrCa_BrCa <- FH_PrCa_BrCa %>%
  dplyr::mutate(FH_PrCa_BrCa = if_else(FH_Prostate_cancer == "1" | FH_Breast_cancer == "1" , 1L, 0L, missing = 0L),
  )

# Step 5: Sanity check - how many participants with FH of breast cancer or prostate cancer? Does it match FH_PrCa_BrCa column?
## Should both be 90201

FH_PrCa_BrCa %>%
  filter(
    FH_Prostate_cancer == "1" | FH_Breast_cancer == "1" 
  ) %>%
  summarise(n = n())

FH_PrCa_BrCa %>%
  filter(
    FH_PrCa_BrCa == "1"
  ) %>%
  summarise(n = n())


detach("package:tidyverse", unload = TRUE)

#############################################
# Important - remove withdrawn participants #
#############################################

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

FH_PrCa_BrCa <- exclude_withdrawn(FH_PrCa_BrCa)

write.csv(FH_PrCa_BrCa, file = "FH_PrCa_BrCa.csv")
system(paste("dx upload", "FH_PrCa_BrCa.csv"))