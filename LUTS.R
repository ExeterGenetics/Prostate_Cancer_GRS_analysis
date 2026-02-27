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

system("dx download Callum/LUTS/Green2022supplementarytable1.csv") # This is the supplementary table 1 from Harry's LUTS paper, converted to .csv format using excel. Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC9553867/
Green2022supplementarytable1 <- read.csv("Green2022supplementarytable1.csv")

LUTS_codes <- Green2022supplementarytable1$read_3 |> as.character()
print(LUTS_codes)

## Using these codes, select corresponding GP records

# Method 1) Exactly as codes are written, all full stops and ellipses included

LUTS_gp_records1 <- read_GP(c(
  "1A27.", "K16y8", "XaNFc", "X30Ni", "1A1Z.", "1A11.", "1A1..", "R084.", "1A12.", "R084z", "R0840",
  "1A1..", "1A1..", "1A34.", "1A34.", "R083z", "R083.", "1AZ6.", "1AZ60", "1AZ61", "1AZ62", "1A2..",
  "1A2Z.", "R086z", "8D7..", "R08..", "R08zz", "66K3.", "1A4..", "1A...", "1A...", "1AZ..", "1AZZ.",
  "1AH1.", "R086.", "R08z.", "Kz...", "Ryu4.", "XaB9O", "XaXHi", "XaXHj", "XaXHk", "Xa96j", "1A13.",
  "R0842", "1A33.", "1A31.", "1A3..", "1A3Z.", "1A3..", "R0861", "R0863", "R0860", "317C.", "1A37.",
  "1A36.", "XaD2w", "X77SF", "X76Y0", "R082.", "R0824", "1A32.", "K196.", "R0820", "1A32.", "R0822",
  "1A25.", "R0862", "1A25.", "R15y0", "B7C20", "14270", "ZV104", "B834.", "1J08.", "B58y5", "B8340",
  "B46..", "XaC0j", "Xa3fu", "XaXGk", "XaFwo", "XaKyV"
))

LUTS_gp_records1 <- LUTS_gp_records1 %>%                                          ## Change all blanks to NA
  mutate(across(all_of(c("read_2", "read_3")),
                ~ na_if(str_trim(as.character(.)), "")))


LUTS_gp_records1$read_code <- coalesce(LUTS_gp_records1$read_2, LUTS_gp_records1$read_3)



LUTS_gp_records1_over40 <- LUTS_gp_records1 %>%
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

LUTS_gp_records2 <- read_GP(c(
  "1A27", "K16y8", "XaNFc", "X30Ni", "1A1Z", "1A11", "1A1", "R084", "1A12", "R084z", "R0840",
  "1A1", "1A1", "1A34", "1A34", "R083z", "R083", "1AZ6", "1AZ60", "1AZ61", "1AZ62", "1A2",
  "1A2Z", "R086z", "8D7", "R08", "R08zz", "66K3", "1A4", "1A", "1A", "1AZ", "1AZZ",
  "1AH1", "R086", "R08z", "Kz", "Ryu4", "XaB9O", "XaXHi", "XaXHj", "XaXHk", "Xa96j", "1A13",
  "R0842", "1A33", "1A31", "1A3", "1A3Z", "1A3", "R0861", "R0863", "R0860", "317C", "1A37",
  "1A36", "XaD2w", "X77SF", "X76Y0", "R082", "R0824", "1A32", "K196", "R0820", "1A32", "R0822",
  "1A25", "R0862", "1A25", "R15y0", "B7C20", "14270", "ZV104", "B834", "1J08", "B58y5", "B8340",
  "B46", "XaC0j", "Xa3fu", "XaXGk", "XaFwo", "XaKyV"
))

LUTS_gp_records2 <- LUTS_gp_records2 %>%                                          ## Change all blanks to NA
  mutate(across(all_of(c("read_2", "read_3")),
                ~ na_if(str_trim(as.character(.)), "")))


LUTS_gp_records2$read_code <- coalesce(LUTS_gp_records2$read_2, LUTS_gp_records2$read_3)



starts_with_any <- function(x, prefixes) {
  # returns a logical vector same length as x
  Reduce(`|`, lapply(prefixes, function(p) startsWith(x, p)))
}


LUTS_gp_records2_over40 <- LUTS_gp_records2 %>%
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

