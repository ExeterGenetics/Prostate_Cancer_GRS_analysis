###############################################################################
# Visually inspecting how well top 10% of GRS stratifies prostate cancer risk #
###############################################################################

library(tidyverse)
library(extrafont)
library(dplyr)
library(ggplot2)
library(showtext)

library(devtools)
source_url("https://raw.githubusercontent.com/hdg204/exeteR/main/ex_theme.R")

## Step 1: Run "Logistic_Regressions_testing_Conti_GRS.R" script all the way
##         down to running a model

## Step 2: Set variables and labels

dataframe <- model$data %>%
  dplyr::filter(ethnicity_group_narrow == "White")                              # Change to view specific groups, or comment out filter to see all participants
outcome <- "PrCa_10yrs"                                                         # Set to same outcome as model
GRS <- "multiethnicGRS"
top_segment <- 10                                                               # 10 means top 10% will be treated as a distinct group

xlabel <- "Group"
ylabel <- "% with prostate cancer (10 yrs)"
title_choice <- "Prostate Cancer in Top 10% of GRS vs. Bottom 90% (White participants)"

use_general_top_group = TRUE                                                    # If TRUE, this is the top 10% of GRS among ALL participants, not just the subgroup tested


## Step 3: Visual inspection

bottom_cutoff <- 0.01*(100-top_segment)
predef_top <- paste0("top", top_segment, "_all_", GRS)                                        # Change number (multiples of 10 only) to see other quantiles

# 1) Summarize: percentage of outcome == 1 by top10_column1

df_bar <- dataframe %>%
  {
    if (use_general_top_group) {
      # Use precomputed flag (assumed logical TRUE/FALSE)
      dplyr::mutate(., top = .data[[predef_top]])
    } else {
      # Recompute top using the subset 'dataframe'
      thr <- stats::quantile(.[[GRS]], probs = bottom_cutoff, na.rm = TRUE, names = FALSE)
      dplyr::mutate(., top = .data[[GRS]] >= thr)
    }
  } %>%
  # Optional: enforce logical type if your predefined flag is 0/1
  # mutate(top = as.logical(top)) %>%
  filter(!is.na(top)) %>%                           # drop NA flags
  group_by(top) %>%
  summarise(
    # robust to numeric 0/1 or factor/character "0"/"1"
    n1  = sum(as.character(.data[[outcome]]) == "1", na.rm = TRUE),
    n   = sum(!is.na(.data[[outcome]])),
    pct = dplyr::if_else(n > 0, n1 / n, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(scales::percent(pct, accuracy = 0.1),
                   " (", scales::comma(n1), "/", scales::comma(n), ")")
  )



# 2) Plot 
ggplot(df_bar, aes(x = factor(top, levels = c(FALSE, TRUE)),
                   y = pct)) +
  geom_col(width = 0.6, fill = ex_bright_green) +
  geom_text(aes(label = label), vjust = -0.4, size = 3.5) +
  scale_y_continuous(limits = c(0, 0.25)) +
  scale_x_discrete(labels = c(`FALSE` = paste0("Bottom ", 100-top_segment, "%"), `TRUE` = paste0("Top ", top_segment, "%"))) +
  labs(
    x = xlabel,
    y = ylabel,
    title = title_choice
  ) +
  scale_fill_manual(values=exeter_full_palette)+scale_colour_manual(values=exeter_full_palette)+
  ex_green_theme +
  theme_minimal()

