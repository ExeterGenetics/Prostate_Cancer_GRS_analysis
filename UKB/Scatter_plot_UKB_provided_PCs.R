####### Scatter plot for visualising UKB Principal Components #########
# Requires first running "Logistic_Regressions_Conti_GRS.R" to create PCa_iv_covariates_GRS_clean

library(ggplot2)

test_population <- PCa_iv_covariates_GRS_clean %>%
dplyr::filter(
#ethnicity_group_wide == "White" | ethnicity_group_narrow == "Black" |
Ethnicity == "White and Black African" | Ethnicity == "White and Black Caribbean"
)

ggplot(test_population, aes(x = p22009_a1, y = p22009_a2, color = ethnicity_group_narrow)) +
  geom_point(size = 2) +
  scale_color_discrete(drop = TRUE) +
  scale_y_continuous(limits = c(-300, 100), expand = expansion(mult = 0)) +
  scale_color_discrete(drop = TRUE) +
  theme_minimal() +
  labs(
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Self-reported Ethnicity"
  )

test_population %>%
  filter(Genomic_ancestry == "European ancestry (EUR)") %>%
  summarise(mean_value = mean(p22009_a1, na.rm = TRUE))
