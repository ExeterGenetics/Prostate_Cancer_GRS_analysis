##################################################################################################
#---------------- Turning Our Future Health summary statistics tables into plots ----------------#
##################################################################################################

### This script is intended to be run AFTER creating the summary statistics tables in "Logistic_Regressions.ipynb"
### and "Age_at_Diagnosis.ipynb" in the Our Future Health TRE on DNAnexus. After retrieving those files from the
### TRE via an airlock request, this script creates plots from the summary statistics tables.

###############################################################################################
# Step 0: set working directory to the folder where the summary statistics tables are stored, #
#         and read in the summary statistics tables                                           #
###############################################################################################

setwd("C:/Users/cw1245/Documents/PhD/LocalCoding/OFH_Airlock_Downloads") 
library(dplyr)
library(ggplot2)

# Read in summary statistics table from "Logistic_Regressions.ipynb"
bulk_results_filtered <- read.csv("bulk_results_filtered.csv", header = TRUE) 

# Read in summary statistics table from "Age_at_Diagnosis.ipynb"
bulk_GRS_vs_age_at_diagnosis <- read.csv("bulk_GRS_vs_age_at_diagnosis.csv", header = TRUE)

# Set colour palette for plots
study_palette <- c(
  "#00c896", "#6ab3e7", "#93272c", "#022020", "#702081", "#007d69", "#ffc62c", "#e78699", "#003c3c",  
  "#00a87e", "#b46a55", "#e60000", "#f9423a",  "#9569be","#fc4c02", "#ff7f41", "#f3d54e", "#250e62", 
  "#00dca5", "#f4c3cd", "#898b8d"
)

########################################################################
#========= Step 1: Create plots for bulk logistic regressions =========#
########################################################################

create_forest_plot <- function(data = bulk_results_filtered,
                               outcome = NULL,
                               screening_date = NULL,
                               populations = c("White", "Black", "EUR", "AFR"),
                               predictors = c("age_asd", "ContiAFRGRS", "ContiEURGRS", "ContiMultiethnicGRS", "WangAFRGRS", "WangEURGRS", "WangMultiethnicGRS"),
                               covariates = c("None", "age_asd", "FHx_PrCa_BrCa", "age_asd + FHx_PrCa_BrCa"),
                               metric = NULL
                               ) {

## Warning messages to guide users (note: do not stop, just print warning messages)

if (!is.null(outcome) & !is.null(screening_date) & !is.null(metric)) {
  print(paste0("Creating forest plot for outcome '", outcome, "' at screening date '", screening_date, "' using metric '", metric, "'"))
} else {
  print("Warning: Please specify an outcome, screening date, and metric to create the forest plot.")
}

if (!outcome %in% unique(data$Outcome)) {
  print(paste0("Warning: Outcome '", outcome, "' not found in the data. Please check the available outcomes. These are typically: 'PrCa_post_assessment', 'PrCa_2yrs', 'PrCa_5yrs', 'PrCa_10yrs'"))
}

if (!screening_date %in% unique(data$Screening_Date)) {
  print(paste0("Warning: Screening date '", screening_date, "' not found in the data. Please check the available screening dates. These are typically: 'artificial_screening_date', 'date40', 'date50', 'date55', 'date60', 'date65'"))
}

if (!all(populations %in% unique(data$Population))) {
  unindentified_populations <- populations[!populations %in% unique(data$Population)]
  print(paste0("Warning: One or more populations not found in the data: ", paste(unindentified_populations, collapse = ", "), ". Please check the available populations. These are typically: 'White', 'Black', 'EUR', 'AFR'"))
}

if (!all(predictors %in% unique(data$Predictor))) {
  print(paste0("Warning: One or more predictors not found in the data. Please check the available predictors. These are typically: 'age_asd', 'ContiAFRGRS', 'ContiEURGRS', 'ContiMultiethnicGRS', 'WangAFRGRS', 'WangEURGRS', 'WangMultiethnicGRS'"))
}

if (!all(covariates %in% unique(data$Covariates))) {
  print(paste0("Warning: One or more covariates not found in the data. Please check the available covariates. These are typically: 'None', 'age_asd', 'FHx_PrCa_BrCa', 'age_asd + FHx_PrCa_BrCa'"))
}

if (!metric %in% c("ROC_AUC", "OR_per_1SD_global", "OR_per_1SD_subpopulation")) {
  print(paste0("Warning: Metric '", metric, "' not found in the data. Please check the available metrics. These are typically: 'ROC_AUC', 'OR_per_1SD_global', 'OR_per_1SD_subpopulation'"))
}

forestplot <- data %>%
  dplyr::filter(
    Population %in% populations,
    Predictor %in% predictors,
    Covariates %in% covariates,
    Outcome == outcome,
    Screening_Date == screening_date
  ) %>%
  dplyr::filter(
    !is.na(.data[[metric]]),
    !is.na(.data[[paste0(metric, "_CI_Lower")]]),
    !is.na(.data[[paste0(metric, "_CI_Upper")]])
  ) %>%
  dplyr::mutate(
    Population = factor(Population, levels = rev(populations)),
    Model = paste(Predictor, Covariates, sep = " | ")
  )

plot_limits <- if (metric == "ROC_AUC") {
    c(0.5, 1.0)
  } else if (metric %in% c("OR_per_1SD_global", "OR_per_1SD_subpopulation")) {
    c(1, 8)
  } else {
    NULL
  }

outcome_label <- dplyr::case_when(
  outcome == "PrCa_2yrs" ~ "2-year",
  outcome == "PrCa_5yrs" ~ "5-year",
  outcome == "PrCa_10yrs" ~ "10-year",
  outcome == "PrCa_post_assessment" ~ "Post-assessment",
  TRUE ~ outcome
)

metric_label <- dplyr::case_when(
  metric == "ROC_AUC" ~ "ROC AUC",
  metric == "OR_per_1SD_global" ~ "OR per 1SD (Global)",
  metric == "OR_per_1SD_subpopulation" ~ "OR per 1SD (Subpopulation)",
  TRUE ~ metric
)

metric_axis_label <- paste0(outcome_label, " ", metric_label, " [95% CI]")

ggplot(
  forestplot,
  aes(x = Population, y = .data[[metric]], color = Model)
) +
  geom_hline(yintercept = ifelse(metric == "ROC_AUC", 0.5, 1), linetype = "dashed", color = "#898b8d") +
  geom_errorbar(
    aes(ymin = .data[[paste0(metric, "_CI_Lower")]], ymax = .data[[paste0(metric, "_CI_Upper")]]),
    position = position_dodge(width = 0.6),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 2.5
  ) +
  coord_flip(ylim = plot_limits) +
  labs(
    title = paste0("Forest plot for outcome '", outcome, "' at screening date '", screening_date, "'"),
    x = "Population",
    y = metric_axis_label,
    color = "Model"
  ) +
  scale_color_manual(values = study_palette[1:length(unique(forestplot$Model))]) +
  guides(color = guide_legend(reverse = TRUE)) +
  theme_minimal()

}


## Plotting time

### Screening date 40

#### No point in 2 year prediction horizon, barely any models for screening date 40

#### No point in 5 year prediction horizon, because only White/EUR men included

#### 10 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date40", outcome = "PrCa_10yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date40", outcome = "PrCa_10yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date40", outcome = "PrCa_10yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

### Screening date 50

#### No point in 2 year prediction horizon, because only White/EUR men included

#### 5 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date50", outcome = "PrCa_5yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date50", outcome = "PrCa_5yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date50", outcome = "PrCa_5yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

#### 10 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date50", outcome = "PrCa_10yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date50", outcome = "PrCa_10yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date50", outcome = "PrCa_10yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

### Screening date 55

#### 2 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date55", outcome = "PrCa_2yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date55", outcome = "PrCa_2yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date55", outcome = "PrCa_2yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

#### 5 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date55", outcome = "PrCa_5yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date55", outcome = "PrCa_5yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date55", outcome = "PrCa_5yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

#### 10 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date55", outcome = "PrCa_10yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date55", outcome = "PrCa_10yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date55", outcome = "PrCa_10yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

### Screening date 60

#### 2 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date60", outcome = "PrCa_2yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date60", outcome = "PrCa_2yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date60", outcome = "PrCa_2yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

#### 5 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date60", outcome = "PrCa_5yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date60", outcome = "PrCa_5yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date60", outcome = "PrCa_5yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

#### 10 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date60", outcome = "PrCa_10yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date60", outcome = "PrCa_10yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date60", outcome = "PrCa_10yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

### Screening date 65

#### 2 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date65", outcome = "PrCa_2yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date65", outcome = "PrCa_2yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date65", outcome = "PrCa_2yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

#### 5 year prediction horizon
create_forest_plot(data = bulk_results_filtered, screening_date = "date65", outcome = "PrCa_5yrs", metric = "ROC_AUC", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date65", outcome = "PrCa_5yrs", metric = "OR_per_1SD_global", populations = c("White", "Black"))
create_forest_plot(data = bulk_results_filtered, screening_date = "date65", outcome = "PrCa_5yrs", metric = "OR_per_1SD_subpopulation", populations = c("White", "Black"))

#### 10 year prediction horizon is removed for date65, so not included here

##############################################################################
# Step 2 - Create plots for linear regressions (predicting Age at Diagnosis) #
##############################################################################

## Since Our Future Health does not like scatter plots being exported, as the points represent individual participants,
## we will plot the trendline of the linear regression instead

library(tidyr)

create_trendline_plot <- function(data = bulk_GRS_vs_age_at_diagnosis,
                                  populations = c("White", "Black", "EUR", "AFR"),
                                  predictors = c("ContiAFRGRS", "ContiEURGRS", "ContiMultiethnicGRS",
                                                 "WangAFRGRS", "WangEURGRS", "WangMultiethnicGRS"),
                                  normalisation = c("Subpopulation", "Global"),
                                  covariates = NULL,
                                  grs_range = c(-3, 3),
                                  grs_step = 0.1) {
  
  plot_data <- data %>%
    dplyr::filter(
      Population %in% populations,
      Predictor %in% predictors,
      Normalisation %in% normalisation
    )
  
  if (!is.null(covariates)) {
    plot_data <- plot_data %>%
      dplyr::filter(Covariates %in% covariates)
  }
  
  plot_data <- plot_data %>%
    dplyr::mutate(
      Model = paste(Predictor, Normalisation, sep = " | "),
      GRS_CI_Lower = GRS.coefficient - 1.96 * GRS.SE,
      GRS_CI_Upper = GRS.coefficient + 1.96 * GRS.SE
    )
  
  prediction_grid <- tidyr::crossing(
    plot_data,
    GRS = seq(grs_range[1], grs_range[2], by = grs_step)
  ) %>%
    dplyr::mutate(
      Predicted_Age_at_Diagnosis = Intercept + GRS.coefficient * GRS,
      Predicted_Age_CI_1 = Intercept + GRS_CI_Lower * GRS,
      Predicted_Age_CI_2 = Intercept + GRS_CI_Upper * GRS,
      Predicted_Age_Lower = pmin(Predicted_Age_CI_1, Predicted_Age_CI_2),
      Predicted_Age_Upper = pmax(Predicted_Age_CI_1, Predicted_Age_CI_2)
    )
  
  ggplot(
      prediction_grid,
      aes(
        x = GRS,
        y = Predicted_Age_at_Diagnosis,
        color = Model,
        fill = Model
      )
    ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#39393a", linewidth = 0.6) +
      geom_ribbon(
        aes(
          ymin = Predicted_Age_Lower,
          ymax = Predicted_Age_Upper
        ),
        alpha = 0.15,
        color = NA
      ) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ Population) +
      labs(
        title = paste("Predicted age at diagnosis by genetic risk score in", paste(populations, collapse = ", "), "men"),
        x = "Genetic Risk Score",
        y = "Predicted age at diagnosis",
        color = "Model",
        fill = "Model"
      ) +
      scale_color_manual(values = study_palette[1:length(unique(prediction_grid$Model))]) +
      scale_fill_manual(values = study_palette[1:length(unique(prediction_grid$Model))]) +
      theme_minimal()
}

## Plotting time

create_trendline_plot(data = bulk_GRS_vs_age_at_diagnosis, populations = c("Black"), predictors = c("ContiMultiethnicGRS", "WangMultiethnicGRS"), normalisation = c("Global"), covariates = NULL, grs_range = c(-1, 4), grs_step = 0.1)
create_trendline_plot(data = bulk_GRS_vs_age_at_diagnosis, populations = c("White"), predictors = c("ContiMultiethnicGRS", "WangMultiethnicGRS"), normalisation = c("Global"), covariates = NULL, grs_range = c(-1, 4), grs_step = 0.1)

create_trendline_plot(data = bulk_GRS_vs_age_at_diagnosis, populations = c("Black"), predictors = c("ContiMultiethnicGRS", "WangMultiethnicGRS"), normalisation = c("Subpopulation"), covariates = NULL, grs_range = c(-1, 4), grs_step = 0.1)
create_trendline_plot(data = bulk_GRS_vs_age_at_diagnosis, populations = c("White"), predictors = c("ContiMultiethnicGRS", "WangMultiethnicGRS"), normalisation = c("Subpopulation"), covariates = NULL, grs_range = c(-1, 4), grs_step = 0.1)


###################################################################################
# Step 3 - plot example net benefit curve to contextualise "useful range" results #
###################################################################################

##############################################
# Example net benefit curve for illustration #
##############################################

plot_example_nb_curve <- function(
    cutoff_values = seq(0, 0.80, by = 0.01),
    prevalence = 0.2,
    plot_title = "Example Decision Curve Analysis"
) {
  
  # Threshold probabilities
  pt <- cutoff_values
  
  # Treat none: net benefit is always 0
  net_benefit_treat_none <- rep(0, length(pt))
  
  # Treat all:
  # NB = prevalence - (1 - prevalence) * pt / (1 - pt)
  net_benefit_treat_all <- prevalence - (1 - prevalence) * (pt / (1 - pt))
  
  # Example model net benefit.
  # This is manually constructed so that:
  # - the model is below treat all at low cutoffs
  # - the model becomes better than both treat all and treat none
  # - the model gradually declines
  # - the model eventually crosses below treat none
  model_control_points <- data.frame(
    cutoff = c(0.00, 0.10, 0.20, 0.30, 0.45, 0.60, 0.70, 0.80),
    nb     = c(0.16, 0.12, 0.055, 0.045, 0.032, 0.015, 0.002, -0.012)
  )
  
  net_benefit_model <- approx(
    x = model_control_points$cutoff,
    y = model_control_points$nb,
    xout = pt,
    rule = 2
  )$y
  
  # Combine into a dataframe, mirroring your nb_curve() output style
  example_nb_results <- data.frame(
    pt = pt,
    net_benefit_model = net_benefit_model,
    net_benefit_treat_all = net_benefit_treat_all,
    net_benefit_treat_none = net_benefit_treat_none
  )
  
  # Plot in the same colour/style as your nb_curve() function
  plot(
    example_nb_results$pt,
    example_nb_results$net_benefit_model,
    type = "l",
    col = "blue",
    lwd = 2,
    xlab = "Threshold Probability",
    ylab = "Net Benefit",
    main = plot_title,
    ylim = c(-0.05, max(prevalence, net_benefit_model, na.rm = TRUE) + 0.03)
  )
  
  lines(
    example_nb_results$pt,
    example_nb_results$net_benefit_treat_all,
    col = "red",
    lwd = 2
  )
  
  lines(
    example_nb_results$pt,
    example_nb_results$net_benefit_treat_none,
    col = "green",
    lwd = 2
  )
  
  legend(
    "topright",
    legend = c("Model", "Treat All", "Treat None"),
    col = c("blue", "red", "green"),
    lwd = 2
  )
  
  return(example_nb_results)
}

plot_example_nb_curve()
