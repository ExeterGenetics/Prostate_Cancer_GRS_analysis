# This script sets up the logistic regression, confusion matrix, and net benefit decision curve functions used in 
# Logistic_Regressions.ipynb

## Helper: attempt installation without stopping immediately if installation fails

try_install_package <- function(pkg) {
  warnings_seen <- character()

  ok <- tryCatch(
    {
      withCallingHandlers(
        {
          install.packages(pkg)
        },
        warning = function(w) {
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      requireNamespace(pkg, quietly = TRUE)
    },
    error = function(e) {
      message(
        "Could not install package '", pkg, "': ",
        conditionMessage(e)
      )
      FALSE
    }
  )

  if (length(warnings_seen) > 0) {
    message(
      "Warnings while installing package '", pkg, "':\n",
      paste("-", warnings_seen, collapse = "\n")
    )
  }

  ok
}

## Required packages

required_packages <- c(
  "dplyr",
  "lubridate"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Required package '", pkg, "' is not available. Attempting installation...")
    try_install_package(pkg)
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required package '", pkg, "' is not available and could not be installed. ",
      "Please request that it is added to the TRE environment."
    )
  }

  library(pkg, character.only = TRUE)
  message("Required package '", pkg, "' loaded successfully.")
}

## Optional packages

optional_packages <- c(
  "pROC",
  "PRROC",
  "epiR"
)

for (pkg in optional_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Optional package '", pkg, "' is not available. Attempting installation...")
    try_install_package(pkg)
  }

  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    message("Optional package '", pkg, "' loaded successfully.")
  } else {
    message(
      "Optional package '", pkg, "' is not available. ",
      "Using fallback behaviour where possible."
    )
  }
}

## Subset controls to match the number of cases

subset_controls_to_case_count <- function(df, case_col) {
  cases    <- df[!is.na(df[[case_col]]) & df[[case_col]] == 1L, , drop = FALSE]
  controls <- df[!is.na(df[[case_col]]) & df[[case_col]] == 0L, , drop = FALSE]

  ncases    <- nrow(cases)
  ncontrols <- nrow(controls)

  if (ncases == 0L || ncontrols == 0L) {
    return(df)
  }

  if (ncontrols < ncases) {
    warning("Fewer controls than cases; returning all controls and all cases.")
  }

  n_to_sample     <- min(ncases, ncontrols)
  sampled_controls <- controls[sample.int(ncontrols, size = n_to_sample, replace = FALSE), , drop = FALSE]

  dplyr::bind_rows(cases, sampled_controls)
}

## Skip models that use age_asd at fixed screening dates (date40, date50 etc.)

age_asd_skip_message <- "Age is not an appropriate covariate for this screening date. Skipping this model."

fixed_age_screening_dates <- c("date40",
                               "date45",
                               "date50",
                               "date55",
                               "date60",
                               "date65"
                          ) 

## Skip models that use age_asd at fixed screening ages

age_asd_skip_message <- "Age is not an appropriate covariate for this screening date. Skipping model"

fixed_age_screening_dates <- c(
  "date40",
  "date45",
  "date50",
  "date55",
  "date60",
  "date65"
)

uses_age_asd_inappropriately <- function(screening_date,
                                         predictor,
                                         covariates = NULL) {
  !is.null(screening_date) &&
    length(screening_date) == 1L &&
    !is.na(screening_date) &&
    screening_date %in% fixed_age_screening_dates &&
    (
      identical(predictor, "age_asd") ||
        "age_asd" %in% covariates
    )
}

## Helper function to calculate calibration metrics and optionally produce a calibration plot

assess_logistic_calibration <- function(
    data,
    outcome,
    pred_col = "pred",
    n_groups = 10L,
    plot_calibration = FALSE,
    plot_title = "Calibration plot",
    show_output = TRUE
) {

  ## Validate inputs

  if (!outcome %in% names(data)) {
    stop(
      "Outcome column not found in calibration data: ",
      outcome
    )
  }

  if (!pred_col %in% names(data)) {
    stop(
      "Prediction column not found in calibration data: ",
      pred_col
    )
  }

  if (
    length(n_groups) != 1L ||
    is.na(n_groups) ||
    n_groups < 2L
  ) {
    stop("n_groups must be a single integer of at least 2.")
  }

  n_groups <- as.integer(n_groups)

  ## Extract and safely convert the outcome

  observed_raw <- data[[outcome]]

  if (is.logical(observed_raw)) {

    observed <- as.integer(observed_raw)

  } else if (
    is.factor(observed_raw) ||
    is.character(observed_raw)
  ) {

    observed <- suppressWarnings(
      as.numeric(as.character(observed_raw))
    )

  } else {

    observed <- as.numeric(observed_raw)
  }

  predicted <- as.numeric(data[[pred_col]])

  ## Retain valid observations only

  valid_rows <- (
    !is.na(observed) &
      !is.na(predicted) &
      is.finite(observed) &
      is.finite(predicted)
  )

  calibration_individual <- data.frame(
    observed = observed[valid_rows],
    pred = predicted[valid_rows],
    stringsAsFactors = FALSE
  )

  if (nrow(calibration_individual) == 0L) {
    stop(
      "No complete observations are available for calibration."
    )
  }

  if (!all(calibration_individual$observed %in% c(0, 1))) {
    stop(
      "The outcome must be coded numerically as 0 and 1."
    )
  }

  if (
    any(
      calibration_individual$pred < 0 |
        calibration_individual$pred > 1
    )
  ) {
    stop(
      "All predicted probabilities must be between 0 and 1."
    )
  }

  ## Prevent infinite logits for predictions exactly equal to 0 or 1

  calibration_individual$pred_clipped <- pmin(
    pmax(calibration_individual$pred, 1e-15),
    1 - 1e-15
  )

  calibration_individual$logit_pred <- qlogis(
    calibration_individual$pred_clipped
  )

  ## Basic calibration metrics

  mean_predicted_risk <- mean(
    calibration_individual$pred
  )

  observed_event_rate <- mean(
    calibration_individual$observed
  )

  brier_score <- mean(
    (
      calibration_individual$observed -
        calibration_individual$pred
    )^2
  )

  n_cases <- sum(
    calibration_individual$observed == 1
  )

  n_controls <- sum(
    calibration_individual$observed == 0
  )

  calibration_in_the_large <- NA_real_
  calibration_slope <- NA_real_

  citl_model <- NULL
  calibration_slope_model <- NULL

  ## Calibration regression requires both outcome classes

  if (n_cases > 0L && n_controls > 0L) {

    citl_model <- tryCatch(
      glm(
        observed ~ 1 + offset(logit_pred),
        data = calibration_individual,
        family = binomial()
      ),
      error = function(e) {
        warning(
          "Calibration-in-the-large model failed: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(citl_model)) {
      calibration_in_the_large <- unname(
        coef(citl_model)[1]
      )
    }

    calibration_slope_model <- tryCatch(
      glm(
        observed ~ logit_pred,
        data = calibration_individual,
        family = binomial()
      ),
      error = function(e) {
        warning(
          "Calibration-slope model failed: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(calibration_slope_model)) {
      calibration_slope <- unname(
        coef(calibration_slope_model)["logit_pred"]
      )
    }

  } else {

    warning(
      "Calibration intercept and slope could not be calculated ",
      "because the calibration data do not contain both cases ",
      "and controls."
    )
  }

  ## Store metrics

  calibration_metrics <- data.frame(
    N = nrow(calibration_individual),
    N_Cases = n_cases,
    N_Controls = n_controls,
    Mean_predicted_risk = mean_predicted_risk,
    Observed_event_rate = observed_event_rate,
    Calibration_in_the_large =
      calibration_in_the_large,
    Calibration_slope = calibration_slope,
    Brier_score = brier_score,
    stringsAsFactors = FALSE
  )

  ## Create quantile-grouped calibration points

  number_of_groups <- min(
    n_groups,
    nrow(calibration_individual)
  )

  calibration_data <- calibration_individual %>%
    dplyr::mutate(
      calibration_group = dplyr::ntile(
        pred,
        number_of_groups
      )
    ) %>%
    dplyr::group_by(calibration_group) %>%
    dplyr::summarise(
      minimum_predicted = min(pred),
      maximum_predicted = max(pred),
      mean_predicted = mean(pred),
      observed_proportion = mean(observed),
      n = dplyr::n(),
      events = sum(observed),
      .groups = "drop"
    )

  ## Format statistics for the plot

  calibration_subtitle <- paste0(
    "Mean predicted risk = ",
    sprintf("%.3f", mean_predicted_risk),
    "    |    Observed event rate = ",
    sprintf("%.3f", observed_event_rate),
    "\nCalibration-in-the-large = ",
    sprintf("%.3f", calibration_in_the_large),
    "    |    Calibration slope = ",
    sprintf("%.3f", calibration_slope),
    "    |    Brier score = ",
    sprintf("%.3f", brier_score),
    "\nBlue: LOESS calibration curve; black points: ",
    number_of_groups,
    " risk groups"
  )

  ## Optionally construct and print the plot

  calibration_plot <- NULL

  if (isTRUE(plot_calibration)) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning(
        "Calibration plot was not produced because ",
        "ggplot2 is unavailable."
      )
    } else {

      calibration_plot <- ggplot2::ggplot() +
        ggplot2::geom_abline(
          intercept = 0,
          slope = 1,
          linetype = "dashed",
          color = "red"
        ) +
        ggplot2::geom_smooth(
          data = calibration_individual,
          ggplot2::aes(
            x = pred,
            y = observed
          ),
          method = "loess",
          se = FALSE,
          color = "blue"
        ) +
        ggplot2::geom_point(
          data = calibration_data,
          ggplot2::aes(
            x = mean_predicted,
            y = observed_proportion
          ),
          size = 2.5,
          color = "black"
        ) +
        ggplot2::geom_rug(
          data = calibration_individual,
          ggplot2::aes(x = pred),
          sides = "b",
          alpha = 0.10
        ) +
        ggplot2::coord_cartesian(
          xlim = c(0, 1),
          ylim = c(0, 1)
        ) +
        ggplot2::labs(
          title = plot_title,
          subtitle = calibration_subtitle,
          x = "Predicted probability",
          y = "Observed proportion diagnosed"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold",
            size = 14
          ),
          plot.subtitle = ggplot2::element_text(
            size = 10,
            lineheight = 1.15,
            margin = ggplot2::margin(b = 10)
          )
        )

      print(calibration_plot)
    }
  }

  if (isTRUE(show_output)) {
    cat("\nCalibration metrics:\n")
    print(calibration_metrics, row.names = FALSE)
  }

  return(list(
    metrics = calibration_metrics,
    grouped_data = calibration_data,
    individual_data = calibration_individual,
    plot = calibration_plot,
    citl_model = citl_model,
    slope_model = calibration_slope_model
  ))
}

## Run logistic regressions + ROC curve (produces ROC curve from same population as the logistic regression model)

run_logreg <- function(data,
                       outcome,
                       predictor,
                       covariates = NULL,
                       screening_date = NULL,
                       subset_controls = FALSE,
                       plot_roc = TRUE,
                       plot_pr = FALSE,
                       plot_calibration = FALSE,
                       calibration_groups = 10L,
                       show_output = TRUE,
                       remove_dead_controls = TRUE,
                       remove_too_young_controls = TRUE,
                       suppress_warnings = FALSE) {

  manual_roc_fallback <- function(data, outcome, show_output = TRUE) {
    y <- data[[outcome]]
    x <- data$pred

    ok <- complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]

    if (is.factor(y)) {
      y <- as.numeric(y == levels(y)[2])
    }

    if (!all(y %in% c(0, 1))) {
      stop("Outcome must be binary (0/1 or factor with 2 levels) for manual ROC fallback.")
    }

    n_cases <- sum(y == 1)
    n_controls <- sum(y == 0)

    if (n_cases == 0 || n_controls == 0) {
      stop("Manual ROC fallback requires at least one case and one control.")
    }

    ord <- order(x, decreasing = TRUE)
    x <- x[ord]
    y <- y[ord]

    tp <- cumsum(y)
    fp <- cumsum(1 - y)

    P <- sum(y)
    N <- length(y) - P

    tpr <- c(0, tp / P, 1)
    fpr <- c(0, fp / N, 1)

    auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)

    cases_scores <- x[y == 1]
    controls_scores <- x[y == 0]

    V10 <- sapply(cases_scores, function(xi) mean(ifelse(xi > controls_scores, 1, ifelse(xi == controls_scores, 0.5, 0))))
    V01 <- sapply(controls_scores, function(yj) mean(ifelse(cases_scores > yj, 1, ifelse(cases_scores == yj, 0.5, 0))))

    var_auc <- var(V10) / n_cases + var(V01) / n_controls
    se_auc <- sqrt(var_auc)

    ci_lower <- auc - 1.96 * se_auc
    ci_upper <- auc + 1.96 * se_auc

    if (show_output) {
      cat(sprintf("Cases: %d | Controls: %d\n", n_cases, n_controls))
      cat(sprintf("AUC: %.4f (95%% CI: %.4f - %.4f)\n", auc, ci_lower, ci_upper))

      plot(fpr, tpr, type = "l", lwd = 2,
           xlab = "False Positive Rate (1 - Specificity)",
           ylab = "True Positive Rate (Sensitivity)",
           main = sprintf("ROC Curve\nAUC = %.4f (95%% CI: %.4f - %.4f)\nCases = %d  |  Controls = %d",
                          auc, ci_lower, ci_upper, n_cases, n_controls))
      abline(0, 1, lty = 2, col = "grey")
    }

    list(
      method = "manual",
      auc = auc,
      ci = c(lower = ci_lower, upper = ci_upper),
      fpr = fpr,
      tpr = tpr,
      n_cases = n_cases,
      n_controls = n_controls
    )
  }

  ## Skip age_asd models for fixed screening dates

  if (uses_age_asd_inappropriately(
    screening_date = screening_date,
    predictor = predictor,
    covariates = covariates
  )) {
    if (!suppress_warnings) message(age_asd_skip_message)

    return(list(
      skipped = TRUE,
      skip_reason = age_asd_skip_message,
      model = NULL,
      data = NULL,
      outcome = outcome,
      predictor = predictor,
      covariates = covariates,
      roc = NULL,
      pr = NULL,
      calibration = NULL,
      converged = NA,
      subset_controls = subset_controls
    ))
  }

  if ((isTRUE(remove_dead_controls) || isTRUE(remove_too_young_controls)) && is.null(screening_date)) {
    stop("screening_date must be supplied when remove_dead_controls or remove_too_young_controls is TRUE.")
  }


  if (isTRUE(remove_dead_controls)) {
      if ((screening_date == "artificial_screening_date") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_asd == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "artificial_screening_date" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_asd_plus_2yrs == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "artificial_screening_date" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_asd_plus_5yrs == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "artificial_screening_date" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_asd_plus_10yrs == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date40") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_40 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_42 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_45 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_50 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date45") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_45 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_47 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_50 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_55 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date50") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_50 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_52 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_55 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_60 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date55") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_55 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_57 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_60 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_65 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date60") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_60 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_62 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_65 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_70 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date65") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_65 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_67 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_70 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_75 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else {
          if (!suppress_warnings) warning("Outcome not recognised for remove_dead; no participants removed based on death within prediction horizon.")
      }
  }

  if (isTRUE(remove_too_young_controls)) {
      if ((screening_date == "artificial_screening_date")) {
          data <- data %>% dplyr::filter(over_20_at_asd == TRUE)
      } else if((screening_date == "date40") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_40 == TRUE)
      } else if(screening_date == "date40" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_42 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_45 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_50 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date45") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_45 == TRUE)
      } else if(screening_date == "date45" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_47 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_50 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_55 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date50") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_50 == TRUE)
      } else if(screening_date == "date50" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_52 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_55 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_60 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date55") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_55 == TRUE)
      } else if(screening_date == "date55" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_57 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_60 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_65 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date60") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_60 == TRUE)
      } else if(screening_date == "date60" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_62 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_65 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_70 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date65") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_65 == TRUE)
      } else if(screening_date == "date65" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_67 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_70 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_75 == TRUE | .data[[outcome]] == 1)
      } else {
          if (!suppress_warnings) warning("Outcome not recognised for remove_too_young_controls; no control participants removed based on age within prediction horizon.")
      }
  }

    if (isTRUE(subset_controls)) {
      data <- subset_controls_to_case_count(data, case_col = outcome)
    }

  if (is.null(covariates)) {
    formula <- as.formula(paste(outcome, "~", predictor))
  } else {
    formula <- as.formula(
      paste(outcome, "~", predictor, "+", paste(covariates, collapse = " + "))
    )
  }
  
  ## Complete-case check

  # Split interaction/transformed terms (e.g. "A*B") into underlying column names
  covariate_vars <- unique(unlist(lapply(covariates, function(x) all.vars(stats::as.formula(paste("~", x))))))
  model_vars <- unique(c(outcome, all.vars(stats::as.formula(paste("~", predictor))), covariate_vars))

  missing_model_vars <- setdiff(model_vars, names(data))

  if (length(missing_model_vars) > 0) {
    if (!suppress_warnings) warning(
      "Missing model variable(s): ",
      paste(missing_model_vars, collapse = ", ")
    )
  }

  model_data_ok <- stats::complete.cases(data[, model_vars, drop = FALSE])

  model_data <- data[model_data_ok, , drop = FALSE]

  n_cases_model <- sum(model_data[[outcome]] == 1, na.rm = TRUE)
  n_controls_model <- sum(model_data[[outcome]] == 0, na.rm = TRUE)

  if (n_cases_model == 0 || n_controls_model == 0) {
    stop("Model requires at least one case and one control among complete cases.")
  }

  if (n_cases_model < 10) {
    if (!suppress_warnings) warning("Fewer than 10 complete-case cases; model estimates may be unstable.")
  }

  # Guard against predictors/covariates that collapse to one observed value in complete-case data.
  predictor_terms <- setdiff(model_vars, outcome)
  for (term in predictor_terms) {
    x <- model_data[[term]]

    if (is.character(x)) {
      x <- factor(x)
    }

    if (is.factor(x)) {
      x <- droplevels(x)
      if (nlevels(x) < 2) {
        stop(sprintf("Model term '%s' has fewer than 2 observed levels after filtering and complete-case selection. Remove this term for this subgroup/outcome or ensure both levels are present.", term))
      }
    }
  }

  ## Run logistic regression model

  logreg <- glm(formula, data = model_data, family = binomial)

  if (show_output) {
    print(summary(logreg))
  }
  data$pred <- predict(logreg, data, type = "response")
  data$residual <- data[[outcome]] - data$pred
  data$predtop10 <- data$pred >= quantile(data$pred, probs = 0.9, na.rm = TRUE, names = FALSE)
  data$predtop20 <- data$pred >= quantile(data$pred, probs = 0.8, na.rm = TRUE, names = FALSE)
  data$predtop30 <- data$pred >= quantile(data$pred, probs = 0.7, na.rm = TRUE, names = FALSE)
  data$predtop40 <- data$pred >= quantile(data$pred, probs = 0.6, na.rm = TRUE, names = FALSE)
  data$predtop50 <- data$pred >= quantile(data$pred, probs = 0.5, na.rm = TRUE, names = FALSE)
  data$predtop60 <- data$pred >= quantile(data$pred, probs = 0.4, na.rm = TRUE, names = FALSE)
  data$predtop70 <- data$pred >= quantile(data$pred, probs = 0.3, na.rm = TRUE, names = FALSE)
  data$predtop80 <- data$pred >= quantile(data$pred, probs = 0.2, na.rm = TRUE, names = FALSE)
  data$predtop90 <- data$pred >= quantile(data$pred, probs = 0.1, na.rm = TRUE, names = FALSE)

  ## Assess apparent calibration in the model-fitting cohort

  calibration_obj <- tryCatch(
    assess_logistic_calibration(
      data = data,
      outcome = outcome,
      pred_col = "pred",
      n_groups = calibration_groups,
      plot_calibration = plot_calibration && show_output,
      plot_title = paste0(
        "Calibration plot: ",
        predictor,
        if (!is.null(covariates)) {
          paste0(
            " + ",
            paste(covariates, collapse = " + ")
          )
        } else {
          ""
        },
        "\nOutcome: ",
        outcome
      ),
      show_output = show_output
    ),
    error = function(e) {
      if (!suppress_warnings) warning(
        "Calibration assessment failed: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  proc_available <- tryCatch(
    requireNamespace("pROC", quietly = TRUE),
    error = function(e) FALSE
  )

  prroc_available <- tryCatch(
    requireNamespace("PRROC", quietly = TRUE),
    error = function(e) FALSE
  )

  roc_obj <- NULL
  if (proc_available) {
    if (plot_roc && show_output) {
      roc_obj <- pROC::roc(
        response = data[[outcome]],
        predictor = data$pred,
        levels = c(0, 1),
        direction = "<",
        plot = TRUE,
        print.auc = TRUE,
        ci = TRUE
      )

    } else {
      
      roc_obj <- pROC::roc(
        response = data[[outcome]],
        predictor = data$pred,
        levels = c(0, 1),
        direction = "<",
        ci = TRUE
      )

    }
    if (show_output) {
      print(roc_obj)
    }
  } else {
    roc_obj <- tryCatch(
      manual_roc_fallback(data = data, outcome = outcome, show_output = (plot_roc && show_output)),
      error = function(e) {
        if (!suppress_warnings) warning(paste0("Manual ROC fallback failed: ", e$message))
        NULL
      }
    )
    if (is.null(roc_obj)) {
      if (!suppress_warnings) warning("pROC is unavailable and manual ROC fallback failed; ROC metrics set to NA.")
    }
  }
  
  pr_obj <- NULL
  if (plot_pr) {
    labels <- as.numeric(data[[outcome]] == 1)
    scores <- data$pred
    keep <- !is.na(labels) & !is.na(scores)
    labels_clean <- labels[keep]
    scores_clean <- scores[keep]
    prevalence <- mean(labels_clean == 1)

    if (!prroc_available) {
      if (!suppress_warnings) warning("PR curve skipped: PRROC is unavailable or failed to load in this session.")
    } else {
      pr_obj <- tryCatch(
        PRROC::pr.curve(
          scores.class0 = scores_clean[labels_clean == 1],
          scores.class1 = scores_clean[labels_clean == 0],
          curve = TRUE
        ),
        error = function(e) {
          if (!suppress_warnings) warning(paste0("PR curve skipped due to PRROC error: ", e$message))
          NULL
        }
      )

      if (!is.null(pr_obj) && plot_pr && show_output) {
        model_label <- if (is.null(covariates)) {
          predictor
        } else {
          paste(predictor, "+", paste(covariates, collapse = " + "))
        }
        plot(pr_obj, main = paste("Precision-Recall Curve:", model_label),
             xlab = "Recall", ylab = "Precision")
        abline(h = prevalence, col = "red", lty = 2, lwd = 2)
        legend("bottomleft",
               legend = paste0("Outcome prevalence = ", sprintf("%.3f", prevalence)),
               col = "red", lty = 2, lwd = 2, bty = "n")
      }
    }
    
    if (show_output) {
      print(pr_obj)
    }
  }
  
  if (!isTRUE(logreg$converged)) {
    if (!suppress_warnings) warning("Logistic regression did not converge.")
  }
  
  return(list(
    model = logreg,
    data = data,
    outcome = outcome,
    predictor = predictor,
    covariates = covariates,
    roc = roc_obj,
    pr = pr_obj,
    calibration = calibration_obj,
    converged = logreg$converged,
    subset_controls = subset_controls
  ))

}


## Run logistic regression model on one population then use to predict outcomes in another population 
## (same as above but uses model_population and eval_population arguments instead of just "data")

run_logreg_advanced <- function(model_population,
                                eval_population,
                                outcome,
                                predictor,
                                covariates = NULL,
                                screening_date = NULL,
                                subset_controls = FALSE,
                                plot_roc = TRUE,
                                plot_pr = FALSE,
                                plot_calibration = TRUE,
                                calibration_groups = 10L,
                                show_output = TRUE,
                                remove_dead_controls = TRUE,
                                remove_too_young_controls = TRUE,
                                suppress_warnings = FALSE) {

  manual_roc_fallback <- function(data, outcome, show_output = TRUE) {
    y <- data[[outcome]]
    x <- data$pred

    ok <- complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]

    if (is.factor(y)) {
      y <- as.numeric(y == levels(y)[2])
    }

    if (!all(y %in% c(0, 1))) {
      stop("Outcome must be binary (0/1 or factor with 2 levels) for manual ROC fallback.")
    }

    n_cases <- sum(y == 1)
    n_controls <- sum(y == 0)

    if (n_cases == 0 || n_controls == 0) {
      stop("Manual ROC fallback requires at least one case and one control.")
    }

    ord <- order(x, decreasing = TRUE)
    x <- x[ord]
    y <- y[ord]

    tp <- cumsum(y)
    fp <- cumsum(1 - y)

    P <- sum(y)
    N <- length(y) - P

    tpr <- c(0, tp / P, 1)
    fpr <- c(0, fp / N, 1)

    auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)

    cases_scores <- x[y == 1]
    controls_scores <- x[y == 0]

    V10 <- sapply(cases_scores, function(xi) mean(ifelse(xi > controls_scores, 1, ifelse(xi == controls_scores, 0.5, 0))))
    V01 <- sapply(controls_scores, function(yj) mean(ifelse(cases_scores > yj, 1, ifelse(cases_scores == yj, 0.5, 0))))

    var_auc <- var(V10) / n_cases + var(V01) / n_controls
    se_auc <- sqrt(var_auc)

    ci_lower <- auc - 1.96 * se_auc
    ci_upper <- auc + 1.96 * se_auc

    if (show_output) {
      cat(sprintf("Cases: %d | Controls: %d\n", n_cases, n_controls))
      cat(sprintf("AUC: %.4f (95%% CI: %.4f - %.4f)\n", auc, ci_lower, ci_upper))

      plot(fpr, tpr, type = "l", lwd = 2,
           xlab = "False Positive Rate (1 - Specificity)",
           ylab = "True Positive Rate (Sensitivity)",
           main = sprintf("ROC Curve\nAUC = %.4f (95%% CI: %.4f - %.4f)\nCases = %d  |  Controls = %d",
                          auc, ci_lower, ci_upper, n_cases, n_controls))
      abline(0, 1, lty = 2, col = "grey")
    }

    list(
      method = "manual",
      auc = auc,
      ci = c(lower = ci_lower, upper = ci_upper),
      fpr = fpr,
      tpr = tpr,
      n_cases = n_cases,
      n_controls = n_controls
    )
  }

  ## Skip age_asd models for fixed screening dates

  if (uses_age_asd_inappropriately(
    screening_date = screening_date,
    predictor = predictor,
    covariates = covariates
  )) {
    if (!suppress_warnings) message(age_asd_skip_message)

    return(list(
      skipped = TRUE,
      skip_reason = age_asd_skip_message,
      model = NULL,
      model_population = NULL,
      eval_population = NULL,
      outcome = outcome,
      predictor = predictor,
      covariates = covariates,
      roc = NULL,
      pr = NULL,
      calibration = NULL,
      converged = NA,
      subset_controls = subset_controls
    ))
  }

  if ((isTRUE(remove_dead_controls) || isTRUE(remove_too_young_controls)) && is.null(screening_date)) {
    stop("screening_date must be supplied when remove_dead_controls or remove_too_young_controls is TRUE.")
  }

  ## Remove controls who died within the prediction horizon before the screening date

  filter_dead_controls <- function(data, screening_date, outcome) {
      if ((screening_date == "artificial_screening_date") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_asd == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "artificial_screening_date" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_asd_plus_2yrs == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "artificial_screening_date" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_asd_plus_5yrs == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "artificial_screening_date" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_asd_plus_10yrs == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date40") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_40 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_42 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_45 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_50 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date45") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_45 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_47 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_50 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_55 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date50") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_50 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_52 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_55 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_60 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date55") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_55 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_57 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_60 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_65 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date60") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_60 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_62 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_65 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_70 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else if((screening_date == "date65") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter((died_by_65 == FALSE & PrCa != 1) | PrCa == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter((died_by_67 == FALSE & PrCa_2yrs != 1) | PrCa_2yrs == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter((died_by_70 == FALSE & PrCa_5yrs != 1) | PrCa_5yrs == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter((died_by_75 == FALSE & PrCa_10yrs != 1) | PrCa_10yrs == 1)
      } else {
          if (!suppress_warnings) warning("Outcome not recognised for remove_dead; no participants removed based on death within prediction horizon.")
      }
      data
  }

  ## Remove controls who have not yet reached the prediction horizon (too young / insufficient follow-up)

  filter_too_young_controls <- function(data, screening_date, outcome) {
      if ((screening_date == "artificial_screening_date")) {
          data <- data %>% dplyr::filter(over_20_at_asd == TRUE)
      } else if((screening_date == "date40") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_40 == TRUE)
      } else if(screening_date == "date40" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_42 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_45 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date40" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_50 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date45") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_45 == TRUE)
      } else if(screening_date == "date45" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_47 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_50 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date45" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_55 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date50") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_50 == TRUE)
      } else if(screening_date == "date50" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_52 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_55 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date50" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_60 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date55") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_55 == TRUE)
      } else if(screening_date == "date55" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_57 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_60 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date55" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_65 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date60") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_60 == TRUE)
      } else if(screening_date == "date60" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_62 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_65 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date60" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_70 == TRUE | .data[[outcome]] == 1)
      } else if((screening_date == "date65") & (outcome == "PrCa" | outcome == "PrCa_post_assessment")) {
          data <- data %>% dplyr::filter(over_65 == TRUE)
      } else if(screening_date == "date65" & outcome == "PrCa_2yrs") {
          data <- data %>% dplyr::filter(over_67 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_5yrs") {
          data <- data %>% dplyr::filter(over_70 == TRUE | .data[[outcome]] == 1)
      } else if(screening_date == "date65" & outcome == "PrCa_10yrs") {
          data <- data %>% dplyr::filter(over_75 == TRUE | .data[[outcome]] == 1)
      } else {
          if (!suppress_warnings) warning("Outcome not recognised for remove_too_young_controls; no control participants removed based on age within prediction horizon.")
      }
      data
  }

  ## Apply control-removal filters independently to the model-fitting and evaluation cohorts

  if (isTRUE(remove_dead_controls)) {
    model_population <- filter_dead_controls(model_population, screening_date, outcome)
    eval_population <- filter_dead_controls(eval_population, screening_date, outcome)
  }

  if (isTRUE(remove_too_young_controls)) {
    model_population <- filter_too_young_controls(model_population, screening_date, outcome)
    eval_population <- filter_too_young_controls(eval_population, screening_date, outcome)
  }

  ## subset_controls only balances the model-fitting cohort; the evaluation cohort is left untouched

  if (isTRUE(subset_controls)) {
    model_population <- subset_controls_to_case_count(model_population, case_col = outcome)
  }

  if (is.null(covariates)) {
    formula <- as.formula(paste(outcome, "~", predictor))
  } else {
    formula <- as.formula(
      paste(outcome, "~", predictor, "+", paste(covariates, collapse = " + "))
    )
  }
  
  ## Complete-case check (model-fitting cohort)

  # Split interaction/transformed terms (e.g. "A*B") into underlying column names
  covariate_vars <- unique(unlist(lapply(covariates, function(x) all.vars(stats::as.formula(paste("~", x))))))
  model_vars <- unique(c(outcome, all.vars(stats::as.formula(paste("~", predictor))), covariate_vars))

  missing_model_vars <- setdiff(model_vars, names(model_population))

  if (length(missing_model_vars) > 0) {
    if (!suppress_warnings) warning(
      "Missing model variable(s) in model_population: ",
      paste(missing_model_vars, collapse = ", ")
    )
  }

  missing_eval_vars <- setdiff(model_vars, names(eval_population))

  if (length(missing_eval_vars) > 0) {
    if (!suppress_warnings) warning(
      "Missing model variable(s) in eval_population: ",
      paste(missing_eval_vars, collapse = ", ")
    )
  }

  model_data_ok <- stats::complete.cases(model_population[, model_vars, drop = FALSE])

  model_data <- model_population[model_data_ok, , drop = FALSE]

  n_cases_model <- sum(model_data[[outcome]] == 1, na.rm = TRUE)
  n_controls_model <- sum(model_data[[outcome]] == 0, na.rm = TRUE)

  if (n_cases_model == 0 || n_controls_model == 0) {
    if (!suppress_warnings) warning("Model requires at least one case and one control among complete cases.")
  }

  if (n_cases_model < 10) {
    if (!suppress_warnings) warning("Fewer than 10 complete-case cases; model estimates may be unstable.")
  }

  # Guard against predictors/covariates that collapse to one observed value in complete-case data.
  predictor_terms <- setdiff(model_vars, outcome)
  for (term in predictor_terms) {
    x <- model_data[[term]]

    if (is.character(x)) {
      x <- factor(x)
    }

    if (is.factor(x)) {
      x <- droplevels(x)
      if (nlevels(x) < 2) {
        if (!suppress_warnings) warning(sprintf("Model term '%s' has fewer than 2 observed levels after filtering and complete-case selection. Remove this term for this subgroup/outcome or ensure both levels are present.", term))
      }
    }
  }

  ## Fit logistic regression model on the model-fitting cohort

  logreg <- glm(formula, data = model_data, family = binomial)

  if (show_output) {
    print(summary(logreg))
  }

  model_data$residual <- residuals(logreg, type = "response")

  ## Predict outcomes in the (independent) evaluation cohort

  eval_population$pred <- predict(logreg, eval_population, type = "response")
  eval_population$residual <- eval_population[[outcome]] - eval_population$pred
  eval_population$predtop10 <- eval_population$pred >= quantile(eval_population$pred, probs = 0.9, na.rm = TRUE, names = FALSE)
  eval_population$predtop20 <- eval_population$pred >= quantile(eval_population$pred, probs = 0.8, na.rm = TRUE, names = FALSE)
  eval_population$predtop30 <- eval_population$pred >= quantile(eval_population$pred, probs = 0.7, na.rm = TRUE, names = FALSE)
  eval_population$predtop40 <- eval_population$pred >= quantile(eval_population$pred, probs = 0.6, na.rm = TRUE, names = FALSE)
  eval_population$predtop50 <- eval_population$pred >= quantile(eval_population$pred, probs = 0.5, na.rm = TRUE, names = FALSE)
  eval_population$predtop60 <- eval_population$pred >= quantile(eval_population$pred, probs = 0.4, na.rm = TRUE, names = FALSE)
  eval_population$predtop70 <- eval_population$pred >= quantile(eval_population$pred, probs = 0.3, na.rm = TRUE, names = FALSE)
  eval_population$predtop80 <- eval_population$pred >= quantile(eval_population$pred, probs = 0.2, na.rm = TRUE, names = FALSE)
  eval_population$predtop90 <- eval_population$pred >= quantile(eval_population$pred, probs = 0.1, na.rm = TRUE, names = FALSE)

  ## Assess calibration in the independent evaluation cohort

  calibration_obj <- tryCatch(
    assess_logistic_calibration(
      data = eval_population,
      outcome = outcome,
      pred_col = "pred",
      n_groups = calibration_groups,
      plot_calibration = plot_calibration && show_output,
      plot_title = paste0(
        "External calibration plot: ",
        predictor,
        if (!is.null(covariates)) {
          paste0(
            " + ",
            paste(covariates, collapse = " + ")
          )
        } else {
          ""
        },
        "\nOutcome: ",
        outcome
      ),
      show_output = show_output
    ),
    error = function(e) {
      if (!suppress_warnings) warning(
        "Calibration assessment failed in evaluation cohort: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  proc_available <- tryCatch(
    requireNamespace("pROC", quietly = TRUE),
    error = function(e) FALSE
  )

  prroc_available <- tryCatch(
    requireNamespace("PRROC", quietly = TRUE),
    error = function(e) FALSE
  )

  roc_obj <- NULL
  if (proc_available) {
    if (plot_roc && show_output) {
      roc_obj <- pROC::roc(
        response = eval_population[[outcome]],
        predictor = eval_population$pred,
        levels = c(0, 1),
        direction = "<",
        plot = TRUE,
        print.auc = TRUE,
        ci = TRUE
      )

    } else {
      
      roc_obj <- pROC::roc(
        response = eval_population[[outcome]],
        predictor = eval_population$pred,
        levels = c(0, 1),
        direction = "<",
        ci = TRUE
      )

    }
    if (show_output) {
      print(roc_obj)
    }
  } else {
    roc_obj <- tryCatch(
      manual_roc_fallback(data = eval_population, outcome = outcome, show_output = (plot_roc && show_output)),
      error = function(e) {
        if (!suppress_warnings) warning(paste0("Manual ROC fallback failed: ", e$message))
        NULL
      }
    )
    if (is.null(roc_obj)) {
      if (!suppress_warnings) warning("pROC is unavailable and manual ROC fallback failed; ROC metrics set to NA.")
    }
  }
  
  pr_obj <- NULL
  if (plot_pr) {
    labels <- as.numeric(eval_population[[outcome]] == 1)
    scores <- eval_population$pred
    keep <- !is.na(labels) & !is.na(scores)
    labels_clean <- labels[keep]
    scores_clean <- scores[keep]
    prevalence <- mean(labels_clean == 1)

    if (!prroc_available) {
      if (!suppress_warnings) warning("PR curve skipped: PRROC is unavailable or failed to load in this session.")
    } else {
      pr_obj <- tryCatch(
        PRROC::pr.curve(
          scores.class0 = scores_clean[labels_clean == 1],
          scores.class1 = scores_clean[labels_clean == 0],
          curve = TRUE
        ),
        error = function(e) {
          if (!suppress_warnings) warning(paste0("PR curve skipped due to PRROC error: ", e$message))
          NULL
        }
      )

      if (!is.null(pr_obj) && plot_pr && show_output) {
        model_label <- if (is.null(covariates)) {
          predictor
        } else {
          paste(predictor, "+", paste(covariates, collapse = " + "))
        }
        plot(pr_obj, main = paste("Precision-Recall Curve:", model_label),
             xlab = "Recall", ylab = "Precision")
        abline(h = prevalence, col = "red", lty = 2, lwd = 2)
        legend("bottomleft",
               legend = paste0("Outcome prevalence = ", sprintf("%.3f", prevalence)),
               col = "red", lty = 2, lwd = 2, bty = "n")
      }
    }
    
    if (show_output) {
      print(pr_obj)
    }
  }
  
  if (!isTRUE(logreg$converged)) {
    if (!suppress_warnings) warning("Logistic regression did not converge.")
  }
  
  return(list(
    model = logreg,
    model_population = model_data,
    eval_population = eval_population,
    outcome = outcome,
    predictor = predictor,
    covariates = covariates,
    roc = roc_obj,
    pr = pr_obj,
    calibration = calibration_obj,
    converged = logreg$converged,
    subset_controls = subset_controls
  ))

}

## Confusion matrix

confusion_matrix <- function(data,
                             outcome,
                             cutoff_value = 0.5,
                             positive_level = 1,
                             negative_level = 0,
                             print_epi = TRUE,
                             print_cm = TRUE,
                             print_net_benefit = TRUE,
                             warn_missing_epi = TRUE) {

  if (cutoff_value <= 0 || cutoff_value >= 1) {
    stop("cutoff_value must be between 0 and 1 (exclusive).")
  }
  if (!("pred" %in% names(data))) {
    stop("data must contain a column named 'pred' with predicted probabilities.")
  }

  if (!outcome %in% names(data)) {
    stop(paste0("Outcome column not found in data: ", outcome))
  }

  valid_rows <- !is.na(data[[outcome]]) & !is.na(data$pred)
  data <- data[valid_rows, , drop = FALSE]

  if (nrow(data) == 0) {
    stop("No complete rows with non-missing outcome and predicted probabilities.")
  }

  # Predicted class at threshold (treat/positive if pred >= pt)
  data$predbin <- factor(
    data$pred >= cutoff_value,
    levels = c(TRUE, FALSE)
  )

  # Actual outcome (ensure positive then negative ordering)
  data$outcome_actual <- factor(
    data[[outcome]],
    levels = c(positive_level, negative_level)
  )

  # Confusion matrix: rows = predicted (TRUE/FALSE), cols = actual (pos/neg)
  cm <- table(
    pred_outcome = data$predbin,
    outcome_actual = data$outcome_actual
  )

  if (print_cm) print(cm)

  get_cell <- function(mat, r, c) {
    r <- as.character(r); c <- as.character(c)
    if (r %in% rownames(mat) && c %in% colnames(mat)) mat[r, c] else 0
  }

  TP <- get_cell(cm, TRUE,  positive_level)
  FP <- get_cell(cm, TRUE,  negative_level)
  FN <- get_cell(cm, FALSE, positive_level)
  TN <- get_cell(cm, FALSE, negative_level)
  N  <- sum(cm)

  epi <- NULL

  if (isTRUE(print_epi)) {

    epi_available <- tryCatch(
      requireNamespace("epiR", quietly = TRUE),
      error = function(e) FALSE
    )

    if (epi_available) {

      epi <- epiR::epi.tests(cm)
      print(epi)

    } else {

      if (isTRUE(warn_missing_epi)) {
        warning(
          "epi.tests is unavailable; returning manual diagnostic metrics."
        )
      }

      safe_div <- function(num, den) {
        ifelse(den == 0, NA_real_, num / den)
      }

      sensitivity <- safe_div(TP, TP + FN)
      specificity <- safe_div(TN, TN + FP)
      ppv <- safe_div(TP, TP + FP)
      npv <- safe_div(TN, TN + FN)
      accuracy <- safe_div(TP + TN, N)
      f1_score <- safe_div(
        2 * ppv * sensitivity,
        ppv + sensitivity
      )
      lr_pos <- safe_div(
        sensitivity,
        1 - specificity
      )
      lr_neg <- safe_div(
        1 - sensitivity,
        specificity
      )

      epi <- list(
        method = "manual",
        confusion_matrix = cm,
        counts = list(
          TP = TP,
          FP = FP,
          FN = FN,
          TN = TN,
          N = N
        ),
        metrics = list(
          sensitivity = sensitivity,
          specificity = specificity,
          ppv = ppv,
          npv = npv,
          accuracy = accuracy,
          f1_score = f1_score,
          lr_pos = lr_pos,
          lr_neg = lr_neg
        )
      )

      cat("\nManual diagnostic metrics:\n")
      print(epi$counts)
      print(epi$metrics)
    }
  }

  # ---- Net benefit calculations ----
  pt <- cutoff_value
  w  <- pt / (1 - pt)  # threshold odds / exchange rate  [1](https://www.bmj.com/content/352/bmj.i6)[2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6777022/)

  # Standard net benefit: TP/N - FP/N * w  [2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6777022/)[1](https://www.bmj.com/content/352/bmj.i6)
  net_benefit <- (TP / N) - (FP / N) * w

  # Reference strategies for decision curves
  prevalence <- sum(data$outcome_actual == as.character(positive_level), na.rm = TRUE) / N

  # Treat none: NB = 0 by definition (no positives, no false positives) [2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6777022/)[1](https://www.bmj.com/content/352/bmj.i6)
  nb_treat_none <- 0

  # Treat all: TP/N = prevalence; FP/N = (1 - prevalence) [2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6777022/)[1](https://www.bmj.com/content/352/bmj.i6)
  nb_treat_all <- prevalence - (1 - prevalence) * w

  if (print_net_benefit) {
    cat("\nThreshold (pt) =", pt,
        "\nNet benefit (model)     =", net_benefit,
        "\nNet benefit (treat all) =", nb_treat_all,
        "\nNet benefit (treat none)=", nb_treat_none,
        "\n")
  }

  return(list(
    confusion_matrix = cm,
    epi_tests = epi,
    net_benefit = net_benefit,
    net_benefit_treat_all = nb_treat_all,
    net_benefit_treat_none = nb_treat_none,
    components = list(
      TP = as.numeric(TP),
      FP = as.numeric(FP),
      FN = as.numeric(FN),
      TN = as.numeric(TN),
      N = as.numeric(N),
      N_cases = as.numeric(TP + FN),
      N_controls = as.numeric(TN + FP),
      pt = as.numeric(pt),
      w = as.numeric(w),
      prevalence = as.numeric(prevalence)
    ),
    data = data
  ))
}

## Net benefit decision curve function

nb_curve <- function(data,
                     outcome,
                     cutoff_values = seq(0.01, 0.99, by = 0.01),
                     positive_level = 1,
                     negative_level = 0,
                     plot_curve = TRUE,
                     plot_title = "Decision Curve Analysis") {

  if (!("pred" %in% names(data))) {
    stop("data must contain a column named 'pred' with predicted probabilities.")
  }

  
  if (!outcome %in% names(data)) {
    stop(paste0("Outcome column not found in data: ", outcome))
  }

  # validate that cutoff values are unique and in ascending order

  if (
    !is.numeric(cutoff_values) ||
    length(cutoff_values) == 0L ||
    any(!is.finite(cutoff_values))
  ) {
    stop(
      "cutoff_values must be a non-empty numeric vector ",
      "containing only finite values."
    )
  }

  if (any(cutoff_values <= 0 | cutoff_values >= 1)) {
    stop(
      "All cutoff_values must be between 0 and 1, exclusive."
    )
  }

  if (is.unsorted(cutoff_values, strictly = TRUE)) {
    stop(
      "cutoff_values must be unique and supplied in ",
      "strictly increasing order."
    )
  }

  # Keep only rows with non-missing outcome and predicted probability
  valid_rows <- !is.na(data[[outcome]]) & !is.na(data$pred)
  data <- data[valid_rows, , drop = FALSE]

  if (nrow(data) == 0) {
    stop("No complete rows with non-missing outcome and predicted probabilities.")
  }

  # Check that outcome contains both positive and negative classes
  n_cases <- sum(data[[outcome]] == positive_level, na.rm = TRUE)
  n_controls <- sum(data[[outcome]] == negative_level, na.rm = TRUE)

  if (n_cases == 0 || n_controls == 0) {
    stop("Decision curve requires at least one case and one control among valid rows.")
  }

  ##

  nb_tolerance <- 1e-12

  results <- lapply(cutoff_values, function(pt) {

    cm_result <- confusion_matrix(
      data = data,
      outcome = outcome,
      cutoff_value = pt,
      positive_level = positive_level,
      negative_level = negative_level,
      print_epi = FALSE,
      print_cm = FALSE,
      print_net_benefit = FALSE,
      warn_missing_epi = FALSE
    )

    TP <- as.numeric(cm_result$components$TP)
    FP <- as.numeric(cm_result$components$FP)
    FN <- as.numeric(cm_result$components$FN)
    TN <- as.numeric(cm_result$components$TN)

    safe_div <- function(num, den) {
      if (den == 0) NA_real_ else num / den
    }

    sensitivity <- safe_div(TP, TP + FN)
    specificity <- safe_div(TN, TN + FP)

    data.frame(
      pt = as.numeric(pt),

      N = as.numeric(TP + FP + FN + TN),
      N_cases = as.numeric(TP + FN),
      N_controls = as.numeric(TN + FP),

      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,

      sensitivity = sensitivity,
      specificity = specificity,
      ppv = safe_div(TP, TP + FP),
      npv = safe_div(TN, TN + FN),

      # Additional cutoff metrics
      treatment_equality = safe_div(FN, FP),
      accuracy = safe_div(TP + TN, TP + FP + FN + TN),
      balanced_accuracy = safe_div(
        sensitivity + specificity,
        2
      ),

      net_benefit_model = as.numeric(
        cm_result$net_benefit
      ),
      net_benefit_treat_all = as.numeric(
        cm_result$net_benefit_treat_all
      ),
      net_benefit_treat_none = as.numeric(
        cm_result$net_benefit_treat_none
      ),

      stringsAsFactors = FALSE
    )
  })

  results_df <- dplyr::bind_rows(results)

  stopifnot(
  all(results_df$TP + results_df$FN == results_df$N_cases),
  all(results_df$TN + results_df$FP == results_df$N_controls),
  all(
    results_df$TP +
      results_df$FP +
      results_df$FN +
      results_df$TN ==
      results_df$N
  )
)

  if (plot_curve) {
    ## Explicitly include zero in the y axis range so that the treat-none line is always visible
    y_values <- c(
      results_df$net_benefit_model,
      results_df$net_benefit_treat_all,
      results_df$net_benefit_treat_none,
      0
    )

    y_values <- y_values[is.finite(y_values)]

    if (length(y_values) == 0L) {
      stop("No finite net-benefit values are available for plotting.")
    }

    ## Upper limit based on the largest observed net benefit
    y_upper <- max(y_values)

    ## Add a small amount of space above the highest curve
    upper_padding <- max(
      0.05 * abs(y_upper),
      0.002
    )

    y_upper_limit <- y_upper + upper_padding

    ## Show only a small amount below zero
    negative_margin <- max(
      0.10 * abs(y_upper_limit),
      0.005
    )

    y_lower_limit <- -negative_margin

    y_limits <- c(
      y_lower_limit,
      y_upper_limit
    )

    plot(
      results_df$pt,
      results_df$net_benefit_model,
      type = "l",
      col = "blue",
      lwd = 2,
      ylim = y_limits,
      xlab = "Threshold Probability",
      ylab = "Net Benefit",
      main = plot_title
    )

    lines(
      results_df$pt,
      results_df$net_benefit_treat_all,
      col = "red",
      lwd = 2
    )

    ## Treat none is always zero, so draw it as a horizontal reference line
    
    abline(
      h = 0,
      col = "green",
      lwd = 2
    )

    legend(
      "topright",
      legend = c("Model", "Treat All", "Treat None"),
      col = c("blue", "red", "green"),
      lwd = 2,
      bty = "n"
    )

    # Assign key thresholds for reporting

    ## A useful threshold is one at which the model has greater
    ## net benefit than both treat-all and treat-none

    model_useful_idx <- which(
      results_df$net_benefit_model >
        results_df$net_benefit_treat_all + nb_tolerance &
      results_df$net_benefit_model >
        results_df$net_benefit_treat_none + nb_tolerance
    )

    ## C1: first threshold where the model beats both alternatives

    idx_start_useful <- if (length(model_useful_idx) > 0L) {
      model_useful_idx[1]
    } else {
      NA_integer_
    }

    ## C2: threshold that maximises Youden's J

    youden_results <- data.frame(
      idx = seq_len(nrow(results_df)),
      youden_j = results_df$sensitivity +
        results_df$specificity - 1,
      stringsAsFactors = FALSE
    )

    valid_youden_idx <- which(
      is.finite(youden_results$youden_j)
    )

    idx_youden <- if (length(valid_youden_idx) > 0L) {
      valid_youden_idx[
        which.max(
          youden_results$youden_j[valid_youden_idx]
        )
      ]
    } else {
      NA_integer_
    }

    ## C3: final candidate threshold before the model first falls
    ## below treat-none after C1

    idx_first_not_better_than_none <- NA_integer_
    idx_end_useful <- NA_integer_

    if (
      !is.na(idx_start_useful) &&
      idx_start_useful < nrow(results_df)
    ) {

      candidate_idx <- seq.int(
        from = idx_start_useful + 1L,
        to = nrow(results_df)
      )

      model_not_better_than_none <- (
        results_df$net_benefit_model[candidate_idx] <=
          results_df$net_benefit_treat_none[candidate_idx] +
            nb_tolerance
      )

      first_boundary_position <- which(
        model_not_better_than_none
      )[1]

      if (
        length(first_boundary_position) == 1L &&
        !is.na(first_boundary_position)
      ) {
        idx_first_not_better_than_none <-
          candidate_idx[first_boundary_position]

        ## C3 is the immediately preceding candidate threshold
        idx_end_useful <- idx_first_not_better_than_none - 1L
      }
    }

    ## Create summary table of key cutoff diagnostics

    key_idx <- c(
      idx_start_useful,
      idx_youden,
      idx_end_useful
    )

    key_labels <- c(
      "Start of useful range",
      "Maximum Youden J",
      "End of useful range"
    )

    valid_key_rows <- !is.na(key_idx)

    key_idx <- key_idx[valid_key_rows]
    key_labels <- key_labels[valid_key_rows]


    key_summary <- dplyr::bind_rows(
      lapply(seq_along(key_idx), function(i) {

        idx <- key_idx[i]

        data.frame(
          Scenario = key_labels[i],
          Cutoff = as.numeric(results_df$pt[idx]),

          N = as.numeric(results_df$N[idx]),
          N_Cases = as.numeric(results_df$N_cases[idx]),
          N_Controls = as.numeric(results_df$N_controls[idx]),

          TP = as.numeric(results_df$TP[idx]),
          FP = as.numeric(results_df$FP[idx]),
          FN = as.numeric(results_df$FN[idx]),
          TN = as.numeric(results_df$TN[idx]),

          Sensitivity = as.numeric(results_df$sensitivity[idx]),
          Specificity = as.numeric(results_df$specificity[idx]),
          PPV = as.numeric(results_df$ppv[idx]),
          NPV = as.numeric(results_df$npv[idx]),

          Treatment_Equality = as.numeric(results_df$treatment_equality[idx]),
          Accuracy = as.numeric(results_df$accuracy[idx]),
          Balanced_Accuracy = as.numeric(results_df$balanced_accuracy[idx]),

          Youden_J = if (
            is.finite(youden_results$youden_j[idx])
          ) {
            as.numeric(youden_results$youden_j[idx])
          } else {
            NA_real_
          },

          Net_Benefit_Model = as.numeric(
            results_df$net_benefit_model[idx]
          ),
          Net_Benefit_Treat_All = as.numeric(
            results_df$net_benefit_treat_all[idx]
          ),
          Net_Benefit_Treat_None = as.numeric(
            results_df$net_benefit_treat_none[idx]
          ),

          stringsAsFactors = FALSE
        )

  if (length(custom_columns) > 0L) {

    custom_df <- as.data.frame(
      custom_columns,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    summary_row <- cbind(
      summary_row,
      custom_df
    )
  }

  summary_row

      })
    )

    cat("\nKey cutoff diagnostics for plotted decision curve:\n")
    print(key_summary, row.names = FALSE)
  }

  return(results_df)
}

###########################################
#============= Bulk functions ============#
###########################################

## helper function to build populations for a given screening date, optionally requiring genetic data

build_populations_for_screening_date <- function(screening_date,
                                                 require_genetic_data = TRUE,
                                                 genetic_data_marker = "ContiMultiethnicGRS") {

  base_data <- PrCa_prediagnosis_males %>%
    dplyr::filter(!is.na(.data[[screening_date]]))

  if (isTRUE(require_genetic_data)) {
    base_data <- base_data %>%
      dplyr::filter(!is.na(.data[[genetic_data_marker]]))
  }

  PrCa_predhorizon_bulk <- base_data %>%
    dplyr::mutate(
      PrCa_post_assessment = as.integer(
        !is.na(earliest_PrCa_date) &
          earliest_PrCa_date >= .data[[screening_date]]
      ),
      PrCa_2yrs = as.integer(
        !is.na(earliest_PrCa_date) &
          earliest_PrCa_date >= .data[[screening_date]] &
          earliest_PrCa_date <= .data[[screening_date]] %m+% lubridate::years(2)
      ),
      PrCa_5yrs = as.integer(
        !is.na(earliest_PrCa_date) &
          earliest_PrCa_date >= .data[[screening_date]] &
          earliest_PrCa_date <= .data[[screening_date]] %m+% lubridate::years(5)
      ),
      PrCa_10yrs = as.integer(
        !is.na(earliest_PrCa_date) &
          earliest_PrCa_date >= .data[[screening_date]] &
          earliest_PrCa_date <= .data[[screening_date]] %m+% lubridate::years(10)
      )
    )

  if (screening_date == "artificial_screening_date") {
    PrCa_predhorizon_allmen_bulk <- PrCa_predhorizon_bulk %>%
      dplyr::filter(
        demog_sex == 1,
        HES_only == FALSE,
        death_only == FALSE,
        pre_diagnosed_asd == FALSE,
        over_20_at_asd == TRUE,
        died_by_asd == FALSE
      )
  } else if (screening_date == "date40") {
    PrCa_predhorizon_allmen_bulk <- PrCa_predhorizon_bulk %>%
      dplyr::filter(
        demog_sex == 1,
        HES_only == FALSE,
        death_only == FALSE,
        pre_diagnosed_date40 == FALSE,
        over_40 == TRUE,
        died_by_40 == FALSE
      )
  } else if (screening_date == "date45") {
    PrCa_predhorizon_allmen_bulk <- PrCa_predhorizon_bulk %>%
      dplyr::filter(
        demog_sex == 1,
        HES_only == FALSE,
        death_only == FALSE,
        pre_diagnosed_date45 == FALSE,
        over_45 == TRUE,
        died_by_45 == FALSE
      )
  } else if (screening_date == "date50") {
    PrCa_predhorizon_allmen_bulk <- PrCa_predhorizon_bulk %>%
      dplyr::filter(
        demog_sex == 1,
        HES_only == FALSE,
        death_only == FALSE,
        pre_diagnosed_date50 == FALSE,
        over_50 == TRUE,
        died_by_50 == FALSE
      )
  } else if (screening_date == "date55") {
    PrCa_predhorizon_allmen_bulk <- PrCa_predhorizon_bulk %>%
      dplyr::filter(
        demog_sex == 1,
        HES_only == FALSE,
        death_only == FALSE,
        pre_diagnosed_date55 == FALSE,
        over_55 == TRUE,
        died_by_55 == FALSE
      )
  } else if (screening_date == "date60") {
    PrCa_predhorizon_allmen_bulk <- PrCa_predhorizon_bulk %>%
      dplyr::filter(
        demog_sex == 1,
        HES_only == FALSE,
        death_only == FALSE,
        pre_diagnosed_date60 == FALSE,
        over_60 == TRUE,
        died_by_60 == FALSE
      )
  } else if (screening_date == "date65") {
    PrCa_predhorizon_allmen_bulk <- PrCa_predhorizon_bulk %>%
      dplyr::filter(
        demog_sex == 1,
        HES_only == FALSE,
        death_only == FALSE,
        pre_diagnosed_date65 == FALSE,
        over_65 == TRUE,
        died_by_65 == FALSE
      )
  } else {
    stop(paste0("Unsupported screening_date: ", screening_date))
  }

  list(
    "All" = PrCa_predhorizon_allmen_bulk,
    "White" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(ethnicity_group_narrow == "White"),
    "Black" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(ethnicity_group_narrow == "Black"),
    "Mixed" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(ethnicity_group_narrow == "Mixed White and Black"),
    "EUR" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(genetic_similarity_group == "EUR"),
    "AFR" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(genetic_similarity_group == "AFR"),
    "EAS" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(genetic_similarity_group == "EAS"),
    "SAS" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(genetic_similarity_group == "SAS"),
    "MID" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(genetic_similarity_group == "MID"),
    "AMR" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(genetic_similarity_group == "AMR"),
    "SouthAsian" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(ethnicity_group_narrow == "South Asian"),
    "EastAsian" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(ethnicity_group_narrow == "Chinese"),
    "Arab" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(ethnicity_group_narrow == "Arab"),
    "Other" = PrCa_predhorizon_allmen_bulk %>%
      dplyr::filter(ethnicity_group_narrow == "Other")
  )
}

## Bulk run logistic regressions for selected screening dates

logreg_table <- function(
    screening_dates = c("artificial_screening_date", "date40", "date45", "date50", "date55", "date60", "date65"),
    populations = c("All", "White", "Black", "Mixed", "EUR", "AFR", "EAS", "SAS", "MID", "AMR", "SouthAsian", "EastAsian", "Arab", "Other"),
    predictors = NULL,
    outcomes = c("PrCa_post_assessment", "PrCa_2yrs", "PrCa_5yrs", "PrCa_10yrs"),
    covariates_list = list(NULL),
    subset_controls = FALSE,
    plot_roc = FALSE,
    verbose = TRUE,
    adjust_for_PCs = FALSE
) {

  # Default GRS list if not specified
  if (is.null(predictors)) {
    predictors <- c(
      # Conti GRSs
      "ContiMultiethnicGRS",
      "ContiEURGRS",
      "ContiAFRGRS",
      "ContiEASGRS",
      "ContiHISGRS",
      # Wang GRSs
      "WangMultiethnicGRS",
      "WangEURGRS",
      "WangAFRGRS",
      "WangEASGRS",
      "WangHISGRS"
    )
  }

  pc_covariates <- paste0("PC", 1:10)

  extract_roc_metrics <- function(roc_obj) {
    if (is.null(roc_obj)) {
      return(list(auc = NA_real_, lower = NA_real_, upper = NA_real_))
    }

    auc_val <- suppressWarnings(as.numeric(roc_obj$auc)[1])
    ci_vals <- suppressWarnings(as.numeric(roc_obj$ci))

    if (length(ci_vals) >= 3) {
      lower <- ci_vals[1]
      upper <- ci_vals[3]
    } else if (length(ci_vals) == 2) {
      lower <- ci_vals[1]
      upper <- ci_vals[2]
    } else {
      lower <- NA_real_
      upper <- NA_real_
    }

    if (!is.finite(auc_val)) auc_val <- NA_real_
    if (!is.finite(lower)) lower <- NA_real_
    if (!is.finite(upper)) upper <- NA_real_

    list(auc = auc_val, lower = lower, upper = upper)
  }

  parse_covariate_terms <- function(covariate_spec) {
    if (is.null(covariate_spec)) {
      return(character(0))
    }

    terms <- trimws(unlist(strsplit(as.character(covariate_spec), "\\+")))
    terms <- terms[nzchar(terms)]
    unique(terms)
  }

  resolve_covariates <- function(base_covariates, predictor) {
    covs <- parse_covariate_terms(base_covariates)

    if (isTRUE(adjust_for_PCs)) {
      covs <- c(covs, pc_covariates)
    }

    covs <- unique(covs)
    covs <- covs[covs != predictor]

    if (length(covs) == 0) {
      return(NULL)
    }

    covs
  }

  covariate_label <- function(covs) {
    if (is.null(covs) || length(covs) == 0) {
      return("None")
    }

    paste(covs, collapse = " + ")
  }

  result_rows <- list()
  result_count <- 0L

  total_combos <- length(screening_dates) *
    length(populations) *
    length(outcomes) *
    length(predictors) *
    length(covariates_list)

  combo_count <- 0

  seen_model_specs <- new.env(parent = emptyenv())

  for (sd in screening_dates) {
    if (isTRUE(verbose)) {
      cat(sprintf("\n===== Screening date: %s =====\n", sd))
      flush.console()
    }

    
      pops <- build_populations_for_screening_date(
        screening_date = sd,
        require_genetic_data = TRUE,
        genetic_data_marker = "ContiMultiethnicGRS"
      )


    missing_keys <- setdiff(populations, names(pops))
    if (length(missing_keys) > 0) {
      stop(sprintf("Unknown population key(s): %s", paste(missing_keys, collapse = ", ")))
    }

    pops_selected <- pops[populations]

    for (current_outcome in outcomes) {
      for (pop_name in names(pops_selected)) {
        pop_data <- pops_selected[[pop_name]]

        if (is.null(pop_data) || nrow(pop_data) == 0) {
          if (isTRUE(verbose)) {
            cat(sprintf(
              "Skipping empty population: %s at %s\n",
              pop_name,
              sd
            ))
            flush.console()
          }
          next
        }

        for (predictor in predictors) {

          if (!predictor %in% colnames(pop_data)) {
            if (isTRUE(verbose)) {
              cat(sprintf(
                "Skipping %s (not in %s, %s)\n",
                predictor,
                pop_name,
                sd
              ))
              flush.console()
            }
            next
          }

          for (cov in covariates_list) {
            combo_count <- combo_count + 1

            effective_covariates <- resolve_covariates(cov, predictor = predictor)
            cov_label <- covariate_label(effective_covariates)
            
            if (uses_age_asd_inappropriately(
              screening_date = sd,
              predictor = predictor,
              covariates = effective_covariates
        )) {

              message(age_asd_skip_message)
              next

            }

            model_spec_key <- paste(
              sd,
              current_outcome,
              pop_name,
              predictor,
              cov_label,
              sep = "__"
            )

            if (isTRUE(verbose)) {
              cat(sprintf(
                "[%d/%d] %s | %s | %s | %s | %s\n",
                combo_count,
                total_combos,
                sd,
                current_outcome,
                pop_name,
                predictor,
                cov_label
              ))
              flush.console()
            }

            missing_covariates <- setdiff(effective_covariates, colnames(pop_data))
            if (length(missing_covariates) > 0) {
              if (isTRUE(verbose)) {
                cat(sprintf(
                  "  Skipping (missing covariates in %s at %s: %s)\n",
                  pop_name,
                  sd,
                  paste(missing_covariates, collapse = ", ")
                ))
              }
              next
            }
  
            if (exists(model_spec_key, envir = seen_model_specs, inherits = FALSE)) {
                          if (isTRUE(verbose)) {
                            cat(sprintf(
                              "  Skipping duplicate resolved model: %s | %s | %s | %s | %s\n",
                              sd,
                              current_outcome,
                              pop_name,
                              predictor,
                              cov_label
                            ))
                          }
                          next
                        }

            assign(model_spec_key, TRUE, envir = seen_model_specs)

            tryCatch({
              model <- run_logreg(
                data = pop_data,
                outcome = current_outcome,
                predictor = predictor,
                covariates = effective_covariates,
                screening_date = sd,
                subset_controls = subset_controls,
                plot_roc = plot_roc,
                plot_pr = TRUE,
                show_output = FALSE
              )

              roc_obj <- model$roc
              roc_metrics <- extract_roc_metrics(roc_obj)

              auc_val <- roc_metrics$auc
              auc_lower <- roc_metrics$lower
              auc_upper <- roc_metrics$upper

              pr_auc <- if (!is.null(model$pr) && !is.null(model$pr$auc.integral)) {
                as.numeric(model$pr$auc.integral)
              } else {
                NA_real_
              }

              valid_rows <- !is.na(model$data[[current_outcome]]) & !is.na(model$data$pred)
              prevalence <- mean(model$data[[current_outcome]][valid_rows] == 1)
              mean_residual <- mean(model$data$residual[valid_rows], na.rm = TRUE)

              n_cases <- sum(model$data[[current_outcome]][valid_rows] == 1, na.rm = TRUE)
              n_controls <- sum(model$data[[current_outcome]][valid_rows] == 0, na.rm = TRUE)

              coef_tbl <- summary(model$model)$coefficients

              ## Calculate SD of the globally z-transformed predictor within the current model population (cases and controls included)
              ## (this allows calculation of subpopulation-specific OR per 1 SD increase in predictor, as well as OR per global SD)
              
              model_vars_for_sd <- c(current_outcome, predictor, effective_covariates)

              complete_rows_for_sd <- stats::complete.cases(
                model$data[, model_vars_for_sd, drop = FALSE]
              )

              rows_for_sd <- complete_rows_for_sd &
                !is.na(model$data[[current_outcome]])

              subpop_sd <- stats::sd(
                model$data[[predictor]][rows_for_sd],
                na.rm = TRUE
              )

              if (!is.finite(subpop_sd) || subpop_sd <= 0) {
                subpop_sd <- NA_real_
              }

              if (!grepl("GRS$", predictor)) { ## Ensure that non-GRS predictors (e.g., age) are not rescaled to the subpopulation
                  subpop_sd <- NA_real_
              }

              ## Calculate OR per 1 SD increase in predictor (z-transformed to global male population and subpopulation-specific)

              if (predictor %in% rownames(coef_tbl)) {
                beta <- coef_tbl[predictor, "Estimate"]
                se <- coef_tbl[predictor, "Std. Error"]
                p_val <- coef_tbl[predictor, "Pr(>|z|)"]

                ## OR per 1 global-population SD.
                ## (assuming the predictor has already been globally z-transformed)
                or_global <- exp(beta)
                or_global_ci_lower <- exp(beta - 1.96 * se)
                or_global_ci_upper <- exp(beta + 1.96 * se)

                ## OR per 1 subpopulation SD.
                ## Since the predictor is globally z-transformed, subpop_sd is usually not exactly 1.
                if (is.finite(subpop_sd)) {
                  or_subpop <- exp(beta * subpop_sd)
                  or_subpop_ci_lower <- exp((beta - 1.96 * se) * subpop_sd)
                  or_subpop_ci_upper <- exp((beta + 1.96 * se) * subpop_sd)
                  p_subpop <- p_val
                } else {
                  or_subpop <- NA_real_
                  or_subpop_ci_lower <- NA_real_
                  or_subpop_ci_upper <- NA_real_
                  p_subpop <- NA_real_
                }

              } else {
                or_global <- NA_real_
                or_global_ci_lower <- NA_real_
                or_global_ci_upper <- NA_real_
                p_val <- NA_real_

                or_subpop <- NA_real_
                or_subpop_ci_lower <- NA_real_
                or_subpop_ci_upper <- NA_real_
                p_subpop <- NA_real_
              }

            result_count <- result_count + 1L

            result_rows[[result_count]] <- data.frame(
              Screening_Date = sd,
              Outcome = current_outcome,
              Population = pop_name,
              Predictor = predictor,
              Covariates = cov_label,
              N_Cases = n_cases,
              N_Controls = n_controls,
              Prevalence = prevalence,
              mean_residual = mean_residual,
              PR_AUC = pr_auc,
              ROC_AUC = auc_val,
              ROC_AUC_CI_Lower = auc_lower,
              ROC_AUC_CI_Upper = auc_upper,

              OR_per_1SD_global = or_global,
              OR_per_1SD_global_CI_Lower = or_global_ci_lower,
              OR_per_1SD_global_CI_Upper = or_global_ci_upper,
              OR_per_1SD_global_p = p_val,

              Subpopulation_SD = subpop_sd,
              OR_per_1SD_subpopulation = or_subpop,
              OR_per_1SD_subpopulation_CI_Lower = or_subpop_ci_lower,
              OR_per_1SD_subpopulation_CI_Upper = or_subpop_ci_upper,
              OR_per_1SD_subpopulation_p = p_subpop,

              Converged = model$converged,
              stringsAsFactors = FALSE
            )

            if (isTRUE(verbose)) {
              cat("  Completed successfully.\n")
              flush.console()
            }

            rm(
              model,
              roc_obj,
              roc_metrics,
              coef_tbl,
              valid_rows,
              complete_rows_for_sd,
              rows_for_sd
            )

            invisible(gc(verbose = FALSE))

            }, error = function(e) {
              if (isTRUE(verbose)) {
                cat(sprintf(
                  "  ERROR: %s | %s | %s | %s | %s: %s\n",
                  sd,
                  current_outcome,
                  pop_name,
                  predictor,
                  cov_label,
                  conditionMessage(e)
                ))
                flush.console()
              }

              invisible(gc(verbose = FALSE))
            })
          }
        }
      }
    }
  }

  all_results <- if (length(result_rows) > 0L) {
    dplyr::bind_rows(result_rows)
  } else {
    data.frame(
      Screening_Date = character(),
      Outcome = character(),
      Population = character(),
      Predictor = character(),
      Covariates = character(),
      N_Cases = integer(),
      N_Controls = integer(),
      Prevalence = numeric(),
      PR_AUC = numeric(),
      ROC_AUC = numeric(),
      ROC_AUC_CI_Lower = numeric(),
      ROC_AUC_CI_Upper = numeric(),
      OR_per_1SD_global = numeric(),
      OR_per_1SD_global_CI_Lower = numeric(),
      OR_per_1SD_global_CI_Upper = numeric(),
      OR_per_1SD_global_p = numeric(),
      Subpopulation_SD = numeric(),
      OR_per_1SD_subpopulation = numeric(),
      OR_per_1SD_subpopulation_CI_Lower = numeric(),
      OR_per_1SD_subpopulation_CI_Upper = numeric(),
      OR_per_1SD_subpopulation_p = numeric(),
      Converged = logical(),
      stringsAsFactors = FALSE
    )
  }

  return(all_results)
}

## Bulk run net benefit curves for selected screening dates, populations, and predictors

bulk_nb_curve <- function(populations = c("All", "White", "Black", "Mixed", "EUR", "AFR", "EAS", "SAS", "MID", "AMR"),
                          outcomes = c("PrCa_post_assessment", "PrCa_2yrs", "PrCa_5yrs", "PrCa_10yrs"),
                          predictors,
                          covariates_list = list(NULL),
                          screening_dates = c("artificial_screening_date", "date40", "date45", "date50", "date55", "date60", "date65"),
                          cutoff_values = seq(0.01, 0.99, by = 0.01),
                          positive_level = 1,
                          negative_level = 0,
                          subset_controls = FALSE,
                          print_all_curves = FALSE,
                          summary_table = TRUE,
                          output_csv = "bulk_nb_summary_table.csv",
                          verbose = TRUE,
                          retain_model_runs = FALSE,
                          suppress_counts = TRUE,         # Set to TRUE for rows with TP/TN/FP/FN counts less than 10 to have algorithmic fairness metrics suppressed
                          custom_thresholds = NULL        # supply as a set (e.g., custom_thresholds = c(0.1, 0.2, 0.3))
                          ) {

  if (length(populations) == 0) stop("Provide at least one population.")
  if (length(outcomes) == 0) stop("Provide at least one outcome.")
  if (length(predictors) == 0) stop("Provide at least one predictor.")

    if (is.null(custom_thresholds)) {

    custom_thresholds <- numeric(0)

  } else {

    if (!is.numeric(custom_thresholds)) {
      stop("custom_thresholds must be NULL or a numeric vector.")
    }

    if (
      anyNA(custom_thresholds) ||
      any(!is.finite(custom_thresholds))
    ) {
      stop("custom_thresholds cannot contain NA, NaN, or infinite values.")
    }

    if (any(custom_thresholds <= 0 | custom_thresholds >= 1)) {
      stop("All custom_thresholds must be greater than 0 and less than 1.")
    }

    custom_thresholds <- sort(unique(as.numeric(custom_thresholds)))
  }

  if (
    anyNA(cutoff_values) ||
    any(!is.finite(cutoff_values)) ||
    any(cutoff_values <= 0 | cutoff_values >= 1)
  ) {
    stop("All cutoff_values must be finite values greater than 0 and less than 1.")
  }

  cutoff_values <- sort(unique(as.numeric(cutoff_values)))

  ## The original cutoff_values will still be used to identify C1, C2, and C3.
  evaluated_cutoff_values <- sort(unique(c(
    cutoff_values,
    custom_thresholds
  )))

  pop_alias <- c(
    all = "All", allmen = "All",
    white = "White", whitemen = "White",
    black = "Black", blackmen = "Black",
    mixed = "Mixed", mixedmen = "Mixed",
    southasian = "SouthAsian", southasianmen = "SouthAsian",
    eastasian = "EastAsian", eastasianmen = "EastAsian",
    arab = "Arab", arabmen = "Arab",
    other = "Other", othermen = "Other",
    eur = "EUR", afri = "AFR", afr = "AFR", eas = "EAS", sas = "SAS", mid = "MID", amr = "AMR"
  )

  pop_input <- tolower(populations)
  bad_pops <- setdiff(pop_input, names(pop_alias))
  if (length(bad_pops) > 0) {
    stop(paste0("Unknown population(s): ", paste(bad_pops, collapse = ", ")))
  }
  selected_pops <- unique(unname(pop_alias[pop_input]))

  total_combos <- length(screening_dates) *
    length(selected_pops) *
    length(outcomes) *
    length(predictors) *
    length(covariates_list)

  combo_count <- 0L

  results <- list()
  summary_rows <- list()

  seen_model_specs <- new.env(parent = emptyenv())

  covariate_label <- function(covs) {
    if (is.null(covs) || length(covs) == 0) {
      return("None")
    }
    paste(covs, collapse = " + ")
  }

    parse_covariate_terms <- function(covariate_spec) {
    if (is.null(covariate_spec)) {
      return(character(0))
    }

    terms <- trimws(unlist(strsplit(as.character(covariate_spec), "\\+")))
    terms <- terms[nzchar(terms)]
    unique(terms)
  }

  resolve_covariates <- function(base_covariates, predictor) {
    covs <- parse_covariate_terms(base_covariates)

    covs <- unique(covs)
    covs <- covs[covs != predictor]

    if (length(covs) == 0) {
      return(NULL)
    }

    covs
  }


summarise_nb_cutoffs <- function(nb_results,
                                 outcome,
                                 screening_date,
                                 population,
                                 predictor,
                                 suppress_counts = suppress_counts,
                                 standard_thresholds,
                                 custom_thresholds = NULL
                                 ) {

  nb_tolerance <- 1e-12

  ## Retain all evaluated results so that custom thresholds can be extracted.
  all_nb_results <- nb_results

  threshold_matches <- function(x, targets, tolerance = 1e-12) {

    if (length(targets) == 0L) {
      return(rep(FALSE, length(x)))
    }

    vapply(
      x,
      function(value) {
        any(abs(value - targets) <= tolerance)
      },
      logical(1)
    )
  }

  ## Restrict the existing C1/C2/C3 selection to the original cutoff grid.
  standard_idx <- which(threshold_matches(
    x = all_nb_results$pt,
    targets = standard_thresholds,
    tolerance = nb_tolerance
  ))

  if (length(standard_idx) == 0L) {
    stop("None of the standard cutoff_values were found in nb_results$pt.")
  }

  nb_results <- all_nb_results[standard_idx, , drop = FALSE]

  ## Ensure thresholds are in ascending order before identifying C1 and C3.
  nb_results <- nb_results[
    order(nb_results$pt),
    ,
    drop = FALSE
  ]

# Useful range: thresholds where model beats both treat all and treat none

useful_idx <- which(
  nb_results$net_benefit_model >
    nb_results$net_benefit_treat_all + nb_tolerance &
  nb_results$net_benefit_model >
    nb_results$net_benefit_treat_none + nb_tolerance
)

  ## Select the cutoff that maximises Youden's J across all candidate cutoffs
  youden_results <- data.frame(
    idx = seq_len(nrow(nb_results)),
    youden_j = nb_results$sensitivity +
      nb_results$specificity - 1,
    stringsAsFactors = FALSE
  )


  valid_youden_idx <- which(is.finite(youden_results$youden_j))
  idx_c2 <- if (length(valid_youden_idx) > 0) {
    valid_youden_idx[which.max(youden_results$youden_j[valid_youden_idx])]
  } else {
    NA_integer_
  }

  ## C1: first threshold where the model beats both treat-all and treat-none (start of useful range)

  idx_c1 <- if (length(useful_idx) > 0L) {
    useful_idx[1]
  } else {
    NA_integer_
  }

  ## C3: final candidate threshold before the model first falls
  ## below treat-none after C1

  idx_first_not_better_than_none <- NA_integer_
  idx_c3 <- NA_integer_

  if (
    !is.na(idx_c1) &&
    idx_c1 < nrow(nb_results)
  ) {

    candidate_idx <- seq.int(
      from = idx_c1 + 1L,
      to = nrow(nb_results)
    )

    model_not_better_than_none <- (
      nb_results$net_benefit_model[candidate_idx] <=
        nb_results$net_benefit_treat_none[candidate_idx] +
          nb_tolerance
    )

    first_boundary_position <- which(
      model_not_better_than_none
    )[1]

    if (
      length(first_boundary_position) == 1L &&
      !is.na(first_boundary_position)
    ) {
      idx_first_not_better_than_none <-
        candidate_idx[first_boundary_position]

      ## Final evaluated threshold before the model is no longer
      ## better than treat-none
      idx_c3 <- idx_first_not_better_than_none - 1L
    }
  }


  ##

  get_nb_row <- function(idx, results_data = nb_results) {

    if (
      length(idx) != 1L ||
      is.na(idx) ||
      idx < 1L ||
      idx > nrow(results_data)
    ) {
      return(list(
        cutoff = NA_real_,

        N = NA_real_,
        N_cases = NA_real_,
        N_controls = NA_real_,

        TP = NA_real_,
        FP = NA_real_,
        FN = NA_real_,
        TN = NA_real_,

        sensitivity = NA_real_,
        specificity = NA_real_,
        ppv = NA_real_,
        npv = NA_real_,

        treatment_equality = NA_real_,
        accuracy = NA_real_,
        balanced_accuracy = NA_real_,

        net_benefit_model = NA_real_,
        net_benefit_treat_all = NA_real_,
        net_benefit_treat_none = NA_real_
      ))
    }

    list(
      cutoff = as.numeric(results_data$pt[idx]),

      N = as.numeric(results_data$N[idx]),
      N_cases = as.numeric(results_data$N_cases[idx]),
      N_controls = as.numeric(results_data$N_controls[idx]),

      TP = as.numeric(results_data$TP[idx]),
      FP = as.numeric(results_data$FP[idx]),
      FN = as.numeric(results_data$FN[idx]),
      TN = as.numeric(results_data$TN[idx]),

      sensitivity = as.numeric(results_data$sensitivity[idx]),
      specificity = as.numeric(results_data$specificity[idx]),
      ppv = as.numeric(results_data$ppv[idx]),
      npv = as.numeric(results_data$npv[idx]),

      treatment_equality =
        as.numeric(results_data$treatment_equality[idx]),
      accuracy = as.numeric(results_data$accuracy[idx]),
      balanced_accuracy =
        as.numeric(results_data$balanced_accuracy[idx]),

      net_benefit_model =
        as.numeric(results_data$net_benefit_model[idx]),
      net_benefit_treat_all =
        as.numeric(results_data$net_benefit_treat_all[idx]),
      net_benefit_treat_none =
        as.numeric(results_data$net_benefit_treat_none[idx])
    )
  }

  ## Helper function - find custom thresholds 

  get_custom_threshold_idx <- function(
      threshold,
      tolerance = nb_tolerance
  ) {

    matching_idx <- which(
      abs(all_nb_results$pt - threshold) <= tolerance
    )

    if (length(matching_idx) == 0L) {
      return(NA_integer_)
    }

    ## If duplicated threshold rows somehow exist, use the closest one
    matching_idx[
      which.min(
        abs(all_nb_results$pt[matching_idx] - threshold)
      )
    ]
  }

  ## Helper function: identify whether a cutoff contains any reportable count from 1 to 9
  ##                  (missing cutoffs not treated as violations)

  ## Identify whether a cutoff contains any reportable count from 1 to 9.
## Missing cutoffs are not treated as violations because they were not evaluated
## as C1, C2, or C3.

check_under_10_rule <- function(
    diag,
    cutoff_found,
    minimum_reportable_count = 10L
) {

  if (!isTRUE(cutoff_found)) {
    return(FALSE)
  }

  counts_to_check <- c(
    N_cases = diag$N_cases,
    N_controls = diag$N_controls,
    TP = diag$TP,
    FP = diag$FP,
    FN = diag$FN,
    TN = diag$TN
  )

  ## Under this definition:
  ## 1 to 9 violate the rule;
  ## zero, NA, and counts of 10 or more do not.
  any(
    counts_to_check > 0 &
      counts_to_check < minimum_reportable_count,
    na.rm = TRUE
  )
}

  ## end of helper function

  c1_diag <- get_nb_row(idx_c1)
  c2_diag <- get_nb_row(idx_c2)
  c3_diag <- get_nb_row(idx_c3)

  ## Check each available cutoff for counts from 1 to 9

  c1_breaks_under_10 <- check_under_10_rule(
    diag = c1_diag,
    cutoff_found = !is.na(idx_c1)
  )

  c2_breaks_under_10 <- check_under_10_rule(
    diag = c2_diag,
    cutoff_found = !is.na(idx_c2)
  )

  c3_breaks_under_10 <- check_under_10_rule(
    diag = c3_diag,
    cutoff_found = !is.na(idx_c3)
  )

  suppress_diag_metrics <- function(diag) {

    diag$sensitivity <- NA_real_
    diag$specificity <- NA_real_
    diag$ppv <- NA_real_
    diag$npv <- NA_real_

    diag$treatment_equality <- NA_real_
    diag$accuracy <- NA_real_
    diag$balanced_accuracy <- NA_real_

    diag
  }

    custom_diags <- list()
  custom_breaks_under_10 <- logical(0)
  custom_threshold_found <- logical(0)

  if (length(custom_thresholds) > 0L) {

    custom_diags <- lapply(
      custom_thresholds,
      function(threshold) {

        idx <- get_custom_threshold_idx(threshold)

        get_nb_row(
          idx = idx,
          results_data = all_nb_results
        )
      }
    )

    custom_threshold_found <- vapply(
      custom_thresholds,
      function(threshold) {
        !is.na(get_custom_threshold_idx(threshold))
      },
      logical(1)
    )

    custom_breaks_under_10 <- vapply(
      seq_along(custom_diags),
      function(i) {
        check_under_10_rule(
          diag = custom_diags[[i]],
          cutoff_found = custom_threshold_found[i]
        )
      },
      logical(1)
    )

    if (isTRUE(suppress_counts)) {

      for (i in seq_along(custom_diags)) {

        if (custom_breaks_under_10[i]) {
          custom_diags[[i]] <-
            suppress_diag_metrics(custom_diags[[i]])
        }
      }
    }
  }

  if (isTRUE(suppress_counts)) {

    if (c1_breaks_under_10) {
      c1_diag <- suppress_diag_metrics(c1_diag)
    }

    if (c2_breaks_under_10) {
      c2_diag <- suppress_diag_metrics(c2_diag)
    }

    if (c3_breaks_under_10) {
      c3_diag <- suppress_diag_metrics(c3_diag)
    }
  }

  ## Combine affected cutoff labels into one field

  custom_threshold_labels <- if (
    length(custom_thresholds) == 0L
  ) {
    character(0)
  } else {
    paste0(
      "Custom_",
      format(
        custom_thresholds,
        scientific = FALSE,
        trim = TRUE
      )
    )
  }

  cutoffs_breaking_under_10 <- c(
    if (c1_breaks_under_10) "C1" else character(0),
    if (c2_breaks_under_10) "C2" else character(0),
    if (c3_breaks_under_10) "C3" else character(0),

    if (length(custom_breaks_under_10) > 0L) {
      custom_threshold_labels[custom_breaks_under_10]
    } else {
      character(0)
    }
  )

  breaks_under_10_rule <- if (
    length(cutoffs_breaking_under_10) == 0L
  ) {
    "None"
  } else {
    paste(cutoffs_breaking_under_10, collapse = ", ")
  }

  c1_cutoff <- c1_diag$cutoff
  c2_cutoff <- c2_diag$cutoff
  c3_cutoff <- c3_diag$cutoff

  fmt_cutoff <- function(x) {
    if (length(x) != 1L || is.na(x)) {
      "NA"
    } else {
      sprintf("%.2f", x)
    }
  }

  peak_nb_useful_range <- paste0(
    fmt_cutoff(c2_cutoff),
    " (",
    fmt_cutoff(c1_cutoff),
    " - ",
    fmt_cutoff(c3_cutoff),
    ")"
  )

  c2_youden_j <- if (is.na(idx_c2)) {
    NA_real_
  } else {
    as.numeric(youden_results$youden_j[idx_c2])
  }

  if (isTRUE(suppress_counts) && c2_breaks_under_10) {
    c2_youden_j <- NA_real_
  }

  make_threshold_prefix <- function(threshold) {

    threshold_text <- format(
      threshold,
      scientific = FALSE,
      trim = TRUE,
      digits = 15
    )

    ## Convert values such as 0.1 to a syntactically simple label:
    ## Custom_0p1
    threshold_text <- gsub("\\.", "p", threshold_text)
    threshold_text <- gsub("-", "m", threshold_text)

    paste0("Custom_", threshold_text)
  }

  custom_columns <- list()

  if (length(custom_thresholds) > 0L) {

    for (i in seq_along(custom_thresholds)) {

      prefix <- make_threshold_prefix(custom_thresholds[i])
      diag <- custom_diags[[i]]

      custom_columns[[paste0(prefix, "_found")]] <-
        custom_threshold_found[i]

      custom_columns[[paste0(prefix, "_cutoff")]] <-
        diag$cutoff

      custom_columns[[paste0(prefix, "_sensitivity")]] <-
        diag$sensitivity

      custom_columns[[paste0(prefix, "_specificity")]] <-
        diag$specificity

      custom_columns[[paste0(prefix, "_PPV")]] <-
        diag$ppv

      custom_columns[[paste0(prefix, "_NPV")]] <-
        diag$npv

      custom_columns[[paste0(
        prefix,
        "_Treatment_Equality"
      )]] <- diag$treatment_equality

      custom_columns[[paste0(prefix, "_Accuracy")]] <-
        diag$accuracy

      custom_columns[[paste0(
        prefix,
        "_Balanced_Accuracy"
      )]] <- diag$balanced_accuracy
    }
  }

  summary_row <- data.frame(
    Screening_Date = screening_date,
    Outcome = outcome,
    Population = population,
    Predictor = predictor,

    "Youden J cutoff (C1-C3)" = peak_nb_useful_range,

    breaks_under_10_rule = breaks_under_10_rule,

    C1_found = !is.na(idx_c1),
    C3_boundary_found = !is.na(idx_first_not_better_than_none),

    C1_cutoff = c1_cutoff,
    C1_sensitivity = c1_diag$sensitivity,
    C1_specificity = c1_diag$specificity,
    C1_PPV = c1_diag$ppv,
    C1_NPV = c1_diag$npv,
    C1_Treatment_Equality = c1_diag$treatment_equality,
    C1_Accuracy = c1_diag$accuracy,
    C1_Balanced_Accuracy = c1_diag$balanced_accuracy,

    C2_cutoff = c2_cutoff,
    C2_sensitivity = c2_diag$sensitivity,
    C2_specificity = c2_diag$specificity,
    C2_PPV = c2_diag$ppv,
    C2_NPV = c2_diag$npv,
    C2_Treatment_Equality = c2_diag$treatment_equality,
    C2_Accuracy = c2_diag$accuracy,
    C2_Balanced_Accuracy = c2_diag$balanced_accuracy,

    C2_Youden_J = c2_youden_j,

    C3_cutoff = c3_cutoff,
    C3_sensitivity = c3_diag$sensitivity,
    C3_specificity = c3_diag$specificity,
    C3_PPV = c3_diag$ppv,
    C3_NPV = c3_diag$npv,
    C3_Treatment_Equality = c3_diag$treatment_equality,
    C3_Accuracy = c3_diag$accuracy,
    C3_Balanced_Accuracy = c3_diag$balanced_accuracy,

    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ## Append custom-threshold columns, if requested
  if (length(custom_columns) > 0L) {

    custom_df <- as.data.frame(
      custom_columns,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    summary_row <- cbind(
      summary_row,
      custom_df
    )
  }

  ## Explicitly return the completed summary row
  summary_row
}

for (screening_date in screening_dates) {

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\n===== Screening date: %s =====\n",
      screening_date
    ))
    flush.console()
  }


  pop_data <- build_populations_for_screening_date(
    screening_date = screening_date,
    require_genetic_data = TRUE,
    genetic_data_marker = "ContiMultiethnicGRS"
  )

  missing_pops <- setdiff(selected_pops, names(pop_data))
  if (length(missing_pops) > 0) {
    stop(paste0("Unknown population(s): ", paste(missing_pops, collapse = ", ")))
  }

  

    for (pop in selected_pops) {
      d <- pop_data[[pop]]
      if (is.null(d) || nrow(d) == 0) {
        if (isTRUE(verbose)) {
          cat(sprintf(
            "Skipping empty population: %s at %s\n",
            pop,
            screening_date
          ))
        }
        next
      }

      for (current_outcome in outcomes) {
        for (pred in predictors) {

        if (!pred %in% colnames(d)) {
          if (isTRUE(verbose)) {
            cat(sprintf(
              "Skipping %s: predictor not found in %s at %s for %s\n",
              pred,
              pop,
              screening_date,
              current_outcome
            ))
          }
          next
        }

          for (cov in covariates_list) {
            effective_covariates <- resolve_covariates(cov, predictor = pred)
            cov_label <- covariate_label(effective_covariates)

            combo_count <- combo_count + 1L

            if (isTRUE(verbose)) {
              cat(sprintf(
                "[%d/%d] %s | %s | %s | %s | %s\n",
                combo_count,
                total_combos,
                screening_date,
                current_outcome,
                pop,
                pred,
                cov_label
              ))
              flush.console()
            }

            if (uses_age_asd_inappropriately(
              screening_date = screening_date,
              predictor = pred,
              covariates = effective_covariates
            )) {
              if (isTRUE(verbose)) {
                cat("  Skipping: ", age_asd_skip_message, "\n", sep = "")
              }
              next
            }

            model_spec_key <- paste(
              screening_date,
              current_outcome,
              pop,
              pred,
              cov_label,
              sep = "__"
            )

            missing_covariates <- setdiff(
              effective_covariates,
              colnames(d)
            )

            if (length(missing_covariates) > 0) {
              if (isTRUE(verbose)) {
                cat(sprintf(
                  "  Skipping: missing covariates: %s\n",
                  paste(missing_covariates, collapse = ", ")
                ))
              }
              next
            }

            if (exists(
              model_spec_key,
              envir = seen_model_specs,
              inherits = FALSE
            )) {
              if (isTRUE(verbose)) {
                cat("  Skipping duplicate resolved model.\n")
              }
              next
            }

              assign(model_spec_key, TRUE, envir = seen_model_specs)

            tryCatch({

              model_fit <- run_logreg(
                data = d,
                outcome = current_outcome,
                predictor = pred,
                covariates = effective_covariates,
                screening_date = screening_date,
                subset_controls = subset_controls,
                plot_roc = FALSE,
                plot_pr = FALSE,
                show_output = FALSE
              )

              nb_results <- nb_curve(
                data = model_fit$data,
                outcome = current_outcome,
                cutoff_values = evaluated_cutoff_values,
                positive_level = positive_level,
                negative_level = negative_level,
                plot_curve = print_all_curves,
                
                plot_title = paste0(
                  "DCA: ", pop, " | ", pred,
                  "\nCovariates: ", cov_label,
                  "\nOutcome: ", current_outcome, " | Screen: ", screening_date
                )

              )

              if (isTRUE(summary_table)) {
                summary_key <- paste(screening_date, current_outcome, pop, pred, cov_label, sep = "__")

                summary_row <- summarise_nb_cutoffs(
                  nb_results = nb_results,
                  outcome = current_outcome,
                  screening_date = screening_date,
                  population = pop,
                  predictor = pred,
                  standard_thresholds = cutoff_values,
                  custom_thresholds = custom_thresholds,
                  suppress_counts = suppress_counts
                )
                summary_row$Covariates <- cov_label

                summary_row <- summary_row %>%
                dplyr::relocate(Covariates, .after = Predictor)
                
                summary_rows[[summary_key]] <- summary_row
              }

              key <- paste(screening_date, current_outcome, pop, pred, cov_label, sep = "__")

              if (isTRUE(retain_model_runs)) {
                results[[key]] <- list(
                  screening_date = screening_date,
                  outcome = current_outcome,
                  population = pop,
                  predictor = pred,
                  covariates = effective_covariates,
                  covariate_label = cov_label,
                  model_results = model_fit,
                  nb_results = nb_results
                )
              }

              if (isTRUE(verbose)) {
                cat("  Completed successfully.\n")
                flush.console()
              }

              if (!isTRUE(retain_model_runs)) {
                rm(model_fit, nb_results)
                invisible(gc(verbose = FALSE))
              }

            }, error = function(e) {
              warning(sprintf(
                "Skipping DCA: %s | %s | %s | %s | %s: %s",
                screening_date,
                current_outcome,
                pop,
                pred,
                cov_label,
                e$message
              ))
            })
          }
        }
      }
    }
  }

  if (isTRUE(summary_table)) {
    summary_df <- if (length(summary_rows) > 0) {
      do.call(rbind, summary_rows)
    } else {
      data.frame(stringsAsFactors = FALSE)
    }

    ## Write to csv
    write.csv(summary_df, file = output_csv, row.names = FALSE)

    message(sprintf("Summary table created! It is saved as '%s'.", output_csv))

    return(list(
      model_runs = results,
      summary_table = summary_df
    ))
  }

  return(results)
}


###############################################################################################################################################
# Simulate participant for a given model to derive estimated absolute risk, based on their GRS and any other covariates included in the model #
###############################################################################################################################################

simulate_participant_risk <- function(logreg_model,
                                      participant_GRS_percentile = 0.5, # default to median GRS
                                      participant_FHx_PrCa = 0, # default to no family history of prostate cancer
                                      participant_FHx_BrCa = 0, # default to no family history of breast cancer
                                      participant_FHx_PrCa_BrCa = 0, # default to no family history of prostate or breast cancer
                                      participant_age_asd = 50, # default to age 50 at artifical screening date (ASD)
                                      GRS_standardisation = "global") { # default to global GRS standardisation (standardised across all men). 
                                                                        # other option would be "model" which would standardise GRS within the model's population

  # warnings

  if (!is.null(logreg_model$subset_controls) && isTRUE(logreg_model$subset_controls)) {
  warning(
    "This model was fit after subsetting controls. Predicted probabilities may not represent absolute risk."
    )
  }

  if (participant_GRS_percentile < 0 || participant_GRS_percentile > 1) {
    stop("participant_GRS_percentile must be between 0 and 1.")
  }

  if (participant_age_asd < 20 || participant_age_asd > 100) {
    stop("participant_age_asd must be between 20 and 100.")
  }
  
  # Extract model and data from logreg_model input

  model <- logreg_model$model
  data <- logreg_model$data

  # Create a new data frame for the participant with the same structure as the model's data

  if (GRS_standardisation == "global") {
    # Use global standardisation: calculate the quantiles based on the entire dataset

    participant_data <- data.frame(
      ContiMultiethnicGRS = quantile(PrCa_prediagnosis_males$ContiMultiethnicGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      ContiEURGRS = quantile(PrCa_prediagnosis_males$ContiEURGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      ContiAFRGRS = quantile(PrCa_prediagnosis_males$ContiAFRGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      ContiEASGRS = quantile(PrCa_prediagnosis_males$ContiEASGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      ContiHISGRS = quantile(PrCa_prediagnosis_males$ContiHISGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangMultiethnicGRS = quantile(PrCa_prediagnosis_males$WangMultiethnicGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangEURGRS = quantile(PrCa_prediagnosis_males$WangEURGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangAFRGRS = quantile(PrCa_prediagnosis_males$WangAFRGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangEASGRS = quantile(PrCa_prediagnosis_males$WangEASGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangHISGRS = quantile(PrCa_prediagnosis_males$WangHISGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      FHx_PrCa = participant_FHx_PrCa,
      FHx_BrCa = participant_FHx_BrCa,
      FHx_PrCa_BrCa = participant_FHx_PrCa_BrCa,
      age_asd = participant_age_asd
    )
  
  } else if (GRS_standardisation == "model") {
    # Use model standardisation: calculate the quantiles based on the model's data

    participant_data <- data.frame(
      ContiMultiethnicGRS = quantile(data$ContiMultiethnicGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      ContiEURGRS = quantile(data$ContiEURGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      ContiAFRGRS = quantile(data$ContiAFRGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      ContiEASGRS = quantile(data$ContiEASGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      ContiHISGRS = quantile(data$ContiHISGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangMultiethnicGRS = quantile(data$WangMultiethnicGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangEURGRS = quantile(data$WangEURGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangAFRGRS = quantile(data$WangAFRGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangEASGRS = quantile(data$WangEASGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      WangHISGRS = quantile(data$WangHISGRS, probs = participant_GRS_percentile, na.rm = TRUE, names = FALSE),
      FHx_PrCa = participant_FHx_PrCa,
      FHx_BrCa = participant_FHx_BrCa,
      FHx_PrCa_BrCa = participant_FHx_PrCa_BrCa,
      age_asd = participant_age_asd
    )
  } else {
    stop("Invalid GRS_standardisation option. Choose 'global' or 'model'.")
  }
  
  # Ensure that the participant_data includes the same covariates as the model (doesn't matter if it has more, but all model covariates must be present)
  model_vars <- all.vars(delete.response(terms(model)))
  missing_vars <- setdiff(model_vars, names(participant_data))

  if (length(missing_vars) > 0) {
    stop(paste("Missing covariates in participant_data:", paste(missing_vars, collapse = ", ")))
  }

  # Predict the probability for the participant

  probability <- predict(model, newdata = participant_data, type = "response")

  return(probability)
}      


#################################################################################################################
# Statistical tests - Paired and Unpaired DeLong tests for comparing ROC AUCs between models and/or populations #
#################################################################################################################

## Paired DeLong test: compare two models in the same population, same screening date, same outcome, different predictor and/or covariates

paired_DeLong_test <- function(populations,
                               outcome,
                               predictor_a,
                               predictor_b,
                               covariates_a = NULL,
                               covariates_b = NULL,
                               screening_date,
                               subset_controls = FALSE,
                               id_col = "pid",
                               min_n = 20,
                               adjust_method = "BH") {

  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required for paired DeLong tests.")
  }

  if (missing(predictor_b) || is.null(predictor_b)) {
    stop("predictor_b must be supplied.")
  }

  if (missing(screening_date) || is.null(screening_date)) {
    stop("screening_date must be supplied.")
  }

  if (isTRUE(subset_controls)) {
    warning(
      "subset_controls = TRUE may undermine paired comparison because controls may be sampled differently between model A and model B. ",
      "For paired DeLong tests, subset_controls = FALSE is recommended."
    )
  }

  model_a_label <- if (is.null(covariates_a) || length(covariates_a) == 0) {
    predictor_a
  } else {
    paste(predictor_a, "+", paste(covariates_a, collapse = " + "))
  }

  model_b_label <- if (is.null(covariates_b) || length(covariates_b) == 0) {
    predictor_b
  } else {
    paste(predictor_b, "+", paste(covariates_b, collapse = " + "))
  }

  out <- lapply(names(populations), function(pop_name) {

    dat <- populations[[pop_name]]

    model_a <- run_logreg(
      data = dat,
      outcome = outcome,
      predictor = predictor_a,
      covariates = covariates_a,
      screening_date = screening_date,
      subset_controls = subset_controls,
      plot_roc = FALSE,
      show_output = FALSE
    )

    model_b <- run_logreg(
      data = dat,
      outcome = outcome,
      predictor = predictor_b,
      covariates = covariates_b,
      screening_date = screening_date,
      subset_controls = subset_controls,
      plot_roc = FALSE,
      show_output = FALSE
    )

    df_a <- model_a$data
    df_b <- model_b$data

    use_id_col <- !is.null(id_col) &&
      id_col %in% names(df_a) &&
      id_col %in% names(df_b)

    if (use_id_col) {
      common <- dplyr::inner_join(
        df_a %>%
          dplyr::select(
            join_id = dplyr::all_of(id_col),
            outcome_value = dplyr::all_of(outcome),
            pred_a = pred
          ),
        df_b %>%
          dplyr::select(
            join_id = dplyr::all_of(id_col),
            outcome_check = dplyr::all_of(outcome),
            pred_b = pred
          ),
        by = "join_id"
      )
    } else {
      warning(
        "id_col not found in both model datasets for population '",
        pop_name,
        "'. Falling back to row order. This is less robust."
      )

      df_a <- df_a %>% dplyr::mutate(join_id = dplyr::row_number())
      df_b <- df_b %>% dplyr::mutate(join_id = dplyr::row_number())

      common <- dplyr::inner_join(
        df_a %>%
          dplyr::select(
            join_id,
            outcome_value = dplyr::all_of(outcome),
            pred_a = pred
          ),
        df_b %>%
          dplyr::select(
            join_id,
            outcome_check = dplyr::all_of(outcome),
            pred_b = pred
          ),
        by = "join_id"
      )
    }

    common <- common %>%
      dplyr::filter(
        !is.na(outcome_value),
        !is.na(outcome_check),
        !is.na(pred_a),
        !is.na(pred_b)
      )

    if (!all(common$outcome_value == common$outcome_check, na.rm = TRUE)) {
      stop(
        "Outcome values differ between model A and model B for population '",
        pop_name,
        "'. This function is only for same-population, same-screening-date, same-outcome paired comparisons."
      )
    }

    n_cases <- sum(common$outcome_value == 1, na.rm = TRUE)
    n_controls <- sum(common$outcome_value == 0, na.rm = TRUE)

    if (
      nrow(common) < min_n ||
      n_cases == 0 ||
      n_controls == 0
    ) {
      return(data.frame(
        Population = pop_name,
        Outcome = outcome,
        Screening_Date = screening_date,
        Model_A = model_a_label,
        Model_B = model_b_label,
        AUC_A = NA_real_,
        AUC_B = NA_real_,
        Delta_AUC = NA_real_,
        p_value = NA_real_,
        N_usable = nrow(common),
        N_Cases = n_cases,
        N_Controls = n_controls,
        stringsAsFactors = FALSE
      ))
    }

    roc_a <- pROC::roc(
      response = common$outcome_value,
      predictor = common$pred_a,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )

    roc_b <- pROC::roc(
      response = common$outcome_value,
      predictor = common$pred_b,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )

    test_obj <- tryCatch(
      pROC::roc.test(
        roc_a,
        roc_b,
        paired = TRUE,
        method = "delong"
      ),
      error = function(e) {
        warning(paste0(
          "Paired DeLong test failed for population '",
          pop_name,
          "': ",
          e$message
        ))
        list(p.value = NA_real_)
      }
    )

    data.frame(
      Population = pop_name,
      Outcome = outcome,
      Screening_Date = screening_date,
      Model_A = model_a_label,
      Model_B = model_b_label,
      AUC_A = as.numeric(pROC::auc(roc_a)),
      AUC_B = as.numeric(pROC::auc(roc_b)),
      Delta_AUC = as.numeric(pROC::auc(roc_a)) - as.numeric(pROC::auc(roc_b)),
      p_value = as.numeric(test_obj$p.value),
      N_usable = nrow(common),
      N_Cases = n_cases,
      N_Controls = n_controls,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(out) %>%
    dplyr::mutate(p_value_adj = p.adjust(p_value, method = adjust_method))
}


## Unpaired DeLong's test function to compare AUCs of two logistic regression models on different populations

unpaired_DeLong_test <- function(populations,
                                 outcome,
                                 predictor,
                                 covariates = NULL,
                                 screening_date = NULL,
                                 subset_controls = FALSE,
                                 control_group = NULL,
                                 pairings = NULL,
                                 min_n = 20,
                                 adjust_method = "BH") {

  model_label <- if (is.null(covariates) || length(covariates) == 0) {
    predictor
  } else {
    paste(predictor, "+", paste(covariates, collapse = " + "))
  }

  pred_by_pop <- lapply(names(populations), function(pop_name) {
    dat <- populations[[pop_name]]

    model_obj <- run_logreg(
      data = dat,
      outcome = outcome,
      predictor = predictor,
      covariates = covariates,
      screening_date = screening_date,
      subset_controls = subset_controls,
      plot_roc = FALSE,
      show_output = FALSE
    )

    model_obj$data %>%
      dplyr::select(outcome_value = dplyr::all_of(outcome), pred = pred) %>%
      dplyr::filter(!is.na(outcome_value), !is.na(pred)) %>%
      dplyr::mutate(Population = pop_name)
  }) %>%
    dplyr::bind_rows()

  pop_names <- names(populations)

  if (!is.null(control_group)) {
    if (!(control_group %in% pop_names)) {
      stop("control_group must match one of names(populations)")
    }
    other_pops <- setdiff(pop_names, control_group)
    pairings <- lapply(other_pops, function(pop_name) c(control_group, pop_name))
  } else if (is.null(pairings)) {
    pairings <- combn(pop_names, 2, simplify = FALSE)
  }

  out <- lapply(pairings, function(pair_vec) {
  pop_a <- pair_vec[1]
  pop_b <- pair_vec[2]

  dat_a <- pred_by_pop %>% dplyr::filter(Population == pop_a)
  dat_b <- pred_by_pop %>% dplyr::filter(Population == pop_b)

  if (nrow(dat_a) < min_n || nrow(dat_b) < min_n ||
      length(unique(dat_a$outcome_value)) < 2 ||
      length(unique(dat_b$outcome_value)) < 2) {
    return(data.frame(
      Population_A = pop_a,
      Population_B = pop_b,
      Outcome = outcome,
      Screening_Date = screening_date,
      Model = model_label,
      AUC_A = NA_real_,
      AUC_B = NA_real_,
      Delta_AUC = NA_real_,
      p_value = NA_real_,
      N_A = nrow(dat_a),
      N_B = nrow(dat_b),
      stringsAsFactors = FALSE
    ))
  }

  roc_a <- pROC::roc(
    response = dat_a$outcome_value,
    predictor = dat_a$pred,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  roc_b <- pROC::roc(
    response = dat_b$outcome_value,
    predictor = dat_b$pred,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  test_obj <- tryCatch(
    pROC::roc.test(roc_a, roc_b, paired = FALSE, method = "delong"),
    error = function(e) {
      warning(paste0(
        "DeLong unpaired test failed for populations '",
        pop_a, "' vs '", pop_b, "': ", e$message
      ))
      list(p.value = NA_real_)
    }
  )

  data.frame(
    Population_A = pop_a,
    Population_B = pop_b,
    Outcome = outcome,
    Screening_Date = screening_date,
    Model = model_label,
    AUC_A = as.numeric(pROC::auc(roc_a)),
    AUC_B = as.numeric(pROC::auc(roc_b)),
    Delta_AUC = as.numeric(pROC::auc(roc_a)) - as.numeric(pROC::auc(roc_b)),
    p_value = as.numeric(test_obj$p.value),
    N_A = nrow(dat_a),
    N_B = nrow(dat_b),
    stringsAsFactors = FALSE
  )
})




  dplyr::bind_rows(out) %>%
    dplyr::mutate(p_value_adj = p.adjust(p_value, method = adjust_method))

}


## DeLong's test for comparing AUCs of Logistic Regression models in the same population on different screening dates (unpaired test)

screening_date_DeLong_test <- function(populations,
                                       populations_b,
                                       outcome,
                                       predictor,
                                       covariates = NULL,
                                       screening_date_a,
                                       screening_date_b,
                                       subset_controls = FALSE,
                                       id_col = "pid",
                                       min_n = 20,
                                       adjust_method = "BH") {

  common_pop_names <- intersect(names(populations), names(populations_b))

  if (length(common_pop_names) == 0) {
    stop("No matching population names found between 'populations' and 'populations_b'.")
  }

  model_label <- if (is.null(covariates) || length(covariates) == 0) {
    predictor
  } else {
    paste(predictor, "+", paste(covariates, collapse = " + "))
  }

  safe_roc <- function(df, outcome, pred_col = "pred") {
    valid <- !is.na(df[[outcome]]) & !is.na(df[[pred_col]])
    d <- df[valid, , drop = FALSE]

    n_cases <- sum(d[[outcome]] == 1, na.rm = TRUE)
    n_controls <- sum(d[[outcome]] == 0, na.rm = TRUE)

    if (
      nrow(d) < min_n ||
      n_cases == 0 ||
      n_controls == 0
    ) {
      return(list(
        roc = NULL,
        auc = NA_real_,
        n = nrow(d),
        n_cases = n_cases,
        n_controls = n_controls
      ))
    }

    roc_obj <- pROC::roc(
      response = d[[outcome]],
      predictor = d[[pred_col]],
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )

    list(
      roc = roc_obj,
      auc = as.numeric(pROC::auc(roc_obj)),
      n = nrow(d),
      n_cases = n_cases,
      n_controls = n_controls
    )
  }

  out <- lapply(common_pop_names, function(pop_name) {
    dat_a <- populations[[pop_name]]
    dat_b <- populations_b[[pop_name]]

    model_a <- run_logreg(
      data = dat_a,
      outcome = outcome,
      predictor = predictor,
      covariates = covariates,
      screening_date = screening_date_a,
      subset_controls = subset_controls,
      plot_roc = FALSE,
      show_output = FALSE
    )

    model_b <- run_logreg(
      data = dat_b,
      outcome = outcome,
      predictor = predictor,
      covariates = covariates,
      screening_date = screening_date_b,
      subset_controls = subset_controls,
      plot_roc = FALSE,
      show_output = FALSE
    )

    roc_a_info <- safe_roc(model_a$data, outcome)
    roc_b_info <- safe_roc(model_b$data, outcome)

    p_val <- NA_real_

    if (!is.null(roc_a_info$roc) && !is.null(roc_b_info$roc)) {
      test_obj <- tryCatch(
        pROC::roc.test(
          roc_a_info$roc,
          roc_b_info$roc,
          paired = FALSE,
          method = "delong"
        ),
        error = function(e) {
          warning(paste0(
            "DeLong test failed for population '",
            pop_name,
            "': ",
            e$message
          ))
          list(p.value = NA_real_)
        }
      )

      p_val <- as.numeric(test_obj$p.value)
    }

    common_n <- NA_integer_
    common_outcome_a_cases <- NA_integer_
    common_outcome_b_cases <- NA_integer_
    common_outcome_agreement <- NA_real_

    if (!is.null(id_col) &&
        id_col %in% names(model_a$data) &&
        id_col %in% names(model_b$data)) {

      common <- dplyr::inner_join(
        model_a$data %>%
          dplyr::select(
            join_id = dplyr::all_of(id_col),
            outcome_a = dplyr::all_of(outcome),
            pred_a = pred
          ),
        model_b$data %>%
          dplyr::select(
            join_id = dplyr::all_of(id_col),
            outcome_b = dplyr::all_of(outcome),
            pred_b = pred
          ),
        by = "join_id"
      ) %>%
        dplyr::filter(
          !is.na(outcome_a),
          !is.na(outcome_b),
          !is.na(pred_a),
          !is.na(pred_b)
        )

      common_n <- nrow(common)
      common_outcome_a_cases <- sum(common$outcome_a == 1, na.rm = TRUE)
      common_outcome_b_cases <- sum(common$outcome_b == 1, na.rm = TRUE)
      common_outcome_agreement <- mean(common$outcome_a == common$outcome_b, na.rm = TRUE)
    }

    data.frame(
      Population = pop_name,
      Outcome = outcome,
      Model = model_label,
      Screening_Date_A = screening_date_a,
      Screening_Date_B = screening_date_b,

      AUC_A = roc_a_info$auc,
      AUC_B = roc_b_info$auc,
      Delta_AUC = roc_a_info$auc - roc_b_info$auc,
      p_value = p_val,

      N_A = roc_a_info$n,
      N_A_Cases = roc_a_info$n_cases,
      N_A_Controls = roc_a_info$n_controls,

      N_B = roc_b_info$n,
      N_B_Cases = roc_b_info$n_cases,
      N_B_Controls = roc_b_info$n_controls,

      #N_Common = common_n,
      #Common_Outcome_A_Cases = common_outcome_a_cases,
      #Common_Outcome_B_Cases = common_outcome_b_cases,       # Commenting out because of disclosure concerns with small counts
      #Common_Outcome_Agreement = common_outcome_agreement,

      Test_Type = "Unpaired DeLong on date-specific outcomes",
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(out) %>%
    dplyr::mutate(p_value_adj = p.adjust(p_value, method = adjust_method))
}

## Plot principal components vs. residuals from a fitted model

plot_pcs_vs_residuals <- function(model,
                                  pcs = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), 
                                  outcome = NULL,
                                  colour_by = "ethnicity_group_narrow"
                                ) {

  if (is.null(model$data)) {
    stop("model must be the result returned by run_logreg(), with participant data in model$data.")
  }
  data <- model$data

  if (!colour_by %in% names(data)) {
    stop(sprintf("Colour column '%s' was not found in model$data.", colour_by))
  }

  if (!"residual" %in% names(data)) {
    if (is.null(outcome) || !outcome %in% names(data) || !"pred" %in% names(data)) {
      stop("model$data must contain 'residual', or both the outcome and 'pred' columns.")
    }
    data$residual <- data[[outcome]] - data$pred
  }

  groups <- factor(data[[colour_by]])
  group_colours <- stats::setNames(
    grDevices::hcl.colors(nlevels(groups), palette = "Dark 3"),
    levels(groups)
  )

  for (pc in pcs) {
    if (pc %in% names(data)) {
      keep <- stats::complete.cases(data[[pc]], data$residual, groups)
      plot_groups <- droplevels(groups[keep])
      plot_colours <- unname(group_colours[as.character(plot_groups)])

      plot(data[[pc]][keep], data$residual[keep],
           xlab = pc,
           ylab = "Residuals",
           main = paste("Residuals vs", pc),
           col = plot_colours,
           pch = 16)
      legend("topright",
             legend = levels(plot_groups),
             col = unname(group_colours[levels(plot_groups)]),
             pch = 16,
             title = colour_by,
             bty = "n")
    }
  }
}