from pathlib import Path
path = Path(r"c:\Users\cw1245\Documents\PhD\Callum_Scripts\LUTS.R")
text = path.read_text(encoding='utf-8')
lines = text.splitlines(keepends=True)
marker = "# Step 9 - Model Logistic Regression, Confusion Matrix, and OR/RR tables #\n"
idx = None
for i, line in enumerate(lines):
    if line == marker:
        idx = i
        break
if idx is None:
    raise SystemExit("Marker not found")
insert = '''# NEW FUNCTION: logreg_table() - Generate summary table for all logistic regressions
# (Population x Outcome x GRS x Covariates)
logreg_table <- function(
    populations = NULL,
    grs_list = NULL,
    outcomes = c("PrCa", "PrCa_2yrs", "PrCa_5yrs", "PrCa_10yrs",
                 "PrCa_actionable", "PrCa_actionable_2yrs", "PrCa_actionable_5yrs", "PrCa_actionable_10yrs",
                 "PrCa_severe", "PrCa_severe_2yrs", "PrCa_severe_5yrs", "PrCa_severe_10yrs"),
    covariates_list = list(NULL, "event_age"),  # event_age is the relevant covariate here
    plot_roc = FALSE,
    verbose = TRUE
) {
  
  # Default populations if not specified
  if (is.null(populations)) {
    populations <- list(
      "All" = PCa_iv_covariates_GRS_predhorizon,
      "White" = PCa_iv_covariates_GRS_predhorizon_WhiteOnly,
      "Black" = PCa_iv_covariates_GRS_predhorizon_BlackOnly,
      "Mixed" = PCa_iv_covariates_GRS_predhorizon_Mixed,
      "Black+Mixed" = PCa_iv_covariates_GRS_predhorizon_BlackMixed,
      "EUR" = PCa_iv_covariates_GRS_predhorizon_EUROnly,
      "AFR" = PCa_iv_covariates_GRS_predhorizon_AFROnly,
      "EAS" = PCa_iv_covariates_GRS_predhorizon_EASOnly,
      "CSA" = PCa_iv_covariates_GRS_predhorizon_CSAOnly,
      "MID" = PCa_iv_covariates_GRS_predhorizon_MIDOnly,
      "AMR" = PCa_iv_covariates_GRS_predhorizon_AMROnly
    )
  }
  
  # Default GRS list if not specified
  if (is.null(grs_list)) {
    grs_list <- c(
      "multiethnicGRS",
      "EuropeanGRS",
      "AfricanGRS",
      "East_AsianGRS",
      "HispanicGRS",
      "adjustedGRS",
      "WangmultiethnicGRS",
      "WangEuropeanGRS",
      "WangAfricanGRS",
      "WangEast_AsianGRS",
      "WangHispanicGRS"
    )
  }
  
  # Initialize results dataframe
  results <- data.frame(
    Outcome = character(),
    Population = character(),
    GRS = character(),
    Covariates = character(),
    N_Cases = integer(),
    N_Controls = integer(),
    ROC_AUC = numeric(),
    ROC_AUC_CI_Lower = numeric(),
    ROC_AUC_CI_Upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Total combinations to process
  total_combos <- length(outcomes) * length(populations) * length(grs_list) * length(covariates_list)
  combo_count <- 0
  
  # Loop over all combinations
  for (current_outcome in outcomes) {
    for (pop_name in names(populations)) {
      pop_data <- populations[[pop_name]]
      
      for (grs_pred in grs_list) {
        
        # Check if GRS exists in data
        if (!grs_pred %in% colnames(pop_data)) {
          if (verbose) cat(sprintf("Skipping %s (not in %s)\n", grs_pred, pop_name))
          next
        }
        
        for (cov in covariates_list) {
          combo_count <- combo_count + 1
          
          # Progress indicator
          if (verbose) {
            cov_label <- ifelse(is.null(cov), "None", cov)
            cat(sprintf("[%d/%d] %s | %s | %s | %s\n", 
                        combo_count, total_combos, current_outcome, pop_name, grs_pred, cov_label))
          }
          
          # Run logistic regression
          tryCatch({
            model <- run_logreg(
              data = pop_data,
              outcome = current_outcome,
              predictor = grs_pred,
              covariates = cov,
              plot_roc = plot_roc
            )
            
            # Extract ROC AUC and CI
            roc_obj <- model$roc
            auc_val <- as.numeric(roc_obj$auc)
            auc_ci <- as.numeric(roc_obj$ci)  # Returns [lower, AUC, upper]
            auc_lower <- auc_ci[1]
            auc_upper <- auc_ci[3]
            
            # Count cases and controls in the dataset used
            n_cases <- sum(model$data[[current_outcome]] == 1, na.rm = TRUE)
            n_controls <- sum(model$data[[current_outcome]] == 0, na.rm = TRUE)
            
            # Covariate label
            cov_label <- ifelse(is.null(cov), "None", cov)
            
            # Append to results
            results <- rbind(results, data.frame(
              Outcome = current_outcome,
              Population = pop_name,
              GRS = grs_pred,
              Covariates = cov_label,
              N_Cases = n_cases,
              N_Controls = n_controls,
              ROC_AUC = auc_val,
              ROC_AUC_CI_Lower = auc_lower,
              ROC_AUC_CI_Upper = auc_upper,
              stringsAsFactors = FALSE
            ))
            
          }, error = function(e) {
            if (verbose) cat(sprintf("  ERROR: %s\n", e$message))
          })
        }
      }
    }
  }
  
  return(results)
}


'''

lines.insert(idx, insert)
path.write_text(''.join(lines), encoding='utf-8')
print('Inserted logreg_table at line', idx+1)
