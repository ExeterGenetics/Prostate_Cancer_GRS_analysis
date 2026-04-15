## This script sets up the packages and functions used in this repository's scripts

# Validate rlang before loading broader package stack.
if (!requireNamespace("rlang", quietly = TRUE)) {
  install.packages("rlang", dependencies = TRUE)
}
if (packageVersion("rlang") < "1.1.6") {
  stop("rlang >= 1.1.6 is required. Restart R, run install.packages('rlang'), then knit again.")
}

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

packages_needed <- c(
  "pROC",
  "epiR",
  "RMySQL",
  "readstata13",
  "survminer",
  "tidyverse",
  "DiagrammeR",
  "extrafont",
  "showtext",
  "readxl",
  "remotes",
  "caret",
  "bigsnpr",
  "kableExtra",
  "PRROC",
  "randomForest"
)

installed_pkgs <- rownames(installed.packages())

new_pkgs <- packages_needed[!packages_needed %in% installed_pkgs]
if (length(new_pkgs) > 0) {
  # Install PRROC without dependencies to prevent it triggering a rlang
  # reinstall attempt that fails because rlang is already loaded.
  if ("PRROC" %in% new_pkgs) {
    install.packages("PRROC", dependencies = FALSE)
    new_pkgs <- setdiff(new_pkgs, "PRROC")
  }
  if (length(new_pkgs) > 0) {
    install.packages(new_pkgs, dependencies = TRUE)
  }
}

# After installation, check whether rlang (or any other critical package) was
# upgraded on disk to a version newer than what is currently loaded in this
# R session.  When that happens, calling library() on any package that requires
# the newer version triggers an unloadNamespace() attempt on rlang, which fails
# because almost everything imports it.  Detect this early and stop with a
# clear, actionable message instead of a confusing backtrace.
.rlang_loaded_ver <- tryCatch(
  numeric_version(getNamespaceVersion("rlang")),
  error = function(e) NULL
)
.rlang_disk_ver <- tryCatch(packageVersion("rlang"), error = function(e) NULL)
if (!is.null(.rlang_loaded_ver) && !is.null(.rlang_disk_ver) &&
    .rlang_disk_ver > .rlang_loaded_ver) {
  stop(
    "rlang was upgraded from ", .rlang_loaded_ver, " to ", .rlang_disk_ver,
    " during package installation.\n",
    "All packages are now installed — please restart R and knit again.\n",
    "This will not happen on the next knit."
  )
}
rm(.rlang_loaded_ver, .rlang_disk_ver)

invisible(lapply(packages_needed, library, character.only = TRUE))

# Helper: randomly subset controls to match the number of cases for a given outcome column
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

# Function to run a logistic regression and compute a ROC AUC curve with 95% CIs

run_logreg <- function(data,
                       outcome,
                       predictor,
                       covariates = NULL,
                       subset_controls = FALSE,
                       plot_roc = TRUE,
                       plot_pr = FALSE,
                       show_output = TRUE) {

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
  
  logreg <- glm(formula, data = data, family = binomial)
  if (show_output) {
    print(summary(logreg))
  }
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
  
  if (plot_roc && show_output) {
    roc_obj <- roc(data[[outcome]] ~ data$pred,
                   plot = TRUE,
                   print.auc = TRUE,
                   ci = TRUE)
  } else {
    roc_obj <- roc(data[[outcome]] ~ data$pred,
                   ci = TRUE)
  }
  
  if (show_output) {
    print(roc_obj)
  }
  
  pr_obj <- NULL
  if (plot_pr) {
    labels <- as.numeric(data[[outcome]] == 1)
    scores <- data$pred
    keep <- !is.na(labels) & !is.na(scores)
    labels_clean <- labels[keep]
    scores_clean <- scores[keep]
    prevalence <- mean(labels_clean == 1)

    prroc_available <- tryCatch(
      requireNamespace("PRROC", quietly = TRUE),
      error = function(e) FALSE
    )

    if (!prroc_available) {
      warning("PR curve skipped: PRROC is unavailable or failed to load in this session.")
    } else {
      pr_obj <- tryCatch(
        PRROC::pr.curve(
          scores.class0 = scores_clean[labels_clean == 1],
          scores.class1 = scores_clean[labels_clean == 0],
          curve = TRUE
        ),
        error = function(e) {
          warning(paste0("PR curve skipped due to PRROC error: ", e$message))
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
  
  return(list(
    model = logreg,
    data = data,
    outcome = outcome,
    predictor = predictor,
    covariates = covariates,
    roc = roc_obj,
    pr = pr_obj
  ))
}

# Function to run a confusion matrix

confusion_matrix <- function(data,
                             outcome,
                             cutoff_value = 0.5,
                             positive_level = 1,
                             negative_level = 0,
                             print_epi = TRUE) {

  if (cutoff_value <= 0 || cutoff_value >= 1) {
    stop("cutoff_value must be between 0 and 1 (exclusive).")
  }
  if (!("pred" %in% names(data))) {
    stop("data must contain a column named 'pred' with predicted probabilities.")
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

  print(cm)

  epi <- epi.tests(cm)
  if (print_epi) print(epi)

  # ---- Net benefit calculations ----
  get_cell <- function(mat, r, c) {
    r <- as.character(r); c <- as.character(c)
    if (r %in% rownames(mat) && c %in% colnames(mat)) mat[r, c] else 0
  }

  TP <- get_cell(cm, TRUE,  positive_level)
  FP <- get_cell(cm, TRUE,  negative_level)
  N  <- sum(cm)

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

  cat("\nThreshold (pt) =", pt,
      "\nNet benefit (model)     =", net_benefit,
      "\nNet benefit (treat all) =", nb_treat_all,
      "\nNet benefit (treat none)=", nb_treat_none,
      "\n")

  return(list(
    confusion_matrix = cm,
    epi_tests = epi,
    net_benefit = net_benefit,
    net_benefit_treat_all = nb_treat_all,
    net_benefit_treat_none = nb_treat_none,
    components = list(TP = TP, FP = FP, N = N, pt = pt, w = w, prevalence = prevalence),
    data = data
  ))
}

## Function to run an Odds Ratio table

ORtable <- function(data,
                    outcome,
                    group_col = "Group",
                    positive_level = 1,
                    pred_var = "pred",
                    bins = c(10, 20, 30, 40, 60, 70, 80, 90, 100),
                    use_existing_cols = FALSE,
                    digits = 2) {
  # Basic checks
  stopifnot(outcome %in% names(data))
  stopifnot(group_col %in% names(data))
  if (!pred_var %in% names(data)) {
    stop("`pred_var` not found in `data`.")
  }
  if (use_existing_cols) {
    message("`use_existing_cols` is ignored in this ORtable version; percentile bands are computed from `pred_var`.")
  }
  
  # Make Group a factor (preserve order if already set)
  data[[group_col]] <- as.factor(data[[group_col]])
  
  # Percentile cut points for pred_var
  q <- stats::quantile(data[[pred_var]], probs = seq(0, 1, 0.1), na.rm = TRUE, names = FALSE)
  
  # Define requested percentile bands and the reference band
  band_defs <- list(
    c(0.0, 0.1),
    c(0.1, 0.2),
    c(0.2, 0.3),
    c(0.3, 0.4),
    c(0.4, 0.6),
    c(0.6, 0.7),
    c(0.7, 0.8),
    c(0.8, 0.9),
    c(0.9, 1.0)
  )
  band_labels <- c(
    "0-10th percentile",
    "10-20th percentile",
    "20-30th percentile",
    "30-40th percentile",
    "40-60th percentile",
    "60-70th percentile",
    "70-80th percentile",
    "80-90th percentile",
    "90-100th percentile"
  )
  ref_label <- "40-60th percentile"
  
  # Build a membership vector for each percentile band
  make_band_vector <- function(lo_prob, hi_prob) {
    lo_idx <- as.integer(lo_prob * 10) + 1
    hi_idx <- as.integer(hi_prob * 10) + 1
    lo_val <- q[lo_idx]
    hi_val <- q[hi_idx]
    x <- data[[pred_var]]
    if (hi_prob < 1.0) {
      return(x >= lo_val & x < hi_val)
    }
    return(x >= lo_val & x <= hi_val)
  }
  band_vectors <- lapply(band_defs, function(b) make_band_vector(b[1], b[2]))
  names(band_vectors) <- band_labels
  ref_vec <- band_vectors[[ref_label]]
  
  # Containers
  groups <- levels(data[[group_col]])
  
  # Long results accumulator
  res_long <- list()
  
  # Iterate over requested percentile bands and groups
  for (i in seq_along(band_labels)) {
    bin_label <- band_labels[i]
    inbin_vec <- band_vectors[[bin_label]]
    
    for (g in groups) {
      # Subset to group g and non-missing outcomes/predictions
      keep <- !is.na(data[[group_col]]) &
        !is.na(data[[outcome]]) &
        !is.na(data[[pred_var]]) &
        (data[[group_col]] == g)
      
      if (!any(keep)) {
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
      refbin <- ref_vec[keep]
      y_pos <- y == positive_level
      
      # Current band counts
      a <- sum(inbin & y_pos, na.rm = TRUE)
      b <- sum(inbin & !y_pos, na.rm = TRUE)
      n_in <- a + b
      
      # Reference (40-60th percentile) counts
      c <- sum(refbin & y_pos, na.rm = TRUE)
      d <- sum(refbin & !y_pos, na.rm = TRUE)
      n_ref <- c + d
      
      # Force reference band OR to 1 by definition
      if (identical(bin_label, ref_label)) {
        res_long[[length(res_long) + 1]] <- data.frame(
          bin = bin_label,
          group = g,
          n_in = n_in,
          events_in = a,
          rate_in = ifelse(n_in > 0, a / n_in, NA_real_),
          n_out = n_ref,
          events_out = c,
          rate_out = ifelse(n_ref > 0, c / n_ref, NA_real_),
          OR_pos = 1,
          CI_low = 1,
          CI_high = 1,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }
      
      # OR undefined if either current band or reference band has no observations
      if (n_in == 0L || n_ref == 0L) {
        res_long[[length(res_long) + 1]] <- data.frame(
          bin = bin_label,
          group = g,
          n_in = n_in,
          events_in = a,
          rate_in = ifelse(n_in > 0, a / n_in, NA_real_),
          n_out = n_ref,
          events_out = c,
          rate_out = ifelse(n_ref > 0, c / n_ref, NA_real_),
          OR_pos = NA_real_,
          CI_low = NA_real_,
          CI_high = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }
      
      # Haldane-Anscombe correction if any zero cells
      zero_cell <- any(c(a, b, c, d) == 0L)
      a2 <- if (zero_cell) a + 0.5 else a
      b2 <- if (zero_cell) b + 0.5 else b
      c2 <- if (zero_cell) c + 0.5 else c
      d2 <- if (zero_cell) d + 0.5 else d
      
      # OR (positive outcome): odds in current band vs odds in reference band
      OR <- (a2 * d2) / (b2 * c2)
      
      # Wald CI on log OR
      se_logOR <- sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
      CI_low <- exp(log(OR) - 1.96 * se_logOR)
      CI_high <- exp(log(OR) + 1.96 * se_logOR)
      
      # Fisher's exact p-value for current band vs reference band
      pval <- tryCatch(
        stats::fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
        error = function(e) NA_real_
      )
      
      res_long[[length(res_long) + 1]] <- data.frame(
        bin = bin_label, group = g,
        n_in = n_in, events_in = a, rate_in = a / n_in,
        n_out = n_ref, events_out = c, rate_out = c / n_ref,
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
  bins_order <- band_labels
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
    " in each percentile band of `", pred_var, "` against the 40-60th percentile reference band. ",
    "Reference band OR is fixed at 1 [1, 1]. 95% CI via Wald on log-OR; p-value via Fisher's exact test. ",
    "Haldane-Anscombe +0.5 applied when a zero cell occurs."
  )
  
  list(
    long = res_long,          # one row per (bin, group) with counts, rates, OR, CI, p
    wide_formatted = wide_format,  # display table: OR [L, U]; p=
    wide_OR = wide_OR         # numeric ORs only
  )
}

## Function to run a Risk Ratio table 

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
      
      # Haldane-Anscombe if any zero cell
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
    "Haldane-Anscombe +0.5 applied if any zero cell occurs.",
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


## Function to calculate the Net Reclassification Improvement (NRI)
## Compares a new prediction model against a reference model.
## Method: continuous NRI (Pencina et al. 2008, Stat Med).
##   - Cases:    NRI_cases    = P(up | event)   - P(down | event)
##   - Controls: NRI_controls = P(down | non-event) - P(up | non-event)
##   - Overall:  NRI          = NRI_cases + NRI_controls
## SE from Pencina 2008; 95% CIs and z-test p-values from normal approximation.

nri <- function(data,
                outcome,
                pred_new = "pred2",
                pred_old = "pred",
                positive_level = 1,
                digits = 3,
                show_output = TRUE) {
  
  stopifnot(outcome  %in% names(data))
  stopifnot(pred_new %in% names(data))
  stopifnot(pred_old %in% names(data))
  
  # Drop rows with any missing values in the three key columns
  complete <- !is.na(data[[outcome]]) &
    !is.na(data[[pred_new]]) &
    !is.na(data[[pred_old]])
  d <- data[complete, ]
  
  is_case <- d[[outcome]] == positive_level
  
  delta <- d[[pred_new]] - d[[pred_old]]
  
  # -- Cases ------------------------------------------------------------------
  n_cases         <- sum(is_case)
  n_up_cases      <- sum(delta[is_case]  > 0)
  n_down_cases    <- sum(delta[is_case]  < 0)
  p_up_cases      <- n_up_cases   / n_cases
  p_down_cases    <- n_down_cases / n_cases
  NRI_cases       <- p_up_cases - p_down_cases
  SE_NRI_cases    <- sqrt((p_up_cases + p_down_cases -
                             (p_up_cases - p_down_cases)^2) / n_cases)
  CI_low_cases    <- NRI_cases - 1.96 * SE_NRI_cases
  CI_high_cases   <- NRI_cases + 1.96 * SE_NRI_cases
  z_cases         <- NRI_cases / SE_NRI_cases
  p_cases         <- 2 * pnorm(-abs(z_cases))
  
  # -- Controls ---------------------------------------------------------------
  n_controls      <- sum(!is_case)
  n_up_controls   <- sum(delta[!is_case] > 0)
  n_down_controls <- sum(delta[!is_case] < 0)
  p_up_controls   <- n_up_controls   / n_controls
  p_down_controls <- n_down_controls / n_controls
  NRI_controls    <- p_down_controls - p_up_controls   # down is good for controls
  SE_NRI_controls <- sqrt((p_up_controls + p_down_controls -
                             (p_up_controls - p_down_controls)^2) / n_controls)
  CI_low_controls  <- NRI_controls - 1.96 * SE_NRI_controls
  CI_high_controls <- NRI_controls + 1.96 * SE_NRI_controls
  z_controls       <- NRI_controls / SE_NRI_controls
  p_controls       <- 2 * pnorm(-abs(z_controls))
  
  # -- Overall NRI ------------------------------------------------------------
  NRI             <- NRI_cases + NRI_controls
  SE_NRI          <- sqrt(SE_NRI_cases^2 + SE_NRI_controls^2)
  CI_low_NRI      <- NRI - 1.96 * SE_NRI
  CI_high_NRI     <- NRI + 1.96 * SE_NRI
  z_NRI           <- NRI / SE_NRI
  p_NRI           <- 2 * pnorm(-abs(z_NRI))
  
  # -- Formatting helpers -----------------------------------------------------
  fmt  <- function(x) formatC(x, format = "f", digits = digits)
  fmtp <- function(x) formatC(x, format = "g", digits = 3)
  pct  <- function(x) paste0(formatC(x * 100, format = "f", digits = 1), "%")
  
  ci_str <- function(est, lo, hi)
    paste0(fmt(est), " [", fmt(lo), ", ", fmt(hi), "]")
  
  # -- Results table ----------------------------------------------------------
  results <- data.frame(
    Component = c(
      "Cases: % reclassified up",
      "Cases: % reclassified down",
      "NRI (Cases)",
      "",
      "Controls: % reclassified up",
      "Controls: % reclassified down",
      "NRI (Controls)",
      "",
      "Overall NRI"
    ),
    N = c(
      n_cases, n_cases, n_cases,
      NA,
      n_controls, n_controls, n_controls,
      NA,
      n_cases + n_controls
    ),
    Estimate_95CI = c(
      pct(p_up_cases),
      pct(p_down_cases),
      ci_str(NRI_cases, CI_low_cases, CI_high_cases),
      "",
      pct(p_up_controls),
      pct(p_down_controls),
      ci_str(NRI_controls, CI_low_controls, CI_high_controls),
      "",
      ci_str(NRI, CI_low_NRI, CI_high_NRI)
    ),
    p_value = c(
      NA, NA, fmtp(p_cases),
      NA,
      NA, NA, fmtp(p_controls),
      NA,
      fmtp(p_NRI)
    ),
    stringsAsFactors = FALSE
  )
  
  if (show_output) {
    print(results, row.names = FALSE, na.print = "")
  }
  
  invisible(list(
    NRI          = NRI,
    CI_low_NRI   = CI_low_NRI,
    CI_high_NRI  = CI_high_NRI,
    p_NRI        = p_NRI,
    NRI_cases    = NRI_cases,
    CI_low_cases = CI_low_cases,
    CI_high_cases = CI_high_cases,
    p_cases      = p_cases,
    p_up_cases   = p_up_cases,
    p_down_cases = p_down_cases,
    NRI_controls    = NRI_controls,
    CI_low_controls = CI_low_controls,
    CI_high_controls = CI_high_controls,
    p_controls       = p_controls,
    p_up_controls    = p_up_controls,
    p_down_controls  = p_down_controls,
    n_cases      = n_cases,
    n_controls   = n_controls,
    table        = results
  ))
}

## Function to compute NRI for each GRS added to an Age-only reference model ---------------------
## 
## How it works:
##   For each combination of (outcome, population):
##     1. Fit an Age-only model (reference) ? generates pred column
##     2. For each GRS in grs_list:
##        - Fit Age + GRS model ? generates pred2 column
##        - Calculate NRI(Age+GRS vs Age only)
##        - Store in results dataframe
##
## Returns: dataframe with columns:
##   - Outcome, Population, Predictor
##   - N_Cases, N_Controls
##   - NRI, NRI_CI_Lower, NRI_CI_Upper, p_NRI
##   - NRI_Cases, NRI_Controls (component NRIs)
##   - Cases_Percent_Reclassified_Up, Cases_Percent_Reclassified_Down
##   - Controls_Percent_Reclassified_Up, Controls_Percent_Reclassified_Down

nri_table <- function(
    populations = NULL,
    grs_list = NULL,
    outcomes = c("PrCa", "PrCa_2yrs", "PrCa_5yrs", "PrCa_10yrs",
                 "PrCa_actionable", "PrCa_actionable_2yrs", "PrCa_actionable_5yrs", "PrCa_actionable_10yrs",
                 "PrCa_severe", "PrCa_severe_2yrs", "PrCa_severe_5yrs", "PrCa_severe_10yrs"),
    covariates_list = list("Age", "FH_PrCa_BrCa", "Age + FH_PrCa_BrCa"),
    subset_controls = FALSE,
    verbose = TRUE,
  version = "Asymptomatic Screening",
  adjust_for_PCs = FALSE
) {
  valid_versions <- c("Asymptomatic Screening", "Symptomatic Triage")
  if (!(version %in% valid_versions)) {
    stop("`version` must be one of: 'Asymptomatic Screening' or 'Symptomatic Triage'.")
  }
  
  # Default populations if not specified
  if (is.null(populations)) {
    if (version == "Asymptomatic Screening") {
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
    } else if (version == "Symptomatic Triage") {
      populations <- list(
        "All" = PCa_iv_covariates_GRS_predhorizon2,
        "White" = PCa_iv_covariates_GRS_predhorizon_WhiteOnly2,
        "Black" = PCa_iv_covariates_GRS_predhorizon_BlackOnly2,
        "Mixed" = PCa_iv_covariates_GRS_predhorizon_Mixed2,
        "Black+Mixed" = PCa_iv_covariates_GRS_predhorizon_BlackMixed2,
        "EUR" = PCa_iv_covariates_GRS_predhorizon_EUROnly2,
        "AFR" = PCa_iv_covariates_GRS_predhorizon_AFROnly2,
        "EAS" = PCa_iv_covariates_GRS_predhorizon_EASOnly2,
        "CSA" = PCa_iv_covariates_GRS_predhorizon_CSAOnly2,
        "MID" = PCa_iv_covariates_GRS_predhorizon_MIDOnly2,
        "AMR" = PCa_iv_covariates_GRS_predhorizon_AMROnly2
      )
    }
  }
  
  # Default GRS list if not specified
  if (is.null(grs_list)) {
    grs_list <- c(
      # Conti GRSs (Conti_script.R version)
      "ContimultiethnicGRS",
      "ContiEuropeanGRS",
      "ContiAfricanGRS",
      "ContiEast_AsianGRS",
      "ContiHispanicGRS",
      "ContiadjustedGRS",
      # Conti GRSs (Conti_GRS_267.R version)
      "ContimultiethnicGRS267",
      "ContiEuropeanGRS265",
      "ContiAfricanGRS246",
      "ContiEast_AsianGRS222",
      "ContiHispanicGRS253",
      "ContiORadjustedGRS",
      # Wang GRSs (Conti_script.R version)
      "WangmultiethnicGRS",
      "WangEuropeanGRS",
      "WangAfricanGRS",
      "WangEast_AsianGRS",
      "WangHispanicGRS",
      # Wang GRSs (Conti_GRS_267.R version)
      "WangmultiethnicGRS450",
      "WangEuropeanGRS445",
      "WangAfricanGRS444",
      "WangEast_AsianGRS379",
      "WangHispanicGRS446",
      # Schumacher and BARCODE1 GRSs
      "SchumacherGRS",
      "BARCODE1GRS",
      "SchumacherGRS145",
      "BARCODE1GRS129",
      # Seibert and Pagadala GRSs
      "SeibertGRS",
      "PagadalaGRS",
      "SeibertGRS52",
      "PagadalaGRS285",
      # Genomics PLC GRS
      "GenomicsPLC_PRS"
    )
  }
  
  grs_list <- setdiff(grs_list, "Age")
  pc_covariates <- paste0("PC", 1:10)

  parse_covariate_terms <- function(covariate_spec) {
    if (is.null(covariate_spec)) {
      return(character(0))
    }

    terms <- trimws(unlist(strsplit(as.character(covariate_spec), "\\+")))
    terms <- terms[nzchar(terms)]
    unique(terms)
  }

  resolve_covariates <- function(base_covariates, predictor = NULL) {
    covs <- parse_covariate_terms(base_covariates)

    if (isTRUE(adjust_for_PCs)) {
      covs <- c(covs, pc_covariates)
    }

    covs <- unique(covs)
    if (!is.null(predictor)) {
      covs <- covs[covs != predictor]
    }

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

  # Initialize results dataframe
  results <- data.frame(
    Outcome = character(),
    Population = character(),
    Predictor = character(),
    Covariates = character(),
    N_Cases = integer(),
    N_Controls = integer(),
    NRI = numeric(),
    NRI_CI_Lower = numeric(),
    NRI_CI_Upper = numeric(),
    p_NRI = numeric(),
    NRI_Cases = numeric(),
    NRI_Controls = numeric(),
    Cases_Percent_Reclassified_Up = numeric(),
    Cases_Percent_Reclassified_Down = numeric(),
    Controls_Percent_Reclassified_Up = numeric(),
    Controls_Percent_Reclassified_Down = numeric(),
    stringsAsFactors = FALSE
  )

  # Total combinations to process
  total_combos <- length(outcomes) * length(populations) * length(covariates_list) * length(grs_list)
  combo_count <- 0

  # Loop over all combinations
  for (current_outcome in outcomes) {
    for (pop_name in names(populations)) {
      pop_data <- populations[[pop_name]]

      for (cov in covariates_list) {
        reference_covariates <- resolve_covariates(cov)
        cov_label <- covariate_label(reference_covariates)

        if (verbose) {
          cat(sprintf("\n=== %s | %s | %s ===\n", current_outcome, pop_name, cov_label))
        }

        if (is.null(reference_covariates) || length(reference_covariates) == 0) {
          if (verbose) cat("Skipping (reference covariate set is empty)\n")
          combo_count <- combo_count + length(grs_list)
          next
        }

        missing_reference_covariates <- setdiff(reference_covariates, colnames(pop_data))
        if (length(missing_reference_covariates) > 0) {
          if (verbose) {
            cat(sprintf("Skipping (missing covariates in %s: %s)\n",
                        pop_name, paste(missing_reference_covariates, collapse = ", ")))
          }
          combo_count <- combo_count + length(grs_list)
          next
        }

        reference_predictor <- reference_covariates[1]
        reference_other_covariates <- if (length(reference_covariates) > 1) reference_covariates[-1] else NULL

        tryCatch({
          reference_model <- run_logreg(
            data = pop_data,
            outcome = current_outcome,
            predictor = reference_predictor,
            covariates = reference_other_covariates,
            subset_controls = subset_controls,
            plot_roc = FALSE,
            show_output = FALSE
          )

          reference_data <- reference_model$data

          for (grs_pred in grs_list) {
            combo_count <- combo_count + 1

            if (!grs_pred %in% colnames(pop_data)) {
              if (verbose) cat(sprintf("[%d/%d] Skipping %s (not in %s)\n",
                                       combo_count, total_combos, grs_pred, pop_name))
              next
            }

            if (verbose) {
              cat(sprintf("[%d/%d] %s + %s | %s\n",
                          combo_count, total_combos, current_outcome, grs_pred, cov_label))
            }

            grs_covariates <- resolve_covariates(cov, predictor = grs_pred)
            missing_grs_covariates <- setdiff(grs_covariates, colnames(reference_data))
            if (length(missing_grs_covariates) > 0) {
              if (verbose) {
                cat(sprintf("  Skipping (missing covariates in %s: %s)\n",
                            pop_name, paste(missing_grs_covariates, collapse = ", ")))
              }
              next
            }

            tryCatch({
              grs_model <- run_logreg(
                data = reference_data,
                outcome = current_outcome,
                predictor = grs_pred,
                covariates = grs_covariates,
                subset_controls = subset_controls,
                plot_roc = FALSE,
                show_output = FALSE
              )

              nri_data <- grs_model$data
              names(nri_data)[names(nri_data) == "pred"] <- "pred2"
              nri_data$pred <- reference_data$pred

              nri_result <- nri(
                data = nri_data,
                outcome = current_outcome,
                pred_new = "pred2",
                pred_old = "pred",
                digits = 3,
                show_output = FALSE
              )

              results <- rbind(results, data.frame(
                Outcome = current_outcome,
                Population = pop_name,
                Predictor = grs_pred,
                Covariates = cov_label,
                N_Cases = nri_result$n_cases,
                N_Controls = nri_result$n_controls,
                NRI = nri_result$NRI,
                NRI_CI_Lower = nri_result$CI_low_NRI,
                NRI_CI_Upper = nri_result$CI_high_NRI,
                p_NRI = nri_result$p_NRI,
                NRI_Cases = nri_result$NRI_cases,
                NRI_Controls = nri_result$NRI_controls,
                Cases_Percent_Reclassified_Up = nri_result$p_up_cases * 100,
                Cases_Percent_Reclassified_Down = nri_result$p_down_cases * 100,
                Controls_Percent_Reclassified_Up = nri_result$p_up_controls * 100,
                Controls_Percent_Reclassified_Down = nri_result$p_down_controls * 100,
                stringsAsFactors = FALSE
              ))
            }, error = function(e) {
              if (verbose) cat(sprintf("  ERROR: %s\n", e$message))
            })
          }
        }, error = function(e) {
          if (verbose) cat(sprintf("  ERROR fitting reference model: %s\n", e$message))
        })
      }
    }
  }

  return(results)
}

## Function to compute logistic regression for all outcome x GRS x population x covariate combinations ---
## 
## Version 1 - for Logistic_Regressions_testing_Conti_GRS.R
##
## How it works:
##   For each combination of (outcome, population, GRS, covariates):
##     - Fit a logistic regression model
##     - Compute ROC AUC with 95% CI
##     - Store outcome, population, GRS, covariates, sample sizes, AUC, and 95% CI in results dataframe
##
## Optionally includes Age-only baseline models for comparison.
##
## Returns: dataframe with columns:
##   - Outcome, Population, Predictor, Covariates
##   - N_Cases, N_Controls
##   - Prevalence, PR_AUC
##   - ROC_AUC, ROC_AUC_CI_Lower, ROC_AUC_CI_Upper
##   - OR_per_1SD, OR_per_1SD_CI_Lower, OR_per_1SD_CI_Upper, OR_per_1SD_p

logreg_table <- function(
    populations = NULL,
    grs_list = NULL,
    outcomes = c("PrCa", "PrCa_2yrs", "PrCa_5yrs", "PrCa_10yrs",
                 "PrCa_actionable", "PrCa_actionable_2yrs", "PrCa_actionable_5yrs", "PrCa_actionable_10yrs",
                 "PrCa_severe", "PrCa_severe_2yrs", "PrCa_severe_5yrs", "PrCa_severe_10yrs"),
    covariates_list = list(NULL, "Age", "FH_PrCa_BrCa", "Age + FH_PrCa_BrCa"),  # uses list() so NULL is preserved as a distinct option
    include_age_only = TRUE,
    subset_controls = FALSE,
    plot_roc = FALSE,
    verbose = TRUE,
    version = "Asymptomatic Screening",
    adjust_for_PCs = FALSE
) {
  valid_versions <- c("Asymptomatic Screening", "Symptomatic Triage")
  if (!(version %in% valid_versions)) {
    stop("`version` must be one of: 'Asymptomatic Screening' or 'Symptomatic Triage'.")
  }
  
  # Default populations if not specified
  
  if (is.null(populations)) {
    if (version == "Asymptomatic Screening") {
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
    } else if (version == "Symptomatic Triage") {
      populations <- list(
        "All" = PCa_iv_covariates_GRS_predhorizon2,
        "White" = PCa_iv_covariates_GRS_predhorizon_WhiteOnly2,
        "Black" = PCa_iv_covariates_GRS_predhorizon_BlackOnly2,
        "Mixed" = PCa_iv_covariates_GRS_predhorizon_Mixed2,
        "Black+Mixed" = PCa_iv_covariates_GRS_predhorizon_BlackMixed2,
        "EUR" = PCa_iv_covariates_GRS_predhorizon_EUROnly2,
        "AFR" = PCa_iv_covariates_GRS_predhorizon_AFROnly2,
        "EAS" = PCa_iv_covariates_GRS_predhorizon_EASOnly2,
        "CSA" = PCa_iv_covariates_GRS_predhorizon_CSAOnly2,
        "MID" = PCa_iv_covariates_GRS_predhorizon_MIDOnly2,
        "AMR" = PCa_iv_covariates_GRS_predhorizon_AMROnly2
      )
    }
  }
  
  # Default GRS list if not specified
  if (is.null(grs_list)) {
    grs_list <- c(
      # Conti GRSs (Conti_script.R version)
      "ContimultiethnicGRS",
      "ContiEuropeanGRS",
      "ContiAfricanGRS",
      "ContiEast_AsianGRS",
      "ContiHispanicGRS",
      "ContiadjustedGRS",
      # Conti GRSs (Conti_GRS_267.R version)
      "ContimultiethnicGRS267",
      "ContiEuropeanGRS265",
      "ContiAfricanGRS246",
      "ContiEast_AsianGRS222",
      "ContiHispanicGRS253",
      "ContiORadjustedGRS",
      # Wang GRSs (Conti_script.R version)
      "WangmultiethnicGRS",
      "WangEuropeanGRS",
      "WangAfricanGRS",
      "WangEast_AsianGRS",
      "WangHispanicGRS",
      # Wang GRSs (Conti_GRS_267.R version)
      "WangmultiethnicGRS450",
      "WangEuropeanGRS445",
      "WangAfricanGRS444",
      "WangEast_AsianGRS379",
      "WangHispanicGRS446",
      # Schumacher and BARCODE1 GRSs
      "SchumacherGRS",
      "BARCODE1GRS",
      "SchumacherGRS145",
      "BARCODE1GRS129",
      # Seibert and Pagadala GRSs
      "SeibertGRS",
      "PagadalaGRS",
      "SeibertGRS52",
      "PagadalaGRS285",
      # Genomics PLC GRS
      "GenomicsPLC_PRS"
    )
  }
  
  grs_list <- setdiff(grs_list, "Age")
  pc_covariates <- paste0("PC", 1:10)

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
  
  # Initialize results dataframe
  results <- data.frame(
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
    OR_per_1SD = numeric(),
    OR_per_1SD_CI_Lower = numeric(),
    OR_per_1SD_CI_Upper = numeric(),
    OR_per_1SD_p = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Total combinations to process
  total_combos <- length(outcomes) * length(populations) * length(grs_list) * length(covariates_list)
  if (include_age_only) {
    total_combos <- total_combos + length(outcomes) * length(populations)
  }
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
          effective_covariates <- resolve_covariates(cov, predictor = grs_pred)
          cov_label <- covariate_label(effective_covariates)
          
          # Progress indicator
          if (verbose) {
            cat(sprintf("[%d/%d] %s | %s | %s | %s\n", 
                        combo_count, total_combos, current_outcome, pop_name, grs_pred, cov_label))
          }

          missing_covariates <- setdiff(effective_covariates, colnames(pop_data))
          if (length(missing_covariates) > 0) {
            if (verbose) {
              cat(sprintf("  Skipping (missing covariates in %s: %s)\n",
                          pop_name, paste(missing_covariates, collapse = ", ")))
            }
            next
          }
          
          # Run logistic regression
          tryCatch({
            model <- run_logreg(
              data = pop_data,
              outcome = current_outcome,
              predictor = grs_pred,
              covariates = effective_covariates,
              subset_controls = subset_controls,
              plot_roc = plot_roc,
              plot_pr = TRUE,
              show_output = FALSE
            )
            
            # Extract ROC AUC and CI
            roc_obj <- model$roc
            auc_val <- as.numeric(roc_obj$auc)
            auc_ci <- as.numeric(roc_obj$ci)  # Returns [lower, AUC, upper]
            auc_lower <- auc_ci[1]
            auc_upper <- auc_ci[3]

            # Extract PR AUC and prevalence from rows with non-missing outcome/predictions
            pr_auc <- if (!is.null(model$pr) && !is.null(model$pr$auc.integral)) {
              as.numeric(model$pr$auc.integral)
            } else {
              NA_real_
            }
            valid_rows <- !is.na(model$data[[current_outcome]]) & !is.na(model$data$pred)
            prevalence <- mean(model$data[[current_outcome]][valid_rows] == 1)
            
            # Count cases and controls in the dataset used
            n_cases <- sum(model$data[[current_outcome]] == 1, na.rm = TRUE)
            n_controls <- sum(model$data[[current_outcome]] == 0, na.rm = TRUE)

            # OR per 1 SD increase in predictor (predictors are z-transformed upstream)
            coef_tbl <- summary(model$model)$coefficients
            if (grs_pred %in% rownames(coef_tbl)) {
              beta <- coef_tbl[grs_pred, "Estimate"]
              se <- coef_tbl[grs_pred, "Std. Error"]
              p_val <- coef_tbl[grs_pred, "Pr(>|z|)"]
              or_val <- exp(beta)
              or_ci_lower <- exp(beta - 1.96 * se)
              or_ci_upper <- exp(beta + 1.96 * se)
            } else {
              or_val <- NA_real_
              or_ci_lower <- NA_real_
              or_ci_upper <- NA_real_
              p_val <- NA_real_
            }
            
            # Append to results
            results <- rbind(results, data.frame(
              Outcome = current_outcome,
              Population = pop_name,
              Predictor = grs_pred,
              Covariates = cov_label,
              N_Cases = n_cases,
              N_Controls = n_controls,
              Prevalence = prevalence,
              PR_AUC = pr_auc,
              ROC_AUC = auc_val,
              ROC_AUC_CI_Lower = auc_lower,
              ROC_AUC_CI_Upper = auc_upper,
              OR_per_1SD = or_val,
              OR_per_1SD_CI_Lower = or_ci_lower,
              OR_per_1SD_CI_Upper = or_ci_upper,
              OR_per_1SD_p = p_val,
              stringsAsFactors = FALSE
            ))
            
          }, error = function(e) {
            if (verbose) cat(sprintf("  ERROR: %s\n", e$message))
          })
        }
      }
      
      if (include_age_only) {
        combo_count <- combo_count + 1
        age_covariates <- resolve_covariates(NULL, predictor = "Age")
        age_cov_label <- covariate_label(age_covariates)
        
        if (verbose) {
          cat(sprintf("[%d/%d] %s | %s | %s | %s\n",
                      combo_count, total_combos, current_outcome, pop_name, "Age", age_cov_label))
        }
        
        if (!"Age" %in% colnames(pop_data)) {
          if (verbose) cat(sprintf("Skipping Age (not in %s)\n", pop_name))
        } else {
          missing_age_covariates <- setdiff(age_covariates, colnames(pop_data))
          if (length(missing_age_covariates) > 0) {
            if (verbose) {
              cat(sprintf("  Skipping Age model (missing covariates in %s: %s)\n",
                          pop_name, paste(missing_age_covariates, collapse = ", ")))
            }
            next
          }

          tryCatch({
            model <- run_logreg(
              data = pop_data,
              outcome = current_outcome,
              predictor = "Age",
              covariates = age_covariates,
              subset_controls = subset_controls,
              plot_roc = plot_roc,
              plot_pr = TRUE,
              show_output = FALSE
            )
            
            roc_obj <- model$roc
            auc_val <- as.numeric(roc_obj$auc)
            auc_ci <- as.numeric(roc_obj$ci)
            auc_lower <- auc_ci[1]
            auc_upper <- auc_ci[3]

            pr_auc <- if (!is.null(model$pr) && !is.null(model$pr$auc.integral)) {
              as.numeric(model$pr$auc.integral)
            } else {
              NA_real_
            }
            valid_rows <- !is.na(model$data[[current_outcome]]) & !is.na(model$data$pred)
            prevalence <- mean(model$data[[current_outcome]][valid_rows] == 1)
            
            n_cases <- sum(model$data[[current_outcome]] == 1, na.rm = TRUE)
            n_controls <- sum(model$data[[current_outcome]] == 0, na.rm = TRUE)
            
            results <- rbind(results, data.frame(
              Outcome = current_outcome,
              Population = pop_name,
              Predictor = "Age",
              Covariates = age_cov_label,
              N_Cases = n_cases,
              N_Controls = n_controls,
              Prevalence = prevalence,
              PR_AUC = pr_auc,
              ROC_AUC = auc_val,
              ROC_AUC_CI_Lower = auc_lower,
              ROC_AUC_CI_Upper = auc_upper,
              OR_per_1SD = NA_real_,
              OR_per_1SD_CI_Lower = NA_real_,
              OR_per_1SD_CI_Upper = NA_real_,
              OR_per_1SD_p = NA_real_,
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

## Function to quantify feature importance in integrated GRS+Age models ---------
##
## How it works:
##   For each (outcome, population, GRS):
##     1. Build a common analysis dataset with complete rows for outcome + GRS + Age
##     2. Fit full model: outcome ~ GRS + Age
##     3. Fit reduced models: outcome ~ Age and outcome ~ GRS
##     4. Quantify feature importance with:
##        - Likelihood-ratio test (LRT) drop-in-fit statistics
##        - Delta ROC AUC (full minus reduced)
##     5. Report adjusted ORs for GRS and Age from the full model

feature_importance_table <- function(
    populations = NULL,
    grs_list = NULL,
    outcomes = c("PrCa", "PrCa_2yrs", "PrCa_5yrs", "PrCa_10yrs",
                 "PrCa_actionable", "PrCa_actionable_2yrs", "PrCa_actionable_5yrs", "PrCa_actionable_10yrs",
                 "PrCa_severe", "PrCa_severe_2yrs", "PrCa_severe_5yrs", "PrCa_severe_10yrs"),
    covariates_list = list("Age", "FH_PrCa_BrCa", "Age + FH_PrCa_BrCa"),
    subset_controls = FALSE,
    verbose = TRUE,
  version = "Asymptomatic Screening",
  adjust_for_PCs = FALSE
) {
  valid_versions <- c("Asymptomatic Screening", "Symptomatic Triage")
  if (!(version %in% valid_versions)) {
    stop("`version` must be one of: 'Asymptomatic Screening' or 'Symptomatic Triage'.")
  }

  # Default populations if not specified
  if (is.null(populations)) {
    if (version == "Asymptomatic Screening") {
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
    } else if (version == "Symptomatic Triage") {
      populations <- list(
        "All" = PCa_iv_covariates_GRS_predhorizon2,
        "White" = PCa_iv_covariates_GRS_predhorizon_WhiteOnly2,
        "Black" = PCa_iv_covariates_GRS_predhorizon_BlackOnly2,
        "Mixed" = PCa_iv_covariates_GRS_predhorizon_Mixed2,
        "Black+Mixed" = PCa_iv_covariates_GRS_predhorizon_BlackMixed2,
        "EUR" = PCa_iv_covariates_GRS_predhorizon_EUROnly2,
        "AFR" = PCa_iv_covariates_GRS_predhorizon_AFROnly2,
        "EAS" = PCa_iv_covariates_GRS_predhorizon_EASOnly2,
        "CSA" = PCa_iv_covariates_GRS_predhorizon_CSAOnly2,
        "MID" = PCa_iv_covariates_GRS_predhorizon_MIDOnly2,
        "AMR" = PCa_iv_covariates_GRS_predhorizon_AMROnly2
      )
    }
  }

  # Default GRS list if not specified
  if (is.null(grs_list)) {
    grs_list <- c(
      "ContimultiethnicGRS", "ContiEuropeanGRS", "ContiAfricanGRS", "ContiEast_AsianGRS", "ContiHispanicGRS",
      "ContiadjustedGRS", "ContimultiethnicGRS267", "ContiEuropeanGRS265", "ContiAfricanGRS246",
      "ContiEast_AsianGRS222", "ContiHispanicGRS253", "ContiORadjustedGRS", "WangmultiethnicGRS",
      "WangEuropeanGRS", "WangAfricanGRS", "WangEast_AsianGRS", "WangHispanicGRS", "WangmultiethnicGRS450",
      "WangEuropeanGRS445", "WangAfricanGRS444", "WangEast_AsianGRS379", "WangHispanicGRS446",
      "SchumacherGRS", "BARCODE1GRS", "SchumacherGRS145", "BARCODE1GRS129", "SeibertGRS", "PagadalaGRS",
      "SeibertGRS52", "PagadalaGRS285", "GenomicsPLC_PRS"
    )
  }

  grs_list <- setdiff(grs_list, "Age")
  pc_covariates <- paste0("PC", 1:10)

  parse_covariate_terms <- function(covariate_spec) {
    if (is.null(covariate_spec)) {
      return(character(0))
    }

    terms <- trimws(unlist(strsplit(as.character(covariate_spec), "\\+")))
    terms <- terms[nzchar(terms)]
    unique(terms)
  }

  resolve_covariates <- function(base_covariates, predictor = NULL) {
    covs <- parse_covariate_terms(base_covariates)

    if (isTRUE(adjust_for_PCs)) {
      covs <- c(covs, pc_covariates)
    }

    covs <- unique(covs)
    if (!is.null(predictor)) {
      covs <- covs[covs != predictor]
    }

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

  results <- data.frame(
    Outcome = character(),
    Population = character(),
    Predictor = character(),
    Covariates = character(),
    N = integer(),
    N_Cases = integer(),
    N_Controls = integer(),
    Full_ROC_AUC = numeric(),
    Age_only_ROC_AUC = numeric(),
    GRS_only_ROC_AUC = numeric(),
    Delta_AUC_drop_GRS = numeric(),
    Delta_AUC_drop_Age = numeric(),
    LRT_drop_GRS_ChiSq = numeric(),
    LRT_drop_GRS_p = numeric(),
    LRT_drop_Age_ChiSq = numeric(),
    LRT_drop_Age_p = numeric(),
    GRS_OR_adj = numeric(),
    GRS_OR_adj_CI_Lower = numeric(),
    GRS_OR_adj_CI_Upper = numeric(),
    GRS_OR_adj_p = numeric(),
    Age_OR_adj = numeric(),
    Age_OR_adj_CI_Lower = numeric(),
    Age_OR_adj_CI_Upper = numeric(),
    Age_OR_adj_p = numeric(),
    stringsAsFactors = FALSE
  )

  total_combos <- length(outcomes) * length(populations) * length(covariates_list) * length(grs_list)
  combo_count <- 0

  for (current_outcome in outcomes) {
    for (pop_name in names(populations)) {
      pop_data <- populations[[pop_name]]

      for (cov in covariates_list) {
        base_covariates <- resolve_covariates(cov)
        cov_label <- covariate_label(base_covariates)

        if (is.null(base_covariates) || length(base_covariates) == 0) {
          if (verbose) cat(sprintf("Skipping empty covariate set for %s\n", pop_name))
          combo_count <- combo_count + length(grs_list)
          next
        }

        for (grs_pred in grs_list) {
          combo_count <- combo_count + 1

          if (verbose) {
            cat(sprintf("[%d/%d] %s | %s | %s | %s\n",
                        combo_count, total_combos, current_outcome, pop_name, grs_pred, cov_label))
          }

          if (!grs_pred %in% colnames(pop_data)) {
            if (verbose) cat(sprintf("  Skipping %s (not in %s)\n", grs_pred, pop_name))
            next
          }

          full_covariates <- resolve_covariates(cov, predictor = grs_pred)
          cov_only_terms <- base_covariates
          grs_only_covariates <- resolve_covariates(NULL, predictor = grs_pred)

          needed <- unique(c(current_outcome, grs_pred, full_covariates, cov_only_terms, grs_only_covariates))
          missing_needed <- setdiff(needed, colnames(pop_data))
          if (length(missing_needed) > 0) {
            if (verbose) {
              cat(sprintf("  Skipping (missing columns in %s: %s)\n",
                          pop_name, paste(missing_needed, collapse = ", ")))
            }
            next
          }

          analysis_data <- pop_data[stats::complete.cases(pop_data[, needed]), , drop = FALSE]

          if (nrow(analysis_data) == 0) {
            if (verbose) cat("  Skipping (no complete rows for this model)\n")
            next
          }

          if (isTRUE(subset_controls)) {
            analysis_data <- subset_controls_to_case_count(analysis_data, case_col = current_outcome)
          }

          if (length(unique(analysis_data[[current_outcome]])) < 2) {
            if (verbose) cat("  Skipping (outcome has <2 classes after filtering)\n")
            next
          }

          cov_only_predictor <- cov_only_terms[1]
          cov_only_covariates <- if (length(cov_only_terms) > 1) cov_only_terms[-1] else NULL

          tryCatch({
            full_model <- run_logreg(
              data = analysis_data,
              outcome = current_outcome,
              predictor = grs_pred,
              covariates = full_covariates,
              subset_controls = FALSE,
              plot_roc = FALSE,
              plot_pr = FALSE,
              show_output = FALSE
            )

            cov_only_model <- run_logreg(
              data = analysis_data,
              outcome = current_outcome,
              predictor = cov_only_predictor,
              covariates = cov_only_covariates,
              subset_controls = FALSE,
              plot_roc = FALSE,
              plot_pr = FALSE,
              show_output = FALSE
            )

            grs_only_model <- run_logreg(
              data = analysis_data,
              outcome = current_outcome,
              predictor = grs_pred,
              covariates = grs_only_covariates,
              subset_controls = FALSE,
              plot_roc = FALSE,
              plot_pr = FALSE,
              show_output = FALSE
            )

            full_auc <- as.numeric(full_model$roc$auc)
            cov_only_auc <- as.numeric(cov_only_model$roc$auc)
            grs_auc <- as.numeric(grs_only_model$roc$auc)

            lrt_drop_grs <- anova(cov_only_model$model, full_model$model, test = "LRT")
            lrt_drop_age <- anova(grs_only_model$model, full_model$model, test = "LRT")

            lrt_drop_grs_chisq <- as.numeric(lrt_drop_grs$Deviance[2])
            lrt_drop_grs_p <- as.numeric(lrt_drop_grs$`Pr(>Chi)`[2])
            lrt_drop_age_chisq <- as.numeric(lrt_drop_age$Deviance[2])
            lrt_drop_age_p <- as.numeric(lrt_drop_age$`Pr(>Chi)`[2])

            coef_tbl <- summary(full_model$model)$coefficients

            if (grs_pred %in% rownames(coef_tbl)) {
              grs_beta <- coef_tbl[grs_pred, "Estimate"]
              grs_se <- coef_tbl[grs_pred, "Std. Error"]
              grs_or <- exp(grs_beta)
              grs_or_ci_lower <- exp(grs_beta - 1.96 * grs_se)
              grs_or_ci_upper <- exp(grs_beta + 1.96 * grs_se)
              grs_p <- coef_tbl[grs_pred, "Pr(>|z|)"]
            } else {
              grs_or <- NA_real_
              grs_or_ci_lower <- NA_real_
              grs_or_ci_upper <- NA_real_
              grs_p <- NA_real_
            }

            if ("Age" %in% rownames(coef_tbl)) {
              age_beta <- coef_tbl["Age", "Estimate"]
              age_se <- coef_tbl["Age", "Std. Error"]
              age_or <- exp(age_beta)
              age_or_ci_lower <- exp(age_beta - 1.96 * age_se)
              age_or_ci_upper <- exp(age_beta + 1.96 * age_se)
              age_p <- coef_tbl["Age", "Pr(>|z|)"]
            } else {
              age_or <- NA_real_
              age_or_ci_lower <- NA_real_
              age_or_ci_upper <- NA_real_
              age_p <- NA_real_
            }

            n_cases <- sum(analysis_data[[current_outcome]] == 1, na.rm = TRUE)
            n_controls <- sum(analysis_data[[current_outcome]] == 0, na.rm = TRUE)

            results <- rbind(results, data.frame(
              Outcome = current_outcome,
              Population = pop_name,
              Predictor = grs_pred,
              Covariates = cov_label,
              N = nrow(analysis_data),
              N_Cases = n_cases,
              N_Controls = n_controls,
              Full_ROC_AUC = full_auc,
              Age_only_ROC_AUC = cov_only_auc,
              GRS_only_ROC_AUC = grs_auc,
              Delta_AUC_drop_GRS = full_auc - cov_only_auc,
              Delta_AUC_drop_Age = full_auc - grs_auc,
              LRT_drop_GRS_ChiSq = lrt_drop_grs_chisq,
              LRT_drop_GRS_p = lrt_drop_grs_p,
              LRT_drop_Age_ChiSq = lrt_drop_age_chisq,
              LRT_drop_Age_p = lrt_drop_age_p,
              GRS_OR_adj = grs_or,
              GRS_OR_adj_CI_Lower = grs_or_ci_lower,
              GRS_OR_adj_CI_Upper = grs_or_ci_upper,
              GRS_OR_adj_p = grs_p,
              Age_OR_adj = age_or,
              Age_OR_adj_CI_Lower = age_or_ci_lower,
              Age_OR_adj_CI_Upper = age_or_ci_upper,
              Age_OR_adj_p = age_p,
              stringsAsFactors = FALSE
            ))
          }, error = function(e) {
            if (verbose) cat(sprintf("  ERROR: %s\n", e$message))
          })
        }
      }
    }
  }

  results
}

## Function to compute Random Forest feature importance for integrated GRS+Age models ---
##
## How it works:
##   For each (outcome, population, GRS):
##     1. Build a complete-case dataset for outcome + GRS + Age
##     2. Fit a Random Forest: outcome ~ GRS + Age (classification)
##     3. Extract permutation importance (Mean Decrease Accuracy) for GRS and Age
##        - This measures how much OOB accuracy drops when each feature's values
##          are randomly shuffled, i.e. how much the model relies on it
##     4. Compute OOB ROC AUC from OOB predicted probabilities
##
## Returns: dataframe with columns:
##   - Outcome, Population, Predictor, N, N_Cases, N_Controls
##   - OOB_ROC_AUC
##   - GRS_MDA, Age_MDA        (Mean Decrease Accuracy - permutation importance)
##   - GRS_MDG, Age_MDG        (Mean Decrease Gini - impurity-based, less preferred)
##   - GRS_Pct_Importance      (GRS_MDA as % of total MDA; relative contribution)
##   - Age_Pct_Importance

rf_feature_importance_table <- function(
    populations = NULL,
    grs_list = NULL,
    outcomes = c("PrCa", "PrCa_2yrs", "PrCa_5yrs", "PrCa_10yrs",
                 "PrCa_actionable", "PrCa_actionable_2yrs", "PrCa_actionable_5yrs", "PrCa_actionable_10yrs",
                 "PrCa_severe", "PrCa_severe_2yrs", "PrCa_severe_5yrs", "PrCa_severe_10yrs"),
    covariates_list = list("Age", "FH_PrCa_BrCa", "Age + FH_PrCa_BrCa"),
    ntree = 500,
    subset_controls = FALSE,
    verbose = TRUE,
  version = "Asymptomatic Screening",
  adjust_for_PCs = FALSE
) {
  valid_versions <- c("Asymptomatic Screening", "Symptomatic Triage")
  if (!(version %in% valid_versions)) {
    stop("`version` must be one of: 'Asymptomatic Screening' or 'Symptomatic Triage'.")
  }

  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required. Install it with install.packages('randomForest').")
  }

  # Default populations if not specified
  if (is.null(populations)) {
    if (version == "Asymptomatic Screening") {
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
    } else if (version == "Symptomatic Triage") {
      populations <- list(
        "All" = PCa_iv_covariates_GRS_predhorizon2,
        "White" = PCa_iv_covariates_GRS_predhorizon_WhiteOnly2,
        "Black" = PCa_iv_covariates_GRS_predhorizon_BlackOnly2,
        "Mixed" = PCa_iv_covariates_GRS_predhorizon_Mixed2,
        "Black+Mixed" = PCa_iv_covariates_GRS_predhorizon_BlackMixed2,
        "EUR" = PCa_iv_covariates_GRS_predhorizon_EUROnly2,
        "AFR" = PCa_iv_covariates_GRS_predhorizon_AFROnly2,
        "EAS" = PCa_iv_covariates_GRS_predhorizon_EASOnly2,
        "CSA" = PCa_iv_covariates_GRS_predhorizon_CSAOnly2,
        "MID" = PCa_iv_covariates_GRS_predhorizon_MIDOnly2,
        "AMR" = PCa_iv_covariates_GRS_predhorizon_AMROnly2
      )
    }
  }

  # Default GRS list if not specified
  if (is.null(grs_list)) {
    grs_list <- c(
      "ContimultiethnicGRS", "ContiEuropeanGRS", "ContiAfricanGRS", "ContiEast_AsianGRS", "ContiHispanicGRS",
      "ContiadjustedGRS", "ContimultiethnicGRS267", "ContiEuropeanGRS265", "ContiAfricanGRS246",
      "ContiEast_AsianGRS222", "ContiHispanicGRS253", "ContiORadjustedGRS", "WangmultiethnicGRS",
      "WangEuropeanGRS", "WangAfricanGRS", "WangEast_AsianGRS", "WangHispanicGRS", "WangmultiethnicGRS450",
      "WangEuropeanGRS445", "WangAfricanGRS444", "WangEast_AsianGRS379", "WangHispanicGRS446",
      "SchumacherGRS", "BARCODE1GRS", "SchumacherGRS145", "BARCODE1GRS129", "SeibertGRS", "PagadalaGRS",
      "SeibertGRS52", "PagadalaGRS285", "GenomicsPLC_PRS"
    )
  }

  grs_list <- setdiff(grs_list, "Age")
  pc_covariates <- paste0("PC", 1:10)

  parse_covariate_terms <- function(covariate_spec) {
    if (is.null(covariate_spec)) {
      return(character(0))
    }

    terms <- trimws(unlist(strsplit(as.character(covariate_spec), "\\+")))
    terms <- terms[nzchar(terms)]
    unique(terms)
  }

  resolve_covariates <- function(base_covariates, predictor = NULL) {
    covs <- parse_covariate_terms(base_covariates)

    if (isTRUE(adjust_for_PCs)) {
      covs <- c(covs, pc_covariates)
    }

    covs <- unique(covs)
    if (!is.null(predictor)) {
      covs <- covs[covs != predictor]
    }

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

  results <- data.frame(
    Outcome = character(),
    Population = character(),
    Predictor = character(),
    Covariates = character(),
    N = integer(),
    N_Cases = integer(),
    N_Controls = integer(),
    OOB_ROC_AUC = numeric(),
    GRS_MDA = numeric(),
    Age_MDA = numeric(),
    GRS_MDG = numeric(),
    Age_MDG = numeric(),
    GRS_Pct_Importance = numeric(),
    Age_Pct_Importance = numeric(),
    stringsAsFactors = FALSE
  )

  total_combos <- length(outcomes) * length(populations) * length(covariates_list) * length(grs_list)
  combo_count <- 0

  for (current_outcome in outcomes) {
    for (pop_name in names(populations)) {
      pop_data <- populations[[pop_name]]

      for (cov in covariates_list) {
        base_covariates <- resolve_covariates(cov)
        cov_label <- covariate_label(base_covariates)

        if (is.null(base_covariates) || length(base_covariates) == 0) {
          if (verbose) cat(sprintf("Skipping empty covariate set for %s\n", pop_name))
          combo_count <- combo_count + length(grs_list)
          next
        }

        for (grs_pred in grs_list) {
          combo_count <- combo_count + 1

          if (verbose) {
            cat(sprintf("[%d/%d] %s | %s | %s | %s\n",
                        combo_count, total_combos, current_outcome, pop_name, grs_pred, cov_label))
          }

          if (!grs_pred %in% colnames(pop_data)) {
            if (verbose) cat(sprintf("  Skipping %s (not in %s)\n", grs_pred, pop_name))
            next
          }

          rf_terms <- unique(c(grs_pred, base_covariates))
          needed <- unique(c(current_outcome, rf_terms))
          missing_needed <- setdiff(needed, colnames(pop_data))
          if (length(missing_needed) > 0) {
            if (verbose) {
              cat(sprintf("  Skipping (missing columns in %s: %s)\n",
                          pop_name, paste(missing_needed, collapse = ", ")))
            }
            next
          }
          analysis_data <- pop_data[stats::complete.cases(pop_data[, needed]), needed, drop = FALSE]

          if (nrow(analysis_data) == 0) {
            if (verbose) cat("  Skipping (no complete rows)\n")
            next
          }

          if (isTRUE(subset_controls)) {
            analysis_data <- subset_controls_to_case_count(analysis_data, case_col = current_outcome)
          }

          # outcome must be a factor with exactly 2 levels for RF classification
          analysis_data[[current_outcome]] <- factor(analysis_data[[current_outcome]])

          if (length(levels(analysis_data[[current_outcome]])) < 2) {
            if (verbose) cat("  Skipping (outcome has <2 classes after filtering)\n")
            next
          }

          tryCatch({
            rf_formula <- stats::as.formula(paste(current_outcome, "~", paste(rf_terms, collapse = " + ")))

            rf_model <- randomForest::randomForest(
              formula = rf_formula,
              data = analysis_data,
              ntree = ntree,
              importance = TRUE,  # enables permutation importance
              keep.inbag = FALSE
            )

            # Permutation importance (Mean Decrease Accuracy) - type = 1
            imp <- randomForest::importance(rf_model, type = 1, scale = TRUE)
            grs_mda <- imp[grs_pred, "MeanDecreaseAccuracy"]
            age_mda <- if ("Age" %in% rownames(imp)) imp["Age", "MeanDecreaseAccuracy"] else NA_real_

            # Gini importance (Mean Decrease Gini) - type = 2
            imp_gini <- randomForest::importance(rf_model, type = 2)
            grs_mdg <- imp_gini[grs_pred, "MeanDecreaseGini"]
            age_mdg <- if ("Age" %in% rownames(imp_gini)) imp_gini["Age", "MeanDecreaseGini"] else NA_real_

            # Relative importance (% of total MDA, treating negative MDA as 0 for the ratio)
            total_mda <- max(grs_mda, 0) + ifelse(is.na(age_mda), 0, max(age_mda, 0))
            grs_pct <- if (total_mda > 0) max(grs_mda, 0) / total_mda * 100 else NA_real_
            age_pct <- if (total_mda > 0 && !is.na(age_mda)) max(age_mda, 0) / total_mda * 100 else NA_real_

            # OOB ROC AUC using OOB vote probabilities for the positive class ("1")
            oob_votes <- rf_model$votes
            pos_level <- "1"
            if (!pos_level %in% colnames(oob_votes)) {
              pos_level <- levels(analysis_data[[current_outcome]])[2]
            }
            oob_probs <- oob_votes[, pos_level]
            true_labels <- as.integer(as.character(analysis_data[[current_outcome]]))
            oob_roc <- pROC::roc(true_labels ~ oob_probs, quiet = TRUE)
            oob_auc <- as.numeric(oob_roc$auc)

            n_cases <- sum(as.character(analysis_data[[current_outcome]]) == "1", na.rm = TRUE)
            n_controls <- sum(as.character(analysis_data[[current_outcome]]) == "0", na.rm = TRUE)

            results <- rbind(results, data.frame(
              Outcome = current_outcome,
              Population = pop_name,
              Predictor = grs_pred,
              Covariates = cov_label,
              N = nrow(analysis_data),
              N_Cases = n_cases,
              N_Controls = n_controls,
              OOB_ROC_AUC = oob_auc,
              GRS_MDA = grs_mda,
              Age_MDA = age_mda,
              GRS_MDG = grs_mdg,
              Age_MDG = age_mdg,
              GRS_Pct_Importance = grs_pct,
              Age_Pct_Importance = age_pct,
              stringsAsFactors = FALSE
            ))
          }, error = function(e) {
            if (verbose) cat(sprintf("  ERROR: %s\n", e$message))
          })
        }
      }
    }
  }

  results
}

# Consistent HTML table renderer for report/slides
pretty_print_table <- function(df, caption = NULL) {
  knitr::kable(df, format = "html", escape = FALSE, caption = caption) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = TRUE,
      font_size = 12,
      position = "left"
    ) %>%
    kableExtra::scroll_box(width = "100%", height = "420px")
}