## This script sets up the functions used in this repository's scripts

# Function to run a logistic regression and compute a ROC AUC curve with 95% CIs

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
    roc_obj <- roc(data[[outcome]] ~ data$pred,
                   ci = TRUE)
  }
  
  print(roc_obj)
  
  return(list(
    model = logreg,
    data = data,
    outcome = outcome,
    predictor = predictor,
    covariates = covariates,
    roc = roc_obj
  ))
}

# Function to run a confusion matrix

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
                digits = 3) {
  
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
  
  # ── Cases ──────────────────────────────────────────────────────────────────
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
  
  # ── Controls ───────────────────────────────────────────────────────────────
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
  
  # ── Overall NRI ────────────────────────────────────────────────────────────
  NRI             <- NRI_cases + NRI_controls
  SE_NRI          <- sqrt(SE_NRI_cases^2 + SE_NRI_controls^2)
  CI_low_NRI      <- NRI - 1.96 * SE_NRI
  CI_high_NRI     <- NRI + 1.96 * SE_NRI
  z_NRI           <- NRI / SE_NRI
  p_NRI           <- 2 * pnorm(-abs(z_NRI))
  
  # ── Formatting helpers ─────────────────────────────────────────────────────
  fmt  <- function(x) formatC(x, format = "f", digits = digits)
  fmtp <- function(x) formatC(x, format = "g", digits = 3)
  pct  <- function(x) paste0(formatC(x * 100, format = "f", digits = 1), "%")
  
  ci_str <- function(est, lo, hi)
    paste0(fmt(est), " [", fmt(lo), ", ", fmt(hi), "]")
  
  # ── Results table ──────────────────────────────────────────────────────────
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
  
  print(results, row.names = FALSE, na.print = "")
  
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
