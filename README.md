## This repository contains Callum's scripts for assessing Prostate Cancer GRS performance in Black Men in the UK Biobank (UKB) and Our Future Health (OFH).

If you are here to reproduce the analyses from the paper "Current Prostate Cancer Genetic Risk Scores Add Limited Predictive Value for Black Men: a prospective cohort study in the UK Biobank and Our Future Health", take the following steps:

### Step 1: Reproduce UK Biobank Results

1. Run `AncestryProbability2calculation.R`. This script requires downloads from online sources, detailed at the relevant place in the script.
2. Run `Conti_GRS_267.r`, `Create_adjustedGRS_weighting_Conti_odds_ratios_by_AncestryProbability2.R`, and `Conti_script.R` to generate Prostate Cancer GRS using 3 different methods. Produce all GRS detailed at the top of `SummaryOfLogisticRegressions.Rmd`.
3. Run `Extracting_rs72725854_carriers.R` to generate a dataframe of rs72725854 carriers (or comment out relevant dx download calls for these files in `SummaryOfLogisticRegressions.Rmd`, since it's never actually used).
4. Run `Family_history.R` to generate a dataframe of prostate/breast cancer family history.
5. Using the UKB-RAP's cohort browser, create the datasets described at the top of `SummaryOfLogisticRegressions.Rmd`.
6. Knit `SummaryOfLogisticRegressions.Rmd` to html and save to persistent storage.

### Step 2: Reproduce Our Future Health Results

1. Run `GRS_script.ipynb`. This requires that you upload supplementary table 4 from both Conti et al. (2021) (link: https://www.nature.com/articles/s41588-020-00748-0) and supplementary table 4 from Wang et al. (2023) (https://www.nature.com/articles/s41588-023-01534-4) as .csv files with the first ~4 rows trimmed, at that you use a polygenic scoring script from lcpilling. If you do not have access to this script, please email L.Pilling@exeter.ac.uk.
2. Run `Logistic_Regressions.ipynb`. Tip: use a higher instance for this (at least 32GB RAM).
3. Run `Age_at_Diagnosis.ipynb`.

### Step 3: Combine results

1. Export results via the OFH (and presumably the UKB, when access resumes) airlock to an unrestricted project.
2. Download the results files to a local machine or server.
3. Run `PrCa_paper.Rmd` on the local machine or server, altering file paths as required

#### Note: this may not work immediately when access to the UK Biobank resumes, as we do not yet know the nature of the airlock. For example, rendered html files may not be permitted for export. These instructions and scripts will be adjusted accordingly once exports resume.
