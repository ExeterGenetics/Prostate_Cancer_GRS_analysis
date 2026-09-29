**This repository includes Callum's scripts for assessing Prostate Cancer GRS performance in Black Men in the UK Biobank.**

Key:

- "**SummaryOfLogisticRegressions.R**" is the main analysis pipeline summarised as an Rmd, which includes instructions for how to derive all relevant datasets using the UKB-RAP cohort browser and/or other scripts in the repository.

- The folder "**AncestryProbability1_bigsnpr**" contains 1 shell script (ideally to run on Slade) and 2 subsequent R scripts for:
  
  i) Calculating Principal Components for the Human Genome Diversity Project (HGDP) + 1000 Genomes project participants

  ii) Projecting UK Biobank participants into the Principal Components space using dispensed PLINK files (code: ukb22418), using the R package bigsnpR

  iii) Using a Random Forest model to calculate each participant's probability of falling into the major ancestry-based clusters of HGDP+1KG

  iv) Saving this output as "AncestryProbability", which is used in the main analysis pipeline

- The folder "**AncestryProbability2_plink**" contains an R script for doing something similar. If in doubt, use this one over **AncestryProbability1_bigsnpr** as it is more reproducible and more faithful to the Pan-UKB study which first derived the "Genomic Ancestry" (Genetic Similarity) variable in the UK Biobank (https://github.com/atgu/ukbb_pan_ancestry/tree/master):

  i) Downloading pre-computed Principal Components and the corresponding Allele Frequencies and Loadings for the Human Genome Diversity Project (HGDP) + 1000 Genomes project participants

  ii) Projecting UK Biobank participants into the Principal Components space using imputed genotypes based on these pre-calculated PCs (code: ukb22828), using plink/plink2

  iii) Using a Random Forest model to calculate each participant's probability of falling into the major ancestry-based clusters of HGDP+1KG

  iv) Saving this output as "AncestryProbability2", which is used in the main analysis pipeline

- Conti_script.R calculates Conti's Genetic Risk Score for UKB participants, using lcpilling's ukbrapR::create.pgs() function. This script was produced by @hdg204 (https://github.com/hdg204).

- Create_adjustedGRS_weighting_Conti_odds_ratios_by_AncestryProbability2.R uses the Ancestry Probability (calculated using method 2) variable for computing an adjusted GRS, where odds ratios of SNPs in the GRS are altered depending on probabilities of fitting into the ancestry groups of the HGDP+1KG reference panel.

- Extracting_rs72725854_carriers.R identifies all participants in the biobank with a recorded rs72725854 SNP (A>G or A>T) and adds them to a dataframe. This script was produced by @hdg204 (https://github.com/hdg204).

- Family_history.R is an adapted version of @hdg204's Family History Script (https://github.com/hdg204/UKBB/blob/main/Family_History_Script.R) which collates Prostate and Breast cancer family history into a more useful binary format that can be applied in models.

- Scatter_plot_HGDP1KG_projected_PCs.R and Scatter_plot_UKB_provided_PCs.R are simple scatter plots for visualising Principal Components from the HGDP+1KG reference panel and the UK Biobank built-in Principal Components.

- Stratifying_risk_by_GRS_percentile.R and Stratifying_risk_by_IRM_pred_percentile.R produce plots to visualise how certain quantiles of GRSs and Integrated Risk Models differ in their prevalence of Prostate Cancer diagnosis.
