####################################################################
# Calculating HGDP+1KG principal components from phased haplotypes #
####################################################################

##### Note - this requires the working directory to have the publicly available phased haplotypes for HGDP+1KG. These can be downloaded at:
##### https://console.cloud.google.com/storage/browser/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2;tab=objects?prefix=&forceOnObjectsSortingFiltering=false 

## Step 1: convert bcf to plink2 files and merge all together

set -euo pipefail
THREADS=${THREADS:-8}   # default to 8 if not set

# Create a temporary PGEN with unique IDs

for chr in {1..22}; do
  inbcf="resources_hgdp_1kg_phased_haplotypes_v2_hgdp1kgp_chr${chr}.filtered.SNV_INDEL.phased.shapeit5.bcf"
  ./plink2 \
    --threads "${THREADS}" \
    --bcf "$inbcf" \
    --maf 0.01 \
    --snps-only just-acgt \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 50 truncate \
    --make-pgen \
    --out "tmp_chr${chr}"
done


# LD-prune per chromosome with PLINK 2

for chr in {1..22}; do
  ./plink2 \
    --threads "${THREADS}" \
    --pfile "tmp_chr${chr}" \
    --indep-pairwise 100 10 0.1 \
    --out "chr${chr}_pruned"
done


# Extract pruned variants to final per-chr PGENs
 
for chr in {1..22}; do
  ./plink2 \
    --threads "${THREADS}" \
    --pfile "tmp_chr${chr}" \
    --extract "chr${chr}_pruned.prune.in" \
    --make-pgen \
    --out "HGDP_1KGG_chr${chr}_pruned"
done


# Merge all autosomes in PGEN space (multi-allelic-safe)

: > pmerge_list.txt
for chr in {2..22}; do
  echo "HGDP_1KGG_chr${chr}_pruned" >> pmerge_list.txt
done

./plink2 \
  --threads ${THREADS} \
  --pfile HGDP_1KGG_chr1_pruned \
  --pmerge-list pmerge_list.txt \
  --make-pgen \
  --out HGDP_1KGG_autosomes_pruned

# (Optional) Export to PLINK 1 binary (BED) for downstream tools

./plink2 \
  --threads ${THREADS} \
  --pfile HGDP_1KGG_autosomes_pruned \
  --snps-only just-acgt \
  --max-alleles 2 \
  --make-bed \
  --out HGDP_1KGG_autosomes_pruned_p1




## Step 2: make a long-range LD mask, remove close relatives/duplicates


cat > longrange_ld_hg38.bed << 'EOF'
chr1  47761740 51761740
chr1  125169943 125170022
chr1  144106678 144106709
chr1  181955019 181955047
chr2  85919365 100517106
chr2  87416141 87416186
chr2  87417804 87417863
chr2  87418924 87418981
chr2  89917298 89917322
chr2  135275091 135275210
chr2  182427027 189427029
chr2  207609786 207609808
chr3  47483505 49987563
chr3  83368158 86868160
chr5  44464140 51168409
chr5  129636407 132636409
chr6  25391792 33424245
chr6  26726947 26726981
chr6  57788603 58453888
chr6  61109122 61357029
chr6  61424410 61424451
chr6  139637169 142137170
chr7  54964812 66897578
chr7  62182500 62277073
chr8  8105067 12105082
chr8  43025699 48924888
chr8  47303500 47317337
chr8  110918594 113918595
chr9  40365644 40365693
chr9  64198500 64200392
chr9  88958735 88959017
chr10 36671065 43184546
chr10 41693521 41885273
chr11 88127183 91127184
chr12 32955798 41319931
chr12 34639034 34639084
chr14 87391719 87391996
chr14 94658026 94658080
chr17 43159541 43159574
chr20 4031884 4032441
chr20 33948532 36438183
chr22 30060084 30060162
chr22 42980497 42980522
EOF


./plink2 \
  --pfile HGDP_1KGG_autosomes_pruned \
  --exclude bed1 longrange_ld_hg38.bed \
  --make-pgen \
  --out HGDP_1KGG_autosomes_pruned_masked

./plink2 \
  --pfile HGDP_1KGG_autosomes_pruned_masked \
  --king-cutoff 0.0884 \
  --make-pgen \
  --out HGDP_1KGG_autosomes_pruned_masked_norelatives

# Export to PLINK 1 binary (BED) for downstream tools (this is then used in Step 2 PCA projection code in R)

./plink2 \
  --threads ${THREADS} \
  --pfile HGDP_1KGG_autosomes_pruned_masked_norelatives \
  --snps-only just-acgt \
  --max-alleles 2 \
  --make-bed \
  --out HGDP_1KGG_autosomes_pruned_masked_norelatives_p1 # Final Plink .bed .bim .fam files used in step 2 (PCA projection) - upload these to UKB project space

# You can stop here if you want to then follow the rest of AncestryProbability1_bigsnpr


## Step 3: Compute Principal Components

./plink2 \
  --pfile HGDP_1KGG_autosomes_pruned_masked_norelatives \
  --pca 20 scols=sid \
  --out HGDP_1KGG_PCA_sid



## Step 4: Add ancestry labels


set -euo pipefail

# Inputs

PCA_EIG="HGDP_1KGG_PCA_sid.eigenvec"
META_TSV="release_3.1_secondary_analyses_hgdp_1kg_v2_metadata_and_qc_gnomad_meta_updated.tsv"

# Extract IID + PCs (IID PC1..PC20)

awk 'BEGIN{OFS=" "} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' \
    "$PCA_EIG" > pca.iid_pcs.txt

# Build (s, GROUP6) from population_inference.pop; collapse to Pan-UKB six groups; drop OCE

awk -F'\t' 'BEGIN{OFS=" "}
  NR==1{
    for(i=1;i<=NF;i++) h[$i]=i
    if(!("s" in h) || !("population_inference.pop" in h)){ print "ERROR: missing s or population_inference.pop" > "/dev/stderr"; exit 1 }
    sidx=h["s"]; pidx=h["population_inference.pop"]; next
  }
  {
    sub(/\r$/,"")
    id=$(sidx); tl=tolower($(pidx)); g="NA"
    if      (tl=="afr") g="AFR"
    else if (tl=="amr" || tl=="ami") g="AMR"
    else if (tl=="eas") g="EAS"
    else if (tl=="sas") g="CSA"
    else if (tl=="mid") g="MID"
    else if (tl=="nfe" || tl=="fin" || tl=="asj") g="EUR"
    else if (tl=="oce") g="OCE"
    if (g!="OCE") print id, g
  }' "$META_TSV" > meta.s_group6_raw.txt

# Merge by IID (hash-join; no sorting required); keep labeled samples only

awk 'NR==FNR{lab[$1]=$2; next} ($1 in lab){print $0, lab[$1]}' \
    meta.s_group6_raw.txt pca.iid_pcs.txt > HGDP_1KGG_PCA_labeled_six.txt

# Quick sanity: counts per group + spot-check a known European sample (HG00096)

awk '{c[$NF]++} END{for(k in c) print k, c[k]}' HGDP_1KGG_PCA_labeled_six.txt | sort
grep -m1 '^HG00096 ' HGDP_1KGG_PCA_labeled_six.txt || true

#### NOTE: The final file containing HGDP+1KG Principal Components plus ancestry labels is called "HGDP_1KGG_PCA_labeled_six.txt"