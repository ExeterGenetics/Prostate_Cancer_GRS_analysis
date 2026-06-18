#!/bin/bash

#### In this script I am going to be running a Step 1 for the Regenie GWAS pipe line using swiss army knife. 
##### This is going to be run for Haemochromatosis



#### Luke's test line to run code

# bash test_script.sh -output_prefix="hh_bin_GWAS" -pheno_file="/users/luke_s/GWAS/haemochromatosis/phenotype/hh_diagnosis.pheno" -covar_file="/data/genetics/GWAS/GWAS_covariates_imputed_noquote.white.covar" -rint="no" -catcovars="region_at_reg" -contcovars="age,age_sqd,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" -output_dir="/users/luke_s/GWAS/test_step1_dir/" -bin_or_cont="BIN" -bfile="/data/genetics/golden/ofh_snv.v9.golden.pruned"

#### Callum's test line to run code

# bash regenie_step_1.sh -output_prefix="PrCa_EUR_GWAS" -pheno_file="/users/callum/GWAS/PrCa_EUR_regenie.txt" - covar_file="/users/callum/GWAS/ofh_covars_regenie_males_EUR.txt" -rint="no" -catcovars="region_at_reg" -contcovars="age,age_sqd,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" -output_dir="/users/callum/GWAS/step_1_dir/" -bin_or_cont="BIN" - bfile="/data/genetics/golden/ofh_snv.v10.golden_EUR"

for argument in "$@"
do
    case $argument in 
        -output_prefix=*)
            OUTPUT_PREFIX_FULL="${argument#*=}"
            ;;
        -pheno_file=*)
            PHENO_FILE="${argument#*=}"
            ;;
        -catcovars=*)
            CAT_COVARS="${argument#*=}"
            ;;
        -contcovars=*)
            CONT_COVARS="${argument#*=}"
            ;;
        -covar_file=*)
            COVAR_FILE="${argument#*=}"
            ;;
        -output_dir=*)
            OUTPUT_DIRECTORY="${argument#*=}"
            ;;
        -bin_or_cont=*)
            BIN_OR_CONT="${argument#*=}"
            ;;
        -rint=*)
            RINT="${argument#*=}"
            ;;
        -bfile=*)
            BFILE_INPUT="${argument#*=}"
            ;;
        *)
            echo "Unknown argument: $argument"
            ;;
    esac
done


##### Show all inputs to user.
echo "Output prefix = ${OUTPUT_PREFIX_FULL}"
echo "Phenotype file path = ${PHENO_FILE}"
echo "Covariate file path = ${COVAR_FILE}"
echo "Rank Inverse Normal Transformation = ${RINT}"
echo "Output directory path = ${OUTPUT_DIRECTORY}"
echo "Binary or Continuous = ${BIN_OR_CONT}"
echo "Bed file for step 1 = ${BFILE_INPUT}"
echo "Categorical Covariates = ${CAT_COVARS}"
echo "Continuous Covariates = ${CONT_COVARS}"



##### Ask user if inputs are correct and continue with script if they answer yes.
read -p "Are the values you inputted correct? (y/n): " answer

if [[ "$answer" == "y" || "$answer" == "Y" ]]; then
  echo "Continuing script..."
elif [[ "$answer" == "n" || "$answer" == "N" ]]; then
    echo "Exiting script"
    exit 1
else
  echo "INVALID INPUT. ONLY PUT y OR n"
  exit 1
fi

######## THERE IS A PROBLEM WITH THIS AND IM NOT SURE WHAT IT IS !!!!!!
##### Creating a variable with relevant flag depending on BIN or CONT
BIN_OR_CONT_VAR=""

if [[ "$BIN_OR_CONT" == "BIN" ]]; then
    BIN_OR_CONT_VAR="--bt"
elif [[ "$BIN_OR_CONT" == "CONT" ]]; then
    if [["$RINT" == "yes"]]; then
        BIN_OR_CONT_VAR="--apply-rint"
    else
        BIN_OR_CONT_VAR=""
    fi
else
    echo "INVALID INPUT MUST BE BIN OR CONT FOR $BIN_OR_CONT"
    exit 1
fi


##### Get file ID for the phenotype file so can mount to SAK. 
PHENO_FILE_FILE_NAME=$(dx ls -l "$PHENO_FILE" | awk '{print $7}' | sed 's/(//g' | sed 's/)//g')
echo "Phenotype File Name"
echo "$PHENO_FILE_FILE_NAME"

##### Get file ID for the covariate file so can mount to SAK.
COVAR_FILE_FILE_NAME=$(dx ls -l "$COVAR_FILE" | awk '{print $7}' | sed 's/(//g' | sed 's/)//g')
echo "Covariate File Name"
echo "$COVAR_FILE_FILE_NAME"

##### Get file IDs for the bed, bim and fam files for the pruned high quality directly sequenced array bases. 
BFILE="${BFILE_INPUT}" 
BFILE_FILE_NAME_BED=$(dx ls -l "${BFILE}.bed" | awk '{print $7}' | sed 's/(//g' | sed 's/)//g')
BFILE_FILE_NAME_BIM=$(dx ls -l "${BFILE}.bim" | awk '{print $7}' | sed 's/(//g' | sed 's/)//g')
BFILE_FILE_NAME_FAM=$(dx ls -l "${BFILE}.fam" | awk '{print $7}' | sed 's/(//g' | sed 's/)//g')
echo "Golden SNP variant files"
echo "$BFILE_FILE_NAME_BED"
echo "$BFILE_FILE_NAME_BIM"
echo "$BFILE_FILE_NAME_FAM"


##### $CMD will be passed into the SAK.
##### I will cover what happens in this section here split by paragraph:

##### Paragraph 1:
##### Resave each of the inputs in the command so they have similar names in the SAK.

##### Paragraph 2:
##### As files have been mounted into the instance they don't use the full file path, so I removed the first
##### part of the file path.

##### Paragraph 3:
##### This is the regenie command that runs step 1.
##### PROBLEM: I NEED TO MAKE IT SO THEY HAVE TO ENTER THE COVARIATES PRESENT IN THE COVARIATE FILE.
##### (This shouldn't effect us too much as we are basically using the same covariate files) 

##### Paragraph 4
##### This sed command removes the file path that is before the loco file in the pred.list file. If it contains this file path
##### It cannot find the loco file in step 2.

### This seems to work but I don't fully understand why this worked and the otherone didnt. 

CMD='
PHENO_FILE="'"${PHENO_FILE}"'"
OUT_PREFIX_FULL="'"${OUTPUT_PREFIX_FULL}"'"
COVAR_FILE="'"${COVAR_FILE}"'"
BFILE="'"${BFILE_INPUT}"'"
BIN_STEP1="'"${BIN_OR_CONT_VAR}"'"
CATS="'"${CAT_COVARS}"'"
CONTS="'"${CONT_COVARS}"'"

PHENO_FILE=$(echo "${PHENO_FILE##*/}")
COVAR_FILE=$(echo "${COVAR_FILE##*/}")
BFILE=$(echo "${BFILE##*/}")

regenie \
--step 1 \
--bed "${BFILE}" \
--phenoFile "${PHENO_FILE}" \
--covarFile "${COVAR_FILE}" ${BIN_STEP1} \
--covarColList "${CONTS}" \
--catCovarList "${CATS}" \
--maxCatLevels 30 \
--lowmem \
--bsize 1000 \
--out "${OUT_PREFIX_FULL}"
sed -i 's@/home/dnanexus/out/out/@@g' "${OUT_PREFIX_FULL}_pred.list"
'

echo "$CMD"

##### Below is the Swiss Army Knife (SAK) command. I input the CMD I wrote above and mount the phenotype file,
##### the covariate file, and the bed, bin and fam files for the golden variants.

dx run swiss-army-knife \
        -icmd="$CMD" \
        -iin="$PHENO_FILE_FILE_NAME" \
        -iin="$COVAR_FILE_FILE_NAME" \
        -iin="$BFILE_FILE_NAME_BED" \
        -iin="$BFILE_FILE_NAME_BIM" \
        -iin="$BFILE_FILE_NAME_FAM" \
        -imount_inputs=true \
        --tag regenie_step1 \
        --instance-type="azure:mem2_ssd1_v2_x8" \
        --destination "${OUTPUT_DIRECTORY}" \
        --priority high \
        --yes



