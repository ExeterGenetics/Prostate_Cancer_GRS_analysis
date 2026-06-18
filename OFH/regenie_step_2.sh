
#!/bin/bash


##### Final changes I need to add :
### Make it so you have to enter the covariate file headers as continuous or binary
### Make it so you add in your own extract file with the variants you want to look at. 

#### Luke's test line to run code

# bash test_script_step2.sh -pheno_file="/users/luke_s/GWAS/haemochromatosis/phenotype/hh_diagnosis.pheno" -output_prefix="hh_bin_GWAS_step2" -catcovars="sex,region_at_reg" -contcovars="age_at_assessment,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" -covar_file="/data/genetics/GWAS/GWAS_covariates_imputed_noquote.white.covar" -info="0.3" -maf="0.001" -rint="no" -output_dir="/users/luke_s/GWAS/test_step2_dir/" -bin_or_cont="BIN" -loco="/users/luke_s/GWAS/test_step1_dir/hh_bin_GWAS_1.loco" -predlist="/users/luke_s/GWAS/test_step1_dir/hh_bin_GWAS_pred.list"

#### Callum's test line to run code

# bash regenie_step_2.sh -pheno_file="/users/callum/GWAS/PrCa_EUR_regenie.txt" -output_prefix="PrCa_EUR_GWAS_step2" -catcovars="region_at_reg" -contcovars="age,age_sqd,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" -covar_file="/users/callum/GWAS/ofh_covars_regenie_males_EUR.txt" -info="0.3" -maf="0.001" -rint="no" -output_dir="/users/callum/GWAS/step_2_dir/" -bin_or_cont="BIN" -loco="/users/callum/GWAS/step_1_dir/PrCa_EUR_GWAS_1.loco" -predlist="/users/callum/GWAS/step_1_dir/PrCa_EUR_GWAS_pred.list"

for argument in "$@"
do
    case $argument in 
        -output_prefix=*)
            PREFIX_STEP2="${argument#*=}"
            ;;
        -pheno_file=*)
            FILE="${argument#*=}"
            ;;
        -catcovars=*)
            CAT_COVARS="${argument#*=}"
            ;;
        -contcovars=*)
            CONT_COVARS="${argument#*=}"
            ;;
        -info=*)
            INFO="${argument#*=}"
            ;;
        -maf=*)
            MAF="${argument#*=}"
            ;;
        -rint=*)
            RINT="${argument#*=}"
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
        -loco=*)
            STEP1_LOCO="${argument#*=}"
            ;;
        -predlist=*)
            STEP1_PRED_LIST="${argument#*=}"
            ;;
        *)
            echo "Unknown argument: $argument"
            ;;
    esac
done


##### Show all inputs to user.
echo "Output prefix = ${PREFIX_STEP2}"
echo "Phenotype file path = ${FILE}"
echo "Categorical Covariates = ${CAT_COVARS}"
echo "Continuous Covariates = ${CONT_COVARS}"
echo "Minimum INFO = ${INFO}"
echo "Minimum MAF = ${MAF}"
echo "Rank Inverse Normal Transformation = ${RINT}"
echo "Covariate file path = ${COVAR_FILE}"
echo "Output directory path = ${OUTPUT_DIRECTORY}"
echo "Binary or Continuous phenotype = ${BIN_OR_CONT}"
echo "Loco file from step 1 = ${STEP1_LOCO}"
echo "Pred.list file from step 1 = ${STEP1_PRED_LIST}"



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

echo "${BIN_OR_CONT_VAR[@]}"

#### Here I am going to convert the MAF into a MAC using the number of individuals in the covariate file.

#### Should I be using the covariate file or the phenotype file or should I be finding out how many of teh peole are similar between teh two files.
### I think covariates might be good enough you know.
NUM_INV=$(tail -n+2 "/mnt/project/${COVAR_FILE}" | wc -l )
echo "${NUM_INV}"
MAC=$(awk -v n="$NUM_INV" -v maf="$MAF" 'BEGIN { print n * 2 * maf }')
echo "Raw MAC = ${MAC}"
MAC=$(printf "%.0f" "$MAC")
echo "Rounded MAC = ${MAC}"

###### This section find the file ID from the supplied pred.list file
STEP1_PRED_LIST_FILE_NAME=$(dx ls -l "$STEP1_PRED_LIST" | awk '{print $7}' | sed 's/(//g' | sed 's/)//g')
echo "PRED LIST FILE: ${STEP1_PRED_LIST_FILE_NAME}"

###### This section find the file ID from the supplied LOCO file
STEP1_LOCO_FILE_NAME=$(dx ls -l "$STEP1_LOCO" | awk '{print $7}' | sed 's/(//g' | sed 's/)//g')
echo "LOCO FILE: ${STEP1_LOCO_FILE_NAME}"

###### This section find the file ID from the supplied phenotype file
FILE_FILE_NAME=$(dx ls -l "$FILE" | awk '{print $7}' | sed 's/(//g' | sed 's/)//g')
echo "PHENOTYPE FILE: ${FILE_FILE_NAME}"

###### This section find the file ID from the supplied covariate file
COVAR_FILE_FILE_NAME=$(dx ls -l "$COVAR_FILE" | awk '{print $7}' | sed 's/(//g' | sed 's/)//g')
echo "COVAR FILE: ${COVAR_FILE_FILE_NAME}"

### In this section I am going to create the list of bgen files then I am going to group them into groups of 15.
dx ls -l /imputed_bgen/ | awk '{print $6,$7}'| grep '.bgen ' | sort | awk '{print $2}' | sed 's/(//g' | sed 's/)//g' | split -l 15 --numeric-suffixes=1 --additional-suffix=.txt - bgen_group_id_

### In this setion I am going to pull out a list of the the index files and group them into groups of 15.
dx ls -l /imputed_bgen/ | awk '{print $6,$7}'| grep '.bgen.bgi' | sort | awk '{print $2}' | sed 's/(//g' | sed 's/)//g' | split -l 15 --numeric-suffixes=1 --additional-suffix=.txt - bgi_group_id_

#### In this section I am going to remove all bit 1 of teh bgen and bgi groups
mv bgen_group_id_01.txt test_group1
rm bgen_group_id_*
mv test_group1 bgen_group_id_01.txt
#
mv bgi_group_id_01.txt test_bgi_group1
rm bgi_group_id_*
mv test_bgi_group1 bgi_group_id_01.txt


##### In this section I will be making the regenie command that will be going into the Swiss army knife.
##### I will describe what each paragraph does here.
##### Paragraph 1:
##### Here I save the file paths for the key files and important information to variables.

##### Paragraph 2:
##### Here I remove the start of the file paths. As they are mounted to the SAK instance including the file path will break it.

##### Paragraph 3:
##### This is a for loop that loop through the the 15 bgen files in the SAK.

##### Paragraph 4:
##### Here I print the name of the bgen that is about to be analysed. Then I create a variable of teh bgen name so I can save
##### it to the output.

##### Paragraph 5:
##### This is the regenie command that runs step 2 of the GWAS.



######## THERE IS A PROBLEM WITH BIN_STEP2 AS IT IS ONLY ADDING --bt AND NOT ALL THE OTHER STUFF.

#### NEED TO ADD IN XARGS

CMD='

STEP1_PRED_LIST="'"${STEP1_PRED_LIST}"'"
STEP1_LOCO="'"${STEP1_LOCO}"'" 
TEST="additive"
FILE="'"${FILE}"'"
COVAR_FILE="'"${COVAR_FILE}"'"
PREFIX_STEP2="'"${PREFIX_STEP2}"'"
BIN_STEP2=("'"${BIN_OR_CONT_VAR[@]}"'")
MAC="'"${MAC}"'"
INFO="'"${INFO}"'"
CATS="'"${CAT_COVARS}"'"
CONTS="'"${CONT_COVARS}"'"

echo ${MAC}
echo ${INFO}

STEP1_PRED_LIST=$(echo "${STEP1_PRED_LIST##*/}")
STEP1_LOCO=$(echo "${STEP1_LOCO##*/}")
FILE=$(echo "${FILE##*/}")
COVAR_FILE=$(echo "${COVAR_FILE##*/}")

for line in *.bgen 
do 
    echo "$line"

    LOOP_FILE_NAME=$(echo "$line" | sed 's/.bgen//g')

    echo "$LOOP_FILE_NAME"

    regenie \
	  --step 2 \
	  --bgen ${line} \
      --bgi "${line}.bgi" \
	  --ref-first \
      ${BIN_STEP2[@]} \
      --minINFO ${INFO} \
      --minMAC ${MAC} \
	  --sample ofh_imputed.v6.chr1-b0001.sample \
	  --phenoFile "${FILE}" \
	  --covarFile "${COVAR_FILE}" \
	  --pred "${STEP1_PRED_LIST}" \
	  --covarColList "${CONTS}" \
	  --catCovarList "${CATS}" \
	  --maxCatLevels 30 \
	  --test "${TEST}" \
	  --bsize 400 \
	  --out ${PREFIX_STEP2}_${LOOP_FILE_NAME}

done
'
echo $CMD

#### END OF CMD


##### Here I am going to put the Swiss army knife command. 

#### file-J8PKFP12ykG8K0151xz6k33G is the b0001 sample file. I have checked and I think that all of these are the same. 

### Getting the bgis in might be a little bit more complicated then I previously thought. 
for bgen_group_loop in bgen_group_id_*
do
    echo "$bgen_group_loop"

    bgi_group_loop=$(echo "$bgen_group_loop" | sed 's/bgen/bgi/g')
    echo "${bgi_group_loop}"


    dx run swiss-army-knife \
        $(cat "$bgen_group_loop" | sed 's/^/-iin=/') \
        $(cat "$bgi_group_loop" | sed 's/^/-iin=/') \
        -iin=file-J8PKFP12ykG8K0151xz6k33G \
        -iin="$STEP1_PRED_LIST_FILE_NAME" \
        -iin="$FILE_FILE_NAME" \
        -iin="$COVAR_FILE_FILE_NAME" \
        -iin="$STEP1_LOCO_FILE_NAME" \
        -imount_inputs=true \
        -icmd="$CMD" \
        --tag regenie_step2 \
        --instance-type="azure:mem2_ssd1_v2_x8" \
        --priority high \
        --destination "${OUTPUT_DIRECTORY}" \
        --brief \
        --yes

done

  