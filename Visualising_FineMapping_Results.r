## Taking a look at fine mapping results

# Note: these fine mapping results were created using Robin Beaumont's "Finemapping_WGS” applet. Preparation of GWAS summary statistics for use with this tool can be found in the script "Conti_summary_stats_MASTER_into_Fine_Mapping_ready_version"

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

dxdownload("Callum/FineMappingResults/FineMap_rs72725854_test.cs")
dxdownload("Callum/FineMappingResults/FineMap_rs72725854_test.snp")
dxdownload("Callum/FineMappingResults/FineMap_rs72725854_test.summary")

dxdownload("Callum/FineMappingResults/FineMap_rs72725854_AFR.cs")
dxdownload("Callum/FineMappingResults/FineMap_rs72725854_AFR.snp")
dxdownload("Callum/FineMappingResults/FineMap_rs72725854_AFR.summary")

cs_generalGWAS <- read.delim("FineMap_rs72725854_test.cs")
snp_generalGWAS <- read.delim("FineMap_rs72725854_test.snp")
summary_generalGWAS <- read.delim("~/FineMap_rs72725854_test.summary", comment.char="#")

cs_AFRGWAS <- read.delim("FineMap_rs72725854_AFR.cs")
snp_AFRGWAS <- read.delim("FineMap_rs72725854_AFR.snp")
summary_AFRGWAS <- read.delim("~/FineMap_rs72725854_AFR.summary", comment.char="#")
