## Taking a look at fine mapping results

# Note: these fine mapping results were created using Robin Beaumont's "Finemapping_WGS” applet. Preparation of GWAS summary statistics for use with this tool can be found in the script "Conti_summary_stats_into_Fine_Mapping_ready_version"

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AllMen.cs")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AllMen.snp")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AllMen.summary")

dxdownload("Callum/FineMappingResults/UsingContiGWAS/Conti_FineMap_rs72725854_AllMen.cs")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/Conti_FineMap_rs72725854_AllMen.snp")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/Conti_FineMap_rs72725854_AllMen.summary")


dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_BlackOnly.cs")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_BlackOnly.snp")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_BlackOnly.summary")

dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AFR.cs")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AFR.snp")
dxdownload("Callum/FineMappingResults/UsingContiGWAS/FineMap_rs72725854_AFR.summary")

dxdownload("Callum/FineMappingResults/UsingWangGWAS/Wang_FineMap_rs72725854_BlackOnly.cs")
dxdownload("Callum/FineMappingResults/UsingWangGWAS/Wang_FineMap_rs72725854_BlackOnly.snp")
dxdownload("Callum/FineMappingResults/UsingWangGWAS/Wang_FineMap_rs72725854_BlackOnly.summary")

cs_generalGWAS <- read.delim("FineMap_rs72725854_AllMen.cs")
snp_generalGWAS <- read.delim("FineMap_rs72725854_AllMen.snp")
summary_generalGWAS <- read.delim("~/FineMap_rs72725854_AllMen.summary", comment.char="#")

cs_generalGWAS_repeated <- read.delim("Conti_FineMap_rs72725854_AllMen.cs")
snp_generalGWAS_repeated <- read.delim("Conti_FineMap_rs72725854_AllMen.snp")
summary_generalGWAS_repeated <- read.delim("~/Conti_FineMap_rs72725854_AllMen.summary", comment.char="#")

cs_BlackOnlyGWAS <- read.delim("FineMap_rs72725854_BlackOnly.cs")
snp_BlackOnlyGWAS <- read.delim("FineMap_rs72725854_BlackOnly.snp")
summary_BlackOnlyGWAS <- read.delim("~/FineMap_rs72725854_BlackOnly.summary", comment.char="#")

cs_AFRGWAS <- read.delim("FineMap_rs72725854_AFR.cs")
snp_AFRGWAS <- read.delim("FineMap_rs72725854_AFR.snp")
summary_AFRGWAS <- read.delim("~/FineMap_rs72725854_AFR.summary", comment.char="#")

cs_Wang_BlackOnlyGWAS <- read.delim("Wang_FineMap_rs72725854_BlackOnly.cs")
snp_Wang_BlackOnlyGWAS <- read.delim("Wang_FineMap_rs72725854_BlackOnly.snp")
summary_Wang_BlackOnlyGWAS <- read.delim("~/Wang_FineMap_rs72725854_BlackOnly.summary", comment.char="#")
