### YOU WILL NEED TO RUN THESE TWO LINES SEPARETELY (AND INDIVIDUALLY), THEN THE REST OF THE SCRIPT
remotes::install_github("lcpilling/ukbrapR@v0.3.10")

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')

library(dplyr)
library(readxl)
library(stringr)
library(tidyr)

## GRSs trained to predict general prostate cancer diagnosis

system('dx download Callum/ContiGWAS/Conti2021supplementarytables.xlsx') # This is the supplementary table file from Conti et al. (2021), downloadable at: https://www.nature.com/articles/s41588-020-00748-0
system('dx download Callum/WangGWAS/Wang2023supplementarytables.xlsx') # This is the supplementary table file from Wang et al. (2023), downloadable at: https://pmc.ncbi.nlm.nih.gov/articles/PMC10841479/
system('dx download Callum/SchumacherGWAS/Schumacher.txt') # This is the list of SNPs and weights from Schumacher et al. (2018), downloadable at: https://www.pgscatalog.org/publication/PGP000019/ 
system('dx download Callum/SchumacherGWAS/BARCODE1.txt') # This is the list of SNPs and weights from BARCODE1 (2021), downloadable at: https://www.pgscatalog.org/publication/PGP000726/ 

## GRSs trained to predict aggressive prostate cancer diagnosis

system('dx download Callum/SeibertGWAS/Seibert.txt') # This is the list of SNPs and weights from Seibert et al. (2018), downloadable at: https://www.pgscatalog.org/publication/PGP000047/ 
system('dx download Callum/SeibertGWAS/Pagadala.txt') # This is the list of SNPs and weights from Pagadala et al. (2022), downloadable at: https://www.pgscatalog.org/publication/PGP000400/

############################
# Define withdrawal filter #
############################

exclude_withdrawn=function(df){
  system('dx download Callum/Withdrawals/withdrawn_20260310.csv --overwrite') ## This file is a list of participants who withdrew from the Biobank up to the date 10th March 2026. This was sent from the UK Biobank team via email to members of approved applications
  df2 = df %>% left_join(
    read_csv("withdrawn_20260310.csv", col_names = FALSE, show_col_types = FALSE) %>%
      mutate(w=1) %>%
      rename(eid = X1),
    by='eid'
  ) %>%
    filter(is.na(w))
  return(df2)
}

#####################################
# Option #1 - Conti Multiethnic GRS #
#####################################


# The below block of code reads the conti xlsx file, and puts it into a form that ukbrapR can use to make a GRS (the outputted tsv)

#If you change the effect_weight column it should be easy to run a different GRS

Conti <- read_excel("Conti2021supplementarytables.xlsx", sheet = "S4", skip = 3, na = "NA")
Conti=Conti[1:269,]
Conti=arrange(Conti,Chromosome,Position)
Conti=Conti%>%rename(
  rsID=`rs*`,
  CHR=Chromosome,
  POS=Position,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`Multiethnic Analysis`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Conti2=Conti%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Conti2$CHR[is.na(Conti2$CHR)]="X"
write.table(Conti2,'Conti.tsv',quote=FALSE,sep='\t',row.names = FALSE)


conti_out=ukbrapR:::create_pgs(
  in_file='Conti.tsv',
  out_file='Conti_multi_ethnic.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Conti_multi_ethnic.pgs.bim')
#outscore=read.table('Conti_multi_ethnic.pgs.profile',header=T)

GRS <- read.table("Conti_multi_ethnic.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "Conti_multi_ethnic.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Conti_multi_ethnic.pgs.tsv"))



###########################################
# Option #2 - Conti European-specific GRS #
###########################################

Conti <- read_excel("Conti2021supplementarytables.xlsx", sheet = "S4", skip = 3, na = "NA")
Conti=Conti[1:269,]
Conti=arrange(Conti,Chromosome,Position)
Conti=Conti%>%rename(
  rsID=`rs*`,
  CHR=Chromosome,
  POS=Position,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`European...16`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Conti2=Conti%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Conti2$CHR[is.na(Conti2$CHR)]="X"
write.table(Conti2,'Conti.tsv',quote=FALSE,sep='\t',row.names = FALSE)


conti_out=ukbrapR:::create_pgs(
  in_file='Conti.tsv',
  out_file='Conti_European.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Conti_European.pgs.bim')
#outscore=read.table('Conti_European.pgs.profile',header=T)

GRS <- read.table("Conti_European.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "Conti_European.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Conti_European.pgs.tsv"))



##########################################
# Option #3 - Conti African-specific GRS #
##########################################

Conti <- read_excel("Conti2021supplementarytables.xlsx", sheet = "S4", skip = 3, na = "NA")
Conti=Conti[1:269,]
Conti=arrange(Conti,Chromosome,Position)
Conti=Conti%>%rename(
  rsID=`rs*`,
  CHR=Chromosome,
  POS=Position,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`African...19`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Conti2=Conti%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Conti2$CHR[is.na(Conti2$CHR)]="X"
write.table(Conti2,'Conti.tsv',quote=FALSE,sep='\t',row.names = FALSE)


conti_out=ukbrapR:::create_pgs(
  in_file='Conti.tsv',
  out_file='Conti_African.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Conti_African.pgs.bim')
#outscore=read.table('Conti_African.pgs.profile',header=T)

GRS <- read.table("Conti_African.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "Conti_African.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Conti_African.pgs.tsv"))




#############################################
# Option #4 - Conti East Asian-specific GRS #
#############################################

Conti <- read_excel("Conti2021supplementarytables.xlsx", sheet = "S4", skip = 3, na = "NA")
Conti=Conti[1:269,]
Conti=arrange(Conti,Chromosome,Position)
Conti=Conti%>%rename(
  rsID=`rs*`,
  CHR=Chromosome,
  POS=Position,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`East Asian...22`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Conti2=Conti%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Conti2$CHR[is.na(Conti2$CHR)]="X"
write.table(Conti2,'Conti.tsv',quote=FALSE,sep='\t',row.names = FALSE)


conti_out=ukbrapR:::create_pgs(
  in_file='Conti.tsv',
  out_file='Conti_East_Asian.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Conti_East_Asian.pgs.bim')
#outscore=read.table('Conti_East_Asian.pgs.profile',header=T)

GRS <- read.table("Conti_East_Asian.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "Conti_East_Asian.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Conti_East_Asian.pgs.tsv"))


###########################################
# Option #5 - Conti Hispanic-specific GRS #
###########################################

Conti <- read_excel("Conti2021supplementarytables.xlsx", sheet = "S4", skip = 3, na = "NA")
Conti=Conti[1:269,]
Conti=arrange(Conti,Chromosome,Position)
Conti=Conti%>%rename(
  rsID=`rs*`,
  CHR=Chromosome,
  POS=Position,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`Hispanic...25`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Conti2=Conti%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Conti2$CHR[is.na(Conti2$CHR)]="X"
write.table(Conti2,'Conti.tsv',quote=FALSE,sep='\t',row.names = FALSE)


conti_out=ukbrapR:::create_pgs(
  in_file='Conti.tsv',
  out_file='Conti_Hispanic.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Conti_Hispanic.pgs.bim')
#outscore=read.table('Conti_Hispanic.pgs.profile',header=T)

GRS <- read.table("Conti_Hispanic.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "Conti_Hispanic.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Conti_Hispanic.pgs.tsv"))



#=====================================#
# -- New section: Wang (2023) GRSs -- #
#=====================================#

# First, we need to rename the columns and ensure they are numeric

Wang <- read_excel("Wang2023supplementarytables.xlsx", sheet = "S4", skip = 3, na = "NA") 

Wang <- Wang %>%
  dplyr::rename(
    EUR_Rsquared = European,
    AFR_Rsquared = African,
    EAS_Rsquared = Asian,          
    HIS_Rsquared = Hispanic,       
    OR_Multiethnic_Marginal = `OR...15`,
    CI95_Multiethnic_Marginal = `95%CI...16`,
    Pvalue_Multiethnic_Marginal = `P-value...17`,
    OR_Multiethnic_Conditional = `OR...18`,
    CI95_Multiethnic_Conditional = `95%CI...19`,
    Pvalue_Multiethnic_Conditional = `P-value...20`,
    RAF_EUR = `RAF...21`,
    OR_EUR = `OR...22`,
    CI95_EUR = `95%CI...23`,
    Pvalue_EUR = `P-value...24`,
    RAF_AFR = `RAF...25`,
    OR_AFR = `OR...26`,
    CI95_AFR = `95%CI...27`,
    Pvalue_AFR = `P-value...28`,
    RAF_EAS = `RAF...29`,
    OR_EAS = `OR...30`,
    CI95_EAS = `95%CI...31`,
    Pvalue_EAS = `P-value...32`,
    RAF_HIS = `RAF...33`,
    OR_HIS = `OR...34`,
    CI95_HIS = `95%CI...35`,
    Pvalue_HIS = `P-value...36`
  ) %>%
  # Clean P-values like "<1e-5" -> "1e-5" then numeric
  mutate(
    across(
      c(Pvalue_Multiethnic_Marginal, Pvalue_Multiethnic_Conditional,
        Pvalue_EUR, Pvalue_AFR, Pvalue_EAS, Pvalue_HIS),
      ~ .x |> as.character() |> str_replace("^\\s*<\\s*", "") |> as.numeric()
    ),
    across(
      c(EUR_Rsquared, AFR_Rsquared, EAS_Rsquared, HIS_Rsquared,
        OR_Multiethnic_Marginal, OR_Multiethnic_Conditional,
        RAF_EUR, OR_EUR, RAF_AFR, OR_AFR, RAF_EAS, OR_EAS, RAF_HIS, OR_HIS),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) 


####################################
# Option #6 - Wang Multiethnic GRS #
####################################


## Filter to Wang 451-SNP GRS-specific variants, and (optionally) to variants with p value below 0.05

Wang2 <- Wang %>%
  dplyr::filter(!is.na(`Wang et al., 451 SNPs`)) #%>%
#dplyr::filter(Pvalue_Multiethnic_Marginal <= 0.05)

## Prepare for bgenix GRS calculation

Wang2=arrange(Wang2,Chromosome,`Position (GRCh37)`)
Wang2=Wang2%>%rename(
  rsID=rsID,
  CHR=Chromosome,
  POS=`Position (GRCh37)`,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`OR_Multiethnic_Marginal`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Wang3=Wang2%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Wang3$CHR[is.na(Wang3$CHR)]="X"
write.table(Wang3,'Wang.tsv',quote=FALSE,sep='\t',row.names = FALSE)


GRS_out=ukbrapR:::create_pgs(
  in_file='Wang.tsv',
  out_file='Wang_multi_ethnic.pgs',
  pgs_name='Wang',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Wang_multi_ethnic.pgs.bim')
#outscore=read.table('Wang_multi_ethnic.pgs.profile',header=T)

GRS <- read.table("Wang_multi_ethnic.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Wang"))

write.table(GRS, "Wang_multi_ethnic.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Wang_multi_ethnic.pgs.tsv"))


#################################
# Option #7 - Wang European GRS #
#################################


## Filter to Wang 451-SNP GRS-specific variants, and (optionally) to variants with p value below 0.05

Wang2 <- Wang %>%
  dplyr::filter(!is.na(`Wang et al., 451 SNPs`)) #%>%
#dplyr::filter(Pvalue_EUR <= 0.05)

## Prepare for bgenix GRS calculation

Wang2=arrange(Wang2,Chromosome,`Position (GRCh37)`)
Wang2=Wang2%>%rename(
  rsID=rsID,
  CHR=Chromosome,
  POS=`Position (GRCh37)`,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`OR_EUR`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Wang3=Wang2%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Wang3$CHR[is.na(Wang3$CHR)]="X"
write.table(Wang3,'Wang.tsv',quote=FALSE,sep='\t',row.names = FALSE)


GRS_out=ukbrapR:::create_pgs(
  in_file='Wang.tsv',
  out_file='Wang_European.pgs',
  pgs_name='Wang',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Wang_European.pgs.bim')
#outscore=read.table('Wang_European.pgs.profile',header=T)

GRS <- read.table("Wang_European.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Wang"))

write.table(GRS, "Wang_European.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Wang_European.pgs.tsv"))



################################
# Option #8 - Wang African GRS #
################################


## Filter to Wang 451-SNP GRS-specific variants, and (optionally) to variants with p value below 0.05

Wang2 <- Wang %>%
  dplyr::filter(!is.na(`Wang et al., 451 SNPs`)) #%>%
#dplyr::filter(Pvalue_AFR <= 0.05)

## Prepare for bgenix GRS calculation

Wang2=arrange(Wang2,Chromosome,`Position (GRCh37)`)
Wang2=Wang2%>%rename(
  rsID=rsID,
  CHR=Chromosome,
  POS=`Position (GRCh37)`,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`OR_AFR`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Wang3=Wang2%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Wang3$CHR[is.na(Wang3$CHR)]="X"
write.table(Wang3,'Wang.tsv',quote=FALSE,sep='\t',row.names = FALSE)


GRS_out=ukbrapR:::create_pgs(
  in_file='Wang.tsv',
  out_file='Wang_African.pgs',
  pgs_name='Wang',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Wang_African.pgs.bim')
#outscore=read.table('Wang_African.pgs.profile',header=T)

GRS <- read.table("Wang_African.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Wang"))

write.table(GRS, "Wang_African.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Wang_African.pgs.tsv"))


###################################
# Option #9 - Wang East Asian GRS #
###################################


## Filter to Wang 451-SNP GRS-specific variants, and (optionally) to variants with p value below 0.05

Wang2 <- Wang %>%
  dplyr::filter(!is.na(`Wang et al., 451 SNPs`)) #%>%
#dplyr::filter(Pvalue_EAS <= 0.05)

## Prepare for bgenix GRS calculation

Wang2=arrange(Wang2,Chromosome,`Position (GRCh37)`)
Wang2=Wang2%>%rename(
  rsID=rsID,
  CHR=Chromosome,
  POS=`Position (GRCh37)`,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`OR_EAS`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Wang3=Wang2%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Wang3$CHR[is.na(Wang3$CHR)]="X"
write.table(Wang3,'Wang.tsv',quote=FALSE,sep='\t',row.names = FALSE)


GRS_out=ukbrapR:::create_pgs(
  in_file='Wang.tsv',
  out_file='Wang_East_Asian.pgs',
  pgs_name='Wang',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Wang_East_Asian.pgs.bim')
#outscore=read.table('Wang_East_Asian.pgs.profile',header=T)

GRS <- read.table("Wang_East_Asian.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Wang"))

write.table(GRS, "Wang_East_Asian.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Wang_East_Asian.pgs.tsv"))



##################################
# Option #10 - Wang Hispanic GRS #
##################################


## Filter to Wang 451-SNP GRS-specific variants, and (optionally) to variants with p value below 0.05

Wang2 <- Wang %>%
  dplyr::filter(!is.na(`Wang et al., 451 SNPs`)) #%>%
#dplyr::filter(Pvalue_HIS <= 0.05)

## Prepare for bgenix GRS calculation

Wang2=arrange(Wang2,Chromosome,`Position (GRCh37)`)
Wang2=Wang2%>%rename(
  rsID=rsID,
  CHR=Chromosome,
  POS=`Position (GRCh37)`,
  effect_allele=`Risk Allele`,
  other_allele=`Reference Allele`,
  effect_weight=`OR_HIS`
)%>%mutate(
  effect_weight=log(effect_weight)
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Wang3=Wang2%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Wang3$CHR[is.na(Wang3$CHR)]="X"
write.table(Wang3,'Wang.tsv',quote=FALSE,sep='\t',row.names = FALSE)


GRS_out=ukbrapR:::create_pgs(
  in_file='Wang.tsv',
  out_file='Wang_Hispanic.pgs',
  pgs_name='Wang',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Wang_Hispanic.pgs.bim')
#outscore=read.table('Wang_Hispanic.pgs.profile',header=T)

GRS <- read.table("Wang_Hispanic.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Wang"))

write.table(GRS, "Wang_Hispanic.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Wang_Hispanic.pgs.tsv"))





#------------------------------#
# -- Schumacher (2018) GRSs -- #
#------------------------------#

########################################
# Option #11 - Schumacher (147-SNP) GRS #
########################################

Schumacher <- read.delim("~/Schumacher.txt", comment.char="#")
Schumacher=arrange(Schumacher,chr_name,chr_position)
Schumacher=Schumacher%>%rename(
  rsID=`rsID`,
  CHR=chr_name,
  POS=chr_position,
  effect_allele=`effect_allele`,
  effect_weight=`effect_weight`
)%>%mutate(
  other_allele=NA_character_  # Reference alleles will be looked up from reference genome
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Schumacher2=Schumacher%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Schumacher2$CHR[is.na(Schumacher2$CHR)]="X"
write.table(Schumacher2,'Schumacher.tsv',quote=FALSE,sep='\t',row.names = FALSE)


conti_out=ukbrapR:::create_pgs(
  in_file='Schumacher.tsv',
  out_file='Schumacher.pgs',
  pgs_name='Schumacher',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Schumacher.pgs.bim')
#outscore=read.table('Schumacher.pgs.profile',header=T)

GRS <- read.table("Schumacher.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Schumacher"))

write.table(GRS, "Schumacher.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Schumacher.pgs.tsv"))


#######################################################
# Option #12 - Schumacher (130-SNP) GRS from BARCODE1 #
#######################################################

BARCODE1 <- read.delim("~/BARCODE1.txt", comment.char="#")
BARCODE1=arrange(BARCODE1,chr_name,chr_position)
BARCODE1=BARCODE1%>%rename(
  rsID=`rsID`,
  CHR=chr_name,
  POS=chr_position,
  effect_allele=`effect_allele`,
  other_allele=`other_allele`,
  effect_weight=`effect_weight`
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

BARCODE1_2=BARCODE1%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

BARCODE1_2$CHR[is.na(BARCODE1_2$CHR)]="X"
write.table(BARCODE1_2,'BARCODE1.tsv',quote=FALSE,sep='\t',row.names = FALSE)


conti_out=ukbrapR:::create_pgs(
  in_file='BARCODE1.tsv',
  out_file='BARCODE1.pgs',
  pgs_name='BARCODE1',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('BARCODE1.pgs.bim')
#outscore=read.table('BARCODE1.pgs.profile',header=T)

GRS <- read.table("BARCODE1.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "BARCODE1"))

write.table(GRS, "BARCODE1.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "BARCODE1.pgs.tsv"))



#------------------------------------#
# -- Seibert (2018 and 2022) GRSs -- #
#------------------------------------#

#####################################
# Option #13 - Seibert (54-SNP) GRS #
#####################################

Seibert <- read.delim("~/Seibert.txt", comment.char="#")
Seibert=arrange(Seibert,chr_name,chr_position)
Seibert=Seibert%>%rename(
  rsID=`rsID`,
  CHR=chr_name,
  POS=chr_position,
  effect_allele=`effect_allele`,
  other_allele=`other_allele`,
  effect_weight=`effect_weight`
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Seibert2=Seibert%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Seibert2$CHR[is.na(Seibert2$CHR)]="X"
write.table(Seibert2,'Seibert.tsv',quote=FALSE,sep='\t',row.names = FALSE)


conti_out=ukbrapR:::create_pgs(
  in_file='Seibert.tsv',
  out_file='Seibert.pgs',
  pgs_name='Seibert',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Seibert.pgs.bim')
#outscore=read.table('Seibert.pgs.profile',header=T)

GRS <- read.table("Seibert.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Seibert"))

write.table(GRS, "Seibert.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Seibert.pgs.tsv"))


#######################################
# Option #14 - Pagadala (290-SNP) GRS #
#######################################

Pagadala <- read.delim("~/Pagadala.txt", comment.char="#")
Pagadala=arrange(Pagadala,chr_name,chr_position)
Pagadala=Pagadala%>%rename(
  rsID=`rsID`,
  CHR=chr_name,
  POS=chr_position,
  effect_allele=`effect_allele`,
  other_allele=`other_allele`,
  effect_weight=`effect_weight`
)%>%select(
  rsID,CHR,POS,effect_allele,other_allele,effect_weight
)

Pagadala2=Pagadala%>%mutate(CHR=as.numeric(CHR))%>%arrange(CHR,POS)

Pagadala2$CHR[is.na(Pagadala2$CHR)]="X"
write.table(Pagadala2,'Pagadala.tsv',quote=FALSE,sep='\t',row.names = FALSE)


conti_out=ukbrapR:::create_pgs(
  in_file='Pagadala.tsv',
  out_file='Pagadala.pgs',
  pgs_name='Pagadala',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Pagadala.pgs.bim')
#outscore=read.table('Pagadala.pgs.profile',header=T)

GRS <- read.table("Pagadala.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Pagadala"))

write.table(GRS, "Pagadala.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Pagadala.pgs.tsv"))
