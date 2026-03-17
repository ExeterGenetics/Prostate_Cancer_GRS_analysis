### YOU WILL NEED TO RUN THESE TWO LINES SEPARETELY (AND INDIVIDUALLY), THEN THE REST OF THE SCRIPT
remotes::install_github("lcpilling/ukbrapR@v0.3.10")
install.packages('readxl')

source('https://raw.githubusercontent.com/ExeterGenetics/ukbextractR/main/session_setup.R')


library(dplyr)
library(readxl)
system('dx download file-J4BBx88Jj59xJPB35Kv3f114') # This is the supplementary table file from Conti et al. (2021), downloadable at: https://www.nature.com/articles/s41588-020-00748-0

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

###############################
# Option #1 - Multiethnic GRS #
###############################


# The below block of code reads the conti xlsx file, and puts it into a form that ukbrapR can use to make a grs (the outputted tsv)

#If you change the effect_weight column it should be easy to run a different GRS

Conti <- read_excel("Conti.xlsx", sheet = "S4", skip = 3, na = "NA")
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
  out_file='multi_ethnic.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('multi_ethnic.pgs.bim')
outscore=read.table('multi_ethnic.pgs.profile',header=T)

GRS <- read.table("multi_ethnic.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "multi_ethnic.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "multi_ethnic.pgs.tsv"))



#####################################
# Option #2 - European-specific GRS #
#####################################

Conti <- read_excel("Conti.xlsx", sheet = "S4", skip = 3, na = "NA")
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
  out_file='European.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('European.pgs.bim')
#outscore=read.table('European.pgs.profile',header=T)

GRS <- read.table("European.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "European.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "European.pgs.tsv"))



#####################################
# Option #3 - African-specific GRS #
#####################################

Conti <- read_excel("Conti.xlsx", sheet = "S4", skip = 3, na = "NA")
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
  out_file='African.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('African.pgs.bim')
#outscore=read.table('African.pgs.profile',header=T)

GRS <- read.table("African.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "African.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "African.pgs.tsv"))




#######################################
# Option #4 - East Asian-specific GRS #
#######################################

Conti <- read_excel("Conti.xlsx", sheet = "S4", skip = 3, na = "NA")
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
  out_file='East_Asian.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('East_Asian.pgs.bim')
#outscore=read.table('East_Asian.pgs.profile',header=T)

GRS <- read.table("East_Asian.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "East_Asian.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "East_Asian.pgs.tsv"))


#####################################
# Option #5 - Hispanic-specific GRS #
#####################################

Conti <- read_excel("Conti.xlsx", sheet = "S4", skip = 3, na = "NA")
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
  out_file='Hispanic.pgs',
  pgs_name='Conti',
  use_imp_pos=TRUE,
  very_verbose=TRUE, # can probably remove
  overwrite=TRUE # overwrites files with same name
)

outbim=read.table('Hispanic.pgs.bim')
#outscore=read.table('Hispanic.pgs.profile',header=T)

GRS <- read.table("Hispanic.pgs.tsv", header = TRUE)

GRS <- exclude_withdrawn(GRS) %>%
  dplyr::select(c("eid", "Conti"))

write.table(GRS, "Hispanic.pgs.tsv", quote=FALSE, sep='\t',row.names = FALSE)

## Upload to project

system(paste("dx upload", "Hispanic.pgs.tsv"))