### YOU WILL NEED TO RUN THESE TWO LINES SEPARETELY (AND INDIVIDUALLY), THEN THE REST OF THE SCRIPT
remotes::install_github("lcpilling/ukbrapR@v0.3.10")
install.packages('readxl')



library('dplyr')
library(readxl)
system('dx download file-J4BBx88Jj59xJPB35Kv3f114') # This is the supplementary table file from Conti et al. (2021), downloadable at: https://www.nature.com/articles/s41588-020-00748-0






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

system(paste("dx upload", "multi_ethnic.pgs.tsv"))






