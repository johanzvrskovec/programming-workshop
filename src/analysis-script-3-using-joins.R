#The analysis converted into using joins
#install.packages('tictoc')
library(tictoc)

set.seed(123) #this is for the randomised generation of a data subset

tic()
dfVariantReference<-read.table(gzfile("/Users/jakz/Documents/work_rstudio/programming-workshop/data/1kgp3.bX.eur.l2.jz2023.gz"),
              sep = "\t",
              header = TRUE,
              fill = TRUE,
              blank.lines.skip = TRUE,
              comment.char = ""
) #264.74 sec
toc()

nrow(dfVariantReference)
colnames(dfVariantReference)


tic()
dfGWAS<-read.table(gzfile("/Users/jakz/Documents/work_rstudio/programming-workshop/data/GCST90027164_buildGRCh37.tsv.gz"),
               sep = "\t",
               header = TRUE,
               fill = TRUE,
               blank.lines.skip = TRUE,
               comment.char = ""
) #98.514 sec
toc()

nrow(dfGWAS)
colnames(dfGWAS)

#let's pretend we are interested in this subset
dfGWAS$rand<-rnorm(n = nrow(dfGWAS), mean = 0, sd = 1)
#lSelectedVariants<-dfGWAS[dfGWAS$rand>0.5,]$rsid
#length(lSelectedVariants)
dfGWAS$selected<-dfGWAS$rand>0.5

#filter dfVariantReference
dfVariantReference.filtered<-dfVariantReference[dfVariantReference$MAF<0.25,]
nrow(dfVariantReference.filtered)

#filter dfGWAS
dfGWAS.filtered<-dfGWAS[dfGWAS$selected,]
nrow(dfGWAS.filtered)

#prepare join variables, because standard R can't merge on two variables
dfVariantReference.filtered$mergeVar<-paste0(dfVariantReference.filtered$SNP,"_",dfVariantReference.filtered$CHR,"_",dfVariantReference.filtered$BP)
dfGWAS.filtered$mergeVar<-paste0(dfGWAS.filtered$rsid,"_",dfGWAS.filtered$chromosome,"_",dfGWAS.filtered$base_pair_location)

#merge
dfOutput<-merge(dfVariantReference.filtered,dfGWAS.filtered,by = "mergeVar", all = FALSE)

#the analysis
dfOutput$Z<-dfOutput$beta/dfOutput$standard_error

#trim output
dfOutput.final<-dfOutput[,c("SNP","Z")]


