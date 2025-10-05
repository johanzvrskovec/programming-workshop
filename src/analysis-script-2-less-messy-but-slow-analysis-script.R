#A slightly less messy but still slow analysis script - you probably don't want to run this
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
lSelectedVariants<-dfGWAS[dfGWAS$rand>0.5,]$rsid
length(lSelectedVariants)

#improved output/metadata
dfOutput <- as.data.frame(matrix(NA,ncol=0,nrow=0))

#the analyses are now contained in a for-loop
for(iVariant in 1:length(lSelectedVariants)) {
  #the 'analysis'
  #iVariant<-1
  
  cVariantRsId<-lSelectedVariants[iVariant]
  CHR1 <- dfVariantReference[dfVariantReference$SNP==cVariantRsId,c("CHR")]
  CHR2 <- dfGWAS[dfGWAS$rsid==cVariantRsId,c("chromosome")]
  MAF <- dfVariantReference[dfVariantReference$SNP==cVariantRsId,c("MAF")]
  BP1 <- dfVariantReference[dfVariantReference$SNP==cVariantRsId,c("BP")]
  BP2 <- dfGWAS[dfGWAS$rsid==cVariantRsId,c("base_pair_location")]
  if(CHR1==CHR2 & BP1==BP2 & MAF<0.25) {
    
    BETA<-dfGWAS[dfGWAS$rsid==cVariantRsId,c("beta")]
    SE<-dfGWAS[dfGWAS$rsid==cVariantRsId,c("standard_error")]
    Z=BETA/SE

    dfOutput[cVariantRsId,c("Z")]<-Z
  }
  
  if(iVariant %% 100 ==0) cat(".")
  
}

