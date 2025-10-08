#A slightly less messy but still slow analysis script - you probably don't want to run this
#install.packages('tictoc')
library(tictoc)

set.seed(123) #this is for the randomised generation of a data subset

#shared variables
folderpath.project<-file.path("/Users","jakz","Documents","work_rstudio","programming-workshop") #platform independent file/folder path

tic()
dfVariantReference<-read.table(gzfile(file.path(folderpath.project,"data","reference.1000G.maf.0.005.txt.gz")),
              #sep = "\t",
              header = TRUE,
              fill = TRUE,
              blank.lines.skip = TRUE,
              comment.char = ""
) #66.785 sec
toc()

nrow(dfVariantReference)
colnames(dfVariantReference)


tic()
dfGWAS<-read.table(gzfile(file.path(folderpath.project,"data","GCST90027164_buildGRCh37.tsv.gz")),
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
  if(length(CHR1)>0 & length(CHR2)>0 & length(BP1)>0 & length(BP2)>0 & length(MAF)>0 & CHR1==CHR2 & BP1==BP2 & MAF<0.25) {
    
    BETA<-dfGWAS[dfGWAS$rsid==cVariantRsId,c("beta")]
    SE<-dfGWAS[dfGWAS$rsid==cVariantRsId,c("standard_error")]
    Z=BETA/SE

    dfOutput[cVariantRsId,c("Z")]<-Z
  }
  
  if(iVariant %% 10 == 0) cat(".")
  
}

# If you run this loop, you can abort it by selecting Session/Interrupt R

