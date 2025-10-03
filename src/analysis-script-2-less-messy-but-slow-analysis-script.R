#The 'messy and slow' analysis script - you probably don't want to run this
#install.packages('tictoc')
library(tictoc)

#First run

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
) #253.414 sec
toc()

nrow(dfGWAS)
colnames(dfGWAS)


lSelectedVariants <- c('rs11240777','rs113656530','rs2096536')

saved<-c()
saved.z<-c()

for(iVariant in 1:length())


#the 'analysis'
#rs11240777
CHR1<-d[d$SNP=='rs11240777',c("CHR")]
CHR2<-d2[d2$rsid=='rs11240777',c("chromosome")]
MAF<-d[d$SNP=='rs11240777',c("MAF")]
BP1<-d[d$SNP=='rs11240777',c("BP")]
BP2<-d2[d2$rsid=='rs11240777',c("base_pair_location")]
if(CHR1==CHR2 & BP1==BP2 & MAF<0.25) {
  saved[length(saved)+1]<-'rs11240777' 
  BETA<-d2[d2$rsid=='rs11240777',c("beta")]
  SE<-d2[d2$rsid=='rs11240777',c("standard_error")]
  Z=BETA/SE
  saved.z[length(saved.z)+1]<-Z
}


#rs113656530
CHR1<-d[d$SNP=='rs113656530',c("CHR")]
CHR2<-d2[d2$rsid=='rs113656530',c("chromosome")]
MAF<-d[d$SNP=='rs113656530',c("MAF")]
BP1<-d[d$SNP=='rs113656530',c("BP")]
BP2<-d2[d2$rsid=='rs113656530',c("base_pair_location")]
if(CHR1==CHR2 & BP1==BP2 & MAF<0.25) {
  saved[length(saved)+1]<-'rs113656530' 
  BETA<-d2[d2$rsid=='rs113656530',c("beta")]
  SE<-d2[d2$rsid=='rs113656530',c("standard_error")]
  Z=BETA/SE
  saved.z[length(saved.z)+1]<-Z
}


#.
#.
#.

#rs2096536
CHR1<-d[d$SNP=='rs2096536',c("CHR")]
CHR2<-d2[d2$rsid=='rs2096536',c("chromosome")]
MAF<-d[d$SNP=='rs2096536',c("MAF")]
BP1<-d[d$SNP=='rs2096536',c("BP")]
BP2<-d2[d2$rsid=='rs2096536',c("base_pair_location")]
if(CHR1==CHR2 & BP1==BP2 & MAF<0.25) {
  saved[length(saved)+1]<-'rs2096536' 
  BETA<-d2[d2$rsid=='rs2096536',c("beta")]
  SE<-d2[d2$rsid=='rs2096536',c("standard_error")]
  Z=BETA/SE
  saved.z[length(saved.z)+1]<-Z
}



