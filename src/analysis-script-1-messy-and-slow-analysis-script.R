#The 'messy and slow' analysis script.
#The script is actually so rudimentary that it has not become slow yet (apart from reading the data files).
#install.packages('tictoc')
library(tictoc)

tic()
d<-read.table(gzfile("/Users/jakz/Documents/work_rstudio/programming-workshop/data/reference.1000G.maf.0.005.txt.gz"),
              #sep = "\t",
              header = TRUE,
              fill = TRUE,
              blank.lines.skip = TRUE,
              comment.char = ""
) #66.785 sec
toc()

nrow(d)
colnames(d)


tic()
d2<-read.table(gzfile("/Users/jakz/Documents/work_rstudio/programming-workshop/data/GCST90027164_buildGRCh37.tsv.gz"),
               sep = "\t",
               header = TRUE,
               fill = TRUE,
               blank.lines.skip = TRUE,
               comment.char = ""
) #98.514 sec
toc()

nrow(d2)
colnames(d2)



saved<-c()
saved.z<-c()


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
# Pretend there is an additional number of manually entered analyses for other rs-ids here as well.
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



