#A less 'messy' but still slow analysis script, applying control structures
#install.packages('tictoc')
library(tictoc)



#First run

tic()
d<-read.table(gzfile("/Users/jakz/Documents/work_rstudio/programming-workshop/data/1kgp3.bX.eur.l2.jz2023.gz"),
              sep = "\t",
              header = TRUE,
              fill = TRUE,
              blank.lines.skip = TRUE,
              comment.char = ""
) #264.74 sec
toc()

nrow(d)
colnames(d)


tic()
d2<-read.table(gzfile("/Users/jakz/Documents/work_rstudio/programming-workshop/data/reference.1000G.maf.0.005.txt.gz"),
              #sep = "\t",
              header = TRUE,
              fill = TRUE,
              blank.lines.skip = TRUE,
              comment.char = ""
) #64.366 sec
toc()

nrow(d2)
colnames(d2)


tic()
d3<-read.table(gzfile("/Users/jakz/Documents/work_rstudio/programming-workshop/data/sumstats_neuroticism_ctg_format.txt.gz"),
               sep = "\t",
               header = TRUE,
               fill = TRUE,
               blank.lines.skip = TRUE,
               comment.char = ""
) #253.414 sec
toc()

nrow(d3)
colnames(d3)


##merge on id/SNP
tic()
d.merged <- merge(d, d2, by = "SNP", all.x = TRUE, all.y = FALSE)
nrow(d.merged)
colnames(d.merged)
toc()

##subset based on condition
#d.merged.overlapping <- d.merged[!is.na(d.merged$CHR.y),]
#nrow(d.merged.overlapping)

tic()
d.merged.matching.coordinates <- d.merged[!is.na(d.merged$CHR.x) & !is.na(d.merged$CHR.y) & !is.na(d.merged$BP.x) & !is.na(d.merged$BP.y) & d.merged$CHR.x == d.merged$CHR.y & d.merged$BP.x == d.merged$BP.y, ]
nrow(d.merged.matching.coordinates)

d.merged.matching.coordinates2<-d.merged.matching.coordinates

d.merged.matching.coordinates2$CHR<-d.merged.matching.coordinates2$CHR.x
d.merged.matching.coordinates2$BP<-d.merged.matching.coordinates2$BP.x
d.merged.matching.coordinates2$A1<-d.merged.matching.coordinates2$A1.x
d.merged.matching.coordinates2$A2<-d.merged.matching.coordinates2$A2.x
d.merged.matching.coordinates2$MAF<-d.merged.matching.coordinates2$MAF.x
d.merged.matching.coordinates2<-d.merged.matching.coordinates2[,c("SNP","CHR","BP","CM","A1","A2","MAF")]
toc()

##loop operation with multiple datasets
cat("Loop\n")
for(iD in 1:nrow(d.merged.matching.coordinates2)){
  if(iD %% 100000 ==0) cat(".")
}






