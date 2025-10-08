#The analysis converted into using data.table
#install.packages('tictoc')
#install.packages('data.table')
library(tictoc)
library(data.table)

set.seed(123) #this is for the randomised generation of a data subset

#shared variables
folderpath.project<-file.path("/Users","jakz","Documents","work_rstudio","programming-workshop") #platform independent file/folder path

tic()
dtVariantReference<-fread(file.path(folderpath.project,"data","reference.1000G.maf.0.005.txt.gz"),
              #sep = "\t",
              header = TRUE,
              fill = TRUE,
              blank.lines.skip = TRUE,
              nThread = 4,
              showProgress = FALSE
) #14.38 sec
toc()

nrow(dtVariantReference)
colnames(dtVariantReference)



tic()
dtGWAS<-fread(file.path(folderpath.project,"data","GCST90027164_buildGRCh37.tsv.gz"),
              sep = "\t",
              header = TRUE,
              fill = TRUE,
              blank.lines.skip = TRUE,
              nThread = 4,
              showProgress = FALSE
) #19.585 sec
toc()

nrow(dtGWAS)
colnames(dtGWAS)

tic()
#let's pretend we are interested in this subset
dtGWAS$rand<-rnorm(n = nrow(dtGWAS), mean = 0, sd = 1) #we use standard data.frame assignment here to keep the behaviour of rnorm from the previous scripts.
dtGWAS[,selected:=rand>0.5] #this is data.table assignment.

#filter dtVariantReference
dtVariantReference.filtered<-dtVariantReference[MAF<0.25,]
nrow(dtVariantReference.filtered)

#filter dtGWAS
dtGWAS.filtered<-dtGWAS[selected==TRUE,]
nrow(dtGWAS.filtered)

setkeyv(dtVariantReference.filtered, cols = c("SNP","CHR","BP")) #we index the columns on which we will access data and join
setkeyv(dtGWAS.filtered, cols = c("rsid","chromosome","base_pair_location")) #we index the columns on which we will access data and join

#merge with an inner join, chaining the analysis, and trimming the output
dtOutput<-dtVariantReference.filtered[dtGWAS.filtered, on=list(SNP=rsid,CHR=chromosome,BP=base_pair_location), nomatch = NULL][,Z:=beta/standard_error][,c("SNP","Z")]

toc() #22.703 sec
