#The analysis using fast(?) loops
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
dtPreOutput<-dtVariantReference.filtered[dtGWAS.filtered, on=list(SNP=rsid,CHR=chromosome,BP=base_pair_location), nomatch = NULL]

#What if there were parts of the analysis that we needed/wanted to keep in a loop?

#using data.table row-wise operation
tic()
for(iVariant in 1L:nrow(dtPreOutput)){
  set(x = dtPreOutput, i = iRow, j = "Z", value = dtPreOutput[iVariant,.(beta)]/dtPreOutput[iVariant,.(standard_error)])
  if(iVariant %% 1000 == 0) cat(".")
}
toc() 


#using a Reference Class

myAnalysisClassHandle <- setRefClass(
  "myAnalysisClass",
  fields = list(
    myDF = "ANY"
  ),
  methods = list
  (
    initialize=function(myDFSet)
    {
      myDF <<- as.data.frame(myDFSet)
    }
  )
)

myAnalysisClassHandle$methods(
  rowWiseOperation=function(index){
    myDF[index,c("Z")] <<- myDF[index,c("beta")]/myDF[index,c("standard_error")]
  }
)

myAnalysisClassHandle$methods(
  loopInFunction=function(){
    for(iVariant in 1L:nrow(analysisObj$myDF)){
      myDF[iVariant,c("Z")] <<- myDF[iVariant,c("beta")]/myDF[iVariant,c("standard_error")]
      if(iVariant %% 1000 == 0) cat(".")
    }
  }
)


tic()
analysisObj <- myAnalysisClassHandle(dtPreOutput)
for(iVariant in 1L:nrow(analysisObj$myDF)){
  analysisObj$rowWiseOperation(iVariant)
  if(iVariant %% 1000 == 0) cat(".")
}
toc() 


tic()
analysisObj <- myAnalysisClassHandle(dtPreOutput)
analysisObj$loopInFunction()
toc() 