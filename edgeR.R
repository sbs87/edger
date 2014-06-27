rm(list=ls())
source("http://bioconductor.org/biocLite.R")
source("http://stevenbsmith.net/source/load_R_enviornment_vars.R")
#biocLite("edgeR")
library("edgeR")
#edgeRUsersGuide()
setwd(paste("/Users/",user,"/bin/Qualifiers/edgeR",sep=""))
datafile =system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
pasillaCountTable =read.table( datafile, header=TRUE, row.names="gene_id" )
head(pasillaCountTable)

pasillaDesign = data.frame(
  row.names = colnames( pasillaCountTable ),
  condition = c( "untreated", "untreated", "untreated",
                 "untreated", "treated", "treated", "treated" ),
  libType = c( "single-end", "single-end", "paired-end",
               "paired-end", "single-end", "paired-end", "paired-end" ) )
pasillaDesign

pairedSamples = pasillaDesign$libType == "paired-end"
countTable = pasillaCountTable[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]

head(countTable)

group<-factor(c(1,1,2,2))



y <- DGEList(counts=countTable,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
topTags(et)


