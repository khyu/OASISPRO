# writeRNAFeatures.R
#
# output RNA-seq feature weights to files
#
# Kun-Hsing Yu
# 2014.11.27

rm(list=ls())
library(FSelector)
setwd("~/Dropbox/CS273AFinalProject/analysis/StageOverall/")

trainFraction<-0.7

# read files
rnaSeqFile<-read.table("../../data/RNASeq/RNAseq.csv", sep=",")
rnaSeqID<-read.table("../../data/RNASeq/RNAseqTcgaID.csv", sep=",")
clinicalFile<-read.csv("../../data/ClinicalInfo/ClinicalAnnot.tsv", sep="\t", header=T)
clinicalID<-read.table("../../data/ClinicalInfo/ClinicalTcgaID.csv")

IDintersection<-intersect(rnaSeqID[,1], clinicalID[,1])
rnaSeq<-rnaSeqFile[(rnaSeqID[,1] %in% IDintersection),]
clinical<-clinicalFile[(clinicalID[,1] %in% IDintersection),]

YFile<-clinical[,40]

Y<-YFile[YFile>(-1)]
Y[Y<20]<-(-1)
Y[Y>=20]<-1
Y<-as.factor(Y)
X<-rnaSeq[YFile>(-1),]

# training, test set: random 0.7, 0.3
set.seed(30)
trainingSet<-sample(length(Y), floor(length(Y)*trainFraction), replace=F)
testSet<-which(!(1:length(Y) %in% trainingSet))

## feature selection
YTrain<-Y[trainingSet]
XTune<-cbind(X[trainingSet,],YTrain)
#numFeatures=20

if (FALSE){
# information gain
weights <- information.gain(YTrain~., XTune)
write.table(weights,"RNAinfog.csv",quote=F,sep=",",row.names=F,col.names=F)

# gain ratio  
weights <- gain.ratio(YTrain~., XTune)
write.table(weights,"RNAgainr.csv",quote=F,sep=",",row.names=F,col.names=F)

# symmetrical uncertainty
weights <- symmetrical.uncertainty(YTrain~., XTune)
write.table(weights,"RNAsymu.csv",quote=F,sep=",",row.names=F,col.names=F)
}

# consistency
weights <- consistency(YTrain~., XTune)
write.table(weights,"RNAcons.csv",quote=F,sep=",",row.names=F,col.names=F)

if (FALSE){
# chi-squared
weights <- chi.squared(YTrain~., XTune)
write.table(weights,"RNAchi.csv",quote=F,sep=",",row.names=F,col.names=F)
}

if (FALSE){
  # cfs
  figureFileName<-paste("rocCFSModels.png",numFeatures,"Features.png",sep="")
  weights <- cfs(YTrain~., XTune)
  subsetFeatures<-cutoff.k(weights, numFeatures)
  X<-X[,weights]
  # random forest importance
  figureFileName<-paste("rocRFIModels.png",numFeatures,"Features.png",sep="")
  weights <- random.forest.importance(YTrain~., XTune, importance.type = 1)
  subsetFeatures<-cutoff.k(weights, numFeatures)
  X<-X[,subsetFeatures]
}
