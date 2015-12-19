# writeProteomicsFeatures.R
#
# output proteomics feature weights to files
#
# Kun-Hsing Yu
# 2014.11.28

rm(list=ls())
library(FSelector)
setwd("~/Dropbox/CS273AFinalProject/analysis/StageOverall/")

trainFraction<-0.7

# read files
proteomicsFile<-read.table("../../data/Proteomics/proteomics.csv", sep=",")
proteomicsNameFile<-read.table("../../data/Proteomics/proteomics.txt", sep="\t")
proteomicsID<-proteomicsNameFile[2:dim(proteomicsNameFile)[1],1]
proteinNames<-proteomicsNameFile[1,2:dim(proteomicsNameFile)[2]]

clinicalFile<-read.csv("../../data/ClinicalInfo/ClinicalAnnot.tsv", sep="\t", header=T)
clinicalID<-read.table("../../data/ClinicalInfo/ClinicalTcgaID.csv")

IDintersection<-intersect(proteomicsID, clinicalID[,1])
proteomics<-proteomicsFile[(proteomicsID %in% IDintersection),]
clinical<-clinicalFile[(clinicalID[,1] %in% IDintersection),]

YFile<-clinical[,40]

Y<-YFile[YFile>(-1)]
Y[Y<20]<-(-1)
Y[Y>=20]<-1
Y<-as.factor(Y)
X<-proteomics[YFile>(-1),]

# training, test set: random 0.7, 0.3
set.seed(30)
trainingSet<-sample(length(Y), floor(length(Y)*trainFraction), replace=F)
testSet<-which(!(1:length(Y) %in% trainingSet))

## feature selection
YTrain<-Y[trainingSet]
XTune<-cbind(X[trainingSet,],YTrain)

if (FALSE){
# information gain
weights <- information.gain(YTrain~., XTune)
write.table(weights,"Proteomicsinfog.csv",quote=F,sep=",",row.names=F,col.names=F)

# gain ratio  
weights <- gain.ratio(YTrain~., XTune)
write.table(weights,"Proteomicsgainr.csv",quote=F,sep=",",row.names=F,col.names=F)

# symmetrical uncertainty
weights <- symmetrical.uncertainty(YTrain~., XTune)
write.table(weights,"Proteomicssymu.csv",quote=F,sep=",",row.names=F,col.names=F)
}

if (FALSE){
# consistency
weights <- consistency(YTrain~., XTune)
write.table(weights,"Proteomicscons.csv",quote=F,sep=",",row.names=F,col.names=F)
}

# chi-squared
weights <- chi.squared(YTrain~., XTune)
write.table(weights,"Proteomicschi.csv",quote=F,sep=",",row.names=F,col.names=F)

if (FALSE){
  # cfs
  weights <- cfs(YTrain~., XTune)
  subsetFeatures<-cutoff.k(weights, numFeatures)
  X<-X[,weights]
  # random forest importance
  weights <- random.forest.importance(YTrain~., XTune, importance.type = 1)
  subsetFeatures<-cutoff.k(weights, numFeatures)
  X<-X[,subsetFeatures]
}
