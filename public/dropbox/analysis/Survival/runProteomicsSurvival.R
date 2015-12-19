# runProteomicsSurvival.R
#
# Proteomics
# LASSO Cox Survival
#
# Kun-Hsing Yu
# 2014.11.19

rm(list=ls())
library(survival)
library(glmnet)
library(rms)
library(ggplot2)
library(GGally)

setwd("public/dropbox/analysis/Survival")
sysargs<-commandArgs(trailingOnly=TRUE)
numFeatures<-as.numeric(sysargs[1])

# read files
proteomicsFile<-read.table("../../data/Proteomics/proteomics.csv", sep=",")
proteomicsNameFile<-read.table("../../data/Proteomics/proteomics.txt", sep="\t")
proteomicsID<-proteomicsNameFile[2:dim(proteomicsNameFile)[1],1]
proteinNames<-proteomicsNameFile[1,2:dim(proteomicsNameFile)[2]]
#clinicalFile<-read.csv("../../data/ClinicalInfo/ClinicalAnnot.tsv", sep="\t", header=T)
clinicalID<-read.table("../../data/ClinicalInfo/ClinicalTcgaID.csv")
survivedDays<-read.table("../../data/ClinicalInfo/survivedDays.csv",sep="")
event<-read.table("../../data/ClinicalInfo/event.csv",sep="")

IDintersection<-intersect(proteomicsID, clinicalID[,1])
proteomics<-proteomicsFile[(proteomicsID %in% IDintersection),]
survivedDays<-survivedDays[(clinicalID[,1] %in% IDintersection),]
event<-event[(clinicalID[,1] %in% IDintersection),]
#clinical<-clinicalFile[(clinicalID[,1] %in% IDintersection),]

#X<-cbind(rnaSeq,clinical[,c(7,8,13,15,17,18,19,20,26,27,28,29,30,44,45,52)])
X<-proteomics

X<-X[survivedDays!=-1,]
event<-event[survivedDays!=-1]
survivedDays<-survivedDays[survivedDays!=-1]

Ymatrix<-Surv(survivedDays,event)
Ymatrix[,1]<-Ymatrix[,1]+1

# training, test set 1: 1 to 256, 257 to 366
#trainingSet<-1:256
#testSet<-257:length(survivedDays)
# training, test set 2: random 256, 110
set.seed(30)
trainingSet<-sample(length(survivedDays), 256, replace=F)
testSet<-which(!(1:length(survivedDays) %in% trainingSet))

# cox
#cv.tr<-cv.glmnet(as.matrix(X[trainingSet,]),Ymatrix[trainingSet],family='cox',nfolds=10)#nrow(X[trainingSet,]))
#plot(cv.tr)
#fit.tr<-glmnet(as.matrix(X[trainingSet,]),Ymatrix[trainingSet],family='cox',lambda=cv.tr$lambda.min)
#fit.tr<-glmnet(as.matrix(X[trainingSet,]),Ymatrix[trainingSet],family='cox',lambda=cv.tr$lambda.1se)
#fit.tr<-glmnet(as.matrix(X[trainingSet,]),Ymatrix[trainingSet],family='cox',lambda=seq(0.0001,0.07,0.0001))
fit.tr<-glmnet(as.matrix(X[trainingSet,]),Ymatrix[trainingSet],family='cox',lambda=seq(0.0001,0.07,0.0001))

# Make predictions
predAllM<-predict(fit.tr,as.matrix(X),type='response')
predTrainM<-predict(fit.tr,as.matrix(X[trainingSet,]),type='response')
predTestM<-predict(fit.tr,as.matrix(X[testSet,]),type='response')

matrixIndex<-floor(median(which(fit.tr$df==numFeatures)))

predAll<-predAllM[,matrixIndex]
predTrain<-predTrainM[,matrixIndex]
predTest<-predTestM[,matrixIndex]


thres<-median(predTrain)

plotAll<-data.frame(cbind(predAll,survivedDays,event))
for (i in 1:dim(plotAll)[1]){
  if (plotAll[i,1]>=thres){
    plotAll[i,1]<-"Poor"
  } else {
    plotAll[i,1]<-"Good"
  }  
}

#plotAll.surv <- survfit(Surv(survivedDays, event) ~ s0, data = plotAll) 
plotAll.surv <- survfit(Surv(survivedDays, event) ~ predAll, data = plotAll) 
plot(plotAll.surv, lty = 2:3) 


plotTrain<-data.frame(cbind(predTrain,survivedDays[trainingSet],event[trainingSet]))
for (i in 1:dim(plotTrain)[1]){
  if (plotTrain[i,1]>=thres){
    plotTrain[i,1]<-"Poor"
  } else {
    plotTrain[i,1]<-"Good"
  }  
}
#plotTrain.surv <- survfit(Surv(V2, V3) ~ s0, data = plotTrain) 
plotTrain.surv <- survfit(Surv(V2, V3) ~ predTrain, data = plotTrain) 
plot(plotTrain.surv, lty = 2:3) 


errorThres<-0.5#.8
plotTest<-data.frame(cbind(predTest,survivedDays[testSet],event[testSet]))
for (i in 1:dim(plotTest)[1]){
  if (plotTest[i,1]>=(thres+0)){
    plotTest[i,1]<-"Poor"
  } else if (plotTest[i,1]<(thres-0)){
    plotTest[i,1]<-"Good"
  } else {
    plotTest[i,1]<-"Unknown"
  }
}
plotTest<-plotTest[plotTest[,1]!="Unknown",]
#plotTest.surv <- survfit(Surv(V2, V3) ~ s0, data = plotTest)
plotTest.surv <- survfit(Surv(V2, V3) ~ predTest, data = plotTest)
plot(plotTest.surv, lty = 2:3) 

# plot by ggsurv
plotTestGGsurv<-plotTest
plotTestGGsurv[,2]<-plotTestGGsurv[,2]/12
plotTest.surv <- survfit(Surv(V2, V3) ~ predTest, data = plotTestGGsurv)

# re-define ggsurv
source("ggsurv.R")

png(filename="survivalOutput.png", width=1000, height=500)
ggsurv(plotTest.surv) + 
  guides(linetype = F) + 
  xlab("Months") + ggtitle("Kaplan-Meier Curve") +
  theme(plot.title = element_text(size=24,face="bold"),
        legend.text=element_text(size=24), 
        legend.title=element_text(size=24), 
        axis.text=element_text(size=18), 
        axis.title=element_text(size=24,face="bold")) +
  scale_colour_discrete(name = 'Survival Groups', labels=c('Good Prognostic Group', 'Poor Prognostic Group'))

# Close and save the PNG file.
dev.off()
