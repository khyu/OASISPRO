# elasticNetCox.R
#
# train and test binary classifiers
# plot ROC for classifiers
#
# param 1: tumor type
# param 2: data type (RNAseq, proteomics, etc)
# param 3: partition type (random, batch)
# param 4: partition seed (random seed, name of batch file)
# param 5: percentage of data in the training set
# param 6: lower bound of alpha 
# param 7: upper bound of alpha
# param 8: lower bound of lambda
# param 9: upper bound of lambda
# param 10: hand-picked clinical variables (file)
#
# Kun-Hsing Yu
# Dec. 8, 2015

rm(list=ls())
library(FSelector)
library(ggplot2)
library(survival)
library(glmnet)
library(rms)
library(GGally)
library(dplyr)
#setwd("public/dropbox/analysis/StageOverall/")

# sysargs<-commandArgs(trailingOnly=TRUE)
# print (sysargs)
# tumorType<-sysargs[1]
# dataType<-sysargs[2]
# partitionType<-sysargs[3]
# if (partitionType == "random") {
# 	partitionSeed<-sysargs[4]
# 	trainingPercentage<-sysargs[5]
# } else {
# 	partitionSeed<-(-1)
# 	selectedBatchesFile<-read.table(sysargs[4], sep="")
# 	selectedBatches<-selectedBatchesFile[,1]
# 	trainingPercentage<-(-1)
# } 
# minAlpha<-as.numeric(sysargs[6])
# maxAlpha<-as.numeric(sysargs[7])
# minLambda<-as.numeric(sysargs[8])
# maxLambda<-as.numeric(sysargs[9])
# if (sysargs[10] != "-1") {
#   clinicalVariablesFile<-read.table(sysargs[10], sep="")
#   clinicalVariables<-clinicalVariablesFile[,1]
# }

tumorTypesFile<-read.table("../../data/tumorTypes.txt", stringsAsFactors = F)
tumorTypes<-tumorTypesFile[,1]

#for (i in 4:length(tumorTypes)) {
  
  tumorType<-tumorTypes[5]
  dataType<-"proteomics"
  partitionType<-"random"
  if (partitionType == "random") {
    partitionSeed<-1
    trainingPercentage<-0.7
  } else {
    partitionSeed<-(-1)
    selectedBatchesFile<-read.table(sysargs[4], sep="")
    selectedBatches<-selectedBatchesFile[,1]
    trainingPercentage<-(-1)
  } 
  minAlpha<--1
  maxAlpha<--1
  minLambda<--1
  maxLambda<--1
  
  print (maxAlpha)
  print (tumorType)
  #print (clinicalVariables)
  
  # read files
  omicsFile<-read.table(paste("../../data/", tumorType, "_", dataType, ".txt", sep=""), stringsAsFactors=F, sep=",")
  print("Finished reading omics file")
  omicsNameFile<-read.table(paste("../../data/", tumorType, "_", dataType, "_ids.txt", sep=""), stringsAsFactors=F, sep="")
  omics<-omicsFile[substr(omicsNameFile[,1],14,15)=="01",]
  omicsName<-omicsNameFile[substr(omicsNameFile[,1],14,15)=="01",]
  omicsID<-substr(omicsName, 1, 12)
  
  clinicalFile<-read.table(paste("../../data/nationwidechildrens.org_clinical_patient_", tumorType, ".txt", sep=""), stringsAsFactors = F, quote="", sep="\t")
  clinical<-clinicalFile[4:dim(clinicalFile)[1],]
  colnames(clinical)<-clinicalFile[2,]
  clinicalID<-clinical[,2]
  
  event<-clinical[,"vital_status"]
  survivedDays_alive<-clinical[,"days_to_last_followup"]
  survivedDays_dead<-clinical[,"days_to_death"]
  survivedDays<-ifelse(event == "alive", survivedDays_alive, survivedDays_dead)
  
  IDintersection<-intersect(omicsID, clinicalID)
  omics<-omics[(omicsID %in% IDintersection),]
  clinical<-clinical[(clinicalID %in% IDintersection),]
  intersectIDs<-omicsID[(omicsID %in% IDintersection)]
  
  survivedDays<-survivedDays[(clinicalID %in% IDintersection)]
  event<-event[(clinicalID %in% IDintersection)]
  
  X<-omics
  
  # if (sysargs[10] == "-1") {
  #   X<-omics
  # } else {
  #   #X<-cbind(rnaSeq,clinical[,c(7,8,13,15,17,18,19,20,26,27,28,29,30,44,45,52)])
  #   clinicalVarColNums<-which(colnames(clinical) %in% clinicalVariables)
  #   X<-cbind(omics, clinical[,clinicalVarColNums])
  # }
  
  survivedDays<-ifelse(substr(survivedDays,1,2)=="[N", -1, survivedDays)
  survivedDays<-as.numeric(survivedDays)
  event<-ifelse(event=="Alive", 0, 1)
  X<-X[survivedDays!=-1,]
  event<-event[survivedDays!=-1]
  iDs<-intersectIDs[survivedDays!=-1]
  survivedDays<-survivedDays[survivedDays!=-1]
  
  Ymatrix<-Surv(survivedDays,event)
  Ymatrix[,1]<-Ymatrix[,1]+1
  
  print("Ymatrix done")
  
  if (partitionType == "random") {
    # set.seed(randSeed)
    trainingSet<-sample(length(survivedDays), floor(length(survivedDays)*trainingPercentage), replace=F)
    testSet<-which(!(1:length(survivedDays) %in% trainingSet))
  } else {
    trainingSet<-which(substr(iDs, 6, 7) %in% selectedBatches)
    testSet<-which(!(1:length(survivedDays) %in% trainingSet))
  }
  
  # cox
  #Find optimal alpha
  if (minAlpha == -1) {
    minAlpha<-0
    maxAlpha<-1
  }
  alphas<-seq(minAlpha,maxAlpha,by=0.1)
  print("alphas done")
  if (minLambda == -1) {
    elasticnet<-lapply(alphas, function(alpha){cv.glmnet(as.matrix(X[trainingSet,]), Ymatrix[trainingSet,], alpha=alpha, family="cox")})
  } else {
    elasticnet<-lapply(alphas, function(alpha){cv.glmnet(as.matrix(X[trainingSet,]), Ymatrix[trainingSet,], lambda=seq(minLambda,maxLambda,0.0001), alpha=alpha, family="cox")})
  }
  print("cv done")
  minCVMs<-rep(0,length(alphas))
  for (i in 1:length(alphas)) {
    minCVMs[i]<-min(elasticnet[[i]]$cvm[2:length(elasticnet[[i]]$cvm)])
  }
  optimalAlphaIndices<-which(minCVMs==min(minCVMs))
  optimalAlphaIndex<-optimalAlphaIndices[1]
  cv.tr<-cv.glmnet(as.matrix(X[trainingSet,]),Ymatrix[trainingSet,],family='cox',alpha=alphas[optimalAlphaIndex])
  if (cv.tr$lambda.min==cv.tr$lambda[1]){
    lambdaFit<-cv.tr$lambda[2]
  } else {
    lambdaFit<-cv.tr$lambda.min
  }
  predAll<-predict(cv.tr,as.matrix(X),s=lambdaFit,type='response')
  predTrain<-predict(cv.tr,as.matrix(X[trainingSet,]),s=lambdaFit,type='response')
  predTest<-predict(cv.tr,as.matrix(X[testSet,]),s=lambdaFit,type='response')
  
  print("prediction done")
  
  # #Make predictions. Use cross validation to find optimal lambda if user did not specify lambda range.
  # if (minLambda == -1) {
  #   cv.tr<-cv.glmnet(X[trainingSet,],Ymatrix[trainingSet,],family='cox',alpha=alphas[optimalAlphaIndex],nfolds=10)
  #   predAll<-predict(cv.tr,as.matrix(X),s=cv.tr$lambda.min,type='response')
  #   predTrain<-predict(cv.tr,as.matrix(X[trainingSet,]),s=cv.tr$lambda.min,type='response')
  #   predTest<-predict(cv.tr,as.matrix(X[testSet,]),s=cv.tr$lambda.min,type='response')
  # } else {
  #   fit.tr<-glmnet(as.matrix(X[trainingSet,]),Ymatrix[trainingSet],family='cox',lambda=seq(minLambda,maxLambda,0.0001), alpha=alphas[optimalAlphaIndex])
  #   predAll<-predict(fit.tr,as.matrix(X),type='response')
  #   predTrain<-predict(fit.tr,as.matrix(X[trainingSet,]),type='response')
  #   predTest<-predict(fit.tr,as.matrix(X[testSet,]),type='response')
  # }
  
  # # Make predictions
  # predAllM<-predict(fit.tr,as.matrix(X),type='response')
  # predTrainM<-predict(fit.tr,as.matrix(X[trainingSet,]),type='response')
  # predTestM<-predict(fit.tr,as.matrix(X[testSet,]),type='response')
  # 
  # matrixIndex<-floor(median(which(fit.tr$df==numFeatures)))
  # 
  # predAll<-predAllM[,matrixIndex]
  # predTrain<-predTrainM[,matrixIndex]
  # predTest<-predTestM[,matrixIndex]
  
  
  thres<-median(predTrain)
  
  plotTrain<-data.frame(cbind(predTrain,survivedDays[trainingSet],event[trainingSet]))
  for (i in 1:dim(plotTrain)[1]){
    if (plotTrain[i,1]>=thres){
      plotTrain[i,1]<-"Poor"
    } else {
      plotTrain[i,1]<-"Good"
    }  
  }
  plotTrain.surv <- survfit(Surv(V2, V3) ~ X1, data = plotTrain) 
  plot(plotTrain.surv, lty = 2:3) 
  
  
  errorThres<-0.5
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
  plotTest.surv <- survfit(Surv(V2, V3) ~ X1, data = plotTest)
  plot(plotTest.surv, lty = 2:3) 
  
  # plot by ggsurv
  plotTestGGsurv<-plotTest
  plotTestGGsurv[,2]<-plotTestGGsurv[,2]/12
  plotTest.surv <- survfit(Surv(V2, V3) ~ X1, data = plotTestGGsurv)
  
  print("everything but ggsurv done")
  
  # re-define ggsurv
  source("../public/dropbox/analysis/Survival/ggsurv.R")
  
  png(filename=paste("./survivalOutput_", tumorType, ".png", sep=""), width=1000, height=500)
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
  
  print("all done")
  
#}
  
  
