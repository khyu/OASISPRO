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
# param 11: session ID
# param 12: user's e-mail
#
# Kun-Hsing Yu
# March 28, 2016

rm(list=ls())
suppressMessages(library(FSelector))
suppressMessages(library(ggplot2))
suppressMessages(library(survival))
suppressMessages(library(glmnet))
suppressMessages(library(rms))
suppressMessages(library(GGally))
suppressMessages(library(dplyr))
suppressMessages(library(cvTools))


sysargs<-commandArgs(trailingOnly=TRUE)
print (sysargs)
sessionID<-sysargs[11]
tumorType<-sysargs[1]
dataType<-sysargs[2]
partitionType<-sysargs[3]
if (partitionType == "random") { # random
	partitionSeed<-sysargs[4]
	trainingPercentage<-as.numeric(sysargs[5])
} else if (partitionType == "batch") { # batch
  partitionSeed<-(-1)
  selectedBatchesFile<-read.table(paste("public/sessions/",sessionID,"/",sysargs[4],sep=""), sep="", colClasses="character")
  selectedBatches<-selectedBatchesFile[,1]
  trainingPercentage<-(-1)
} else if (partitionType == "kfold") {
  nFolds<-as.numeric(sysargs[4])
}

minAlpha<-as.numeric(sysargs[6])
maxAlpha<-as.numeric(sysargs[7])
minLambda<-as.numeric(sysargs[8])
maxLambda<-as.numeric(sysargs[9])
if (sysargs[10] != "-1") {
  clinicalVariablesFile<-read.table(paste("public/sessions/",sessionID,"/",sysargs[10],sep=""), sep="")
  clinicalVariables<-clinicalVariablesFile[,1]
}
userEmail<-sysargs[12]
if (is.na(userEmail)){
  userEmail<-""
}
osType<-system("uname -a", intern=T)
if (substr(osType,1,5) == "Linux"){ ## google cloud, get ip address
  ipAddress<-system(paste("curl -H 'Metadata-Flavor: Google' http://169.254.169.254/computeMetadata/v1/instance/network-interfaces/0/access-configs/0/external-ip 2>public/sessions/",sessionID,"/ipOutput.txt",sep=""), intern=T)
} else { ## apple, test on local host
  ipAddress<-"localhost:3000"
}

argsFileName<-paste("public/sessions/",sessionID,"/args.txt",sep="")
write(sysargs,argsFileName)

print (maxAlpha)
print (tumorType)
#print (clinicalVariables)
milestonesFileName<-paste("public/sessions/",sessionID,"/milestones.txt",sep="")
percentageFinish<-0
print(paste("Initializing...",percentageFinish,sep=","))
write(paste("Initializing...",percentageFinish,sep=","),milestonesFileName)


# read files
percentageFinish<-1
print(paste("Reading omics files...",percentageFinish,sep=","))
write(paste("Reading omics files...",percentageFinish,sep=","),milestonesFileName)

omicsFileName<-paste("../data/", tumorType, "_", dataType, ".RData", sep="")
if (file.exists(omicsFileName)){
  #try(omicsFile<-read.table(omicsFileName, stringsAsFactors=F, sep=","), TRUE)
  try(load(omicsFileName), TRUE)
  if (!exists("omicsFile")) {
    stop ("Omics type not available for this tumor. Please reselect.")
  }
} else {
  stop ("Omics type not available for this tumor. Please reselect.")
}

omicsNameFile<-t(read.table(paste("../data/", tumorType, "_", dataType, "_elemid.txt", sep=""), stringsAsFactors=F, sep=","))
percentageFinish<-2
print(paste("Finished reading omics file",percentageFinish,sep=","))
write(paste("Finished reading omics file",percentageFinish,sep=","),milestonesFileName,append=T)

omicsIDFile<-read.table(paste("../data/", tumorType, "_", dataType, "_ids.txt", sep=""), stringsAsFactors=F, sep="")
omics<-omicsFile[substr(omicsIDFile[,1],14,15)=="01",]
omicsName<-omicsIDFile[substr(omicsIDFile[,1],14,15)=="01",]
omicsID<-substr(omicsName, 1, 12)

clinicalFile<-read.table(paste("../data/nationwidechildrens.org_clinical_patient_", tumorType, ".txt", sep=""), stringsAsFactors = F, quote="", sep="\t")
clinical<-clinicalFile[4:dim(clinicalFile)[1],]
colnames(clinical)<-clinicalFile[2,]
clinicalID<-clinical[,2]

event<-clinical[,"vital_status"]
survivedDays_alive<-clinical[,"days_to_last_followup"]
survivedDays_dead<-clinical[,"days_to_death"]
survivedDays<-ifelse(event == "Alive", survivedDays_alive, survivedDays_dead)

IDintersection<-intersect(omicsID, clinicalID)
omics<-omics[(omicsID %in% IDintersection),]
clinical<-clinical[(clinicalID %in% IDintersection),]
intersectIDs<-omicsID[(omicsID %in% IDintersection)]

survivedDays<-survivedDays[(clinicalID %in% IDintersection)]
event<-event[(clinicalID %in% IDintersection)]

#X<-omics

if (sysargs[10] == "-1") {
  X<-omics
} else {
  #X<-cbind(rnaSeq,clinical[,c(7,8,13,15,17,18,19,20,26,27,28,29,30,44,45,52)])
  clinicalVarColNums<-which(colnames(clinical) %in% clinicalVariables)
  X<-cbind(omics, clinical[,clinicalVarColNums])
}

survivedDays<-ifelse(substr(survivedDays,1,1)=="[", -1, survivedDays)
survivedDays<-as.numeric(survivedDays)
event<-ifelse(event=="Alive", 0, 1)
X<-X[survivedDays>=0,]
event<-event[survivedDays>=0]
iDs<-intersectIDs[survivedDays>=0]
survivedDays<-survivedDays[survivedDays>=0]

Ymatrix<-Surv(survivedDays,event)
Ymatrix[,1]<-Ymatrix[,1]+1

percentageFinish<-5
print(paste("Finished building survival matrix",percentageFinish,sep=","))
write(paste("Finished building survival matrix",percentageFinish,sep=","),milestonesFileName,append=T)

if (partitionType == "random"){ #random
  # set.seed(randSeed)
  trainingSet<-sample(length(survivedDays), floor(length(survivedDays)*trainingPercentage), replace=F)
  testSet<-which(!(1:length(survivedDays) %in% trainingSet))
  nFolds<-1
} else if (partitionType == "kfold"){ # k fold
  folds<-cvFolds(length(survivedDays), K=nFolds)
} else if (partitionType == "LOOCV"){ # LOOCV
  nFolds<-length(survivedDays)
  folds<-cvFolds(length(survivedDays), K=length(survivedDays))
} else { #batch
  trainingSet<-which(substr(iDs, 6, 7) %in% selectedBatches)
  testSet<-which(!(1:length(survivedDays) %in% trainingSet))
  nFolds<-1
}

write.table(c(length(trainingSet),length(testSet)),paste("public/sessions/",sessionID,"/nSamples.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F)

# cox
#Find optimal alpha
if (minAlpha == -1) {
  minAlpha<-0
  maxAlpha<-1
}
alphas<-seq(minAlpha,maxAlpha,by=0.1)
print("alphas done")


percentageStep<-90/nFolds
plotTest<-as.data.frame(cbind(0,survivedDays,event))
for (fold in 1:nFolds) {
  if (nFolds>1){ # kfold or loocv
    trainingSet<-folds$subsets[,1][folds$which != fold]
    testSet<-which(!(1:length(survivedDays) %in% trainingSet))
    print(paste("Fold: ",fold,sep=""))
  }

  if (minLambda == -1) {
    elasticnet<-lapply(alphas, function(alpha){cv.glmnet(as.matrix(X[trainingSet,]), Ymatrix[trainingSet,], alpha=alpha, family="cox")})
  } else {
    elasticnet<-lapply(alphas, function(alpha){cv.glmnet(as.matrix(X[trainingSet,]), Ymatrix[trainingSet,], lambda=seq(minLambda,maxLambda,0.0001), alpha=alpha, family="cox")})
  }

  percentageFinish<-percentageFinish+(percentageStep/2)
  if (nFolds>1){
    print(paste(paste("Fold ",fold,": Finished regularization",sep=""),round(percentageFinish,0),sep=","))
    write(paste(paste("Fold ",fold,": Finished regularization",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
  } else {
    print(paste("Finished regularization",round(percentageFinish,0),sep=","))
    write(paste("Finished regularization",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
  }

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

  percentageFinish<-percentageFinish+(percentageStep/2)
  if (nFolds>1){
    print(paste(paste("Fold ",fold,": Finished survival prediction",sep=""),round(percentageFinish,0),sep=","))
    write(paste(paste("Fold ",fold,": Finished survival prediction",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
  } else {
    print(paste("Finished survival prediction",round(percentageFinish,0),sep=","))
    write(paste("Finished survival prediction",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
  }

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

  plotTrain<-as.data.frame(cbind(predTrain,survivedDays[trainingSet],event[trainingSet]))
  for (i in 1:dim(plotTrain)[1]){
    if (plotTrain[i,1]>=thres){
      plotTrain[i,1]<-"Poor"
    } else {
      plotTrain[i,1]<-"Good"
    }  
  }
  plotTrain.surv <- survfit(Surv(plotTrain[,2], plotTrain[,3]) ~ plotTrain[,1], data = plotTrain)
  #plot(plotTrain.surv, lty = 2:3) 

  plotTest[testSet,1]<-predTest
  for (i in testSet){
    if (plotTest[i,1]>=(thres)){
      plotTest[i,1]<-"Poor"
    } else if (plotTest[i,1]<(thres)){
      plotTest[i,1]<-"Good"
    }
  }
  print(plotTest)
  #plotTest.surv <- survfit(Surv(V2, V3) ~ X1, data = plotTest)
  #plot(plotTest.surv, lty = 2:3)
}


# plot by ggsurv
if (nFolds>1){ # kfold, LOOCV
  plotTestGGsurv<-plotTest  
} else { # random, batch
  plotTestGGsurv<-plotTest[testSet,]  
}
plotTestGGsurv[,2]<-plotTestGGsurv[,2]/12
plotTest.surv <- survfit(Surv(plotTestGGsurv[,2], plotTestGGsurv[,3]) ~ plotTestGGsurv[,1], data = plotTestGGsurv)


survdiffTrain<-survdiff(Surv(plotTrain[,2], plotTrain[,3]) ~ plotTrain[,1], plotTrain)
pTrain<-1-pchisq(survdiffTrain$chisq, length(survdiffTrain$n)-1)

if (length(unique(plotTestGGsurv[,1]))>1){
  survdiffTest<-survdiff(Surv(plotTestGGsurv[,2], plotTestGGsurv[,3]) ~ plotTestGGsurv[,1], plotTestGGsurv)
  pTest<-1-pchisq(survdiffTest$chisq, length(survdiffTest$n)-1)
} else {
  pTest<-"NA"
}
print(pTest)
write(pTest,paste("public/sessions/",sessionID,"/pTest.txt",sep=""))

percentageFinish<-percentageFinish+2
print(paste("Finished building survival groups",percentageFinish,sep=","))
write(paste("Finished building survival groups",percentageFinish,sep=","),milestonesFileName,append=T)



# output feature weights
coefFit<-as.data.frame(as.matrix(coef(cv.tr, s=lambdaFit)))
coefFit<-cbind(0,coefFit)
rownames(coefFit)<-make.names(omicsNameFile, unique=T)
coefFitNonZero<-coefFit[coefFit[,2]!=0,]
coefFitNonZero<-coefFitNonZero[order(-coefFitNonZero[,2]),]
coefFitOutput<-rbind(coefFitNonZero,coefFit[coefFit[,2]==0,])
coefFitOutput[,1]<-rownames(coefFitOutput)

write.table(coefFitOutput,paste("public/sessions/",sessionID,"/featureWeights.txt",sep=""),quote=F, sep=",",row.names=F, col.names=F)



# re-define ggsurv
#source("../public/dropbox/analysis/Survival/ggsurv.R")
source("public/RCodes/ggsurv.R")

png(filename=paste("public/sessions/",sessionID,"/survivalOutput.png", sep=""), width=1000, height=500)
suppressMessages(ggsurv(plotTest.surv) + 
  guides(linetype = F) + 
  xlab("Months") + ylab("Probability of Survival") + ggtitle("Kaplan-Meier Curves") +
  theme(plot.title = element_text(size=24,face="bold"),
    legend.text=element_text(size=24),
    legend.title=element_text(size=24),
    axis.text=element_text(size=18),
    axis.title=element_text(size=24,face="bold")) +
  scale_colour_discrete(name = 'Survival Groups', labels=c('Good Prognostic Group', 'Poor Prognostic Group'))
)
# Close and save the PNG file.
dev.off()

percentageFinish<-100
print(paste("Completed!",percentageFinish,sep=","))
write(paste("Completed!",percentageFinish,sep=","),milestonesFileName,append=T)




# send notification e-mail
emailBody<-paste("Your requested analysis is complete. Please visit http://",ipAddress,
                 ":3000/analysis/survival?done=1&tumor_type=",tumorType,
                 "&data_source=",dataType,
                 "&prediction_target=survival",
                 "&partition=",partitionType,
                 "&var1=",minAlpha,
                 "&var2=",maxAlpha,
                 "&var3=",minLambda,
                 "&var4=",maxLambda,
                 "&session_id=",sessionID,
                 "#results to see the detailed results.",sep="")
emailBodyFileName<-paste("public/sessions/",sessionID,"/emailBody.txt",sep="")
write(emailBody,emailBodyFileName)

if (userEmail!=""){
  system(paste("cat ",emailBodyFileName," | mail -s 'Your OASISPRO Analysis is Complete' ",userEmail,sep=""))
}


