# binaryClassification.R
#
# train and test binary classifiers
# plot ROC for classifiers
#
# param 1: tumor type
# param 2: data type (RNAseq, proteomics, etc)
# param 3: prediction target
# param 4: partition type (random, batch)
# param 5: partition seed (random seed, name of batch file)
# param 6: percentage of data in the training set
# param 7: feature selection method
# param 8: number of features
# param 9: session ID
# param 10: user's e-mail
#
# Kun-Hsing Yu
# April 15, 2016

rm(list=ls())
suppressMessages(library(FSelector))
suppressMessages(library(cvTools))
suppressMessages(library(dplyr))
suppressMessages(library(rpart))
suppressMessages(library(party))
suppressMessages(library(ipred))
suppressMessages(library(e1071))
suppressMessages(library(randomForest))
suppressMessages(library(ROCR))
suppressMessages(library(ggplot2))

sysargs<-commandArgs(trailingOnly=TRUE)
print (sysargs)
tumorType<-sysargs[1]
dataType<-sysargs[2]
predictionTarget<-sysargs[3]
partitionType<-sysargs[4]
if (partitionType == "random") { # random
	partitionSeed<-as.numeric(sysargs[5])
	trainingPercentage<-as.numeric(sysargs[6])
	nRandGroups<-20
} else if (partitionType == "batch") { # batch
	partitionSeed<-(-1)
	selectedBatchesFile<-read.table(paste("public/sessions/",sysargs[9],"/",sysargs[5],sep=""), sep="")
	selectedBatches<-selectedBatchesFile[,1]
	trainingPercentage<-(-1)
	nRandGroups<-1
} else if (partitionType == "kfold") {
  nFolds<-as.numeric(sysargs[5])
} 
featureSelectionMethod<-sysargs[7]
numFeatures<-as.numeric(sysargs[8])
sessionID<-sysargs[9]
nCores<-as.numeric(sysargs[10])
userEmail<-sysargs[11]
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

mlParameters<-read.table(paste("public/sessions/",sessionID,"/mlparameters.txt",sep=""), stringsAsFactors=F)

AUCs<-rep(0,10)

print (featureSelectionMethod)
print (numFeatures)
milestonesFileName<-paste("public/sessions/",sessionID,"/milestones.txt",sep="")
percentageFinish<-0
print(paste("Initializing...",percentageFinish,sep=","))
write(paste("Initializing...",percentageFinish,sep=","),milestonesFileName)


# read files
percentageFinish<-1
print(paste("Reading omics files...",percentageFinish,sep=","))
write(paste("Reading omics files...",percentageFinish,sep=","),milestonesFileName,append=T)

omicsFileName<-paste("../data/", tumorType, "_", dataType, ".RData", sep="")
if (file.exists(omicsFileName)) {
  #try(omicsFile<-read.table(omicsFileName, stringsAsFactors=F, sep=","), TRUE)
  try(load(omicsFileName), TRUE)
  if (!exists("omicsFile")) {
    stop ("Omics type not available for this tumor. Please reselect.")
  }
} else {
  stop ("Omics type not available for this tumor. Please reselect.")
}


percentageFinish<-2
print(paste("Finished reading omics file",percentageFinish,sep=","))
write(paste("Finished reading omics file",percentageFinish,sep=","),milestonesFileName,append=T)


omicsIDFile<-read.table(paste("../data/", tumorType, "_", dataType, "_ids.txt", sep=""), stringsAsFactors=F, sep="")
omicsFile<-omicsFile[substr(omicsIDFile[,1],14,15)=="01",]
omicsIDFile<-omicsIDFile[substr(omicsIDFile[,1],14,15)=="01",]
omicsID<-substr(omicsIDFile, 1, 12)

omicsNameFile<-t(read.table(paste("../data/", tumorType, "_", dataType, "_elemid.txt", sep=""), stringsAsFactors=F, sep=","))
#colnames(omicsFile)<-paste(omicsNameFile,"_",colnames(omicsFile),sep="")
colnames(omicsFile)<-paste(colnames(omicsFile),"_",omicsNameFile,sep="")
colnames(omicsFile)<-gsub("-","_",colnames(omicsFile))

clinicalFile<-read.table(paste("../data/", "nationwidechildrens.org_clinical_patient_", tumorType, ".txt", sep=""), stringsAsFactors = F, sep="\t", quote="")
clinical<-clinicalFile[4:dim(clinicalFile)[1],]
colnames(clinical)<-clinicalFile[2,]
clinicalID<-clinical[,2]

IDintersection<-intersect(omicsID, clinicalID)
omics<-omicsFile[(omicsID %in% IDintersection),]
clinical<-clinical[(clinicalID %in% IDintersection),]
intersectIDs<-omicsID[(omicsID %in% IDintersection)]

YFile<-clinical[,predictionTarget]
uniqueValues<-unique(as.factor(YFile))
YFile[substr(YFile,1,2)=="[N"]<-NA
YFile[substr(YFile,nchar(YFile),nchar(YFile))=="X"]<-NA


#threshold<-median(as.numeric(as.factor(YFile)),na.rm=T)
#Y<-YFile[!is.na(YFile)]
#Y<-(as.numeric(as.factor(Y))>threshold)
#Y<-ifelse(Y=="TRUE", 1, 0)
#Y<-as.factor(Y)
Ygroup1<-read.table(paste("public/sessions/",sessionID,"/outcomeLabel1.txt",sep=""),sep="\n",stringsAsFactors=F)
Ygroup2<-read.table(paste("public/sessions/",sessionID,"/outcomeLabel2.txt",sep=""),sep="\n",stringsAsFactors=F)

Y<-YFile[!is.na(YFile)]
useIndex<-which((Y %in% Ygroup1[,1]) | (Y %in% Ygroup2[,1]))
Y[Y %in% Ygroup1[,1]]<-0
Y[Y %in% Ygroup2[,1]]<-1

Y<-Y[useIndex]
Y<-as.factor(Y)
Xall<-omics[!is.na(YFile),]
Xall<-Xall[useIndex,]
iDs<-intersectIDs[!is.na(YFile)]
iDs<-iDs[useIndex]


percentageFinish<-5
print(paste("Finished building clinical matrix",percentageFinish,sep=","))
write(paste("Finished building clinical matrix",percentageFinish,sep=","),milestonesFileName,append=T)



if (partitionType == "random"){ # random
  trainingSet<-sample(length(Y), floor(length(Y)*trainingPercentage), replace=F)
  testSet<-which(!(1:length(Y) %in% trainingSet))
  nFolds<-1
} else if (partitionType == "kfold"){ # k fold
  folds<-cvFolds(length(Y), K=nFolds)
  trainingSet<-1:length(Y)
  testSet<-1:length(Y)
} else if (partitionType == "LOOCV"){ # LOOCV
  nFolds<-length(Y)
  folds<-cvFolds(length(Y), K=length(Y))
  trainingSet<-1:length(Y)
  testSet<-1:length(Y)
} else { # batch
  trainingSet<-which(substr(iDs, 6, 7) %in% selectedBatches)
  testSet<-which(!(1:length(Y) %in% trainingSet))
  nFolds<-1
}

write.table(c(sum(Y[trainingSet]==0),sum(Y[trainingSet]==1)),paste("public/sessions/",sessionID,"/nSamplesTraining.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F)
write.table(c(sum(Y[testSet]==0),sum(Y[testSet]==1)),paste("public/sessions/",sessionID,"/nSamplesTest.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F)

percentageStep<-90/nFolds
x.nb.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.rp.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.ct.prob<-as.data.frame(matrix(ncol=1, nrow=length(Y)))
x.ip.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.rf.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.cf.prob<-as.data.frame(matrix(ncol=1, nrow=length(Y)))
x.svm.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.svm.l.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.svm.p.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.svm.s.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))



if (nFolds==1 | nCores==1){
  i=1
  source("public/RCodes/parallelBinaryClassificationHelper.R")
} else {
  save.image(file=paste("public/sessions/",sessionID,"/parallelInput.RData",sep=""))
  
  pids<-NULL
  for (i in 1:nFolds){
    print (paste("UsedCores:",length(pids)))
    print (pids)
    while (length(pids)>=nCores) { # wait for some process to complete
      for (processIndex in 1:length(pids)) {
        if (processIndex>length(pids)) { # after removal of pid, the length(pids) in the for loop will not be updated
          break
        }
        doneFlag<-system(paste("if [ $(ps -p ",pids[processIndex]," -o pid= | wc -l) -eq 0 ]; then echo 'done'; else echo 'not_done'; fi",sep=""),intern=T)
        if (doneFlag=="done") {
          print (paste("Done:",pids[processIndex]))
          pids<-pids[-processIndex]
        } else {
          system("sleep 1")
        }
      }
    }
    system(paste("Rscript public/RCodes/parallelBinaryClassificationHelper.R ",nCores," ",sessionID," ",i,"& echo $! > ","public/sessions/",sessionID,"/parallelPS.txt",sep=""))
    pid<-read.table(paste("public/sessions/",sessionID,"/parallelPS.txt",sep=""),stringsAsFactors=F)
    pids<-c(pids,pid[1,1])
    print (pids)
    percentageFinish<-percentageFinish+(percentageStep/2)
    print(paste(paste("Fold ",i,": Building machine learning models (parallel model)",sep=""),round(percentageFinish,0),sep=","))
    write(paste(paste("Fold ",i,": Building machine learning models (parallel model)",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
  }
  #system("wait")
  for (i in 1:nFolds){
    waitFlag<-1
    while (waitFlag==1){
      if (system(paste("[ -f ","public/sessions/",sessionID,"/parallelOutputFold",i,".RData ] && echo '1' || echo '0'",sep=""),intern=T)=='1'){
        load(file=paste("public/sessions/",sessionID,"/parallelOutputFold",i,".RData",sep=""))
        #system(paste("rm public/sessions/",sessionID,"/parallelOutputFold",i,".RData",sep=""))
        testSet<-parallelOutput[[1]]
        weightsOutput<-parallelOutput[[2]]
        x.nb.prob[testSet,]<-parallelOutput[[3]][testSet,]
        x.rp.prob[testSet,]<-parallelOutput[[4]][testSet,]
        x.ct.prob[testSet,]<-parallelOutput[[5]][testSet,]
        x.ip.prob[testSet,]<-parallelOutput[[6]][testSet,]
        x.rf.prob[testSet,]<-parallelOutput[[7]][testSet,]
        x.cf.prob[testSet,]<-parallelOutput[[8]][testSet,]
        x.svm.prob[testSet,]<-parallelOutput[[9]][testSet,]
        x.svm.l.prob[testSet,]<-parallelOutput[[10]][testSet,]
        x.svm.p.prob[testSet,]<-parallelOutput[[11]][testSet,]
        x.svm.s.prob[testSet,]<-parallelOutput[[12]][testSet,]
        x.rp<-parallelOutput[[13]]
        #percentageFinish<-parallelOutput[[14]]
        waitFlag<-0
        if (nFolds>1){
          write.table(c(paste("Fold: ",i,sep="")),paste("public/sessions/",sessionID,"/featureWeightsAll.txt",sep=""),quote=F, sep=",",row.names=F, col.names=F,append=T)
        }
        write.table(weightsOutput,paste("public/sessions/",sessionID,"/featureWeightsAll.txt",sep=""),quote=F, sep=",",row.names=F, col.names=F,append=T)
        
        percentageFinish<-percentageFinish+(percentageStep/2)
        print(paste(paste("Fold ",i,": Finished prediction",sep=""),round(percentageFinish,0),sep=","))
        write(paste(paste("Fold ",i,": Finished prediction",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
      } else {
        system("sleep 1")
      }
    }
    colnames(x.svm.prob)<-colnames(parallelOutput[[9]][testSet,])
    colnames(x.svm.l.prob)<-colnames(parallelOutput[[10]][testSet,])
    colnames(x.svm.p.prob)<-colnames(parallelOutput[[11]][testSet,])
    colnames(x.svm.s.prob)<-colnames(parallelOutput[[12]][testSet,])
  }
}

write.table(weightsOutput,paste("public/sessions/",sessionID,"/featureWeights.txt",sep=""),quote=F, sep=",",row.names=F, col.names=F)

## plot ROC curves
if (nFolds>1){ # kfold, LOOCV
  testSet=1:length(Y)
}

if (length(testSet)<5){
  stop("Too few observations in the test set. Please reselect training / test partition.")
}


rocPlot<-ggplot()
# naive Bayes
if (mlParameters[1,1]=="on"){
  x.nb.prob.rocr <- prediction(x.nb.prob[testSet,2], Y[testSet])
  x.nb.perf <- performance(x.nb.prob.rocr, "tpr","fpr")
  performance(x.nb.prob.rocr,"auc")@y.values
  AUCs[1]<-unlist(performance(x.nb.prob.rocr,"auc")@y.values)
  d1 <- data.frame(x1=x.nb.perf@x.values[[1]], y1=x.nb.perf@y.values[[1]], Methods="Naive Bayes")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d1)
}

# rpart
if (mlParameters[2,1]=="on"){
  x.rp.prob.rocr <- prediction(x.rp.prob[testSet,2], Y[testSet])
  x.rp.perf <- performance(x.rp.prob.rocr, "tpr","fpr")
  performance(x.rp.prob.rocr,"auc")@y.values
  AUCs[2]<-unlist(performance(x.rp.prob.rocr,"auc")@y.values)
  d2 <- data.frame(x1=x.rp.perf@x.values[[1]], y1=x.rp.perf@y.values[[1]], Methods="Recursive Partitioning Trees")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d2)
}

# ctree
if (mlParameters[3,1]=="on"){
  x.ct.prob.rocr <- prediction(x.ct.prob[testSet,], Y[testSet])
  x.ct.perf <- performance(x.ct.prob.rocr, "tpr","fpr")
  performance(x.ct.prob.rocr,"auc")@y.values
  AUCs[3]<-unlist(performance(x.ct.prob.rocr,"auc")@y.values)
  d3 <- data.frame(x1=x.ct.perf@x.values[[1]], y1=x.ct.perf@y.values[[1]], Methods="Conditional Inference Trees (CITs)")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d3)
}

# bagging
if (mlParameters[4,1]=="on"){
  x.ip.prob.rocr <- prediction(x.ip.prob[testSet,2], Y[testSet])
  x.ip.perf <- performance(x.ip.prob.rocr, "tpr","fpr")
  performance(x.ip.prob.rocr,"auc")@y.values
  AUCs[4]<-unlist(performance(x.ip.prob.rocr,"auc")@y.values)
  d4 <- data.frame(x1=x.ip.perf@x.values[[1]], y1=x.ip.perf@y.values[[1]], Methods="Bagging")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d4)
}

# random forest
if (mlParameters[5,1]=="on"){
  x.rf.prob.rocr <- prediction(x.rf.prob[testSet,2], Y[testSet])
  x.rf.perf <- performance(x.rf.prob.rocr, "tpr","fpr")
  performance(x.rf.prob.rocr,"auc")@y.values
  AUCs[5]<-unlist(performance(x.rf.prob.rocr,"auc")@y.values)
  d5 <- data.frame(x1=x.rf.perf@x.values[[1]], y1=x.rf.perf@y.values[[1]], Methods="Random Forest")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d5)
}

# cforest
if (mlParameters[6,1]=="on"){
  x.cf.prob.rocr <- prediction(x.cf.prob[testSet,], Y[testSet])
  x.cf.perf <- performance(x.cf.prob.rocr, "tpr","fpr")
  performance(x.cf.prob.rocr,"auc")@y.values
  AUCs[6]<-unlist(performance(x.cf.prob.rocr,"auc")@y.values)
  d6 <- data.frame(x1=x.cf.perf@x.values[[1]], y1=x.cf.perf@y.values[[1]], Methods="Random Forest with CITs")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d6)
}

# svm
if (mlParameters[7,1]=="on"){
  x.svm.prob.rocr <- prediction(x.svm.prob[testSet,"1"], Y[testSet])
  x.svm.perf <- performance(x.svm.prob.rocr, "tpr","fpr")
  performance(x.svm.prob.rocr,"auc")@y.values
  AUCs[7]<-unlist(performance(x.svm.prob.rocr,"auc")@y.values)
  d7 <- data.frame(x1=x.svm.perf@x.values[[1]], y1=x.svm.perf@y.values[[1]], Methods="SVMs with Gaussian Kernel")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d7)
}

# svm, linear
if (mlParameters[8,1]=="on"){
  x.svm.l.prob.rocr <- prediction(x.svm.l.prob[testSet,"1"], Y[testSet])
  x.svm.l.perf <- performance(x.svm.l.prob.rocr, "tpr","fpr")
  performance(x.svm.l.prob.rocr,"auc")@y.values
  AUCs[8]<-unlist(performance(x.svm.l.prob.rocr,"auc")@y.values)
  d8 <- data.frame(x1=x.svm.l.perf@x.values[[1]], y1=x.svm.l.perf@y.values[[1]], Methods="SVMs with Linear Kernel")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d8)
}

# svm, polynomial
if (mlParameters[9,1]=="on"){
  x.svm.p.prob.rocr <- prediction(x.svm.p.prob[testSet,"1"], Y[testSet])
  x.svm.p.perf <- performance(x.svm.p.prob.rocr, "tpr","fpr")
  performance(x.svm.p.prob.rocr,"auc")@y.values
  AUCs[9]<-unlist(performance(x.svm.p.prob.rocr,"auc")@y.values)
  d9 <- data.frame(x1=x.svm.p.perf@x.values[[1]], y1=x.svm.p.perf@y.values[[1]], Methods="SVMs with Polynomial Kernel")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d9)
}

# svm, sigmoid
if (mlParameters[10,1]=="on"){
  x.svm.s.prob.rocr <- prediction(x.svm.s.prob[testSet,"1"], Y[testSet])
  x.svm.s.perf <- performance(x.svm.s.prob.rocr, "tpr","fpr")
  performance(x.svm.s.prob.rocr,"auc")@y.values
  AUCs[10]<-unlist(performance(x.svm.s.prob.rocr,"auc")@y.values)
  d10 <- data.frame(x1=x.svm.s.perf@x.values[[1]], y1=x.svm.s.perf@y.values[[1]], Methods="SVMs with Sigmoid Kernel")
  rocPlot<-rocPlot + geom_path(aes(x1, y1, colour=Methods), d10)
}


percentageFinish<-percentageFinish+2
print(paste("Finished model evaluation",round(percentageFinish,0),sep=","))
write(paste("Finished model evaluation",round(percentageFinish,0),sep=","),milestonesFileName,append=T)


## output decision tree
figureFileName<-paste("public/sessions/",sessionID,"/decisionTree.png",sep="")
png(filename=figureFileName, width=1000, height=500)
if (dim(x.rp$frame)[1]>1) {
  plot(x.rp, uniform=TRUE)
  text(x.rp, use.n=TRUE, all=TRUE, cex=.8)
} else {
  plot(0,type='n',axes=FALSE,xlab="The decision tree is just a root.",ylab="")
}
dev.off()

rocPlot<-rocPlot + geom_abline(intercept=0, slope=1)

rocPlot<-rocPlot + 
  theme(plot.title = element_text(size=24,face="bold"),legend.text=element_text(size=16), legend.title=element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=24,face="bold")) +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  ggtitle("ROC Curves for Machine Learning Models")


## output ROC curves
figureFileName<-paste("public/sessions/",sessionID,"/roc.png",sep="")
png(filename=figureFileName, width=1000, height=500)
print(rocPlot)
# Close and save the PNG file.
dev.off()

# output AUC
write.table(cbind(c("Naive Bayes","Recursive Partitioning Trees","Conditional Inference Trees (CITs)","Bagging","Breiman's Random Forest","Random Forest with CITs","SVMs with Gaussian Kernel","SVMs with Linear Kernel","SVMs with Polynomial Kernel","SVMs with Sigmoid Kernel"),AUCs),paste("public/sessions/",sessionID,"/AUCs.txt",sep=""),quote=F, sep=",",row.names=F, col.names=F)


percentageFinish<-100
print(paste("Completed!",percentageFinish,sep=","))
write(paste("Completed!",percentageFinish,sep=","),milestonesFileName,append=T)


# send notification e-mail
emailBody<-paste("Your requested analysis is complete. Please visit http://",ipAddress,
                 "/analysis/stage?done=1&tumor_type=",tumorType,
                 "&data_source=",dataType,
                 "&prediction_target=",predictionTarget,
                 "&partition=",partitionType,
                 "&var1=",featureSelectionMethod,
                 "&var2=",numFeatures,
                 "&var3=",userEmail,
                 "&var4=0&session_id=",sessionID,
                 "#results to see the detailed results.",sep="")
emailBodyFileName<-paste("public/sessions/",sessionID,"/emailBody.txt",sep="")
write(emailBody,emailBodyFileName)

if (userEmail!=""){
  system(paste("cat ",emailBodyFileName," | mail -s 'Your OASISPRO Analysis is Complete' ",userEmail,sep=""))
}















