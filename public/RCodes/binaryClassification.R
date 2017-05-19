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
userEmail<-sysargs[10]
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


AUCs<-rep(0,12)

print (featureSelectionMethod)
print (numFeatures)
milestonesFileName<-paste("public/sessions/",sessionID,"/milestones.txt",sep="")
percentageFinish<-0
print(paste("Initializing...",percentageFinish,sep=","))
write(paste("Initializing...",percentageFinish,sep=","),milestonesFileName)


# read files
percentageFinish<-1
print(paste("Reading omics files...",percentageFinish,sep=","))
write(paste("Reading omics files...",percentageFinish,sep=","),milestonesFileName)

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
colnames(omicsFile)<-paste(omicsNameFile,"_",colnames(omicsFile),sep="")


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
x.rp.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.ct.prob<-as.data.frame(matrix(ncol=1, nrow=length(Y)))
x.cf.prob<-as.data.frame(matrix(ncol=1, nrow=length(Y)))
x.ip.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.svm.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.svm.l.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.svm.p.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.svm.s.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.dt.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.rf.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.nb.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))
x.nb.l.prob<-as.data.frame(matrix(ncol=2, nrow=length(Y)))

for (i in 1:nFolds) {
  if (nFolds>1){ # kfold or loocv
    trainingSet<-folds$subsets[,1][folds$which != i]
    testSet<-which(!(1:length(Y) %in% trainingSet))
    print(paste("Fold: ",i,sep=""))
  }

  ## feature selection
  YTrain<-Y[trainingSet]
  XTune<-cbind(Xall[trainingSet,],YTrain)
  
  
  if (sum(YTrain==0)<5){
    stop("Too few observations in Group 1 (n < 5). Please reselect.")
  } else if (sum(YTrain==1)<5){
    stop("Too few observations in Group 2 (n < 5). Please reselect.")
  }
  
  if (featureSelectionMethod == "custom") { # custom
    customFile<-read.table(featureSelectionMethod, sep="")
    subsetFeatures<-customFile[,1]
    X<-Xall[,subsetFeatures]
  } else {
    weights<-matrix(data=0, nrow=dim(XTune)[2]-1, ncol=1)
    rownames(weights)<-colnames(XTune)[1:(dim(XTune)[2]-1)]
    featureSelectionStep<-2000
    nFeatureSelectionParts<-length(seq(1,dim(XTune)[2]-1,featureSelectionStep))
    
    for (featureStart in seq(1,dim(XTune)[2]-1,featureSelectionStep)){
      featureEnd<-min((featureStart+featureSelectionStep-1),dim(XTune)[2]-1)
      XTunePart<-XTune[,featureStart:featureEnd]
      print(featureEnd)
      
      if (featureSelectionMethod == "infog"){ # information gain
        weightsPart <- information.gain(YTrain~., XTunePart)
      } else if (featureSelectionMethod == "gainr"){ # gain ratio
        weightsPart <- gain.ratio(YTrain~., XTunePart)
      } else if (featureSelectionMethod == "symu"){ # symmetrical uncertainty
        weightsPart <- symmetrical.uncertainty(YTrain~., XTunePart)
      } else if (featureSelectionMethod == "randf"){ # random forest importance
        weightsPart <- random.forest.importance(YTrain~., XTunePart, importance.type = 1)
      }
      weights[featureStart:featureEnd,1]<-weightsPart[1:dim(weightsPart)[1],1]
      
      
      
      percentageFinish<-percentageFinish+((percentageStep/2)/nFeatureSelectionParts)
      if (nFolds>1){
        print(paste(paste("Fold ",i,": Running feature selection",sep=""),round(percentageFinish,0),sep=","))
        write(paste(paste("Fold ",i,": Running feature selection",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
      } else {
        print(paste("Running feature selection",round(percentageFinish,0),sep=","))
        write(paste("Running feature selection",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
      }
    }
    
    weightOrder<-order(weights,decreasing=T)
    weightsOutput<-cbind(omicsNameFile[weightOrder,1],weights[weightOrder,1])
    subsetFeatures<-cutoff.k(weights, numFeatures)
    X<-Xall[,subsetFeatures]
  }
  
  write.table(weightsOutput,paste("public/sessions/",sessionID,"/featureWeights.txt",sep=""),quote=F, sep=",",row.names=F, col.names=F)


  if (nFolds>1){
    print(paste(paste("Fold ",i,": Finished feature selection",sep=""),round(percentageFinish,0),sep=","))
    write(paste(paste("Fold ",i,": Finished feature selection",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
  } else {
    print(paste("Finished feature selection",round(percentageFinish,0),sep=","))
    write(paste("Finished feature selection",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
  }
  

  # fit recursive partitioning tree
  x.rp <- rpart(Y[trainingSet]~., data=X[trainingSet,])
  x.rp.pred <- predict(x.rp, X[testSet,], type="class")
  x.rp.prob[testSet,] <- predict(x.rp, X[testSet,], type="prob")
  
  # fit conditional inference trees
  x.ct <- ctree(Y[trainingSet]~., data=X[trainingSet,])
  x.ct.pred <- predict(x.ct, X[testSet,])
  x.ct.prob[testSet,] <- 1- unlist(treeresponse(x.ct, X[testSet,]), use.names=F)[seq(1,nrow(X[testSet,])*2,2)]
  
  # fit random forest with conditional inference trees
  x.cf <- cforest(Y[trainingSet]~., data=X[trainingSet,], control = cforest_unbiased(mtry = ncol(X)-2))
  x.cf.pred <- predict(x.cf, newdata=X[testSet,])
  x.cf.prob[testSet,] <- 1- unlist(treeresponse(x.cf, X[testSet,]), use.names=F)[seq(1,nrow(X[testSet,])*2,2)]

  # fit bagging (bootstrap aggregating)
  x.ip <- bagging(Y[trainingSet]~., data=X[trainingSet,])
  x.ip.prob[testSet,] <- predict(x.ip, type="prob", newdata=X[testSet,])

  # fit svm (support vector machine), radial
  XTune<-cbind(X,Y)
  x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                     ranges = list(gamma = 2^(-8:1), cost = 2^(0:4)),
                     tunecontrol = tune.control(sampling = "fix"))
  x.svm.tune
  x.svm <- svm(Y[trainingSet]~., data = X[trainingSet,], cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
  x.svm.prob.tmp <- attr(predict(x.svm, type="prob", newdata=X[testSet,], probability = TRUE), "probabilities")
  x.svm.prob[testSet,] <- x.svm.prob.tmp
  colnames(x.svm.prob) <- colnames(x.svm.prob.tmp)

  # fit svm (support vector machine), linear
  XTune<-cbind(X,Y)
  x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                     ranges = list(gamma = 2^(-8:1), cost = 2^(0:4)),
                     tunecontrol = tune.control(sampling = "fix"), kernel = "linear")
  x.svm.tune
  x.svm.l <- svm(Y[trainingSet]~., data = X[trainingSet,], kernel = "linear", cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
  x.svm.l.prob.tmp <- attr(predict(x.svm.l, type="prob", newdata=X[testSet,], probability = TRUE), "probabilities")
  x.svm.l.prob[testSet,] <- x.svm.l.prob.tmp
  colnames(x.svm.l.prob) <- colnames(x.svm.l.prob.tmp)

  # fit svm (support vector machine), polynomial
  XTune<-cbind(X,Y)
  x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                     ranges = list(gamma = 2^(-8:1), cost = 2^(0:4)),
                     tunecontrol = tune.control(sampling = "fix"), kernel = "polynomial")
  x.svm.tune
  x.svm.p <- svm(Y[trainingSet]~., data = X[trainingSet,], kernel = "polynomial", cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
  x.svm.p.prob.tmp <- attr(predict(x.svm.p, type="prob", newdata=X[testSet,], probability = TRUE), "probabilities")
  x.svm.p.prob[testSet,] <- x.svm.p.prob.tmp
  colnames(x.svm.p.prob) <- colnames(x.svm.p.prob.tmp)

  # fit svm (support vector machine), sigmoid
  XTune<-cbind(X,Y)
  x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                     ranges = list(gamma = 2^(-8:1), cost = 2^(0:4)),
                     tunecontrol = tune.control(sampling = "fix"), kernel = "linear")
  x.svm.tune
  x.svm.s <- svm(Y[trainingSet]~., data = X[trainingSet,], kernel = "linear", cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
  x.svm.s.prob.tmp <- attr(predict(x.svm.l, type="prob", newdata=X[testSet,], probability = TRUE), "probabilities")
  x.svm.s.prob[testSet,] <- x.svm.s.prob.tmp
  colnames(x.svm.s.prob) <- colnames(x.svm.s.prob.tmp)

  # shallow decision tree
  dep=10
  x.dt<-rpart(Y[trainingSet]~.,X[trainingSet,], control=rpart.control(minsplit=0, minbucket=0,cp=-1, maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
  x.dt.prob[testSet,] <- predict(x.dt, type="prob", newdata=X[testSet,], probability = TRUE)

  # fit random forest
  x.rf<-randomForest(X[trainingSet,],Y[trainingSet])
  x.rf.prob[testSet,] <- predict(x.rf, type="prob", newdata=X[testSet,], probability = TRUE)

  # fit naive Bayes
  x.nb<-naiveBayes(X[trainingSet,],Y[trainingSet])
  x.nb.prob[testSet,] <- predict(x.nb, type="raw", newdata=X[testSet,], probability = TRUE)

  # fit naive Bayes with laplace smoothing
  x.nb.l<-naiveBayes(X[trainingSet,],Y[trainingSet], laplace = 3)
  x.nb.l.prob[testSet,] <- predict(x.nb.l, type="raw", newdata=X[testSet,], probability = TRUE)

  percentageFinish<-percentageFinish+(percentageStep/2)
  if (nFolds>1){
    print(paste(paste("Fold ",i,": Finished prediction",sep=""),round(percentageFinish,0),sep=","))
    write(paste(paste("Fold ",i,": Finished prediction",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
  } else {
    print(paste("Finished prediction",round(percentageFinish,0),sep=","))
    write(paste("Finished prediction",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
  }
}

## plot ROC curves
if (nFolds>1){ # kfold, LOOCV
  testSet=1:length(Y)
}

if (length(testSet)<5){
  stop("Too few observations in the test set. Please reselect training / test partition.")
}

# rpart
x.rp.prob.rocr <- prediction(x.rp.prob[testSet,2], Y[testSet])
x.rp.perf <- performance(x.rp.prob.rocr, "tpr","fpr")
performance(x.rp.prob.rocr,"auc")@y.values
AUCs[1]<-unlist(performance(x.rp.prob.rocr,"auc")@y.values)
d1 <- data.frame(x1=x.rp.perf@x.values[[1]], y1=x.rp.perf@y.values[[1]], Methods="Recursive Partitioning Trees")

# ctree
x.ct.prob.rocr <- prediction(x.ct.prob[testSet,], Y[testSet])
x.ct.perf <- performance(x.ct.prob.rocr, "tpr","fpr")
performance(x.ct.prob.rocr,"auc")@y.values
AUCs[2]<-unlist(performance(x.ct.prob.rocr,"auc")@y.values)
d2 <- data.frame(x1=x.ct.perf@x.values[[1]], y1=x.ct.perf@y.values[[1]], Methods="CITs")

# cforest
x.cf.prob.rocr <- prediction(x.cf.prob[testSet,], Y[testSet])
x.cf.perf <- performance(x.cf.prob.rocr, "tpr","fpr")
performance(x.cf.prob.rocr,"auc")@y.values
AUCs[3]<-unlist(performance(x.cf.prob.rocr,"auc")@y.values)
d3 <- data.frame(x1=x.cf.perf@x.values[[1]], y1=x.cf.perf@y.values[[1]], Methods="Random Forest with CITs")

# bagging
x.ip.prob.rocr <- prediction(x.ip.prob[testSet,2], Y[testSet])
x.ip.perf <- performance(x.ip.prob.rocr, "tpr","fpr")
performance(x.ip.prob.rocr,"auc")@y.values
AUCs[4]<-unlist(performance(x.ip.prob.rocr,"auc")@y.values)
d4 <- data.frame(x1=x.ip.perf@x.values[[1]], y1=x.ip.perf@y.values[[1]], Methods="Bagging")

# svm
x.svm.prob.rocr <- prediction(x.svm.prob[testSet,"1"], Y[testSet])
x.svm.perf <- performance(x.svm.prob.rocr, "tpr","fpr")
performance(x.svm.prob.rocr,"auc")@y.values
AUCs[5]<-unlist(performance(x.svm.prob.rocr,"auc")@y.values)
d5 <- data.frame(x1=x.svm.perf@x.values[[1]], y1=x.svm.perf@y.values[[1]], Methods="SVMs with Gaussian Kernel")

# svm, linear
x.svm.l.prob.rocr <- prediction(x.svm.l.prob[testSet,"1"], Y[testSet])
x.svm.l.perf <- performance(x.svm.l.prob.rocr, "tpr","fpr")
performance(x.svm.l.prob.rocr,"auc")@y.values
AUCs[6]<-unlist(performance(x.svm.l.prob.rocr,"auc")@y.values)
d6 <- data.frame(x1=x.svm.l.perf@x.values[[1]], y1=x.svm.l.perf@y.values[[1]], Methods="SVMs with Linear Kernel")

# svm, polynomial
x.svm.p.prob.rocr <- prediction(x.svm.p.prob[testSet,"1"], Y[testSet])
x.svm.p.perf <- performance(x.svm.p.prob.rocr, "tpr","fpr")
performance(x.svm.p.prob.rocr,"auc")@y.values
AUCs[7]<-unlist(performance(x.svm.p.prob.rocr,"auc")@y.values)
d7 <- data.frame(x1=x.svm.p.perf@x.values[[1]], y1=x.svm.p.perf@y.values[[1]], Methods="SVMs with Polynomial Kernel")

# svm, sigmoid
x.svm.s.prob.rocr <- prediction(x.svm.s.prob[testSet,"1"], Y[testSet])
x.svm.s.perf <- performance(x.svm.s.prob.rocr, "tpr","fpr")
performance(x.svm.s.prob.rocr,"auc")@y.values
AUCs[8]<-unlist(performance(x.svm.s.prob.rocr,"auc")@y.values)
d8 <- data.frame(x1=x.svm.s.perf@x.values[[1]], y1=x.svm.s.perf@y.values[[1]], Methods="SVMs with Sigmoid Kernel")

# shallow decision trees
x.dt.prob.rocr <- prediction(x.dt.prob[testSet,2], Y[testSet])
x.dt.perf <- performance(x.dt.prob.rocr, "tpr","fpr")
performance(x.dt.prob.rocr,"auc")@y.values
AUCs[9]<-unlist(performance(x.dt.prob.rocr,"auc")@y.values)
d9 <- data.frame(x1=x.dt.perf@x.values[[1]], y1=x.dt.perf@y.values[[1]], Methods="Shallow Decision Trees")

# random forest
x.rf.prob.rocr <- prediction(x.rf.prob[testSet,2], Y[testSet])
x.rf.perf <- performance(x.rf.prob.rocr, "tpr","fpr")
performance(x.rf.prob.rocr,"auc")@y.values
AUCs[10]<-unlist(performance(x.rf.prob.rocr,"auc")@y.values)
d10 <- data.frame(x1=x.rf.perf@x.values[[1]], y1=x.rf.perf@y.values[[1]], Methods="Random Forest")

# naive Bayes
x.nb.prob.rocr <- prediction(x.nb.prob[testSet,2], Y[testSet])
x.nb.perf <- performance(x.nb.prob.rocr, "tpr","fpr")
performance(x.nb.prob.rocr,"auc")@y.values
AUCs[11]<-unlist(performance(x.nb.prob.rocr,"auc")@y.values)
d11 <- data.frame(x1=x.nb.perf@x.values[[1]], y1=x.nb.perf@y.values[[1]], Methods="NB")

# naive Bayes with Laplace smoothing
x.nb.l.prob.rocr <- prediction(x.nb.l.prob[testSet,2], Y[testSet])
x.nb.l.perf <- performance(x.nb.l.prob.rocr, "tpr","fpr")
performance(x.nb.l.prob.rocr,"auc")@y.values
AUCs[12]<-unlist(performance(x.nb.l.prob.rocr,"auc")@y.values)
d12 <- data.frame(x1=x.nb.l.perf@x.values[[1]], y1=x.nb.l.perf@y.values[[1]], Methods="NB with Laplace Smoothing")


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

## output ROC curves
figureFileName<-paste("public/sessions/",sessionID,"/roc.png",sep="")
png(filename=figureFileName, width=1000, height=500)
ggplot() +
  geom_path(aes(x1, y1, colour=Methods), d1) +
  geom_path(aes(x1, y1, colour=Methods), d2) +
  geom_line(aes(x1, y1, colour=Methods), d3) +
  geom_line(aes(x1, y1, colour=Methods), d4) +
  geom_line(aes(x1, y1, colour=Methods), d5) +
  geom_line(aes(x1, y1, colour=Methods), d6) +
  geom_line(aes(x1, y1, colour=Methods), d7) +
  geom_line(aes(x1, y1, colour=Methods), d8) +
  geom_line(aes(x1, y1, colour=Methods), d9) +
  geom_line(aes(x1, y1, colour=Methods), d10) +
  geom_line(aes(x1, y1, colour=Methods), d11) +
  geom_line(aes(x1, y1, colour=Methods), d12) +
  geom_line(aes(x1, x1), colour="#000000", d12) +
  theme(plot.title = element_text(size=24,face="bold"),legend.text=element_text(size=16), legend.title=element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=24,face="bold")) +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  ggtitle("ROC Curves for Machine Learning Models")

# Close and save the PNG file.
dev.off()

# output AUC
write.table(cbind(c('Recursive Partitioning Trees','Conditional Inference Trees (CITs)','Random Forest with CITs','Bagging','SVMs with Gaussian Kernel','SVMs with Linear Kernel','SVMs with Polynomial Kernel','SVMs with Sigmoid Kernel','Shallow Decision Trees','Random Forest','Naive Bayes','Naive Bayes with Laplace Smoothing'),AUCs),paste("public/sessions/",sessionID,"/AUCs.txt",sep=""),quote=F, sep=",",row.names=F, col.names=F)


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















