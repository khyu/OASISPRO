# parallelBinaryClassificationHelper.R
#
# parallelize binary classification tasks
# plot ROC for classifiers
#
# param 1: number of CPU cores
# param 2: sessionID
# param 3: case index
#
# Kun-Hsing Yu
# April 15, 2016

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


# check if this run is a part of parallel sessions
if (!exists("nCores")){
  sysargs<-commandArgs(trailingOnly=TRUE)
  print (sysargs)
  nCores<-as.numeric(sysargs[1])
  sessionID<-sysargs[2]
  i<-as.numeric(sysargs[3]) # fold_i
  load(file=paste("public/sessions/",sessionID,"/parallelInput.RData",sep=""))
  iterStart<-i
  iterEnd<-i
} else {
  iterStart<-1
  iterEnd<-nFolds
}

for (i in iterStart:iterEnd){
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
      #print(featureEnd)
      
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
      if (nCores==1){
        if (nFolds>1){
          print(paste(paste("Fold ",i,": Running feature selection",sep=""),round(percentageFinish,0),sep=","))
          write(paste(paste("Fold ",i,": Running feature selection",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
        } else {
          print(paste("Running feature selection",round(percentageFinish,0),sep=","))
          write(paste("Running feature selection",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
        }
      }
    }
    
    weightOrder<-order(weights,decreasing=T)
    weightsOutput<-cbind(omicsNameFile[weightOrder,1],weights[weightOrder,1])
    subsetFeatures<-cutoff.k(weights, numFeatures)
    X<-Xall[,subsetFeatures]
  }
  if (nCores==1){
    if (nFolds>1){
      write.table(c(paste("Fold: ",i,sep="")),paste("public/sessions/",sessionID,"/featureWeightsAll.txt",sep=""),quote=F, sep=",",row.names=F, col.names=F,append=T)
    }
    write.table(weightsOutput,paste("public/sessions/",sessionID,"/featureWeightsAll.txt",sep=""),quote=F, sep=",",row.names=F, col.names=F,append=T)
  }
  
  
  if (nCores==1){
    if (nFolds>1){
      print(paste(paste("Fold ",i,": Finished feature selection",sep=""),round(percentageFinish,0),sep=","))
      write(paste(paste("Fold ",i,": Finished feature selection",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
    } else {
      print(paste("Finished feature selection",round(percentageFinish,0),sep=","))
      write(paste("Finished feature selection",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
    }
  }
  
  # fit naive Bayes with/without laplace smoothing
  if (mlParameters[1,1]=="on"){
    x.nb<-naiveBayes(X[trainingSet,],Y[trainingSet], laplace = as.numeric(mlParameters[11,1]))
    x.nb.prob[testSet,] <- predict(x.nb, type="raw", newdata=X[testSet,], probability = TRUE)
  }
  
  # fit recursive partitioning tree
  if (mlParameters[2,1]=="on"){
    x.rp <- rpart(Y[trainingSet]~., data=X[trainingSet,], control=rpart.control(maxdepth=as.numeric(mlParameters[12,1]), cp=as.numeric(mlParameters[13,1])))
    x.rp.pred <- predict(x.rp, X[testSet,], type="class")
    x.rp.prob[testSet,] <- predict(x.rp, X[testSet,], type="prob")
  }
  
  # fit conditional inference trees
  if (mlParameters[3,1]=="on"){
    x.ct <- ctree(Y[trainingSet]~., data=X[trainingSet,], control=ctree_control(maxdepth=as.numeric(mlParameters[14,1])))
    x.ct.pred <- predict(x.ct, X[testSet,])
    x.ct.prob[testSet,] <- 1- unlist(treeresponse(x.ct, X[testSet,]), use.names=F)[seq(1,nrow(X[testSet,])*2,2)]
  }
  
  # fit bagging (bootstrap aggregating)
  if (mlParameters[4,1]=="on"){
    x.ip <- bagging(Y[trainingSet]~., data=X[trainingSet,], nbagg=as.numeric(mlParameters[15,1]))
    x.ip.prob[testSet,] <- predict(x.ip, type="prob", newdata=X[testSet,])
  }
  
  # fit random forest
  if (mlParameters[5,1]=="on"){
    x.rf<-randomForest(X[trainingSet,],Y[trainingSet],ntree=as.numeric(mlParameters[16,1]))
    x.rf.prob[testSet,] <- predict(x.rf, type="prob", newdata=X[testSet,], probability = TRUE)
  }
  
  # fit random forest with conditional inference trees
  if (mlParameters[6,1]=="on"){
    x.cf <- cforest(Y[trainingSet]~., data=X[trainingSet,], control = cforest_unbiased(ntree=as.numeric(mlParameters[17,1])))
    x.cf.pred <- predict(x.cf, newdata=X[testSet,])
    x.cf.prob[testSet,] <- 1- unlist(treeresponse(x.cf, X[testSet,]), use.names=F)[seq(1,nrow(X[testSet,])*2,2)]
  }
  
  # fit svm (support vector machine), radial
  if (mlParameters[7,1]=="on"){
    XTune<-cbind(X,Y)
    x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                       ranges = list(gamma = 2^(-8:1), cost = 2^(as.numeric(mlParameters[18,1]):as.numeric(mlParameters[19,1]))),
                       tunecontrol = tune.control(sampling = "fix"))
    #x.svm.tune
    x.svm <- svm(Y[trainingSet]~., data = X[trainingSet,], cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
    x.svm.prob.tmp <- attr(predict(x.svm, type="prob", newdata=X[testSet,], probability = TRUE), "probabilities")
    x.svm.prob[testSet,] <- x.svm.prob.tmp
    colnames(x.svm.prob) <- colnames(x.svm.prob.tmp)
  }
  
  # fit svm (support vector machine), linear
  if (mlParameters[8,1]=="on"){
    XTune<-cbind(X,Y)
    x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                       ranges = list(gamma = 2^(-8:1), cost = 2^(as.numeric(mlParameters[20,1]):as.numeric(mlParameters[21,1]))),
                       tunecontrol = tune.control(sampling = "fix"), kernel = "linear")
    #x.svm.tune
    x.svm.l <- svm(Y[trainingSet]~., data = X[trainingSet,], kernel = "linear", cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
    x.svm.l.prob.tmp <- attr(predict(x.svm.l, type="prob", newdata=X[testSet,], probability = TRUE), "probabilities")
    x.svm.l.prob[testSet,] <- x.svm.l.prob.tmp
    colnames(x.svm.l.prob) <- colnames(x.svm.l.prob.tmp)
  }
  
  # fit svm (support vector machine), polynomial
  if (mlParameters[9,1]=="on"){
    XTune<-cbind(X,Y)
    x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                       ranges = list(gamma = 2^(-8:1), cost = 2^(as.numeric(mlParameters[22,1]):as.numeric(mlParameters[23,1]))),
                       tunecontrol = tune.control(sampling = "fix"), kernel = "polynomial")
    #x.svm.tune
    x.svm.p <- svm(Y[trainingSet]~., data = X[trainingSet,], kernel = "polynomial", cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
    x.svm.p.prob.tmp <- attr(predict(x.svm.p, type="prob", newdata=X[testSet,], probability = TRUE), "probabilities")
    x.svm.p.prob[testSet,] <- x.svm.p.prob.tmp
    colnames(x.svm.p.prob) <- colnames(x.svm.p.prob.tmp)
  }
  
  # fit svm (support vector machine), sigmoid
  if (mlParameters[10,1]=="on"){
    XTune<-cbind(X,Y)
    x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                       ranges = list(gamma = 2^(-8:1), cost = 2^(as.numeric(mlParameters[24,1]):as.numeric(mlParameters[25,1]))),
                       tunecontrol = tune.control(sampling = "fix"), kernel = "linear")
    #x.svm.tune
    x.svm.s <- svm(Y[trainingSet]~., data = X[trainingSet,], kernel = "linear", cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
    x.svm.s.prob.tmp <- attr(predict(x.svm.l, type="prob", newdata=X[testSet,], probability = TRUE), "probabilities")
    x.svm.s.prob[testSet,] <- x.svm.s.prob.tmp
    colnames(x.svm.s.prob) <- colnames(x.svm.s.prob.tmp)
  }
  
  
  percentageFinish<-percentageFinish+(percentageStep/2)
  if (nCores==1){
    if (nFolds>1){
      print(paste(paste("Fold ",i,": Finished prediction",sep=""),round(percentageFinish,0),sep=","))
      write(paste(paste("Fold ",i,": Finished prediction",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
    } else {
      print(paste("Finished prediction",round(percentageFinish,0),sep=","))
      write(paste("Finished prediction",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
    }
  }
  
  if (nCores>1){
    parallelOutput<-list(testSet,weightsOutput,x.nb.prob,x.rp.prob,x.ct.prob,x.ip.prob,x.rf.prob,x.cf.prob,x.svm.prob,x.svm.l.prob,x.svm.p.prob,x.svm.s.prob,x.rp,percentageFinish)
    save(parallelOutput,file=paste("public/sessions/",sessionID,"/parallelOutputFold",i,".RData",sep=""))
  }
}





