# parallelElasticNetCoxHelper.R
#
# train and test binary classifiers
# plot ROC for classifiers
#
# param 1: number of CPU cores
# param 2: sessionID
# param 3: case index
#
# Kun-Hsing Yu
# March 28, 2016


suppressMessages(library(FSelector))
suppressMessages(library(ggplot2))
suppressMessages(library(survival))
suppressMessages(library(glmnet))
suppressMessages(library(rms))
suppressMessages(library(GGally))
suppressMessages(library(dplyr))
suppressMessages(library(cvTools))


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

for (fold in iterStart:iterEnd){
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
  if (nCores==1){
    if (nFolds>1){
      print(paste(paste("Fold ",fold,": Finished regularization",sep=""),round(percentageFinish,0),sep=","))
      write(paste(paste("Fold ",fold,": Finished regularization",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
    } else {
      print(paste("Finished regularization",round(percentageFinish,0),sep=","))
      write(paste("Finished regularization",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
    }
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
  if (nCores==1){
    if (nFolds>1){
      print(paste(paste("Fold ",fold,": Finished survival prediction",sep=""),round(percentageFinish,0),sep=","))
      write(paste(paste("Fold ",fold,": Finished survival prediction",sep=""),round(percentageFinish,0),sep=","),milestonesFileName,append=T)
    } else {
      print(paste("Finished survival prediction",round(percentageFinish,0),sep=","))
      write(paste("Finished survival prediction",round(percentageFinish,0),sep=","),milestonesFileName,append=T)
    }
  }
  
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
  
  plotTest[testSet,1]<-predTest
  for (i in testSet){
    if (plotTest[i,1]>=(thres)){
      plotTest[i,1]<-"Poor"
    } else {
      plotTest[i,1]<-"Good"
    }
  }
  print(plotTest)
  
  # output feature weights
  coefFit<-as.data.frame(as.matrix(coef(cv.tr, s=lambdaFit)))
  coefFit<-cbind(0,coefFit)
  coefFit[,1]<-rownames(coefFit)
  coefFit[1:length(omicsNameFile),1]<-omicsNameFile[,1]
  coefFitNonZero<-coefFit[coefFit[,2]!=0,]
  coefFitNonZero<-coefFitNonZero[order(-coefFitNonZero[,2]),]
  coefFitOutput<-rbind(coefFitNonZero,coefFit[coefFit[,2]==0,])
  
  parallelOutput<-list(testSet,plotTest,cv.tr,lambdaFit,coefFitOutput,percentageFinish)
  save(parallelOutput,file=paste("public/sessions/",sessionID,"/parallelOutputFold",fold,".RData",sep=""))
}
