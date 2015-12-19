# proteomicsStageOverall.R
#
# use information gain / gain ratio / symmetric uncertainty to select features
# plot ROC for classifiers
#
# Kun-Hsing Yu
# 2014.11.19

rm(list=ls())
library(FSelector)
library(ggplot2)
setwd("~/Dropbox/CS273AFinalProject/analysis/StageOverall/")

nRandGroups<-20
maxFeatures<-15
numClassifiers<-12
trainFraction<-0.7
AUCs<-array(0,dim=c(nRandGroups+1,maxFeatures,numClassifiers))

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

for (randGroup in 1:nRandGroups){
  for (numFeatures in 2:maxFeatures){
    # training, test set 1: 1 to 321, 322 to 459
    #trainingSet<-1:floor(length(Y)*trainFraction)
    #testSet<-(floor(length(Y)*trainFraction)+1):length(Y)
    # training, test set 2: random 0.7, 0.3
    trainingSet<-sample(length(Y), floor(length(Y)*trainFraction), replace=F)
    testSet<-which(!(1:length(Y) %in% trainingSet))
    
    ## feature selection
    YTrain<-Y[trainingSet]
    XTune<-cbind(X[trainingSet,],YTrain)
    #numFeatures=20
    if (FALSE){
      # information gain
      figureFileName<-paste("rocInfoGainModels.png",numFeatures,"Features.png",sep="")
      weights <- information.gain(YTrain~., XTune)
        #X<-X[,weights>0.155]
      subsetFeatures<-cutoff.k(weights, numFeatures)
      X<-X[,subsetFeatures]
    }
      # gain ratio
      figureFileName<-paste("rocGainRatioModels.png",numFeatures,"Features.png",sep="")
      weights <- gain.ratio(YTrain~., XTune)
        #X<-X[,(weights>0.22 & is.na(weights)==F)]
      subsetFeatures<-cutoff.k(weights, numFeatures)
      X<-X[,subsetFeatures]
    if (FALSE){
      # symmetrical uncertainty
      figureFileName<-paste("rocSymUncModels.png",numFeatures,"Features.png",sep="")
      weights <- symmetrical.uncertainty(YTrain~., XTune)
        #X<-X[,weights>0.21]
      subsetFeatures<-cutoff.k(weights, numFeatures)
      X<-X[,subsetFeatures]
    }
    if (FALSE){
      # consistency
      figureFileName<-paste("rocConsistModels.png",numFeatures,"Features.png",sep="")
      weights <- consistency(YTrain~., XTune)
      X<-X[,weights]
    
      # chi-squared
      figureFileName<-paste("rocChiSqModels.png",numFeatures,"Features.png",sep="")
      weights <- chi.squared(YTrain~., XTune)
      subsetFeatures<-cutoff.k(weights, numFeatures)
      X<-X[,subsetFeatures]
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
      # forward feature selection from matlab
      figureFileName<-"rocForwardModels.png"
      X<-X[,c(15339,3003,2405,3025,5092,6668,1215,13913,7357,14037,12315,2498,15470,3560,3169,14077,12462,6,15,4367)]
    }
    # create model using recursive partitioning on the training data set
    library(rpart)
    x.rp <- rpart(Y[trainingSet]~., data=X[trainingSet,])
    # predict classes for the evaluation data set
    x.rp.pred <- predict(x.rp, X[testSet,], type="class")
    # score the evaluation data set (extract the probabilities)
    x.rp.prob <- predict(x.rp, X[testSet,], type="prob")
    # To view the decision tree, uncomment this line.
    # plot(x.rp, main="Decision tree created using rpart")
    
    # create model using conditional inference trees
    require(party)
    x.ct <- ctree(Y[trainingSet]~., data=X[trainingSet,])
    x.ct.pred <- predict(x.ct, X[testSet,])
    x.ct.prob <-  1- unlist(treeresponse(x.ct, X[testSet,]), use.names=F)[seq(1,nrow(X[testSet,])*2,2)]
    # To view the decision tree, uncomment this line.
    # plot(x.ct, main="Decision tree created using condition inference trees")
    
    # create model using random forest and bagging ensemble using conditional inference trees
    x.cf <- cforest(Y[trainingSet]~., data=X[trainingSet,], control = cforest_unbiased(mtry = ncol(X)-2))
    x.cf.pred <- predict(x.cf, newdata=X[testSet,])
    x.cf.prob <-  1- unlist(treeresponse(x.cf, X[testSet,]), use.names=F)[seq(1,nrow(X[testSet,])*2,2)]
    
    # create model using bagging (bootstrap aggregating)
    require(ipred)
    x.ip <- bagging(Y[trainingSet]~., data=X[trainingSet,])
    x.ip.prob <- predict(x.ip, type="prob", newdata=X[testSet,])
    
    # create model using svm (support vector machine), radial
    require(e1071)
    # svm requires tuning
    XTune<-cbind(X,Y)
    x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                       ranges = list(gamma = 2^(-8:1), cost = 2^(0:4)),
                       tunecontrol = tune.control(sampling = "fix"))
    # display the tuning results (in text format)
    x.svm.tune
    # If the tuning results are on the margin of the parameters (e.g., gamma = 2^-8), 
    # then widen the parameters.
    # I manually copied the cost and gamma from console messages above to parameters below.
    x.svm <- svm(Y[trainingSet]~., data = X[trainingSet,], cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
    x.svm.prob <- predict(x.svm, type="prob", newdata=X[testSet,], probability = TRUE)
    
    
    # create model using svm (support vector machine), linear
    # svm requires tuning
    XTune<-cbind(X,Y)
    x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                       ranges = list(gamma = 2^(-8:1), cost = 2^(0:4)),
                       tunecontrol = tune.control(sampling = "fix"), kernel = "linear")
    # display the tuning results (in text format)
    x.svm.tune
    # If the tuning results are on the margin of the parameters (e.g., gamma = 2^-8), 
    # then widen the parameters.
    # I manually copied the cost and gamma from console messages above to parameters below.
    x.svm.l <- svm(Y[trainingSet]~., data = X[trainingSet,], kernel = "linear", cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
    x.svm.l.prob <- predict(x.svm.l, type="prob", newdata=X[testSet,], probability = TRUE)
    
    
    # create model using svm (support vector machine), polynomial
    # svm requires tuning
    XTune<-cbind(X,Y)
    x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                       ranges = list(gamma = 2^(-8:1), cost = 2^(0:4)),
                       tunecontrol = tune.control(sampling = "fix"), kernel = "polynomial")
    # display the tuning results (in text format)
    x.svm.tune
    # If the tuning results are on the margin of the parameters (e.g., gamma = 2^-8), 
    # then widen the parameters.
    # I manually copied the cost and gamma from console messages above to parameters below.
    x.svm.p <- svm(Y[trainingSet]~., data = X[trainingSet,], kernel = "polynomial", cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
    x.svm.p.prob <- predict(x.svm.p, type="prob", newdata=X[testSet,], probability = TRUE)
    
    
    # create model using svm (support vector machine), sigmoid
    # svm requires tuning
    XTune<-cbind(X,Y)
    x.svm.tune <- tune(svm, Y~., data = XTune[trainingSet,],
                       ranges = list(gamma = 2^(-8:1), cost = 2^(0:4)),
                       tunecontrol = tune.control(sampling = "fix"), kernel = "sigmoid")
    # display the tuning results (in text format)
    x.svm.tune
    # If the tuning results are on the margin of the parameters (e.g., gamma = 2^-8), 
    # then widen the parameters.
    # I manually copied the cost and gamma from console messages above to parameters below.
    x.svm.s <- svm(Y[trainingSet]~., data = X[trainingSet,], kernel = "sigmoid", cost=x.svm.tune$best.parameters[2], gamma=x.svm.tune$best.parameters[1], probability = TRUE)
    x.svm.s.prob <- predict(x.svm.s, type="prob", newdata=X[testSet,], probability = TRUE)
    
    # decision tree
    dep=10
    x.dt<-rpart(Y[trainingSet]~.,X[trainingSet,], control=rpart.control(minsplit=0, minbucket=0,cp=-1, maxcompete=0, maxsurrogate=0, usesurrogate=0, xval=0,maxdepth=dep))
    x.dt.prob <- predict(x.dt, type="prob", newdata=X[testSet,], probability = TRUE)
    
    # random forest
    library(randomForest)
    x.rf<-randomForest(X[trainingSet,],Y[trainingSet])
    x.rf.prob <- predict(x.rf, type="prob", newdata=X[testSet,], probability = TRUE)
    
    # naive Bayes
    #library(e1071)
    x.nb<-naiveBayes(X[trainingSet,],Y[trainingSet])
    x.nb.prob <- predict(x.nb, type="raw", newdata=X[testSet,], probability = TRUE)
    
    # naive Bayes with laplace smoothing
    #library(e1071)
    x.nb.l<-naiveBayes(X[trainingSet,],Y[trainingSet], laplace = 3)
    x.nb.l.prob <- predict(x.nb.l, type="raw", newdata=X[testSet,], probability = TRUE)
    
    ##
    ## plot ROC curves to compare the performance of the individual classifiers
    ##
    
    # Output the plot to a PNG file for display on web.  To draw to the screen, 
    # comment this line out.
    # png(filename=figureFileName, width=1000, height=500)
    
    # load the ROCR package which draws the ROC curves
    require(ROCR)
    
    # create an ROCR prediction object from rpart() probabilities
    x.rp.prob.rocr <- prediction(x.rp.prob[,2], Y[testSet])
    # prepare an ROCR performance object for ROC curve (tpr=true positive rate, fpr=false positive rate)
    x.rp.perf <- performance(x.rp.prob.rocr, "tpr","fpr")
    performance(x.rp.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,1]<-unlist(performance(x.rp.prob.rocr,"auc")@y.values)
    # plot it
    #plot(x.rp.perf, col=2, main="ROC Curves for Machine Learning Models")
    d1 <- data.frame(x1=x.rp.perf@x.values[[1]], y1=x.rp.perf@y.values[[1]], Methods="Recursive Partitioning Trees")
    
    # Draw a legend.
    #legend(0.6, 0.6, c('rpart', 'ctree', 'cforest','bagging','svm'), 2:6)
    
    # ctree
    x.ct.prob.rocr <- prediction(x.ct.prob, Y[testSet])
    x.ct.perf <- performance(x.ct.prob.rocr, "tpr","fpr")
    performance(x.ct.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,2]<-unlist(performance(x.ct.prob.rocr,"auc")@y.values)
    # add=TRUE draws on the existing chart 
    #plot(x.ct.perf, col=3, add=TRUE)
    d2 <- data.frame(x1=x.ct.perf@x.values[[1]], y1=x.ct.perf@y.values[[1]], Methods="CITs")
    
    # cforest
    x.cf.prob.rocr <- prediction(x.cf.prob, Y[testSet])
    x.cf.perf <- performance(x.cf.prob.rocr, "tpr","fpr")
    performance(x.cf.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,3]<-unlist(performance(x.cf.prob.rocr,"auc")@y.values)
    #plot(x.cf.perf, col=4, add=TRUE)
    d3 <- data.frame(x1=x.cf.perf@x.values[[1]], y1=x.cf.perf@y.values[[1]], Methods="Random Forest with CITs")
    
    # bagging
    x.ip.prob.rocr <- prediction(x.ip.prob[,2], Y[testSet])
    x.ip.perf <- performance(x.ip.prob.rocr, "tpr","fpr")
    performance(x.ip.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,4]<-unlist(performance(x.ip.prob.rocr,"auc")@y.values)
    #plot(x.ip.perf, col=5, add=TRUE)
    d4 <- data.frame(x1=x.ip.perf@x.values[[1]], y1=x.ip.perf@y.values[[1]], Methods="Bagging")
    
    # svm
    x.svm.prob.rocr <- prediction(attr(x.svm.prob, "probabilities")[,"1"], Y[testSet])
    x.svm.perf <- performance(x.svm.prob.rocr, "tpr","fpr")
    performance(x.svm.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,5]<-unlist(performance(x.svm.prob.rocr,"auc")@y.values)
    #plot(x.svm.perf, col=6, add=TRUE)
    d5 <- data.frame(x1=x.svm.perf@x.values[[1]], y1=x.svm.perf@y.values[[1]], Methods="SVMs with Gaussian Kernel")
    
    # svm, linear
    x.svm.l.prob.rocr <- prediction(attr(x.svm.l.prob, "probabilities")[,"1"], Y[testSet])
    x.svm.l.perf <- performance(x.svm.l.prob.rocr, "tpr","fpr")
    performance(x.svm.l.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,6]<-unlist(performance(x.svm.l.prob.rocr,"auc")@y.values)
    #plot(x.svm.perf, col=7, add=TRUE)
    d6 <- data.frame(x1=x.svm.l.perf@x.values[[1]], y1=x.svm.l.perf@y.values[[1]], Methods="SVMs with Linear Kernel")
    
    # svm, polynomial
    x.svm.p.prob.rocr <- prediction(attr(x.svm.p.prob, "probabilities")[,"1"], Y[testSet])
    x.svm.p.perf <- performance(x.svm.p.prob.rocr, "tpr","fpr")
    performance(x.svm.p.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,7]<-unlist(performance(x.svm.p.prob.rocr,"auc")@y.values)
    #plot(x.svm.perf, col=8, add=TRUE)
    d7 <- data.frame(x1=x.svm.p.perf@x.values[[1]], y1=x.svm.p.perf@y.values[[1]], Methods="SVMs with Polynomial Kernel")
    
    # svm, sigmoid
    x.svm.s.prob.rocr <- prediction(attr(x.svm.s.prob, "probabilities")[,"1"], Y[testSet])
    x.svm.s.perf <- performance(x.svm.s.prob.rocr, "tpr","fpr")
    performance(x.svm.s.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,8]<-unlist(performance(x.svm.s.prob.rocr,"auc")@y.values)
    #plot(x.svm.perf, col=9, add=TRUE)
    d8 <- data.frame(x1=x.svm.s.perf@x.values[[1]], y1=x.svm.s.perf@y.values[[1]], Methods="SVMs with Sigmoid Kernel")
    
    # decision trees
    x.dt.prob.rocr <- prediction(x.dt.prob[,2], Y[testSet])
    x.dt.perf <- performance(x.dt.prob.rocr, "tpr","fpr")
    performance(x.dt.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,9]<-unlist(performance(x.dt.prob.rocr,"auc")@y.values)
    #plot(x.rf.perf, col=10, add=TRUE)
    d9 <- data.frame(x1=x.dt.perf@x.values[[1]], y1=x.dt.perf@y.values[[1]], Methods="Decision Trees")
    
    # random forest
    x.rf.prob.rocr <- prediction(x.rf.prob[,2], Y[testSet])
    x.rf.perf <- performance(x.rf.prob.rocr, "tpr","fpr")
    performance(x.rf.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,10]<-unlist(performance(x.rf.prob.rocr,"auc")@y.values)
    #plot(x.rf.perf, col=11, add=TRUE)
    d10 <- data.frame(x1=x.rf.perf@x.values[[1]], y1=x.rf.perf@y.values[[1]], Methods="Random Forest")
    
    # naive Bayes
    x.nb.prob.rocr <- prediction(x.nb.prob[,2], Y[testSet])
    x.nb.perf <- performance(x.nb.prob.rocr, "tpr","fpr")
    performance(x.nb.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,11]<-unlist(performance(x.nb.prob.rocr,"auc")@y.values)
    #plot(x.svm.perf, col=12, add=TRUE)
    d11 <- data.frame(x1=x.nb.perf@x.values[[1]], y1=x.nb.perf@y.values[[1]], Methods="NB")
    
    # naive Bayes with Laplace smoothing
    x.nb.l.prob.rocr <- prediction(x.nb.l.prob[,2], Y[testSet])
    x.nb.l.perf <- performance(x.nb.l.prob.rocr, "tpr","fpr")
    performance(x.nb.l.prob.rocr,"auc")@y.values
    AUCs[randGroup,numFeatures,12]<-unlist(performance(x.nb.l.prob.rocr,"auc")@y.values)
    #plot(x.svm.perf, col=13, add=TRUE)
    d12 <- data.frame(x1=x.nb.l.perf@x.values[[1]], y1=x.nb.l.perf@y.values[[1]], Methods="NB with Laplace Smoothing")
    
    
    #
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
    # dev.off()
  }
}


for (classifier in 1:numClassifiers){
  for (feature in 1:maxFeatures){
    AUCs[(nRandGroups+1),feature,classifier]<-mean(AUCs[1:nRandGroups,feature,classifier])
  }
}
for (classifier in 1:numClassifiers){
  write(paste("Classifier",classifier,sep=""), "proteomicsSymuAUCs321R1.csv", append=T)
  write.table(AUCs[,,classifier], "proteomicsSymuAUCs321R1.csv", append=T, quote=F, row.names=F, col.names=F, sep=",")
}
