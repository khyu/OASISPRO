# omicsVisualization.R
#
# generate boxplots to visualize quantitative omics data
#
# param 1: tumor type
# param 2: data type (RNAseq, proteomics, etc)
# param 3: gene name
# param 4: clinical categories
# param 5: session ID
#
# Kun-Hsing Yu
# April 3, 2016

rm(list=ls())
library(ggplot2)
#setwd("~/Movies/OASISPRO/oasispro/")

sysargs<-commandArgs(trailingOnly=TRUE)
print (sysargs)
tumorType<-sysargs[1]
dataType<-sysargs[2]
geneName<-sysargs[3]
clinicalCat<-sysargs[4]
sessionID<-sysargs[5]
argsFileName<-paste("public/sessions/",sessionID,"/args.txt",sep="")
write(sysargs,argsFileName)

nContinuousSplit<-5
#tumorType<-"acc"
#dataType<-"proteomics"
#geneName<-"14-3-3_beta-R-V"
#geneName<-"TFRC-R-V"
#clinicalCat<-"gender"
#clinicalCat<-"pathologic_T"
#sessionID<-"810110794241357"
milestonesFileName<-paste("public/sessions/",sessionID,"/milestones.txt",sep="")


# read files
#omicsFile<-read.table(paste("../data/", tumorType, "_", dataType, ".txt", sep=""), stringsAsFactors=F, sep=",")
load(paste("../data/", tumorType, "_", dataType, ".RData", sep=""))
omicsName<-read.table(paste("../data/", tumorType, "_", dataType, "_elemid.txt", sep=""), stringsAsFactors=F, sep=",")
print("Finished reading omics file")
write("Finished reading omics file",milestonesFileName)


omicsIDFile<-read.table(paste("../data/", tumorType, "_", dataType, "_ids.txt", sep=""), stringsAsFactors=F, sep="")
if (tumorType=="laml"){
  tumorPartID="03"
} else {
  tumorPartID="01"
}
omicsFile<-omicsFile[substr(omicsIDFile[,1],14,15)==tumorPartID,]
omicsIDFile<-omicsIDFile[substr(omicsIDFile[,1],14,15)==tumorPartID,]
omicsID<-substr(omicsIDFile, 1, 12)

if (clinicalCat!="NULL"){
  clinicalFile<-read.table(paste("../data/", "nationwidechildrens.org_clinical_patient_", tumorType, ".txt", sep=""), stringsAsFactors = F, sep="\t", quote="")

  clinicalTitlesUniqueAnnotated<-read.table(paste("../data/clinicalTitlesUniqueAnnotated.txt", sep=""), stringsAsFactors = F, sep="\t", quote="")
  clinicalTitlesUniqueAnnotated<-t(clinicalTitlesUniqueAnnotated)
  clinicalDataType<-clinicalTitlesUniqueAnnotated[which(clinicalTitlesUniqueAnnotated[,1]==clinicalCat),2]
  
  print("Finished reading clinical file")
  write("Finished reading clinical file",milestonesFileName)

  clinical<-clinicalFile[4:dim(clinicalFile)[1],]
  colnames(clinical)<-clinicalFile[2,]
  clinicalID<-clinical[,"bcr_patient_barcode"]

  IDintersection<-intersect(omicsID, clinicalID)
  omics<-omicsFile[(omicsID %in% IDintersection),]
  clinical<-clinical[(clinicalID %in% IDintersection),]
  intersectIDs<-omicsID[(omicsID %in% IDintersection)]
  
} else {
  omics<-omicsFile
  intersectIDs<-omicsID
}

omicsColNum<-which(omicsName[1,]==geneName)[1]
## plot boxplots
figureFileName<-paste("public/sessions/",sessionID,"/boxplot.png",sep="")
png(filename=figureFileName, width=1000, height=500)
if (clinicalCat!="NULL"){
  if (clinicalDataType=="c") {
    clinicalNArow<-substr(clinical[,clinicalCat],1,1)=="["
    if (length(clinicalNArow)>0){
      clinical[clinicalNArow,clinicalCat]<-NA
    }
    clinical[,clinicalCat]<-as.numeric(clinical[,clinicalCat])
    clinicalSplitGroup<-cut(clinical[,clinicalCat], nContinuousSplit)
    qplot(clinicalSplitGroup,omics[,omicsColNum], geom = "boxplot", xlab=geneName, ylab="Value")
  }
  else { # clinicalDataType=="d"
    qplot(clinical[,clinicalCat],omics[,omicsColNum], geom = "boxplot", xlab=geneName, ylab="Value")
  }
} else { # clinicalCat=="NULL"
  qplot(1,omics[,1], geom = "boxplot", xlab=geneName, ylab="Value")
}
dev.off()

print("Finished building boxplots")
write("Finished building boxplots",milestonesFileName,append=T)


## perform statistical tests
if (clinicalCat!="NULL"){
  if (length(unique(clinical[,clinicalCat]))==1){ # No test
    pValue<-"Not enough observations"
  } else if (length(unique(clinical[,clinicalCat]))==2){ # Wilcoxon test
    uniqueClinicalCats<-unique(clinical[,clinicalCat])
    pValue<-tryCatch({
      wilcox.test(omics[clinical[,clinicalCat]==uniqueClinicalCats[1],omicsColNum],omics[clinical[,clinicalCat]==uniqueClinicalCats[2],omicsColNum])$p.value
    }, error = function(err) {
      pValue<-"Not enough observations"
    })
  } else { # ANOVA
    uniqueClinicalCats<-unique(clinical[,clinicalCat])
    fit <- aov(omics[,omicsColNum] ~ as.factor(clinical[,clinicalCat]))
    pValue<-summary(fit)[[1]][["Pr(>F)"]][1]
  }
} else {
  pValue<-""
}
print(pValue)
write(pValue,paste("public/sessions/",sessionID,"/pValue.txt",sep=""))


# write histogram tables
if (FALSE){
  if (clinicalCat!="NULL"){
    h<-hist(omics[,1])
    write.table(h$breaks,paste("public/sessions/",sessionID,"/HistogramsBreaks.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F)
    write.table(h$counts,paste("public/sessions/",sessionID,"/HistogramsCounts.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F)
  } else {
    clinicalCatLevels<-unique(clinical[,clinicalCat])
    for (i in 1:length(clinicalCatLevels)){
      print(i)
      hist(omics[clinical[,clinicalCat]==clinicalCatLevels[i],1])
      if (i==1){
        write.table(h$breaks,paste("public/sessions/",sessionID,"/HistogramsBreaks.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F)
        write.table(h$counts,paste("public/sessions/",sessionID,"/HistogramsCounts.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F)
      } else { # append
        write.table(h$breaks,paste("public/sessions/",sessionID,"/HistogramsBreaks.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F, append=T)
        write.table(h$counts,paste("public/sessions/",sessionID,"/HistogramsCounts.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F, append=T)
      }
    } 
  }
}

print("Completed!")
write("Completed!",milestonesFileName,append=T)

