# omicsVisualization.R
#
# generate boxplots to visualize quantitative omics data
#
# param 1: tumor type
# param 2: data type (RNAseq, proteomics, etc)
# param 3: gene name
# param 4: clinical categories
#
# Kun-Hsing Yu
# April 3, 2016

rm(list=ls())
library(ggplot2)
setwd("~/Movies/OASISPRO/oasispro/")

sysargs<-commandArgs(trailingOnly=TRUE)
print (sysargs)
tumorType<-sysargs[1]
dataType<-sysargs[2]
geneName<-sysargs[3]
clinicalCat<-sysargs[4]


#tumorType<-"acc"
#dataType<-"proteomics"
#geneName<-"ACTN1"
#clinicalCat<-"gender"

#milestonesFileName<-paste("public/sessions/",sessionID,"/milestones.txt",sep="")


# read files
omicsFile<-read.table(paste("../data/", tumorType, "_", dataType, ".txt", sep=""), stringsAsFactors=F, sep=",")
print("Finished reading omics file")
#write("Finished reading omics file",milestonesFileName)


omicsIDFile<-read.table(paste("../data/", tumorType, "_", dataType, "_ids.txt", sep=""), stringsAsFactors=F, sep="")
omicsFile<-omicsFile[substr(omicsIDFile[,1],14,15)=="01",]
omicsIDFile<-omicsIDFile[substr(omicsIDFile[,1],14,15)=="01",]
omicsID<-substr(omicsIDFile, 1, 12)

if (clinicalCat!="NULL"){
  clinicalFile<-read.table(paste("../data/", "nationwidechildrens.org_clinical_patient_", tumorType, ".txt", sep=""), stringsAsFactors = F, sep="\t", quote="")

  print("Finished reading clinical file")
  write("Finished reading clinical file",milestonesFileName)

  clinical<-clinicalFile[4:dim(clinicalFile)[1],]
  colnames(clinical)<-clinicalFile[2,]
  clinicalID<-clinical[,2]

  IDintersection<-intersect(omicsID, clinicalID)
  omics<-omicsFile[(omicsID %in% IDintersection),]
  clinical<-clinical[(clinicalID %in% IDintersection),]
  intersectIDs<-omicsID[(omicsID %in% IDintersection)]
} else {
  omics<-omicsFile
  intersectIDs<-omicsID
}


## plot boxplots
figureFileName<-paste("public/sessions/",sessionID,"/decisionTree.png",sep="")
png(filename=figureFileName, width=1000, height=500)
if (clinicalCat!="NULL"){
  qplot(clinical[,clinicalCat],omics[,1], geom = "boxplot", xlab=geneName, ylab="Value")
} else {
  qplot(1,omics[,1], geom = "boxplot", xlab=geneName, ylab="Value")
}
dev.off()

print("Finished building boxplots")
write("Finished building boxplots",milestonesFileName,append=T)



# write histogram tables
if (FALSE){
  if (clinicalCat!="NULL"){
    h<-hist(omics[,1])
    write.table(AUCs,paste("public/sessions/",sessionID,"/Histograms.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F)
  } else {
    clinicalCatLevels<-unique(clinical[,clinicalCat])
    for (i in 1:length(clinicalCatLevels)){
      print(i)
      hist(omics[clinical[,clinicalCat]==clinicalCatLevels[i],1])
      if (i==1){
        write.table(AUCs,paste("public/sessions/",sessionID,"/Histograms.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F)
      } else {
        write.table(AUCs,paste("public/sessions/",sessionID,"/Histograms.txt",sep=""),quote=F, sep="\t",row.names=F, col.names=F, append=T)
      }
    } 
  }
}

print("Completed!")
write("Completed!",milestonesFileName,append=T)

