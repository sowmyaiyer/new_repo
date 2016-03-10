inputbam <- commandArgs(TRUE)[1]
sampleName <- commandArgs(TRUE)[2]
myCAGEset <- new("CAGEset", genomeName="BSgenome.Hsapiens.UCSC.hg19", inputFiles=inputbam, inputFilesType="bam",sampleLabels=c(sampleName))
getCTSS(myCAGEset)
ctss <- CTSStagCount(myCAGEset)
write.table(ctss, file=paste("/home/si14w/gnearline/flu/txt/",sampleName,".txt",sep=""), quote=FALSE, row.names=FALSE)
