input <- commandArgs(TRUE)[1]
output <- commandArgs(TRUE)[2]
reads_df <- as.data.frame(read.table(input))
colnames(reads_df) <- c("chr","start","end","gene_name","score","strand","readCount")

scalingFactor <- sum(as.numeric(reads_df$readCount))/1000000
reads_df$rpm <- reads_df$readCount/scalingFactor
write.table(reads_df, file=output, quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")
cat(" check", sum(reads_df$rpm),"\n")
