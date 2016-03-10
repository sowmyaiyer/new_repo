#motifs <- as.data.frame(read.table("/home/si14w/gnearline/hdc/txt/motif_info.txt", sep="\t", header=TRUE))
motif_info_file <- commandArgs(TRUE)[1]
motif_info_file_name <- basename(motif_info_file)
motifs <- as.data.frame(read.table(motif_info_file, sep="\t"))
pdf(paste("/home/si14w/gnearline/hdc/results/merged_E70_E72_nucl_signal.",motif_info_file_name,".pdf",sep=""))
par(mfrow=c(3,3))
for (i in 1:nrow(motifs))
{
	cat(i,"\n")
	motif <- as.character(motifs[i,1])
	motif_name <- as.character(motifs[i,2])
	cat(motif_name,"\n")
	frac_in_expressed_promoters <- round(as.numeric(motifs[i,3])/as.numeric(motifs[i,4]), digits=2)
	motif_length <- as.numeric(motifs[i,5])
	tryCatch({
		df_0h <- as.data.frame(read.table(paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72.0h.",motif,".all_signals.txt",sep=""), row.names=1))
		df_30min <- as.data.frame(read.table(paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72.30min_LPS.",motif,".all_signals.txt",sep=""), row.names=1))  
		df_2h <- as.data.frame(read.table(paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72.2h_LPS.",motif,".all_signals.txt",sep=""), row.names=1))  
		df_4h <- as.data.frame(read.table(paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72.4h_LPS.",motif,".all_signals.txt",sep=""), row.names=1)) 
		df_24h <- as.data.frame(read.table(paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72.24h_LPS.",motif,".all_signals.txt",sep=""), row.names=1)) 
		filter_indices <- which(rowSums(df_0h > 0) & rowSums(df_30min) > 0 & rowSums(df_2h) > 0 & rowSums(df_4h) > 0 & rowSums(df_24h) > 0)
#		df_0h <- df_0h[which(rowSums(df_2h) > 1),]
#		df_30min <- df_30min[which(rowSums(df_2h) > 1),]
#		df_2h <- df_2h[which(rowSums(df_2h) > 1),]
#		df_4h <- df_4h[which(rowSums(df_2h) > 1),]
#		df_24h <- df_24h[which(rowSums(df_2h) > 1),]
		cmeans_0h <- colMeans(df_0h[filter_indices,])
		cmeans_30min <- colMeans(df_30min[filter_indices,])
		cmeans_2h <- colMeans(df_2h[filter_indices,])
		cmeans_4h <- colMeans(df_4h[filter_indices,])
		cmeans_24h <- colMeans(df_24h[filter_indices,])
		ymax <- max(c(cmeans_0h,cmeans_30min,cmeans_2h,cmeans_4h,cmeans_24h))
		write(c(motif,cmeans_0h),file=paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72_nucl_signal.colmeans.",motif,".0h.txt",sep=""),sep="\t", ncolumns=length(cmeans_0h)+1)
		write(c(cmeans_30min),file=paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72_nucl_signal.colmeans.",motif,".30min_LPS.txt",sep=""),sep="\t", ncolumns=length(cmeans_30min)+1)
		write(c(cmeans_2h),file=paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72_nucl_signal.colmeans.",motif,".2h_LPS.txt",sep=""),sep="\t", ncolumns=length(cmeans_2h)+1)
		write(c(cmeans_4h),file=paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72_nucl_signal.colmeans.",motif,".4h_LPS.txt",sep=""),sep="\t", ncolumns=length(cmeans_4h)+1)
		write(c(cmeans_24h),file=paste("/farline/umw_flustore_schiffer/hdc/merged_E70_E72_nucl_signal.colmeans.",motif,".24h_LPS.txt",sep=""),sep="\t", ncolumns=length(cmeans_24h)+1)
	#	left <- -(200+motif_length-1)/2
	#	right <- (200+motif_length-1)/2
		plot(cmeans_0h, type="l", col="magenta", ylim=c(0,max(c(cmeans_0h,cmeans_30min,cmeans_2h,cmeans_4h))), xlab="", ylab="", bty="n", main=substr(motif_name,1,20), cex.main=0.65, xaxt="n", bty="n", cex.lab=0.3, cex.axis=0.5)
		lines(cmeans_30min, col="red")
		lines(cmeans_2h, col="orange")
		lines(cmeans_4h, col="yellow")
		lines(cmeans_24h, col="black")
		legend("bottom",c("0h","30min","2h","4h","24h"),col=c("magenta","red","orange","yellow","black"), lty="solid", bty="n",cex=0.6)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
dev.off()
