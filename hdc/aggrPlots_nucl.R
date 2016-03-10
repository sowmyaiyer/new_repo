f_0h <- commandArgs(TRUE)[1]
f_30min <- commandArgs(TRUE)[2]
f_2h <- commandArgs(TRUE)[3]
f_4h <- commandArgs(TRUE)[4]
f_24h <- commandArgs(TRUE)[5]
outfile <- commandArgs(TRUE)[6]
title <- commandArgs(TRUE)[7]
motif_length <- as.numeric(commandArgs(TRUE)[7])
df_0h <- as.data.frame(read.table(f_0h, row.names=1, sep="\t"))
cat("read\n")
df_30min <- as.data.frame(read.table(f_30min, row.names=1, sep="\t"))
cat("read\n")
df_2h <- as.data.frame(read.table(f_2h, row.names=1, sep="\t"))
cat("read\n")
df_4h <- as.data.frame(read.table(f_4h, row.names=1, sep="\t"))
cat("read\n")
cmeans_0h <- colMeans(df_0h)
cmeans_30min <- colMeans(df_30min)
cmeans_2h <- colMeans(df_2h)
cmeans_4h <- colMeans(df_4h)
pdf(outfile)
left <- -(200+motif_length-1)/2
right <- (200+motif_length-1)/2
cat(length(seq(left, right,1)),"\n")
plot(seq(left, right,1), cmeans_0h, type="l", col="magenta", ylim=c(0,max(c(cmeans_0h,cmeans_30min,cmeans_2h,cmeans_4h))), xlab="distance from motif", ylab="mean signal", bty="n", main=title, cex.main=0.65)
lines(seq(left, right,1),cmeans_30min, col="red")
lines(seq(left, right,1),cmeans_2h, col="orange")
lines(seq(left, right,1),cmeans_4h, col="yellow")
lines(seq(left, right,1),cmeans_24h, col="brown")
legend("topright",c(paste("0h(n =",nrow(df_0h),")"), paste("30min(n =",nrow(df_30min),")"), paste("2h(n =",nrow(df_2h),")"),paste("4h(n =",nrow(df_4h),")"), paste("24h(n =",nrow(df_24h),")")), col=c("magenta","red","orange","yellow","brown"), lty="solid", bty="n",cex=0.6)
dev.off()
