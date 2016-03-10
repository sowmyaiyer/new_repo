f_high <- commandArgs(TRUE)[1]
f_med <- commandArgs(TRUE)[2]
f_low <- commandArgs(TRUE)[3]
f_verylow <- commandArgs(TRUE)[4]
outfile <- commandArgs(TRUE)[5]
title <- commandArgs(TRUE)[6]
df_high <- as.data.frame(read.table(f_high, row.names=1, sep="\t"))
cat("read\n")
df_med <- as.data.frame(read.table(f_med, row.names=1, sep="\t"))
cat("read\n")
df_low <- as.data.frame(read.table(f_low, row.names=1, sep="\t"))
cat("read\n")
df_verylow <- as.data.frame(read.table(f_verylow, row.names=1, sep="\t"))
cat("read\n")
cmeans_high <- colMeans(df_high)
cmeans_med <- colMeans(df_med)
cmeans_low <- colMeans(df_low)
cmeans_verylow <- colMeans(df_verylow)
pdf(outfile)
plot(seq(-995,995,10), cmeans_high, type="l", col="magenta", ylim=c(0,max(c(cmeans_high,cmeans_med,cmeans_low,cmeans_verylow))), xlab="distance to TSS", ylab="mean signal", bty="n", main=title, cex.main=0.65)
lines(seq(-995,995,10),cmeans_med, col="red")
lines(seq(-995,995,10),cmeans_low, col="orange")
lines(seq(-995,995,10),cmeans_verylow, col="yellow")
legend("topright",c(paste("high(n =",nrow(df_high),")"), paste("med(n =",nrow(df_med),")"), paste("low(n =",nrow(df_low),")"),paste("verylow(n =",nrow(df_verylow),")")), col=c("magenta","red","orange","yellow"), lty="solid", bty="n",cex=0.6)
dev.off()
