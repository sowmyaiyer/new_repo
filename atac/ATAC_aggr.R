require("ggplot2")
sites_in_peaks <- as.data.frame(read.table(commandArgs(TRUE)[1],stringsAsFactors=FALSE))
sites_notin_peaks <- as.data.frame(read.table(commandArgs(TRUE)[2],stringsAsFactors=FALSE))
cat(max(sites_in_peaks[,2]),"\n")
cat(max(sites_notin_peaks[,2]),"\n")
sites_in_peaks$position <- "in_peaks"
sites_notin_peaks$position <- "not_in_peaks"
sites <- rbind(sites_in_peaks,sites_notin_peaks)
colnames(sites) <- c("window","mean_score","position")
pdf("../results/atac_aggr.pdf", onefile=TRUE)
plot(seq(-107,107,1),as.numeric(sites_in_peaks[,2]), type="l", col="blue")
lines(seq(-107,107,1),as.numeric(sites_notin_peaks[,2]), col="red")
dev.off()
