pdf("/home/si14w/gnearline/flu/results/multiplicity_snatches_paralogs_removed.pdf", onefile=TRUE)
for (time in c("A","B","C","D","E","F","G","H"))
{
	cat(time,"\n")
	multiplicity <- as.numeric(scan(paste("/home/si14w/gnearline/flu/txt/histogram_multireads_",time,".paralogs_removed.txt",sep=""), what="numeric"))
	multiplicity[multiplicity > 20] <- 20
	h <- hist(multiplicity, main=time, cex.main=0.75, breaks=1:20)
	legend("topright", c(paste("total snatch reads mapped (0 mismatches) =", sum(h$counts)),paste("fraction of reads mapping uniquely =", round(h$counts[1]/sum(h$counts), digits=2))), bty="n", cex=0.75)
}
dev.off()
