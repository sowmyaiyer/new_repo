library("ggplot2")
library(gridExtra)
setwd("~/gnearline/flu/scripts/")
pdf("/home/si14w/gnearline/flu/results/tss_scores.pdf", onefile=TRUE)
plots <- list()
for (time in c("A","B","C","D","E","F","G","H"))
{
	exons=data.frame("score"=as.numeric(scan(paste("../txt/tss_tpm_exons_",time,".txt",sep=""))), "label"="exons")
	introns <- data.frame("score"=as.numeric(scan(paste("../txt/tss_tpm_introns_",time,".txt",sep=""))), "label"="intron")
	intergenic <- data.frame("score"=as.numeric(scan(paste("../txt/tss_tpm_intergenic_",time,".txt",sep=""))), "label"="intergenic")
	all_scores <- rbind(exons,introns,intergenic)
	print(ggplot(all_scores, aes(x=label,y=score,fill=label)) + geom_boxplot() +  coord_cartesian(ylim = quantile(all_scores$score, c(0.1, 0.9))) + ggtitle(time) + xlab(""))
}
dev.off()
