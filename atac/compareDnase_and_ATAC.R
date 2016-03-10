require("ggplot2")
corr_atac_stam <- numeric(10)
corr_atac_crawford <- numeric(10)
corr_stam_crawford <- numeric(10)

z_score <- function(vec)
{
	return ((vec - mean(vec))/sd(vec))
}

for (i in 1:10)
{
	atac_scores <- as.numeric(scan(paste("~/scratch/atac_scores_",i,".txt",sep="")))
	dnase_stam_scores <- as.numeric(scan(paste("~/scratch/dnase_stam_scores_",i,".txt",sep="")))
	dnase_crawford_scores <- as.numeric(scan(paste("~/scratch/dnase_crawford_scores_",i,".txt",sep="")))
	
	corr_atac_stam[i] <- cor(z_score(atac_scores),z_score(dnase_stam_scores))
	corr_atac_crawford[i] <- cor(z_score(atac_scores),z_score(dnase_crawford_scores))
	corr_stam_crawford[i] <- cor(z_score(dnase_stam_scores),z_score(dnase_crawford_scores))
}
df_corr_atac_stam <- data.frame("corr"=corr_atac_stam, "types_compared"="atac_vs_dnaseStam")
df_corr_atac_crawford <- data.frame("corr"=corr_atac_crawford, "types_compared"="atac_vs_dnaseCrawford")
df_corr_stam_crawford <- data.frame("corr"=corr_stam_crawford, "types_compared"="dnaseStam_vs_dnaseCrawford")

df_all <- rbind(df_corr_atac_stam,df_corr_atac_crawford,df_corr_stam_crawford)
pdf("../results/atac_dnase_z_score_comparisons.pdf",onefile=TRUE)
ggplot(df_all, aes(x=types_compared, y=corr, fill=types_compared)) + geom_boxplot() + theme(axis.title.x = element_blank()) + ylab("correlation") + theme(axis.text.x=element_blank())
dev.off()
