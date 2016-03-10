eigens <- as.data.frame(read.table("/home/si14w/gnearline/hdc/txt/YLJLDCHiC-hgDCE05-LPSminus-R2__hg19__genome__C-100000-iced__chrX.zScore.compartments", header=TRUE))
eigen1_v <- eigens$eigen1
eigen2_v <- eigens$eigen2
gene_density <- eigens$geneDensity
print(summary(gene_density))
print(quantile(gene_density, probs=seq(0,1,0.25)))
top25 <- quantile(gene_density, probs=seq(0,1,0.1))[["90%"]]
bottom25 <- quantile(gene_density, probs=seq(0,1,0.1))[["50%"]]
non_nan_indices <- !is.nan(eigen1_v) & !is.nan(eigen2_v)
pdf("/home/si14w/gnearline/hdc/results/compartments.pdf")
plot(eigen1_v[non_nan_indices], eigen2_v[non_nan_indices], col=ifelse(gene_density > top25, "red", ifelse(gene_density < bottom25, "blue","black")))
dev.off()
