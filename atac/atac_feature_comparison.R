require("ggplot2")
require("gridExtra")
output_dir_name <- commandArgs(TRUE)[1]
cat(paste("../results/",output_dir_name,"/H1_ATAC_trainingset_withhismods.txt",sep=""),"\n")
atac_features <- as.data.frame(read.table(paste("../results/",output_dir_name,"/H1_ATAC_trainingset_withhismods.txt",sep=""), sep="\t",header=TRUE))
atac_features[which(atac_features$label == 0),"label"] <- "non-intersecting"
atac_features[which(atac_features$label == 1),"label"] <- "intersecting"

atac_features$label <- as.character(atac_features$label)
length_intersecting <- length(which(atac_features$label == "intersecting"))
length_nonintersecting <- length(which(atac_features$label == "non-intersecting"))
cat(length_intersecting,length_nonintersecting,"\n")
cat("data loaded\n")
pdf(paste("../results/",output_dir_name,"/H1_feature_comparison_int_vs_nonint.with_hismods.pdf", sep=""),onefile=TRUE)
features <- c("atac_score","phylop_cons","dist_to_gene","gc","H3K4me1","H3K4me3","H3K4me2","H3K9ac","H3K9me3","H3K27me3","H3K27ac")
for (feature in features)
{
	cat(feature,"\n")
	p_val <- wilcox.test(atac_features[,as.character(feature)],atac_features[,as.character(feature)])$p.value
	gg_atac <- ggplot(atac_features, aes_string(x="label", y=feature, fill="label")) + geom_boxplot(outlier.shape=NA) + labs(title=paste(feature,"int =",length_intersecting,"non-int =",length_nonintersecting),x="",y=feature) 
	#gg_atac_rescaled = gg_atac + coord_cartesian(ylim = quantile(atac_features[,as.character(feature)], c(0.1, 0.9))) + scale_y_continuous(limits = quantile(atac_features[,as.character(feature)], c(0.1, 0.9)))
	gg_atac_rescaled = gg_atac + scale_y_continuous(limits = quantile(atac_features[,as.character(feature)], c(0.1, 0.9))) 
	#+ coord_cartesian(ylim = quantile(atac_features[,as.character(feature)], c(0.1, 0.9)))
	grid.arrange(gg_atac_rescaled)
}
#features_to_compare <- c("phylop_cons","dist_to_gene","gc","H3K4me1","H3K4me3","H3K4me2","H3K9ac","H3K9me3","H3K27me3","H3K27ac")
#for (f in features_to_compare)
#{
#	cat(f,"\n")
#	gg <- ggplot(atac_features, aes_string(x="atac_score",y=f)) + geom_point(aes(color=factor(label)))
#	grid.arrange(gg)
#}

dev.off()
