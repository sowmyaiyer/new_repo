require("ggplot2")
require("reshape2")
output_dir_name <- commandArgs(TRUE)[1]
tf_intersection_data <- as.data.frame(read.table(paste("../results/",output_dir_name,"/H1hesc_TF_intersection_stats.txt",sep=""), header=TRUE))
tf_intersection_data$fraction_intersecting_atac_enhancers <- (tf_intersection_data$peaks_intersecting_atac_enhancers - tf_intersection_data$peaks_intersecting_atac_both_enh_and_nonenh)/tf_intersection_data$total_peaks
tf_intersection_data$fraction_intersecting_atac_nonenhancers <- (tf_intersection_data$peaks_intersecting_atac_nonenhancers - tf_intersection_data$peaks_intersecting_atac_both_enh_and_nonenh)/tf_intersection_data$total_peaks
tf_intersection_data$fraction_intersecting_atac_enh_and_nonenh <- (tf_intersection_data$peaks_intersecting_atac_both_enh_and_nonenh)/tf_intersection_data$total_peaks
tf_intersection_data$fraction_not_intersecting_atac <- tf_intersection_data$peaks_not_intersecting_atac/tf_intersection_data$total_peaks

tf_frac <- data.frame(paste(tf_intersection_data$tfname,"(",tf_intersection_data$total_peaks,")",sep=""),tf_intersection_data$fraction_intersecting_atac_nonenhancers,tf_intersection_data$fraction_intersecting_atac_enhancers,tf_intersection_data$fraction_intersecting_atac_enh_and_nonenh,tf_intersection_data$fraction_not_intersecting_atac)
print(rowSums(tf_frac[,c(2,3,4,5)]))
colnames(tf_frac) <- c("TF(num_peaks)","in_ATAC_nonenhancers","in_ATAC_enhancers","in_ATAC_enh_and_non-enh","not_in_ATAC")
tf_frac$TF <- reorder(tf_frac$TF, tf_frac$in_ATAC_nonenhancers)

tf_melted <- melt(tf_frac)
print(tf_melted)

enhancer_intersection_data <- as.data.frame(read.table(paste("../results/",output_dir_name,"/enhancer_intersecting_tf.txt",sep=""), header=TRUE))
nonenhancer_intersection_data <- as.data.frame(read.table(paste("../results/",output_dir_name,"/nonenhancer_intersecting_tf.txt",sep=""), header=TRUE))
enhancer_intersection_data$fraction <- (enhancer_intersection_data$enhancers_bound_to_tf/enhancer_intersection_data$total_enhancers)*10000/enhancer_intersection_data$totalTFPeaks
nonenhancer_intersection_data$fraction <- (nonenhancer_intersection_data$nonenhancers_bound_to_tf/nonenhancer_intersection_data$total_nonenhancers)*10000/nonenhancer_intersection_data$totalTFPeaks


tf_atac_int_score <- as.data.frame(read.table(paste("../results/",output_dir_name,"/TFs_and_ATAC_scores.txt",sep="")))
colnames(tf_atac_int_score) <- c("TF","ATAC_score")


tf_scores_in_enhancer_atac <- as.data.frame(read.table(paste("../results/",output_dir_name,"/H1_TFPeaks_in_enhancer_atac_peaks.txt",sep="")))
tf_scores_in_nonenhancer_atac <- as.data.frame(read.table(paste("../results/",output_dir_name,"/H1_TFPeaks_in_nonenhancer_atac_peaks.txt",sep="")))
tf_scores_in_non_atac <- as.data.frame(read.table(paste("../results/",output_dir_name,"/H1_TFPeaks_not_in_atac_peaks.txt",sep="")))

tf_scores_by_category <- rbind(tf_scores_in_enhancer_atac,tf_scores_in_nonenhancer_atac,tf_scores_in_non_atac)
colnames(tf_scores_by_category) <- c("TF", "score", "type")

pdf(paste("../results/",output_dir_name,"/H1hesc_tf_intersect_atac.pdf",sep=""), width=15, height=15)
m <- t(as.matrix(tf_frac[,c(2,3,4,5)]))
colnames(m) <- tf_frac[,1]
bp <- barplot(m, horiz=TRUE, axes=FALSE)
axis(1, at=seq(0,1,0.05), seq(0,1,0.05), las=2, cex=0.3)
axis(2, at=bp, colnames(m),las=1)
ggplot(tf_melted, aes(x=TF,y=value,fill=variable)) + geom_bar(position="fill", stat="identity") + geom_hline(yintercept=seq(0,1,0.05), linetype = "longdash") +  geom_hline(yintercept=c(0.25,0.5,0.75)) + ylab("fraction of ChIP-seq peaks") + coord_flip() + xlab("TF(number of peaks)") 
ggplot(enhancer_intersection_data, aes(x=reorder(TF,fraction), y=fraction)) + geom_bar(stat="identity") + xlab("TF") + ylab("Fraction of ATAC enhancers bound by TF") + coord_flip()
ggplot(nonenhancer_intersection_data, aes(x=reorder(TF,fraction), y=fraction)) + geom_bar(stat="identity") + xlab("TF") + ylab("Fraction of ATAC non-enhancers bound by TF") + coord_flip()
ggplot(tf_atac_int_score, aes(x=reorder(TF,ATAC_score,FUN=median), y=ATAC_score)) + geom_boxplot() + coord_flip()
dev.off()
ggplot(tf_scores_by_category, aes(x=TF, y=score, fill=type)) + geom_boxplot() + coord_cartesian(ylim=c(0,750))
ggsave(file=paste("../results/",output_dir_name,"/H1hesc_tf_intersect_atac_peak_scores_by_category.svg",sep=""),width=75, height=10, limitsize=FALSE)
