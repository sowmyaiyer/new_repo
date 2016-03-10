require("ggplot2")
require("gridExtra")
#atac_distal_fg <- as.data.frame(read.table("../data/H1_ATAC.distal.bed"))
#atac_distal_fg <- as.data.frame(read.table("../txt/H7Hesc_DNase_Stam_signal_in_ATAC_peaks.bed"))
atac_distal_fg <- as.data.frame(read.table("../txt/H7Hesc_ATAC_signal_in_ATAC_peaks.bed"))
atac_distal_bg <- as.data.frame(read.table("../txt/H7Hesc_ATAC_signal_in_random_regions.bed"))


atac_peak_lengths <- as.numeric(atac_distal_fg[,3]) - as.numeric(atac_distal_fg[,2])
atac_scores_fg_scores <- as.numeric(atac_distal_fg[,5])
atac_scores_bg_scores <- as.numeric(atac_distal_bg[,5])
atac_distal_fg$type="atac_peaks"
atac_distal_bg$type="bg"

atac_all <- rbind(atac_distal_fg,atac_distal_bg)

pdf("../results/atac_scores_peak_lengths.pdf", onefile=TRUE)
atac_peak_lengths[atac_peak_lengths > 1000] <- 1000
hist(atac_peak_lengths, breaks=seq(1,1001,25), xlab="length of ATAC peaks")
atac_all[atac_all[,"V5"] > 50,"V5"] <- 50
#atac_score_plot <- ggplot(atac_all, aes(x=V5)) + geom_histogram(data=subset(atac_all,type == "atac_peaks"), fill="red",binwidth = 1, alpha=0.2) + geom_histogram(data=subset(atac_all,type =="bg"), fill="blue",binwidth = 1, alpha=0.2) + xlab("ATAC scores")
atac_score_plot2 <- ggplot(atac_all, aes(x=V5, fill=type)) + geom_histogram(position="identity",binwidth = 1, alpha=0.5) + xlab("ATAC scores")
#grid.arrange(atac_score_plot)
grid.arrange(atac_score_plot2)
dev.off()
