#library("ggplot2")
motifs <- as.data.frame(read.table("/home/si14w/gnearline/hdc/txt/motif_info.txt", sep="\t", header=TRUE))
fold_change_df <- list()
#for (i in 1:nrow(motifs))
for (i in 1:5)
{
	motif_id <- as.character(motifs[i,"motif_id"])
	cat(motif_id,"\n")
	scores_0h <- as.numeric(scan(paste("/farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_0h.",motif_id,".mean_scores.txt",sep=""),quiet=TRUE)) + 0.01
	scores_2h <- as.numeric(scan(paste("/farline/umw_flustore_schiffer/hdc/nucleoatac_out_E70_E72_2h_LPS.",motif_id,".mean_scores.txt",sep=""),quiet=TRUE)) + 0.01
	fold_change_df[[motif_id]] <- log2(scores_2h/scores_0h)
}
pdf("/home/si14w/gnearline/hdc/results/test.pdf")
boxplot(fold_change_df,range=0)
dev.off()
