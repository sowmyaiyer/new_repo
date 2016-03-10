library("ggplot2")
library("gridExtra")

pdf("/home/si14w/gnearline/ATAC_analysis/results/DNase_seq_plots.pdf", width=12, height=7, onefile=TRUE)

df1 <- as.data.frame(read.table("../txt/dnase_scores_percentiles.txt", sep=" ", header=TRUE))
df1$fraction <- df1$numbers/df1$total
g10 <- ggplot(df1, aes(x=reorder(bin,-1*as.numeric(bin)), y=fraction)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=6)) + xlab("score at each percentile") + ylab("fraction of peaks in bin that intersect H3K27ac peaks") + ggtitle("DNase peaks")
grid.arrange(g10)

#df2 <- as.data.frame(read.table("~/scratch/dnase_scores_percentiles.distal.txt", sep=" ", header=TRUE))
#df2$fraction <- df2$numbers/df2$total
#g12 <- ggplot(df2, aes(x=reorder(bin,-1*as.numeric(bin)), y=fraction)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=6)) + xlab("score at each percentile") + ylab("fraction of peaks in bin that intersect H3K27ac peaks") + ggtitle("Distal ATAC peaks")

#df3 <- as.data.frame(read.table("~/scratch/dnase_scores_percentiles.proximal.txt", sep=" ", header=TRUE))
#df3$fraction <- df3$numbers/df3$total
#g13 <- ggplot(df3, aes(x=reorder(bin,-1*as.numeric(bin)), y=fraction)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=6)) + xlab("score at each percentile") + ylab("fraction of peaks in bin that intersect H3K27ac peaks") + ggtitle("Proximal ATAC peaks")
#grid.arrange(g12)
#grid.arrange(g13)
dev.off()
