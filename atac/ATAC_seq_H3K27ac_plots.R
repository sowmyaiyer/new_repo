library("ggplot2")
library("gridExtra")

pdf("/home/si14w/scratch/ATAC_seq_plots.pdf", width=12, height=7, onefile=TRUE)
df <- as.data.frame(read.table("~/scratch/ATAC_numbers.txt", sep=" ", header=TRUE))
fraction_intersecting_H3K27ac_peaks <- df$numbers/df$rank_cutoff
g1 <- ggplot(df,aes(x=rank_cutoff,y=fraction_intersecting_H3K27ac_peaks)) + geom_point() + xlab("Rank cutoff for ATAC score") + ylab("fraction of ATAC peaks above rank cutoff that intersect with H3K27ac peaks")
scores_int <- as.numeric(scan("~/TR/ATAC_scores_intersecting.txt"))
scores_nonint <- as.numeric(scan("~/TR/ATAC_scores_non_intersecting.txt"))

scores_nonint_df <- data.frame("ATAC_scores"=scores_nonint,"type"="not intersecting with H3K27ac")
scores_int_df <- data.frame("ATAC_scores"=scores_int,"type"="intersecting with H3K27ac")
scores_all <- rbind(scores_int_df, scores_nonint_df)
p_val <- wilcox.test(scores_int,scores_nonint)$p.value

#g2 <- ggplot(scores_int_df, aes(x=type, y=ATAC_scores, fill=type)) + geom_boxplot(data=scores_int_df)+geom_boxplot(data=scores_nonint_df) + annotate("text", label = paste("wilcox p-value =", p_val), x = 1.5, y = 3000, size = 3, colour = "black") + theme(axis.title.x = element_blank()) + ylab("ATAC seq score") + annotate("text", label = paste("n =",length(scores_int)), x = 1, y = 3000, size=3, colour="black", fontface="bold") + annotate("text", label = paste("n =",length(scores_nonint)), x = 2, y = 3000, size=3, colour="black", fontface="bold")
g2 <- ggplot(scores_all, aes(x=type, y=ATAC_scores, fill=type)) + geom_boxplot()+ annotate("text", label = paste("wilcox p-value =", p_val), x = 1.5, y = 100, size = 3, colour = "black") + theme(axis.title.x = element_blank()) + ylab("ATAC seq score")  + annotate("text", label = paste("n =",length(scores_int)), x = 1, y = 100, size=3, colour="black", fontface="bold") + annotate("text", label = paste("n =",length(scores_nonint)), x = 2, y = 100, size=3, colour="black", fontface="bold")
#+ scale_fill_discrete(name="H3K27ac peak intersection", breaks=c("intersecting", "non-intersecting"), labels=c("ATAC scores for ATAC peaks intersecting H3K27ac peaks","ATAC scores for ATAC peaks not intersecting H3K27ac peaks")) + theme(axis.title.x = element_blank()) + ylab("ATAC scores")

df <- as.data.frame(read.table("~/scratch/ATAC_tfchip_numbers.txt", sep=" ", header=TRUE))
fraction_intersecting_tfchip_peaks <- df$numbers/df$rank_cutoff
g3 <- ggplot(df,aes(x=rank_cutoff,y=fraction_intersecting_tfchip_peaks)) + geom_point()

scores_tfchip_int <- as.numeric(scan("~/TR/ATAC_scores_intersecting_tfchip.txt"))
scores_tfchip_nonint <- as.numeric(scan("~/TR/ATAC_scores_non_intersecting_tfchip.txt"))
scores_tfchip_int_df <- data.frame("ATAC_scores"=scores_tfchip_int,"type"="intersecting with TF")
scores_tfchip_nonint_df <- data.frame("ATAC_scores"=scores_tfchip_nonint,"type"="not intersecting with TF")
p_val_tfchip <- wilcox.test(scores_tfchip_int,scores_tfchip_nonint)$p.value
g4 <- ggplot(scores_tfchip_int_df, aes(x=type, y=ATAC_scores, fill=type)) + geom_boxplot(data=scores_tfchip_int_df)+geom_boxplot(data=scores_tfchip_nonint_df) + annotate("text", label = paste("wilcox p-value =", p_val_tfchip), x = 1.5, y = 100, size = 3, colour = "black") + theme(axis.title.x = element_blank()) + ylab("ATAC seq score") + annotate("text", label = paste("n =",length(scores_tfchip_int)), x = 1, y = 100, size=3, colour="black", fontface="bold") + annotate("text", label = paste("n =",length(scores_tfchip_nonint)), x = 2, y = 100, size=3, colour="black", fontface="bold")
#+ scale_fill_discrete(name="TF ChIP peak intersection",breaks=c("intersecting", "non-intersecting"),labels=c("ATAC scores for ATAC peaks intersecting TF peaks","ATAC scores for ATAC peaks not intersecting TF peaks")) +  theme(axis.title.x = element_blank()) + ylab("ATAC scores")


#df <- as.data.frame(read.table("ATAC_intersect_tfChip.bed", header=TRUE))
#g5 <- ggplot(df, aes(x=reorder(tf,ATAC_score,FUN=median),y=ATAC_score))+ geom_boxplot()+ theme(axis.title.y = element_blank())+ ylab("ATAC seq score") + coord_flip()


df <- as.data.frame(read.table("~/scratch/tf_numbers_inK27ac_and_ATAC_peaks.txt", sep=" ", header=TRUE))
df$fraction_intersecting_withK27ac_peaks <- df$num_in_K27ac/df$total
df$atac_fraction <- df$num_in_atac/df$total
g6 <- ggplot(df, aes(x=reorder(tf,fraction_intersecting_withK27ac_peaks), y=fraction_intersecting_withK27ac_peaks))+geom_bar(stat="identity")+ scale_y_continuous(limits=c(0, 1))+ theme(axis.title.y = element_blank()) + xlab("Fraction of peaks intersecting H3K27ac seq peaks") + coord_flip() 
#g7 <- ggplot(df, aes(x=reorder(tf,atac_fraction), y=atac_fraction))+geom_bar(stat="identity")+ scale_y_continuous(limits=c(0, 1)) + theme(axis.title.y = element_blank()) + xlab("Fraction of peaks intersecting ATAC seq peaks") + coord_flip()

#df <- as.data.frame(read.table("tfs_in_atac_peaks.txt", header=TRUE))
#df$fraction <- df$number/150142
#g8 <- ggplot(df, aes(x=reorder(tf,fraction), y=fraction)) + geom_bar(stat="identity") + theme(axis.title.y = element_blank()) + xlab("Fraction of ATAC peaks intersecting with TF") + coord_flip()

df <- as.data.frame(read.table("~/scratch/ATAC_numbers_by_bin_binnumbers_bin1000.txt",sep="\t", header=TRUE))
df$fraction_intersecting <- df$numbers/1000
g9 <- ggplot(df, aes(x=reorder(bin,as.numeric(bin)), y=fraction_intersecting)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7,color="black")) + xlab("ATAC seq peaks binned by rank(in 1000s)") + ylab("fraction of ATAC peaks in bin intersecting H3K27ac peaks") + scale_x_discrete(breaks=seq(1,83500,2))

#df <- as.data.frame(read.table("ATAC_scores_by_intersection_status.txt", sep="\t"))
#rank_p_val <- wilcox.test(df[which(df[,2] == "non_int"),1], df[which(df[,2] == "int"),1])$p.value
#g10 <- ggplot(df, aes(x=V2, y=rank(V1), fill=V2)) + geom_boxplot() + annotate("text", label = paste("wilcox p-value =", rank_p_val), x = 1.5, y = 150000, size = 3, colour = "black") + theme(axis.title.x = element_blank()) + ylab("ATAC seq rank") + annotate("text", label = paste("n =",length(scores_int)), x = 1, y = 150000, size=3, colour="black", fontface="bold") + annotate("text", label = paste("n =",length(scores_nonint)), x = 2, y = 150000, size=3, colour="black", fontface="bold")
df1 <- as.data.frame(read.table("~/scratch/atac_scores_percentiles.txt", sep=" ", header=TRUE))
df1$fraction <- df1$numbers/df1$total
#g11 <- ggplot(df1, aes(x=reorder(bin,-1*as.numeric(bin)), y=fraction)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=6)) + xlab("score at each percentile") + ylab("fraction of peaks in bin that intersect H3K27ac peaks")


df_int_proximal <- data.frame("ATAC_score"=as.numeric(scan("~/TR/ATAC_scores_intersecting_H3K27ac.proximal.txt")),"type"="intersecting","distance"="proximal")
df_int_distal <- data.frame("ATAC_score"=as.numeric(scan("~/TR/ATAC_scores_intersecting_H3K27ac.distal.txt")), "type"="intersecting","distance"="distal")
df_nonint_proximal <- data.frame("ATAC_score"=as.numeric(scan("~/TR/ATAC_scores_nonintersecting_H3K27ac.proximal.txt")),type="nonintersecting","distance"="proximal")
df_nonint_distal <- data.frame("ATAC_score"=as.numeric(scan("~/TR/ATAC_scores_nonintersecting_H3K27ac.distal.txt")), type="nonintersecting", "distance"="distal")
df_all <- rbind(df_int_proximal,df_int_distal,df_nonint_proximal,df_nonint_distal)
g14 <- ggplot(df_all, aes(x=distance, y=ATAC_score, fill=type)) + geom_boxplot() + annotate("text", label=paste("n =",nrow(df_int_proximal)), x=0.75, y=100, size=3) + annotate("text", label=paste("n =",nrow(df_nonint_proximal)),x=1.25, y=100, size=3)+ annotate("text", label=paste("n =",nrow(df_int_distal)),x=1.75, y=100, size=3) + annotate("text", label=paste("n =",nrow(df_nonint_distal)),x=2.25, y=100, size=3)

df2 <- as.data.frame(read.table("~/scratch/atac_scores_percentiles.distal.txt", sep=" ", header=TRUE))
df2$fraction <- df2$numbers/df2$total
g12 <- ggplot(df2, aes(x=reorder(bin,-1*as.numeric(bin)), y=fraction)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=6)) + xlab("score at each percentile") + ylab("fraction of peaks in bin that intersect H3K27ac peaks") + ggtitle("Distal ATAC peaks")

df3 <- as.data.frame(read.table("~/scratch/atac_scores_percentiles.proximal.txt", sep=" ", header=TRUE))
df3$fraction <- df3$numbers/df3$total
g13 <- ggplot(df3, aes(x=reorder(bin,-1*as.numeric(bin)), y=fraction)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=6)) + xlab("score at each percentile") + ylab("fraction of peaks in bin that intersect H3K27ac peaks") + ggtitle("Proximal ATAC peaks")
#grid.arrange(g1)
grid.arrange(g9)
#grid.arrange(g10)
#grid.arrange(g11)
grid.arrange(g2)

#grid.arrange(g14)
#grid.arrange(g12)
#grid.arrange(g13)
#grid.arrange(g3)
#grid.arrange(g4)
#grid.arrange(g5)
#grid.arrange(g6)
#grid.arrange(g7)
#grid.arrange(g8)
dev.off()
