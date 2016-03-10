require("ggplot2")
tf_table <- as.data.frame(read.table("../txt/reads_intersect_tf.txt", header=TRUE))
tf_table[which(tf_table[,2] > 1000),2] <- 1000
means <- aggregate(fragment_size~tf, tf_table, mean)
cat("read data\n")
pdf("../results/reads_intersect_tf.pdf", onefile=TRUE)
ggplot(tf_table, aes(x=reorder(tf,fragment_size,FUN=median), y=fragment_size)) + geom_boxplot() + stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show_guide = FALSE) +  xlab("TF") + coord_flip()
dev.off()
