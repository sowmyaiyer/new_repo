#times <- c("6h","12h","24h","48h","mock12h")
times <- c("24h","mock12h")
combos <- combn(times,2)
print(combos)
for (i in 1:ncol(combos))
{
	cat(combos[1,i],combos[2,i],"\n")
	time1_abundance <- scan(paste("/home/si14w/gnearline/flu/txt/04SEP15.A549_B59_totalRNA_",combos[1,i],".transcript_abundances.txt",sep=""),sep="\n", what="character")
	time1_abundance_list <- strsplit(time1_abundance, "\t")
	names(time1_abundance_list) <- lapply(time1_abundance_list, "[[", index=1)
	
	time2_abundance <- scan(paste("/home/si14w/gnearline/flu/txt/04SEP15.A549_B59_totalRNA_",combos[2,i],".transcript_abundances.txt",sep=""),sep="\n",  what="character")
	time2_abundance_list <- strsplit(time2_abundance, "\t")
	names(time2_abundance_list) <- lapply(time1_abundance_list, "[[", index=1)
	df <- data.frame(names(time1_abundance_list),0.0)
	for (gene in names(time1_abundance_list))
	{
		time1_abundance_values <- as.numeric(time1_abundance_list[[gene]][2:length(time1_abundance_list[[gene]])])
		time2_abundance_values <- as.numeric(time2_abundance_list[[gene]][2:length(time2_abundance_list[[gene]])])
		if ((mean(time1_abundance_values) > 100 & mean(time2_abundance_values) > 100 ) & length(time1_abundance_values) > 1 & length(time1_abundance_values) < 5 & length(time1_abundance_values) == length(time2_abundance_values))
		{
			if (length(time1_abundance_values) == 2 & time1_abundance_values[1] == time1_abundance_values[2])
			{
				time1_abundance_values[2] <- time1_abundance_values[1] + 0.01
			}	
			if (length(time2_abundance_values) == 2 & time2_abundance_values[1] == time2_abundance_values[2]) 
			{
				time2_abundance_values[2] <- time2_abundance_values[1] + 0.01
			}
			cat(gene,"\t")
			p <- chisq.test(time1_abundance_values, time2_abundance_values, rescale.p=TRUE)$p.value
			df[which(df[,1] == gene),2] <- p
			cat(p,"\n")
		}
	}
	write.table(df, file=paste("/home/si14w/gnearline/flu/txt/transcript_dist_differences_",combos[1,i],"_",combos[2,i],".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
