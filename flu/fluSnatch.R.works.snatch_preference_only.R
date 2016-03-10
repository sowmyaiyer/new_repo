readCounts_file <- commandArgs(TRUE)[1]
multimapper_file <- commandArgs(TRUE)[2]

df_readCounts <- as.data.frame(read.table(readCounts_file, header=TRUE))

genes_t_df <- data.frame(row.names=df_readCounts$gene,snatchReadCounts=df_readCounts$snatchReadCounts,capReadCounts=df_readCounts$capReadCounts)
genes_t_df$snatch_preference <- genes_t_df$snatchReadCounts/sum(genes_t_df$snatchReadCounts)
genes_t_df$snatch_contribution_unnormalized <- genes_t_df$snatchReadCounts/(genes_t_df$snatchReadCounts+genes_t_df$capReadCounts)
genes_t_df$snatch_contribution_normalized <- genes_t_df$snatch_contribution_unnormalized/sum(genes_t_df$snatch_contribution_unnormalized)
print(genes_t_df)
genes_t_plus_1_df <- data.frame(row.names=rownames(genes_t_df),
				snatchReadCounts=genes_t_df$snatchReadCounts*genes_t_df$snatch_preference,
				capReadCounts=genes_t_df$capReadCounts, 
				snatch_preference=genes_t_df$snatch_preference,
				snatch_contribution_unnormalized=genes_t_df$snatch_contribution_unnormalized,
				snatch_contribution_normalized=genes_t_df$snatch_contribution_normalized)
log_l_t <- sum(genes_t_df$snatchReadCounts * (log10(genes_t_df$snatch_preference)))
log_l_t_plus_1 <- log_l_t + 0.01
multimap_readCounts <- as.data.frame(read.table(multimapper_file, header=TRUE))
pdf("/home/si14w/gnearline/flu/scripts/log_l_t_plus_1.pdf")
v_log_l_t_plus_1 <- numeric()
while ((log_l_t_plus_1 - log_l_t) > 0.00001) 
{
	for (i in 1:nrow(multimap_readCounts))
	{
		genes <- as.character(multimap_readCounts[i,"genenames"])
		multiplicity <- multimap_readCounts[i,"multiplicity"]
		genes_vector <- unlist(strsplit(genes, split=","))
		sum_snatch_preference <- sum(genes_t_df[genes_vector,"snatch_preference"])
		sum_snatch_contribution <- sum(genes_t_df[genes_vector,"snatch_contribution_normalized"])
		for (gene in genes_vector)
		{
			#genes_t_plus_1_df[gene,"snatchReadCounts"] <- genes_t_df[gene,"snatchReadCounts"] + (multiplicity * (genes_t_df[gene,"snatch_preference"]/sum_snatch_preference) *  (genes_t_df[gene,"snatch_contribution_normalized"]/sum_snatch_contribution))
			genes_t_plus_1_df[gene,"snatchReadCounts"] <- genes_t_df[gene,"snatchReadCounts"]*genes_t_df[gene,"snatch_preference"] + (multiplicity * (genes_t_df[gene,"snatch_preference"]/sum_snatch_preference)) 
		}
	}
	genes_t_plus_1_df$snatch_preference <- genes_t_plus_1_df$snatchReadCounts/sum(genes_t_plus_1_df$snatchReadCounts)
	genes_t_plus_1_df$snatch_contribution_unnormalized <- genes_t_plus_1_df$snatchReadCounts/(genes_t_plus_1_df$snatchReadCounts+genes_t_plus_1_df$capReadCounts)
	genes_t_plus_1_df$snatch_contribution_normalized <- genes_t_plus_1_df$snatch_contribution_unnormalized/sum(genes_t_plus_1_df$snatch_contribution_unnormalized)
	log_l_t <- sum(genes_t_df$snatchReadCounts * (log10(genes_t_df$snatch_preference)))
	log_l_t_plus_1 <- sum(genes_t_plus_1_df$snatchReadCounts * (log10(genes_t_plus_1_df$snatch_preference))) 
	genes_t_df <- genes_t_plus_1_df
	v_log_l_t_plus_1 <- c(v_log_l_t_plus_1, log_l_t_plus_1)
}
cat(v_log_l_t_plus_1,"\n")
plot(v_log_l_t_plus_1)
print(genes_t_plus_1_df)
dev.off()
