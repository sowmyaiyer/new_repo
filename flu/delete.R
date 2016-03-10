print("here")

{
        print("here")
        for (i in 1:10)
        {
                genes <- "A,B"
                multiplicity <- 100
                genes_vector <- unlist(strsplit(genes, split=","))
                for (gene in genes_vector)
                {
                        genes_t_plus_1_df[gene,"snatchReadCounts"] <- genes_t_plus_1_df[gene,"snatchReadCounts"] #+                                                                      multiplicity * genes_t_plus_1_df[gene,"snatch_preference"] *  genes_t_plus_1_df[gene,"snatch_contribution_normalized]
                }
        }
#        genes_t_plus_1_df$snatch_preference <- genes_t_plus_1_df$snatchReadCounts/sum(genes_t_plus_1_df$snatchReadCounts)
#        genes_t_plus_1_df$snatch_contribution_unnormalized <- genes_t_plus_1_df$snatchReadCounts/(genes_t_plus_1_df$snatchReadCounts+genes_t_plus_1_df$capReadCounts)
#        genes_t_plus_1_df$snatch_contribution_normalized <- genes_t_plus_1_df$snatch_contribution_unnormalized/sum(genes_t_plus_1_df$snatch_contribution_unnormalized)
#
#        log_l_t <- sum(genes_t_df$snatchReadCounts * (log10(genes_t_df$snatch_preference) + log10(genes_t_df$snatch_contribution_normalized)))
#        log_l_t_plus_1 <- sum(genes_t_plus_1_df$snatchReadCounts * (log10(genes_t_plus_1_df$snatch_preference) + log10(genes_t_plus_1_df$snatch_contribution_normalized)))
}
