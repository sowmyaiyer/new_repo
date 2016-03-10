
atac_scores_all <- as.numeric(scan(commandArgs(TRUE)[1])) #"~/scratch/ATAC_scores_all.txt"
outfile <- commandArgs(TRUE)[2]
q <- quantile(atac_scores_all, prob=seq(0,1,0.01))
quantile_values <- rev(as.numeric(q))
write(quantile_values,file=outfile, sep="\n")
