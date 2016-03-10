tads <- as.data.frame(read.table(commandArgs(TRUE)[1])
tad_lengths <- as.numeric(tads[,3]) -as.numeric(tads[,2])
tad_lengths[tad_lengths > 3000000] <- 3000000
pdf("../results/tad_length_hist.pdf")
hist(tad_lengths, breaks=seq(0,3000000,40000))
dev.off()
