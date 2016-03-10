cov_input <- as.numeric(scan("/home/si14w/gnearline/hdc/txt/E37_input_coverage_across_genome.txt", what="numeric"))
cov_input[cov_input > 0.1] <- 0.1
cov_IgG <- as.numeric(scan("/home/si14w/gnearline/hdc/txt/E37_IgG_coverage_across_genome.txt", what="numeric"))
cov_IgG[cov_IgG > 0.1] <- 0.1
cov_H3K4me3 <- as.numeric(scan("/home/si14w/gnearline/hdc/txt/E37_H3K4me3_coverage_across_genome.txt", what="numeric"))
cov_H3K4me3[cov_H3K4me3 > 0.1] <- 0.1
cov_H3K4me1 <- as.numeric(scan("/home/si14w/gnearline/hdc/txt/E37_H3K4me1_coverage_across_genome.txt", what="numeric"))
cov_H3K4me1[cov_H3K4me1 > 0.1] <- 0.1

pdf("/home/si14w/gnearline/hdc/results/coveragePlots.pdf")
plot(density(cov_input), main="input", xlab="coverage", ylab="frequency", col="red", xlim=c(0,0.5), ylim=c(0,50))
lines(density(cov_IgG), main="IgG", col="blue")
lines(density(cov_H3K4me3), main="H3K4me3", col="green")
lines(density(cov_H3K4me1), main="H3K4me1", col="orange")
legend("topright", c("input","IgG","H3K4me3","H3K4me1"), col=c("red","blue","green","orange"),bty="n", lty="solid")
plot(cov_input, main="input", ylab="coverage", pch=".")
dev.off()
