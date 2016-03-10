require("ggplot2")
output_dir_name <- commandArgs(TRUE)[1]
fpkm_enhancers_closest_gene <- as.data.frame(read.table(paste("../results/",output_dir_name,"/fpkm_genes_closest_to_atac_enhancers.txt",sep="")))
fpkm_nonenhancers_closest_gene <- as.data.frame(read.table(paste("../results/",output_dir_name,"/fpkm_genes_closest_to_atac_nonenhancers.txt",sep="")))
fpkm_random_closest_gene <- as.data.frame(read.table(paste("../results/",output_dir_name,"/fpkm_genes_closest_to_random_distal_regions.txt",sep="")))
fpkm_all_closest <- rbind(fpkm_enhancers_closest_gene,fpkm_nonenhancers_closest_gene,fpkm_random_closest_gene)
colnames(fpkm_all_closest) <- c("type","fpkm")
lab_enhancer_all <- paste("enhancer(",nrow(fpkm_enhancers_closest_gene),")",sep="")
lab_nonenhancer_all <- paste("nonenhancer(",nrow(fpkm_nonenhancers_closest_gene),")",sep="")
lab_random_all <- paste("random(",nrow(fpkm_random_closest_gene),")",sep="")


fpkm_enhancers_closest_gene_within5k <- as.data.frame(read.table(paste("../results/",output_dir_name,"/fpkm_genes_closest_to_atac_enhancers_within5k.txt",sep="")))
fpkm_nonenhancers_closest_gene_within5k <- as.data.frame(read.table(paste("../results/",output_dir_name,"/fpkm_genes_closest_to_atac_nonenhancers_within5k.txt",sep="")))
fpkm_random_closest_gene_within5k <- as.data.frame(read.table(paste("../results/",output_dir_name,"/fpkm_genes_closest_to_random_distal_regions_within5k.txt",sep="")))

fpkm_all_closest_5k <- rbind(fpkm_enhancers_closest_gene_within5k, fpkm_nonenhancers_closest_gene_within5k,fpkm_random_closest_gene_within5k)
colnames(fpkm_all_closest_5k) <-  c("type","fpkm")
lab_enhancer_5k <- paste("enhancer(",nrow(fpkm_enhancers_closest_gene_within5k),")",sep="")
lab_nonenhancer_5k <- paste("nonenhancer(",nrow(fpkm_nonenhancers_closest_gene_within5k),")",sep="")
lab_random_5k <- paste("random(",nrow(fpkm_random_closest_gene_within5k),")",sep="")


fpkm_enhancers_closest_gene_within10k <- as.data.frame(read.table(paste("../results/",output_dir_name,"/fpkm_genes_closest_to_atac_enhancers_within10k.txt",sep="")))
fpkm_nonenhancers_closest_gene_within10k <- as.data.frame(read.table(paste("../results/",output_dir_name,"/fpkm_genes_closest_to_atac_nonenhancers_within10k.txt",sep="")))
fpkm_random_closest_gene_within10k <- as.data.frame(read.table(paste("../results/",output_dir_name,"/fpkm_genes_closest_to_random_distal_regions_within10k.txt",sep="")))

fpkm_all_closest_10k <- rbind(fpkm_enhancers_closest_gene_within10k, fpkm_nonenhancers_closest_gene_within10k,fpkm_random_closest_gene_within10k)
colnames(fpkm_all_closest_10k) <-  c("type","fpkm")
lab_enhancer_10k <- paste("enhancer(",nrow(fpkm_enhancers_closest_gene_within10k),")",sep="")
lab_nonenhancer_10k <- paste("nonenhancer(",nrow(fpkm_nonenhancers_closest_gene_within10k),")",sep="")
lab_random_10k <- paste("random(",nrow(fpkm_random_closest_gene_within10k),")",sep="")


pdf(paste("../results/",output_dir_name,"/expression_closest_genes.pdf",sep=""), onefile=TRUE)
ggplot(fpkm_all_closest, aes(x=type, y=fpkm, fill=type)) + geom_boxplot() + coord_cartesian(ylim = quantile(fpkm_all_closest$fpkm, c(0.1, 0.9))) + ylab("fpkm of closest transcript") + scale_x_discrete(labels=c(lab_enhancer_all,lab_nonenhancer_all,lab_random_all))
#ggplot(fpkm_all_5k, aes(x=type, y=fpkm, fill=type)) + geom_boxplot() + coord_cartesian(ylim = quantile(fpkm_all_5k$fpkm, c(0.1, 0.9))) + ylab("fpkm of transcript within 5K of ATAC peak")
ggplot(fpkm_all_closest_10k, aes(x=type, y=fpkm, fill=type)) + geom_boxplot() + coord_cartesian(ylim = quantile(fpkm_all_closest_10k$fpkm, c(0.1, 0.9))) + ylab("fpkm of closest transcript within 10K") + scale_x_discrete(labels=c(lab_enhancer_10k,lab_nonenhancer_10k,lab_random_10k))
ggplot(fpkm_all_closest_5k, aes(x=type, y=fpkm, fill=type)) + geom_boxplot() + coord_cartesian(ylim = quantile(fpkm_all_closest_5k$fpkm, c(0.1, 0.9))) + ylab("fpkm of closest transcript within 5K") + scale_x_discrete(labels=c(lab_enhancer_5k,lab_nonenhancer_5k,lab_random_5k))
dev.off()
