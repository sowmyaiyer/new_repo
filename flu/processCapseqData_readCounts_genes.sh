for time in {"A","B","C","D","E","F","G","H"}
do
	echo """
#	samtools view -b -F 4  ../bowtie_out/112213flucap_${time}.sorted.bam | bedtools bamtobed -i stdin > ../bowtie_out/112213flucap_${time}.sorted.bamtobed.bed
	bedtools intersect -nonamecheck -a ../txt/hg19_genes_150.bed -b ../bowtie_out/112213flucap_${time}.sorted.bamtobed.bed -c -F 1 -s > ../txt/112213flucap_${time}.readCounts_150bp.allRNA.genes.bed
	Rscript normalizeReadCounts.R $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.allRNA.genes.bed $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.rpm.allRNA.genes.bed
	""" > ../bsubFiles/capseq_readcounts_flu_${time}.bsub
done
