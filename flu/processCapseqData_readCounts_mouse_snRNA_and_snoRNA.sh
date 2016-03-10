awk '{ if ($3 == "-") {start=$5-150;end=$5} else if ($3 == "+") {start=$4;end=$4+150} printf("%s\t%d\t%d\t%s_%s_%s\t0\t%s\n", $2,start,end,$1,$6,$8,$3) }' mm9_snRNA_and_snoRNA-2.txt > mm9_snRNA_and_snoRNA.bed
for time in {}
do
	echo """
	samtools view -b -F 4  <bamFile> | bedtools bamtobed -i stdin > <bamFile>bamtobed.bed
	bedtools coverage -nonamecheck -a ../txt/hg19_ucsc_genes_and_flu.5prime_plus150.allRNA.bed -b ../bowtie_out/112213flucap_${time}.sorted.bamtobed.bed -counts -F 1 -s > ../txt/112213flucap_${time}.readCounts_150bp.allRNA.bed
	Rscript normalizeReadCounts.R $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.allRNA.bed $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.rpm.allRNA.bed
	""" > ../bsubFiles/capseq_readcounts_flu_${time}.bsub
done
