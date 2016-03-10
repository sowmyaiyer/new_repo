awk '{ if (($3 -$2) >= 200) print }'  ../txt/hg19_ucsc_genes_and_flu.bed >  ../txt/hg19_ucsc_genes_and_flu.longer_than_200.bed
awk '{ if ($6 == "-") {start=$3-150;end=$3} else if ($6 == "+") {start=$2;end=$2+150} printf("%s\t%d\t%d\t%s\t%s\t%s\n", $1,start,end,$4,$5,$6) }' ../txt/hg19_ucsc_genes_and_flu.longer_than_200.bed > ../txt/hg19_ucsc_genes_and_flu.5prime_plus150.bed
for time in {"A","B","C","D","E","F","G","H"}
do
	echo """
#	samtools view -b -F 4  ../bowtie_out/112213flucap_${time}.sorted.bam | bedtools bamtobed -i stdin > ../bowtie_out/112213flucap_${time}.sorted.bamtobed.bed
	bedtools coverage -nonamecheck -a ../txt/hg19_ucsc_genes_and_flu.5prime_plus150.bed -b ../bowtie_out/112213flucap_${time}.sorted.bamtobed.bed -counts -F 1 -s > ../txt/112213flucap_${time}.readCounts_150bp.bed
	Rscript normalizeReadCounts.R $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.bed $HOME/gnearline/flu/txt/112213flucap_${time}.readCounts_150bp.rpm.bed
	""" > ../bsubFiles/capseq_readcounts_flu_${time}.bsub
done
